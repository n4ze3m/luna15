import csv
import json
import os
import subprocess
import sys

import numpy as np
from scipy import signal

FS = 22050
C = 299792458.0
F_RF = 922.8e6
LAM = C / F_RF
PIECES = [(198.4, 232.85), (234.6, 252.45), (253.8, 277.4)]
POST_LOS = (278.5, 305.0)
CLEAN_FROM = 260.5
ORBIT = [(88, 130, 2050, 2180, '20 Jul, 40th orbit, just after emerging from behind the Moon'),
         (132, 190, 2000, 2120, '21 Jul, start of the 52nd orbit')]
TIME_CALL = {'tape_s': [258.4, 260.6], 'said': '15:50, 20 seconds'}
SITES = {'17N 60E': (17, 60), '16N 57E': (16, 57), '12N 60E': (12, 60)}
SUB_OBS = (3.548, 353.060)


def load(path):
    raw = subprocess.run(['ffmpeg', '-v', 'quiet', '-i', path, '-f', 'f32le', '-ac', '1', '-ar', str(FS), '-'],
                         capture_output=True, check=True).stdout
    return np.frombuffer(raw, dtype=np.float32).astype(np.float64)


def quiet_gaps(x, thr_db=-50, min_s=0.15):
    env = np.sqrt(np.maximum(signal.sosfiltfilt(signal.butter(2, 20, fs=FS, output='sos'), x * x), 1e-12))
    lev = 20 * np.log10(env[::int(FS * 0.05)])
    q = lev < thr_db
    out, i = [], 0
    while i < len(q):
        if q[i]:
            j = i
            while j < len(q) and q[j]:
                j += 1
            if (j - i) * 0.05 >= min_s:
                out.append((round(i * 0.05, 2), round(j * 0.05, 2)))
            i = j
        else:
            i += 1
    return out


def track(x, a, b, lo=600, hi=1700, nper=8192, hop=441, win=25):
    f, t, Z = signal.stft(x[int(a * FS):int(b * FS)], FS, window='hann', nperseg=nper, noverlap=nper - hop,
                          boundary=None, padded=False)
    P = np.abs(Z) ** 2
    t = t + a
    df = f[1] - f[0]
    band = (f > lo) & (f < hi)
    i0 = np.argmax(P[band].max(axis=0))
    f_start = f[np.flatnonzero(band)[np.argmax(P[band, i0])]]
    fr = np.full(len(t), np.nan)
    snr = np.full(len(t), np.nan)
    for step in (1, -1):
        fc, i = f_start, i0
        while 0 <= i < len(t):
            m = np.abs(f - fc) < win
            k = np.flatnonzero(m)[np.argmax(P[m, i])]
            y0, y1, y2 = np.log(P[k - 1:k + 2, i] + 1e-30)
            fr[i] = f[k] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * df
            snr[i] = 10 * np.log10(P[k, i] / np.median(P[band, i]))
            if snr[i] < 6:
                break
            fc = fr[i]
            i += step
    return t, fr, snr


def fit_piece(t, fr, snr, a, b):
    g = (t >= a) & (t <= b) & (snr > 8) & np.isfinite(fr)
    p = np.polyfit(t[g] - a, fr[g], 3)
    res = fr[g] - np.polyval(p, t[g] - a)
    return p, float(res.std())


def model_freq(p, a, b, ts):
    dp = np.polyder(p)
    return (np.polyval(p, np.clip(ts, a, b) - a) + (ts > b) * np.polyval(dp, b - a) * (ts - b)
            + (ts < a) * np.polyval(dp, 0) * (ts - a))


def dechirp(x, p, a, b, t0, t1):
    seg = x[int(t0 * FS):int(t1 * FS)]
    ts = np.arange(len(seg)) / FS + t0
    return signal.hilbert(seg) * np.exp(-2j * np.pi * np.cumsum(model_freq(p, a, b, ts)) / FS), ts


def carrier_level(x, p, a, b, t0, t1, bw=8.0):
    z, ts = dechirp(x, p, a, b, t0, t1)
    sos = signal.butter(4, bw, fs=FS, output='sos')
    lp = lambda u: signal.sosfiltfilt(sos, u.real) + 1j * signal.sosfiltfilt(sos, u.imag)
    car = np.abs(lp(z))
    ref = np.abs(lp(z * np.exp(-2j * np.pi * 150 * (ts - t0))))
    return ts[::110], 20 * np.log10(car[::110] / np.median(ref))


def mirror_spacing(x, p, a, b):
    f, t, Z = signal.stft(x[int(a * FS):int(b * FS)], FS, window='hann', nperseg=16384, noverlap=16384 - 441,
                          boundary=None, padded=False)
    t = t + a
    L = np.log10(np.abs(Z) ** 2 / np.median(np.abs(Z) ** 2, axis=0, keepdims=True) + 1e-12)
    fm = np.polyval(p, t - a)
    df = f[1] - f[0]
    Fs = np.arange(780, 792, 0.05)
    score = []
    for F in Fs:
        s = []
        for k in (3, 4, 5, 6):
            kk = np.round((k * F - fm) / df).astype(int)
            ok = (kk > 2) & (kk < len(f) - 2)
            s.append(L[kk[ok], np.flatnonzero(ok)].mean())
        score.append(np.mean(s))
    score = np.array(score)
    return float(Fs[np.argmax(score)]), float((score.max() - np.median(score)) / score.std())


def band_test(x, p, a, b, F_img, noise_ref=None, step=1.0):
    nper, hop = 2048, 441
    f, t, Z = signal.stft(x[int(a * FS):int(b * FS)], FS, window='hann', nperseg=nper, noverlap=nper - hop,
                          boundary=None, padded=False)
    P = np.abs(Z) ** 2
    t = t + a
    fm = model_freq(p, a, b, t)
    df = f[1] - f[0]
    if noise_ref is None:
        far = np.abs(f[:, None] - fm[None, :]) > 350
        N = np.array([np.median(P[i, far[i]]) if far[i].sum() > 20 else np.median(P[i]) for i in range(len(f))])
    else:
        fn, tn, Zn = signal.stft(x[int(noise_ref[0] * FS):int(noise_ref[1] * FS)], FS, window='hann', nperseg=nper,
                                 noverlap=nper - hop, boundary=None, padded=False)
        N = np.median(np.abs(Zn) ** 2, axis=1)
    offs = np.arange(-320, 321, df)
    G = np.full((len(offs), len(t)), np.nan)
    for j in range(len(t)):
        k = np.round((fm[j] + offs) / df).astype(int)
        ok = (k > 0) & (k < len(f))
        G[ok, j] = P[k[ok], j] / N[k[ok]]
        for kk in range(1, 8):
            D = kk * F_img - 2 * fm[j]
            G[(np.abs(offs - D) < 45) | (np.abs(offs + D) < 45), j] = np.nan
    pos = (offs > 40) & (offs < 300)
    neg = (offs < -40) & (offs > -300)
    rows = []
    for tb in np.arange(a, b - step + 1e-9, step):
        jj = (t >= tb) & (t < tb + step)
        gp, gn = G[pos][:, jj], G[neg][:, jj]
        if np.isfinite(gp).sum() < 40 or np.isfinite(gn).sum() < 40:
            continue
        ratio = 10 * np.log10(np.nanmean(gp) / np.nanmean(gn))
        ex = np.nanmean(G[pos][:, jj], axis=1) - np.nanmean(G[neg][:, jj], axis=1)[::-1][:pos.sum()]
        ex = np.clip(np.nan_to_num(ex), 0, None)
        q = [np.nan] * 3
        if ex.sum() > 0:
            c = np.cumsum(ex) / ex.sum()
            q = [float(offs[pos][min(np.searchsorted(c, v), pos.sum() - 1)]) for v in (0.1, 0.5, 0.9)]
        rows.append({'t_tape': round(tb + step / 2, 2), 'above_vs_below_dB': round(ratio, 2),
                     'above_dB_re_noise': round(10 * np.log10(np.nanmean(gp)), 2),
                     'below_dB_re_noise': round(10 * np.log10(np.nanmean(gn)), 2),
                     'offset_q10': q[0], 'offset_q50': q[1], 'offset_q90': q[2]})
    return rows, offs, G, t


def orbit_drift(x, a, b, lo, hi):
    t, fr, snr = track(x, a, b, lo=lo, hi=hi, nper=16384, hop=2205, win=15)
    g = snr > 10
    p = np.polyfit(t[g], fr[g], 2)
    return {'f_hz': round(float(np.median(fr[g])), 1), 'drift_hz_s': round(float(np.polyfit(t[g], fr[g], 1)[0]), 3),
            'curvature_hz_s2': round(float(2 * p[0]), 4)}


def hum(x, a, b, f0=98.2):
    seg = x[int(a * FS):int(b * FS)]
    f, P = signal.welch(seg, FS, nperseg=65536)
    m = (f > f0 - 2) & (f < f0 + 2)
    return float(f[m][np.argmax(P[m])])


def earth_elevation(lat, lon):
    la, lo = np.radians([lat, lon])
    sa, so = np.radians(SUB_OBS)
    cosd = np.sin(la) * np.sin(sa) + np.cos(la) * np.cos(sa) * np.cos(lo - so)
    return float(90 - np.degrees(np.arccos(cosd)))


def main():
    x = load(sys.argv[1])
    outdir = sys.argv[2] if len(sys.argv) > 2 else 'out'
    os.makedirs(outdir, exist_ok=True)
    summary = {'audio_s': round(len(x) / FS, 2), 'quiet_gaps_tape_s': quiet_gaps(x), 'wavelength_m': round(LAM, 5)}

    fits, rows_car, pieces = [], [], []
    for n, (a, b) in enumerate(PIECES):
        t, fr, snr = track(x, a - 0.3, b + 1.5 if n == 2 else b + 0.2)
        p, rms = fit_piece(t, fr, snr, a, b)
        fits.append(p)
        dp = np.polyder(p)
        sl = [float(np.polyval(dp, u)) for u in np.linspace(0, b - a, 5)]
        pieces.append({'tape_s': [a, b], 'f_start_hz': round(float(np.polyval(p, 0)), 1),
                       'f_end_hz': round(float(np.polyval(p, b - a)), 1),
                       'slope_hz_s_at_0_25_50_75_100pc': [round(v, 2) for v in sl],
                       'los_accel_m_s2_at_0_25_50_75_100pc': [round(-v * LAM, 2) for v in sl],
                       'fit_rms_hz': round(rms, 2)})
        for ti, fi, si in zip(t, fr, snr):
            if np.isfinite(fi) and a <= ti <= b + 0.3:
                rows_car.append([n + 1, round(float(ti), 3), round(float(fi), 2), round(float(si), 1)])
    summary['pieces'] = pieces

    a, b = PIECES[2]
    tl, lev = carrier_level(x, fits[2], a, b, 275.0, 279.0)
    above = lev > 6
    i_last = np.flatnonzero(above & (tl < 278.5))[-1]
    t_los = float(tl[i_last])
    pre = float(np.median(lev[(tl > 275.5) & (tl < 276.5)]))
    summary['loss_of_signal'] = {
        'tape_s': round(t_los, 3), 'tape_mmss': f'{int(t_los // 60)}:{t_los % 60:06.3f}',
        'level_db_275.5_276.5': round(pre, 1),
        'level_db_last_0.5s': round(float(np.median(lev[(tl > t_los - 0.5) & (tl <= t_los)])), 1),
        'slope_hz_s_at_end': pieces[2]['slope_hz_s_at_0_25_50_75_100pc'][-1],
        'time_call': TIME_CALL,
        'tape_s_from_time_call': [round(t_los - TIME_CALL['tape_s'][1], 1), round(t_los - TIME_CALL['tape_s'][0], 1)],
    }

    summary['image_line_spacing_hz'] = [round(mirror_spacing(x, fits[n], *PIECES[n])[0], 2) for n in range(3)]

    band_rows = []
    summary['band'] = {}
    for n, (a, b) in enumerate(PIECES):
        rows, _, _, _ = band_test(x, fits[n], a, b, summary['image_line_spacing_hz'][n],
                                  noise_ref=POST_LOS if n == 2 else None)
        r = np.array([q['above_vs_below_dB'] for q in rows])
        summary['band'][f'piece{n + 1}'] = {'windows': len(r), 'mean_above_vs_below_dB': round(float(r.mean()), 2),
                                            'windows_over_3dB': int((r > 3).sum()), 'windows_below_0dB': int((r < 0).sum())}
        if n == 2:
            clean = [q for q in rows if q['t_tape'] > CLEAN_FROM]
            rc = np.array([q['above_vs_below_dB'] for q in clean])
            tq = np.array([q['t_tape'] for q in clean])
            q50 = np.array([q['offset_q50'] for q in clean])
            k = np.polyfit(tq, q50, 1)
            summary['band']['piece3_clean'] = {
                'from_tape_s': CLEAN_FROM, 'windows': len(rc), 'mean_above_vs_below_dB': round(float(rc.mean()), 2),
                'windows_over_3dB': int((rc > 3).sum()), 'windows_below_0dB': int((rc < 0).sum()),
                'centre_offset_hz_fit_at_start_end': [round(float(np.polyval(k, tq[0])), 0), round(float(np.polyval(k, tq[-1])), 0)]}
        for q in rows:
            band_rows.append({'piece': n + 1, **q})
    a, b = PIECES[2]
    ctrl = x.copy()
    shift = int((POST_LOS[0] + 0.5 - a) * FS)
    ctrl[int(a * FS):int(b * FS)] = x[int(a * FS) + shift:int(b * FS) + shift]
    rows, _, _, _ = band_test(ctrl, fits[2], a, b, summary['image_line_spacing_hz'][2], noise_ref=POST_LOS)
    r = np.array([q['above_vs_below_dB'] for q in rows])
    summary['band']['control_no_carrier'] = {'windows': len(r), 'mean_above_vs_below_dB': round(float(r.mean()), 2)}

    elev = {k: round(earth_elevation(*v), 1) for k, v in SITES.items()}
    summary['earth_elevation_deg_at_site'] = elev
    s_lo, s_hi = (np.sin(np.radians(min(elev.values()))), np.sin(np.radians(max(elev.values()))))
    summary['closing_speed_m_s_per_100hz'] = [round(100 * LAM / (2 * s_hi), 1), round(100 * LAM / (2 * s_lo), 1)]

    summary['orbit_segments'] = [{'tape_s': [a, b], 'what': w, **orbit_drift(x, a, b, lo, hi)} for a, b, lo, hi, w in ORBIT]
    summary['hum_hz'] = {'descent': round(hum(x, 191.8, 306.7), 2), 'orbit_20jul': round(hum(x, 88, 130), 2)}

    with open(os.path.join(outdir, 'carrier.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['piece', 't_tape_s', 'f_hz', 'snr_db'])
        w.writerows(rows_car)
    with open(os.path.join(outdir, 'band.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(band_rows[0].keys()))
        w.writeheader()
        w.writerows(band_rows)
    json.dump(summary, open(os.path.join(outdir, 'summary.json'), 'w'), indent=1, ensure_ascii=False)
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == '__main__':
    main()
