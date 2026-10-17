import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import signal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import luna15 as L

x = L.load(sys.argv[1])
out = sys.argv[2]
os.makedirs(out, exist_ok=True)
TR, CO = '/System/Library/Fonts/Supplemental/Times New Roman.ttf', '/System/Library/Fonts/Supplemental/Courier New.ttf'
F, C = (lambda s: ImageFont.truetype(TR, s)), (lambda s: ImageFont.truetype(CO, s))
INK, GREY, RED, BLUE = (0, 0, 0), (110, 110, 110), (190, 30, 30), (30, 80, 190)

fits = []
for n, (a, b) in enumerate(L.PIECES):
    t, fr, snr = L.track(x, a - 0.3, b + 1.5 if n == 2 else b + 0.2)
    fits.append(L.fit_piece(t, fr, snr, a, b)[0])
F_IMG = [L.mirror_spacing(x, fits[n], *L.PIECES[n])[0] for n in range(3)]


def mmss(s):
    return f'{int(s // 60)}:{s % 60:04.1f}' if s % 1 else f'{int(s // 60)}:{int(s % 60):02d}'


def grey_spec(P, lo, hi):
    v = np.clip((P - lo) / (hi - lo), 0, 1)
    return Image.fromarray((255 * (1 - v) ** 1.3).astype(np.uint8)[::-1]).convert('RGB')


t0, t1, f0, f1 = 189.0, 286.0, 400.0, 2700.0
f, t, Z = signal.stft(x[int(t0 * L.FS):int(t1 * L.FS)], L.FS, window='hann', nperseg=4096, noverlap=4096 - 441,
                      boundary=None, padded=False)
m = (f > f0) & (f < f1)
P = 10 * np.log10(np.abs(Z[m]) ** 2 + 1e-12)
P -= np.median(P, axis=1, keepdims=True)
W, H, X0, Y0 = 1300, 460, 70, 60
im = Image.new('RGB', (W + X0 + 30, H + Y0 + 110), (255, 255, 255))
im.paste(grey_spec(P, 3, 24).resize((W, H)), (X0, Y0))
d = ImageDraw.Draw(im)
d.rectangle([X0, Y0, X0 + W, Y0 + H], outline=INK)
X = lambda s: X0 + int((s - t0) / (t1 - t0) * W)
Yf = lambda fr: Y0 + H - int((fr - f0) / (f1 - f0) * H)
for s in range(190, 286, 10):
    d.line([(X(s), Y0 + H), (X(s), Y0 + H + 6)], fill=INK)
    d.text((X(s) - 14, Y0 + H + 8), mmss(s), font=F(15), fill=INK)
for fr in range(500, 2701, 500):
    d.line([(X0 - 6, Yf(fr)), (X0, Yf(fr))], fill=INK)
    d.text((X0 - 48, Yf(fr) - 8), f'{fr}', font=F(15), fill=INK)
d.text((8, Y0 - 22), 'Hz', font=F(15), fill=INK)
for n, (a, b) in enumerate(L.PIECES):
    d.text((X((a + b) / 2) - 22, Y0 + 8), f'piece {n + 1}', font=F(17), fill=BLUE)
for g0, g1 in [(233.0, 234.45), (252.55, 253.65)]:
    d.text((X(g0) - 12, Y0 + H + 28), 'edit', font=F(14), fill=GREY)
d.line([(X(198.1), Y0), (X(198.1), Y0 + 30)], fill=GREY)
d.text((X(190.5), Y0 + 32), 'orbit tone', font=F(14), fill=GREY)
d.line([(X(259.5), Y0 + H - 40), (X(259.5), Y0 + H)], fill=RED)
d.text((X(259.5) - 60, Y0 + H - 60), '"Time is 15:50, 20 seconds"', font=F(15), fill=RED)
d.line([(X(277.6), Y0 + 40), (X(277.6), Y0 + H)], fill=RED)
d.text((X(277.6) + 4, Y0 + 44), 'signal stops', font=F(15), fill=RED)
d.text((X(277.6) + 4, Y0 + 62), '4:37.6', font=F(15), fill=RED)
d.text((X0, 12), 'Jodrell Bank tape, 21 July 1969: Luna 15\'s 922.8 MHz carrier during the descent (audio beat note). '
       'Faint rising lines are receiver images.', font=F(18), fill=INK)
d.text((X0, Y0 + H + 52), 'Tape time (min:s). The falling line is the carrier; it falls because the braking engine changes '
       'the probe\'s speed along our line of sight.', font=F(15), fill=GREY)
im.save(f'{out}/descent.png')

W2, H2, X0, Y0 = 900, 300, 80, 50
im = Image.new('RGB', (W2 + X0 + 40, H2 + Y0 + 90), (255, 255, 255))
d = ImageDraw.Draw(im)
ta, tb, alo, ahi = 196.0, 280.0, 0.0, 4.5
X = lambda s: X0 + int((s - ta) / (tb - ta) * W2)
Ya = lambda v: Y0 + H2 - int((v - alo) / (ahi - alo) * H2)
d.rectangle([X0, Y0, X0 + W2, Y0 + H2], outline=INK)
for v in np.arange(0, 4.51, 1.0):
    d.line([(X0 - 6, Ya(v)), (X0, Ya(v))], fill=INK)
    d.text((X0 - 30, Ya(v) - 8), f'{v:.0f}', font=F(15), fill=INK)
for s in range(200, 281, 10):
    d.line([(X(s), Y0 + H2), (X(s), Y0 + H2 + 6)], fill=INK)
    d.text((X(s) - 14, Y0 + H2 + 8), mmss(s), font=F(15), fill=INK)
g_los = 1.62 * np.sin(np.radians(23.1))
d.line([(X0, Ya(g_los)), (X0 + W2, Ya(g_los))], fill=GREY, width=1)
d.text((X(200), Ya(g_los) - 20), f'gravity alone (no engine): {g_los:.2f} m/s²', font=F(15), fill=GREY)
for n, (a, b) in enumerate(L.PIECES):
    dp = np.polyder(fits[n])
    ss = np.linspace(a, b, 200)
    pts = [(X(s), Ya(-np.polyval(dp, s - a) * L.LAM)) for s in ss]
    d.line(pts, fill=RED, width=3)
d.text((X0, 12), 'Acceleration along the line of sight, from the slope of the carrier (m/s²; slope × 0.325 m)',
       font=F(18), fill=INK)
d.text((X0, Y0 + H2 + 34), 'Tape time. Gaps are edits in the tape. If the tape plays 1.8% slow (see the hum), add 3.6%.',
       font=F(15), fill=GREY)
im.save(f'{out}/acceleration.png')

a, b = L.PIECES[2]
p = fits[2]
t0, t1 = 253.8, 280.5
z, ts = L.dechirp(x, p, a, b, t0, t1)
f, t, Z = signal.stft(z, L.FS, window='hann', nperseg=2048, noverlap=2048 - 110, return_onesided=False,
                      boundary=None, padded=False)
f = np.fft.fftshift(f)
Z = np.fft.fftshift(Z, axes=0)
t = t + t0
SPAN = 400
m = np.abs(f) <= SPAN
P = 10 * np.log10(np.abs(Z[m]) ** 2 + 1e-12)
fn, tn, Zn = signal.stft(x[int(L.POST_LOS[0] * L.FS):int(L.POST_LOS[1] * L.FS)], L.FS, window='hann', nperseg=2048,
                         noverlap=2048 - 441, boundary=None, padded=False)
P -= np.median(P[:, t > 279.0], axis=1, keepdims=True)
P = signal.convolve2d(P, np.ones((3, 9)) / 27, mode='same')
W3, H3, X0, Y0 = 1100, 420, 70, 50
im = Image.new('RGB', (W3 + X0 + 160, H3 + Y0 + 80), (255, 255, 255))
im.paste(grey_spec(P, 1.5, 13).resize((W3, H3)), (X0, Y0))
d = ImageDraw.Draw(im)
d.rectangle([X0, Y0, X0 + W3, Y0 + H3], outline=INK)
X = lambda s: X0 + int((s - t0) / (t1 - t0) * W3)
Yo = lambda o: Y0 + H3 // 2 - int(o / SPAN * H3 / 2)
for s in range(255, 281, 5):
    d.line([(X(s), Y0 + H3), (X(s), Y0 + H3 + 6)], fill=INK)
    d.text((X(s) - 14, Y0 + H3 + 8), mmss(s), font=F(15), fill=INK)
for o in range(-400, 401, 100):
    d.line([(X0 - 6, Yo(o)), (X0, Yo(o))], fill=INK)
    d.text((X0 - 44, Yo(o) - 8), f'{o:+d}' if o else '0', font=F(15), fill=INK)
d.text((X0 + W3 + 8, Yo(0) - 9), 'carrier', font=F(16), fill=INK)
d.text((X0 + W3 + 8, Yo(170) - 9), 'one-sided', font=F(16), fill=RED)
d.text((X0 + W3 + 8, Yo(170) + 9), 'band', font=F(16), fill=RED)
d.text((X0 + W3 + 8, Yo(-170) - 9), 'no band below', font=F(16), fill=GREY)
d.text((X0 + W3 + 8, Yo(-170) + 9), '(thin lines: hum)', font=F(16), fill=GREY)
tt = np.linspace(t0, b, 800)
fm = L.model_freq(p, a, b, tt)
for k in (2, 3):
    D = k * F_IMG[2] - 2 * fm
    pts = [(X(s), Yo(o)) for s, o in zip(tt, D) if abs(o) < SPAN - 5]
    for i in range(0, len(pts) - 6, 12):
        d.line(pts[i:i + 6], fill=BLUE, width=2)
d.text((X(256.2), Yo(330)), 'receiver image line', font=F(15), fill=BLUE)
d.line([(X(L.CLEAN_FROM), Y0 + H3), (X(L.CLEAN_FROM), Y0 + H3 + 20)], fill=BLUE)
d.text((X(L.CLEAN_FROM) + 4, Y0 + H3 + 20), 'image-free from here', font=F(14), fill=BLUE)
d.line([(X(277.6), Y0), (X(277.6), Y0 + H3)], fill=RED)
d.text((X(277.6) + 4, Y0 + 4), 'signal stops', font=F(15), fill=RED)
d.text((X0, 12), 'The last 24 seconds with the carrier moved to 0 Hz (offset in Hz from the carrier). '
       'Shading: power over the receiver noise measured after the signal stopped.', font=F(17), fill=INK)
d.text((X0, Y0 + H3 + 34), 'Tape time. Dotted blue: predicted receiver image lines (masked in the measurements).',
       font=F(15), fill=GREY)
im.save(f'{out}/last-piece.png')

rows = []
for n, (a, b) in enumerate(L.PIECES):
    r, _, _, _ = L.band_test(x, fits[n], a, b, F_IMG[n], noise_ref=L.POST_LOS if n == 2 else None)
    rows.append(r)
W4, H4, X0, Y0 = 900, 260, 80, 50
im = Image.new('RGB', (W4 + X0 + 40, H4 + Y0 + 90), (255, 255, 255))
d = ImageDraw.Draw(im)
ta, tb, vlo, vhi = 196.0, 280.0, -6.0, 12.0
X = lambda s: X0 + int((s - ta) / (tb - ta) * W4)
Yv = lambda v: Y0 + H4 - int((v - vlo) / (vhi - vlo) * H4)
d.rectangle([X0, Y0, X0 + W4, Y0 + H4], outline=INK)
d.line([(X0, Yv(0)), (X0 + W4, Yv(0))], fill=GREY)
for v in range(-6, 13, 3):
    d.line([(X0 - 6, Yv(v)), (X0, Yv(v))], fill=INK)
    d.text((X0 - 34, Yv(v) - 8), f'{v:+d}' if v else '0', font=F(15), fill=INK)
for s in range(200, 281, 10):
    d.line([(X(s), Y0 + H4), (X(s), Y0 + H4 + 6)], fill=INK)
    d.text((X(s) - 14, Y0 + H4 + 8), mmss(s), font=F(15), fill=INK)
for r in rows:
    for q in r:
        cx, cy = X(q['t_tape']), Yv(np.clip(q['above_vs_below_dB'], vlo, vhi))
        d.ellipse([cx - 4, cy - 4, cx + 4, cy + 4], fill=RED if q['t_tape'] > 253 else INK)
d.text((X0, 12), 'Power 40-300 Hz above the carrier minus 40-300 Hz below it (dB, one dot per second)', font=F(18),
       fill=INK)
d.text((X0, Y0 + H4 + 34), 'Tape time. Pieces 1-2 sit at about 0 dB; the last piece (red) sits 4-10 dB above.',
       font=F(15), fill=GREY)
im.save(f'{out}/above-below.png')

tl, lev = L.carrier_level(x, p, a, b, 274.6, 278.6)
W5, H5, X0, Y0 = 900, 240, 80, 50
im = Image.new('RGB', (W5 + X0 + 40, H5 + Y0 + 90), (255, 255, 255))
d = ImageDraw.Draw(im)
ta, tb, vlo, vhi = 274.6, 278.6, -10.0, 35.0
X = lambda s: X0 + int((s - ta) / (tb - ta) * W5)
Yv = lambda v: Y0 + H5 - int((np.clip(v, vlo, vhi) - vlo) / (vhi - vlo) * H5)
d.rectangle([X0, Y0, X0 + W5, Y0 + H5], outline=INK)
for v in range(-10, 36, 10):
    d.line([(X0 - 6, Yv(v)), (X0, Yv(v))], fill=INK)
    d.text((X0 - 34, Yv(v) - 8), f'{v}', font=F(15), fill=INK)
for s in np.arange(275.0, 278.61, 0.5):
    d.line([(X(s), Y0 + H5), (X(s), Y0 + H5 + 6)], fill=INK)
    d.text((X(s) - 16, Y0 + H5 + 8), mmss(s), font=F(15), fill=INK)
d.line([(X(s_), Yv(v_)) for s_, v_ in zip(tl, lev)], fill=INK, width=2)
d.line([(X(277.6), Y0), (X(277.6), Y0 + H5)], fill=RED)
d.text((X(277.6) + 4, Y0 + 6), '4:37.6', font=F(15), fill=RED)
d.text((X0, 12), 'Carrier level in an 8 Hz band around its track (dB over the noise beside it)', font=F(18), fill=INK)
d.text((X0, Y0 + H5 + 34), 'Tape time. The tone fades and flutters for about a second, then stops within a few hundredths '
       'of a second.', font=F(15), fill=GREY)
im.save(f'{out}/last-seconds.png')
print('figures saved to', out)
