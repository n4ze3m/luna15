Luna 15 descent on Jodrell Bank's "Luna 15 and Apollo 11" tape (21 July 1969): extracted data
Produced by luna15.py from Luna15-Apollo11.mp3 (https://www.jb.man.ac.uk/multimedia/audio/luna15-apollo11.html,
(c) Jodrell Bank Centre for Astrophysics, University of Manchester). The audio itself is not included.
All times are seconds into that MP3 (tape time). Audio is resampled to 22050 Hz.

carrier.csv   the carrier (audio beat note) through the three descent pieces
  piece       1, 2 or 3 (tape 198.4-232.85, 234.6-252.45, 253.8-277.4 s)
  t_tape_s    centre of the spectrum (0.37 s window, 20 ms steps)
  f_hz        carrier frequency (peak interpolation)
  snr_db      peak over the median of 600-1700 Hz in that spectrum

band.csv      one row per second: is there more power above the carrier than below it?
  piece                one of the three pieces
  t_tape               centre of the 1 s window
  above_vs_below_dB    mean power 40..300 Hz above the carrier minus 40..300 Hz below it, after dividing every
                       frequency by the receiver noise (piece 3: noise measured after the signal stops,
                       278.5-305 s; pieces 1-2: each bin's median while the carrier is >350 Hz away).
                       Receiver image lines (k * ~786 Hz - carrier, k = 1..7) are masked within 45 Hz.
  above_dB_re_noise    mean above, dB over that noise
  below_dB_re_noise    mean below, dB over that noise
  offset_q10/q50/q90   10/50/90 % points of the excess above the carrier (Hz above the carrier)

summary.json  piece fits (frequency, slope in Hz/s and the line-of-sight acceleration = slope x 0.3249 m at
              0/25/50/75/100 % of each piece), the loss of signal (tape time, level, slope at the end, the
              spoken time check), receiver image-line spacing, the band test (with the image-free part of
              piece 3 and a carrier-free control), Earth elevation at the published crash sites, closing speed
              per 100 Hz for flat ground, in-orbit drifts from earlier on the tape, and the hum frequency.
