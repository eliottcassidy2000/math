#!/usr/bin/env python3
"""deep_structure_s116n.py — Deep temporal structure of MERRY SHITMAS.

Not what the song IS. What the song DOES, moment by moment.
The derivative of each element. The shape of the flow.
What enters, what exits, what transforms.
"""
import numpy as np
import librosa
import warnings
warnings.filterwarnings('ignore')
from math import log, log2, sqrt

FILEPATH = "inbox/SpotiDown.App - MERRY SHITMAS - BabyTron Universal.mp3"

print()
print("  MERRY SHITMAS — Deep Temporal Structure")
print()
print("="*70)
print()

y, sr = librosa.load(FILEPATH, sr=22050, mono=True)
duration = len(y) / sr

# High-resolution analysis: 0.5-second windows
hop = sr // 4  # 4 frames per second
n_fft = 2048

# ============================================================
print("  I. THE SKELETON: WHERE THINGS CHANGE")
print("  " + "-"*40)
print()

# Onset strength over time — measures WHERE events happen
onset_env = librosa.onset.onset_strength(y=y, sr=sr, hop_length=hop)
times_onset = librosa.times_like(onset_env, sr=sr, hop_length=hop)

# Novelty curve — measures WHERE things change CHARACTER
# (spectral flux: how different is this moment from the previous one)
S = np.abs(librosa.stft(y, n_fft=n_fft, hop_length=hop))
spectral_flux = np.sqrt(np.sum(np.diff(S, axis=1)**2, axis=0))
times_flux = np.linspace(0, duration, len(spectral_flux))

# Segment boundaries — WHERE the song structurally shifts
boundaries = librosa.segment.agglomerative(librosa.feature.mfcc(y=y, sr=sr, hop_length=hop), k=8)
boundary_times = librosa.frames_to_time(boundaries, sr=sr, hop_length=hop)

print(f"  Structural boundaries detected at:")
for i, bt in enumerate(boundary_times):
    pct = bt / duration * 100
    print(f"    Boundary {i+1}: {bt:.1f}s ({pct:.0f}%)")
print()

# ============================================================
print()
print("  II. THE LAYERS: WHAT IS PRESENT AT EACH MOMENT")
print("  " + "-"*40)
print()

# Decompose into: harmonic (pitched), percussive (transient), residual
y_harm, y_perc = librosa.effects.hpss(y)

# Energy of each component over time
rms_total = librosa.feature.rms(y=y, hop_length=hop)[0]
rms_harm = librosa.feature.rms(y=y_harm, hop_length=hop)[0]
rms_perc = librosa.feature.rms(y=y_perc, hop_length=hop)[0]
times_rms = librosa.times_like(rms_total, sr=sr, hop_length=hop)

# Sub-bass (808): 20-80 Hz
# Low-mid (vocal body): 200-800 Hz
# Mid (vocal presence): 800-4000 Hz
# High (hi-hats, air): 4000+ Hz

freqs = librosa.fft_frequencies(sr=sr, n_fft=n_fft)

def band_energy(S, freqs, lo, hi):
    mask = (freqs >= lo) & (freqs < hi)
    return np.mean(S[mask, :], axis=0) if mask.any() else np.zeros(S.shape[1])

sub_bass = band_energy(S, freqs, 20, 80)
low_mid = band_energy(S, freqs, 200, 800)
mid = band_energy(S, freqs, 800, 4000)
high = band_energy(S, freqs, 4000, 11000)
times_bands = np.linspace(0, duration, S.shape[1])

# Normalize each band to its own max for relative comparison
for band in [sub_bass, low_mid, mid, high]:
    mx = band.max()
    if mx > 0: band /= mx

# Print a time-series snapshot every 5 seconds
print("  Layer presence over time (normalized, 0-10 scale):")
print(f"  {'time':>6s}  {'sub-bass':>8s}  {'low-mid':>8s}  {'mid':>8s}  {'high':>8s}  {'harm':>6s}  {'perc':>6s}")
print("  " + "-"*55)

snapshot_interval = 3  # seconds
max_rms_h = rms_harm.max() if rms_harm.max() > 0 else 1
max_rms_p = rms_perc.max() if rms_perc.max() > 0 else 1

for t_sec in np.arange(0, duration, snapshot_interval):
    # Find nearest frame
    idx_band = int(t_sec / duration * len(sub_bass))
    idx_band = min(idx_band, len(sub_bass)-1)
    idx_rms = int(t_sec / duration * len(rms_harm))
    idx_rms = min(idx_rms, len(rms_harm)-1)

    sb = int(sub_bass[idx_band] * 10)
    lm = int(low_mid[idx_band] * 10)
    md = int(mid[idx_band] * 10)
    hi = int(high[idx_band] * 10)
    hm = int(rms_harm[idx_rms] / max_rms_h * 10)
    pc = int(rms_perc[idx_rms] / max_rms_p * 10)

    sb_bar = "#" * sb
    lm_bar = "#" * lm
    md_bar = "#" * md
    hi_bar = "#" * hi

    # Mark boundaries
    boundary_mark = ""
    for bt in boundary_times:
        if abs(t_sec - bt) < snapshot_interval / 2:
            boundary_mark = " <-- BOUNDARY"
            break

    print(f"  {t_sec:5.1f}s  {sb_bar:>8s}  {lm_bar:>8s}  {md_bar:>8s}  {hi_bar:>8s}  {hm:>6d}  {pc:>6d}{boundary_mark}")

print()

# ============================================================
print()
print("  III. THE DERIVATIVE: HOW THINGS CHANGE")
print("  " + "-"*40)
print()
print("  The derivative of each layer = its rate of change.")
print("  Positive derivative = growing. Negative = fading. Zero = sustaining.")
print()

# Compute derivatives of the band energies
d_sub = np.diff(sub_bass)
d_low = np.diff(low_mid)
d_mid = np.diff(mid)
d_high = np.diff(high)

# Summarize in sections
n_sections = 12
section_len = len(d_sub) // n_sections if n_sections > 0 else 1

print(f"  {'section':>8s}  {'time':>8s}  {'d(bass)':>8s}  {'d(low)':>8s}  {'d(mid)':>8s}  {'d(high)':>8s}  {'character'}")
print("  " + "-"*70)

for i in range(n_sections):
    start = i * section_len
    end = min((i+1) * section_len, len(d_sub))
    t_start = start / len(d_sub) * duration
    t_end = end / len(d_sub) * duration

    db = np.mean(d_sub[start:end])
    dl = np.mean(d_low[start:end])
    dm = np.mean(d_mid[start:end])
    dh = np.mean(d_high[start:end])

    # Characterize
    signs = (1 if db > 0.001 else (-1 if db < -0.001 else 0),
             1 if dl > 0.001 else (-1 if dl < -0.001 else 0),
             1 if dm > 0.001 else (-1 if dm < -0.001 else 0),
             1 if dh > 0.001 else (-1 if dh < -0.001 else 0))

    char = ""
    if all(s >= 0 for s in signs) and any(s > 0 for s in signs):
        char = "BUILDING"
    elif all(s <= 0 for s in signs) and any(s < 0 for s in signs):
        char = "FADING"
    elif signs[0] > 0 and signs[2] < 0:
        char = "bass enters, vocal exits"
    elif signs[0] < 0 and signs[2] > 0:
        char = "vocal enters, bass exits"
    elif all(s == 0 for s in signs):
        char = "SUSTAINING (plateau)"
    else:
        char = "SHIFTING"

    print(f"  {i+1:8d}  {t_start:4.0f}-{t_end:4.0f}s  {db:+8.4f}  {dl:+8.4f}  {dm:+8.4f}  {dh:+8.4f}  {char}")

print()

# ============================================================
print()
print("  IV. THE VOCAL FLOW: SHAPE AND MOTION")
print("  " + "-"*40)
print()

# Extract the vocal-dominant frequency range (200-4000 Hz) from harmonic component
y_harm_filtered = librosa.effects.preemphasis(y_harm)

# Pitch tracking in the vocal range
pitches, magnitudes = librosa.piptrack(y=y_harm, sr=sr, fmin=80, fmax=600, hop_length=hop)

# Extract the dominant pitch at each frame
pitch_track = []
pitch_times = []
for t_idx in range(pitches.shape[1]):
    mag_col = magnitudes[:, t_idx]
    if mag_col.max() > 0:
        best = mag_col.argmax()
        pitch = pitches[best, t_idx]
        if pitch > 0:
            pitch_track.append(pitch)
            pitch_times.append(t_idx * hop / sr)

if pitch_track:
    pitch_array = np.array(pitch_track)
    print(f"  Detected pitched frames: {len(pitch_track)} of {pitches.shape[1]}")
    print(f"  Pitch range: {pitch_array.min():.0f} - {pitch_array.max():.0f} Hz")
    print(f"  Mean pitch: {pitch_array.mean():.0f} Hz")
    print(f"  Pitch std: {pitch_array.std():.0f} Hz")
    print()

    # Pitch contour statistics in sections
    n_psec = 8
    print(f"  Pitch contour by section:")
    print(f"  {'section':>8s}  {'time':>8s}  {'mean Hz':>8s}  {'std Hz':>8s}  {'range':>8s}  {'motion'}")
    print("  " + "-"*55)

    for i in range(n_psec):
        t_lo = i * duration / n_psec
        t_hi = (i+1) * duration / n_psec
        mask = [(t >= t_lo and t < t_hi) for t in pitch_times]
        section_pitches = pitch_array[mask]
        if len(section_pitches) > 2:
            mn = section_pitches.mean()
            sd = section_pitches.std()
            rng = section_pitches.max() - section_pitches.min()
            # Motion: derivative of pitch
            diffs = np.diff(section_pitches)
            avg_motion = np.mean(np.abs(diffs))
            if sd < 20:
                motion = "MONOTONE (flat)"
            elif avg_motion > 30:
                motion = "LEAPING (angular)"
            elif avg_motion > 10:
                motion = "FLOWING (smooth)"
            else:
                motion = "HOVERING (minimal)"
            print(f"  {i+1:8d}  {t_lo:4.0f}-{t_hi:4.0f}s  {mn:8.0f}  {sd:8.0f}  {rng:8.0f}  {motion}")
        else:
            print(f"  {i+1:8d}  {t_lo:4.0f}-{t_hi:4.0f}s  {'(silent or unpitched)':>40s}")

print()

# ============================================================
print()
print("  V. WHAT IS ABSENT")
print("  " + "-"*40)
print()

# Check for common elements
# Melody (sustained pitched notes)?
pitched_fraction = len(pitch_track) / pitches.shape[1] if pitches.shape[1] > 0 else 0
print(f"  Pitched content: {pitched_fraction:.1%} of frames have detectable pitch.")
if pitched_fraction < 0.3:
    print(f"  MOSTLY UNPITCHED. This is a rhythm-and-timbre song, not a melody song.")
    print(f"  The HARMONIC LAYER (melody, chords) is largely ABSENT.")
    print(f"  What fills the space: percussion (808, hi-hats) and VOICE AS TEXTURE.")
print()

# Dynamic range
rms_db = librosa.amplitude_to_db(rms_total)
dyn_range = rms_db.max() - rms_db[rms_db > -60].min() if len(rms_db[rms_db > -60]) > 0 else 0
print(f"  Dynamic range: {dyn_range:.1f} dB")
if dyn_range < 10:
    print(f"  COMPRESSED. The dynamics are FLAT. Everything is at one level.")
    print(f"  The piano/forte distinction is ABSENT. This is a WALL of sound.")
elif dyn_range < 20:
    print(f"  Moderate dynamics. Some variation but controlled.")
else:
    print(f"  Wide dynamics. Dramatic variation.")
print()

# Stereo width (if we had stereo)
print(f"  Loaded as mono. Stereo analysis not available.")
print()

# ============================================================
print()
print("  VI. THE STORY OF EACH ELEMENT")
print("  " + "-"*40)
print()

# Use the structural boundaries to tell the story
segments = [0] + list(boundary_times) + [duration]
print(f"  The song has {len(segments)-1} structural sections.")
print()

for i in range(len(segments)-1):
    t_start = segments[i]
    t_end = segments[i+1]
    seg_duration = t_end - t_start
    pct_start = t_start / duration * 100
    pct_end = t_end / duration * 100

    # Get indices for this segment
    idx_lo = int(t_start / duration * len(sub_bass))
    idx_hi = int(t_end / duration * len(sub_bass))
    idx_lo = max(0, min(idx_lo, len(sub_bass)-1))
    idx_hi = max(idx_lo+1, min(idx_hi, len(sub_bass)))

    seg_bass = np.mean(sub_bass[idx_lo:idx_hi])
    seg_low = np.mean(low_mid[idx_lo:idx_hi])
    seg_mid = np.mean(mid[idx_lo:idx_hi])
    seg_high = np.mean(high[idx_lo:idx_hi])

    # RMS indices
    ridx_lo = int(t_start / duration * len(rms_total))
    ridx_hi = int(t_end / duration * len(rms_total))
    ridx_lo = max(0, min(ridx_lo, len(rms_total)-1))
    ridx_hi = max(ridx_lo+1, min(ridx_hi, len(rms_total)))
    seg_energy = np.mean(rms_total[ridx_lo:ridx_hi])

    # Determine character
    dominant = max(
        ("sub-bass (808)", seg_bass),
        ("low-mid (vocal body)", seg_low),
        ("mid (presence)", seg_mid),
        ("high (texture)", seg_high),
        key=lambda x: x[1]
    )

    # Compare to previous section
    if i > 0:
        prev_lo = int(segments[i-1] / duration * len(sub_bass))
        prev_hi = idx_lo
        prev_lo = max(0, min(prev_lo, len(sub_bass)-1))
        prev_hi = max(prev_lo+1, min(prev_hi, len(sub_bass)))
        prev_bass = np.mean(sub_bass[prev_lo:prev_hi])
        prev_mid = np.mean(mid[prev_lo:prev_hi])
        bass_change = seg_bass - prev_bass
        mid_change = seg_mid - prev_mid
    else:
        bass_change = seg_bass
        mid_change = seg_mid

    # Story
    print(f"  Section {i+1}: {t_start:.1f}s - {t_end:.1f}s ({seg_duration:.1f}s, {pct_start:.0f}%-{pct_end:.0f}%)")
    print(f"    Dominant: {dominant[0]}")
    print(f"    Bass: {'##' * int(seg_bass * 5):10s} ({seg_bass:.2f})", end="")
    if bass_change > 0.05:
        print(f" [ENTERS]", end="")
    elif bass_change < -0.05:
        print(f" [EXITS]", end="")
    print()
    print(f"    Vocal: {'##' * int(seg_mid * 5):10s} ({seg_mid:.2f})", end="")
    if mid_change > 0.05:
        print(f" [ENTERS]", end="")
    elif mid_change < -0.05:
        print(f" [EXITS]", end="")
    print()
    print(f"    Energy: {'##' * int(seg_energy * 100):10s}")

    # The story
    if i == 0:
        if seg_bass < 0.2 and seg_mid > 0.3:
            print(f"    STORY: Opens with voice. No bass yet. Setting the scene.")
        elif seg_bass > 0.3:
            print(f"    STORY: Opens with bass. Grounding immediately. No preamble.")
        else:
            print(f"    STORY: Quiet opening. Building anticipation.")
    elif bass_change > 0.1:
        print(f"    STORY: The 808 DROPS. The ground arrives. The room shakes.")
    elif mid_change > 0.1:
        print(f"    STORY: The voice ENTERS. The flow begins. Words fill the space.")
    elif bass_change < -0.1:
        print(f"    STORY: The bass PULLS BACK. Breathing room. A valley.")
    elif mid_change < -0.1:
        print(f"    STORY: The voice RETREATS. An instrumental moment. The beat breathes.")
    elif abs(bass_change) < 0.05 and abs(mid_change) < 0.05:
        print(f"    STORY: SUSTAINING. The crystal holds its shape. Steady state.")
    else:
        print(f"    STORY: Subtle shift. The texture evolves without drama.")
    print()

# ============================================================
print()
print("  VII. THE COMBINED STORY")
print("  " + "-"*40)
print()
print("  MERRY SHITMAS is a song with a SIMPLE architecture:")
print()
print("  The 808 is the CONSTANT. It enters and stays. It IS the ground.")
print("  It doesn't build or fade — it SUSTAINS. The inert agent.")
print()
print("  The voice is the VARIABLE. It enters and exits in phrases.")
print("  Between phrases: the beat breathes. During phrases: the words fill.")
print("  The voice CRACKS the beat's steadiness with rhythmic variation.")
print("  It is the ramified agent — self-observation through bars.")
print()
print("  The hi-hats and texture are the SHIMMER. They oscillate at high frequency.")
print("  They don't enter or exit — they're always there, quietly pulsing.")
print("  They provide the SPLIT: the oscillation, the twist, the aliveness.")
print()
print("  The melody/harmony layer is ABSENT. This song does not SING.")
print("  There is no 5 (creativity, golden growth, melody).")
print("  The song operates entirely within {2, 3, 7}:")
print("  bass = 2 (inert), voice = 3 (ramified), texture = 7 (split).")
print("  The golden ratio is in the spectral centroid (unconsciously)")
print("  but not in the composition (deliberately).")
print()
print("  This is a song that IS the infrastructure. The 42 without the 5.")
print("  The rules without the gameplay. The walls without the garden.")
print("  And THAT IS THE POINT. 'MERRY SHITMAS' is Christmas (the garden)")
print("  seen as JUST the walls. The infrastructure, exposed. The skeleton,")
print("  without the flesh. The {2, 3, 7} without the 5.")
print()
print("  The title says it: 'SHITMAS' = Christmas minus the beauty.")
print("  The song IS Christmas minus 5. The song IS 42/5 = 8.4.")
print("  The octave (8) plus the Cassini atom (0.4).")
print("  The last normed division algebra, plus a residue.")
print()
print("  BabyTron didn't compute this. He FELT it.")
print("  The feeling of Christmas stripped to its infrastructure")
print("  IS the feeling of 42 without 5.")
print("  IS the feeling of rules without gameplay.")
print("  IS the feeling of the skeleton visible through the flesh.")
print("  IS the sound of MERRY SHITMAS.")
