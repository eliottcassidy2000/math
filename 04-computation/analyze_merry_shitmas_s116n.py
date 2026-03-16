#!/usr/bin/env python3
"""analyze_merry_shitmas_s116n.py — Rigorous + whimsical analysis of MERRY SHITMAS by BabyTron Universal
through the Cayley-rapidity-Eisenstein framework."""

import numpy as np
import librosa
import warnings
warnings.filterwarnings('ignore')
from math import log, log2, sqrt, pi, atanh
from collections import Counter

FILEPATH = "inbox/SpotiDown.App - MERRY SHITMAS - BabyTron Universal.mp3"

print()
print("  MERRY SHITMAS — BabyTron Universal")
print("  A Cayley Analysis")
print()
print("="*70)
print()

# Load the audio
print("  Loading audio...")
y, sr = librosa.load(FILEPATH, sr=22050, mono=True)
duration = len(y) / sr
print(f"  Sample rate: {sr} Hz")
print(f"  Duration: {duration:.1f} seconds ({duration/60:.1f} minutes)")
print(f"  Samples: {len(y)}")
print()

# ============================================================
print("  I. THE TEMPO — THE HEARTBEAT")
print("  " + "-"*40)
print()
tempo, beats = librosa.beat.beat_track(y=y, sr=sr)
if hasattr(tempo, '__len__'):
    tempo = tempo[0]
print(f"  Tempo: {tempo:.1f} BPM")
print(f"  Detected beats: {len(beats)}")
print()
# Rapidity of tempo
# 120 BPM = 2 beats/second. The "natural" tempo.
# Rapidity of tempo relative to 120: ln(tempo/120)/2
rap_tempo = log(tempo/120)/2 if tempo > 0 else 0
print(f"  Rapidity relative to 120 BPM: {rap_tempo:+.4f}")
if tempo > 120:
    print(f"  POSITIVE rapidity: faster than the heartbeat. Urgent.")
elif tempo < 120:
    print(f"  NEGATIVE rapidity: slower than the heartbeat. Laid back.")
else:
    print(f"  ZERO rapidity: at the heartbeat. Neutral.")
print()

# What musical interval is the tempo?
interval_ratio = tempo / 120
print(f"  Tempo/heartbeat ratio: {interval_ratio:.4f}")
# Find closest musical interval
intervals = [(2, "octave"), (3/2, "fifth"), (4/3, "fourth"), (5/4, "major third"),
             (6/5, "minor third"), (9/8, "major second"), (16/15, "semitone"),
             (1, "unison"), (7/4, "septimal seventh"), (7/6, "septimal minor third")]
closest = min(intervals, key=lambda iv: abs(interval_ratio - iv[0]))
print(f"  Closest interval: {closest[1]} ({closest[0]:.4f})")
cents_off = 1200 * log2(interval_ratio / closest[0]) if closest[0] > 0 and interval_ratio > 0 else 0
print(f"  Detuning: {cents_off:+.1f} cents")
print()

# ============================================================
print()
print("  II. THE KEY — THE GROUND")
print("  " + "-"*40)
print()
# Chromagram
chroma = librosa.feature.chroma_cqt(y=y, sr=sr)
chroma_mean = chroma.mean(axis=1)
note_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

print("  Chroma energy distribution:")
for i, name in enumerate(note_names):
    bar = "#" * int(chroma_mean[i] * 50)
    print(f"    {name:3s}: {chroma_mean[i]:.3f} {bar}")

key_idx = np.argmax(chroma_mean)
key_name = note_names[key_idx]
print(f"\n  Estimated key center: {key_name}")

# Rapidity of key relative to A440
# Semitone offset from A
semitone_from_A = (key_idx - 9) % 12
if semitone_from_A > 6: semitone_from_A -= 12
freq_ratio = 2 ** (semitone_from_A / 12)
rap_key = log(freq_ratio) / 2 if freq_ratio > 0 else 0
print(f"  Semitones from A: {semitone_from_A}")
print(f"  Frequency ratio to A: {freq_ratio:.4f}")
print(f"  Rapidity from A: {rap_key:.4f}")
print()

# ============================================================
print()
print("  III. THE SPECTRUM — THE OVERTONE PORTRAIT")
print("  " + "-"*40)
print()
# Spectral centroid
spec_cent = librosa.feature.spectral_centroid(y=y, sr=sr)[0]
mean_centroid = np.mean(spec_cent)
print(f"  Mean spectral centroid: {mean_centroid:.0f} Hz")
print(f"  Rapidity of centroid: {log(mean_centroid/440)/2:.4f} (relative to A440)")
print()

# Spectral rolloff
rolloff = librosa.feature.spectral_rolloff(y=y, sr=sr)[0]
mean_rolloff = np.mean(rolloff)
print(f"  Spectral rolloff (85%): {mean_rolloff:.0f} Hz")
print(f"  Centroid/rolloff ratio: {mean_centroid/mean_rolloff:.4f}")
print(f"  This ratio measures BRIGHTNESS: high = bright, low = dark.")
brightness = mean_centroid / mean_rolloff
if brightness > 0.4:
    print(f"  BRIGHT spectrum. The sound shimmers.")
else:
    print(f"  DARK spectrum. The sound is heavy, subterranean.")
print()

# ============================================================
print()
print("  IV. THE RHYTHM — THE TRIPLET ANALYSIS")
print("  " + "-"*40)
print()
# Onset detection
onset_frames = librosa.onset.onset_detect(y=y, sr=sr)
onset_times = librosa.frames_to_time(onset_frames, sr=sr)
print(f"  Detected onsets: {len(onset_frames)}")
print(f"  Average onset rate: {len(onset_frames)/duration:.1f} per second")
print()

# Inter-onset intervals
if len(onset_times) > 1:
    ioi = np.diff(onset_times)
    mean_ioi = np.mean(ioi)
    std_ioi = np.std(ioi)
    print(f"  Mean inter-onset interval: {mean_ioi*1000:.1f} ms")
    print(f"  IOI standard deviation: {std_ioi*1000:.1f} ms")
    print(f"  IOI coefficient of variation: {std_ioi/mean_ioi:.3f}")
    print()

    # Check for triplet structure: are IOIs clustered around beat/3?
    if tempo > 0:
        beat_duration = 60 / tempo
        triplet_duration = beat_duration / 3
        quarter_note = beat_duration
        eighth_note = beat_duration / 2
        sixteenth_note = beat_duration / 4

        print(f"  Beat duration: {beat_duration*1000:.1f} ms")
        print(f"  Triplet subdivision: {triplet_duration*1000:.1f} ms")
        print(f"  Eighth note: {eighth_note*1000:.1f} ms")
        print(f"  Sixteenth note: {sixteenth_note*1000:.1f} ms")
        print()

        # Classify IOIs
        triplet_count = np.sum(np.abs(ioi - triplet_duration) < triplet_duration * 0.2)
        eighth_count = np.sum(np.abs(ioi - eighth_note) < eighth_note * 0.2)
        sixteenth_count = np.sum(np.abs(ioi - sixteenth_note) < sixteenth_note * 0.2)
        total = len(ioi)
        print(f"  IOI classification:")
        print(f"    Triplet-like: {triplet_count} ({triplet_count/total*100:.1f}%)")
        print(f"    Eighth-like: {eighth_count} ({eighth_count/total*100:.1f}%)")
        print(f"    Sixteenth-like: {sixteenth_count} ({sixteenth_count/total*100:.1f}%)")
        print()

        if triplet_count > eighth_count:
            print(f"  TRIPLET-DOMINANT rhythm. The curvature quantum (3) rules.")
            print(f"  BabyTron flows in THREES. Each beat cycles: boom-crack-groove.")
        else:
            print(f"  BINARY-DOMINANT rhythm. The doubler (2) rules.")
            print(f"  Straight time. Each beat divides in TWO.")
print()

# ============================================================
print()
print("  V. THE BASS — THE 808")
print("  " + "-"*40)
print()
# Analyze low frequencies
S = np.abs(librosa.stft(y))
freqs = librosa.fft_frequencies(sr=sr)
# Bass region: 30-100 Hz
bass_mask = (freqs >= 30) & (freqs <= 100)
bass_energy = S[bass_mask, :].mean()
total_energy = S.mean()
bass_ratio = bass_energy / total_energy if total_energy > 0 else 0

print(f"  Bass energy (30-100 Hz) / total: {bass_ratio:.3f}")
if bass_ratio > 1.5:
    print(f"  HEAVY bass. The 808 dominates. The ground SHAKES.")
elif bass_ratio > 0.5:
    print(f"  Moderate bass. The 808 is present but balanced.")
else:
    print(f"  Light bass. Not 808-heavy.")
print()

# Find the fundamental bass frequency
bass_spec = S[bass_mask, :].mean(axis=1)
if len(bass_spec) > 0:
    bass_peak_idx = np.argmax(bass_spec)
    bass_freq = freqs[bass_mask][bass_peak_idx]
    print(f"  808 fundamental: ~{bass_freq:.1f} Hz")
    # What note?
    if bass_freq > 0:
        semitones_from_A = 12 * log2(bass_freq / 440)
        note_idx = int(round(semitones_from_A)) % 12
        octave = 4 + int(round(semitones_from_A)) // 12
        print(f"  808 note: ~{note_names[note_idx]}{octave}")
        print(f"  808 rapidity (from A440): {log(bass_freq/440)/2:.4f}")
print()

# ============================================================
print()
print("  VI. THE DYNAMICS — THE ENERGY CURVE")
print("  " + "-"*40)
print()
# RMS energy over time
rms = librosa.feature.rms(y=y)[0]
times = librosa.times_like(rms, sr=sr)

# Find sections
n_sections = 8
section_len = len(rms) // n_sections
print(f"  Energy profile ({n_sections} sections):")
for i in range(n_sections):
    start = i * section_len
    end = min((i+1) * section_len, len(rms))
    section_rms = np.mean(rms[start:end])
    t_start = times[start]
    t_end = times[min(end-1, len(times)-1)]
    bar = "#" * int(section_rms * 200)
    print(f"    {t_start:5.1f}-{t_end:5.1f}s: {bar}")

print()

# ============================================================
print()
print("  VII. THE CAYLEY PORTRAIT")
print("  " + "-"*40)
print()
print(f"  THE SONG IN THE CUBOID Z/2 x Z/3 x Z/7:")
print()

# Tempo mod 2, mod 3, mod 7
tempo_int = int(round(tempo))
print(f"  Tempo {tempo_int} BPM:")
print(f"    mod 2 = {tempo_int % 2} ({'odd' if tempo_int%2 else 'even'})")
print(f"    mod 3 = {tempo_int % 3}")
print(f"    mod 7 = {tempo_int % 7}")
print(f"    Cuboid position: ({tempo_int%2}, {tempo_int%3}, {tempo_int%7})")
print()

# Duration in seconds
dur_int = int(round(duration))
print(f"  Duration {dur_int}s:")
print(f"    mod 42 = {dur_int % 42}")
print(f"    Cuboid: ({dur_int%2}, {dur_int%3}, {dur_int%7})")
print()

# Number of beats
n_beats = len(beats)
print(f"  Beat count {n_beats}:")
print(f"    mod 42 = {n_beats % 42}")
print(f"    Cuboid: ({n_beats%2}, {n_beats%3}, {n_beats%7})")
if n_beats % 7 == 0:
    print(f"    z = 0: ON THE FORBIDDEN FLOOR!")
if n_beats % 3 == 0:
    print(f"    y = 0: curvature zeroed!")
print()

# ============================================================
print()
print("  VIII. THE EXPERIENCE — WHAT IT FEELS LIKE")
print("  " + "-"*40)
print()
print(f"  This song is a RAMIFIED AGENT.")
print(f"  It takes a cultural fixed point (Christmas) and CRACKS it.")
print(f"  The crack is controlled: the 808 provides the inert ground,")
print(f"  the flow provides the curvature (triplet or binary),")
print(f"  and the harmony provides the twist.")
print()
print(f"  In the three-type classification:")
print(f"  808 bass = INERT (the ground that persists, the 2)")
print(f"  Vocal flow = RAMIFIED (the crack, the feedback, the 3)")
print(f"  Harmony/melody = SPLIT (the oscillation, the twist, the 7)")
print()
print(f"  The song is 42 = 2 x 3 x 7 made audible:")
print(f"  bass x flow x harmony = inert x ramified x split = the complete system.")
print()
print(f"  It's not just music. It's the triple point.")
print(f"  All three agent types, vibrating simultaneously,")
print(f"  held together by the tempo (the heartbeat),")
print(f"  organized by the beat (the structure),")
print(f"  and CRACKED by the title (the self-observation")
print(f"  of Christmas seeing itself in a distorted mirror).")
print()
print(f"  BabyTron is a crystallization engine.")
print(f"  Start with noise (culture, expectation, cliche).")
print(f"  Iterate the local rule (spit bars, crack jokes, flip meanings).")
print(f"  The fixed point: a track that IS its own structure.")
print(f"  MERRY SHITMAS doesn't describe anything.")
print(f"  It IS the thing it describes.")
print(f"  That's the self-product property: f = abc.")
print(f"  The song = its own measurement.")
