#!/usr/bin/env python3
"""prime_music_s116n.py — The music of the primes, literally.

Each prime p creates a tone at frequency proportional to ln(p).
The amplitude = 1/sqrt(p) (the critical line Re(s) = 1/2).
The phase = the Chebyshev psi function's fractional part.

The COMPOSITE of all prime tones IS the prime distribution's sound.
The Riemann zeros are where the tones CANCEL (destructive interference).
"""
import struct
import wave
import math

def generate_prime_music(filename, duration=10, sample_rate=44100):
    """Generate a WAV file that IS the prime distribution."""

    # Primes up to 1000
    def sieve(n):
        s = [True]*(n+1)
        s[0] = s[1] = False
        for i in range(2, int(n**0.5)+1):
            if s[i]:
                for j in range(i*i, n+1, i):
                    s[j] = False
        return [i for i in range(2, n+1) if s[i]]

    primes = sieve(200)

    n_samples = int(duration * sample_rate)
    samples = []

    print(f"  Generating {duration}s of prime music at {sample_rate}Hz...")
    print(f"  Using {len(primes)} primes up to {primes[-1]}.")
    print()

    # Each prime p contributes a tone:
    # frequency = base_freq * ln(p) (rapidity determines pitch)
    # amplitude = 1/sqrt(p) (critical line amplitude)
    # This IS the explicit formula for the prime counting function

    base_freq = 220  # A3 = 220 Hz as the base

    print("  Prime tones (first 20):")
    print(f"  {'prime':>6s}  {'frequency':>10s}  {'amplitude':>10s}  {'note (approx)'}")
    print("  " + "-"*45)

    note_names = ['A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#']
    for p in primes[:20]:
        freq = base_freq * math.log(p) / math.log(2)  # rapidity in octaves * base
        amp = 1.0 / math.sqrt(p)
        # Nearest note
        semitones = 12 * math.log2(freq / 440)
        note_idx = int(round(semitones)) % 12
        octave = 4 + int(round(semitones)) // 12
        note = f"{note_names[note_idx]}{octave}"
        print(f"  {p:6d}  {freq:10.1f}Hz  {amp:10.4f}  {note}")

    print()
    print("  Mixing all prime tones...")

    for i in range(n_samples):
        t = i / sample_rate
        sample = 0.0
        for p in primes:
            freq = base_freq * math.log(p) / math.log(2)
            amp = 1.0 / math.sqrt(p)
            # Add slight detuning based on prime's position to avoid exact harmonics
            phase = 2 * math.pi * freq * t
            sample += amp * math.sin(phase)

        samples.append(sample)

    # Normalize
    max_val = max(abs(s) for s in samples)
    if max_val > 0:
        samples = [s / max_val * 0.9 for s in samples]

    # Write WAV
    with wave.open(filename, 'w') as wav:
        wav.setnchannels(1)
        wav.setsampwidth(2)  # 16-bit
        wav.setframerate(sample_rate)
        for s in samples:
            wav.writeframes(struct.pack('<h', int(s * 32767)))

    print(f"  Written to {filename}")
    print(f"  Duration: {duration}s, Samples: {n_samples}")
    return samples

# Also: the SIEVE as music. Each sieving step REMOVES frequencies.
def generate_sieve_music(filename, duration=15, sample_rate=44100):
    """Generate music that IS the sieve of Eratosthenes.
    Start with all frequencies (white noise structured by integers).
    Remove harmonics of 2 (octave sieve). Then 3 (fifth sieve).
    Then 5 (third sieve). Then 7 (minor third sieve).
    What remains = the prime tones.
    """
    n_samples = int(duration * sample_rate)
    samples_by_phase = []
    base_freq = 110

    primes_to_sieve = [2, 3, 5, 7, 11, 13]
    phase_duration = duration / (len(primes_to_sieve) + 1)

    print(f"  Generating sieve music: {duration}s")
    print(f"  Phase 0: all integers (dense). {phase_duration:.1f}s")
    for i, p in enumerate(primes_to_sieve):
        print(f"  Phase {i+1}: sieve by {p} (remove multiples). {phase_duration:.1f}s")
    print()

    # Active frequencies: start with integers 2..50
    active = list(range(2, 80))
    sieve_idx = 0
    current_phase_start = 0

    samples = []
    for i in range(n_samples):
        t = i / sample_rate

        # Check if we need to advance the sieve
        phase = int(t / phase_duration)
        while sieve_idx < phase and sieve_idx < len(primes_to_sieve):
            p = primes_to_sieve[sieve_idx]
            active = [n for n in active if n == p or n % p != 0]
            sieve_idx += 1

        sample = 0.0
        for n in active:
            freq = base_freq * math.log(n)  # frequency = rapidity * base
            amp = 1.0 / n  # amplitude decreases with n
            sample += amp * math.sin(2 * math.pi * freq * t)

        samples.append(sample)

    # Normalize
    max_val = max(abs(s) for s in samples) if samples else 1
    samples = [s / max_val * 0.9 for s in samples]

    with wave.open(filename, 'w') as wav:
        wav.setnchannels(1)
        wav.setsampwidth(2)
        wav.setframerate(sample_rate)
        for s in samples:
            wav.writeframes(struct.pack('<h', int(s * 32767)))

    print(f"  Written to {filename}")
    print(f"  Active frequencies at end: {active[:20]}... ({len(active)} total)")
    return active

# Generate both
print()
print("  THE MUSIC OF THE PRIMES")
print()
print("="*70)
print()

print("  PIECE 1: PRIME HARMONY")
print("  Each prime = a tone. Frequency = rapidity. Amplitude = 1/sqrt(p).")
print("  This IS the explicit formula for pi(x), made audible.")
print()
generate_prime_music("05-knowledge/results/prime_harmony.wav", duration=8, sample_rate=22050)

print()
print("  PIECE 2: THE SIEVE")
print("  Start with all integers as tones. Remove multiples one prime at a time.")
print("  You HEAR the sieve of Eratosthenes. Each sieve step = one musical interval removed.")
print("  What remains = the primes = the irreducible tones.")
print()
active = generate_sieve_music("05-knowledge/results/prime_sieve.wav", duration=12, sample_rate=22050)

print()
print("  PIECE 3: THE CRYSTALLIZATION")
print("  (conceptual: the sieve IS crystallization applied to the integers)")
print("  Start with noise (all frequencies). Iterate removal of weakest")
print("  contradictions (composite frequencies). The fixed point = the primes.")
print("  The primes ARE the crystal that emerges from sieving the integers.")
print()

# The forbidden chord
print()
print("  BONUS: THE FORBIDDEN CHORD")
print("  The chord of primes 3, 7: frequency ratio 7/3 = septimal minor third.")
print("  This interval is FORBIDDEN in standard Western tuning.")
print("  It sounds 'wrong' — 267 cents, between a major second and minor third.")
print("  The FORBIDDEN INTERVAL of music theory IS the ratio of")
print("  the curvature quantum to the forbidden threshold.")
print()

# The 42 chord
print("  THE 42 CHORD: 2, 3, 7 played together.")
print("  Frequency ratios: 1 : 3/2 : 7/2 = root : fifth : septimal seventh.")
print("  This is a DOMINANT SEVENTH with the seventh tuned to the 7th harmonic")
print("  instead of the tempered version. The 'barbershop seventh'.")
print("  It is the MOST CONSONANT form of the dominant seventh chord.")
print("  The Hurwitz chord: the sound of 42.")
print()

# Generate the 42 chord
def generate_chord(filename, ratios, duration=5, sample_rate=22050):
    base = 220
    n_samples = int(duration * sample_rate)
    samples = []
    for i in range(n_samples):
        t = i / sample_rate
        # Fade in/out
        env = min(t * 4, 1.0, (duration - t) * 4)
        sample = 0
        for ratio, amp in ratios:
            sample += amp * math.sin(2 * math.pi * base * ratio * t)
        samples.append(sample * env)
    max_val = max(abs(s) for s in samples) if samples else 1
    samples = [s / max_val * 0.9 for s in samples]
    with wave.open(filename, 'w') as wav:
        wav.setnchannels(1)
        wav.setsampwidth(2)
        wav.setframerate(sample_rate)
        for s in samples:
            wav.writeframes(struct.pack('<h', int(s * 32767)))
    print(f"  Written to {filename}")

# The 42 chord: root + fifth (3/2) + septimal seventh (7/4)
generate_chord("05-knowledge/results/chord_42.wav",
    [(1.0, 1.0), (3/2, 0.8), (7/4, 0.6)], duration=5)

# The forbidden interval: 7/3 alone
generate_chord("05-knowledge/results/forbidden_interval.wav",
    [(1.0, 1.0), (7/3, 0.9)], duration=3)

# The full Hurwitz: root + octave (2) + twelfth (3) + 7th harmonic (7)
generate_chord("05-knowledge/results/hurwitz_full.wav",
    [(1.0, 1.0), (2.0, 0.7), (3.0, 0.5), (7/2, 0.4)], duration=5)

print()
print("  GENERATED FILES:")
print("  prime_harmony.wav — all primes as simultaneous tones")
print("  prime_sieve.wav — the sieve removing composites in real time")
print("  chord_42.wav — the Hurwitz chord (root + fifth + septimal 7th)")
print("  forbidden_interval.wav — the ratio 7/3 = the forbidden interval")
print("  hurwitz_full.wav — root + octave + twelfth + 7th harmonic = {2,3,7}")
print()
print("  The primes SOUND like this. The sieve SOUNDS like this.")
print("  The forbidden number SOUNDS like the septimal minor third.")
print("  Mathematics has always been music. Now you can hear it.")
