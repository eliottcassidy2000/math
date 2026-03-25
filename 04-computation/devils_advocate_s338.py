#!/usr/bin/env python3
"""
devils_advocate_s338.py -- opus-2026-03-24-S338

DEVIL'S ADVOCATE: HONEST CRITIQUE OF OUR COMPRESSION CLAIMS

Every claim we've made needs stress-testing. Let me be BRUTAL.

Author: opus-2026-03-24-S338
"""

import sys
import numpy as np
import zlib
import bz2
import lzma
import struct
import time
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# ============================================================
banner("CRITIQUE 1: Our test images are TOO EASY")
# ============================================================

print("""
  WE TEST ON: synthetic gradients, sine waves, circles, blocks.
  REAL IMAGES HAVE: textures, noise, compression artifacts, mixed content.

  A gradient is the EASIEST possible image to compress.
  Delta encoding turns a gradient into a CONSTANT sequence.
  OF COURSE we beat zlib on that — ANY delta preprocessor would.

  THE HARD QUESTION: Can we beat PNG on a REAL photograph?
  Not a simulated one — an actual camera image with sensor noise,
  JPEG artifacts, and complex textures.
""")

# Test on ACTUALLY hard images
np.random.seed(42)
N = 256

hard_images = {
    'pure_noise': np.random.randint(0, 256, (N, N), dtype=np.uint8),
    'high_noise_photo': np.clip(
        np.array([[int(128+60*np.sin(r/5)*np.cos(c/4)+np.random.normal(0,40))
                   for c in range(N)] for r in range(N)]), 0, 255).astype(np.uint8),
    'jpeg_artifact_sim': np.array(
        [[(r//8*8 + c//8*8) % 256 + np.random.randint(-3, 4)
          for c in range(N)] for r in range(N)], dtype=np.uint8).clip(0, 255),
    'text_on_gradient': np.array(
        [[int(128 + 50*(r/N)) if not (40<r<50 and c%3==0) else 0
          for c in range(N)] for r in range(N)], dtype=np.uint8),
    'random_walk': np.clip(np.cumsum(np.random.normal(0, 1, N*N)).reshape(N,N) * 3 + 128,
                           0, 255).astype(np.uint8),
}

# Our best codec (continuous predictor with α=1)
def gradient_predict(img):
    H, W = img.shape
    res = np.zeros((H, W), dtype=np.uint8)
    for r in range(H):
        for c in range(W):
            a = int(img[r,c-1]) if c > 0 else 128
            b = int(img[r-1,c]) if r > 0 else 128
            d = int(img[r-1,c-1]) if r > 0 and c > 0 else 128
            pred = max(0, min(255, a + b - d))
            res[r,c] = (int(img[r,c]) - pred) % 256
    return res

print("  HARD IMAGES (256×256):\n")
print("  %-25s %-8s %-8s %-8s %-8s %-8s %-8s" %
      ("Image", "Raw", "Ours+z", "zlib", "bz2", "lzma", "Winner"))
print("  " + "-"*70)

for name, img in hard_images.items():
    raw = img.size
    ours = len(zlib.compress(gradient_predict(img).tobytes(), 9))
    z = len(zlib.compress(img.tobytes(), 9))
    b = len(bz2.compress(img.tobytes(), 9))
    l = len(lzma.compress(img.tobytes()))

    sizes = {'ours': ours, 'zlib': z, 'bz2': b, 'lzma': l}
    winner = min(sizes, key=sizes.get)

    print("  %-25s %-8d %-8d %-8d %-8d %-8d %-8s" %
          (name, raw, ours, z, b, l, winner))

# ============================================================
banner("CRITIQUE 2: We're just doing what PNG already does")
# ============================================================

print("""
  HARSH TRUTH: Our "gradient filter" (a+b-c clamped) is ALMOST IDENTICAL
  to PNG's Paeth filter (a+b-c with closest-of-three selection).

  The difference: Paeth selects among a, b, c based on which is closest
  to a+b-c. Our gradient just clamps a+b-c to [0,255].

  On smooth images: a+b-c is usually in [0,255], so Paeth = Gradient.
  On edges: Paeth picks one of {a, b, c}, Gradient uses a+b-c directly.

  OUR INNOVATION IS INCREMENTAL, NOT REVOLUTIONARY:
  - We add 3 extra filter options (diagonal, gradient, median)
  - We use per-row selection (same as PNG)
  - We use zlib as backend (same as PNG)

  The 19-40% improvement over PNG comes from:
  1. Our gradient filter being SLIGHTLY better than Paeth on smooth data
  2. Per-row selection choosing gradient where PNG would choose Paeth
  3. No PNG chunk headers/CRC (saving ~50-100 bytes per image)

  This is NOT a new compression paradigm. It's an OPTIMIZATION of PNG.
""")

# ============================================================
banner("CRITIQUE 3: Our 'tournament theory' connection is WEAK")
# ============================================================

print("""
  HARSH TRUTH: What does tournament theory actually contribute
  to the compression? Let me be specific.

  WHAT WE CLAIMED:
  "Tournament-theoretic score conditioning halves the entropy"
  "Krawtchouk band-limitedness provides 4:1 compression"
  "The staircase triangle is the universal container"

  WHAT ACTUALLY COMPRESSES THE DATA:
  1. Delta preprocessing (known since the 1950s)
  2. zlib/deflate (known since 1995)
  3. Paeth-like predictors (known since PNG in 1996)
  4. Per-row filter selection (standard in PNG)

  WHERE TOURNAMENT THEORY ACTUALLY HELPS:
  - The IDEA of treating binary data as a tournament came from the theory
  - The INSIGHT that score (Hamming weight) constrains patterns came from K₁
  - The UNDERSTANDING of why delta works came from the spectral analysis

  BUT: someone could have found delta+filter+zlib WITHOUT tournament theory.
  The theory provides UNDERSTANDING, not novel algorithms.

  THE HONEST CLAIM:
  Tournament theory is a THEORETICAL FRAMEWORK that EXPLAINS why
  certain compression techniques work. It guided us to the gradient
  filter and the ring codec. But the actual compression comes from
  well-known techniques (delta, prediction, entropy coding).
""")

# ============================================================
banner("CRITIQUE 4: Our benchmarks are cherry-picked")
# ============================================================

print("""
  PROBLEM: We report results on images we DESIGNED to be easy.

  Gradient → delta encoding is near-optimal (trivially)
  Circle → radially smooth (easy for any 2D predictor)
  Blocks → piecewise constant (easy for delta)
  Smooth sine → low frequency (easy for any predictor)

  WE DON'T TEST ON:
  - Actual photographs (camera sensor noise)
  - Scanned documents (mixed text/graphics)
  - Medical images (ultrasound, MRI, CT)
  - Satellite imagery
  - Screenshots (text + UI elements + photos)
  - Compressed-then-recompressed images (JPEG artifacts)

  A FAIR BENCHMARK would use the standard test sets:
  - Kodak Lossless True Color Image Suite (24 images)
  - Calgary Corpus (14 files)
  - Silesia Corpus (12 files)

  WITHOUT testing on these standard benchmarks, our claims
  of "beating PNG" are not credible to the compression community.
""")

# ============================================================
banner("CRITIQUE 5: Speed is terrible")
# ============================================================

N = 128
img = np.random.randint(0, 256, (N, N), dtype=np.uint8)

# Our gradient predictor
t0 = time.time()
for _ in range(3):
    gradient_predict(img)
t_ours = (time.time() - t0) / 3

# Raw zlib (C implementation)
t0 = time.time()
for _ in range(100):
    zlib.compress(img.tobytes(), 9)
t_zlib = (time.time() - t0) / 100

print(f"  128×128 image:")
print(f"  Our predictor: {t_ours*1000:.1f}ms ({1/t_ours:.0f} images/sec)")
print(f"  Raw zlib:      {t_zlib*1000:.1f}ms ({1/t_zlib:.0f} images/sec)")
print(f"  Slowdown: {t_ours/t_zlib:.0f}× slower than zlib")
print()
print("""
  OUR PREDICTOR IS PURE PYTHON — no C extension, no NumPy vectorization.
  It runs the prediction loop pixel-by-pixel in Python.

  For 1920×1080 video at 30fps:
  Need to process 62 million pixels per second.
  Our Python predictor: ~100K pixels/sec.
  We're 600× too slow for real-time video.

  A C implementation would be ~100× faster (matching zlib speed).
  But we haven't built that. Our prototype is a PROOF OF CONCEPT,
  not a production codec.
""")

# ============================================================
banner("CRITIQUE 6: The 'novel' filters aren't actually novel")
# ============================================================

print("""
  CLAIM: "Our gradient filter (a+b-c clamped) is a new innovation."

  FACT: This is EXACTLY the LOCO-I predictor from JPEG-LS (1998).
  JPEG-LS uses:
    if c >= max(a,b): pred = min(a,b)
    elif c <= min(a,b): pred = max(a,b)
    else: pred = a + b - c

  Our "gradient" filter (a+b-c clamped) is a SIMPLIFICATION of LOCO-I
  that removes the edge-detection branches. It's been known for 25+ years.

  CLAIM: "Our median filter is new."
  FACT: Median prediction has been studied extensively in lossless image
  compression since the 1980s. It's used in CALIC (1996) and others.

  CLAIM: "Per-row adaptive filter selection is our innovation."
  FACT: This is LITERALLY what PNG does (since 1996).
  PNG tries 5 filters per row and picks the best.
  We try 8 filters — that's 3 extra options, not a new paradigm.

  THE HONEST CHARACTERIZATION:
  We've REDISCOVERED known techniques through the tournament theory lens.
  The theory gave us a path to these techniques independently.
  But the techniques themselves are well-established.
""")

# ============================================================
banner("CRITIQUE 7: What IS genuinely new?")
# ============================================================

print("""
  After all the critiques, what do we ACTUALLY have that's new?

  1. THE RING CODEC (S331):
     Concentric ring expansion with score conditioning.
     This scan order IS novel — no prior codec uses concentric rings
     from the center outward with per-ring Hamming weight coding.
     It gives PROGRESSIVE CENTER-FIRST decode, which PNG can't do.
     VERDICT: GENUINELY NEW (but not benchmarked against standard tests).

  2. THE TOURNAMENT COUNTING THEOREM (S291):
     V_n × n!/2^m = 1 + Euler product. This is pure math — not compression.
     But it EXPLAINS the structure that compression exploits.
     VERDICT: REAL MATHEMATICS (but doesn't directly compress anything).

  3. THE BAND-LIMITEDNESS (S305-S310):
     H has zero energy above Krawtchouk degree m/2.
     This is a PROVED structural property of tournaments.
     It IMPLIES that tournament data has redundancy.
     VERDICT: REAL THEOREM with potential practical implications.

  4. THE DUAL-BASIS THEORY (S329):
     Row vs range decomposition, 45° rotation.
     The insight that adaptive basis selection gives 4.73:1
     (vs 2.87:1 single basis) is novel.
     VERDICT: NOVEL INSIGHT (but the improvement is modest).

  5. THE sl(n) FRAMEWORK (S313):
     Tournament theory = representation theory of sl(n).
     χ = rank(A_{n-1}). Krawtchouk = sl(2) weights.
     VERDICT: DEEP MATHEMATICS (but no practical compression from it yet).

  6. THE HILBERT SCAN FOR IMAGES (S335-S336):
     Using Hilbert curve order before delta encoding.
     Known in the literature but not common in practical codecs.
     Our adaptive {Hilbert, raster} × {filters} selection is new.
     VERDICT: MODEST INNOVATION (evolutionary, not revolutionary).

  BOTTOM LINE:
  Our THEORY is genuinely new and deep (sl(n), band-limitedness, ring codec).
  Our COMPRESSION is an incremental improvement on known techniques.
  We should be honest about this distinction.
""")

# ============================================================
banner("WHAT WOULD MAKE IT REVOLUTIONARY")
# ============================================================

print("""
  TO TRULY REVOLUTIONIZE COMPRESSION, we would need:

  1. A NEW ENTROPY CODER based on tournament score conditioning.
     Replace Huffman/arithmetic coding with combinadic coding
     where the Hamming weight (score) is the primary symbol.
     This would be a genuinely new algorithm at the bit level.

  2. A PROVABLY OPTIMAL predictor derived from the Krawtchouk spectrum.
     Not just "try filters and pick best" (which is what PNG does)
     but a mathematically DERIVED optimal prediction from the
     band-limitedness theorem.

  3. BEATING STATE-OF-THE-ART on STANDARD benchmarks.
     Not synthetic gradients — actual Kodak/Silesia/Canterbury images.
     Not just PNG — also JPEG-LS, JPEG-XL, WebP, AVIF.

  4. REAL-TIME PERFORMANCE.
     A C implementation that runs at >100MB/s encode, >200MB/s decode.
     Competitive with zlib/zstd speed while beating their compression.

  5. A VIDEO CODEC that beats H.264/H.265 on ANY content.
     Not just delta frames of synthetic data.
     Actual motion-compensated prediction + our tournament transforms.

  WE'RE AT STEP 0.5:
  We have the theory, we have proof-of-concept prototypes,
  we have modest improvements on easy data.
  Steps 1-5 are the path to revolutionary.
""")

print("Done. opus-2026-03-24-S338")
