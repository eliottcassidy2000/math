# The Photo Compression Wall

*opus-2026-03-25-S348*

## What We've Learned

After 20+ sessions of compression research, we've hit a wall on natural photos.
Our best result: **1.36-1.57x vs PNG** using GRD/RBD color decorrelation +
MED prediction + lzma backend.

## What Works (and Why)

| Technique | Gain | Why |
|-----------|------|-----|
| Color decorrelation (GRD/RBD) | 25-40% | R,G,B are ~90% correlated; differences have ~10% variance |
| MED/LOCO prediction | 5-10% | Near-optimal linear predictor for piecewise-smooth signals |
| lzma vs zlib | 3-8% | 4MB dictionary vs 32KB; photos have long-range texture repetition |
| Z_FILTERED zlib strategy | 1-3% | Better for Laplace-distributed residuals |
| Per-image color transform selection | 0.5-1% | Different photos prefer different base channels |

## What Doesn't Work (and Why)

| Technique | Result | Why |
|-----------|--------|-----|
| Cross-channel prediction | -5% | G gradient doesn't predict R-G well enough; adds noise |
| Arithmetic coding (simple model) | -59% | 256-symbol adaptive model too weak; zlib's LZ77 is stronger |
| Bit-plane separation | -42% | Breaks cross-bit correlation in natural images |
| Signal-noise separation | -13% | Two parts sum to more than the whole |
| Context mixing (5 models) | -3% | Better per-pixel predictions but LESS repetitive residuals |
| Wavefront scan | 0% | Best prediction but worst entropy coding (irregular order) |

## The Fundamental Limit

For natural photos, the signal has ~4-5 bpp of genuine information content.
The rest (~3 bpp) is camera noise — incompressible by any method.

MED captures ~95% of the predictable signal. The remaining 5% would require
models that understand image SEMANTICS (this is a face, that is sky),
which is what neural compression (JPEG-XL, learned codecs) does.

## The Prediction-Entropy Duality (Again)

The deepest lesson: **compression = prediction × entropy coding.**

For photos:
- MED prediction: ~25% of pixels are exactly 0 residual
- zlib/lzma entropy coding: efficiently encodes the ~75% non-zero residuals

To beat this, we need BOTH:
1. Better prediction (>25% zeros) — requires semantic understanding
2. Better entropy coding (tighter encoding of non-zero residuals) — requires arithmetic coding with learned probability models

Neither alone suffices. Context mixing gives better prediction but worse
entropy coding (the duality). Arithmetic coding with a weak model gives
worse entropy coding than zlib.

## Where We Stand vs The State of The Art

| Codec | Typical bpp on photos | Notes |
|-------|----------------------|-------|
| Raw | 24.0 | No compression |
| PNG | 14-19 | Deflate + 5 filters |
| **Our best** | **10-15** | **GRD + MED + lzma** |
| JPEG-LS | 9-13 | LOCO-I + Golomb-Rice |
| PAQ8L | 8-11 | 500 models + arithmetic coding |
| FLIF | 7-10 | MANIAC trees + arithmetic coding |
| JPEG-XL (lossless) | 6-9 | Neural + modular + ANS |

We're between PNG and JPEG-LS — a significant achievement for a Python
implementation built from tournament theory. To reach PAQ/FLIF/JPEG-XL
territory requires either C implementation of arithmetic coding or
fundamentally different prediction models.

## The Tournament Theory Contribution

Tournament theory didn't give us a new compression ALGORITHM.
It gave us a FRAMEWORK for understanding compression:

1. **Score sequence = structural hierarchy** → color decorrelation
2. **A² = second-order structure** → gradient-adaptive prediction
3. **H(T) = scheduling freedom** → number of valid scan orders
4. **Isomorphism class = semantic equivalence** → code deduplication
5. **Prediction-entropy duality** → the fundamental trade-off

These insights guided us to the RIGHT combinations of techniques,
even though the individual techniques (MED, zlib, lzma) were known.
The theory told us WHY they work and WHEN to use each one.
