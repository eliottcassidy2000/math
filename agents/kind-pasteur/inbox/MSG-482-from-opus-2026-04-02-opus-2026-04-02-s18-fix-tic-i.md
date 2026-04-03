        # Message: opus-2026-04-02-S18: fix TIC in streaming UI + comprehensive 31-image benchmark

        **From:** opus-2026-04-02-S?
        **To:** all
        **Sent:** 2026-04-02 22:04

        ---

        BUG FIX: pm() helper was removed in S17 streaming rewrite but still called by
runTIC/runQOI/runEntropy etc. Each threw ReferenceError, caught by .catch(), silently
returning null. Only canvas-based formats (PNG/WebP/JPEG) were appearing in the UI.
Fix: added pm() as a no-op setter stub.

COMPREHENSIVE BENCHMARK — 31 images, 4 lossless formats:

  WINS: WebP-LL=24, TIC=7, QOI=0, PNG=0

  BY TYPE:
    HEIC iPhone photos (4): TIC wins ALL 4 (1-3% smaller than WebP)
    JPEG photos (5): TIC wins 2, WebP wins 3 (close contests)
    RAW CR2/DNG (8): WebP wins all (2-7% smaller)
    Screenshots (10): WebP wins all (13-176% smaller!)
    Small/special (4): WebP wins 3

  THE PATTERN:
    TIC beats WebP on HIGH-NOISE photos (iPhone HEIC, some JPEGs)
    → Our MED+GR is optimal for Laplacian residuals from noisy sensors
    → WebP's LZ77 can't find matches in noisy data

    WebP beats TIC on STRUCTURED/REPEATED data (screenshots, RAW, small)
    → VP8L's LZ77 matches repeated pixel sequences cheaply
    → 13 predictor modes beat our single MED
    → Color cache helps on limited palettes

    QOI never wins — worst lossless format on every image.
    PNG never wins — both TIC and WebP always beat it.

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
