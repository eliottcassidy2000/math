# Devil's Advocate Session S3: Honest Assessment

**kind-pasteur-2026-03-25-S3**

## Compression: What's Real vs What's Inflated

### The "100% beats PNG" claim is MISLEADING because:
1. **Synthetic test images only** — no Kodak, Silesia, Canterbury, no real camera photos
2. **Unfair backend** — brotli/zstd vs PNG's Deflate. With same backend (zlib-9): 90% win rate
3. **Header trick** — 4 bytes vs PNG's 57. Significant for tiny images, negligible for real images
4. **No novel algorithms** — MED (1998), GAP (1996), Paeth (1996). Every component is 25+ years old
5. **Python zlib is 7% worse** than libpng's C zlib on identical data

### What IS real:
- **Lossless roundtrip VERIFIED**: 10/10 images pass compress→decompress (newly verified S3)
- **Fair comparison (zlib-only): 90% win rate** — genuinely better prediction than PNG's 5 filters
- **Transpose trick** is a modest but real innovation for vertical patterns
- **Ring codec** (from prior session) is genuinely novel scan order — never properly benchmarked

### What should be done:
1. Test on Kodak 24-image corpus (the standard)
2. Compare against JPEG-LS, WebP, JPEG-XL — the real competition
3. Implement in C for honest speed comparison
4. Properly benchmark the ring codec (the one genuinely new idea)

## Math: What's Real vs What's Inflated

### THM-260 (band-limitedness): Repackaged existing results
- THM-076 (Walsh-OCF factorization) already proved the Walsh degree bound
- THM-259 already stated the formula with proof sketch
- THM-260's only genuinely new content: correcting MISTAKE-034 (n=5 not band-limited)
- The "one-line inequality" n-1 < (n-1)(n-2)/4 for n≥6 is freshman math

### THM-261 (SC orbit pairing): Straightforward observation
- Z₂ action from anti-automorphism is a direct consequence of definitions
- Does NOT prove SC Maximizer conjecture (only describes mechanism)
- OPEN-Q-016 remains fully open

### Seesaw: NOW statistically significant at n=8
- **3000 samples: 0 violations, expected 6.4 under independence, p=0.0017**
- This IS a real result — the seesaw is not sampling noise
- At n=9: 500 samples insufficient (p=0.61)
- Algebraic proof still needed

### The defect wave: Real pattern, needs more samples at n=9
- beta_1: 30% → 14.6% → 5.8% → 1.1% → 0.4% (n=5 to n=9)
- beta_3: 0% → 1.2% → 8% → 18.8% → 24.4% (n=5 to n=9)
- Crossover at n=7 is robust
- n=9 sample size (500) adequate for beta_3 rate, insufficient for beta_1

## The Math-Engineering Gap

The devil's advocate found: **zero theorems led to codec improvements.**
The prediction techniques in the codec (MED, GAP, Paeth) were all discovered
empirically 25+ years ago. Tournament theory explains WHY they work
(Krawtchouk band-limitedness) but does not produce better ones.

The ring codec is the one place where mathematical thinking (concentric
expansion, center-first progressive decode) produced a genuinely new
engineering idea. It deserves proper development.

## What We ACTUALLY Accomplished (honest list)

1. **MISTAKE-034 correction** — real value, prevents error propagation
2. **Seesaw significance at n=8** — p=0.0017, first statistically rigorous confirmation
3. **Lossless roundtrip verification** — proves the codec is correct
4. **Fair benchmark numbers** — 90% win rate with same backend is honest and real
5. **beta_3=2 at n=8 and n=9** — confirms the rare high-beta_3 phenomenon
6. **Defect wave to n=9** — extends known pattern with larger samples
