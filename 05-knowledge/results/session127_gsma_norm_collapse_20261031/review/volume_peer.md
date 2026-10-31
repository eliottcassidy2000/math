# Independent volume peer review

Reviewer lane: recover_lucas_47_triad. Scope: mathematical formulas and proof logic
in `volume/THEORY.md`, not a second execution audit of its checker.
Reviewed file SHA-256: `5e141fed011c9e21d76dfb85d1e356c2b64f7d8073341dde3e38dcb66eaad662`.

**Verdict: mathematically sound within its stated reconstruction; no load-bearing
gap identified. Three small precision edits are recommended below.**

## Checked

- Scaled complex coordinates `(sqrt(2) Re z,sqrt(2) Im z)` give complex norm
  factor `(u²+v²)/2`. The real/imaginary change contributes determinant modulus
  `1/2` per pair and the subsequent scaling contributes2, so the lattice
  covolume is `sqrt(abs(discriminant))`, not `2^(-s)` times that. Positive
  Euclidean Gram determinant is therefore `abs(discriminant)`; it is not the
  indefinite algebraic trace matrix.
- Applying AM-GM to the squared embedding magnitudes yields
  `abs(N alpha)<=(Q(alpha)/d)^(d/2)`, with the stated equality criterion.
  The Euclidean-ball Minkowski constants at degrees3,4,5 simplify to the
  displayed `2/(sqrt(3) pi)`, `2/pi²`, `12/(5sqrt(5)pi²)`.
- The binary Gram matrix `[[k,-1],[-1,S_k]]`, determinant `kS_k-1`, and
  coefficient ellipse area `pi R²/sqrt(kS_k-1)` are correct. Complex-pair
  comparison gives `S_k=Tr(theta²)+4 sum(Im theta_j)²`, so substituting the
  algebraic trace3 for the metric coefficient would indeed be invalid.
- For the tribonacci point, `3*103²-2*103*56=20291`; the remaining coefficient
  is `56²=3136`. The lower bound `S_3>106/25` gives
  `Q_3>839691/25>180²`, correctly demonstrating that norm magnitude1 need not
  mean a short Euclidean vector.
- The full norm-body infinite-volume proof uses two archimedean factors with
  real dimensions w1,w2. The reciprocal scalings preserve both their product
  norm and their total Jacobian. A positive-volume box away from zero factors
  has disjoint geometric translates, giving infinite measure, not merely
  unboundedness. All degree>=3 signatures have at least two factors.
- The binary-area substitution `x=ty`, Jacobian `abs(y)`, gives exactly
  `integral_R abs(f(t))^(-2/d) dt`, without an extra factor2. At a simple
  real root the exponent2/d is integrable for d>2; at infinity the integrand
  decays quadratically. Hence the binary slice may have finite area even
  though the full norm body has infinite volume. Its real-root rays and an
  interior disk rule out convexity. The note properly refuses to deduce
  lattice-point finiteness or effective heights from finite area alone.
- Homogeneity gives threshold scaling `M^(2/d)`; the embedded-plane Jacobian
  is the stated square root of the binary Gram determinant. Degree, rank and
  signature remain distinct throughout.

## Minor precision edits (not mathematical counterexamples)

1. VOL7: replace “with only simple real roots” by **“every real root of which
   is simple”**. The current wording can be read as demanding that all roots
   be real, whereas the intended proof permits nonreal roots and the degree
   3/4/5 examples have them.
2. VOL1: say explicitly that the Gaussian polar calculation is for n>=1;
   the n=0 case is the separately declared convention V0=1. The displayed
   radial integral is not literally valid at n=0.
3. VOL4: if the cited Minkowski version requires volume strictly greater
   than `2^d covol(L)`, add the usual closed-ball limit/discreteness sentence
   to justify the displayed non-strict radius bound. This is a routine valid
   repair, not an additional hypothesis.

The source/morphism architecture and historical `normcollapse` notation are
explicitly reconstructed rather than authenticated. This review does not
change that provenance status and does not cover any later edits to the
hashed version above.
