# Independent audit: characteristic-four norm shadow and convolution units

**Verdict: PASS**, subject to the exact stated carriers and thresholds. No mathematical correction is required. The finite checks support, but do not replace, the dimension-uniform proofs. This audit does not claim a new solution of a major open problem or literature priority.

## Scope and proof review

Reviewed `../arithmetic/ARITHMETIC_SHADOW.md` and `check_arithmetic.py`, and the inherited Session112 `norms/NORM_SIGNS.md` at the isolated Session115 baseline. The inherited dimension obstruction is an explicit five-equation signed-XOR certificate, not an appeal to an unproved extension of Hurwitz's theorem. Normalization forces the negative diagonal and oriented-line form whenever real norm composition holds. The five equations then sum to `0=1`; restriction to a four-dimensional binary subspace covers every d≥4.

1. **Integral defect ideal.** Expanding the bilinear product squared leaves only mixed monomials supported on XOR rectangles. Each coefficient is 0 or ±4. A surviving monomial is isolated by setting its two row and two column coordinates to1 and all others to0: the defect is exactly that coefficient. Hence coefficient and integer-value ideals agree and are either(0) or(4). This does not require associativity or normalized signs.
2. **All commutative coefficient rings.** Formal sufficiency follows from the coefficient ideal. Necessity follows by the same actual sparse evaluation, so finite rings do not introduce polynomial-function exceptions. For an incoherent table the condition is exactly `4·1_R=0`, including characteristic4, characteristic2 and the zero ring. The modulus statement includes m=1.
3. **Parity norm.** Squaring any integral lift modulo4 depends only on its parity. Thus reducing signed multiplication modulo2 and taking Hamming weight modulo4 gives a well-defined multiplicative map. It is not additive and is not octonion multiplication. A modulo8 lift cannot in general be reconstructed from parity, because even entries may square to4 modulo8.
4. **All units.** In the commutative group algebra, every cross term vanishes in characteristic2 and each basis element squares to the identity. Therefore `a*a=aug(a)·δ₀`. Odd augmentation is sufficient for invertibility, and even augmentation makes the element square-zero and noninvertible. Thus the unit group is an elementary abelian2-group of order `2^(2^d−1)`.
5. **Extreme-Walsh subgroup.** For d≥3, the affine-function ideal has square zero: every relevant intersection of two affine equations has an even number of solutions. Translation of affine functions supplies the cross terms. In parameters `f(v,r,c)(x)=c+r·v+r·x+δ_v(x)`, the resulting law adds v,r,c. Distinct defect points cannot differ by an affine function because the difference has weight2 whereas nonzero affine weights are `2^(d−1)` or `2^d`. This establishes uniqueness, group size `2^(2d+1)`, and index `2^(2^d−2d−2)`.
6. **Character.** The affine part takes value c at the defect point. Its weight is divisible by4 for d≥3, so toggling the point changes weight by `(-1)^c` modulo4. Both character fibers have `2^(2d)` members. H₃ exhausts all units; H₄ has index64. This is a convolution group, not the pointwise additive sign-gauge group.
7. **Sharp unit-weight ideals.** At d≤1 units all have weight1. At d2 the only unit weights are1,3; a product of two weight3 units has weight1 by the modulo4 character, giving defect−8. For d≥3, supports `{0,e₁,e₂}` and `{0,e₁,e₃}` multiply to `{e₂,e₁+e₂,e₃,e₁+e₃,e₂+e₃}`, producing defect−4. Embedding preserves this witness in all higher dimensions. Hence the complete ideals are(0),(0),(8),(4),(4),….

## Boundary and information-loss checks

- Below d3, the affine square-zero ideal and unique phase parameters fail. At d2 the two independent coordinate functions convolve to the constant-one function. No negative-exponent index formula is being asserted there.
- The modulo8 **unit-weight** failure starts at d3 (vector-space dimension8). The real signed-XOR **norm-composition** obstruction starts at d4 (dimension16). These are distinct statements, and the manuscript explicitly separates them.
- The divisor-sum statement applies only to the specified seven-prime bank and orientation map. Modulo4 makes the norm-composition predicate insensitive to its syndrome; it does not assert that arithmetic squareclasses themselves cease to exist or that all integer selections have square divisor sums.
- The classification of the extreme-Walsh family as one-point affine defects is inherited. The convolution law, exact unit-group index and weight character are the new structural assertions under audit.

## Independent computation

`check_independent.py` imports no producer implementation. Its independent routes are:

- direct formal polynomial expansion of every normalized dimension4 signed table (512 tables), checking **every** surviving sparse witness, rather than using the producer's rectangle coefficient loop;
- direct formal expansion of all128 dimension8 oriented-line tables and exact integer square testing of the prescribed sigma bank;
- direct reconstruction of the five contradictory line equations from their rectangle coordinates;
- tuple-based convolution of all pairs for d=0,1,2,3, including nonunits, and the entire unit-defect gcd census;
- every generator multiplication on every phase in H₃, H₄ and H₅ (sizes128,512,2048), independent of the producer's full pair-table strategy; associativity then certifies the full group law on these finite cases;
- Walsh peaks in d3,d4; exact embedded hostile supports through d7; the d2 affine-product hostile.

Run from the session scratch directory:

```sh
python3 session116_research/audit_arithmetic/check_independent.py > /tmp/arithmetic-normal.json
python3 -O session116_research/audit_arithmetic/check_independent.py > /tmp/arithmetic-optimized.json
cmp /tmp/arithmetic-normal.json /tmp/arithmetic-optimized.json
```

The saved `normal.json` and `optimized.json` must be byte-identical. `producer_rerun.json` records a separate reproduction of the producer checker and must agree with its saved `results.json`. The independent script uses explicit exception gates, with no `assert` dependence.
