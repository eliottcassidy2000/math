# Independent audit: Carmichael power retractions and axis-precision clocks

**Verdict: PASS, with one wording clarification applied.**

Audit completed 2026-10-20 UTC. Scope is the staged
`session8-carmichael-frobenius/CARMICHAEL_POWER_RETRACTIONS.md` and its two
verifiers. This is not a re-audit of the cited repository canon or an
independent audit of the attached paper's analytic inputs.

Reviewed note SHA-256:
`16b4ac3a34342536769281fffbf0635f51f379a328619147181cf5c49068f508`.

## Findings and resolution

No mathematical defect was found in the stated scopes. Corollary 2.1's
phrase “monomial ideal over Z/p^a” could have included coefficient-torsion
relations such as pT. The author clarified that this flat corollary uses
monomials with unit coefficients and a finite quotient. The general
coefficient-torsion theorem separately covers pT and other staircases.

The odd-prime assumption on the exact arithmetic/axis theorem is essential
and is present. For example, on Z/8, the sixth-power map reaches its constant
Teichmüller target in one iteration for every input, whereas the naive
odd-prime formula `(a-1)/v_2(6)` would give two. Also, under squaring on Z/8,
the nonunit 2 takes two steps `2 -> 4 -> 0`, while odd units take at most one.
The note does not extend its sharp odd-prime arithmetic formula across this
boundary. Its flat nonzero-jet examples at p=2 are valid independent controls.

## Proof audit

### Nilpotent contraction and the fixed scheme

If x-y lies in I^j, the linear binomial term has an additional I factor
because e lies in I; every higher term lies in I^(2j). Hence each power
iteration gains one ideal order. Starting from F(x)-x in I proves that
F^(k-1) is a fixed section when I^k=0. The k=1 boundary is valid: I=0
and the assumed residue identity already says F is identity.

The section is a multiplicative set/monoid map. Its failure to be additive
at arithmetic precision n^k, k>=2, follows from the characteristic
obstruction and is correctly stated.

The fixed polynomial X^e-X is monic and its derivative is a unit since e
is nilpotent. Therefore its fixed locus is finite étale of degree e.
This does not count A-valued points as e, and does not make the original
power map étale. The note maintains these distinctions.

### Sharp arithmetic timing

For odd p dividing a Carmichael n, v_p(n)=1. At a unit residue the
principal-unit error therefore gains exactly one p-adic order per
iteration by LTE. At a zero residue the valuation is multiplied by n
instead. Taking the maximum of local hitting times gives the CRT clock.
The witness 1+n attains k-1 iterations simultaneously at all prime
components, so the uniform arithmetic bound is genuinely sharp.

### Multivariate torsion-staircase theorem

The monotonicity `gamma<=delta => w_gamma>=w_delta` makes the specified
coefficient ideals consistent and gives a genuine monomial normal form.
No cancellation between distinct monomials is silently allowed.

The coefficient estimate follows from the integral identity
`gamma_i c_gamma = E * [T^gamma](P^(E-1) T_i partial_i P)`.
Thus `v_p(c_gamma)>=v_p(E)-v_p(gamma_i)` before passing to the quotient.
For every surviving nonconstant monomial, a surviving p-power axis lies
coordinatewise below it with at least its coefficient precision. This
validates the upper bound using only axis precision.

For sharpness, if an axis requirement exceeds v_p(E), the selected
coefficient `j=p^min(h,v_p(E))` is within the exact binomial-valuation
range, even when r=0. Its valuation is below the available coefficient
precision. The witness 1+T_i therefore survives. The scalar witness 1+p
handles the independent constant-precision constraint. Together they prove
the minimum clock, including rings with no surviving variables.

The nonunit upper bound is valid: if v_p(E)>=a-1, then
`E>=p^(a-1)>=a`, so the E-th power of a p-divisible constant vanishes
modulo p^a. The a=1, r=0 case is also correct.

Consequently the exact formula
`ceil(max(a-1,H)/v_p(e))`
passes the proof audit for the stated odd-prime monomial quotients.
Mixed nonmonomial relations remain outside that theorem.

### Jets, extension fields, and carry normalization

The first surviving local Hasse term at a nonzero point is at degree p,
with coefficient n/p modulo p. At the zero point its first degree is n;
the note correctly keeps that distinction. The nonadditivity example at
jet length p+1 requires n>p and is valid for composite Carmichael n.

For n=561, independent cyclic-group arithmetic confirms:

* F_9: the power map is identity because 8 divides 560;
* F_27: it is a permutation since gcd(561,26)=1, with 3 fixed points;
* F_121: it has `1+120/gcd(561,120)=41` image points and
  `1+gcd(560,120)=41` fixed points.

The least noninjective extension degrees at 3, 11, 17 are 5, 2, 2.
The map's inseparable degree p and separable degree n/p are different
from the separability of X^n-X; the note does not conflate them.

The precision-two carry normalization uses the equality of the two
Teichmüller lifts a^n and a^p modulo p^2. Dividing the evaluated difference
by n contributes the factor `(n/p)^(-1)` modulo p. Its cocycle and
nonintegral formal-polynomial boundary are correct.

## Independent controls

The audit imported no author helper functions.

1. `independent_staircase_check.py`: 232,544 checks across 176 complete
   one-variable coefficient-precision rings. Includes odd primes, p=2
   boundary controls, exponents 2, 3, 4, 5, 6, 9, 21, 45, arithmetic-only
   rings, torsion staircases, and every element in each stated ring.
2. `independent_multivariate_check.py`: 97,710 checks across 62 complete
   multivariate rings. Includes squarefree mixed monomials, disappearing
   mixed terms, cubic p-power axes, three variables, nonflat coefficients,
   and independent measurement of the worst first hitting time.

Total independent checks: **330,254 across 238 rings**. Both audit scripts
pass in normal and optimized Python, with identical JSON outputs. Explicit
exceptions retain all checks under optimization. These finite controls
corroborate the proofs; they do not extend their quantified scope.

Both author verifiers were also rerun during this audit. Their outputs
match `verification.json` and `axis-verification.json` byte for byte.
The first author verifier reports 2,063,548 checks. These replays are not
counted as independent algorithms or added to the independent total.

## Reproduction

From `session8-frobenius-review/`:

`python3 independent_staircase_check.py`

`python3 -O independent_staircase_check.py`

`python3 independent_multivariate_check.py`

`python3 -O independent_multivariate_check.py`

The two `author-*-replay.json` files preserve the audit's reruns. No
checkout or Git files were modified by this reviewer.
