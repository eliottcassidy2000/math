# Independent geometry/arithmetic review

Reviewer: Session292 manuscript lane. Scope: `square_sum/PROOF.md`,
`arithmetic/PROOF.md`, their reproducibility scripts, and independent controls.
This is not a review of the reviewer's own manuscript lane or the whole repo.

**Verdict: PASS WITH CORRECTION**, after the simple-cycle boundary below was
repaired. All infinite assertions reviewed here have hand arguments; finite
controls corroborate their explicit arithmetic and boundary cases.

## Correction found and repaired

The first all-word cycle statement/script accepted the length-one word(4) at
start2,N2 as a simple cycle, and length-two backtracks likewise. Both comparison
implementations shared this convention, so their agreement did not catch it.
The worker added the standard length>=3 condition to the theorem and both
implementations and added explicit hostile controls. No hexagon, generating
function, or first-birth formula changed.

## Square-sum reasoning checked

- Prefix affine formulas and intersection of interval guards are exact.
  Opposite-parity collisions occur at one rational start, deleted only if
  integral; same-parity collisions are identical or impossible.
- For ABCABC the complete cross-parity collision set reduces to
  `0,±p,±q,±(p+q)` in y. The y↔−y identification is exactly the two starts
  with the prescribed oriented colour word on one unoriented cycle.
  Reversal changes ABC to ACB, so no extra factor of two is warranted.
- The birth polynomial uses maximum vertex `(B+C−A+|y|)/2` and retains
  parity; its at-most-three holes are a set, hence p=q is counted once.
- Consecutive-square threshold and the triangle/hexagon parity split check
  algebraically. First hexagon23 and first triangle30 arguments retain their
  stated word-family scope and exhaust their respective bounded square triples.
- Independent enumeration identifies cycles by undirected **edge sets**, not
  the producer's dihedral canonical tuples:9,555 target-triple/size checks
  pass. Normal and optimized producer runs pass15,288 word checks,69,000
  hexagon checks and97 consecutive-square boundary checks.

## Arithmetic reasoning checked

- Redundant29/37 clauses lie inside the5 clause. All active discrete-log
  classes and shared-core/petal clocks agree with modular computation.
- A least period of the one-sided sequence divides the displayed period:
  extend by that period to all integers, then use gcd of periods. A proper
  divisor divides M/q for some prime q|M, so the seven opposing-value
  certificates establish minimality. The k0 witness is shifted into the
  stated positive domain by M.
- The conditional product is legitimate only after conditioning on b mod360.
  In the general statement, g_i|C and gcd(C,q_i)=1 imply that C/g_i is
  invertible modulo q_i; pairwise coprimality of q_i then gives CRT independence.
  Pairwise intersection of arbitrary guard sets alone does not do this.
- Hidden-factor induction is sound: every forced prime is coprime to the
  base; pinning exponent modulo its order at p^(e+1) fixes its valuation at e.
  Unbounded growth on the retained infinite progression forces a new prime.
  Retain nonnegative k so f(k) is an integer. The author incorporated this small domain clarification
  as k_min>=0. Factorization by trial division makes the search
  effective, regardless of efficiency.
- Full61,209,720-residue sieve, generalized-CRT inclusion/exclusion, and the
  factored360-cell computation agree. Normal/optimized full-period output
  agrees; this is not a primality test or global Chan–Pak computation.

## Additional independently derived control

`99*2^34+23 = 149*347*2111*15583`. Their base2 orders are respectively
148,346,1055,1113. Thus the entire progression

    k=34+30670763342539320*t, t>=0

survives the first25 prime filters and has these four distinct prime factors.
This is a concrete check of the hidden-factor mechanism, not its all-order
proof. `independent_geometry_arithmetic.py` checks primality, exact factorization,
orders, and ten huge-exponent modular samples without constructing huge integers.
Normal and optimized outputs agree.

## Reproduction

From repository root:

```sh
python3 05-knowledge/results/session292_triple_bridges_20261113/square_sum/verify.py
python3 05-knowledge/results/session292_triple_bridges_20261113/arithmetic/verify.py --full-period
python3 05-knowledge/results/session292_triple_bridges_20261113/audit/independent_geometry_arithmetic.py
```

Repeat with `python3 -O`. The final coordinator manifest should pin reviewed
files after any later prose or receipt changes; a stale hash is not a final
source audit. No assertion here supplies a physical tournament addition move.
