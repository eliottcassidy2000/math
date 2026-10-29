# Audit, proof dependencies, and exact failure boundaries

2026-10-28. This is a **same-agent adversarial audit with alternate computational
paths**, not an independent-agent audit. No canon theorem ID was reserved or
promoted. All all-size statements have written proofs in this packet;
bounded computations are controls on them, not replacements for them.

## 1. Proof dependency audit

- The general prime-degree fixed-point-gap result depends on **CITED**
  Jones Theorem1.2(k=0), whose proof uses finite-simple-group classification.
  Every projective exception is checked: prime projective degree forces
  d prime; d=2 has even q; odd d has a large fixed subspace; semilinear
  involutions require quadratic descent. The norm normalization and fixed
  representative of a stable projective line are included in the proof.
- The Lucas fixed-word count is a folded path count. Its strict inequality
  is the exact Fibonacci identity, not an asymptotic estimate. Group parity
  is computed modulo4, with3|n excluded only after assuming primality.
- The odd golden reader is **inherited PROVED**, from Session74's independent
  all-period audit. Even words have three zero representatives and are not
  silently put in a field of size L_n. The even-index group theorem uses a
  separately chosen cycle and states that distinction.
- Full holonomy follows from the actual Schreier generators of the marked
  action. The zero-prefix rotations generate the full native wreath group,
  so C intersect <K>=C. The repair word is expanded into literal marking
  triangles and every edge is charged. It is not a native-complement theorem.
- Joint-port counts concern **ordered tuples allowing repetition**. The
  pair-distinct formula subtracts the diagonal. Linked child actions are
  one diagonal A action, never two independent actions. Burnside moments
  use uniform automorphisms, not uniform sources in a different problem.
- The finite equidecomposition theorem specifies singleton pieces and the
  acting group. It is neither equivariant G-set isomorphism nor Euclidean
  scissors congruence. Classical geometric Dehn assertions are not imported.
- Prime-divisor records distinguish order from valuation. The lifting proof
  uses odd r!=5, an unramified quadratic algebra and unit roots. Density2/3
  is imported from Moree for prime **divisors of some term**, not for indices
  producing prime-valued terms. No independence assumption is used.

## 2. Producer universe and active controls

`check.py` reports **311033 active checks**, unchanged under `python3 -O`.
It has no Python assertions.

- All cyclic legal binary words at lengths3 through19. Reflection counts,
  exact golden label bijections at odd prime-valued lengths in that range,
  involution parity, the660-element eleven group and the inherited A29
  prime-cycle witness are checked. Fibonacci identities are tested through
  m=499, and even sign identities at k=3,...,9. Prime examples47 and2207
  are verified by trial division.
- Every normalized projective involution in PGL3(F2) and PGL3(F3), giving
  respectively21 involutions with3 fixed points and117 with5. All120
  semilinear projective maps on P1(F4), including15 involutions fixing1
  point and10 fixing3. General lower-bound polynomials are checked for
  q=2,...,32 and the declared small dimensions; these checks do not prove
  or replace the general linear/semilinear argument.
- Every odd prime<=5000 except5:667 primes. Projective order by powering
  in the quadratic algebra is compared with a separately iterated Fibonacci
  rank. The exact Lucas midpoint/valuation rule is checked for n=1,...,120.
  Of these primes445 have even rank; this finite ratio is not used to prove
  a density theorem. The six-mirror character table is tested throughout,
  with full projective group enumeration at r=5,7,11,13.
- Complete marking orbits at W1,W2,W3:3,27,6561 marks. All Schreier triangles
  from the zero-prefix generators are verified as literal transports; their
  generated holonomy groups have orders1,3,243. Every correction in these
  groups is separately reconstructed by native local rotations and by a
  based marking loop. The longest constructed generator words are0,6,38;
  these are finite maxima for this construction, not optimal loop lengths.
- Exact joint-port moments through r=8 and depth12. Burnside sums over the
  actual retained centralizers agree through depth3; finite differences
  check the predicted polynomial degrees and leading terms. The all-r,
  all-depth polynomial proof is written in `PORT_MOMENTS.md`.
- All76 period-nine states; the F4 projection, the four-way F19 ambiguity,
  and every one of72 nonzero chart/phase reconstructions. A separate
  degree76 full cycle and involution preserve nineteen blocks of size4,
  demonstrating exactly why primitivity cannot be dropped.

Hostile controls include dihedral groups with one fixed point, the degree11
exception, the composite degree76 imprimitive construction, even golden
zero-fibre loss, two nonconjugate projective involutions at11, and the order9
joint mod3 port charge.

## 3. Alternate computational paths

`audit_alternate.py` imports no producer code. It constructs centralizers
directly from the recursive block formula and compares the full permutation
set hashes, not merely their orders, with the loop-generated subgroups.
It also generates legal words as binary tuples and obtains exact group
orders at11 and29 using SymPy1.14.0.

Attempting general Schreier-Sims enumeration at the larger degrees was
stopped without accepting its unfinished output. The replacement is stronger
as a small reviewable certificate:

- At p=47, using lexicographic tuple labels and the edge mirror R,
  R T^12 has cycle lengths1,2,2,2,3,11,26. Its78th power is an11-cycle.
- At p=199, using the golden root62 and R:i->-i,
  R T^31 R T^180 has cycle lengths19,24,37,59,60. Its84360th power is
  a59-cycle, fixing140 points.

Primality, the actual cycle decompositions, label bijections and reflection
parities are explicitly checked. **CITED** Jordan's theorem (Jones
Theorem1.1) then supplies A47 and S199, using a different proof route from
the all-prime fixed-point classification. These finite certificates do not
establish the theorem for every later prime on their own.

## 4. Reproduction hashes and documentation boundary

Run `python3 reproduce.py`. Both scripts reproduce their saved JSON exactly
in normal and optimized modes. Receipt SHA256 values:

    receipt.json
    2b25608beaab5c2a2dd39fe8e195164440874249b098abd96f92f17e17082baf

    audit_receipt.json
    74e0ad718a64bc14754de1f6875bdccc3d3f4777139f56200458a79b8a27ce84

The alternate script requires SymPy; the producer uses only the Python
standard library. Tests print consequence objects: full subgroup hashes,
port counts, exact prime-cycle certificates and ambiguity counts. No Lean
formalization, unrestricted tournament census, or independent-agent audit
is implied. No demonstrated false published claim arose in this packet;
inherited corrected near misses are linked in README rather than overwritten.
