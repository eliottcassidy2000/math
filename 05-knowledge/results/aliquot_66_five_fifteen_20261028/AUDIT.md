# Verification and scope review

2026-10-28. Same-agent proof review plus independent computational
formulations. **No independent-agent audit is claimed.** No scarce theorem
or hypothesis ID was reserved or promoted.

## Mathematical checks

`check.py` uses only the Python standard library and explicit exceptions;
none of its checks disappears under optimization. The normal and `-O`
outputs agree byte-for-byte with `receipts.json`: **68,305 checks**.

Declared universes and independent paths:

- A literal additive divisor sieve covers every1<=N<=200,000. Independent
  factorization verifies σ for every N<=5,000. All2,242 square-sigma
  sources in the larger window are tested under every prime multiplier
  p<=97 by exact factorization of the product, even beyond the source
  window. There are1,388 accepted edges. The universal classification
  remains the inherited corollary of Ljunggren, not a finite inference.
- Both complete exponent1/2 cubes on supports{2,3,11} and{2,5,7} are
  recorded with all sigmas, square classes and balance defects. Their
  three promotion signatures have binary rank3. Literal proper divisors
  independently check the balance formula through1,500.
- All mixed balanced square-sigma sources through200,000 are listed.
  The minimum364 has a complete smaller-source check; the fixed-pair
  parameter identity is checked at every listed source. This is an
  inherited result, not new priority.
- HNF enumeration checks all content strata, total σ and primitive ψ
  for every N<=500. Two-stage lattice chains are enumerated at p2,3,5,7
  to check every endpoint multiplicity independently of the product law.
- Lagrangian planes are independently enumerated as isotropic spans over
  F2 and F3, giving15 and40. The shear witness supplies the all-prime
  action obstruction in the proof, not just those two finite counts.
- The sigma144 fibre is exhaustively enumerated under its proved bound.
  Modular genera are calculated with exact fractions; elliptic-point
  counts are obtained by literal polynomial roots modulo N, separately
  from the Legendre-symbol product formulas used in the proof.
- Primitive reduced forms are exhaustively enumerated at discriminants
  −20,−60,−84,−252,−264,−280,−7392, using a<=sqrt(|D|/3). The actual
  representatives and self-inverse counts are saved. Class-group types
  at−264 and−280 follow from the proved counts and finite-abelian-group
  classification, not from an assumed numerical class-group routine.
- The all-k form formula is checked for1<=k<=1,000. Its actual presence
  and distinct inverse are checked by complete form enumeration for
  every prime q=12k²−1 with k<=30. This checks the algebra; the all-size
  theorem follows from its explicit discriminant, gcd and strict
  inequalities, so no finite bound is used in place of a proof.
  The stronger parent family8t²−6 is checked through t1,000, including
  complete reduced-form checks at its exceptional t1 and small t2.
- All15 marked K6 matching images are checked; all binary vertex-sign
  assignments on Petersen reject a bipartite gauge. The literal odd
  cycle and the corrected height-cycle normalization are also checked.

## Proof audit: obligations and hostile examples

P1 retains the ambient lattice and quotient type. HNF is not an enumeration
of abstract finite groups. The content decomposition has a literal inverse.
P2 uses the genus formula for X₀, not X₀^+, and records the elliptic terms.
P3/P4 use proper classes of the specified **order**, not an unmarked field.
The strict inequalities in P4 avoid every ambiguous reduced-form boundary;
both signs of b are canonical and distinct. For prime q the discriminant
is fundamental with three prime-discriminant factors, giving genus order4.
The non-idoneal conclusion itself does not require that extra prime guard.

P5's obstruction is expressly for a target action factoring through
reduction modulo p. It does not deny arbitrary set bijections or other
actions with extra level data. P6's gauge changes the label domain and
does not assert preservation of positive square-sum or absolute graceful
labels. P7 retains the actual column dependence, not only its zero total.

Positive and hostile controls:

- 3 and22 legitimately combine to66;3 and81 do not combine by the
  coprime rule. The prime3 is the unique zero prime column.
- 364 shows that square sigma and FUS balance coexist; their conjunction
  cannot be globally excluded by the promotion obstruction.
- k5 gives composite q299, so the uniform form theorem survives but the
  prime-based sigma formula fails.70 gives a square-sigma idoneal source
  outside the ray and blocks a universal extension of P4.
- The class numbers along the prime-filtered ray are not monotone:
  h=56 at k6 and h=32 at k7. Only the lower bound is proved.
- The integral shear is trivial modulo p but moves an index-p³ lattice.
  The fibre at pZ² has p+1 intermediate paths; uniform endpoint averaging
  is therefore wrong even at p2.
- The Monster/Baby ratio is not the actual involution-class index. Replacing
  B with2.B changes the small-prime height cycle, as checked explicitly.

## Historical repair and boundaries

The March `five-as-bridge.md` lacked the scope warning already present on
related reflections. This session adds a warning distinguishing the valid
curvature and rational-linear-independence calculations from unsupported
global identifications and algebraic-independence wording. The mechanism
is recorded in `01-canon/MISTAKES.md`; the historical body is preserved.
No other historical source is silently promoted.

The mathematical organizer is the existing guarded-context/weighted-fibre
framework, not a newly claimed equivalence among all listed structures.
The all-ray non-idoneal exclusion is exact and unconditional. Its transfer
to a Collatz positive-current estimate, a geometric six-colour forcing
context, or a planar Jacobian construction is **OPEN** and has no supplied
map preserving the respective acceptance condition. The session's
limited literature searches establish no novelty claim.

Documentation validation (`python3 agents/check_docs.py`) passes; local
Markdown targets resolve against the current files and pinned incoming
Git tree; `git diff --check` passes. Final Git synchronization is verified
after committing and reported in the session close-out.
