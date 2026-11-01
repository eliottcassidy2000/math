# Session 133: executable tournament sequences and sharp H-spectrum certificates

**Correction to checkpoint 1:** [the atlas hostile test confused trailing-zero
padding with distinct polynomials](spectrum/CORRECTION.md). It is repaired;
all graph counts and mathematical inventories are unchanged.

**PROVED program/structural results; FINITE-EXACT specified replays;
CITED classical strong bounds; unrestricted spectrum completeness OPEN.**
No literature-priority or smallest-program claim is made.

## Actual finite programs

There are **sixteen explicit ordered FRACTRAN lists**, not just algorithms
claimed compilable by universality. Their proofs specify valid initial
integers, output prime exponents, first-eligible order and termination or
streaming checkpoints. Detailed execution ranges distinguish native formula,
primitive-register, fraction-exponent and literal-integer runs.

- [Eight arithmetic counters plus two small core counters](fractran/README.md):
  A000568, self-converse, converse-merged, strongly connected, ordered pairs,
  two-colored pairs, two-colored tournaments, and core-floor extremizer counts.
  The A000568 compiler falls from5543 to4294 fractions; source and compiler
  equivalence are tested independently. The core floor has20/14-fraction lists.
- [Three positive multiplicity streams](distribution/README.md):397 fractions
  from order0;67 fractions from order6;43 fractions from order6 using proved
  equal-sequence aliases. Each gives all eleven N_n(h), h=1,3,...,21.
- [Three path/word programs](path_counter/README.md): a general encoded
  tournament evaluator, a32-fraction one-apex word machine, and a square-shell
  counter for the converse-fixed fibre theorem below.
- [Compiler and arithmetic proofs](fractran/THEORY.md), [linear stream
  proof and exact cost](distribution/THEORY.md), [word/adjacency proofs](path_counter/THEORY.md).

## Structural results supporting the programs

1. **Sharp finite local certificate at every threshold.** If r is the least
   integer with3^r>B, every induced-minimal obstruction to H<=B has at most3r
   vertices; r ordered directed triangles attain this bound. Apart from that
   triangle chain, the sharper bound is K(B)+1.
2. **Exact H<=21 basis.** There are56 minimal forbidden induced tournaments:
   29 strong order6,13 strong order7,13 nonstrong order8 and one order9.
   Exact forbidden-basis sets agree at thresholds5/7 and19/21. This recertifies
   the already-known omissions7 and21; equality of cardinalities alone would
   not prove it. [Proof, sequence and explicit adjacencies](spectrum/README.md).
3. **Exact reduced-core budget and extremizers.**
   K(B)=max{3a+4b:3^a5^b<=B}; it suffices to test b=0,1,2.
   In particular K(7)=4 and K(21)=7. The dual minimum g(N) satisfies
   g(N+3)=3g(N) for N>=6, and all equality cores are classified. A capped-DP
   corollary recognizes H<=B on a supplied tournament efficiently in B,
   independently of arbitrarily long transitive padding.
   [Core proofs and precise complexity](core_reduction/THEORY.md).
4. **Converse and square shells in an exact word family.** For pointed strong
   one-apex tournaments, let R_h count all isomorphism classes with H=h.
   Exactly one odd R_h occurs in each shell q²<h<=(q+1)². Its location is
   G(q)=q²+(q mod2)+2f(q)-1, where f is the inherited digital one-apex counter.
   Forgetting the apex requires an explicit correction; H5 is the hostile.
   [Fibre theorem, correction and actual counter](path_counter/THEORY.md).

## Inheritance and the open target

The all-order exclusions7 and21 are inherited from
[THM-1370, exact tournament slug](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
which credits the2016 literature. Fixed-H rational formulas and the digital
inverse are also explicitly inherited. The new work concerns executable
refinements and sharper structural certificates, not rediscovery claims.

**Every other positive odd number occurs as H(T)** remains a conjecture.
The one-apex family itself misses13; finite verified coverage and a general
counting program do not prove universal realization. [Proof bridges and
remaining obligation](SYNTHESIS.md).

[Verification and exact replay coverage](VERIFICATION.md) ·
[Independent audit](independent/README.md) · [Concept board](BOARD.md).
