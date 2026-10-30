# Session 119: odd cycles, odd functions and arithmetic transport

**PROVED scoped classifications and exact synthesis; no major-conjecture resolution.**
Oddness is not itself a transfer principle. A specified action and the obstruction
left by transporting a value around a closed orbit give the common theorem.

- [Affine transport and odd eighths](synthesis/TRANSPORT_AND_EIGHTHS.md): a cycle
  closes exactly when (id-M)t=B. The four odd eighths have an exact character
  expansion; odd-weight Boolean functions and reflection-odd functions are kept distinct.
- [Signed cycles](cycles/SIGNED_CYCLE_INTEGRAL.md) and
  [signed graphs](cycles/SIGNED_GRAPH_INHOMOGENEOUS.md): complete integral Smith
  forms, inhomogeneous solvability and modular fibers. Negative cycle transport
  creates a factor-two obstruction; free compatibility equations remain separate.
- [Reflection-odd multiplier lattices](odd_functions/ODD_LATTICE_CLASSIFICATION.md):
  all divisor strata and signed quotient cycles, complete Smith factors for I+P
  and I-P, and a four-way prime-order valuation criterion.
  [Independent audit](audit_operators/AUDIT.md).
- [Weighted arithmetic cycles](arithmetic/WEIGHTED_CYCLE_OBSTRUCTIONS.md): the
  same cyclic cokernel theorem contains odd-cycle sign equations and closed Collatz
  word equations. The determinant and exponent multiset lose the ordered carry;
  an explicit negative-cycle hostile distinguishes integral from nonintegral words.
  [Independent audit](audit_arithmetic/AUDIT.md).

[Inheritance](INHERITANCE.md) names the actual prior theorems and corrections.
Classical mechanisms and previously established homogeneous formulas are credited,
not presented as new discoveries. No Collatz convergence, positive-cycle exclusion,
Ramsey record, LRC closure or unrestricted circuit lower-bound record is claimed.

## Complete arithmetic tournament families

- [One fixed character](characters/CHARACTER_STRATA.md): conductor-divisibility
  classifies all allowed strata and their signs. Full support forces a prime power
  p^a with p=3 mod4. [All valuation frames](characters/FRAME_ISOMORPHISM.md) form
  one tournament isomorphism class; an explicit wreath subgroup gives Hamiltonian
  divisibility. Ring residues and finite-field Paley constructions are kept distinct.
- [Different characters by stratum](characters/LOCAL_CHARACTER_ESCAPE.md): full
  support exists exactly when every prime factor is3 mod4, with a uniform count.
  At N=21, all sixteen members form exactly four isomorphism classes within this
  family, separated by spectral moments and checked cycle counts.
- [All square-invariant odd symbols](characters/SQUARE_QUOTIENT_CLASSIFICATION.md):
  a complete quotient-function parametrization gives 2^C symbols, where
  C=(product_p(1+2a_p)-1)/2 for N=product_p p^a_p with all p=3 mod4.
  Character formulas suffice exactly up to two distinct prime factors. Three
  factors admit genuinely nonlinear odd quotient functions. The first eligible
  modulus231 has8192 symbols versus4096 locally-character symbols, with a literal
  majority-type example. These are labelled counts, not8192 isomorphism classes.
- [Independent character audit](audit_characters/AUDIT.md) and
  [square-class audit](audit_characters/SQUARE_LIFT_AUDIT.md) verify these boundaries.

[Research synthesis](RESEARCH_SYNTHESIS.md) identifies precise routes toward
extremal and algorithmic advances, and the counterexamples that constrain them.
The extra square-class symbols have not been shown to improve a Ramsey bound.

## Reproduction

Run `python3 verify.py` in this directory or by its repository-relative path.
It checks the sealed file hashes and runs all twelve standard-library Python
checkers normally and with -O, requiring identical outputs and expected receipts.
All twenty-four runs passed during publication. Finite universes are stated in
individual proofs and audits. They corroborate the uniform arguments; they are
not Lean builds or independent human reviews. No external theorem submission or
unrequested external communication was made.
