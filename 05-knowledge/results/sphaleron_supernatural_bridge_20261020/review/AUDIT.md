# Independent audit: sphaleron winding and supernatural completion

Date: 2026-10-20 UTC. **Verdict: PASS within the stated categories and hypotheses.**
No blocking mathematical error was found. No change to the source packet was
made. This is an independent mathematical/source/code audit, not a proof of
new electroweak field equations or a Lean verification.

## Reviewed objects

- `session9-sphaleron/SPHALERON_STEINITZ_BRIDGE.md`, SHA-256
  `f9f9df810e3ca7070c0c286c5c574df88c2f13e98448e2c6ed5483c34e6443a4`.
- `session9-sphaleron/verify_bridge.py`, SHA-256
  `92f786c2fcf91a1ca3568c24456d435c798207d16100cb2d8ba5539c15e9e6ac`.
- The saved verification output and README.

The original verifier was replayed without importing or editing it. Its output
is byte-identical to `verification.json`; the fresh copy is `replay.json`.

## 1. Physical topology and gauge/boundary scope

The note correctly works in the SU(2)-Higgs model with hypercharge omitted,
rather than silently claiming every statement for the full physical model.
The based sphere-at-infinity family and smash-product construction match
[Klinkhamer–Rupp, §§II–IV.A, equations (4.3)–(4.5)](https://arxiv.org/pdf/hep-th/0304167).
The cited map covers the target sphere once, up to orientation. The same
source explicitly distinguishes topological minimax motivation from existence
of an actual solution, and its Fig. 11 places the standard unwrapped
sphaleron at half-integer Chern–Simons coordinate.

[Manton 2019, §4, equations (4)–(5)](https://arxiv.org/pdf/1903.11573)
supports the based mapping-space homotopy description. The note's safer
algebraic hypothesis is a specified boundary homomorphism w with w(g)=1;
it does not need an assertion that every possible configuration-space model
has fundamental group Z. From w(g)=1, infinite order of g follows immediately.

The alternative quotient-by-gauge-group argument is correctly conditional on
a contractible ambient field space and a principal-bundle model. Quotienting
by the identity component versus the full based gauge group is not conflated.
A connecting path before endpoint identification is not silently used as a
based loop. These hypotheses must survive any condensed presentation.

**Scope caution, already respected:** the integer is a loop/boundary winding
or an asymptotic vacuum difference. It is not a blanket quantization statement
for all time-dependent Minkowski gauge fields. Klinkhamer–Rupp §VI.A itself
warns about this distinction. The note makes no such broader claim.

## 2. Profinite completion and supernatural arithmetic

The chain H≅Z -> Z-hat ≅ Gal(F2bar/F2) is correctly typed. The first map is
injective with dense image and is not surjective; the second is a topological
group isomorphism. A supernatural number is then an invariant of the resulting
procyclic group, not a substitute for the group or the original path.

The proof uses every finite quotient Z/nZ. Restricted towers carry their own
supernatural least common multiple; infinitude alone cannot supply all prime
powers. The finite-field realization uses coherent restriction maps and a
common extension of degree lcm(n,m), so choosing normal bases later does not
create a false compatibility requirement.

The inverse-limit and finite-field claims are self-contained. The supporting
normal-basis statement is correctly cited as theorem 5.18 in
[Milne, Fields and Galois Theory](https://www.jmilne.org/math/CourseNotes/FT.pdf).
The stated finite-subfield classification by supernatural degree follows
because finite subdegrees are divisor-closed and lcm-closed. This classifies
subextensions inside a chosen algebraic closure; it says nothing about a
physical field interpretation.

## 3. Frobenius, multiplication, and repeated factors

The two actions are correctly separated:

- primitive-element multiplication has order 2^n−1 and one nonzero state cycle;
- Frobenius has order n and, in a normal basis, is conjugate to an n-coordinate
  cyclic shift, whose minimal polynomial is x^n−1.

Repeated factors in the latter polynomial when n is even are compatible with
separability of the field extension. These are different polynomials with
different roles. The general register formula
lcm_i(r_i 2^ceil(log2 e_i)) is also correct: after the odd root order has been
paid, characteristic-two powering determines the primary thickness needed.
Its distinction between operator order and individual state period is valid.

An independent exhaustive control below checks this formula for every monic
constant-one binary polynomial of degrees 1 through 9, including 170 with
repeated factors. This independently tests a formula not exhaustively checked
by the original field-only verifier.

## 4. Connectedness, midpoint, and character hostiles

All three obstructions are sound:

1. A continuous image of a connected interval or circle in a profinite group
   is a singleton, as seen in every finite discrete quotient.
2. No additive homomorphism Z[1/2] -> Z-hat can take 1 to 1, since reduction
   modulo two would require 0=1. This obstructs a specified additive extension,
   not every arbitrary encoding of the midpoint.
3. Every continuous character Z-hat -> U(1) has finite image. An open subgroup
   inside a sufficiently small identity arc lies in its kernel. Irrational
   characters of discrete Z therefore cannot extend continuously.

For nonzero k, kZ-hat is isomorphic abstractly to Z-hat but has index |k| in
the ambient group. Thus the same supernatural order does forget the chosen
inclusion index. The note states this loss correctly.

## 5. Cover categories and UHF scope

The finite-cover comparison is valid after choosing a circle basepoint,
orientation/generator, and a geometric point over F2. Both resulting fibre
functors give finite sets with permutation monodromy. A connected n-cover
corresponds to Spec(F_(2^n)). Arithmetic versus geometric Frobenius is a
choice of generator convention, which the note fixes explicitly.

[Stacks, tag 03QR](https://stacks.math.columbia.edu/tag/03QR)
gives the required equivalence between étale schemes over a field and
continuous Galois sets; restricting to finite objects yields the claimed
finite category. No equivalence of full physical configuration spaces,
underlying topological spaces, or dynamics follows. The note does not claim one.

The factorial matrix embeddings and the universal UHF invariant agree with
[Ghasemi–Kubiś, §8.1](https://arxiv.org/pdf/1903.10392).
The scaled ordered K0 group is (Q,Q_+,1). The note correctly treats agreement
of supernatural invariants across categories as a dictionary of isomorphism
classes, not an isomorphism of a C*-algebra with a profinite group.

## 6. Independent controls and reproduction

`independent_controls.py` imports none of the original verifier. It performs
2,256 always-active exact checks:

- all 511 monic constant-one binary polynomials of degrees 1–9: factor-order
  formula against direct quotient-ring iteration, plus identity on every
  basis vector at the predicted operator order;
- cyclic shifts of lengths 1–12: exact operator order and every fixed-point
  count for iterates 1–n;
- all 36 ordered degree pairs 1–6: equivariant maps between transitive cyclic
  sets compared with roots of an irreducible degree-m polynomial in F_(2^n).
  The common morphism count is m when m divides n and zero otherwise. This
  checks the direction of the scheme/cover correspondence;
- midpoint equations 2a=1 in every modulus 1–64, with solutions exactly at
  odd moduli;
- an independent degree-one method: target (0,−1,0,0) has the unique preimage
  mu=theta=pi/2, phi=0 modulo 2pi, since
  sin(mu) sin(theta) cos(phi)=1 forces all three factors to equal one.
  The exact local orientation determinant at that point is +1.

Run from the workspace root:

```
python -B session9-sphaleron-review/independent_controls.py
python -B -O session9-sphaleron-review/independent_controls.py
```

The normal and optimized JSON outputs are byte-identical. SymPy emits a
non-fatal deprecation warning about its internal modular-integer ordering
when factoring; this does not affect the arithmetic, stdout, or replay result.
Finite controls support the scoped algebra and conventions. They are not
used to infer an infinite cofinality theorem or a physical isomorphism.

**Conclusion:** the abstract completion/Galois/finite-monodromy bridge may be
published with the existing qualifications intact. The literal path-versus-
supernatural-number identification remains ill-typed; a new physical
intertwining map remains unspecified.
