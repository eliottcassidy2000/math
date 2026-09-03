---
title: "JC2 source-normal weight-twelve row-twelve extinction"
status: >
  INTEGRATED AS THM-4380: PROVED FINITE-ROW RELATIVE TO THM-4308/4315
  and VERIFIED-EXACT + INDEPENDENTLY AUDITED by a 314-check primary plus an
  independent theorem/code/branch-algebra review and triple replay. Scratch broadcast commit
  84590342a2 independently corroborates only the strict Phi*U!=0 branch;
  a full second implementation remains pending. Chart or seam
  entry, higher residual weights, all-row lifting, JC(2), and DC(2) remain
  OPEN.
source: source_seam / JC2 continuation session, 2026-09-03
artifacts:
  - 01-canon/theorems/THM-4380-source-normal-weight-twelve-row-twelve-extinction.md
  - 04-computation/jc2_source_normal_weight12_row12_extinction_thm4380.py
  - 05-knowledge/results/jc2_source_normal_weight12_row12_extinction_thm4380.out
script_sha256: 15cd129452c3da033fa59985dda435077fe4526dbdf9d9df6feae32a8cb0ac6a
output_sha256: 5b37b424137dba93fe7e7c4c621cf96d0068a9247a8bc9a3ff962bbf74279429
hash_basis: raw LF bytes
related:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
  - THM-4366-source-normal-u-zero-row-eleven-hierarchy-selected-extinction
  - THM-4376-source-normal-u-zero-row-eleven-depth-hierarchy-completeness-and-bracket-blindness
---

# The full capped source family is absorbed by row twelve

## 1. Inheritance and live board

The nearest proved mechanism was THM-4308's exact source-normal row-eight
chart, with THM-4315 supplying the full row-nine bracket hypersurface. The
canonical hostile was THM-4376: projected-depth hierarchies can remain
complete on an already bracket-dead counterfactual carrier. It therefore
would have been unsafe to infer later bracket survival from depth survival.
The least-used sidecar was the unexamined `Phi=0` stratum of the full source
equation, rather than another special coordinate slice such as `U=0` or
`beta=0`.

The live board was

```text
full E_9 hypersurface | Phi boundary | row-ten joint-depth mismatch
ratio quotient (r,Y) | sign sheet | row-eleven depth | row-twelve Student cut.
```

## 2. The signal that changed the route

The first exact pull was stronger than expected: on the whole row-nine
hypersurface, both projected-depth systems are compatible and their joint
terminal rank is only three. Depth adds no source equation and leaves an
affine-seven fibre. That made the right next operation a genuine stratum
split, not another hierarchy consumer.

The resulting deterministic tree is

```text
E_9=0
 |
 +-- Phi=0:    two reduced row-ten points -- row-eleven fixed residual --> empty
 |
 +-- Phi!=0:   row-ten curve D=0
                  |
                  +-- row eleven: K(r)=0, fourteen reduced points, depth automatic
                                      |
                                      +-- row twelve: gcd(K,J)=1 --> empty.
```

This closes every source point in the fixed residual-weight-at-most-twelve
universe. It is not just a continuation of the `U=0` calculation: the
fourteen nonzero-`Phi` points all satisfy `U!=0`, proved by the additional
univariate certificate `gcd(K,Ucar)=1`.

## 3. Why the ratio quotient worked

On `Phi!=0`, setting

```text
r=eta/Phi,                  Y=Phi^2
```

turns the row-ten joint-depth equation into a linear graph

```text
A(r)Y=B(r).
```

This quotient removes the simultaneous sign symmetry without erasing its
fibre. The essential sidecar is the equation `Phi^2=Y`: after the degree-seven
ratio eliminant is found, every nonzero `Y` has two sign lifts. Auditing
`gcd(K,A*B)=1` proves both that the ratio graph has no exceptional fibre
and that every sign lift is present. Squarefreeness of `K`, together with
`Y!=0`, upgrades the count from length fourteen to fourteen reduced points.

This is a reusable warning about quotients. The quotient carries seven
addresses; the original source carrier has fourteen points. Dropping the
sign-sheet sidecar would preserve the extinction test but state the wrong
geometry and could corrupt later rank or multiplicity claims.

## 4. Depth and bracket separate cleanly

At row eleven the two projected-depth universes have left nullities 25 and
20. All 45 residuals vanish on the fourteen-point carrier, while the next
bracket does not. The primary explicitly audits that every denominator used
in the Groebner reductions is a monomial in `Phi`; no factor of `A`, `B`,
`K'`, or a source response is silently inverted.

The final cut is especially clean. The row-twelve residual pulls back to a
quadratic `J(r)`, and a literal Bezout identity modulo 29 proves
`gcd_Q(K,J)=1`. Because both leading coefficients survive modulo 29, the
certificate has no bad-prime degree loss. This gives a precise instance of
the phenomenon highlighted by THM-4376: depth completeness and bracket
survival are genuinely different predicates.

## 5. Stochastic-process analogy, used carefully

The staged calculation resembles an absorbing process: the row index is a
filtration, the `Phi` strata are deterministic states, and the first nonzero
residual is an absorbing boundary. The analogy helped organize nonvacuity
controls—one must exhibit a live state immediately before each claimed
absorption—but no probability or averaging enters the proof. The source
points are algebraic, the transitions are exact eliminations, and the
finite-field control is a hostile trajectory rather than a random sample.

The useful transfer from stochastic thinking is therefore procedural:
separate entry, survival, and killing events, and never condition on a state
whose probability—or here, algebraic nonemptiness—has not been proved. The
mathematics itself remains deterministic commutative algebra.

## 6. Audit seams that mattered

Several small checks prevented a plausible but incomplete proof:

- coefficient comparison in degrees `0,...,n` is exhaustive only after
  verifying that the capped bracket difference has degree at most `n`;
- the `Phi=0` formula for `alpha` divides by `eta`, so the nonzero constant of
  its quadratic carrier must be used explicitly;
- `A=0`, `B=0`, roots at infinity, and a forgotten sign sheet are distinct
  possible losses in the ratio parametrization;
- Groebner reduction of rational terminal formulas is valid only after
  recording every localization;
- modular coprimality needs retained degrees at the chosen prime;
- a strict `Phi*U!=0` corroboration does not independently audit the
  `Phi=0` branch or the inherited full source construction.

The primary now checks each item. Scratch broadcast commit `84590342a2`
independently found the same strict-branch septimic, fourteen-point carrier,
45 vanishing depth residuals, and mod-29 terminal cut. Since no checked-in
artifact accompanied it, it remains corroboration rather than a dependency.

## 7. The next genuinely open seam

There is no value in pushing this same capped source carrier to row thirteen:
it is already empty. The next sharp problems are instead:

1. prove or refute entry of the relevant global/seam geometry into the fixed
   THM-4308 source-normal chart;
2. determine whether residual weights above twelve create new source channels
   before the row-twelve cut;
3. find a conceptual module-theoretic reason for the automatic joint depth at
   rows nine and eleven, rather than recomputing its left-null equations;
4. build a clean-room full-family referee, especially reconstructing the
   `Phi=0` branch independently.

Until one of the first two entry problems is solved, THM-4380 is an exact
finite-chart extinction result and not progress from a live Keller pair to
`JC(2)`. That stopping boundary is part of the result, not an omitted final
step.
