---
id: THM-2104
title: "Every constant 3-adic quotient layer has a universal sixth-clock mixed escape"
status: >
  PROVED. Let p,q be independent integer characters and write the guard as p
  and the terminals as c_i=a_i p+n_i q. If all nonzero quotient coefficients
  n_i have the same 3-adic valuation v, the torus point determined by
  p.X=1/2 and q.X=1/(6*3^v) is mixed-safe: every terminal distance is 1/6 or
  1/3. This holds for any number of terminals and arbitrary guard-direction
  lifts a_i. It strictly generalizes the dyadic quotient towers behind all
  six THM-2103 no-go rows. Therefore every mixed cover must use at least two
  3-adic quotient layers in every such integral two-character presentation.
  This closes that family, not arbitrary rank eight or LRC(14).
source: codex-2026-07-22-LRC-sixth-clock-3adic-layer
depends_on:
  - THM-2103
related:
  - THM-1065
  - THM-2069
  - THM-2073
  - THM-2095
  - THM-2099
---

# THM-2104 -- the sixth-clock constant-valuation escape

Let `p,q:T^2->T` be independent integer characters, so the torus map

```text
X |->(p.X,q.X)                                         (1)
```

is surjective. For any finite index set, let

```text
g=p,
c_i=a_i p+n_i q,              a_i in Z, n_i in Z\{0}. (2)
```

Assume there is a common integer `v>=0` such that

```text
nu_3(|n_i|)=v                         for every i.      (3)
```

Then there is an `X in T^2` satisfying

```text
||g.X||=1/2,
||c_i.X|| in {1/6,1/3}                for every i.     (4)
```

In particular `X` is strictly safe for a radius-`1/7` guard and every
radius-`1/14` terminal.

## Proof

Put `b=3^v`. By surjectivity in (1), choose `X` with

```text
p.X=1/2,                       q.X=1/(6b).             (5)
```

Write `m_i=n_i/b`. Hypothesis (3) says that `m_i` is an integer not divisible
by three. Equation (2) gives

```text
c_i.X=a_i/2+m_i/6                         in T.        (6)
```

Modulo six, every integer not divisible by three belongs to

```text
{1,2,4,5}.                                             (7)
```

If `a_i` is even, (6) is one of the four residues in (7), divided by six.
If `a_i` is odd, adding `1/2=3/6` permutes the same set:

```text
1<->4,                       2<->5 mod 6.              (8)
```

The four circle distances are exactly `1/6` and `1/3`. This proves (4).
Since

```text
1/2>1/7,                     1/6>1/14,
```

the point is a strict mixed escape. QED.

## Dyadic and multi-orbit corollaries

The earlier dyadic tower has

```text
n_i=b_0 2^(e_i).                                      (9)
```

All its coefficients have `3`-adic valuation `nu_3(|b_0|)`, so (9) is an
immediate special case. More generally, (3) permits arbitrary odd factors,
arbitrary powers of two, signs, repetitions, and any number of distinct
dyadic orbits. The sixth-clock point depends only on the common `3`-adic
layer, not on those finer multiplicative data.

Thus the six THM-2103 exceptions were not six sporadic geometries and not
merely one dyadic orbit. They lie in a full valuation layer on which the
threshold-labelled residue evaluation is constant-safe.

## Coefficient-plane consequence

In an integral presentation (2), a putative mixed cover must have

```text
#{nu_3(|n_i|):i}=at least 2.                            (10)
```

This is a basis-aware statement. Replacing `q` by another transverse
character changes the quotient coordinates `n_i`; one may use (10) only in
a named integral two-character presentation. The invariant formulation is
that the terminal images in the rank-one quotient lattice by `Zp` cannot all
lie in one `3`-adic valuation shell.

When `Zp` is not saturated in the ambient character lattice, the quotient
has torsion. The theorem still applies exactly as stated whenever the
decompositions (2) exist; passing to a saturated quotient requires retaining
the torsion class as a sidecar rather than silently identifying it with an
integer coordinate.

## Frontier effect

The carrier ladder sharpens to

```text
pair weights
 -> signed affine rank
 -> sixth-clock residue
 -> 3-adic quotient valuation profile.                 (11)
```

THM-2103 found the residue empirically in a finite cube. The present theorem
makes it all-height and shows that the first live quotient obstruction is
not “non-dyadic”; it is **valuation diversity at the prime three**. The next
lawful split is by the ordered valuation layers and their unit residues,
with THM-2095's p-adic deck warning that layer multiplicities, not only prime
support, must be retained.

## Assumption challenge and Tournament Analysis

The challenged assumption is that powers of two cause the sixth-clock
escape. They only keep a coefficient away from zero modulo three. The actual
invariant is constant `3`-adic quotient valuation; arbitrary unit factors are
harmless.

Candidate tournament vertices were terminals, dyadic exponents, valuation
layers, unit residues modulo six, clock states, and proof obligations.
Orienting terminals by `nu_3(n_i)` collapses to a complete tie in this theorem;
any tie-break tournament has artificial scores, cycles/SCCs, edge flips, and
Hamiltonian paths. The faithful carrier is the quotient valuation shell
decorated by the mod-six unit residue and first-coordinate parity. No
computational theorem is required.

This closes only the constant-layer branch (3). It does not bound the number
or separation of valuation layers in a cover, treat quotient torsion without
its sidecar, discharge the vertical branch, or prove LRC(14).
