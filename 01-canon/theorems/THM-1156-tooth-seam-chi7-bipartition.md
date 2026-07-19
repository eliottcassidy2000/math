---
id: THM-1156
title: THE TOOTH-SEAM / chi7 BIPARTITION -- exact radius-1/14 abutments are projectively complementary mod 14 and flip a canonical quadratic sign
status: PROVED exact seam identity, abutment criterion, chi7 bipartition, Fano pair obstruction, and directed strict-seam quantum. The open-cover consequence is qualitative: a zero-overlap two-tooth seam needs third support. No uniform F7 overlap constant or LRC14 closure is claimed
source: codex-2026-07-18-S75
depends_on: [THM-770, THM-965]
related: [THM-856, THM-863, THM-1153, THM-1155, THM-1165, HYP-7678]
script: 04-computation/lrc14_tooth_seam_chi7_bipartition_codex_20260718.py
output: 05-knowledge/results/lrc14_tooth_seam_chi7_bipartition_codex_20260718.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCToothSeamChi7.lean
---

# THM-1156 -- the tooth-seam / chi7 bipartition

The density ledgers in THM-1153 and THM-1155 vanish at seven combs.  This
theorem keeps the information they discard: which labelled tooth endpoint
meets which other endpoint.  It produces a genuine `chi_7` law, but it also
shows precisely why that law alone does not cross the wall.

## 1. Exact directed seam law

On a real lift, write the open radius-`1/14` tooth of speed `a` and index `m`
as

```text
T_(a,m)=((14m-1)/(14a),(14m+1)/(14a)),
L_(a,m)=(14m-1)/(14a),
R_(a,m)=(14m+1)/(14a).
```

For positive integers `a,b`, the signed gap from the right endpoint of the
`a`-tooth to the left endpoint of the `b`-tooth is

```text
Delta(a,b;m,n)=L_(b,n)-R_(a,m)
 = [14(an-bm)-(a+b)]/(14ab).                            (1)
```

Thus `Delta>0` is a gap, `Delta<0` is directed penetration, and `Delta=0` is
exact abutment.  If `g=gcd(a,b)`, Bezout gives

```text
{an-bm:m,n in Z}=gZ.                                   (2)
```

Consequently an exact oppositely oriented seam exists if and only if

```text
14g divides a+b.                                       (3)
```

This is the general determinant coordinate behind THM-1147's specialized
adjacent-label law.  It retains both endpoint owners and both tooth indices.

## 2. The canonical `chi_7` colour

For a positive integer `v`, strip its entire 7-adic factor and define

```text
epsilon(v)=Legendre_7(v/7^nu7(v)) in {+1,-1}.           (4)
```

> **Theorem A (seam bipartition).** If speeds `a,b` admit an exact seam, then
> `epsilon(a)=-epsilon(b)`.  Hence the graph on speeds whose edges are exact
> abutment-compatible pairs is bipartite.  In particular it has no odd cycle
> and no triangle.

**Proof.** Put `a=gA`, `b=gB`, with `(A,B)=1`.  Condition (3) says
`14|(A+B)`.  Therefore `A,B` are both odd.  Neither is divisible by `7`,
since otherwise both would be, contradicting coprimality.  It follows that
`nu7(a)=nu7(b)=nu7(g)`.  Write `g=7^t h` with `7` not dividing `h`.  Modulo
seven, `B=-A`, and the Legendre character is multiplicative, so

```text
epsilon(b)=chi7(hB)=chi7(h)chi7(-A)
          =chi7(-1)chi7(hA)=-epsilon(a),
```

because `chi7(-1)=-1`. ∎

This explains the previously decorative `chi_7` hint: the natural binary
switch is not an orientation of runners by size, but the quadratic colour of
the 7-primitive speed.

## 3. Fano consequence

Label any seven speeds by the seven points of the Fano plane.  If `p` of the
speeds have colour `+1`, the number of same-colour pairs is

```text
C(p,2)+C(7-p,2)>=9,                                    (5)
```

with equality for a `3+4` split.  Every such pair is forbidden from exact
abutment by Theorem A.  Since every Fano line contains three points and every
pair lies on a unique line, every Fano line contains at least one forbidden
pair (and the seven lines collectively carry at least nine forbidden pairs,
counted once each).

This is an exact Fano obstruction, not an overlap lower bound.  The cover is
free to use only cross-colour transitions, repeat owners, form a path rather
than an odd cycle, or replace a pair transition by a multiway event.

## 4. The sharp one-sided seam quantum

There is also an exact metric residue.  Define

```text
r_plus(a,b)=a+b-14g floor((a+b-1)/(14g)),
1<=r_plus(a,b)<=14g.                                   (6)
```

Among all strictly negative values of (1), the least positive directed
penetration is

```text
min_{Delta<0} (-Delta)=r_plus(a,b)/(14ab).              (7)
```

Indeed strict negativity is `14gz<a+b`; the largest admissible integer `z`
is the floor in (6), and (2) makes that value attainable.  When (3) holds,
zero is attainable and the next strictly overlapping seam has
`r_plus=14g`.

Equation (7) must not be promoted to an unconditional pair-intersection
quantum.  If the `b`-centre lies to the right of the `a`-centre and
`d=an-bm`, the penetration equals the actual tooth intersection length only
in the crossing regime

```text
|a-b|<14d<a+b.                                         (8)
```

When `14d<=|a-b|`, one tooth contains the other and the intersection is the
whole narrower tooth.  The seam residue and the pair-overlap measure are
different observables.

## 5. Open-cover consequence: zero seams become triple events

The strict/open convention is decisive.  If two open teeth exactly abut at
`x`, then `x` belongs to neither tooth.  Therefore, inside any closed
protected needle covered by open danger combs, an exact `a`-to-`b` seam must
have a third danger comb containing `x`.  Since that third tooth is open, it
also covers a neighbourhood of `x`.

Thus every interior handoff has one of two forms:

1. a strictly penetrating pair seam, with the quantum (7) when it is a
   genuine crossing; or
2. a third-supported hyperedge at the putative zero seam.

This is the precise local private/triple alternative suggested by THM-1149.
It upgrades the wall carrier from a pair tournament to an endpoint-owner
hypergraph, but supplies no uniform amount of overlap: `1/(14ab)` can tend to
zero, containment changes the metric, and triple-supported seams can bypass a
simple owner cycle.

## 6. Tournament and carrier audit

The pairwise observable is the integer seam numerator

```text
14(an-bm)-(a+b),
```

while the switch/gauge is `epsilon`.  The zero-seam graph is bipartite, not a
tournament.  Orienting it by speed order produces a transitive quotient with
singleton SCCs, no directed triangles, and one Hamiltonian path, but loses
the endpoint indices, metric residue, crossing/containment flag, and third
support.  Fano lines are useful test hyperedges; runners alone are not the
faithful vertices.

The smallest faithful current carrier is

```text
(labelled tooth endpoints, owner speeds, chi7 colours,
 seam numerator, crossing/containment flag, third-support incidence).
```

The exact replay audits (1)--(7), bipartiteness, triangle absence, and the
Fano `>=9` count on finite banks.  The Lean consumer formalizes the arithmetic
seam identity and the finite residue/colour core exposed by the proof.  These
checks support the theorem; the all-range proof is the Bezout/Legendre
argument above.

## Honest frontier

THM-1156 identifies a missing structural coordinate and proves a qualitative
pair-or-triple debt at every zero seam.  It does **not** prove HYP-7678's
quantitative `(F7)` inequality: a tree is bipartite, owner labels can repeat,
and the seam quantum is scale-dependent.  A closing theorem must aggregate
these local debts across the protected needle while controlling reuse of the
same third-support tooth.
