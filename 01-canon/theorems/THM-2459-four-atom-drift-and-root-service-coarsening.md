---
id: THM-2459
title: "Four-atom drift and root-service coarsening"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT. A nonzero aggregate Hilbert-valued linear observable
  and one positive directed atom edge can be retained simultaneously
  by a Boolean union of at most four atoms. The exact identity gives
  norm loss at most 2r-1, where r is the number of atoms outside the
  selected atom and edge. For the 128 THM-2452 atoms, a maximum
  THM-2457 edge yields one lawful union with endpoint drift at least
  D_0/63001 and root-service mass at least M_0/16384. Four atoms and
  the factor 2r-1 are sharp. Nonzero selected-atom drift without
  nonzero aggregate drift is insufficient, already for two
  nonnegative atoms. An arbitrary rational root-constant observer can
  repair non-pointwise cancellation, but is not automatically an LRC
  factor or canonical relation current. No scalar row is excluded.
source: codex-2026-07-26-drift-service-coarsening
depends_on:
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
related:
  - THM-2280-centered-polynomial-grid-avoidance-and-bounded-generic-keller-fibre
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2448-right-endpoint-cospan-transition-atlas
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2456-two-root-replica-uniform-offset-boundary
script: 04-computation/lrc14_four_atom_drift_service_thm2459.py
output: 05-knowledge/results/lrc14_four_atom_drift_service_thm2459.out
script_sha256: 76e24ced413a53e20e349eac4e4cf83f89738f0589c6c128971fc497a16de5db
output_sha256: bfe22300e13fe3dcb7e51a123f947d8e7e44ed21a1fdad099b3665d9fd4fd813
hash_basis: working-tree bytes (LF)
---

# THM-2459 -- four atoms retain drift and root service

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

THM-2457 proves that a positive directed edge in the complete-atom
co-support graph gives root service after Boolean coarsening.  It also
warns that an arbitrary such coarsening can cancel the signed
THM-2452 endpoint drift.  The missing fact is that THM-2452 has more
than one nonzero atom: its **aggregate** table has nonzero drift.

One additional atom then suffices.  The mechanism is the exact
Hilbert-space identity

```text
q_total
 =sum_(j outside S)(q_S+q_j)
  -(#(outside S)-1)q_S.                              (1)
```

It is a degree-one Boolean avoidance argument.  Unlike a generic
polynomial-grid choice, every coefficient stays zero or one, so the
selected object remains an indicator union.

## 1. Abstract Hilbert coarsening lemma

Let `Omega` be a finite set of `N` atoms.  Let `H` be a real or complex
Hilbert space and let

```text
q_omega in H,                    omega in Omega,

q_total=sum_omega q_omega!=0.                       (2)
```

Let

```text
M=(M_(omega,nu))_(omega,nu in Omega)
```

be a nonnegative directed matrix.  A loop is allowed.  Fix:

- a selected atom `kappa`; and
- a positive directed edge `u->v`, so `M_(u,v)>0`.

Write

```text
E={u,v},             with repetitions removed,

S={kappa} union E,

J=Omega\S,

r=|J|.                                                (3)
```

For `I subset Omega`, put

```text
q_I=sum_(omega in I)q_omega,

M(I)=sum_(omega,nu in I)M_(omega,nu).                (4)
```

If `r=0`, take `I=Omega=S`; then already `|I|<=3`.  Suppose
henceforth that `r>=1`.  Define

```text
z_0=q_S,

z_j=q_(S union {j})=q_S+q_j,          j in J.        (5)
```

Then

```text
q_total=sum_(j in J)z_j-(r-1)z_0.                   (6)
```

The triangle inequality gives

```text
max(||z_0||,max_(j in J)||z_j||)
 >=||q_total||/(2r-1).                              (7)
```

Choose `I=S` or `I=S union {j}` attaining the left side.  It satisfies

```text
kappa in I,

|I|<=4,

M(I)>=M_(u,v)>0,

||q_I||>=||q_total||/(2r-1)>0.                      (8)
```

Thus one Boolean union retains the selected atom, the supplied root
edge, and a nonzero aggregate linear observable.

The statement is valid for every Hilbert-valued **linear** observable.
In particular, if `T_omega` are nonnegative matched endpoint tables,
`Q` is a linear orthogonal projection, and

```text
q_omega=Q T_omega,
```

then

```text
D(T_I)=||Q T_I||^2
 >=D(T_total)/(2r-1)^2.                             (9)
```

Nonnegativity is not needed for (6)--(9).  It is needed to say that
the union cannot erase a positive selected whole table merely by
adding other atoms.

### Proof of (6)

There are `r` terms in the first sum, so

```text
sum_(j in J)z_j
 =r q_S+sum_(j in J)q_j.
```

Subtracting `(r-1)q_S` leaves

```text
q_S+sum_(j in J)q_j=q_total.
```

Equation (7) follows because the right side of (6) is a signed sum of
`r+(r-1)=2r-1` vectors, each bounded by the displayed maximum. QED.

## 2. The 128-atom THM-2452/2457 consequence

Use the 128 complete atoms of THM-2452.  Write

```text
T_0=sum_omega T_omega,

D_0=||Q T_0||^2>0.                                  (10)
```

Choose the THM-2452 atom `kappa` with

```text
D(T_kappa)>=D_0/16384.                              (11)
```

Assume the additional THM-2457 service hypotheses: one common
oriented physical `C_13` chart, one disjoint rational semantic pair,
and positive initial co-support mass

```text
M_0=sum_(omega,nu)M_(omega,nu)>0.                   (12)
```

Choose a maximum directed entry.  Then

```text
M_(u,v)>=M_0/128^2=M_0/16384.                       (13)
```

If `(u,v)=(kappa,kappa)`, the single atom already has both properties,
and (11), (13) apply.

Otherwise `|S|>=2`, so

```text
r=128-|S|<=126,

2r-1<=251.                                          (14)
```

Section 1 gives a set `I` of at most four complete atoms such that

```text
D(T_I)>=D_0/251^2=D_0/63001,

M(I)>=M_0/16384.                                    (15)
```

Equation (15) is uniform across the loop, incident-edge, and
disjoint-edge cases: in the exceptional loop at `kappa`, the stronger
`D_0/16384` bound from (11) holds.

THM-2457 now gives, for its normalized common-root coefficients,

```text
sum_(a!=0)|J_I(a)|^2
 >=M_0^2/92001420705792,

max_(a!=0)|J_I(a)|
 >=M_0/33226752,                                    (16)
```

and every one of the twelve nonzero root colours survives.  The
denominators are

```text
92001420705792=16384^2*342732,

33226752=16384*2028.                                (17)
```

Independently, THM-2365 can be applied to the retained drift in (15)
to reselect an eligible target/deep colour and then a fresh exact
address.  Equations (15)--(16) put both nonzero structures in one
Boolean packet.  They do not say that the two later Fourier
selections choose the same phase or address.

## 3. Why the union is a lawful complete-atom filter

For the THM-2452 partition, define

```text
P_I=sum_(omega in I)P_omega.                        (18)
```

The complete atoms are pairwise orthogonal indicators.  Hence

```text
P_I^2=P_I,

T_I=sum_(omega in I)T_omega.                        (19)
```

The same `P_I` is copied to the present and bare endpoint legs.  Thus
the union in Section 2 is a lawful Boolean formula in the already
installed seven physical truth bits, not an integer-weighted formal
combination or an external observer.  Every factor fixed outside the
complete local atom -- including a terminal word when it was fixed
independently before the atom split -- remains fixed.

The union generally no longer names one literal complete atom.  It
also cannot infer a semantic pure/fork word from local truth bits;
THM-2457's strict-profile hostile still applies to that inference.

## 4. Four atoms are necessary

The cardinality bound in (8) is sharp even for nonnegative two-entry
tables.  Let

```text
Q(x,y)=((x-y)/2,(y-x)/2)                            (20)
```

be projection onto the mean-zero line.  Take four atoms

```text
kappa,u,v,j
```

and the tables

```text
T_kappa=(1/4,0),

T_u    =(0,1/4),

T_v    =(1/8,1/8),

T_j    =(1/4,0).                                    (21)
```

They are nonnegative and their coordinatewise sum is at most one.
Let the only positive co-support entry be

```text
M_(u,v)=1.                                          (22)
```

Then

```text
Q(T_kappa+T_u+T_v)=0,

Q(T_kappa+T_u+T_v+T_j)!=0.                          (23)
```

Every service union which retains `kappa` must contain `u,v`.
With at most three atoms it is therefore exactly
`{kappa,u,v}` and has zero drift.  Four atoms are necessary.

The edge matrix (22) has a common-root realization: on one rational
base interval put an `A` singleton at root zero in atom `u`, an `F`
singleton at root one in atom `v`, and no other semantic packet
roots.  The packets are disjoint and only `u->v` fires.

## 5. The factor 2r-1 is sharp

Fix `r>=1` and one nonzero vector `e`.  Arrange

```text
q_S=-e,

q_j=2e                    for every j in J.          (24)
```

Then every candidate in (5) has norm `||e||`, whereas

```text
q_total=(2r-1)e.                                    (25)
```

So (7) is exact.  These vectors can be realized by nonnegative
two-entry tables under (20): add a sufficiently large constant vector
to each table, then scale all tables by one common positive rational
so their total remains at most one.  Constant vectors lie in the
kernel of `Q` and do not change (24).

This sharpness is for the fixed selected atom and fixed service edge.
Additional graph structure can permit a better choice of edge.

## 6. Aggregate drift is the load-bearing sidecar

If (2) is removed, even existence fails with two atoms.  Take

```text
T_kappa=(1/2,0),

T_v=(0,1/2),

M_(kappa,v)=1
```

and all other co-support entries zero.  Then

```text
Q T_kappa!=0,

Q(T_kappa+T_v)=0.                                   (26)
```

The only service union retaining `kappa` contains both atoms, so it
has zero drift.  This is minimal: with one atom, a positive loop and
nonzero selected drift already coexist.

The hostile can be made pointwise in the common base coordinate:
make the two projected table densities opposite at every `y`.  Then
every root-constant gate multiplies the zero sum and still has zero
drift.  Thus

```text
one nonzero selected atom
```

cannot replace

```text
one nonzero aggregate observable.                  (27)
```

THM-2452 supplies exactly the stronger datum in (10).

## 7. Optional post-hoc observer selection

There is a useful but noncanonical complement.  Let `Y` have a finite
rational step partition into positive-measure cells modulo null sets,
let

```text
q:Y->H
```

be a step-function projected-table density, and let

```text
s:Y->R_(>=0)
```

be a service density with positive integral.  An arbitrary rational
root-constant Boolean gate `G(y)` can make

```text
integral Gq!=0,

integral Gs>0                                      (28)
```

if and only if `q` is not identically zero almost everywhere.
Moreover, `G` may be
chosen as the union of at most two cells of the common step
partition.

Indeed, choose a cell `C` on which `s>0`.  If `q` is nonzero on `C`,
take `G=1_C`.  Otherwise choose any cell `D` on which `q` is nonzero
and take

```text
G=1_(C union D).
```

The two contributions cannot cancel because the contribution from
`C` is zero.  Necessity is immediate.

This gate is lawful for the abstract common-filter algebra of
THM-2401.  It is not automatically lawful for the LRC current:
the chosen rational cells need not lie in the Boolean algebra
generated by the installed combs, owner word, and delayed terminal
factor.  Introducing them as new physical factors changes Fourier
support and may erase the canonical owner/address interpretation.
If a prescribed lawful filter algebra contains the two selected
cells, the argument becomes physical; without that separation
sidecar, it is only an observer/coefficient grouping.

The pointwise hostile in Section 6 proves that even an arbitrary
root-constant observer cannot repair every aggregate-zero packet.

## 8. Exact companion

Run

```text
python 04-computation/lrc14_four_atom_drift_service_thm2459.py
python -O 04-computation/lrc14_four_atom_drift_service_thm2459.py
```

The dependency-free `Fraction` companion:

- checks (6) and the sharp `2r-1` ratio for every `1<=r<=127`;
- verifies the exact `128`-atom constants and both root-energy
  denominators;
- checks the four-atom nonnegative-table hostile and its sole edge;
- checks the minimal two-atom aggregate-zero hostile;
- verifies rational nonnegative realizations of the sharp vector
  pattern; and
- exhausts all Boolean unions in both finite hostiles.

Both transcripts must match

```text
05-knowledge/results/lrc14_four_atom_drift_service_thm2459.out
```

byte-for-byte after LF normalization.  Every truth-bearing executable
check uses explicit `require`, including under optimized Python.

## 9. Scope

The theorem removes one precise post-THM-2457 debt:

```text
nonzero aggregate THM-2452 drift
  + one positive THM-2457 atom edge
  -> one at-most-four-atom packet carrying both.     (29)
```

It does not construct THM-2457's common oriented root chart or
disjoint semantic packets, prove their positive same-base
co-support, align the root-service phase with the independently
reselected target/deep phase, or turn the resulting endpoint
coefficients into one canonical `91`-unit relation current.

THM-2456's averaged uniform-offset hostile is also not excluded by
the abstract vector identity; the theorem works after a positive
semantic edge has been supplied.  No scalar profile is excluded, the
ledger remains `165`, and LRC(14) remains open.

QED candidate.
