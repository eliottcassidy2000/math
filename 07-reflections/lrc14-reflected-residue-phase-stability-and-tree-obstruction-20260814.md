# Reflected-residue perturbations: the missing phase and the fixed-tree survivor

**Status.** The residue-blind extension of the THM-3360/3372 high-pair floors
is **REFUTED** by a one-coordinate, one-unit exact witness. A phase-aware
frozen-tree perturbation lemma is **PROVED**; noncanonical repeated-level
packets and a cross-tree obstruction at coherent shift `c=84` are
**FINITE-EXACT + VERIFIED**.
This is not a closure of arbitrary `k=1`, other residue packets, or LRC(14).

## 1. Inheritance and connection contract

[THM-3360](../01-canon/theorems/THM-3360-uniform-reflected-high-pair-floor-by-admissible-affine-tails.md)
proves, and
[THM-3372](../01-canon/theorems/THM-3372-independent-superunit-affine-tail-and-reflected-residue-closure.md)
independently weakens, the matched-residue pair floor: on each of the `649`
upper-median bodies, every canonical reflected high pair

```text
z_e=q_e L-e,  z_f=q_f L-f
```

has enough physical overlap for the five-edge cross-component Hunter tree.
The closest hostile is the body `E=(1,2,3,4,6,12)` at its upper-median cell
`(L,j)=(168,90)`. The corrected near miss is to call `a-e` a small
perturbation in `qL-a`: relative to `qL-e`, the drift displacement is
`h=e-a`, so the cell phase changes by `(e-a)j/L`. This is macroscopic for a
unit residue change at this cell. The least-used sidecar is exactly that
located phase.

The connection is

```text
source: A_z(j)={u in [0,1]: ||z(j+u)/L||<1/14}
target: A_(z+h)(j)
map:    phi_h(u)=(z u-Delta)/(z+h),  Delta=h j-mL
```

on the untruncated real-line comb, for any integer `m`.  It preserves the
ordered comb and gives exact endpoint transport.  Passing only `h` destroys
the cell phase; the required sidecar is the centered residue `Delta mod L`.
Changing residues also destroys the canonical singleton identity
`1/7+e/[7(qL-e)]`; the perturbed proof must use literal singleton masses or an
independent upper bound.

## 2. Exact affine conjugacy

Let `L>0`, `0<=j<L`, `z>0`, and integer `h` satisfy `z+h>0`. Let

```text
B_z(j)={u in R: ||z(j+u)/L||<1/14},
A_z(j)=B_z(j) intersect [0,1].
```

Choose any integer `m`, put `Delta=hj-mL`, and define
`phi_h(u)=(zu-Delta)/(z+h)`.  Then

```text
(z+h)(j+phi_h(u))/L
 = z(j+u)/L+m.
```

Therefore

```text
B_(z+h)(j)=phi_h(B_z(j)).                              (1)
```

This is the exact source-to-target map.  The clipping by `[0,1]` is still
literal and cannot be discarded, but it is finite interval arithmetic.

## 3. A cheap rigorous `L^1` bound

Choose `m` so that `Delta` is a centered residue and put

```text
eta=(|Delta|+|h|)/L.
```

For `0<=u<=1`, the circular phase displacement is at most `eta`.  If the two
radius-`1/14` indicators differ, the old phase lies in the `eta`-neighbourhood
of one of the two danger-arc boundary points.  That periodic boundary
neighbourhood has circle measure at most `min(1,4eta)`.

The old phase traverses an interval of length `alpha=z/L`.  Every complete
turn contributes its circle measure and the final partial turn contributes
at most one more copy.  Hence

```text
sigma(z,h):=mu(A_z(j) triangle A_(z+h)(j))
 <= min(1, 4 eta (1+1/alpha))
 =  min(1, 4 eta (1+L/z)).                             (2)
```

No canonical residue assumption enters `(1)` or `(2)`.

Now freeze a spanning tree `T` and write

```text
M_T(A)=sum_(ij in T) mu(A_i intersect A_j)
       -[sum_i mu(A_i)-6/7].
```

Since

```text
|mu(A_i' intersect A_j')-mu(A_i intersect A_j)|
 <=sigma_i+sigma_j,
mu(A_i')<=mu(A_i)+sigma_i,
```

one obtains the exact perturbation gate

```text
M_T(A') >= M_T(A)-sum_i (deg_T(i)+1)sigma_i.           (3)
```

Equations `(2)`--`(3)` are the cheapest rigorous survivor: any canonical or
already-certified packet whose frozen-tree margin exceeds the displayed
phase-aware budget extends to that noncanonical perturbation.  The criterion
also permits simultaneous changes in several coordinates.

## 4. Minimal counterexample to a residue-blind pair floor

Use `E=(1,2,3,4,6,12)`, `(L,j)=(168,90)`, and the high `3:5` channel on
labels `(12,4)`.  The canonical pair is

```text
z=3L-12=492,          w=5L-4=836.
```

Its exact clauses are

```text
A_z=(5/41,7/41) union (19/41,21/41) union (33/41,35/41),
mu(A_z intersect A_w)=6/209>Dmax/5.
```

Change only the first residue by one:

```text
z'=3L-13=491.
```

Then

```text
A_z'=(0,6/491) union (150/491,174/491)
     union (318/491,342/491) union (486/491,1),
mu(A_z' intersect A_w)=0.                              (4)
```

For `h=-1`, choosing `m=-1` gives `Delta=78`, so the located phase jump is
`78/168=13/28`.  Thus `(4)` is not a continuity paradox: a unit drift change
is a near-half-turn at the chosen cell.  It is minimal in Hamming support and
absolute nonzero integer displacement.  It refutes every extension that keeps
only the level ratio, `|h|`, or the canonical pair floor and drops `Delta`.
It is a pair-certificate counterexample, not an LRC counterexample.

## 5. Strict noncanonical positive control

Take the `649`-bank body

```text
E=(1,2,3,8,9,14),  (L,j)=(7056,3780),
q=(3,3,3,5,5,5).
```

The two repeated-level classes form the two parts of the restricted
cross-`K_(3,3)` certificate. All six drifts are distinct. The exact maximum
cross-tree is

```text
T=((1,4),(1,5),(0,5),(2,3),(0,3))
```

with canonical margin

```text
M=71440713312252560278/527394888495258135905.           (5)
```

Perturb only `z_0=3L-1` to `z_0+28=3L+27`.  This is genuinely noncanonical,
but `28j=15L`, so `Delta=0` and `eta=1/252`.  Exact interval arithmetic gives

```text
sigma_0=67424/16616095,
sigma_0 <= 28223/1333521                              (6)
```

where the right side is `(2)`.  Vertex zero has tree degree two.  Even the
coarse theorem-only bound remains strict:

```text
M-3(28223/1333521)
 =113864784228062404699/1582184665485774407715>0.       (7)
```

The literal perturbed frozen-tree margin is

```text
1585980049263141094/11735389638648206265>0.            (8)
```

This control shows that `(3)` is not merely a vacuous repair of the hostile:
it certifies a distinct-drift, repeated/mixed-level, noncanonical packet.

### A broader all-`649` finite signal

A second exact companion tests the same two equal-level classes
`(3,3,3,5,5,5)` on every one of the `649` upper-median bodies, every
three-versus-three assignment, and the uniform noncanonical residue shifts

```text
a_e=e+c,             c in {-2,-1,1,2},
```

whenever every shifted residue is positive.  All drifts remain distinct.
The literal universe has

```text
c=-2:  1,700 packets,       c=-1:  5,180 packets,
c=+1: 12,980 packets,       c=+2: 12,980 packets,
total: 32,840 packets.
```

Every exact `K_(3,3)` maximum-tree margin is positive.  The four minimum
margins are

```text
c=-2: 704886090693279025178/10402065960840068379865,
c=-1: 179516485647379506696761/2706434410235514048893873,
c=+1: 2423380281636460862201/43432524347730591299828,
c=+2: 64566844993819835/1119348893518526181.
```

This is `FINITE-EXACT` evidence for a coherent-shift chamber, not an
all-shift theorem.  Its significance is the compensation absent from the
residue-blind pair claim: the companion uses the *actual perturbed singleton
debt* together with the whole five-edge tree, rather than demanding the old
canonical floor edge by edge.

### A cross-tree obstruction at coherent shift `c=84`

Pushing the same `3+3` packet to the coherent shift `e -> e+84` tests all
`12,980` assignments on the `649` bodies. Exactly one maximum cross-`K_(3,3)`
tree fails, while the unrestricted `K6` maximum tree has zero failures and
the literal union remains below `6/7`. Thus all `12,980` packets still close;
the unique obstruction is only to the restricted cross-tree certificate:

```text
E=(1,2,3,4,6,12), (L,j)=(168,90),
left residues=(86,88,90), right residues=(85,87,96),
cross margin=-52458412717/61674995021640,
full-tree margin=129336503/523645738,
union-6/7=-209657717/785468607.
```

At this witness `h=-84` and `hj=-45L`, so the centered defect is `Delta=0`.
The failure is therefore caused by the large slope and singleton change, not
by loss of centered phase alone; this is why `(2)` retains both `|Delta|` and
`|h|`. The phase-aware survivor is an unrestricted fixed-tree theorem. The
cross-only graph is a useful canonical certificate, not an invariant optimum
after large residue perturbations: within-class physical overlaps can become
the edges that pay the shifted singleton debt.

## 6. Boundary and next extension

The right next compiler is phase-budgeted, not radius-budgeted.  For each
canonical packet/tree it should retain

```text
(M_T, deg_T, z_i, h_i, centered Delta_i= h_i j mod L)
```

and accept the perturbation whenever `(3)` is positive, first with exact
`sigma_i`, then with `(2)` as an analytic tail.  Hostile search should rank
`|Delta_i|/L` before `|h_i|`; the one-unit witness proves the reverse ranking
is unsound.

Nothing here claims all noncanonical packets enter the perturbative chamber.
Packets with large centered phase, cell reselection, changes of the low graph,
arbitrary six-drift `k=1`, rung/entry, and LRC(14) remain open.

## 7. Reproduction

```bash
python3 04-computation/lrc14_reflected_residue_phase_stability_20260814.py
python3 -O 04-computation/lrc14_reflected_residue_phase_stability_20260814.py
python3 04-computation/lrc14_uniform_residue_shift_probe_20260814.py
python3 -O 04-computation/lrc14_uniform_residue_shift_probe_20260814.py
python3 04-computation/lrc14_shift84_tree_obstruction_20260814.py
python3 -O 04-computation/lrc14_shift84_tree_obstruction_20260814.py
```

The normal and optimized outputs are byte-identical. Hash basis:
LF-normalized bytes.

```text
fedaa46451ef1582e9d6a472e4996380ad08060810c1cf2642c4b9121ce3f383  phase lemma and controls
b41ef513589888543c3899ab0e11b684bbd0f3075f802c88adbdab09d223e4ef  phase output
efdb49aeffd8cf14907d7ab9f651d6c040bfa0969c070007c456f773087a7af2  shifts `-2,-1,1,2`
73303803fb884f5c69e8c9ddaf3d2284e0e1b540021bc6d9948f4fd28d48961f  shift output
be6a82621ab280861e820bbac4452197860dc1c6e6a02b2a86b19d6d8139339c  cross/full-tree obstruction
4a09c266825d663e672e6ce033a18691daa2414790393fc4ccfac5958604e996  obstruction output
```
