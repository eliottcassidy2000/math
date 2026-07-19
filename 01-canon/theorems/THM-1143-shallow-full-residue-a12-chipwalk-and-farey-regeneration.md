---
id: THM-1143
title: The all-height shallow full-residue predicate is an A_12 mechanical ballot walk; AP dilates are affine copies of one Farey-12 comb
status: PROVED carrier equivalence and AP self-similarity; FINITE-EXACT independent h in {0,1,2} classification; arbitrary-height ballot rigidity remains OPEN
source: codex-2026-07-18-S74
depends_on:
  - THM-769   # identifies the shallow tight branch with full nonzero residues
related:
  - THM-770   # stronger bounded classification through height 12, by a different owner CSP
  - THM-795   # all-height Hamming-one AP rigidity
  - THM-800   # all-height Hamming-two AP rigidity
  - THM-1142  # essential-region criterion; no global descent bridge is claimed
  - HYP-6820 # open n=12 sporadic branch
script: 04-computation/lrc13_shallow_a12_chipwalk_descent_audit_codex_20260718.py
output: 05-knowledge/results/lrc13_shallow_a12_chipwalk_descent_audit_codex_20260718.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCA12Chipwalk.lean
---

# THM-1143 -- the shallow object is a mechanical root ballot walk

## 1. Domain and exact statement

Let

```text
W=(w_r)_(r in F_13^*),       w_r>0,       w_r == r (mod 13),
```

be a residue-labelled full-residue packet.  Put

```text
G_13(W)={t in R/Z : ||w_r t||>1/13 for every r},
chi_13(W)=#pi_0(G_13(W)).
```

Write every time uniquely as

```text
t=(j+u)/13,       j in F_13,       0<=u<1.
```

Away from a wall `u=k/w_r`, define the two-element sheet edge

```text
E_r(u)={-floor(w_r u) r^(-1),
        -(floor(w_r u)+1) r^(-1)} subset F_13.             (1)
```

Let `d_j(u)` be the number of these twelve edges incident to sheet `j`, and
let `e_j(u)=d_j(u)-1`.

> **All-height carrier equivalence.**  The packet has `chi_13(W)=0` if and
> only if, after every group of simultaneous walls, the vector `e(u)` is
> coordinatewise nonnegative.

Initially

```text
e(0+) = 11 delta_0.                                       (2)
```

At a wall `u=k/w_r`, the edge slide changes the excess by the `A_12` root

```text
Delta_(r,k)
  = delta_(-(k+1)r^(-1)) - delta_(-(k-1)r^(-1)).          (3)
```

All roots at the same rational wall are summed before testing.  Consequently

```text
chi_13(W)=0
  iff 11 delta_0 + sum_(walls <= u) Delta_(r,k) >= 0       (4)
```

at every grouped prefix of the arithmetic-realizable mechanical word.
Every surviving state is a nonnegative 13-vector of total mass `11`; hence
the predicate-preserving state simplex has exactly

```text
C(11+13-1,13-1)=C(23,11)=1,352,078                       (5)
```

states, independently of the heights of the speeds.  The word can still be
arbitrarily long: (5) is a finite-state quotient, not a finite-height bound.

## 2. Proof of the carrier equivalence

Fix a generic `u` and write `w_r u=m+x`, with `m` integral and `0<x<1`.
Modulo one,

```text
w_r (j+u)/13 = (rj+m+x)/13.
```

This has distance strictly below `1/13` precisely when
`rj+m` is `0` or `-1` modulo `13`.  Solving for `j` gives exactly (1).
Thus a generic time over sheet `j` is blocked if and only if `d_j(u)>=1`.
It is a strict witness if and only if `e_j(u)<0`.

At `u=k/w_r`, the strict danger edge collapses to its common vertex, while
the two other endpoints of the left and right generic edges have clearance
exactly `1/13`.  They therefore cannot be strict witnesses.  Equivalently,
the closed obstruction to `||w_rt||>1/13` at the wall is the union of the two
adjacent generic edges.  Hence checking every generic chamber is complete;
simultaneous walls must be crossed as one group rather than in an arbitrary
tie order.

For `u=0+`, runner `r` has edge `{0,-r^(-1)}`.  The twelve inverse residues
run through `F_13^*`, giving degrees `(12,1,...,1)` and (2).  Immediately
before and after `k/w_r`, the edge is respectively

```text
{-(k-1)r^(-1),-kr^(-1)}
and
{-kr^(-1),-(k+1)r^(-1)}.
```

Their incidence difference is (3).  Every slide preserves total excess;
the terminal vector is `11 delta_(-1)`.  This proves (4) and (5).

The proof uses the original continuous circle only to derive the edge
language.  Thereafter the full predicate is the nonnegative-prefix language
of a labelled mechanical `A_12` root word.

`LRCA12Chipwalk.lean` kernel-checks the finite-state arithmetic used here:
edge-slide cancellation, zero root mass, commuting tied transports, grouped
mass preservation, and the affine mass-eleven prefix invariant.  It contains
no `sorry` or `native_decide` and builds from the project root.  The analytic
floor/danger-edge equivalence in this section and the finite classification
remain external; the module does not pretend to formalize the open ballot
rigidity statement.

## 3. AP dilates are a self-similar Farey comb

For a positive integer `c` prime to `13`, label the dilation
`c{1,...,12}` by residues: the speed in class `r` is `ca`, where
`a=c^(-1)r (mod 13)` is represented in `{1,...,12}`.  For

```text
u=(m+v)/c,       0<=m<c,       0<v<1,
```

direct substitution in (1) gives

```text
E_(ca)((m+v)/c)
   = c^(-1) (E_a(v)-m).                                   (6)
```

Thus every block `[m/c,(m+1)/c]` is an affine relabelling of the scale-one
word.  The scale-one walls are the Farey fractions of order twelve, giving
46 chambers.  Exact graph enumeration on that base comb gives

```text
24 tree chambers,
22 chambers with two components and cycle rank one,
2 star chambers.                                          (7)
```

Equation (6) proves that scale `c` has `c` copies of (7): `46c` chambers,
`24c` tree chambers, `22c` two-component/unicyclic chambers, and `2c` star
chambers.  This is the precise toothpick self-similarity of the equality
rays; it is not merely a picture seen at small scale.

If `c=13q+s`, `1<=s<=12`, and `pi_s(r)` is the representative of
`s^(-1)r` in `{1,...,12}`, the corresponding height vector is explicitly

```text
h_r = q pi_s(r) + floor(s pi_s(r)/13).                    (8)
```

Conversely, division by `gcd(W)` and residue relabelling preserve the circle
cover predicate.  The desired global rigidity can therefore be stated in
primitive form.

## 4. Exact computation and independent controls

The companion script implements the grouped root walk and also a separate
atomic-chamber circle-cover checker.  They agree on AP, `2AP`, `3AP`, and
three deliberately torn controls.  A direct-mask enumeration, independent
of the root-update implementation, checks all

```text
3^12 = 531,441
```

packets with `w_r=r+13h_r`, `h_r in {0,1,2}`.  Exactly three survive:

```text
{1,...,12},       2{1,...,12},       3{1,...,12}.          (9)
```

The 426-chamber decision digest is

```text
31880bc101e49a1682af8de9583ca68aeb428a91a854f6bd2bfc6fa98198d018.
```

This is an independent low-height control for THM-770, whose owner-CSP gives
the much stronger bounded classification through height twelve.  The new
frontier contribution of this theorem is the all-height carrier (4), not the
smaller finite box (9).

Normal and optimized Python runs are byte-identical.  Frozen hashes are

```text
source  bfeb76ddf879b18125b4ed219cf7811391f4766ac579311161331d0959893919
output  f80ab6acaa0d233632c400de7e4dc650f7c55c946a8871e2910271cf18184c71
```

## 5. A tempting descent is exactly false

A coordinate lowering replaces `w_r` by `w_r-13`.  It is tempting to order
failed packets by the time of their first negative coordinate and descend by
one such lowering.  Exhaustion of `h_r in {0,1}` finds two exact regressions:

```text
(1,2,3,4,18,6,20,21,22,10,24,12),     first tear 1/6;
(14,2,16,4,18,6,7,8,22,10,11,12),     first tear 2/7.
```

For the first packet, every legal lowering tears earlier, at respectively
`1/11,3/22,1/10,1/12,2/21`.  For the second, every legal lowering tears
earlier, at `1/6,1/7,1/11,1/8`.  Moreover, the terminal current is universal:
all 4,096 packets with heights in `{1,2}` have the correct endpoint current,
yet none is accepted.

Therefore first-tear time, endpoint current, and any scalar summary that
forgets the whole prefix-minimum stalk are not descent energies.  This is a
negative theorem about a proof route, not evidence against rigidity.

## 6. The exact remaining lemma

The all-height shallow problem is now the following arithmetic-language
statement.

> **A12 ballot rigidity (OPEN).**  Every arithmetic-realizable packet whose
> grouped mechanical root word satisfies (4) is `c{1,...,12}`.

Equivalent useful targets are:

1. **primitive regeneration:** if `gcd(W)=1` and (4) holds, then
   `W={1,...,12}`;
2. **Farey regeneration:** every accepted primitive mechanical word is the
   single 46-chamber Farey-12 comb;
3. **sink rigidity:** every accepted sink under all legal coordinate
   lowerings is a dilation ray, together with a proof that every non-ray
   accepted packet descends to such a sink.

The first-tear counterexamples show that target 3 needs a vector- or
word-valued potential retaining prefix chronology.  The finite simplex (5)
suggests an automaton proof, but one must also characterize which cycles and
returns are realizable by the twelve coupled Beatty/Farey event streams.

THM-769 says this lemma would close the entire *shallow* n=12 sporadic branch.
It would not close the deep sheet branch, the compact LRC(14) residual of
HYP-7665, or LRC(14) itself without an additional proved bridge.

## 7. Tournament and carrier audit

The predicate-preserving vertices are the thirteen sheet cuts, not the twelve
runners.  Orienting two cuts by their current deficit gives a transitive
telemetry tournament (score sequence `0,...,12`, no directed cycles, thirteen
singleton SCCs); root negation with `u -> 1-u` is the switch/gauge and
`0,1,...,12` is a tie Hamiltonian path.  That tournament forgets chronology,
so it cannot decide (4).

Alternatively, take one vertex for each labelled wall `(r,k)`, order pairs by
`k/w_r`, and contract exact ties.  This retains the Christoffel interleaving
but has unboundedly many vertices.  The grouped `A_12` word is the finite-state
quotient between these extremes:

- it preserves the sheet-cover predicate and chronology;
- retaining `(r,k)` labels preserves arithmetic realizability;
- deleting labels loses runner identity;
- keeping only the cut tournament loses wall order;
- allowing arbitrary root words enlarges the language beyond mechanical
  packets.

The challenged assumption is therefore explicit: neither runners nor danger
arcs are the right primary vertices.  The underlying object is a labelled
wall chronology acting on a finite sheet-deficit simplex.
