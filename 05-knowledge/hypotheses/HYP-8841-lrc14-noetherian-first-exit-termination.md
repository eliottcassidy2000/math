---
id: HYP-8841
title: "LRC14 Noetherian first-exit termination"
status: >
  OPEN / exact local no-go and bounded carrier tournament available. THM-2050
  proves that complete period-14 top germs are globally blind. THM-2043 proves
  every audited one-lift alias has a labelled strict exit by q<=42, including
  an infinite Hasse-indistinguishable family. THM-2051 adds a universal finite-
  height higher-relation gate: either there is a positive-measure strict exit or a
  support-three-through-five relation of coefficient height at most 2^21. The target is a
  height-decreasing termination theorem inside that structured branch, ending
  either at a signed Euler boundary certificate or at a strict rational lonely
  phase.
source: codex-2026-07-21-DC2-LRC14-termination
related:
  - THM-523
  - THM-597
  - THM-2043
  - THM-2047
  - THM-2049
  - THM-2050
  - THM-2051
  - THM-2052
  - HYP-2108
  - HYP-3117
  - HYP-3120
  - HYP-8845
script: 04-computation/lrc14_termination_sidecar_codex_20260721.py
output: 05-knowledge/results/lrc14_termination_sidecar_codex_20260721.out
---

# HYP-8841 -- termination, not another local invariant

For a finite speed set `S`, define the first strict exit

```text
q_exit(S)=min{q>=2: exists a, gcd(a,q)=1,
                       min_(v in S)||av/q||>1/14},   (1)
```

with value infinity if the set is empty.

THM-2050 shows why this sidecar is necessary: `AP13` and `12->26` have the
same complete local germs at every unit point `a/14`, but

```text
q_exit(AP13)=infinity,       q_exit(12->26)<=12.     (2)
```

THM-2043 supplies a stronger finite base case: all 156 one-lift aliases in its
scope either are the AP/Goddyn--Wong boundary rows or have an exact labelled
exit by denominator `42`; its infinite `12->96+3444n` family retains the same
exit `(q,a)=(41,17)` at every Hasse depth.

THM-2047 makes the strict branch finite without an experimental cutoff. If
`M(S)>1/14`, choose a maximizing phase `a/q` in lowest terms. Its opposite-sign
active pair gives

```text
q divides v_i+v_j <= 2 max(S),
```

and the same phase is strict. Consequently

```text
M(S)>1/14  ==>  q_exit(S)<=2 max(S).                 (2a)
```

Thus resolved exits are an adaptive finite sidecar. The remaining global
difficulty is exactly the non-strict branch: force a height-`1/14` Euler cell
when no pair-sum vertex lies above that height.

THM-2051 supplies the first universal terminal gate away from structured
circuits. Every row either already has a positive-measure strict lonely set or
obeys one of the finitely many templates

```text
sum_(i in A) k_i v_i=0,
3<=|A|<=5,      0<|k_i|<=2^21.                       (2b)
```

Thus a valid height descent may be localized to a finite union of rational
circuit hyperplanes. Each hyperplane still contains infinitely many speed
rows, and THM-2051 does not orient a decreasing move inside it.

There is nevertheless an exact rank cap. If the same primitive row acquires
twelve independent templates of height at most `H`, maximal minors and
Hadamard give `max v_i<=5^6 H^12`. At that point the problem is a finite exact
enumeration. The honest missing supplier is therefore “new independent
relation or exit,” not merely repeated detection of one relation.

THM-2052 proves that the rank is already high before any new descent is found.
THM-763 bounds a primitive counterexample by `sum v_i<=91^12`; a
pigeonhole/coding argument then supplies at least nine independent support-at-
most-five relations of height `2^20`, and at least eleven independent support-
at-most-three relations of height `91^6`. Thus every counterexample lies in a
finite atlas of rational subspaces of dimension at most two. More explicitly,
either the triple-relation rank is twelve and maximal minors give a finite
speed box, or two anchor coordinates `p,q` can be chosen so that each other
speed satisfies

```text
a_i v_p+b_i v_q+c_i v_i=0,      c_i!=0,
|a_i|,|b_i|,|c_i|<=91^6.                              (2c)
```

The eleven private-coordinate rows are independent. After common dilation is
removed, the entire residual is therefore a finite union of one-projective-
parameter two-anchor star families. Equivalently, bounded pair relations lock
coordinates into scale clusters and genuine cross-cluster triples leave at
most two cluster scales. The unresolved rank jump is from eleven to twelve,
or a direct Euler/strict-exit proof on those one-parameter stars.

The companion audit supplies the essential hostile control for that statement:
every one of its eight rows already has a coefficient-height-one circuit of
support three. The first seven, including AP, GW, `12->26`, `12->36`,
`12->96`, `12->84`, and P10+K33, all share the irrelevant relation

```text
1+2-3=0.                                                   (2d)
```

The genuine Cover14 tax-gain row similarly has `1+11-12=0`. Raw circuit
existence therefore ranks below even the local unit germ in the carrier
tournament. The next theorem must localize a circuit to the peeled speed, an
active endpoint owner, or a relation that changes the signed cyclic wall word;
a relation elsewhere in the row cannot drive deletion.

THM-2052 makes independence load-bearing too. A useful active-owner relation
must lie outside the already harvested rank-at-least-eleven code; otherwise it
does not descend. One new independent row raises the code to rank twelve, where
the maximal-minor bound makes the remaining decision finite. Thus the refined
target is **Euler survival or rank gain**, not merely active incidence.

The proposed LRC14 termination theorem is:

> Every primitive 13-speed row admits a finite sequence of labelled
> phase-height layer deletions/restrictions that strictly decreases a
> Noetherian height tuple and ends either at an AP/GW boundary certificate or
> at a strict exit (1).

A candidate height tuple is

```text
(unpaid q-witness count,
 first-exit denominator or infinity,
 maximum speed in the active magnitude fiber,
 active relation-layer rank,
 deleted-essential-point count).                   (3)
```

The order and even the correct coordinates in (3) are open; in particular,
`q_exit` is a terminal certificate rather than a demonstrated monotone under
deletion. What is fixed by
the DC(2) comparison is the proof architecture:

```text
local layer complex          usually solvable / acyclic
valuation or magnitude      records what localization forgot
strict height descent       prevents an infinite repair series
terminal atom               AP/GW boundary or explicit lonely phase
```

THM-2049 proves that local associated-graded solvability does not imply finite
polynomial closure.  For LRC14, THM-597 and THM-2050 give the corresponding
warning: local safe-component opening and complete top germs do not imply a
global lonely phase.  The missing theorem in both cases is termination.

## Exact carrier tournament and the resulting no-go

The companion audit completes the first concrete experiment on AP, GW,
`12->26`, `12->36`, `12->96`, covering `12->84`, P10+K33, and the incoming
genuine Cover14 tax-gain row. It
computes exact pair-sum maxima, the closed interval/point decomposition of
`G_{1/14}`, every first primitive exit through the complete bound (2a), and
THM-2048's peel tax. The sharp outputs include

```text
row          M          G volume        chi   first exit    tax fires
AP           1/14       0                6     none          no
GW 12->24    1/14       0                6     none          no
12->26       1/12       426/35035       10     (12,1,2)     no
12->36       3/41       1/1260           8     (41,17,1)    no
12->96       8/101      5219/840840     14     (41,17,1)    yes
12->84       7/89       563/105105       8     (41,17,1)    yes
P10+K33      2/27       4/2205          10     (27,8,1)     no
tax Cover14  7/43       589595/5213208  50     (16,7,12)    yes
```

This is table (4).

Among the first seven rows, all but `12->84` share the AP unit-germ digest.
The peel tax is therefore a
valid but selective strict-volume certificate: it catches the deep `96` and
covering `84` lifts, but misses the smallest hostile row, the K33 control, and
P10+K33. It cannot be the scalar termination height by itself.

The final row independently reproduces THM-2048's genuine Cover14 gain. At
peel `v=93` the core has

```text
mu=35517/280280, theta=19801/40040, r=50,
tax excess=2413467317/235670635200>0,                  (4a)
```

so the new tax succeeds exactly where the old uniform tail bound was
inconclusive. This cross-check both validates the interval implementation and
locates the tax honestly: it is strong on some deep Cover14 fibers, not on all
strict rows.

With proof carriers as tournament vertices, the observable is

```text
(cross-route merges, negative number of signature fibers, cost),
```

and the switch/gauge quotients the named bank by a carrier signature. The
tournament is transitive (score histogram `0,...,7`, no directed triangles,
singleton SCCs, one Hamiltonian path). Signed threshold topology wins; the
fused germ/tax/exit carrier is second; raw unlocalized circuit existence is
last. This challenges both the assumption that a deeper local packet should be
sought and the assumption that the THM-2051 circuit can be used without an
incidence sidecar. The missing item is a labelled deletion law that preserves
Euler points when the strict-volume tests do not fire.

## Next decisive target

For a primitive pair-sum phase define the integer signed margin

```text
C_(q,a)(S)=min_(v in S) (14 min(av mod q,q-av mod q)-q).   (5)
```

By THM-2047, LRC14 is equivalent to `max C_(q,a)(S)>=0` over the finite
pair-sum rulers `q | v_i+v_j`; positivity is the strict-exit branch and zero is
the Euler boundary branch. The next decisive theorem should prove a labelled
deletion/restriction alternative for a peel `S=C union {w}`:

1. a THM-2048 inequality is violated, giving positive safe volume; or
2. some pair-sum margin in (5) is positive; or
3. an owner-labelled endpoint of `G_{1/14}(C)` survives deletion by `w`.

The third clause is the only genuinely open Wall-A clause. A scalar residue,
collapse coefficient, volume moment, or unlabelled toric layer is an invalid
vertex because it merges THM-2050's hostile pair or misses AP's six isolated
Euler points.

THM-2051 is complementary rather than a replacement. Its proved
Fejer--BV/Bonferroni alternative shows that absence of a bounded
support-`3..5` integer relation yields positive safe volume. Every hard row
must therefore satisfy all three filters:

```text
all peel taxes pass,
a bounded genuine higher-support relation exists,
no positive pair-sum margin exists.                            (6)
```

The Euler endpoint-survival lemma only needs to be proved on this structured
intersection. THM-2052 further reduces it to a finite list of two-anchor star
templates with one projective parameter, unless rank twelve already gives a
finite box. This is a proved narrowing of the LRC reduction, but it does not
classify the one-parameter stars or orient a Noetherian move inside them. The
first necessary refinement is an **active-owner circuit lemma**: either a
bounded circuit meets the peeled/endpoint-owner set and lies outside the
already harvested rank-eleven code, or a different peel/strict exit is
available.

## Prior-art interface and exact inequality to prove

HYP-2108 already supplies the model active-owner gate in its multiple-speed
branch, while HYP-3117/HYP-3120 explain why endpoint ownership is a required
proof-circuit input. In the notation of the exact deletion identity, peel
`S=C union {w}` and write each positive-length component of
`G_{1/14}(C)` as

```text
I_i=[m_i-l_i/2,m_i+l_i/2].
```

Because a connected interval contained in the open danger set
`D_w={t:||wt||<1/14}` must lie in one danger component, direct endpoint
arithmetic gives

```text
I_i subset D_w
 iff ||w m_i|| + (w/2)l_i < 1/14.                    (7)
```

At equality an endpoint is weak-safe and already proves LRC14. Therefore the
third clause of the deletion trichotomy is exactly

```text
P_w(C):=max_i (||w m_i||+(w/2)l_i-1/14) >= 0.        (8)
```

This is the concrete next decisive target. Let `W` be an eleven-dimensional
bounded relation code supplied by THM-2052:

> On every residual row for which all THM-2048 peel taxes pass and all pair-sum
> margins are nonpositive, either some peel `w` satisfies `P_w(C)>=0`, or the
> active endpoint-owner data supplies a bounded relation `k dot v=0` with
> `k notin W`.

Equations (7)--(8) do not prove the target; they identify the exact missing
input. The raw relation `1+2=3` does not mention `w`, `m_i`, or the endpoint
owners and hence cannot control `P_w`. A successful circuit-elimination step
must manufacture an independent incidence row or emit a strict exit. This
recovers, in the new Fejer/phase-height setting, the older endpoint-cover
circuit rather than inventing another scalar carrier.

HYP-8845 supplies a symmetry dividend on the covering branch. The involution
`t->1-t` acts freely on `G_{1/14}` when an even speed is present, so a surviving
endpoint/component comes with its mirror and `chi` is even. It is enough to
force (8) on one half-circle; the mirror produces the second Euler component.
This halves the owner-word problem but does not supply the rank gain or the
first survivor.
