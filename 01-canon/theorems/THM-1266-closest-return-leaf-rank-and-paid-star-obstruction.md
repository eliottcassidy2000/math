---
id: THM-1266
title: CLOSEST RETURNS HAVE A 6/5 LEAF RANK, AND A CONSECUTIVE TOOTHPICK STAR HAS AT MOST FIVE RUNGS
status: PROVED (closest-return distinct-internal leaf; exact 6/5 metric ascent; height-five fixed-owner and height-42 compact transported ranks; recursive disjoint-packet packing; additive return stalk; exact repeated-low endpoint separation; sharp five-rung consecutive-address star bound; primitive c=140 centered-blocker/two-cycle realization; exact four-fastest-safe-gap tail audit; optimization-safe referee; sorry-free Lean arithmetic consumer). This terminates the local r=1 star and localizes its tail; it does not prove six-comb noncoverage or LRC(14)
source: codex-2026-07-19 return-recursion audit
depends_on: [THM-1252, THM-1253, THM-1256, THM-1262, THM-1264]
related: [THM-841, THM-1196, THM-1199, THM-1233, THM-1254]
script: 04-computation/lrc14_closest_return_leaf_paid_star_thm1266.py
output: 05-knowledge/results/lrc14_closest_return_leaf_paid_star_thm1266.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCClosestReturnLeafPaidStar.lean
script_sha256: 8e6f364e300c1b8c052578d4c644a2d334b0cb9af75266aa317966cf0613b7ab
output_sha256: b1907d7c4ec43379739a46f77820c1dec7b4fbdd10fcbaabc99314133fb5e7de
formalization_sha256: e5f78b09ccdfa6d8fcef09018e8b3498caeb19ef7ed29a206141d8e45a39a8a3
---

# THM-1266 -- closest-return leaf rank and the sharp paid-star tail

## 1. Closest literal returns are short leaves

Use THM-1253's deletion-minimal chronological cover word

```text
(I_0,...,I_(N-1)),          I_q=T(s_q,n_q),
T(s,n)=((14n-1)/(14s),(14n+1)/(14s)).                 (1)
```

Adjacent teeth overlap and have different owners.  Suppose a contiguous
block of the word contains a repeated owner.  Among all equal-owner pairs in
that block choose positions `p<q` of minimum distance and put

```text
a=s_p=s_q,             m=q-p,             R=n_q-n_p. (2)
```

The outer teeth are distinct teeth of one comb, so `R>=1`.  No internal
owner equals `a`, by closestness.  Nor can two internal owners agree, since
they would give a closer equal-owner pair.  There are only six fast owners;
hence

```text
2<=m<=6,                                                 (3)
```

and the `m-1` internal owners are distinct.

THM-1264's exact identity applies because this is a **literal consecutive
subword**, not merely a cycle in the unordered support graph.  If `omega_r`
is the positive overlap of teeth `r,r+1`, then

```text
sum_(r=p)^(q-1) omega_r
 = (1/7)sum_(r=p)^(q-1)1/s_r - R/a >0.                (4)
```

After multiplication by `7a`, (4) gives

```text
sum_(r=p+1)^(q-1) a/s_r > 7R-1.                       (5)
```

Let `b` be the minimum internal speed.  Every summand in (5) is at most
`a/b`, so (3), `R>=1`, and (5) imply

```text
a/b > (7R-1)/(m-1) >= 6/5.                            (6)
```

Thus every repeated block contains a short **return leaf** carrying a strict
metric witness

```text
b -> a,                       a>(6/5)b.                (7)
```

This is the weakest all-length return ascent.  The two-edge toothpick has
the stronger factor `a>6b`; a three-edge return has factor greater than
three.

## 2. The rank is well-founded, but transport is not automatic

Orient every witness as in (7).  Speed strictly increases, so the witness
digraph is acyclic.  On the fixed six owners, any directed path has at most
five edges.  If a future transport construction introduces new owner speeds
but keeps all of them in THM-1233's projective box

```text
c<s<2345c,                                                (8)
```

use the integral rank

```text
rho(s)=max{r in N : (6/5)^r c <=s}.                      (9)
```

Every edge (7) raises `rho` by at least one, and the exact comparison

```text
(6/5)^42<2345<(6/5)^43                                  (10)
```

bounds a transported witness path by 42 edges.

Equation (10) is a conditional recursion height, not yet a recursive proof.
A closest leaf has no repeated owner inside it.  Recursing on the left and
right exterior blocks need not lower `rho`, and THM-1262's third-owner bridge
does not yet say that the next return uses the witness `b` in (7).  This is
the precise gap between a well-founded coordinate and a well-founded
transport.

## 3. Long words contain linearly many paid leaves

There is nevertheless a genuine disjoint packing.  On an original
contiguous word block:

1. stop if its owners are distinct;
2. otherwise remove all vertices of one closest return packet; and
3. recurse on the original left and right blocks.

The selected packets are vertex-disjoint, hence their chronological seam
indices are disjoint.  Each packet uses at most seven teeth.  If `P` packets
are selected, their removal leaves at most `P+1` terminal blocks, each with
at most six teeth.  Consequently

```text
N <=7P+6(P+1)=13P+6,
P >=ceil((N-6)/13).                                    (11)
```

Every selected packet carries (7).  There are only 15 possible strict
ordered speed pairs among six owners, so at least `ceil(P/15)` packet
witnesses use one pair.

This breadth is already paid by THM-1253.  For any family of packets whose
seam occurrences are disjoint,

```text
sum_(packet seams {u,v}) 1/lcm(u,v)
 <=(12/7)(H_fast-1/c),                                (12)

sum_(packet seams {u,v}) 1/lcm(u,v)
 <=16F_6/c.                                           (13)
```

Equations (12)--(13) classify existing chronological mass; they do not
create a second invoice.

There is an exact contraction law, but its state is not a naked owner.  For
an `a`-return packet put

```text
Stalk(P)=(a,R,D),              D=sum_(seams in P)omega. (14)
```

Two consecutive `a`-returns concatenate by

```text
(a,R_1,D_1)*(a,R_2,D_2)=(a,R_1+R_2,D_1+D_2).          (15)
```

The shared `a`-tooth contributes exactly the extra `1/(7a)` needed in the
long return identity.  Replacing a packet only by the symbol `a` is invalid:
it deletes `R`, the paid seam mass, both wall germs, and the interval whose
coverage still has to be transported.  The faithful recursive state must be
at least the stalk (14) plus those boundary germs.

## 4. Consecutive-address toothpick stars terminate at five

The remaining danger is breadth concentrated at one high owner.  Consider
an irredundant subword

```text
H_0,J_0,H_1,J_1,...,J_(p-1),H_p,                      (16)
```

where every `H_i` has speed `h` and address `M+i`.  Thus rung `i` fills the
one complete `h`-safe gap between consecutive `h`-teeth.  Write `j_i` for
the owner of `J_i`.  Each triple

```text
H_i,J_i,H_(i+1)                                        (17)
```

is an `r=1` immediate return, so THM-1256 gives

```text
h>6j_i.                                                (18)
```

Suppose two low slots `i<k` have the same owner `j`, and put `Delta=k-i`.
Their `j`-teeth are distinct, so their left endpoints differ by at least
`1/j`.  Endpoint order and the overlap of `J_k` with `H_k` give

```text
1/j <= alpha(J_k)-alpha(J_i)
    < beta(H_k)-alpha(H_i)
    =(14Delta+2)/(14h)
    =(7Delta+1)/(7h).                                  (19)
```

Combining (18)--(19) yields

```text
Delta>=6.                                              (20)
```

There are only five owners other than `h`.  Six consecutive low slots would
repeat one owner at distance at most five, contradicting (20).  Therefore

> **Five-rung star theorem.**  A six-owner irredundant tooth word contains
> no consecutive-address toothpick star with six rungs.  The sharp bound is
> `p<=5`.

This strengthens THM-1256's `ABAB` exclusion.  That theorem forbids one low
owner from returning in the next star slot; (20) forbids it throughout the
next five slots.

The breadth remains quantitatively paid.  If rung `i` has low owner `j_i`,
the two seams have exact total mass

```text
Omega_i=[h-6j_i]/(7hj_i)
       >=1/[7 lcm(h,j_i)].                             (21)
```

All `2p` seams are pairwise disjoint.  In a full six-cover word, (12) becomes

```text
sum_i 1/lcm(h,j_i) <=(6/7)(H_fast-1/c),               (22)
```

with functional form `sum_i 1/lcm(h,j_i)<=8F_6/c`.

The theorem concerns the `r=1`, consecutive-fast-address star.  General
nonconsecutive returns and turns away from the star remain possible, though
their seams are still paid by THM-1253.

## 5. A primitive exact row attains all five rungs

The bound is sharp in the compact LRC geometry.  Put

```text
c=140,             k=80,
G=[1121/1960,1133/1960],
{d_1,...,d_6}={254,255,256,257,292,1805}.              (23)
```

The following eleven teeth, listed as `(speed,address)`, form a strict
irredundant chain wholly inside `G`:

```text
(1805,1036), (256,147), (1805,1037), (254,146),
(1805,1038), (292,168), (1805,1039), (257,148),
(1805,1040), (255,147), (1805,1041).                  (24)
```

Both endpoint lists in (24) increase strictly, every adjacent pair
overlaps, and

```text
beta(I_i)<alpha(I_(i+2))                               (25)
```

for all nine applicable indices.  The minimum separation in (25) is
`11/1050616`.  The five detunings `1805-6j_i` are

```text
(269,281,53,263,275),                                  (26)
```

so every return is strict.  The ten seam occurrences have total mass

```text
138130160177/393047548696320.                          (27)
```

The backtrack/turn ledger is

```text
B,T,B,T,B,T,B,T,B,              (B,T)=(5,4),          (28)
```

simultaneously attaining both sharp bounds in THM-1256.  The row is
primitive, `1805<2345*140`, and even passes the first harmonic gate:

```text
sum_i 1/d_i -1/140
 =981965948779/78609509739264>0.                       (29)
```

Thus neither primitive reduction, compactification, the scalar harmonic
gate, nor the local word ledger removes it.

## 6. The sharp row also realizes the blocker-cycle package

The six THM-1240 centered phases in (23) are

```text
t_254 =227/394,       t_255 =227/395,
t_256 =19/33,         t_257 =228/397,
t_292 =31/54,         t_1805=1118/1945.               (30)
```

Every phase is strictly contained in a tooth of the one chain (24).  Choose
the following blocker targets, where `@r` is the zero-based word position:

```text
254  ->1805@8,        255 ->254@3,
256  ->257@7,         257 ->256@1,
292  ->256@1,         1805->254@3.                    (31)
```

All target positions are internal, so both THM-1252 wall neighbors occur in
the same word.  The functional graph has two cycles,

```text
254 <->1805,                  256<->257,               (32)
```

with tails `255->254` and `292->256`.

Both cycles realize, rather than evade, THM-1262.  In the first cycle the
ascent and reverse targets have word positions `8,3`,

```text
t_1805-t_254=-1023/766330<0,                           (33)
```

and the corridor-facing third owner at position `7` is `257`.  In the
second cycle the target positions are `7,1`,

```text
t_257-t_256=-19/13101<0,                               (34)
```

and the third owner at position `6` is `1805`.  Thus tooth order and phase
order agree in both cases, the marked targets are disjoint and
nonconsecutive, and the aligned third-owner bridges are literal.

The target `1805@8` has distinct wall owners `257,255`, giving THM-1252's
rank-two fork.  Targets `254@3`, `257@7`, and `256@1` have the same high
owner at both walls and give the exact toothpick branch.  Therefore the
closest-return, fork/toothpick, binary landing, aligned-cycle, and protected
bridge laws are mutually compatible in one primitive exact cell.

## 7. The full-comb sweep leaves exactly four fastest-safe tails

Sweep **all** teeth of the six speeds in (23) that meet `G`, rather than only
the selected eleven.  Exactly 20 teeth meet the gap.  Their strict open
danger union leaves four components:

```text
[14477/25270,14489/25270],
[14491/25270,14503/25270],
[2915/5054,14587/25270],
[14589/25270,14601/25270].                             (35)
```

Each has length

```text
6/12635=6/(7*1805),                                   (36)
```

and the total uncovered mass is `24/12635`.  They are precisely the complete
`1805`-safe gaps after addresses

```text
1034, 1035, 1041, 1042.                               (37)
```

The five slower owners fill the five consecutive fastest-safe gaps between
addresses `1036` and `1041`, which is the sharp star (24), but they miss the
two immediately before and two immediately after.

This row is **not** a six-cover and is not a counterexample to LRC.  Its
importance is diagnostic: all current local return and blocker-cycle laws
hold, and the only failure is global tail completion.  The `r=1` ladder has
terminated exactly where the owner count says it must.

## 8. Remaining implication and the highest-leverage target

THM-1266 rules out an indefinitely self-similar consecutive toothpick
ladder.  It does not show that every hypothetical six-cover contains six
consecutive fastest-gap rungs.  A full word may instead

* jump by a return address `r>1`;
* leave the high-owner star through a nonbacktracking turn;
* reuse a low owner after six or more fastest-gap addresses; or
* cover a tail by a tooth which is not part of the local alternating word.

The next global consumer should therefore prove one of the following.

1. **Star completion:** a six-cover extending a maximal five-rung star must
   fill one adjacent fastest-safe gap, forcing either a forbidden sixth rung
   or a quantified non-star turn.
2. **Tail tax:** both end extensions contribute phase-located mass outside
   the seams already counted in (22), enough to violate the `H`-drift or the
   `j=4` component-span budget.
3. **Transport landing:** a THM-1262 bridge lands on a new return stalk whose
   rank is linked to (7), making the height-five/42 coordinate genuinely
   recursive.

The exact row (23)--(37) is a positive-control test for all three proposals.
A claimed local closure which rejects the row is using an unstated global
hypothesis; a claimed global closure must explain one of its four tails.

## 9. Tournament and alternate-carrier audit

The runner tournament on the six speeds is transitive.  Its score histogram
is `(0,1,2,3,4,5)`, it has no directed cycles, six singleton SCCs, no edge
flips in the speed gauge, and one Hamiltonian path.  It cannot see (24), its
ten seam occurrences, or the four tails.

For closest returns, use **return packets** as vertices.  The pair observable
is the strict metric witness `b->a` in (7); chronology is the tie Hamiltonian
path.  In the sharp row the witness digraph is a five-arrow in-star with
center `1805`, no directed cycle, and depth one.  This is the exact reason
the five ascent factors do not multiply.

For the star bound, use **fastest-gap slots** as vertices.  The switch is
equality of their low owner.  Equation (20) says equal colors cannot occur
within five chronological steps; the conflict graph on any five consecutive
slots is `K_5`.  Five non-high owner colors suffice sharply, while a sixth
slot is impossible.  The packet-order tournament on the five realized rungs
is again transitive, with scores `(0,1,2,3,4)`, no cycles, singleton SCCs,
and one chronological Hamiltonian path.  Collapsing packets to the owner
pair `{j,1805}` destroys their positions; collapsing positions to runners
destroys detuning and address drift.

We explicitly challenged runners, gaps, gap boundaries, selected teeth,
wall events, seam cells, return packets, address stalks, residues, and proof
obligations as vertices.  The smallest faithful state for recursion is

```text
(outer owner; first/last tooth addresses; R; ordered internal owners;
 left/right wall germs; paid seam occurrences; carrier subinterval).      (38)
```

It preserves the cover predicate, exact drift, chronology, and invoice.
Runner tournaments preserve none of those simultaneously.

## 10. Verification and scope

The optimization-safe exact referee has zero Python `assert` nodes.  It
checks `585,936` proper six-owner words, including `583,980` closest-return
words; `86,856` exact ratio rows; `118,096` packet-packing words with
`235,356` selected packets; `75,750` forbidden repeated-low separation rows;
all endpoint, seam, detuning, centered-phase, blocker-map, alignment, and
bridge data of (23)--(34); and the exact 20-tooth full-comb sweep (35)--(37).
Normal and `python -O` outputs are byte-identical and match the stored
certificate.

The sorry-free Lean module kernel-checks the `6/5` consumer, packet count,
rank cutoff, repeated-low separation, six-slot pigeonhole contradiction,
all sharp-chain endpoint comparisons, compact/harmonic checks, centered
blocker membership, two-cycle alignment signs, and the four exact tail-gap
lengths.  It contains no `sorry`, no `admit`, and no `native_decide`.
Minimal-subcover extraction, closest-pair topology, identification of the
full-comb uncovered components, and transport from one return stalk to the
next remain explicit paper/referee providers.

Frozen artifact hashes are

```text
source         8e6f364e300c1b8c052578d4c644a2d334b0cb9af75266aa317966cf0613b7ab
output         b1907d7c4ec43379739a46f77820c1dec7b4fbdd10fcbaabc99314133fb5e7de
formalization  e5f78b09ccdfa6d8fcef09018e8b3498caeb19ef7ed29a206141d8e45a39a8a3
```

THM-1266 proves local star termination and exact tail localization.  It does
not prove uniform six-comb noncoverage, the empty sporadic branch, or
LRC(14).  ∎
