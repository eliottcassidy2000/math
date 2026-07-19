---
id: THM-1248
title: THE CENTERED BLOCKER CYCLE HAS A FINITE RELATIVE-ADDRESS WORD AND EXPORTS A GCD-QUANTIZED WALL SEAM
status: PROVED (exact lower-gap cocycle and translation gauge; least-positive determinant remainder and gcd sheet; uniform 1,174-symbol relative-address bank; binary chronological digits on every target-at-most-source-clock edge; additive and contracting affine cycle holonomies; phase-located fourth-support wall seam on every cycle ascent; slowest-orbit handoff anchor; exact guardrails; dependency-free exact referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S82 continuation with address-compression, cycle-model, and package-audit agents
depends_on: [THM-1233, THM-1237, THM-1240, THM-1244]
related: [THM-1198, THM-1226, THM-1242, THM-1243, THM-1246, THM-1247, THM-1250, HYP-7870]
script: 04-computation/lrc14_centered_blocker_address_compression_thm1248.py
output: 05-knowledge/results/lrc14_centered_blocker_address_compression_thm1248.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCenteredBlockerAddressCompression.lean
script_sha256: d236f1ec6b48a95e7dfc3a7f3beb5a68061aa41a61a2eb01e3b1ba38827af704
output_sha256: ba1b6323d40c34d0b204d75f4658f67fa9d053ced7a9d2b45c6cf32298f2f8bb
formalization_sha256: 54487de2f3886570b2752f044457d2d174edc1b2c70a3067eac397254b06ae35
---

# THM-1248 — centered blocker address compression and wall export

## 1. The lower-gap blocker cocycle

Let

```text
G=G_k(c)=[(14k+1)/(14c),(14k+13)/(14c)],   0<=k<c,   (1)
```

be a complete closed safe gap of the positive integer speed `c`.  Suppose
the six strict danger combs of

```text
c<d1<d2<d3<d4<d5<d6                                  (2)
```

cover `G`.  Put

```text
t0=(k+1/2)/c,             Q_i=c+d_i.                  (3)
```

For every `i`, choose an integer `P_i` nearest `Q_i t0`, set

```text
t_i=P_i/Q_i,
epsilon_i=P_i-Q_i t0 in [-1/2,1/2],                  (4)
```

and choose one THM-1240 blocker `j=b(i)`.  Thus `j!=i` and `d_j` is
strictly dangerous at `t_i`.  There is a unique integer `N_i` such that

```text
beta_i=d_j t_i-N_i,               |beta_i|<1/14.      (5)
```

The raw tooth address `N_i` is not the correct invariant.  Define instead

```text
X_i=cP_i-kQ_i,
r_i=P_i d_j-N_iQ_i=Q_i beta_i,
M_i=k+N_i,
ell_i=P_iQ_j-M_iQ_i,
delta_i=P_j-M_i.                                      (6)
```

Direct expansion using `Q_j=c+d_j` gives the exact cocycle

```text
ell_i=X_i+r_i.                                        (7)
```

THM-1240 says

```text
Q_i/4<X_i<3Q_i/4.                                    (8)
```

Combining (5), (7), and (8) produces the sharp positive central band

```text
5Q_i/28<ell_i<23Q_i/28.                              (9)
```

In particular `0<ell_i<Q_i`; no modular choice is hidden in (9).

## 2. Gauge removal, determinant remainder, and gcd sheet

An integral translation of time by `n` changes the addresses by

```text
P_i ->P_i+nQ_i,       P_j->P_j+nQ_j,
k   ->k+nc,           N_i->N_i+nd_j,
M_i ->M_i+nQ_j.                                      (10)
```

The four quantities `X_i,r_i,ell_i,delta_i` are unchanged.  Thus the
absolute lift `N_i` is gauge, whereas `delta_i` is a genuine relative tooth
address.

Let

```text
D_ij=P_iQ_j-P_jQ_i.                                  (11)
```

Then (6) gives

```text
D_ij=ell_i-delta_iQ_i.                               (12)
```

By (9), `ell_i` is exactly the least positive residue of `D_ij` modulo
`Q_i`, and

```text
M_i=floor(P_iQ_j/Q_i),
gcd(Q_i,Q_j) divides ell_i.                           (13)
```

The second statement follows before any reduction: every common divisor of
the two clocks divides both terms in `P_iQ_j-M_iQ_i`.  Consequently the
faithful sampled edge is not just a directed pair of runners but

```text
(Q_i,Q_j; least-positive ell_i; relative delta_i;
 gcd(Q_i,Q_j)-sheet).                                 (14)
```

It preserves blocker danger, phase chronology, and exact torsion while
discarding only the integral tooth lift.

## 3. The finite relative-address bank

Substitute (3)--(6) and `ct0=k+1/2`.  All absolute addresses cancel:

```text
beta_i=delta_i-1/2+(d_j/Q_i)epsilon_i-epsilon_j,      (15)
epsilon_j=(d_j/Q_i)epsilon_i+delta_i-1/2-beta_i.      (16)
```

Writing `a_i=d_j/Q_i`, the rounding and danger bounds imply

```text
-a_i/2-1/14<delta_i<1+a_i/2+1/14.                    (17)
```

THM-1233 gives `d_j<2345c`, while `Q_i>2c`.  Hence
`a_i<2345/2`; since `delta_i` is integral, (17) yields the uniform bank

```text
-586<=delta_i<=587.                                  (18)
```

There are therefore only `1,174` possible relative address symbols on an
arbitrary blocker edge.  More sharply, if

```text
d_j<=Q_i=c+d_i,                                      (19)
```

then `a_i<=1`, so (17) lies strictly inside `(-1,2)` and

```text
delta_i in {0,1}.                                    (20)
```

Every speed-descent edge has this binary property, but (19) also includes
many ascent edges.

Put `theta_i=ell_i/Q_i`.  Equations (6) and (12) give

```text
Q_j(t_i-t_j)=theta_i-delta_i.                         (21)
```

When (20) holds, the digit records which chronological side of `t_j` the
source phase occupies, and (9) becomes the quantitative separation

```text
5/(28Q_j)<|t_i-t_j|<23/(28Q_j).                      (22)
```

Thus a binary edge is an adjacent-centered-tooth transport, not a bare
orientation.

## 4. Exact cycle invoices

Let

```text
i_1->i_2->...->i_m->i_1,              2<=m<=6,       (23)
```

be a directed cycle of the selected blocker map.  Summing (21) around the
cycle gives the additive invoice

```text
sum_s delta_s/Q_(s+1)=sum_s theta_s/Q_(s+1).          (24)
```

After division by `sum_s 1/Q_(s+1)`, the right side is a weighted average
with every `theta_s in (5/28,23/28)`.  Therefore every cycle contains both

```text
some delta_s<=0,                    some delta_r>=1.  (25)
```

The positive remainders also have a multiplicative form.  In the canonical
gap gauge, `P_i>0`.  Moreover `N_i>0`: since `t_i>1/(14c)` and `d_j>c`,
one has `d_jt_i>1/14`, so (5) excludes `N_i=0`; hence `M_i=k+N_i>0`.
Each identity

```text
P_iQ_j=M_iQ_i+ell_i>M_iQ_i                           (26)
```

can now be multiplied around (23), giving the positive address holonomy

```text
prod_s P_s>prod_s M_s.                               (27)
```

This product statement is canonical-gauge data; (12), (15), and (24) are
the gauge-invariant content.

## 5. The blocker cycle is a uniformly contracting affine local system

Equation (16) transports the centered rounding errors along every edge:

```text
epsilon_(s+1)=a_s epsilon_s+b_s,
a_s=d_(s+1)/(c+d_s),
b_s=delta_s-1/2-beta_s.                              (28)
```

After one circuit,

```text
epsilon_1=A epsilon_1+B,
A=prod_s a_s=prod_s d_s/(c+d_s),
B=sum_s b_s w_s,
w_s=prod_(r>s) a_r.                                  (29)
```

Every factor is positive and less than one.  Since `m>=2` and THM-1233
gives `d_s/(c+d_s)<2345/2346`,

```text
0<A<(2345/2346)^2,
1-A>1-(2345/2346)^2=4691/5503716.                    (30)
```

Thus every fixed choice of slopes and beta/delta data has a unique affine
fixed point.  For the actual spoke errors, `|epsilon_1|<=1/2`, so

```text
|B|<=(1-A)/2.                                        (31)
```

Eliminating the strict blocker errors from (29)--(31) yields the finite-word
exclusion target

```text
|sum_s (delta_s-1/2)w_s|
 <(1/14)sum_s w_s+(1-A)/2.                           (32)
```

The sampled cycle is therefore neither an unbounded address walk nor a free
combinatorial cycle.  It is a length-at-most-six word over (18), carried by
a uniformly contracting affine error system and the exact gcd residues
(13).

## 6. Every cycle exports a phase-located gcd seam

Every directed cycle (23) has a speed-ascent edge `i->j` with `d_i<d_j`;
otherwise speed would decrease strictly around a cycle.  At `t_i`, runner
`j` lies strictly inside one danger tooth `J` by construction.  At `t_j`,
runner `j` has depth greater than `1/4`, so the segment from `t_i` to `t_j`
leaves `J`.  Let `w` be its first boundary point.

The source remains quantitatively safe there.  Since the width of `J` is
`1/(7d_j)`,

```text
|w-t_i|<1/(7d_j),
||d_iw||>1/4-d_i/(7d_j)>3/28>1/14.                   (33)
```

This wall event is in fact protected.  Let `S_i` be the closed `d_i`-safe
component through `t_i` and put

```text
K_i=G intersect S_i.                                 (34)
```

If `Delta_i=||d_it_i||`, the distance from `t_i` to either endpoint of
`S_i` is at least

```text
rho_i=(Delta_i-1/14)/d_i>5/(28d_i).                  (35)
```

The centered-spoke calculation used in THM-1244 shows that the distance to
either endpoint of `G` is still larger than `rho_i`.  Every point of the
whole target tooth `J` is less than its width from `t_i`, and

```text
1/(7d_j)<1/(7d_i)<5/(28d_i)<rho_i.                   (36)
```

Hence the closure of `J` lies in `int(K_i)`.

At `w`, the carrier `c`, source `d_i`, and target `d_j` are all safe (the
target is equality-safe).  The assumed six-cover therefore forces a third
fast label

```text
h notin {i,j}                                         (37)
```

to be strictly dangerous at `w`; equivalently, `h` is a fourth support when
the carrier is counted.  Its open tooth contains `w`, so it overlaps `J` on
a nonempty open interval `Omega_ijh` adjacent to `w`.  By (36), this overlap
lies wholly in `int(K_i)`.  Its endpoints are rational tooth endpoints.  A
positive endpoint difference has an integer numerator divisible by
`gcd(d_j,d_h)`, and therefore

```text
mu(Omega_ijh)>=gcd(d_j,d_h)/(14d_jd_h).               (38)
```

Thus every actual blocker cycle exports at least one phase-located,
gcd-quantized seam.  This is the coverage-sensitive information absent from
the six sampled phases: the functional cycle has become a wall-event
hypergraph with a genuine metric credit.

## 7. The slowest orbit anchors one positioned Hunter credit

Apply the same argument to the selected outgoing edge `d1->d_j`.  It is
automatically a speed ascent.  Here `K_1` in (34) is exactly THM-1244's
protected slowest-spoke component `K`, so (38) is a located overlap between
two of the five active blocker labels inside `K`.

THM-1244 supplies a rank-two forest of two distinct projected handoff edges
inside `K`.  At least one of those edges differs from the wall edge
`{d_j,d_h}`.  Adjoin that edge to `{d_j,d_h}`.  Two distinct edges form a
forest, so THM-1244's forest can be reselected with one prescribed credit:

```text
one of its two positioned Hunter edges is the wall seam (38).             (39)
```

The quantitative scalar and two-unlocated-edge debts of THM-1244 remain
valid, since each selected overlap still contributes at least
`gcd(d2,...,d6)/(14d6^2)`.

Starting the functional orbit at `d1` connects this anchored seam by at most
six sampled edges to its eventual cycle, and (15)--(18) encode the entire
tail.  The initial edge need not itself lie on that cycle.  Conversely, the
cycle seam from Section 6 need not lie in the particular slowest component
`K`.  Equation (39) is therefore a genuine sampled-orbit-to-positioned-debt
bridge, but not yet the missing transport along the tail.

THM-1250 now sharpens (39).  Its connected chronological overlap graph on all
six private-provider labels, together with the prescribed wall edge (38),
contains a five-edge spanning tree including that edge.  Every edge is an
actual interior tooth overlap with quantum `1/[14 lcm(u,v)]`, and repeated
handoffs are charged by its multiplicity-averaged Hunter debt.  THM-1248 was
proved independently of that later statement, but their composition is now
available.  What remains is not tree extraction: it is oriented transport
from this prescribed seam through the affine clocks `d_i -> c+d_i`.

## 8. What address data remains

The quotient has removed the unbounded tooth lift, not the actual speed
packet or its torsion.  Within the sampled centered subsystem, the remaining
unbounded phase-address coordinate can be written

```text
E_i=2c epsilon_i
   =2cP_i-(2k+1)Q_i in [-c,c] intersect Z,
E_i==-(2k+1)Q_i                         (mod 2c).      (40)
```

Thus all spokes share one odd multiplier `2k+1` on the modulus `2c`.  The
finite `delta` word, the contracting errors, the exact `gcd(Q_i,Q_j)` sheet,
the wall seams, and this common residue orbit are the complete retained
sampled object.  Absolute scale and gcd/torsion remain deliberately visible
sidecars; no claim that all arithmetic is finite is intended.

## 9. Two exact guardrails

The finite relative bank cannot be replaced by a bound on absolute tooth
addresses.  For every `c>=27`, set `k=c-1` and take the fast packet

```text
{c+1,c+2,c+3,c+4,2c,3c}.                             (41)
```

At the six centered spokes, the first four labels may point to `2c` and
`2c<->3c` is an exact blocker two-cycle with relative word `(0,1)`, while
the absolute tooth addresses grow linearly with `c`.  The full seven-speed
packet is primitive.  Nevertheless

```text
t*=1-2/(5c),             min_v ||vt*||=1/5.           (42)
```

These packets are globally lonely, not six-comb covers.  They prove that
primitivity, compact ratios, sampled blocker obligations, and cycle
contraction alone cannot control the absolute lift.

Nor can every ascent digit be forced binary.  At `c=1,k=0`, the packet

```text
(1;2,16,17,34,35,2343)                               (43)
```

satisfies all THM-1233 adjacent-ratio bounds and splits into the centered
two-cycles

```text
2<->2343,              16<->17,              34<->35. (44)
```

On `2->2343`,

```text
(Q_i,P_i,N_i,M_i,delta_i,ell_i)=(3,2,1562,1562,-390,2), (45)
```

while the three cycle words are `(-390,1),(1,0),(1,0)`.  Yet `t=1/6` gives
the seven-speed packet depth `1/6`.  This is again a globally lonely
guardrail, not a cover counterexample.  The full alphabet (18), or a
coverage-sensitive replacement for it, is genuinely necessary.

## 10. Tournament and alternate-vertex audit

For the runner tournament, use the pairwise observable

```text
D_ij=P_iQ_j-P_jQ_i=Q_iQ_j(t_i-t_j),                  (46)
```

orient by increasing lexicographic `(t_i,i)`, and use the runner index to
break `D_ij=0` ties.  The tie Hamiltonian path is the same sorted order.  It
is transitive, with score histogram `(0,1,2,3,4,5)`, no directed cycles,
six singleton SCCs, no edge flips after the gauge is fixed, and one
Hamiltonian path.  It loses the central remainder, relative digit, blocker
truth, gcd sheet, and wall ownership.

The proof-bearing directed relation is instead

```text
i -> j  iff d_j is dangerous at the selected centered phase t_i,          (47)
```

decorated by (14).  Its cycle is then lifted to the boundary-obligation
hyperedge `(i,j,h,w,Omega_ijh)` from Section 6.  We challenged runners,
gaps, fixed sections, tooth addresses, phase determinants, residues, gcd
sheets, wall crossings, safe components, overlap seams, private owners, and
proof obligations as vertices.  The smallest faithful carrier found here is

```text
(finite delta word; contracting epsilon transport; positive ell residues;
 gcd sheets; common odd residue orbit; located wall-event hyperedges).     (48)
```

It preserves the exact predicates used in the proof and explicitly destroys
the absolute tooth lift and off-sample cover ownership.  Those destroyed
off-sample incidences, not another orientation of runners, are where the
remaining theorem must act.

## 11. Verification and scope

The dependency-free exact referee replays (7), (9), (12)--(18), (21)--(32),
the common-clock divisibility, positive address holonomy, both guardrail
families, all stated constants, and the tournament methodology using only
integers and exact rational arithmetic.  Normal and optimized outputs are
byte-identical.  It does not computationally certify the continuum wall
argument in Section 6; (33), (35), and (36) are exact paper inequalities.

The Lean module kernel-checks the lower-gap cocycle and translation gauge,
determinant and gcd identities, recentered blocker identity, central band,
global and binary digit bounds, binary phase separation, a three-edge affine
transport, the contraction constant, a three-edge positive product
holonomy, and the wall margin.  Arbitrary-cycle telescoping, nearest-integer
selection, strict danger/safety interpretation, interval continuation, and
the gcd-quantized overlap seam remain explicit paper providers.  There are no
proof placeholders or `native_decide` calls.

Frozen hashes are

```text
source         d236f1ec6b48a95e7dfc3a7f3beb5a68061aa41a61a2eb01e3b1ba38827af704
output         ba1b6323d40c34d0b204d75f4658f67fa9d053ced7a9d2b45c6cf32298f2f8bb
formalization  54487de2f3886570b2752f044457d2d174edc1b2c70a3067eac397254b06ae35
```

THM-1248 proves uniform relative-address compression and a real
cycle-to-wall seam.  It does not prove six-comb noncoverage or LRC(14).  The
principal live target is no longer to bound blocker tooth addresses.  It is
to transport the slowest-orbit seam through the finite tail to the eventual
positive-holonomy cycle, or to exploit a fully located tree plus oriented
private-owner germs to force an off-sample escape/well-founded alternate-gap
descent.  THM-1247 shows that contracted Fano incidence without this oriented
metric continuation is insufficient.
