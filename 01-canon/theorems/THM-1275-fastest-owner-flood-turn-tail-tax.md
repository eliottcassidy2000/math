---
id: THM-1275
title: FASTEST-OWNER FLOOD/TURN TAIL TAX -- every long fastest subsequence pays whole skipped teeth or multi-low return seams
status: PROVED (global deletion-minimal-word theorem; exact eligible-owner regular-run bound; functional private-occurrence refinement; whole-tooth flood invoice added pointwise to the entire chronological seam invoice; exact all-return packet form; weighted and unweighted tail taxes; 7/2-dominated corollary; c=140 six-spoke nonadditivity countermodel; optimization-safe exact referee; sorry-free Lean arithmetic core).  The new lower bound can remain scale-small and does not prove six-comb noncoverage or LRC(14)
source: codex-2026-07-19 multispoke tail-tax audit
depends_on: [THM-1198, THM-1253, THM-1264, THM-1266]
related: [THM-1199, THM-1256, THM-1267]
script: 04-computation/lrc14_fastest_owner_flood_turn_tail_tax_thm1275.py
output: 05-knowledge/results/lrc14_fastest_owner_flood_turn_tail_tax_thm1275.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCFastestFloodTurnTailTax.lean
script_sha256: 622028f9c849e96d221a7bb9a1840a5b66ee458f5c9fdadf7348d551eb56bada
output_sha256: 3eedb36da439609e6e9d24781452de44c6132dce38a607e5d76ecb838f750ea0
formalization_sha256: 56d4fc5aedb1f25605e2ea53d4acece53a85200909b32535c7ceffc9cec31eb0
---

# THM-1275 -- fastest-owner flood/turn tail tax

## 1. The fastest-owner subsequence and its private count

Let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]
```

be a complete `c`-safe gap covered by the strict danger combs of

```text
c<d_1<...<d_5<h,                    h=d_6.             (1)
```

Choose THM-1253's deletion-minimal cover by individual open teeth and order
it chronologically.  Let `K` be the number of selected teeth owned by `h`.
THM-1198's private-region argument already gives

```text
K>=ceil(h/(7c)).                                         (2)
```

There is a useful functional version.  With `Pbar` the exact/BV one-comb
majorant from THM-1198, put

```text
eta_h=1-sum_(i=1)^5 Pbar(6d_i/(7c))>=1/36.              (3)
```

The private set

```text
Q_h=G minus union_(i=1)^5 D_(d_i)                       (4)
```

has six-bin dual mass at least `eta_h` and, under full coverage, is contained
in the `h` comb.  One normalized `h` tooth has length `c/(6h)`, so its dual
mass is at most

```text
(7/6)c/(6h)=7c/(36h).                                   (5)
```

Every `h` tooth meeting positive private mass is indispensable in every tooth
subcover.  Therefore

```text
K>=ceil(36h eta_h/(7c))>=ceil(h/(7c)).                  (6)
```

This is a private-mass count of literal teeth, not a count of abstract runner
labels.

## 2. Regular rungs, floods, and turns

List the selected `h` teeth in chronological order as

```text
H_1,...,H_K,
address(H_r)=M_r,                     M_1<...<M_K.      (7)
```

Between `H_r` and `H_(r+1)` there is a nonempty literal word of low teeth;
write its length as `s_r` and put

```text
R_r=M_(r+1)-M_r>=1.                                    (8)
```

Call transition `r` **regular** when

```text
R_r=1,                         s_r=1.                   (9)
```

Thus it is one consecutive-address toothpick rung `H,J,H`.  Define

```text
e=#{i in {1,...,5}: h>6d_i}.                           (10)
```

THM-1256 says the low owner of every regular rung belongs to this eligible
set.  THM-1266 says two occurrences of one low owner in a regular star have
slot distance at least six.  Since `e<=5`, a run of `e+1` regular transitions
would repeat one of its `e` eligible labels at distance at most `e<=5`.
Hence

> **Eligible-owner run law.** Every contiguous run of regular fastest-owner
> transitions has length at most `e`.

The two nonregular mechanisms are kept separate:

```text
X=sum_(r=1)^(K-1)(R_r-1),                              (11)
T=#{r:R_r=1 and s_r>=2}.                               (12)
```

Here `X` is the number of **skipped fastest addresses**, not merely the number
of transitions that skip.  A transition is nonregular iff it contributes a
positive amount to `X` or is counted by `T`.  If `E` is the number of
nonregular transitions, the regular transitions split into at most `E+1`
runs.  Therefore

```text
K-1-E<=e(E+1),
E>=ceil(K/(e+1))-1,
X+T>=E>=ceil(K/(e+1))-1.                              (13)
```

The ceiling is sharp for the abstract skeleton: concatenate runs of `e`
regular transitions separated by one exception.  At `e=5,K=6`, (13) permits
the exact five-rung star and gives zero forced tail, as it must.

## 3. A skipped address is a whole covered fastest tooth

Suppose consecutive selected fastest teeth have addresses `M` and `M+R` with
`R>=2`.  Deletion minimality gives a witness `x in H_M intersect G` and a
witness `y in H_(M+R) intersect G`.  For every `1<=q<R`, tooth endpoint order
gives

```text
x<right(H_M)<left(H_(M+q))
 <right(H_(M+q))<left(H_(M+R))<y.                    (14)
```

Since `G` is an interval, the **whole** open tooth `H_(M+q)` lies in
`int(G)`.  It was not selected, while the selected low teeth still cover all
of `G`.  Consequently the full six-comb multiplicity is at least two
throughout that skipped tooth.

The skipped fastest teeth are pairwise disjoint.  In normalized slow-gap
coordinates each has length `c/(6h)`.  The six-bin density has floor `3/4`,
so every skipped address pays weighted multiplicity excess

```text
(7c/6)(3/4)(1/(7h))=c/(8h).                           (15)
```

This is the **flood invoice**.  It detects the branch in which a broad low
tooth, or a low-tooth flood, makes an intermediate fastest tooth redundant;
that branch is invisible in a naked `H,J,H` star word.

## 4. Return packets and the two-layer multiplicity invoice

For any successive-fastest transition, its literal word is

```text
H_M,J_1,...,J_s,H_(M+R),                 s>=1.         (16)
```

and write `j_q` for the speed of `J_q`.  THM-1264's exact return polygon
identity gives the total raw seam mass of this packet:

```text
Omega_r
 =sum_(packet seams) overlap
 =(1/7)(sum_(q=1)^s 1/j_q-(7R-1)/h)>0.               (17)
```

For a transition counted by `T`, `R=1,s>=2`, and (17) specializes to

```text
Omega_r=(1/7)(sum_(q=1)^s1/j_q-6/h).                  (18)
```

THM-1253 makes the full family of raw chronological seams pairwise disjoint.
The skipped fastest teeth are also pairwise disjoint.  The two families need
not be disjoint **from each other**, but they may still be added.  Pointwise:

```text
1_(skipped h tooth)+1_(raw seam) <= C(t)-1.           (19)
```

Indeed, away from an intersection each indicator already has multiplicity at
least two.  On an intersection, the raw seam belongs to two distinct selected
teeth.  Neither can be an `h` tooth, because distinct teeth of the `h` comb
are disjoint and the skipped one was not selected.  The full skipped `h`
tooth is therefore a third active owner, so `C(t)>=3` and both units are paid.

This is stronger than proving geometric cross-disjointness: it adds the flood
invoice to **every** seam in THM-1253, including seams in regular packets and
in nonconsecutive returns.  Normalization and `f>=3/4` give

```text
(7c/6)(3/4)|W_a|=(7c/8)|W_a|,
(7c/8)[1/(14 lcm(s_a,s_(a+1)))]
 =c/[16 lcm(s_a,s_(a+1))].                           (20)
```

On the packet (16), (17) contributes the exact lower expression

```text
(7c/8)Omega_r
 =(c/8)(sum_q1/j_q-(7R-1)/h).                        (21)
```

The prefix before the first selected `h` tooth and suffix after the last may
contain further paid seams.  Dropping them gives a safe packet-only form.

## 5. The global functional and harmonic tail taxes

Let

```text
F_6=sum_(i=1)^6 Pbar(6d_i/(7c))-1.                    (22)
```

The actual six-bin weighted coverage excess is at most `F_6`.  The pointwise
two-layer law (19) proves the additive upgrade

```text
F_6 >= cX/(8h)+(7c/8)sum_(a=1)^(N-1)|W_a|
    >= cX/(8h)+(c/16)sum_(a=1)^(N-1)
                          1/lcm(s_a,s_(a+1)).          (23)
```

Using (17) on every successive-fastest packet and dropping prefix/suffix
seams gives the exact return form

```text
F_6 >=(c/8)[X/h
       +sum_(r=1)^(K-1){sum_(q=1)^(s_r)1/j_(r,q)
                         -(7R_r-1)/h}].               (24)
```

The unweighted singleton-discrepancy upper bound used in THM-1253 is

```text
physical coverage excess <=(6/49)(sum_i1/d_i-1/c).    (25)
```

Applying the same two-layer lower invoice before (25) gives the full seam
upgrade

```text
sum_(i=1)^6 1/d_i-1/c
 >=7X/(6h)+(49/6)sum_a|W_a|
 >=7X/(6h)+(7/12)sum_a1/lcm(s_a,s_(a+1)),             (26)
```

and the packet form

```text
sum_i1/d_i-1/c
 >=(7/6)[X/h
       +sum_(r=1)^(K-1){sum_q1/j_(r,q)-(7R_r-1)/h}].   (27)
```

Every packet bracket in (24) and (27) is strictly positive.  The new term
`X/h` is additional to the complete THM-1253 seam invoice; it is not a second
copy of any seam.  Individual packet brackets can be arithmetically small,
so the theorem retains their exact functional form instead of replacing them
by a false scale-free constant.

There is a clean top-dominated branch.  Every turn packet has at least two
internal speeds, each at most `d_5`.  If

```text
h>=(7/2)d_5,                                          (28)
```

then

```text
sum_q1/j_q-6/h >=2/d_5-6/h>=1/h.                     (29)
```

For every skipped transition, the flood term contributes at least `1/h`; for
every consecutive-address turn, (27) does the same.  Discarding the positive
regular and nonconsecutive packet seams, then combining (6), (13), and
(23)--(29), yields

```text
F_6 >= c/(8h)[ceil(K/(e+1))-1]
    >= c/(8h)[ceil(36h eta_h/(7c(e+1)))-1]
    >= c/(8h)[ceil(h/(7c(e+1)))-1],                  (30)

sum_i1/d_i-1/c
 >=7/(6h)[ceil(K/(e+1))-1].                          (31)
```

For `e=5`, the basic extra-tail bound first becomes nonzero at `h>42c`.  It
approaches only a small constant after multiplication by the functional
scale, so (30) is a genuine global completion tax but not by itself a
contradiction.

## 6. Why the six centered spokes do not simply add

The primitive sharp row from THM-1266 is also an exact warning against the
most tempting alternative coupling.  At

```text
c=140, k=80,
(d_1,...,d_6)=(254,255,256,257,292,1805),             (32)
```

the six THM-1240 nearest-integer errors have absolute values

```text
(9/20,1/8,3/10,11/40,2/5,3/8).                       (33)
```

Substitution in THM-1267's normalized protrusion expression

```text
ell_i=max(0,1/2+7rho_i/6-d_i/(2c))                   (34)
```

gives the exact raw values

```text
d=254:   33/280,
d=255:  -89/336,
d=256:   -9/140,
d=257: -163/1680,
d=292:   -8/105,
d=1805: -617/112.                                    (35)
```

Thus only the `254` component protrudes, by `33/280`; the other five are
contained in `G`.  The same row realizes the two blocker cycles and aligned
third-owner bridges of THM-1266.  It is not a six-cover--its full comb sweep
has four holes--so it does not refute a theorem using additional global
coverage.  It does refute any attempt to add six positive protrusion taxes
from centered-spoke geometry and the local blocker package alone.

Its displayed fastest chain has exactly `K=6`, `e=5`, `X=T=0`: the sharp
five-rung star is precisely the zero-**extra-tail** boundary of (13).  Its ten
ordinary seams remain paid in (23), as they should.  Completing the row with
a seventh indispensable fastest occurrence must enter a flood or a turn,
which is the new global information recorded by (30).

## 7. Tournament and carrier audit

Runner speed order remains a transitive tournament and destroys the new
predicate.  The useful vertices are the `K-1` **successive fastest-owner
obligations**.  The binary observable is

```text
regular  versus  flood/turn,                          (36)
```

with fastest-address order as the switch/gauge and chronological order as
the tie Hamiltonian path.  This carrier preserves address jumps, internal-low
counts, and the two-layer multiplicity invoice.  It destroys the phases and exact
low-low clocks inside a turn, so the faithful state retains

```text
(M_r,R_r; internal low tooth word; raw seam intervals).          (37)
```

We challenged runners, gaps, fixed circle sections, high safe components,
high tooth addresses, wall crossings, low flood blocks, return stalks, and
proof obligations as vertices.  The challenged assumption that matters is:
indispensable fastest teeth need not form one consecutive-address star.
Address skips are not a failure of the toothpick method; they are whole-tooth
flood invoices.

## 8. Verification and scope

The dependency-free exact referee exhausts `185,136` admissible
regular/exception words and verifies the sharp ceiling on `96` parameter
rows.  It checks `21,036` eligible-colour words, `15,680` skipped-tooth
endpoint rows, `4,480` turn-corridor disjointness rows, the literal return

```text
H_26,J_5,J_12,H_26,       Omega=41/5460,              (38)
```

`95,950` normalization/private-count rows, `13,135,589` dominated-turn
rows, and every rational in (15), (18), and (29)--(35).  Normal and optimized
Python outputs are byte-identical.

The sorry-free Lean module kernel-checks the regular-run/exception product
bound, its flood/turn consumer, the private-mass product inequality, the
two-layer pointwise multiplicity rule, all normalization coefficients, the
general packet coefficient, the `7/2` simplification, the literal return
constant, and all six `c=140` protrusion values.  Minimal-subcover extraction,
the interval placement and internal-disjointness arguments, and integration
remain explicit paper providers.

Frozen hashes are

```text
source         622028f9c849e96d221a7bb9a1840a5b66ee458f5c9fdadf7348d551eb56bada
output         3eedb36da439609e6e9d24781452de44c6132dce38a607e5d76ecb838f750ea0
formalization  56d4fc5aedb1f25605e2ea53d4acece53a85200909b32535c7ceffc9cec31eb0
```

THM-1275 closes the local-to-global logical gap in the sharp five-rung star:
after five regular rungs, continuation is forced into a quantitatively paid
flood or turn.  More strongly, every whole-tooth flood payment adds pointwise
to the already complete chronological seam invoice.  The resulting extra
lower bound may still decay with scale and does not prove universal six-comb
noncoverage, the empty sporadic branch, or LRC(14).
