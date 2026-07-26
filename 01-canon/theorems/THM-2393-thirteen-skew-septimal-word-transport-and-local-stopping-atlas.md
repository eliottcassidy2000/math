---
id: THM-2393
title: "Thirteen-skew septimal word transport and local stopping atlas"
status: >
  PROVED STRUCTURAL + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.
  In THM-2367's last k=2,(t,b)=(1,0) septimal lane, THM-2388's
  quotient-blocker cage saturates on every generic C_3-safe top bin:
  the two divided low-blocker slices partition the guard slice
  disjointly and every lower q slice lies inside that guard slice.
  Seven-root transfer forces nu_7(H)=nu_7(C_1)=nu_7(C_2); combined with
  THM-2390, this first leaves weight seven or eight. The independently
  audited THM-2391 full-bin period theorem then forces all four lower q
  labels into that layer and excludes weight seven. Multiplication by
  thirteen does not identify the divided and physical blocker words on
  one fibre: with a reduced base z={13y}, it sends W_(13A)(y) to the
  affine reflection floor(13y)-W_A(z) modulo seven. An exact strict W=8
  packet is a genuine local stopping atlas. An exact W=7 word remains a
  useful hostile to one-fibre gluing, but two displayed points elsewhere
  in its 49-point top bin violate respectively the blocker partition and
  lower-q containment required by THM-2391. An infinite triple-shield
  family separately refutes a bare joint-comb anti-shield. The remaining
  obligation is global labelled 13-skew transport, not a same-fibre
  capacity contradiction. No septimal branch, thirteen-adic row, target,
  LRC(14), or ledger entry is closed.
source: codex-2026-07-26-thirteen-skew-septimal-transport
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
related:
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2385-two-top-septimal-blocker-collision-reduction
script: 04-computation/lrc14_thirteen_skew_septimal_word_transport_thm2393.py
output: 05-knowledge/results/lrc14_thirteen_skew_septimal_word_transport_thm2393.out
script_sha256: 2ee29f92566cc2144e3971a3791bd8c18ab1faf3a10ff4cea55f76d5ead8f3ac
output_sha256: 3a3e6ee05c6f6090ec12f96d7eddb9f37be9567a64ab396c8ab9007ca2671d1e
hash_basis: working-tree bytes (LF)
---

# THM-2393 -- thirteen-skew septimal word transport

**PROVED STRUCTURAL + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT
AUDIT.**

THM-2388 and THM-2390 expose two apparently incompatible finite words in
the sole remaining septimal lane:

```text
quotient-blocker word:
  two divided blockers partition the guard on a top-q bin;

terminal heavy word:
  the physical blockers participate in a seven-root partition or
  one-double cover.
```

They do not live at the same scale. If `c_j=13C_j`, multiplication by
thirteen reflects the root labels and also contributes the integer digit
discarded when `13y` is reduced modulo one. Retaining that digit produces
the honest chain

```text
capacity-saturated quotient split
  -> common guard/divided-blocker septimal depth
  -> THM-2390 local words W=7 or W=8
  -> THM-2391 full-bin period forces W=8
  -> affine thirteen-reflection from C_j to c_j
  -> one nonempty live local cross-time atlas
  -> global labelled transport remains open.                       (1)
```

These local packets are stopping controls, not scalar covers.

## 1. Last-lane notation

On `T=R/Z`, put

```text
D_v={x:||vx||<1/14},

E_H={x:||Hx||<1/7}.                                  (2)
```

Retain the primitive scalar cover and the only-`c_3`-dominant branch

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))>0,

nu_7(q_*)=M,

nu_7(H),nu_7(q_i),nu_7(c_1),nu_7(c_2)<M
                                      for q_i!=q_*,

nu_7(c_3)>M,                                         (3)

k=2,                         (t,b)=(1,0).
```

The six labels `H,q_1,...,q_5` are units modulo thirteen and

```text
c_j=13C_j,                         j=1,2,3.           (4)
```

Division by thirteen preserves every septimal valuation. Define

```text
K=1_(E_H)+sum_(i=1)^5 1_(D_(q_i)),

U=1_(E_H)+sum_(q_i!=q_*)1_(D_(q_i)),

B=D_(C_1) union D_(C_2) union D_(C_3).               (5)
```

THM-2388 proves the directed collision cage

```text
{K>=2} subset B                                      (6)
```

almost everywhere. The blockers in (6) are the **divided** blockers at
the physical collision point.

## 2. Capacity saturation forces the quotient split

On `D_(q_*)`, one has

```text
K=1+U.
```

Therefore (6) gives

```text
D_(q_*) intersection {U>0} intersection D_(C_3)^c
 subset D_(C_1) union D_(C_2).                       (7)
```

Put

```text
N=7^(M+1)
```

and take a generic additive `N`-orbit on which `C_3` is safe. The top
word `D_(q_*)` occupies one complete bin of `N/7` points. Inside that
bin:

```text
E_H has                     2N/49 points,

D_(C_1),D_(C_2) have        N/49 points each.        (8)
```

The left side of (7) includes the guard slice. Its size is already the
sum of the two blocker capacities. Consequently every inclusion
saturates:

```text
D_(q_*) intersection D_(C_3)^c intersection E_H

 =

D_(q_*) intersection D_(C_3)^c intersection
  (D_(C_1) disjoint_union D_(C_2))                  (9)
```

almost everywhere. Likewise each lower unit word is caged by that same
guard slice:

```text
D_(q_*) intersection D_(C_3)^c intersection D_(q_i)
 subset E_H,                         q_i!=q_*.        (10)
```

The inherited cover-null set, its finitely many orbit translates, and all
strict endpoints are discarded before the word count. They remain null.
Equation (9) is pointwise on every remaining top bin, not only an
integrated identity.

In particular,

```text
D_(q_*) intersection D_(C_1) intersection D_(C_2)
 subset closure(D_(C_3)).                            (11)
```

Section 7 shows why (11) alone is much weaker than (9).

## 3. The guard and divided blockers have one septimal depth

Put

```text
h=nu_7(H),        r_j=nu_7(C_j),             j=1,2,

e=min(h,r_1,r_2).                                   (12)
```

Since `e<M<nu_7(C_3)`, the top and high masks remain invariant on every
seven-root fibre used below.

For each layer strictly below `e`, all three masks in (9) are divisible
by seven. Sum (9) over the seven roots of multiplication by seven. The
carrier bits are constant, every term contributes seven times its
divided mask, and the same equality descends one layer. Iterate to depth
`e`.

At the first terminating layer, choose a generic base in

```text
D_(q_*/7^(e+1)) minus closure(D_(C_3/7^(e+1))).     (13)
```

This set has positive measure. Both combs have measure `1/7`; if their
difference were null, they would agree almost everywhere, which is
impossible because their frequencies have distinct septimal depths.

On the seven roots over (13), a terminating ordinary blocker contributes
one incidence, a terminating guard contributes two, and an
unterminated mask contributes a multiple of seven. Reducing the root-sum
identity modulo seven gives

```text
2*1_(h=e)
 =
1_(r_1=e)+1_(r_2=e)                    in F_7.       (14)
```

At least one of the three depths equals `e`. The two sides of (14) lie
in `{0,2}` and `{0,1,2}`, respectively. Equality is possible only when
all three terminate:

```text
nu_7(H)=nu_7(C_1)=nu_7(C_2)=e<M.                    (15)
```

Since multiplication by thirteen does not change septimal depth,

```text
nu_7(c_1)=nu_7(c_2)=e.                              (16)
```

This argument uses only the saturated equality (9). THM-2391, now
independently hostile-audited, uses the individual containments (10) on
the whole lifted top bin to obtain the stronger conclusion

```text
nu_7(H)=nu_7(q_i)=nu_7(C_1)=nu_7(C_2)=0,
                                      q_i!=q_*.      (16a)
```

That full-bin conclusion, rather than the one-fibre residue gate (14),
is what removes the weight-seven alternative below.

## 4. Synthesis with the THM-2390 heavy word

THM-2390 says every survivor has a lower layer of weight seven or eight.
That layer contains the guard and at least five of the six lower ordinary
labels

```text
{q_i:q_i!=q_*} union {c_1,c_2}.                     (17)
```

The layer containing the guard is the common layer `e` in (15)--(16), so
both physical blockers already contribute. Its fixed weight is

```text
guard + c_1+c_2 =2+1+1=4.                           (18)
```

Therefore at least three of the four lower `q_i` also have depth `e`.
Before the full-bin theorem is imposed, there are only two local
possibilities:

```text
W=7: guard, c_1,c_2, and exactly three lower q labels;

W=8: guard, c_1,c_2, and all four lower q labels.   (19)
```

On THM-2390's positive family of generic higher-safe fibres:

- the `W=7` word is an exact partition of `F_7`;
- the `W=8` word covers every root once except for one root covered
  twice.

THM-2391 uses the whole lifted top bin in (10), rather than only (14),
to force every lower `q_i` into the common layer and then uses
primitivity to make `e=0`. Hence `W=7` is globally excluded and every
survivor has the `W=8` word. The `W=7` packet in Section 6 is retained
only as a one-fibre hostile to a naive local permutation proof; its
explicit failure of the full-bin hypotheses is audited there. The
`W=8` atlas is the actual local stopping boundary.

## 5. The exact thirteen-skew word law

Scale the common depth `e` out of the three labels in (15). For a
seven-unit integer `v`, width `a in {1,2}`, and a real base `y`, define

```text
W_(v,a)(y)
 ={k in F_7:||v(y+k)/7||<a/14}.                     (20)
```

Suppose

```text
C_j=7^e A_j,                  c_j=7^e(13A_j).
```

For the unreduced real lift `13y`,

```text
W_(13A_j,1)(y)=-W_(A_j,1)(13y).                     (21)
```

Indeed, put `l==13k==-k mod 7`. Then

```text
13A_j(y+k)/7
 =A_j(13y+l)/7                         mod Z,
```

so the physical root labelled `k` is the divided root labelled `l=-k`.

There is an additional digit when both bases are reduced to `[0,1)`.
Write

```text
z={13y},                  m=floor(13y).
```

Changing the base from `13y=z+m` to `z` translates the root labels.
The reduced-base form of (21) is

```text
W_(13A_j,1)(y)
 =m-W_(A_j,1)(z)                         in F_7.     (22)
```

The integer `m mod 7` is not cosmetic. It is the root digit lost by
reducing the base before comparing the two words.

The high-state typing aligns under the same map. Put

```text
Q=q_*/7^(e+1),                 A_3=C_3/7^(e+1).
```

On a terminal seven-root fibre:

```text
q_* is safe at y       iff y notin D_Q,

q_* is dangerous at z  iff z in D_Q,

c_3 is safe at y       iff C_3 is safe at z,        (23)
```

where the last equivalence follows from

```text
(c_3/7^(e+1))y=13A_3y=A_3z                   mod 1.
```

Thus the honest comparison is cross-time:

```text
terminal physical-blocker word at y

versus

quotient-blocker split at z={13y},                  (24)
```

with affine digit `m`. They are not two descriptions of one same-base
word.

## 6. The live atlas and an excluded one-fibre hostile

The following exact packets are not scalar covers. The first is the
live local stopping atlas. The second proves that a one-fibre
permutation argument cannot establish THM-2391's full-bin conclusion:
the local word exists, but the required bin-wide identities fail away
from its displayed seven roots.

### 6.1. The actual `W=8` stopping atlas

Take

```text
y=73/1009,                         z=13y=949/1009,

H=1,                               q_*=7,

(q_2,q_3,q_4,q_5)=(37,82,67,22),

(C_1,C_2)=(15,1),                  (c_1,c_2)=(195,13),

C_3=637=13*7^2,                    c_3=8281=13^2*7^2. (25)
```

Here `m=floor(13y)=0`. The exact seven-root words are

```text
label             physical word at y     divided/same word at z

guard H, width 2        {0,6}                       {0,6}

q_2=37                 {2}                         {0}
q_3=82                 {3}                         {0}
q_4=67                 {4}                         {0}
q_5=22                 {5}                         {0}

c_1=13C_1             {0}             C_1:        {0}
c_2=13C_2             {1}             C_2:        {6}.          (26)
```

At `z`, the divided blockers partition the guard and every lower `q`
chooses one guard endpoint. At `y`, the physical heavy word has
multiplicity vector

```text
(2,1,1,1,1,1,1),                                   (27)
```

with the unique double at root zero.

All state inequalities are strict. With

```text
Delta(v;n)=14||vn||_(1009)-1009,
```

where positive means safety and negative means danger,

```text
Delta(q_*/7;73)=13,              Delta(q_*/7;949)=-169,

Delta(c_3/7;73)=4801,
Delta(C_3/7;949)=4801.                              (28)
```

For the word tests, use denominator `7*1009=7063` and signed endpoint
distance

```text
14||v(n+1009k)||_(7063)-a*7063.
```

The minimum absolute value over every entry in (26) is `840`, so no word
uses an aligned endpoint. The valuation roles are exactly

```text
nu_7(H,q_2,q_3,q_4,q_5,c_1,c_2,q_*,c_3)
 =(0,0,0,0,0,0,0,1,2),

nu_13(c_1,c_2,c_3)=(1,1,2).                         (29)
```

Thus this is a strict repeated-first local model of the last lane.

### 6.2. The excluded `W=7` one-fibre hostile

Take

```text
y=11/1009,                          z=13y=143/1009,

H=8,                                q_*=49,

(q_2,q_3,q_4)=(79,68,53),

q_5=56                              [off the heavy layer],

(C_1,C_2)=(4,8),                    (c_1,c_2)=(52,104),

C_3=4459=13*7^3,                    c_3=57967=13^2*7^3. (30)
```

Again `m=0`. The heavy words are

```text
label             physical word at y     divided/same word at z

guard H, width 2        {0,6}                       {5,6}

q_2=79                 {3}                         {5}
q_3=68                 {4}                         {5}
q_4=53                 {5}                         {5}

c_1=13C_1             {2}             C_1:        {5}
c_2=13C_2             {1}             C_2:        {6}.          (31)
```

The physical word at `y` partitions all seven roots exactly. At `z`,
the divided blockers partition the guard and all three heavy `q` words
choose one endpoint. The off-layer speed `q_5=56` is constant and safe
on both displayed seven-root fibres.

The strict state signs are

```text
Delta(q_*/7;11)=69,                Delta(q_*/7;143)=-897,

Delta(c_3/7;11)=2925,
Delta(C_3/7;143)=2925,

Delta(q_5/7;11)=223,               Delta(q_5/7;143)=881.        (32)
```

The minimum absolute word-endpoint distance is `161`. The septimal
roles are

```text
heavy layer e=0,

nu_7(q_5)=1,             nu_7(q_*)=2,       nu_7(c_3)=3,

nu_13(c_1,c_2,c_3)=(1,1,2).                         (33)
```

The seven roots displayed in (31) are precisely the addresses
`s=7k`, `k in F_7`, in the 49-point top bin

```text
x_s=(1001+1009s)/(49*1009),            0<=s<49.     (33a)
```

Indeed `x_(7k)=(z+k)/7`. Every `x_s` is `q_*=49`
dangerous and `C_3=4459` safe. The full-bin hypotheses nevertheless
fail strictly in two different ways. Define

```text
Gamma(v,a;s)
 =14||v(1001+1009s)||_(49441)-a*49441,              (33b)
```

where `a=2` for the guard and `a=1` otherwise. At `s=6`,

```text
Gamma(q_*,1;6)=-43953,       Gamma(H,2;6)=-896,

Gamma(C_1,1;6)=247653,       Gamma(C_2,1;6)=48545,

Gamma(C_3,1;6)=143325.                              (33c)
```

Thus the top and guard are dangerous while both divided low blockers
and the high blocker are safe, contradicting the partition inclusion
(7). At `s=13`,

```text
Gamma(q_*,1;13)=-43953,      Gamma(q_5,1;13)=-43169,

Gamma(H,2;13)=97986,

Gamma(C_1,1;13)=48993,       Gamma(C_2,1;13)=147427,

Gamma(C_3,1;13)=143325.                             (33d)
```

Here the top and off-layer `q_5=56` are dangerous while the guard and
all divided blockers are safe, contradicting (10). Hence this packet
is not compatible with the proved THM-2391 full-bin hypotheses and is
not a survivor. It is a sharp local hostile only: a correct proof of
the single-layer theorem must retain the whole 49-point bin rather than
glue one seven-root word.

## 7. A bare triple anti-shield is false

The intersection consequence (11) cannot by itself close the lane.

> **Infinite triple-shield family.** Let `m` be a seven-unit and let
> `B` be a seven-unit integer satisfying
>
> ```text
> 50<=B<=91.
> ```
>
> Then
>
> ```text
> D_m intersection D_(Bm) intersection D_(7m)
>   subset D_(49m).                                  (34)
> ```

**Proof.** Pull back by multiplication by `m`; it suffices to take
`m=1`. Since `D_1` is its central interval,

```text
D_1 intersection D_7={x:||x||<1/98}.                (35)
```

For `B<=91`, the nearest noncentral `B`-tooth begins at

```text
13/(14B)>=1/98.
```

Thus only the central `B`-tooth can meet (35). For `B>=50`, its radius
satisfies

```text
1/(14B)<=1/700<1/686,
```

so that central tooth lies strictly inside the central `D_49` tooth.
This proves (34). QED.

The geometric boundaries are sharp:

```text
B=48: central escape,

B=91: endpoint touch,

B=92: first noncentral penetration.                 (36)
```

The endpoint value `B=91` is not itself a seven-unit; it records the
sharp geometric boundary of the family.

This is a hostile to the **triple containment route**, not a hostile to
the full quotient split (9). The distinction is literal. For

```text
B=50,             H=1,             x=49/175,
```

one has

```text
x in D_50 intersection D_7 intersection D_49^c,

x notin D_1 union E_1.                               (37)
```

Thus this shield packet fails the XOR/partition identity. A bare
joint-comb anti-shield is false, while the capacity-saturated quotient
split remains strictly stronger.

## 8. Faithful switching view and the global obligation

On one quotient-split seven-root fibre, the guard supplies two adjacent
root addresses in its own gauge. `C_1,C_2` occupy those addresses
bijectively, and each same-layer lower `q` chooses an endpoint. Swapping
the two guard addresses changes every vertex sign simultaneously. The
gauge-invariant pair relation records whether two labels choose the same
endpoint.

This is a balanced two-class switching word, not a tournament on runners.
It forgets:

- the actual root digit `m` in (22);
- which endpoint belongs to which blocker;
- the labelled slopes of the singleton words;
- the base transition `y -> {13y}`; and
- the top/high/off-layer safety states.

An anchored endpoint and the affine digit restore the lost coordinate.
The W=8 packet proves that the resulting local switching class is
nonempty.

The honest remaining target is therefore:

> **Global thirteen-skew transport obligation.** Exclude a labelled
> family of centred comb words which realizes the quotient split (9)
> on almost every `C_3`-safe top bin, the physical `W=8` one-double word
> on THM-2390's positive higher-safe carrier, and the affine transport
> (22)--(23) between them.

The local packets show that this cannot be replaced by a one-shot
`13 x 7` capacity or permutation check. A proof must use cross-base
coherence: endpoint current, address monodromy, a full lifted-bin period,
or an equivalent global sidecar.

No branch is excluded here. The thirteen-adic ledger remains `165`, and
LRC(14) remains open.

## 9. Exact companion

The dependency-free exact companion:

- checks `207,854` ordinary and guard instances of the affine
  thirteen-skew word law, including every reduced-base lift digit;
- checks the capacity saturation and all `42` labelled two-endpoint
  quotient splits;
- exhausts `216` depth triples and retains exactly the six common-depth
  profiles;
- reconstructs the `W=7/W=8` heavy-role counts after the common-depth
  reduction;
- verifies every word, valuation, state sign, and strict endpoint
  clearance in (25)--(33);
- embeds the excluded `W=7` fibre in its full 49-point top bin, finds
  exactly six blocker-partition failures and five `q_5`-containment
  failures, and checks the strict controls (33c)--(33d);
- proves the `36` admissible members of the shield family and checks the
  sharp `48/91/92` boundaries; and
- verifies explicitly that shielding does not imply the XOR split.

Run

```bash
python3 04-computation/lrc14_thirteen_skew_septimal_word_transport_thm2393.py
python3 -O 04-computation/lrc14_thirteen_skew_septimal_word_transport_thm2393.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_thirteen_skew_septimal_word_transport_thm2393.out
```

after LF normalization. Every executable check uses the
optimization-safe `require` function. No Python `assert` is used.
