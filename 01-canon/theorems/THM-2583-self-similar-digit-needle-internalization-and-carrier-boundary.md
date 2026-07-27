---
id: THM-2583
title: "Self-similar digit-needle internalization and carrier boundary"
status: >
  PROVED + VERIFIED-EXACT.  Let U be any positive rational Boolean step
  carrier on the physical circle which is neutral for a new F_13 endpoint
  action.  For every positive k, L in {1,2}, and every sufficiently large
  N, the future gate d_L(k T^N x+s/13), T(x)=13x mod 1, has a labelled
  boundary x_N in the interior of U.  One can choose it so that the old
  digit floor(13x_N) equals the immediate future digit
  floor(13T^N x_N).  A sufficiently deep physical base-thirteen cylinder
  H_N containing x_N lies inside U, fixes that root equality on the whole
  cylinder, and meets the complete 26*13^N*k target/tooth boundary atlas
  only at x_N.  Thus H_N internalizes THM-2578's absolute boundary needle:
  the lawful positive danger/safe layers have one target-delta handoff and
  every target Abel character is nonzero at every physical offset.  On each
  live row the same geometry can place such a cylinder inside THM-2559's
  source carrier and a nonzero translated c_3 probe, giving exact base
  incidence with the owner source, terminal source word, late owner block,
  and equal old/future physical root.  But that source carrier contains an
  unshifted k_a safety factor and is not neutral.  Its facts do not transport
  with the shifted cylinder orbit; combining them with the nonzero colours
  would repeat MISTAKE-266.  Thus the live cylinder is an absolute-reference
  candidate with base provenance only.  No covariant packet comparison,
  paired-blocker action, THM-2334 current, row exclusion, or LRC(14)
  conclusion follows.
source: common-endpoint-seam-2026-07-28-self-similar-digit-needle
depends_on:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2559-target-informed-chord-and-universal-old-repair-packet
  - THM-2578-target-neutral-boundary-needle-and-all-colour-abel-normal-carrier
related:
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
  - THM-2565-target-active-self-return-and-future-root-overlap-audit
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
  - THM-2581-b-word-depth-five-owner-clock-host-and-reflection-breaking
script: 04-computation/lrc14_self_similar_digit_needle_thm2583.py
output: 05-knowledge/results/lrc14_self_similar_digit_needle_thm2583.out
script_sha256: a68bf5e2d77c9b376c124fd83b1675d5deaa3c0143bd41c4f01ba931f7988de9
output_sha256: eb882d1552679ca811a3d19b4dcfd258305fdc35d340e4315cfc14b5e97cc39d
hash_basis: working-tree bytes (LF)
---

# THM-2583 -- a later tooth grid internalizes the boundary needle

**PROVED + VERIFIED-EXACT.**

THM-2578 constructs a fixed rational interval which isolates one boundary
of a target-translated danger/safe gate.  Its remaining defect is that the
interval was chosen in the ambient endpoint algebra, with no reason for it
to lie inside a live carrier.

The `13`-adic toothpick self-similarity removes that defect for every
positive rational **target-neutral** carrier.  A later gate has a boundary
grid whose mesh tends to zero.  Once one boundary pierces the interior of
the carrier, a still deeper physical digit cylinder isolates it without
leaving the carrier:

```text
positive neutral rational carrier U
  + future tooth grid at speed 13^N k
  -> equal-root boundary point in int(U)
  -> base-13 cylinder H_N inside U with one boundary trace
  -> absolute all-colour Abel-normal reference.              (1)
```

This is a genuine inheritance theorem only for facts separately proved
neutral under the whole action.  For a nonneutral live packet it is instead
a geometric base-incidence theorem.  It does not make a target-informed or
unshifted moving factor neutral and does not transport the paired endpoint
action.

## 1. The fixed-sign boundary atlas is one uniform grid

Put

```text
T(x)=13x mod 1,                  I_h=(h/13,(h+1)/13),

d_L(y)=1_(||y||<L/14),           L in {1,2}.                 (2)
```

Fix a positive integer `k`.  Before the future pullback, the boundaries of
the thirteen target translates are

```text
y_(j,epsilon,s)
 =(j+epsilon L/14-s/13)/k mod 1,

j in {0,...,k-1},       epsilon in {+1,-1},
s in F_13.                                                (3)
```

For a fixed sign `epsilon`, the map

```text
(j,s) -> 13j-s mod 13k                                  (4)
```

is a bijection.  Indeed equality in (4) gives
`13(j-j')=s-s' mod 13k`; the representative difference on the right has
absolute value at most `12`, so it is zero, and then `j=j'`.  Therefore the
`13k` points with fixed sign in (3) are a translate of the uniform grid

```text
(1/(13k)) Z / Z.                                          (5)
```

No point in (3) lies on a physical base-thirteen digit wall.  If
`y_(j,epsilon,s)=h/13 mod 1`, then after clearing denominators one obtains

```text
14(13j-s)+13 epsilon L-14kh in 182k Z.                    (6)
```

Modulo seven, (6) says `13 epsilon L=0`, impossible for `L=1,2`.
Consequently each open digit cell `I_h` contains exactly `k` points of each
fixed-sign grid, hence `2k` points of the complete base boundary atlas.

The same congruence proves the stronger fact needed below: no boundary in
any pulled-back atlas can equal `a/13^D` for any integers `a,D`.  Multiplying
a base-thirteen rational by an integer leaves only a power of `13` in its
denominator, whereas the right side

```text
epsilon L/14-s/13=(13 epsilon L-14s)/182                 (7)
```

retains a factor `7` after reduction.

## 2. Every positive rational carrier is pierced on an equal-root branch

Let

```text
U:T->{0,1}                                                 (8)
```

be a rational step function of positive integral, fixed independently of
the new target label `s`.  Some digit `h` and some open rational interval

```text
J subset {U=1} intersection I_h,

length(J)=lambda>0                                        (9)
```

exist.  By Section 1, choose once and for all one labelled base boundary

```text
y=y_(j_0,epsilon_0,s_0) in I_h.                           (10)
```

For every `N` with

```text
13^N lambda>1,                                            (11)
```

the preimage grid

```text
x_r=(r+y)/13^N,             r=0,...,13^N-1               (12)
```

has mesh `13^(-N)`, so one of its points lies strictly inside `J`.  Call it
`x_N`.  Since `J subset I_h`, and since `T^N x_N=y`, it satisfies the exact
old/future root diagonal

```text
floor(13x_N)=h=floor(13T^N x_N).                          (13)
```

Moreover (10) makes `x_N` a boundary of

```text
P_s^(N)(x)=d_L(k T^N x+s/13)
          =d_L(13^N kx+s/13),                             (14)
```

with target label `s_0`.  Thus every sufficiently large delay works; there
is no finite-clock extrapolation, mixing limit, or congruence subsequence.

## 3. A physical digit cylinder isolates that boundary inside the carrier

At speed

```text
K_N=13^N k,                                                (15)
```

THM-2578 says the complete `26K_N` boundary atlas of (14) is pairwise
separated.  The chosen `x_N` is in the interiors of all three sets

```text
{U=1},                  I_h,
                  (T^N)^(-1)(I_h),                         (16)
```

and, by (7), it is not on a base-thirteen cylinder wall.  Remove the other
`26K_N-1` atlas points from a sufficiently small neighbourhood in (16).

The unique depth-`D` physical digit cylinders containing `x_N` shrink to
`x_N` as `D` tends to infinity.  Hence for some `D>N` the cylinder

```text
C_N=[a/13^D,(a+1)/13^D),

H_N=1_(C_N)                                                (17)
```

obeys, modulo its null endpoints,

```text
C_N subset {U=1},

floor(13x)=floor(13T^N x)=h             for every x in C_N,

C_N intersection Boundary(P_s^(N):s in F_13)={x_N}.       (18)
```

Because every atlas point retains the septimal denominator obstruction,
neither endpoint of `C_N` is a gate boundary.  Thus `H_N` has two-sided
trace one at `x_N`, zero on a neighbourhood of every other atlas point,
and is a fixed target-neutral physical inverse-branch atom.  Literal
inclusion in (18), rather than an averaged correlation, is what inherits
every neutral fact encoded by `U`.

## 4. The internal all-colour boundary normal

Form the lawful whole layers

```text
A_s=H_N P_s^(N),                 C_s=H_N(1-P_s^(N)).        (19)
```

They are nonnegative, disjoint, and both have positive measure at
`s=s_0`.  Their THM-2573 positive handoff measures are

```text
nu_s=delta_(x_N),                s=s_0,

nu_s=0,                          s!=s_0.                    (20)
```

Consequently, for every physical offset `M in Z` and every target character
`q in F_13`, the plus-sign normalized target DFT of the one-sided
logarithmic Abel normal is

```text
Nhat(q;M)
 =exp(2 pi i Mx_N)zeta_13^(q s_0)/(26 pi^2) !=0.           (21)
```

The target profile in (20) has nonzero augmentation.  It is an absolute
charged reference of the kind required by THM-2579, not a target difference.
Unlike THM-2578's ambient interval, however, this reference is now an actual
subatom of the prescribed neutral carrier `U` and fixes the old/future
physical root on its whole support.

## 5. Uniform live-row base incidence, not transported semantics

On every one of the `165` live rows, retain THM-2559 notation

```text
W(x)=F(x)g(x),

F=1_(E_j) 1_(Q_(j,sigma))(13^k x).                         (22)
```

The positive source strata in THM-2559 imply `integral W>0`.  At the
untwisted base, `W` retains literally:

```text
owner source E_j,                  terminal source word sigma,

late owner block g,                j dangerous,

a and c_3 safe,                    all six guard/unit roles safe.       (23)
```

One may also retain a translated deepest probe before applying the theorem.
For every circle point, the thirteen arcs

```text
d_1(c_3x-r/13),                 r in F_13                 (24)
```

have covering multiplicity one or two: their centres have spacing `1/13`,
while each arc has length `1/7`, strictly between `1/13` and `2/13`.
Therefore some `r_0` makes

```text
U=W d_1(c_3x-r_0/13)                                  (25)
```

positive.  Since `W d_1(c_3x)=0`, necessarily `r_0!=0`.  The geometric
argument of Sections 1--3 applies to this positive rational `U` without any
neutrality assumption: it produces a plain physical digit cylinder

```text
H_N subset W d_1(c_3x-r_0/13).                            (26)
```

Thus at the untwisted base the cylinder has exact incidence with all facts
in (23), old `c_3` safety, and a literal nonzero translated-deep danger
probe.

Apply the geometry with

```text
k=k_a,                       L=L_a.                         (27)
```

Choose `N` beyond the chronology threshold and put `k_a` first in the
freely ordered future bank.  The shifted cylinder layers by themselves are
the lawful auxiliary THM-2578 orbit; at `s=s_0` its danger half contains a
literal later `k_a` gate.  The geometric base ledger is

```text
source owner + source word + late owner block
  + nonzero translated deepest probe
  + equal old/future physical root
  + one later target-labelled boundary.                    (28)
```

All choices may depend on the row.  Because there are only `165` rows and
the piercing assertion holds for every sufficiently large `N`, their
thresholds admit one common larger delay if desired.

The word **base** is load-bearing.  `W` contains the unshifted present
`k_a` safety factor.  Therefore it is not neutral under the action which
moves that gate.  Although `H_N` is a legitimate fixed geometric filter and
(26) proves its untwisted provenance, replacing the lawful shifted layers by

```text
W H_N P_s^(N),                 W H_N(1-P_s^(N))             (29)
```

would freeze the old moving factor.  The nonzero target colours of the
auxiliary `H_N` orbit cannot be combined with (23), (26), or (28) as a live
semantic current.  Doing so would repeat exactly MISTAKE-266.  To transport
the provenance, one must co-shift the whole `W` packet or prove a comparison
map from the fixed digit cylinder to a covariant packet orbit.

By contrast, when `U` is separately proved fixed under the whole target
action—for example an appropriately typed positive member of THM-2549's
future-pullback algebra—Sections 3--4 genuinely retain its neutral facts.
That conditional application does not attach the selected old head or the
missing old target charge.

## 6. Relation to the depth-five owner-clock needle

THM-2581 supplies a complementary endogenous reference on the canonical
`sigma={b}`, `r=5` collision fibre.  There

```text
theta=t-2u=floor(2y) in {0,1},                             (30)
```

the signed half-circle label changes sign under reflection, and the central
owner-clock drift has every nonzero collision-root colour.  That is a
same-fibre oriented binary needle.  But target translation sends
`theta -> theta+q`; the binary chart is not a fixed target-neutral physical
cylinder.

The present `H_N` is the opposite side of the bridge: a full physical
base-thirteen atom, fixed in the auxiliary future endpoint orbit and chosen
with exact base incidence inside the source word/owner carrier.  That base
incidence is not a transported word/owner label, and `H_N` is not yet on
THM-2581's collision `y`-fibre.  The precise missing map is therefore a
temporal intertwiner which carries the collision `y,theta` chart to the later
physical digit-cylinder/endpoint chart while transporting, rather than
freezing, the word, owner clock, and deep sidecars.

## 7. Exact stopping boundary

The theorem internalizes the ambient absolute reference on every genuinely
neutral carrier.  On THM-2559 it proves only that the same geometric cylinder
can be chosen with the desired source word/owner/deep **base provenance** and
a future root equal to the old digit.  It does not transport that provenance
to the nonzero target normal or identify the cylinder with THM-2559's
target-informed **head**.  The source gate, head selector, and slope stratum
are target-active data; freezing any of them would repeat MISTAKE-266.

Nor does (28) co-shift the blocker paired with `k_a`, retain both THM-2563
moving dipoles, construct THM-2537 equation (56), or identify (21) with a
THM-2334 relation-residue current.  Those are packet/action comparison
problems, not further boundary-support problems.  No scalar row is excluded;
the live ledger remains `165`, and LRC(14) remains open.

## 8. Exact companion

Run

```bash
python3 04-computation/lrc14_self_similar_digit_needle_thm2583.py
python3 -O 04-computation/lrc14_self_similar_digit_needle_thm2583.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_self_similar_digit_needle_thm2583.out.
```

The dependency-free exact referee checks:

- `1,045,200` base boundary points in `800` fixed-sign grids, their exact
  uniform spacing, their `10,400` root-cell counts, and the persistent
  septimal denominator obstruction;
- `40` equal-root carrier piercings, `1,845,480` complete-atlas cylinder
  traces, and exact old/future digit equality, with isolating digit depth at
  most five in the controls; and
- `997` rational translated-deep-probe controls exercising both exact cover
  multiplicities.

The arbitrary-carrier, all-sufficiently-large-`N`, cylinder-shrinking, and
live-row base-incidence statements are symbolic proofs above, not finite
extrapolations. **QED.**
