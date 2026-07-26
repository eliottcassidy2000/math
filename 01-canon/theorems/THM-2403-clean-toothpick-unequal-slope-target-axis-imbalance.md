---
id: THM-2403
title: "Clean-toothpick unequal-slope target-axis imbalance"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. On a
  complete clean thirteen-root toothpick, delete q_* and lawfully move
  any retained ordinary q_i along one owner-pivot target generator.
  The paired ordinary target-blocker factor is root-constant and can
  suppress at most two of the thirteen shifts. The resulting endpoint
  mass vector has total minus thirteen times its base mass at least the
  retained parent weight; the sharp abstract gap is one. Rational
  cyclotomic rigidity makes all twelve nonzero target-axis characters
  survive, with exact variance and maximum-mode floors. A common
  Poisson--Abel expansion types these as nonzero bare-deletion target
  marginals. This realizes a lawful unequal-slope deformation of the
  literal deletion support, but it does not restore the deleted factor,
  produce an all-91-unit atom, replace the bare endpoint by the
  canonical fully masked endpoint, exclude a row, or prove LRC(14).
source: codex-2026-07-26-clean-target-axis-imbalance
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2384-owner-pivot-primal-dual-obstruction-and-two-probe-repair-corner
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
  - THM-2397-clean-root-same-parent-charged-role-partition
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2366-retained-probe-target-covariance-and-subthirteen-budget-bridge
  - THM-2400-clean-parent-root-gauge-quotient-and-target-slope-boundary
  - THM-2401-common-filter-endpoint-or-first-death-certificate
script: 04-computation/lrc14_clean_target_axis_imbalance_thm2403.py
output: 05-knowledge/results/lrc14_clean_target_axis_imbalance_thm2403.out
script_sha256: 9c26e6372d2468c26b65daf9b898a93480b0ebb694ae8d45e7787d67a033d91f
output_sha256: e64d9ad3ec7eb72ba0c943cbcbe4231d7485ae7fb939bdc50708aaef0ea74c17
hash_basis: working-tree bytes (LF)
---

# THM-2403 -- unequal slopes force target imbalance on a clean toothpick

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2397 identifies the exclusive `q_*` word with literal
single-factor deletion support at one common bare endpoint. THM-2400
records the remaining difficulty: a true target covector does not
translate all six unit factors with one common slope.

For the standard owner-pivot target this unequal-slope branch is not a
hostile. Only one retained unit factor moves. Its paired target blocker
is constant across the inverse-root fibre and can erase at most two
shift values. The exact clean-toothpick incidence then forces a strict
mass imbalance:

```text
clean one-double root cover
  + one retained ordinary-role translation
  + one ordinary target-blocker safe gate
  -> a nonconstant rational C_13 target-twist bank
  -> all twelve nonzero target-axis characters.                 (1)
```

This is a lawful bare-deletion target marginal. The theorem deliberately
stops before canonical fully masked endpoint restoration.

## 1. The abstract clean-cover lemma

Put

```text
G=F_13.
```

Let

```text
Q_*,Q_i,Q_2,Q_3,Q_4 subset G,       |Q_j|=2,

H subset G,                          |H|=4,          (2)
```

and suppose these six labelled masks cover `G`. Their total incidence is

```text
4+5*2=14.                                             (3)
```

Hence exactly one root has multiplicity two and the other twelve have
multiplicity one. There are no uncovered roots and no triple root.

Delete the two ordinary labels `Q_*` and `Q_i` and retain the root core

```text
B=G minus (H union Q_2 union Q_3 union Q_4).          (4)
```

The roots in `B` are precisely those whose remaining incidence lies in
`Q_* union Q_i`. Put

```text
b=|B|,

A_0=B minus Q_i,

a=|A_0|.                                             (5)
```

The set `A_0` is the exclusive `Q_*` support after its safe factor has
been deleted. The unique double root gives exactly three possibilities:

```text
(a,b) in {(1,3),(2,3),(2,4)}.                        (6)
```

Indeed:

- if the double uses `Q_*`, then one `Q_*` root is not exclusive and
  `(a,b)=(1,3)`;
- if the double uses `Q_i` but not `Q_*`, then `(a,b)=(2,3)`; and
- if the double uses neither deleted label, then `(a,b)=(2,4)`.

Let `alpha in F_13^*` and translate only the retained ordinary word:

```text
A_s=B minus (Q_i+alpha s),             s in F_13.    (7)
```

For each fixed `r in B`, exactly two values of `s` put `r` in
`Q_i+alpha s`. Therefore

```text
sum_s |A_s|=11b.                                    (8)
```

This is a complete finite incidence identity, not an averaged
equidistribution assertion.

## 2. An ordinary target-blocker gate cannot flatten the bank

Let

```text
g_s in {0,1},

g_0=1,

#{s:g_s=0}<=2.                                      (9)
```

Since `|A_s|<=b`, equations (8)--(9) give

```text
sum_s g_s|A_s|
 >=11b-2b
 =9b.                                               (10)
```

At the base shift,

```text
g_0|A_0|=a.
```

Thus

```text
sum_s g_s|A_s|-13g_0|A_0|
 >=9b-13a
 >=1.                                               (11)
```

The three typewise floors are

```text
(a,b)=(1,3):       14,

(a,b)=(2,3):        1,

(a,b)=(2,4):       10.                              (12)
```

Equation (11) is sharp. Take

```text
B={0,1,3},             Q_i={2,3},

g_2=g_3=0,             g_s=1 otherwise.
```

Then

```text
(g_s|A_s|)_(s=0)^12
 =(2,2,0,0,3,3,3,3,3,3,2,1,2),                    (13)

sum_s g_s|A_s|-13g_0|A_0|=27-26=1.
```

The two killed shifts are realized by one strict ordinary danger comb:
at phase `z=5/26`,

```text
||z-s/13||<1/14

iff s in {2,3}.                                     (14)
```

The sharpness claim is for the labelled root/gate universe. It does not
assert that the complete word (13) occurs on a canonical scalar row.

## 3. Weighted parent sets and all twelve target characters

Let `Y` be a finite rational interval set. Permit the clean masks in
Section 1 and their common root origin to vary measurably with `y in Y`,
while retaining the same labelled incidence hypotheses. Let

```text
f:Y->[0,infinity)
```

be a nonzero rational step weight independent of the root `r`, and put

```text
rho_f=integral_Y f(y)dy>0.                           (15)
```

For every `y`, let `g_s(y)` satisfy (9), and define the unnormalized
fibre-sum masses

```text
N_s=integral_Y f(y)g_s(y)|A_s(y)|dy.                (16)
```

Integrating (11) gives the strict anchored imbalance

```text
Delta_N:=sum_s N_s-13N_0>=rho_f.                    (17)
```

In particular `N=(N_s)` is not constant. Every `N_s` is rational. With
`zeta=exp(2*pi*i/13)`, put

```text
Ntilde(k)=(1/13)sum_s N_s zeta^(k s).               (18)
```

If one `Ntilde(k)` vanished for `k!=0`, the rational polynomial

```text
P(X)=sum_s N_s X^s
```

would be divisible by `Phi_13(X)=1+X+...+X^12`.
Since `deg P<=12`, all thirteen coefficients would be equal, contrary
to (17). Therefore

```text
Ntilde(k)!=0                 for every k!=0.         (19)
```

The finite variance identity and Cauchy--Schwarz give

```text
sum_(k!=0)|Ntilde(k)|^2
 >=Delta_N^2/2028
 >=rho_f^2/2028,                                    (20)

max_(k!=0)|Ntilde(k)|
 >=Delta_N/156
 >=rho_f/156.                                       (21)
```

Equality in the first inequality of (20), for fixed `N_0` and total
mass, occurs when the other twelve entries are equal. Rationality is
needed for the every-colour conclusion (19), but not for the energy
floor (20).

Under physical circle normalization,

```text
M_s=N_s/13,                     Mtilde=Ntilde/13.    (22)
```

Hence the corresponding floors are

```text
sum_(k!=0)|Mtilde(k)|^2>=rho_f^2/(169*2028),

max_(k!=0)|Mtilde(k)|>=rho_f/2028.                  (23)
```

The factor thirteen in (22) is the Jacobian/discrete-root normalization;
it must not be silently dropped.

## 4. The ordinary blocker gate is exactly the LRC target gate

Retain THM-2392/2397's complete clean-parent set `S`. For `y in S`, use
the inverse roots

```text
x_r=(y+r)/13.
```

The guard and five ordinary unit dangers give exactly the six masks in
(2)--(3). Choose the graft helper for target `a` to be any lower
ordinary label `q_i!=q_*`.

Let `c_a=13C_a` be any ordinary target blocker. Since `S` lies outside
the quotient-blocker cage,

```text
1_(D_(C_a)^c)(y)=1.                                 (24)
```

The THM-2384 owner-pivot dipole, built on THM-2309's exact pivot, is

```text
eta_a=e_a-e_(q_i).                                  (25)
```

At target shift `s`, its two moving safe factors are

```text
1-d(q_i x+s/13),

1-d(c_a x-s/13),                                    (26)
```

where `d(t)=1_(||t||<1/14)`. On the inverse roots,

```text
d(q_i x_r+s/13)
 =d((q_i y+q_i r+s)/13),                            (27)
```

so the `q_i` mask is translated through every root shift because
`q_i mod 13` is a unit. Meanwhile

```text
g_s(y)
 :=1-d(c_a x_r-s/13)
 =1-d(C_a y-s/13)                                   (28)
```

is independent of `r`. Equation (24) gives `g_0(y)=1`. The exact
thirteen-root count, almost everywhere away from strict endpoints, is

```text
sum_s d(C_a y-s/13)
 =2-d(c_a y)
 <=2.                                               (29)
```

Thus (28) satisfies precisely the gate hypotheses (9). Equations
(17)--(23) apply pointwise to the lawful target action (25), not to an
auxiliary duplicate probe.

The other two blocker-safe factors are also root-constant at the base
clean parent. They may be retained in `f`; the moving target blocker is
kept separately as `g_s`.

## 5. Target-neutral delayed words preserve the imbalance

Let `Q` be any nonnegative rational circle word and let `R=13^m`,
`m>=1`. On the inverse-root fibre,

```text
Q(Rx_r)=Q(13^(m-1)y),                               (30)
```

so the delayed word is root-constant and target-neutral modulo
thirteen. Put

```text
f_R(y)=1_S(y)Q(13^(m-1)y).
```

Whenever

```text
rho_R=integral_S f_R(y)dy>0,                        (31)
```

Sections 3--4 give all twelve target-axis characters with `rho_f=rho_R`.
THM-2397 proves (31) for every sufficiently large clock when `Q` has
positive mass, via its two-BV mixing estimate. Thus a fixed positive
delayed terminal word cannot flatten the unequal-slope target bank.

For the unfiltered choice `f=1_S`, THM-2396 gives

```text
rho_f=delta>=1/26754                                (32)
```

uniformly, and `delta>=66/4459` on the common-core chain. The
fibre-sum floors (20)--(21) become

```text
universal:
  energy>=1/1451594774448,
  max mode>=1/4173624;

common core:
  energy>=363/3360173089,
  max mode>=11/115934.                              (33)
```

The physical endpoint-normalized floors from (23) are

```text
universal:
  energy>=1/245319516881712,
  max mode>=1/54257112;

common core:
  energy>=363/567869252041,
  max mode>=11/1507142.                             (34)
```

No denominator-free positive floor is claimed for each individual one
of the twelve nonzero modes.

## 6. Poisson--Abel target-current typing

Let `E_s` be the actual bare-deletion endpoint packet: retain the guard
and the four ordinary safe factors other than `q_*`, move `q_i` and
`c_a` as in (26), retain the other root-constant blocker factors, and
multiply by the delayed word in Section 5. Then

```text
mu(E_s)=M_s.                                        (35)
```

Poisson-smooth all displayed interval factors at one common radius
`0<rho<1`. At finite `rho` their Fourier series are absolutely
convergent, so the product may be expanded and the finite `s`-sum may
be taken termwise. If `n_a,n_i` are the harmonics of the two moving
factors, the target phase is

```text
zeta^(s(n_i-n_a)).
```

Consequently

```text
(1/13)sum_s mu(E_s)zeta^(k s)                       (36)
```

selects exactly the relation terms satisfying

```text
n_a-n_i=k mod 13,

eta_a.r=k mod 13.                                  (37)
```

The remaining target coordinate is not projected in (36). Thus a
nonzero coefficient in (36) is a nonzero target **marginal**:

```text
sum_(q_b in F_13) A(k,q_b)!=0,
```

and therefore some full target vector with first coordinate `k` has a
nonzero coefficient. It is not a claim that the pure-axis fibre
`(k,0)` survives.

Bounded-product `L^1` convergence identifies the `rho->1-` boundary
with (35)--(36). Internal infinite series are never reordered after
the Abel limit. Equations (19) and (36)--(37) therefore give:

> **Lawful bare-deletion target conclusion.** For every
> `k in F_13^*`, the owner-pivot target-axis marginal of the
> `q_*`-deletion packet is nonzero. Hence at least one full nonzero
> target fibre with first coordinate `k` survives.

This is a genuine target-shift-covariant realization of the literal
bare-endpoint deletion support. It is not merely a root-gauge
relabeling: the unequal slopes change `|A_s|`, and that change is the
source of (17).

## 7. Sharp boundaries and hostiles

Four boundaries remain load-bearing.

1. **The fully masked endpoint is zero on the clean cover.** Restoring
   the deleted `q_*`-safe factor annihilates `A_0`. The theorem does not
   identify the bare deletion current with the canonical fully masked
   owner endpoint.
2. **The all-`91`-unit mask is not retained.** Equation (37) fixes one
   mod-thirteen target coordinate but does not control every atomic
   relation coordinate modulo seven or thirteen. THM-2334's
   unrestricted target aggregate and its all-unit aggregate remain
   distinct.
3. **A guard target gate would not satisfy (9).** A width-two guard
   danger can occupy four root shifts. The uniform estimate would then
   be `7b-13a`, which is negative for `(a,b)=(2,3)`. The actual three
   target blockers are ordinary, so this hostile does not occur in
   (25)--(29).
4. **Arbitrary signed or root-dependent filtering can cancel.** The
   proof uses one nonnegative root-constant weight `f` and the exact
   ordinary blocker gate. It does not cover a change of endpoint
   section, signed mixing, or an extra root-sheet-dependent terminal
   filter.

The gap-one example (13)--(14) shows that no stronger universal integer
imbalance follows from the clean incidence and two-shift gate alone.
The prime-rational all-colour step is also sharp in type: over composite
cyclic groups, nonconstant rational mass vectors can have vanishing
nontrivial characters.

## 8. LRC interface

The new proof-graph edge is

```text
uniform clean q_* deletion support
  + one standard owner-pivot target generator
  + any positive root-constant delayed word
  -> all twelve nonzero target-axis marginals
     at the lawful common bare endpoint.                         (38)
```

This removes two former debts:

- unequal factor slopes are no longer merely a deformation hostile on a
  clean toothpick; and
- the literal deletion support now has a lawful target-shift-covariant
  Poisson--Abel realization.

It does **not** remove the final endpoint/owner debt:

```text
bare deletion target marginal
  -> canonical fully masked owner endpoint
     with the required all-91-unit/current-phase typing.         (39)
```

THM-2401 remains complementary. If a common Boolean terminal filtration
of the two charged descendants is constructed, its positive branch
preserves all colours and its first-death branch identifies the failed
factor. The present theorem instead supplies an unconditional lawful
bare-deletion target marginal from the exact clean incidence.

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 9. Exact companion

The dependency-free exact companion:

- verifies canonical clean-cover examples for all three types in (6);
- exhausts all `13,728` residual cores/retained two-root words and all
  `1,084,512` gates with at most two killed nonzero shifts;
- obtains the exact nine-entry sharp gap table and `73,220` distinct
  nonconstant mass vectors;
- checks the gap-one word and exact strict circle gate (13)--(14);
- checks the cyclotomic reduction and variance bound on every distinct
  mass vector;
- verifies all `2,197` target-phase projector cases; and
- checks every fraction in (33)--(34).

Run

```bash
python3 04-computation/lrc14_clean_target_axis_imbalance_thm2403.py
python3 -O 04-computation/lrc14_clean_target_axis_imbalance_thm2403.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_clean_target_axis_imbalance_thm2403.out
```

Every truth-bearing finite check raises explicitly, so optimized mode
executes the same audit.
