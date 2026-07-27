---
id: THM-2548
title: "Seven-step C91 transfer, full-norm separation, and the two-cut arrival criterion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The theorem identifies the exact
  integral kernel and image of the seven-step transfer on THM-2542's C91
  mapping torus and gives a conditional directed-gain/Hall arrival criterion.
  It does not construct the required common-ancestry transition automaton,
  a semantic vertical edge, a Hall-deficient live table, or an LRC(14)
  contradiction.
source: codex-2026-07-27-root-holotopy-transfer
depends_on:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
  - THM-2539-diagonal-cubic-owner-clock-boundary-current
  - THM-2543-augmentation-norm-relative-phase-local-system-dichotomy
related:
  - THM-840-hamming-five-continuation-congruence-boundary
  - THM-1254-coherent-blocker-cycle-chronological-descent
  - THM-2089-persistent-cut-affine-holonomy-normal-form
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
script: 04-computation/lrc14_c91_seven_step_transfer_thm2548.py
output: 05-knowledge/results/lrc14_c91_seven_step_transfer_thm2548.out
script_sha256: 5578734098a94c0ef5a63f867dab95dff17ad3e642c37be2a44a54eecb341293
output_sha256: d15e165f70eac230fe71ef57a72114834f13a026676adb702db8d68675a4db49
hash_basis: LF-normalized bytes
---

# THM-2548 -- the partial transfer remembers the root deck

Let

```text
G=F_7 x F_13,
S_a(k,r)=(k+1,r+a),                  a in F_13^*.            (1)
```

THM-2542 proves that (1) is the abstract mapping torus of its seven root
charts.  It is one `C91` orbit.  This theorem separates two superficially
similar operations on that orbit:

```text
D_a=1+S_a+...+S_a^6,                 seven-step transfer;
N_a=1+S_a+...+S_a^90,                full-orbit norm.         (2)
```

The distinction is exact over both `Z` and every characteristic-zero Fourier
field.  The partial transfer kills only root-uniform clock augmentation and
retains every root-charged character.  Its integer image is primitive: there
is no hidden `13`-adic index.  The full norm, by contrast, kills every
nonconstant character.

This is a cancellation theorem inside horizontal degree.  It does not create
THM-2542's missing semantic vertical arrow.  Combining its gain-cover view
with THM-2545 instead produces an exact two-cut closure architecture: a
directed gain cut must first force target activity, and an independent Hall
cut must then force same-root arrival.

## 1. The integral seven-step exact sequence

Let `Z^G` be the integer tables on `G`; use the pullback action of (1), so a
choice of inverse convention only reverses every displayed character.  For a
table `f`, put

```text
Y_k(f)=sum_(r in F_13) f(k,r).                              (3)
```

Define the two primitive lattices

```text
A_clk={f:f(k,r)=c_k independent of r, sum_k c_k=0},
L_eq ={g:Y_0(g)=Y_1(g)=...=Y_6(g)}.                         (4)
```

Then

```text
0 -> A_clk -> Z^G --D_a--> L_eq -> 0                       (5)
```

is exact.  In particular,

```text
ker D_a=A_clk,                  rank ker D_a=6;
im  D_a=L_eq,                   rank im  D_a=85;             (6)
```

and the Smith form of `D_a` is

```text
diag(1,...,1,0,...,0)=1^85 direct-sum 0^6.                 (7)
```

Thus the image is saturated and `Z^G/im D_a` is free of rank six.
As a module over the orbit group ring, one may identify

```text
coker D_a = Z[z]/(Phi_7(z)) = Z[zeta_7].                    (7a)
```

### Proof of the kernel

The polynomial identity

```text
(S_a-1)D_a=S_a^7-1                                        (8)
```

shows that `D_a f=0` implies `S_a^7 f=f`.  Seven applications of (1) fix the
clock and translate the root by `7a`, a generator of `F_13`.  Hence `f` is
constant in `r` on every clock fibre.  On such a table, `D_a f` is the
constant table `sum_k c_k`.  This proves the first line of (6), including the
converse.

### Proof of the image, including integral surjectivity

For arbitrary `f`, summing `D_a f` over the thirteen roots at clock `k`
visits every clock exactly once and merely permutes each root fibre.  Hence

```text
Y_k(D_a f)=sum_((j,r) in G)f(j,r),                          (9)
```

independent of `k`.  Thus `im D_a` is contained in `L_eq`.

For the converse, enumerate the single `S_a` orbit by `n in Z/91`, so fixed
`n mod 7` is one clock and the thirteen values `n=k+7t` run through all its
roots.  With the convention

```text
(D x)_n=sum_(j=0)^6 x_(n-j),                               (10)
```

the difference equation is

```text
x_n-x_(n-7)=y_n-y_(n-1).                                  (11)
```

For each of the seven residue classes modulo seven, choose one initial
integer and solve (11) around its thirteen-cycle.  The closing condition is
exactly

```text
Y_k(y)-Y_(k-1)(y)=0.                                      (12)
```

When (12) holds, the resulting error `y-Dx` is constant on the whole
91-cycle.  Adding that integer constant to one complete root fibre of `x`
changes every seven-term window by the same amount, removing the error.
Therefore every integer `y in L_eq` has an integer preimage.

Finally, the quotient map

```text
y -> (Y_0-Y_6,...,Y_5-Y_6) in Z^6                           (13)
```

is onto: the six point masses at `(k,0)`, `0<=k<6`, map to the standard
basis.  Its kernel `L_eq` is primitive.  Equations (5)--(7) follow.

### Coprime-prime cycle form

Nothing in the proof uses the special values except coprimality.  On
`C_q x C_p`, with `p,q` distinct primes and a nonzero root translation in
`C_p`, the
`q`-step transfer has

```text
kernel = root-uniform augmentation on the q clocks,
image  = tables whose p-root-fibre sums are clock-independent,
Smith  = 1^(pq-q+1) direct-sum 0^(q-1).                    (14)
```

This is the integral transfer exact sequence of the skew mapping torus.

## 2. Fourier support: all root charge survives

Let `xi` and `zeta` be primitive seventh and thirteenth roots.  On

```text
chi_(beta,alpha)(k,r)=xi^(beta k)zeta^(alpha r),             (15)
```

the seven-step multiplier, up to the harmless global sign convention, is

```text
d_a(alpha,beta)
 =sum_(j=0)^6 (xi^beta zeta^(alpha a))^j.                   (16)
```

If `alpha=0` and `beta!=0`, (16) is the sum of all seventh roots and vanishes.
If `alpha!=0`, then

```text
(xi^beta zeta^(alpha a))^7=zeta^(7 alpha a)!=1,             (17)
```

so the geometric sum cannot vanish.  At `(alpha,beta)=(0,0)` it equals
seven.  Therefore

```text
d_a(alpha,beta)=0 iff alpha=0 and beta!=0.                  (18)
```

The twelve root-only and all seventy-two primitive mixed characters survive.
Equivalently, a seven-chart signed table can cancel under `D_a` only when it
is root-uniform clock augmentation.  Any nonzero root-charged component
survives without selecting a root origin.

This last statement is stronger than a nonzero determinant over a cyclotomic
field: (7) proves that there is no integral torsion invoice hidden behind it.
The apparent componentwise cyclotomic resultants do not describe the integral
regular lattice, whose components are not an integral direct product.

## 3. The full norm erases what the partial transfer keeps

Put

```text
R_a=sum_(t=0)^12 S_a^(7t).                                  (19)
```

Then

```text
N_a=R_a D_a.                                                (20)
```

The first factor in (20) is the complete norm along the thirteen-point root
fibre.  It kills every `alpha!=0` character.  The second factor kills the six
remaining nontrivial clock characters.  Hence

```text
N_a f=(sum_G f) 1_G,                                       (21)
```

and `N_a` kills all ninety nonconstant characters.

Thus a proof which first sums the seven transported charts may retain every
root-deck charge, while a proof which averages the full C91 orbit necessarily
forgets root, clock, and their relative phase.  The right order is therefore

```text
pair with a labelled contragredient for every charge intended to survive;
apply the seven-step transfer;
only then take an invariant or norm if the target permits it.               (22)
```

Equation (22) is a preservation rule, not a construction of the sidecar.
In particular, a contraction which has already summed away `alpha` and leaves
only `beta!=0` lands exactly in (18)'s kernel; `D_a` cannot repair that loss.
THM-2543's future-phase augmentation/norm dichotomy is a different `C7`
module.  Its norm is a multiplicative product of seven phase means, whereas
`N_a` is the additive full-orbit transfer (21).

Three boundaries are sharp.  First, if `a=0`, then the root and clock do not
form one C91 orbit: `D_0` is the seven-clock norm separately on every root and
kills all `13*6=78` modes with `beta!=0`, including mixed modes.  Second, the
two positive root-uniform sheets `1_(k=0)` and `1_(k=1)` have the same image
`1_G`; partial transfer does not recover an absolute clock origin.  Third,
the Fourier wording here is characteristic zero (or characteristic prime to
`7*13`); in characteristic seven the trivial multiplier also degenerates.

## 4. A 1,296-coefficient signed cut amplification

There is one lawful coefficient-level application to the current horizontal
bank.  Before the character contraction in THM-2539 equation (29), retain its
finite labelled direct-sum table and denote the three slope pieces by
`J_i(u,ell,s)`, where `u` is the physical root, `ell` the seven-valued visible
row, and `s` collects the remaining retained labels.  This is a reindexing of
the table already summed there, not a new common-base pairing.  In the
nonconstant branch of THM-2543 use its analogous relative-phase table.  Center
the seven-valued source/relative-phase coordinate by

```text
d_i(u,ell,s)=7J_i(u,ell,s)-sum_(ell')J_i(u,ell',s).          (22a)
```

THM-2535's finite reindexing is linear, so it extends coefficientwise over
the `Q`-vector space of rational step functions and finite labelled direct
sums containing these `J_i`.  No multiplication of independently realized
physical fixtures is being introduced.

Then `sum_ell d_i=0`, while every nontrivial clock coefficient of `J_i` is
simply multiplied by seven.  Thus none of the existing `216=3*6*12`
nonzero signed source/target/root coefficients is lost.

Apply THM-2535's boundary-relative cut transform at slope
`tau=lambda_i` and any cut scale `c in F_7^*`.  Up to a nonzero character
phase, the new multiplier is

```text
sum_(j=0)^6
 (zeta^(-alpha lambda_i) xi^(kappa c^(-1)))^j.              (22b)
```

On the old physical-root diagonal `alpha=lambda_i b`, its base has seventh
power

```text
zeta^(-7 lambda_i^2 b)!=1                                  (22c)
```

for every `b in F_13^*`.  Hence (22b) is nonzero for all six cut scales.
Thus every horizontal bank to which THM-2539 or the corresponding THM-2543
branch supplies its `216` incidences extends over all six cut scales and has

```text
3*6*12*6=1,296                                             (22d)
```

labelled nonzero signed cut/clock/target coefficients.

This is a direct-sum coefficient amplification, not a physical pairing of
independently constructed fixtures.  Its scope is correspondingly narrow.
Centering in (22a) destroys a fixed positive orientation and literal single
owner.  In THM-2543's augmentation lane the clock is relative future phase;
in its norm lane it is THM-2539's visible source row.  The two are not thereby
identified with one semantic coordinate.  Equations (22a)--(22d) produce no
positive Boolean event, common ancestry, or arrival.

## 5. What the transfer proves, and what it cannot prove

Suppose a signed current is already lawfully defined on one common C91
mapping-torus carrier.  If its root-charged projection is nonzero, the sum of
its seven consecutive covariant transports is nonzero by (18).  This removes
one possible cancellation at the horizontal gluing stage.

But `D_a` is a sum of invertible horizontal chart arrows.  It has vertical
semantic degree zero.  It cannot create a target-active head, a later owner,
or a common physical ancestry.  This has an exact THM-2545 control.  Tensor a
single positive mapping-torus atom with the two couplings

```text
aligned={(0,0),(1,1)},       swapped={(0,1),(1,0)}.          (23)
```

After seven-step transfer, both have head and later-root margins `(7,7)`.
The aligned diagonal has mass fourteen and the swapped diagonal has mass
zero.  Every root-charged mapping-torus mode survives in both.  Thus partial
transfer does not recover the joint root pairing discarded by one-point
marginals, exactly as THM-2545 predicts.

## 6. Directed neutral gain is the first arrival cut

The holotopy invoice nevertheless gives a finite conditional test.  Let

```text
Gamma -> C_7                                                   (24)
```

be a **lawful common-ancestry directed transition automaton**.  Its vertices
retain every coordinate needed for composition: terminal word `sigma`,
source/deep sheets, owner role, clock chart, root chart, and relevant carry.
An edge `e` over one clock step carries

```text
c(e) in F_13,                root correction;
chi(e) in {0,1},             target-active flag.             (25)
```

Pairwise positive edges do not suffice: (24) must be built on one common
rational-BV atom refinement, or an equivalent certified fibre product, so
seven edges compose on one actual ancestry.

For one word stratum, define the neutral directed gain spectrum

```text
Z^0_sigma={sum_(k=0)^6 c(e_k):
  (e_0,...,e_6) is a closed sigma-preserving degree-one section,
  chi(e_k)=0 for every k}.                                  (26)
```

THM-2542's horizontal class is `7a`.  A root-trivializing physical section
must therefore satisfy

```text
sum_k c(e_k)=-7a.                                           (27)
```

Among existing root-trivializing sections,

```text
every such section has a target-active edge
  iff -7a notin Z^0_sigma.                                  (28)
```

This is immediate from (26)--(27), but it is the exact first semantic gate.
It is gauge invariant: replacing `c(e)` by
`c(e)+phi(target(e))-phi(source(e))` telescopes on a closed section.

Two useful sufficient certificates are:

1. `c` is a coboundary on the neutral subgraph, so every neutral closed gain
   is zero;
2. the neutral subgraph has no directed section winding once around the
   clock, for example because a lawful chronological height strictly
   increases on every neutral edge.

The exact decision procedure is reachability in the thirteen-fold gain
cover.  Starting from `(v,0)`, follow only neutral edges for seven layers and
ask whether `(v,-7a)` is reachable.  Equivalently, use Boolean group-ring
transition matrices and inspect the coefficient of `z^(-7a)` on the diagonal
of their product.  A failed reachability search yields a finite cut in the
gain cover; every root-trivializing path must cross an active edge.

There is a quantitative min-plus form.  Give edge `e` cost `chi(e)` and let

```text
m_sigma(r)=minimum active-edge count among closed sections of gain r.       (29)
```

If a positive common-base family of root-trivializing sections has mass
`rho` and `m_sigma(-7a)=m>=1`, then the sum of the seven active incidence
masses is at least `m rho`; some fixed clock carries active mass at least

```text
m rho/7.                                                       (30)
```

This integration step is valid only because the sections live on one common
base.

### Why independent neutral branching is anti-coercive

Suppose instead that neutral choices factor independently by clock layer,
with nonempty gain sets `D_0,...,D_6 subset F_13`.  Then

```text
Z^0=D_0+...+D_6.                                             (31)
```

Cauchy--Davenport gives

```text
|Z^0|>=min(13,1+sum_k(|D_k|-1)).                            (32)
```

If the sum in (32) is at least twelve, `Z^0=F_13`, and the neutral system can
pay every holonomy invoice.  For example seven independent three-choice
layers already fill `F_13`.  More spectral support or more local branches is
therefore not automatically helpful: closure needs word/chronological
correlation which removes the neutral sumset freedom.

## 7. Hall alignment is the second, independent cut

The output of (28) has two possible physical types.

1. If an active vertical edge is already a genuine common-ancestry
   selected-head-to-later-arrival 2-cell, then its positive incidence proves
   the desired hit directly.
2. If it records only a genuinely later target-active role and its root `b`,
   (28) removes THM-2545's cemetery alternative but does not identify `b`
   with the selected head root `h`.

In the second case, push the common-base atoms forward by

```text
omega -> (sigma(omega),h(omega),b(omega))                   (33)
```

to obtain THM-2545's exact table `C^sigma_(t,s)`.  A separate word-stratified
Hall deficiency is then necessary and sufficient to force diagonal mass.
The two stages are genuinely independent: the uniform active rule

```text
b=h+1 in F_13                                                (34)
```

has a later active role on every atom, equal uniform root marginals, and zero
diagonal.

The sharp post-THM-2545 architecture is therefore

```text
gain-cover reachability cut       -> target activity;
transportation/Hall min-cut       -> same-root arrival.      (35)
```

Both cuts must be computed on one lawful atomized carrier.  If the first
vertical edge is already diagonal, the second cut collapses.  Ordinary
linear homology or Smith nonmembership is only a sufficient relaxation of
the first gate: it forgets direction and positivity.  The directed gain
cover is exact.

## 8. Exact referee and stopping boundary

Run

```bash
python3 04-computation/lrc14_c91_seven_step_transfer_thm2548.py
python3 -O 04-computation/lrc14_c91_seven_step_transfer_thm2548.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_c91_seven_step_transfer_thm2548.out.
```

The dependency-free referee constructs the `91 x 91` integer transfer
matrix; checks rank `85`, six explicit kernel generators, an explicit
85-vector basis of `L_eq`, and an integral preimage of every basis vector;
verifies the primitive quotient controls giving (7), all `91` instances of
(8), and the positive collision of two distinct full-root clock sheets;
checks all `1,092`
characters for all twelve `a!=0` over `F_547` and all `1,092` factorizations
(20); records `1,008/1,008` root-charged survivors, `72/72` killed clock-only
modes, `1,080/1,080` killed nonconstant full-norm modes, all `78` clock modes
killed by the zero-holonomy hostile, and all `1,296` multipliers in (22d); and
verifies the transferred aligned/swap Hall hostile, one-vertex neutral
flat/hostile/three-choice gain controls, and the active cyclic-shift Hall
boundary.  The one-vertex controls test the gain arithmetic only; they are
not evidence that the typed live automaton (24) has been constructed.

The theorem does not build (24), prove a neutral-gain exclusion on any live
row, construct the later root selector in (33), prove a Hall deficiency,
identify the scheduler/root mapping torus with the inherited source/deep
sheets, produce a positive semantic arrival, remove one of the `165` rows, or
prove LRC(14).

**QED, conditional where (24), (29), and (33) are invoked.**

The independent hostile audit reconstructed the integral exact sequence and
`Z[zeta_7]` cokernel, the Fourier/full-norm split, all `1,296` cut
multipliers, the gain-cover/Hall separation, and the `m rho/7` floor.  It also
checked the zero-holonomy, positive-sheet, transferred Hall, and one-vertex
gain hostiles.  Normal, optimized, and stored executions byte-match the
declared LF hashes; the documentation check passes.
