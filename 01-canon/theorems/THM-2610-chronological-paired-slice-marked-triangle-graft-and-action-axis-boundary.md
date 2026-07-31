---
id: THM-2610
title: "Chronological paired-slice marked-triangle graft and action-axis boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A fixed THM-2349/2537 marked source coefficient can be multiplied, before
  the finite shift transform, by one sufficiently late same-root THM-2599
  target-danger/blocker-safe chamber.  One common delay retains the exact
  old X,Y,m, word, owner, root, and 91-unit deep triangle while producing
  every one of the twelve nonzero future moving-shift characters.  This is
  a genuine fixed-frequency chronological escape from THM-2568's same-time
  full-X annihilation.  The future character is nevertheless an independent
  local action axis: positive Koopman time kills the old C13 root deck and
  rebases C91 to C7.  No relation-residue identification, future terminal
  word, scalar-row exclusion, or LRC(14) conclusion is proved.
source: kind-pasteur-2026-07-28-chronological-holotopy-graft
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2599-rootwise-opposite-shift-paired-slice-law
  - THM-2604-unshifted-future-root-accessibility-and-selector-cross-mixing-boundary
related:
  - THM-2568-full-x-transition-annihilation-and-refined-pair-drift-boundary
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-2607-constant-six-rail-boundary-holonomy-invoice
  - THM-2608-alternative-rail-clock-collapse-and-missing-transition-index
script: 04-computation/lrc14_chronological_marked_triangle_graft_thm2610.py
output: 05-knowledge/results/lrc14_chronological_marked_triangle_graft_thm2610.out
script_sha256: 07be5ffd964c541c01149d2f2a79832ac5c92a9685521bc51a2bd20ee87db91f
output_sha256: af60cc9bf82ee581277e48c3494c03a22a45cca57dc79e8684031be1d37605ba
hash_basis: LF-normalized bytes
---

# THM-2610 -- a marked source triangle survives every later paired colour

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2599 puts the THM-2537 selected source and a later opposite-shift
target-danger/blocker-safe event on one literal orbit.  Its scalar overlap
retains the source word and owner only as provenance.  It does not retain the
old exact marked coefficient.  THM-2478, in a different future-collision
setting, retains an old finite component in every later colour but explicitly
does not retain the same physical frequency `X` or deep multiplier `m`.

The two mechanisms can be made to meet before either of those losses.  Fix the
old marked Fourier coefficient first.  Then use BV mixing on each of the
thirteen nonnegative future-shift packets and only afterward take their finite
Fourier transform.  The zero temporal sideband is the old coefficient times
the nonzero future mean, and the remaining Koopman sidebands are too small to
cancel it at a sufficiently late common clock.

This produces the first exact same-`(X,Y,m)` chronological paired-slice graft.
It also exposes the remaining holotopy defect cleanly: source and future carry
two locally isomorphic `C_13` character groups, but positive time supplies no
connection identifying them.

## 1. The abstract fixed-frequency chronological graft

Put

```text
p=13,                  T(x)=13x mod 1,

zeta=exp(2 pi i/13).                                      (1)
```

Let `f` be a nonzero rational Boolean circle-step function and fix an integer
`X` such that

```text
f_hat(X)!=0.                                               (2)
```

Let `c_s`, `s in F_13`, be nonnegative rational BV step functions.  Define

```text
g_beta(y)
 =1/13 sum_(s in F_13)c_s(y)zeta^(-beta s),

j_beta=integral_T g_beta(y)dy.                             (3)
```

For a delay `L>=1`, form the physical packets first,

```text
F_(s,L)(x)=f(x)c_s(T^Lx)>=0,                               (4)
```

and then take the finite shift transform of their fixed physical-frequency
coefficients:

```text
A_(beta,L)(X)
 =1/13 sum_s (F_(s,L))_hat(X)zeta^(-beta s)

 =integral_T f(x)g_beta(T^Lx)exp(-2 pi i Xx)dx.             (5)
```

Write

```text
V_X=Var(f(x)exp(-2 pi i Xx)),

V_beta=Var(g_beta).                                        (6)
```

The complex BV covariance estimate used in THM-2478 gives

```text
|A_(beta,L)(X)-f_hat(X)j_beta|
 <=V_X V_beta/(12*13^L).                                  (7)
```

Consequently, for any finite set `B` of characters with `j_beta!=0`, every
delay satisfying

```text
13^L
 >max_(beta in B)
    V_X V_beta/[6|f_hat(X)j_beta|]                         (8)
```

obeys

```text
A_(beta,L)(X)!=0                         for all beta in B. (9)
```

Any finitely many already fixed source or word clocks may be added to the
lower bound on `L`.  Thus one common clock works simultaneously for all
twelve nonzero characters.

### Proof

Apply THM-2478's complex BV covariance estimate to

```text
b_X(x)=f(x)exp(-2 pi i Xx)
```

and `g_beta(13^Lx)`.  Its means are `f_hat(X)` and `j_beta`, which gives
(7).  Equation (8) makes the right side of (7) strictly less than half the
modulus of the limiting product.  This proves (9).  Linearity in the finite
`s` variable gives (5), and therefore the construction is Boolean before
DFT.  QED.

There is an exact frequency-graph form.  Cauchy--Schwarz makes the following
sum absolutely convergent:

```text
A_(beta,L)(X)
 =sum_(n in Z) f_hat(X-13^L n)(g_beta)_hat(n).              (10)
```

The `n=0` vertex is precisely `f_hat(X)j_beta`; equation (7) controls the
whole nonzero Koopman graph.  Unlike the full-`X` sum in THM-2568, (10) fixes
the output frequency `X` before summing the temporal sidebands.

## 2. The THM-2599 chamber has all twelve nonzero means

Fix one scalar-cover row.  THM-2604 chooses a pivot-eligible target-active
role `u` and THM-2599 pairs it with its actual `13`-divisible blocker `a`.
Fix a physical root `t` on which the THM-2537 selected source is positive.
THM-2599 gives a positive open chamber

```text
J subset I_t=[t/13,(t+1)/13)                               (11)
```

and an opposite-shift profile

```text
P_s(y)
 =d_(L_u)(u y+s/13) u_1(a y-s/13),

S={s:P_s|_J=1},                     1<=|S|<=4.              (12)
```

Put

```text
c_s=1_J P_s=1_(s in S)1_J.                                (13)
```

Then

```text
g_beta=j'_beta 1_J,

j'_beta=1/13 sum_(s in S)zeta^(-beta s),

j_beta=mu(J)j'_beta.                                      (14)
```

For every `beta!=0`, the sum in (14) is nonzero.  Otherwise the integer
`0/1` polynomial of `S` would be divisible by `Phi_13`; since `S` is nonempty
and proper, that is impossible.  THM-2599 gives the explicit invoice

```text
|j_beta|>=1/[13*4^11*182au].                              (15)
```

Thus Section 1 applies with

```text
B=F_13^*.                                                  (16)
```

## 3. The old marked triangle is literally unchanged

Let `e=1_(E_j)` be the old positive shallow owner carrier.  Refine the
THM-2537 selected source to one positive root piece

```text
0<f<=e,                  support(f) subset I_t,             (17)
```

retaining its terminal word and late owner.  As recorded in THM-2599
Section 8, this is still a rational Boolean subset of the old carrier.
Reapply THM-2349's abstract shallow-carrier lemma to `f<=e`.  For every
prescribed `kappa in F_13^*` it supplies integers `X,Y,m` with

```text
Y=X+m c_3,

gcd(m,91)=1,

X/13=Y/13=kappa                         mod 13,             (18)
```

and

```text
f_hat(X)
(1_(D_(c_3)))_hat(m c_3)
conjugate(e_hat(Y))!=0.                                  (19)
```

Choose one `L` satisfying (8) for the twelve functions (14), beyond all old
source and terminal-word clocks.  Combining (9) and (19) gives, for every
`beta!=0`,

```text
A_(beta,L)(X)
(1_(D_(c_3)))_hat(m c_3)
conjugate(e_hat(Y))!=0.                                  (20)
```

Equation (20) retains **the same**

```text
X, Y, m, kappa, old shallow carrier, deepest-comb edge,
terminal word, late owner, and source root t               (21)
```

while appending every nonzero future moving-shift character `beta`.
Nothing is reselected after `beta` is chosen.

Before the finite transform, each summand in (20) is the coefficient of the
physical nonnegative Boolean packet

```text
f(x)1_J(T^Lx)P_s(T^Lx).                                  (22)
```

On its support,

```text
floor(13x)=t=floor(13T^Lx),                               (23)
```

and the future point has the literal target-danger/blocker-safe event (12).
Thus (20) is a same-root chronological graft on one physical ancestry, not a
product of unrelated coefficient arrays.

## 4. Why this escapes THM-2568

THM-2568 considers a dangerous left factor and its safe right complement at
the **same physical point** under one common lawful twist.  Their product is
pointwise zero, so full physical-frequency recombination annihilates every
coarse target character.

Here the complementary roles, when present, are evaluated at `x` and
`T^Lx`.  There is no pointwise Boolean orthogonality.  The exact finite
control in the companion has

```text
sum_x P(x)(1-P(x))=0,

sum_x P(x)(1-P(Tx))=35>0.                                (24)
```

More importantly, (20) retains one fixed old `X` and sums only along the
Koopman graph (10).  It therefore realizes the first of THM-2568's three
analytic escape mechanisms--a fixed physical frequency before annihilating
full-`X` recombination--in a genuine delayed packet.  It does not yet supply
the complete two-ended atom required by that theorem's closing scope.

It does not contradict THM-2568.  Summing (20) over all old physical
frequencies, or replacing the delayed event by a same-time complementary
endpoint, restores that theorem's hypotheses and its zero.

## 5. The two action axes are not yet connected

The gain in (20) is coexistence of an old marked triangle coefficient and a
future moving-shift character.  Downstream address embedding can place the old
triangle in a relation-current expansion, but (20) is not yet their difference
character.

The reason is exact and factorwise.  For one old atomic coordinate
`y_i={w_i x}`, a lawful target co-shift acts by

```text
y_i -> y_i+theta_i/13.
```

In a factor evaluated at positive time,

```text
T^L(y_i+theta_i/13)
 =T^Ly_i+13^(L-1)theta_i
 =T^Ly_i.                                                  (25)
```

Thus the complete future stalk is neutral under the **old** target action.
No common translation of the physical variable `x` is asserted.  The
parameter `s` in (12) is a separately installed paired future local action,
and `beta` is its Fourier character.  There is no canonical physical map

```text
old relation residue  <-->  future shift character beta.  (26)
```

The same quotient kills the old root deck:

```text
T^L(x+q/13)=T^Lx                         for L>=1.          (27)
```

At the `C_91` scale it retains only the coarse septimal phase:

```text
T^L(x+m/91)=T^Lx+13^(L-1)m/7.                            (28)
```

Equations (25)--(28) explain the holotopy boundary.  The source and future
packets are two endpoint local systems over one temporal edge.  The edge
exists and carries a nonzero section (20), but its Koopman map is the quotient
which forgets precisely the `C_13` fibre needed to compare their charges.
A physical adjacent-clock section, ordered kernel, or other action
intertwiner is still required.

The valid implication is

```text
old fixed marked X,Y,m triangle
 + later same-root paired chamber with full primitive shift spectrum
 + complex BV mixing before finite DFT
 -> same X,Y,m survives in every future shift character.    (29)
```

The invalid promotions are

```text
future beta                  -/-> old relation residue;
same physical root t         -/-> transported C13 deck state;
source word/owner provenance -/-> future terminal endpoint;
nonzero two-time coefficient -/-> completed THM-2334 current;
coexistence                  -/-> scalar-row exclusion or LRC(14). (30)
```

## 6. Sharp controls and exact scope

Three hypotheses are load-bearing.

1. If the future `s` profile is flat, then `j_beta=0` for every
   `beta!=0`; no mixing argument creates a colour.
2. If the delay is not beyond (8), nonzero temporal sidebands in (10) may
   cancel its zero-frequency term.  The theorem asserts existence of one
   sufficiently large common clock, not every clock.
3. If one identifies `beta` with an old residue merely because both are
   elements of `F_13`, equations (25)--(27) give an immediate counterexample:
   the old deck acts trivially on the future factor.

What is closed is the exact-frequency/source-provenance loss in the
THM-2599 composition.  The remaining object is narrower: a same-carrier
identification of the future shift state with the next relation-address
section, including its adjacent clock and complete terminal endpoint.

No scalar profile is removed.  The live ledger remains `165`, and LRC(14)
remains open.

## 7. Exact companion

Run

```bash
python3 04-computation/lrc14_chronological_marked_triangle_graft_thm2610.py
python3 -O 04-computation/lrc14_chronological_marked_triangle_graft_thm2610.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_chronological_marked_triangle_graft_thm2610.out
```

The dependency-free exact companion verifies:

- all `1,092` nonempty proper profiles of sizes one through four and all
  `13,104` primitive cyclotomic character values;
- the `28,561` exponent selectors in the finite Koopman frequency-graph
  identity on `Z/13^2`;
- `104` exact positive-time root-deck erasures and `728` exact
  `C_91 -> C_7` carry rebases;
- a Boolean same-time complementary overlap of zero and chronological
  overlap of `35`.

The BV inequality (7) is the proved analytic estimate inherited from
THM-2478, not a floating-point experiment.  The companion uses no numerical
cyclotomic tolerance.
