---
id: THM-2618
title: "Selected-source product-torus Mobius lift and orbit-atom boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED AFTER SCOPE REPAIR.
  Relative to the retained labelled endpoint-factor presentation, the
  THM-2537 selected-source truth function has a unique Boolean Mobius lift
  to the independent coordinate torus.  The full THM-2350 dipole translation
  co-shifts every target-active factor occurring in that selected-source
  construction and the selector itself.  On every
  positive untwisted selected-source point its thirteen-state source profile
  is nonempty and proper; a positive rational BV chamber has one fixed such
  profile, so all twelve primitive source-translation colours are nonzero.
  Its profile atoms form a free C13 orbit, and their single-copy character
  transform is a lawful regular source-action sidecar.  The associated
  delta autocorrelation is instead a relative two-copy observer and is not a
  THM-2334 endpoint current.  The lift is canonical only relative to the
  labelled factor presentation: the relation-toric kernel gives strictly
  positive off-diagonal lifts with identical physical restriction and
  different target residues.  This is a lawful model and source-side input
  for THM-2573, not a closure of its type gate: the selector frozen there is
  THM-2569's distinct target-informed head, including a failure mask and a
  later occurrence.  The theorem proves no common fixed physical frequency,
  no nonzero Abel jump/tooth coefficient, no populated old present/bare
  Radon diagonal, no row exclusion, and no LRC(14).
source: lrc-source-axis-2026-07-28-product-torus-mobius-lift
depends_on:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2531-prime-necklace-guard-boundary-selector
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2502-endpoint-boolean-newton-carry-tournament-and-dipole-boundary
  - THM-2568-full-x-transition-annihilation-and-refined-pair-drift-boundary
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
  - THM-2573-logarithmic-abel-normal-and-common-endpoint-jump-pairing
  - THM-2574-oriented-tooth-component-holonomy-and-fixed-frequency-descent
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
  - THM-2613-canonical-root-diagonal-opposite-shift-section
  - THM-2615-physical-diagonal-toric-kernel-and-dipole-radon-invoice
script: 04-computation/lrc14_selected_source_product_torus_mobius_lift_thm2618.py
output: 05-knowledge/results/lrc14_selected_source_product_torus_mobius_lift_thm2618.out
script_sha256: 85ea8fffdd93ecbd8d33e3a914cb3b601cf3729f90dd947c0cc5e909e2513944
output_sha256: 2e55977cf0c0d131732f8b837ce5d7adbb12b12df99110b3e0d6c032458e0143
hash_basis: LF-normalized bytes
---

# THM-2618 -- the selected source has a lawful coordinate-torus orbit

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED AFTER SCOPE
REPAIR.**

THM-2573 isolates a precise type gate.  A target colour is lawful only after
**every** target-active factor, including a target-informed selector, has
been co-shifted.  Freezing a selector and moving a duplicate gate creates an
auxiliary colour of the kind excluded by THM-2568 and MISTAKE-266.

For the THM-2537 selected source, its own whole-factor orbit exists.  The
reason is that its selector is a finite Boolean function of the thirteen
complete root truth values.  Boolean Mobius inversion lifts this function to
the labelled independent-coordinate torus, where a THM-2350 dipole is an
honest translation.  One ordinary graft factor then forces the
source-translation profile to have a hole.  The result is a lawful
source-action sidecar at this frontier:

```text
labelled endpoint factors
  -> complete thirteen-root truth table
  -> Boolean Mobius selector lift
  -> full covariant dipole orbit
  -> nonempty proper source profile on a positive chamber
  -> every primitive source-translation colour.                 (1)
```

The construction is deliberately presentation-relative.  Restriction to the
one-dimensional physical curve forgets the relation lattice, so the scalar
truth function alone does not determine its off-diagonal lift.  It is also
selector-relative: THM-2569's later target-informed head is a different
Boolean event, not an alias for the selector lifted here.

Throughout, `p=13`, all equalities ignore the finite collection of strict
interval walls, and additions in root or target-state labels are in
`F_13`.

## 1. The labelled product torus and the physical root word

Let `E` be the nine labelled endpoint coordinates of one canonical live row,
let

```text
X=T^E,                         w=(w_i)_(i in E) in Z^E,

D_w(x)=(w_i x)_(i in E) in X.                              (2)
```

For each present endpoint role take its rational Boolean interval factor

```text
I_i:T->{0,1}.
```

Let `W^#` be the old terminal-word factor and let `g^#` be the sufficiently
late Boolean owner or owner--word factor.  At the live depth both are neutral
under the root translation used below.  They are also neutral under the old
target action: their dilation exponent is divisible by thirteen, so

```text
J(R(Z_i+s/p))=J(RZ_i),                  p|R.                (3)
```

It is convenient to encode this proved factorwise fact abstractly as

```text
W^#(Z+r w/p+s eta/p)=W^#(Z),
g^#(Z+r w/p+s eta/p)=g^#(Z)             for all r,s.        (4)
```

The target-neutral factors may equivalently be carried in an extra product
coordinate with trivial `C_13` action.

Put

```text
B^#(Z)=product_(i in E) I_i(Z_i),

e_r^#(Z)=W^#(Z) B^#(Z+r w/p),             r in F_13.       (5)
```

On the physical curve, with `x=z/13`, equation (5) is exactly the direct
inverse-root word

```text
e_r^#(D_w(z/13))=F((z+r)/13)                           (6)
```

from THM-2537 Section 7.  Thus blockers divisible by thirteen and the old
word are root-constant, while the guard/unit factors retain their labelled
root shifts.  No quotient or reconstructed mask is substituted for the
actual factors.

For a nonconstant Boolean word `e`, let `s_tau(e)` be THM-2531's canonical
occupied tail of the selected oriented `1 -> 0` wall.  For a fixed root
`rho`, define the selected-source truth table

```text
Phi_rho(e)
 =1_(e is nonconstant) 1_(rho=s_tau(e)).                    (7)
```

In particular,

```text
Phi_rho(e)=1  =>  e_rho=1.                                 (8)
```

The product-torus selected-source lift is

```text
S_rho^#(Z)=g^#(Z) Phi_rho((e_r^#(Z))_(r in F_13)).          (9)
```

Equation (6) makes the physical restriction of (9) exactly the `rho`-root
piece of THM-2537's `S^sel_tau`.

## 2. Exact Boolean Mobius lift and coordinate factorization

Every function on the Boolean thirteen-cube has a unique squarefree
multilinear expansion.  For (7), put

```text
c_(rho,A)
 =sum_(B subseteq A)(-1)^(|A|-|B|) Phi_rho(1_B),           (10)
```

where `1_B` is the indicator word of `B subseteq F_13`.  Boolean Mobius
inversion gives the identity

```text
Phi_rho(e)
 =sum_(A subseteq F_13)c_(rho,A) product_(r in A)e_r.      (11)
```

Here `c_(rho,empty)=Phi_rho(0)=0`.  Substituting (5), Boolean idempotence and
root neutrality of `W^#` give, for every nonempty `A`,

```text
product_(r in A)e_r^#(Z)
 =W^#(Z)
  product_(i in E) product_(r in A)I_i(Z_i+r w_i/p).       (12)
```

Thus every summand of (11) is again a product of one BV function in each
labelled coordinate.  Equations (9)--(12) are not an arbitrary extension of
one scalar function: they are the unique Boolean extension relative to the
retained labelled factor truth table.  A target twist therefore acts on
every occurrence of every active factor, including those used to decide the
selector.

The expansion is finite, so the Abel and absolute-rearrangement arguments of
THM-2334 apply term by term whenever the lifted packet is inserted into a
complete endpoint series.  Signed Mobius coefficients do not change the
pointwise Boolean value of (9).

## 3. Exact covariance under the owner-pivot dipole

Choose one THM-2350 target axis

```text
eta=e_a-e_(k_a),                                           (13)
```

where `a` is its blocker and `k_a` is the ordinary unit graft.  On `X` let

```text
T_s Z=Z+s eta/p,

S_(rho,s)^#(Z)=S_rho^#(T_s Z).                             (14)
```

The signs in (13)--(14) are a gauge convention; reversing every `s` gives
THM-2350's other displayed convention.  Directly from the definition,

```text
S_(rho,s)^#(T_t Z)=S_(rho,s+t)^#(Z).                       (15)
```

This is the whole-factor covariance law for the THM-2537 selected source.  It
is stronger than co-shifting only the `a,k_a` factors visible in one product:
the complete root truth table is recomputed after the coordinate translation,
and then the canonical wall is reselected from that translated table.  It is
not the still-missing covariance law for THM-2569's target-informed head.

For `alpha in F_13`, define the normalized source-action character

```text
Psi_(rho,alpha)(Z)
 =1/p sum_(s in F_13)S_(rho,s)^#(Z) zeta^(-alpha s),        (16)
```

where `zeta` is a primitive thirteenth root.  Reindexing in (15) gives

```text
Psi_(rho,alpha)(T_t Z)
 =zeta^(alpha t) Psi_(rho,alpha)(Z).                       (17)
```

Hence (16) is an actual character of the old coordinate action, not a root
character and not the separately installed future character of THM-2610.

## 4. Every positive source point has a nonempty proper profile

Suppose

```text
S_(rho,0)^#(Z)=1.                                         (18)
```

Then the selected root is occupied, so every factor of its complete present
endpoint is one.  In particular, the ordinary graft factor satisfies

```text
u_1(theta)=1,

theta=Z_(k_a)+rho w_(k_a)/p.                              (19)
```

Under (14), that factor becomes `u_1(theta-s/p)`.  A translated
thirteen-grid meets the open danger arc of length `1/7` in one or two points:

```text
#{s:d_1(theta-s/p)=1} in {1,2}.                           (20)
```

Choose one such `s`.  Then the `k_a` factor of the `rho` branch is zero, so

```text
e_rho^#(T_s Z)=0.
```

Equation (8) forces `S_(rho,s)^#(Z)=0`.  Consequently the source-state
profile

```text
L_(rho,Z)={s:S_(rho,s)^#(Z)=1}                            (21)
```

satisfies

```text
0 in L_(rho,Z),                 1<=|L_(rho,Z)|<=12.        (22)
```

No capacity estimate for the other twelve states is used.  The single
ordinary graft hole proves properness even if the translated selector moves
to another root.

THM-2537 proves positive mass of the untwisted selected source.  Decompose
that event first by its selected root `rho` and then by the finite profile
`L_(rho,Z)`.  All factors are rational BV step functions and the selector is
a finite lookup, so each profile atom is a rational finite union of
interval chambers on the physical base.  At least one pair `(rho,L)` has a
positive chamber `C` on which

```text
{s:S_(rho,s)^#=1}=L,

emptyset !=L!=F_13,                    0 in L.              (23)
```

For every `alpha!=0`, cyclotomic irreducibility gives

```text
sum_(s in L)zeta^(-alpha s)!=0.                            (24)
```

Indeed, a degree-at-most-twelve `0/1` polynomial vanishing at a primitive
thirteenth root is a multiple of `1+X+...+X^12`; its coefficients would all
be equal.  Equations (16) and (23)--(24) prove

```text
Psi_(rho,alpha) is not identically zero
                         for every alpha in F_13.           (25)
```

This is complete source-translation spectrum as a family of physical-base
BV functions.  It is not yet simultaneous nonvanishing of their means or of
one prescribed ordinary Fourier coefficient.

## 5. Profile atoms give the exact regular source torsor

Fix the positive profile `L` from (23), and form its complete Boolean atom

```text
K_L(Z)
 =product_(s in L)S_(rho,s)^#(Z)
  product_(s notin L)(1-S_(rho,s)^#(Z)).                   (26)
```

It is nonzero on `C`, and because `0 in L`,

```text
0<=K_L<=S_(rho,0)^#.                                      (27)
```

The physical pullback used here is explicitly

```text
kappa_(rho,L)(z)=K_L(D_w(z/13)).
```

Although its factors with `s!=0` are evaluations off the physical diagonal,
(27) says that this pullback is a positive rational Boolean subset of the old
THM-2537 carrier.  THM-2349's abstract shallow-carrier triangle may be
reapplied to this **base atom**.  This preserves the source word/owner and
supplies its usual marked `X,Y,m` triangle for each prescribed shallow root
colour.  Nothing here says that the same triangle survives every translated
atom.

From (15),

```text
K_L(T_t Z)=K_(L+t)(Z).                                    (28)
```

A nonempty proper subset of the prime cyclic group has trivial translation
stabilizer.  Hence the thirteen atoms in (28) are distinct and pairwise
orthogonal:

```text
K_(L+t)K_(L+u)=0                         for t!=u.          (29)
```

They are the regular `C_13` orbit of one lawful selected-source atom.  Their
character sections

```text
Theta_alpha
 =1/p sum_t K_(L+t) zeta^(-alpha t)                        (30)
```

are nonzero for every `alpha`; on `K_L=1`, equation (30) equals `1/p`.
This is a concrete thirteen-state source-action sidecar, not merely a count
of abstract equivariant identifications.

There is a tempting stronger identity.  Equations (28)--(29) give

```text
K_L(Z)K_L(T_t Z)=1_(t=0)K_L(Z).                            (31)
```

Its state DFT has every colour.  But (31) freezes the first, target-active
copy of `K_L` and translates only the second.  It is therefore a relative
two-copy autocorrelation.  It is **not** the one-copy orbit (28), and it is
not the common-action present/bare diagonal of THM-2334.

This distinction is exactly THM-2502's primitive-idempotent boundary.  At
one physical point the pattern atoms multiply diagonally, but a fixed-`X`
endpoint coefficient is a double integral over two endpoint variables.
Only after summing every physical frequency does Parseval recover a
pointwise product, at which point THM-2568's full-`X` boundary applies.
Thus (31) may be used as a Hilbert-energy or relative-pair observer, not as a
completed target current.

## 6. Comparison with THM-2573 and the exact remaining lift

THM-2573 Section 6 requires a family of complete layers

```text
(L_(s,h),R_(s,h))_(s in F_13)                             (32)
```

formed by co-shifting every target-active factor, including the selector and
the later occurrence, before Abel smoothing.  The frozen selector there is
**not** (7).  In THM-2565/2569 it is the target-informed head

```text
A=T_delta,

A_h(x)=A(x)1_(floor(13x)=h),

w_(N,h)(x)=A_h(x)A_h(T^N x),                              (33)
```

whose marker is selected from the separate `k_a` failure mask and its slope
stratum.  THM-2565 explicitly retains the old canonical word only as a
sidecar; the two selectors can differ.  Therefore equations (9)--(15) do
**not** unfreeze `w_(N,h)` and do not close THM-2573's type gate.

What they supply is an exact construction pattern and a lawful source-side
ingredient.  To close the gate one must give the complete event in (33) a
labelled factor presentation, include the old danger/safe facts and every
factor used to choose `delta`, apply the same Boolean-Mobius lift to that
whole truth table, and translate both `A_h(x)` and its genuinely later
occurrence `A_h(T^N x)` before forming either endpoint layer.  Only after that
construction may THM-2573's logarithmic Abel normal and target-state DFT be
read as lawful target observables.

Even after that missing lift, one would still have to prove the analytic
conclusions which could finish that route:

```text
the total-layer common-jump measure is nonzero;
its target-state profile is nonconstant;
an allowed deep frequency M sees it;
the connected tooth coefficient at r=M mod k is nonzero.
```

Those four clauses come strictly after the whole-head lift.  The immediate
proof target is (33) with full covariance; the bare uniform gate in THM-2574
remains a sharp zero-colour hostile.

## 7. Three sharp boundaries

### 7.1 The lift is factor-presentation-relative

Let

```text
Lambda(w)={a in Z^E:a.w=0}.
```

For `a in Lambda(w)`, the product-torus character
`chi_a(Z)=exp(2 pi i a.Z)` satisfies

```text
chi_a(D_w(x))=1                         for every x.        (34)
```

The THM-2350 owner-pivot data contain relations whose dipole residue is
nonzero modulo thirteen.  For example, with omitted unit `u_0`,

```text
a_0=w_a e_(u_0)-w_(u_0)e_a in Lambda(w),

delta_eta(a_0)=(a_0)_a-(a_0)_(k_a)=-w_(u_0)!=0 mod 13.    (35)
```

For any `epsilon>0`, the two strictly positive functions

```text
F_0(Z)=1,

F_1(Z)=1+epsilon(2-chi_(a_0)(Z)-chi_(-a_0)(Z))
      =1+4epsilon sin^2(pi a_0.Z)                          (36)
```

agree on the entire physical curve by (34), including every inverse root and
every temporal iterate, but have different nonzero dipole support off that
curve.  Therefore no physical truth function determines a unique
product-torus lift.  Equations (9)--(12) are canonical only because the
labelled factor presentation and selector truth table were retained.  This
is the toric ambiguity independently isolated in THM-2615.

### 7.2 Pointwise colours need not give scalar colours

The nonzero functions in (25) can have zero mean or vanish at a prescribed
ordinary frequency.  The finite hostile is a regular thirteen-point base
whose profile atoms `K_(L+t)` each have mass `1/13`.  Every point has one
nonempty proper translated profile and every pointwise source character is
nonzero, but

```text
integral Theta_alpha=0                    for alpha!=0.     (37)
```

Thus neither a common physical frequency `X` nor one nonzero scalar target
coefficient follows from (25).  A fixed-frequency theorem must be proved
after all word, owner, deep and bare factors are inserted.

### 7.3 The future shift is not the old bare endpoint

THM-2613 canonically fixes the physical root-to-**future-local-shift**
section.  It does not turn that future state into the independently
translated old bare endpoint of THM-2334.  The exact left-minus-right Radon
line therefore still requires a lawful present/bare square.

Even a populated nonnegative `13 x 13` square can have a blind diagonal.
The uniform identity permutation square has a positive entry on every
diagonal state but a constant diagonal profile, hence zero in all twelve
primitive Radon colours.  The unit-shift permutation has the same row and
column marginals and zero diagonal.  Separate source and future spectra or a
positive total diagonal do not decide the old target residue.

The valid new implication is only

```text
labelled selected-source factorization
  -> lawful full old-source dipole orbit with all pointwise colours
  -> typed input for a future complete present/bare construction.       (38)
```

No scalar row is excluded.  The live ledger remains `165`, and LRC(14)
remains **OPEN**.

## 8. Exact companion

Run

```text
python3 04-computation/lrc14_selected_source_product_torus_mobius_lift_thm2618.py
python3 -O 04-computation/lrc14_selected_source_product_torus_mobius_lift_thm2618.py
```

The dependency-free referee checks the complete `8,192`-word THM-2531
selector truth table, its Boolean Mobius transform and inverse, root
covariance, all `8,190` nonempty proper profile stabilizers, all `98,280`
primitive cyclotomic profile coefficients over an exact finite field, the
translated ordinary-graft hole count on the complete rational wall atlas,
the regular-orbit atom character law, the equal-mass scalar hostile, and the
two permutation-square diagonal controls.  Normal and optimized executions
must byte-match the stored transcript.

The cyclotomic certificate is computed in
`F_2[X]/(1+X+...+X^12)`.  This is a field because `2` has order `12`
modulo `13`.  Nonvanishing after reduction modulo two is a one-way exact
certificate of nonvanishing of the original integral cyclotomic sum; the
referee never infers complex nonvanishing from a numerical root sample.

The independent hostile audit rederived the physical restriction, Boolean
Mobius signs and idempotence, dipole/character signs, ordinary-graft hole,
cyclotomic function-level spectrum, profile-atom orbit, and autocorrelation
boundary.  It found and repaired one semantic overreach: the initial
candidate had conflated THM-2537's selected source with THM-2569's distinct
target-informed head and therefore incorrectly claimed to close THM-2573's
type gate.  Section 6 records the repaired scope and the exact missing event.

**QED.**
