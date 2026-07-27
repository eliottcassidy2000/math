---
id: THM-2540
title: "Weighted live-event Kakeya flux and transverse gain boundary refinement"
status: >
  PROVED + VERIFIED-EXACT.  For an external nonnegative density invariant
  under the thirteen root translations, the Cayley pairing of a Boolean
  root event is exactly its positively oriented boundary mass and half its
  translation energy.  The density is external: inserting a non-Boolean
  amplitude g inside the event produces |g|^2, not g.  Each algebraic root
  displacement has the sharp run-boundary floor 13/42 against centred
  fibre energy; on the unit-guard carrier it has the stronger sharp 1/10
  boundary-to-event mass floor.  On every one of the 165 live THM-2349
  rows, every one of the twelve boundary events literally retains the
  terminal word and common delayed Boolean owner, has a positive-frequency
  grouped-jump Prony lift bounded by 26J(W_R)-1, and re-enters the carrier
  triangle in every shallow colour with the stated finite invoice.  This is
  a 12 by 12 indexed family, not a family of distinct edges or a joint
  character theorem; only the guard-selected displacement has inherited
  physical orientation.  A separate product-torus statement converts the
  three already-proved THM-2533 gain seeds into one nonnegative rational sum
  of atomwise exclusive transverse boundaries retaining the same three
  72-element Galois orbits.  It does not add a fourth gain, preserve the
  deep anchor, select one common Boolean atom or ordinary-frequency gauge,
  or prove semantic arrival, row exclusion, or LRC(14).
source: codex-2026-07-27-weighted-live-event-kakeya-flux
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2525-unit-guard-collision-floor-and-word-owner-cross-cospan-collapse
  - THM-2529-deep-comb-adjacent-fibre-odd-consumer-and-zero-target-boundary
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2533-owner-weighted-phase-and-mixed-gain-radon-ladders
  - THM-2534-positive-kakeya-boundary-transform-and-crofton-reconstruction
related:
  - THM-2526-affine-skew-orientation-gauge-boundary
  - THM-2530-anchored-deep-gram-cone-and-lossless-skew-target-refinement
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
  - THM-2539-diagonal-cubic-owner-clock-boundary-current
script: 04-computation/lrc14_weighted_live_event_kakeya_flux_thm2540.py
output: 05-knowledge/results/lrc14_weighted_live_event_kakeya_flux_thm2540.out
independent_script: 04-computation/lrc14_transverse_gain_boundary_refinement_thm2540.py
independent_output: 05-knowledge/results/lrc14_transverse_gain_boundary_refinement_thm2540.out
script_sha256: 45dc32a4512efaa787f73b03d7c5474b8a8b8e262abe3e82fe18070228651dc3
output_sha256: 82f21452e4caf201b30f92252d5a031f7605ca37a0e85543e942136a156564db
independent_script_sha256: 864ee2eb212b2af715d7e5bf270a3b040db422963539717f98c857ae0f68cbbb
independent_output_sha256: 7ea5b301804b3b2e0e59a050a7c33316462ce3578616dfb3c2e8ec88295707a5
hash_basis: working-tree bytes (LF)
---

# THM-2540 -- weighted live-event Kakeya flux and transverse gain boundary refinement

**PROVED + VERIFIED-EXACT.**

Canonical THM-2534 already identifies

```text
K_tau(e)=e(1-P_tau e)                                      (1)
```

with the complete directed occupied-to-empty cut, proves its Crofton and
anchored Gram reconstruction, and gives the exact three-vector collapse on
the singleton/adjacent-pair deep cone.  THM-2533 already proves the
eventwise `72/5,184` root/cut bundles and at least three mixed gain orbits,
hence at least `216` mixed phase ladders.  None of those results is reproved
as a new conclusion here.

This strict companion records what survives after two extra operations:

```text
external root-invariant weighting
  -> exact positive Cayley flux and two quantitative boundary floors;

live owner-marked Boolean restriction
  -> twelve boundary subevents, each with its own Prony and carrier invoice;

transverse product translation
  -> one positive rational packing of the three existing gain seeds.    (2)
```

The distinctions between external density and internal amplitude, algebraic
root displacement and physical guard orientation, and spatial boundary and
chronological arrival are part of the theorem.

## 1. External-density Cayley flux

Let `T=R/Z`, and for `tau in F_13` use

```text
(P_tau f)(x)=f(x+tau/13).                                  (3)
```

Let `w in L^infinity(T)` satisfy

```text
w>=0,                    P_u w=w for every u in F_13.       (4)
```

This is an **external density**.  On the real weighted Hilbert space write

```text
<f,h>_w=integral_T w f h.
```

For `tau!=0`, let

```text
C_tau=P_tau-P_tau^2+P_tau^3-...+P_tau^11-P_tau^12.          (5)
```

THM-2532 gives

```text
(I+P_tau)C_tau=P_tau-I.                                    (6)
```

Root invariance of `w` makes `P_tau` unitary in `L^2(w)` and makes
`C_tau` skew-adjoint.  If `e:T->{0,1}` is measurable, put `y=C_tau e`.
Then

```text
<C_tau e,P_tau e>_w
 =<y,(P_tau-I)e>_w
 =<y,(I+P_tau)y>_w
 =1/2||(I+P_tau)y||_w^2
 =1/2||(P_tau-I)e||_w^2.                                  (7)
```

The first equality uses `<C_tau e,e>_w=0`; the third uses that `y` is real
and `P_tau` is unitary.  Booleanity and translation invariance of `w` give

```text
1/2||(P_tau-I)e||_w^2
 =integral_T w e(1-P_tau e)
 =integral_T w K_tau(e).                                   (8)
```

Consequently

```text
<C_tau e,P_tau e>_w
 =integral_T w K_tau(e)
 =1/2||(P_tau-I)e||_w^2 >=0.                               (9)
```

This is the exact weighted flux law.  It requires invariance under the
same translation used in the boundary.  For a noninvariant density, the
two oriented endpoint masses need not match, and (8) can fail.

### Density is not amplitude

The placement of the weight is load-bearing.  In the usual complex inner
product, if `g` is a complex root-invariant amplitude and `q=ge`, linearity
instead gives

```text
<C_tau q,P_tau q>
 =integral_T |g|^2 K_tau(e).                               (10)
```

For real nonnegative `g`, this is `g^2`, not `g`.  Moreover
`q(1-P_tau q)=ge(1-gP_tau e)` is not the Boolean boundary in (1) unless
`g` is Boolean.  Thus (9) has exactly two safe readings used below:

1. `e` is Boolean and `w` is an external root-invariant density; or
2. `g` is a root-invariant **Boolean** owner factor, so `g^2=g` and `ge`
   remains Boolean.

No form of (9) is applied to THM-2533's complex mixed character channel.

## 2. Two sharp directional floors

Disintegrate a root fibre by

```text
q_r(z)=e((z+r)/13),
m(z)=sum_r q_r(z),
p_r(z)=q_r(z)-m(z)/13.                                     (11)
```

For a fixed `tau!=0`, let

```text
k_tau(z)=sum_r q_r(z)(1-q_(r+tau)(z)).                      (12)
```

Canonical THM-2534 proves that `k_tau` is the number of occupied runs in
the `tau`-cyclic order.  Hence on a mixed fibre

```text
k_tau>=1,
sum_r p_r^2=m(13-m)/13<=42/13.                              (13)
```

On a constant fibre both sides below vanish.  Therefore

```text
k_tau >= (13/42)sum_r p_r^2.                               (14)
```

The constant is sharp.  Equality on a mixed fibre holds exactly when

```text
m in {6,7},                    k_tau=1,                     (15)
```

that is, when the occupied roots form one `tau`-consecutive block of size
six or seven.  There are `2*13*12=312` labelled mask--slope equality cases.
This is not THM-2527's different `Psi_tau` cut score and not its
four-necklace equality locus.

Because (4) makes `w` constant on each root fibre, integration gives the
weighted coercivity

```text
integral_T w K_tau(e)
 >=(13/42)||e-E_0e||_w^2.                                  (16)
```

Now assume in addition

```text
support(e) subset C_H={x:||Hx||>1/7},       13 does not divide H. (17)
```

THM-2525 proves `m(z)<=10` on every fibre.  Thus every active fibre obeys

```text
k_tau>=1>=m/10,                                             (18)
```

and consequently every algebraic displacement satisfies

```text
integral_T w K_tau(e) >= (1/10)integral_T w e.              (19)
```

If `wbar(z)` denotes the common value of `w` on the fibre, the exact gap is

```text
integral w K_tau(e)-(1/10)integral w e
 =1/130 integral_T wbar(z)[10k_tau(z)-m(z)]dz.              (20)
```

Hence equality holds exactly when `w`-almost every active fibre has ten
occupied roots in one `tau`-run.  The THM-2525 sharp `H=1` model with roots
`{2,...,11}` attains equality for `tau=1` and `tau=-1`.  Thus `1/10` is
sharp over the guard-supported class, but equality is not asserted for
every fixed displacement on one common extremizer.

## 3. Every live row has twelve owner-marked boundary subevents

Fix any one of the `165` live THM-2349 rows.  In THM-2533's notation, a
sufficiently delayed common Boolean owner--word block gives

```text
W_R(x)=F(x)G(13^(R+1)x),                                   (21)
```

where `F` is the nonzero rational Boolean terminal-word carrier and `W_R`
is a nonzero rational Boolean event with

```text
E_aW_R!=0                         for every a!=0.            (22)
```

THM-2349 has terminal clock `k>=2`.  On the root fibre, THM-2533 equation
(68g) reads

```text
W_R((z+r)/13)=A(z)e_r(z),                                  (23)

A(z)=Q_(j,sigma)(13^(k-1)z)G(13^Rz) in {0,1}.              (24)
```

The factor `A` is independent of `r`.  For every `tau!=0`, Boolean
idempotence therefore gives the literal factorization

```text
K_tau(W_R)((z+r)/13)
 =A(z)e_r(z)(1-e_(r+tau)(z)).                               (25)
```

Thus

```text
0<=K_tau(W_R)<=W_R<=F.                                     (26)
```

Every boundary tail retains the selected terminal word and the same delayed
owner.  Equation (22) implies that the set of nonconstant root fibres has
positive measure.  Since every nonzero `tau` generates `F_13`, (25) has at
least one occupied-to-empty wall on every such fibre.  Hence

```text
K_tau(W_R) is nonzero for every tau!=0.                     (27)
```

The translated head

```text
P_(-tau)K_tau(W_R)=(P_(-tau)W_R)(1-W_R)                    (28)
```

retains the common factors as well, but it marks an empty alternative at
the same horizon.  It is not a later arrival event.

All quantifiers are uniform over the finite live ledger.  For each row one
may use its THM-2533 delay threshold.  If one fixed Boolean block `G` is used
for all rows, the maximum of the `165` finite thresholds is one common
admissible delay.  No row-independent numerical value is claimed.

## 4. Per-boundary Prony and carrier invoices

Products and translations of BV step functions give

```text
J(K_tau(W_R))
 <=J(W_R)+J(P_tau W_R)
 =2J(W_R).                                                  (29)
```

For almost every root base `z`, the vector in (25) is either zero or a
nonempty proper rational Boolean mask.  If one of its nontrivial root
characters vanished, cyclotomic irreducibility of `Phi_13` would force all
thirteen coordinates to agree.  This is impossible on every active fibre.
Together with (27), this proves

```text
E_aK_tau(W_R)!=0                    for every a,tau!=0.      (30)
```

Apply THM-2533's grouped-jump Prony lemma to each boundary event.  For every
`a,tau!=0`, there is a positive integer

```text
n_(tau,a)=a mod 13,

1<=n_(tau,a)
 <=13L_a(K_tau(W_R))-1
 <=13J(K_tau(W_R))-1
 <=26J(W_R)-1,                                              (31)

(K_tau(W_R))hat(n_(tau,a))!=0.                              (32)
```

The positive frequencies and their gauge quotients need not agree across
`a`, `tau`, or rows.

There is a separate geometric invoice.  Reapply THM-2349's abstract
shallow-carrier triangle with

```text
e=F,                       f=K_tau(W_R).                    (33)
```

The hypotheses follow from (26).  For every `tau,kappa!=0`, it supplies a
boundary/terminal/future-owner-marked endpoint on a `91`-unit
deepest-comb carrier edge, with integer edge multiplier `m_(tau,kappa)` and

```text
gcd(m_(tau,kappa),91)=1,

|m_(tau,kappa)|
 <=13J(F)J(K_tau(W_R))+7J(F)-13
 <=26J(F)J(W_R)+7J(F)-13.                                  (34)
```

This gives a `12*12=144` **indexed** `tau`-by-`kappa` family on every live
row.  The statement does not say that its physical edges are distinct,
that their ordinary-frequency gauges agree, or that `tau` is a second
Fourier character jointly coupled to `kappa`.  The twelve `tau` values are
an algebraic root-displacement atlas.  Only THM-2526's guard-selected
`tau_H` has an inherited physical orientation.

The already-canonical pointwise `72/5,184` bundle is not integrated here.
In particular, canonical THM-2534's `132`-mixed-mode hostile forbids using
(30) to infer a general averaged joint-mode theorem.

## 5. One residual deep-consumer identity

No new Gram reconstruction is needed: canonical THM-2534 already proves
the target-anchor inverse and the complete deep path-cone collapse.  There
is, however, one useful scalar bookkeeping form not written there.

Let `A` be total singleton-ray mass, `B` total adjacent-pair-ray mass, and

```text
M=A+B,                    T=12A+22B,                         (35)
```

where `T` is total all-slope boundary mass.  Then

```text
B=(T-12M)/10,
A=(22M-T)/10,                                                (36)

7A+2B=13M-T/2.                                              (37)
```

Equation (37) rewrites THM-2529's positive odd consumer as the complement
of half the Crofton boundary mass.  It is an immediate consequence of the
canonical cone formulas, not a new lossless coordinate or a proof of
target-cell drift.

## 6. Transverse product boundaries positively pack the existing gains

This section is algebraic and deliberately unanchored.  Let

```text
B_(ell,s,alpha):T->{0,1},
ell in F_7,                   s in F_13,                     (38)
```

be a finite rational-BV Boolean atomization of the relevant nonnegative
table, and define its source--target character channels by

```text
W_(kappa,b)
 =1/91 sum_(ell,s,alpha)
   zeta_7^(kappa ell)zeta_13^(b s)B_(ell,s,alpha).           (39)
```

For `tau!=0` and a target co-shift `h`, set

```text
(P_(tau,h)B)_(ell,s,alpha)(x)
 =B_(ell,s+h,alpha)(x+tau/13).                              (40)
```

At ordinary frequency `n congruent a (mod 13)`, reindexing `s` and
translating `x` gives

```text
(P_(tau,h)W_(kappa,b))hat(n)
 =zeta_13^(a tau-bh) What_(kappa,b)(n).                     (41)
```

Thus the characteristic direction for gain `lambda=a/b` is

```text
h=lambda tau.                                               (42)
```

The tangent difference vanishes at that selected mode; for every transverse
choice the difference under `P_(tau,h)-I` is nonzero on that seed.

Use the three distinct gains `lambda_i=a_i/b_0` already proved in
THM-2533, together with their selected nonzero positive-frequency seeds.
For fixed `tau!=0`, exactly three of the twelve nonzero `h` are tangent, so
there are exactly nine common nonzero transverse choices.

For one Boolean atom write `A=P_(tau,h)B`.  Pointwise,

```text
A-B=A(1-B)-B(1-A).                                         (43)
```

The two terms are Boolean and mutually exclusive **within that atom**.
Sum the first and second orientations over the finite atomization, with
the original nonnegative rational coefficients, to obtain nonnegative
rational tables `U,V`.  Let `u_i,v_i` be their three selected linear
Fourier functionals.  Equations (41)--(43) imply

```text
u_i-v_i!=0,                         i=1,2,3.                 (44)
```

For one fixed `i`, the equation `u_i+t v_i=0` excludes at most one scalar
`t`.  Three rows cannot cover four candidates.  Therefore some

```text
t in {1,2,3,4}                                               (45)
```

makes the single nonnegative rational table

```text
Z=U+tV>=0                                                   (46)
```

retain all three seed functionals.

For the invariant nonvanishing statement use the squared projector norm

```text
Gamma_Z(kappa,b,a)=||E_a Z_(kappa,b)||_2^2.                 (47)
```

Rationality gives the exact Galois transport

```text
(kappa,b,a)->(v_7kappa,v_13b,v_13a).                        (48)
```

Each seed therefore fills its already-known `6*12=72` gain orbit, and the
three distinct ratios give the same `216` nonzero incidences as THM-2533.
The advance is not the count: it is the representation of all three seeds
on one nonnegative sum of atomwise exclusive transverse boundaries.

This positive packing cannot in general be replaced by one Boolean atom.
On the Boolean square, the three live functionals

```text
L_1(x)=x_1,             L_2(x)=x_2,
L_3(x)=x_1-x_2                                               (49)
```

have no common nonzero Boolean point, while the nonnegative vector `(2,1)`
makes all three values nonzero.  After summing different atom labels, no
mutual disjointness of `U` and `V` is claimed.

The shift (40) generally changes THM-2449's delta anchor and THM-2530's zero
row.  A target-dependent far square must be transported rather than treated
as common.  The weighted flux law (9) would require invariance under the
**full** product shift and is not asserted for `Z`.  This construction is
therefore not the anchored shift studied in THM-2538 and not a
same-ancestry arrival edge.

## 7. Exact boundary of the advance

The new statements remove four narrow invoices:

```text
external density:     exact positive Cayley flux, with weight type fixed;
guard carrier:        sharp 1/10 boundary/event mass floor;
165 live rows:        per-boundary Prony and indexed carrier certificates;
three mixed gains:    one nonnegative transverse boundary packing.      (50)
```

They do not supply:

```text
a common ordinary-frequency gauge across colours or displacements;
a joint tau/kappa character or 144 distinct physical edges;
a fourth mixed gain or more than the canonical 216 incidences;
one Boolean component retaining all three gains;
preservation of the deep anchor or zero row under product translation;
a later occupied endpoint on the same ancestry trajectory;
a source-to-arrival or blocker-to-owner loop;
a scalar row exclusion or LRC(14).                           (51)
```

In particular the `K_tau` head is a positive spatial absence certificate,
not a temporal handoff.  LRC(14) remains **OPEN**.

## 8. Exact referees

Run

```bash
python3 04-computation/lrc14_weighted_live_event_kakeya_flux_thm2540.py
python3 -O 04-computation/lrc14_weighted_live_event_kakeya_flux_thm2540.py

python3 04-computation/lrc14_transverse_gain_boundary_refinement_thm2540.py
python3 -O 04-computation/lrc14_transverse_gain_boundary_refinement_thm2540.py
```

The first dependency-free referee performs `2,721,357` exact checks.  It
exhausts every nonconstant Boolean root mask and every displacement, checks
the Cayley flux, run-boundary coercivity and its `312` equality cases, the
`1/10` floor on all `97,188` mass-at-most-ten mask--slope cases and its
`156` combinatorial equality cases, the sharp guard model, primitive
boundary colours, target-anchor regression controls, and the deep consumer
formula (37).  The canonical Crofton/Gram checks are retained only as
regression controls; they are not counted as new THM-2540 claims.

The independent product-torus referee performs `5,412,159` exact checks.  It
checks the complete characteristic multiplier, all `220` gain triples and
their nine common transverse choices, exact Boolean boundary splitting,
the four-scalar packing, Galois covariance, the `216` orbit census, the
Boolean-selection hostile, and `4,096` failures of the flux identity for a
noninvariant weight.  It imports THM-2533's live three-gain premise and
explicitly verifies neither anchor preservation nor semantic arrival.

Normal and optimized runs byte-match the stored transcripts.  **QED.**
