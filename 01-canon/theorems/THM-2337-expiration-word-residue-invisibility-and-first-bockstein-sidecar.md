---
id: THM-2337
title: "Expiration-word residue invisibility and the first Bockstein sidecar"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. On every
  positive strict shallow-owner word stratum covered by THM-2327, an
  atomic terminal-word harmonic is multiplied by
  R=13^(lambda_j+1), so it is identically invisible to THM-2309's
  mod-13 target quotient. Retaining the base/word factor split exposes
  the first divided-difference jet beta mod 13. Its ordered two-target
  response is surjective and contains the pure-a, pure-b, and fork
  polarizers. Under THM-2325's primitive sharp septimal-support
  hypotheses, every nonzero target quotient and every such target jet
  contain at least 734,664,038,400,000,000 nonzero all-91-unit-address
  atomic terms. More strongly, the full semantic marked current, not only
  its zero-word-mode slice, has an Abel decomposition by target quotient
  and target jet in which every joint limit exists and at least one is
  nonzero. The surviving target may still be zero and its jet need not
  match the word-support mask. The full relation address alone cannot
  recover the jet: an exact atomic gauge changes beta while fixing the
  address, and it is not weight preserving. No scalar row is excluded and
  LRC(14) remains open.
source: codex-2026-07-25-expiration-word-bockstein
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-2327-two-colour-marked-unit-c3-triangle
  - THM-2331-two-sided-septimal-address-embedding-in-marked-current
  - THM-2333-abel-target-fibre-sum-landing-and-zero-fibre-boundary
related:
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
  - THM-2334-relation-residue-current-and-character-twist-pushforward
script: 04-computation/lrc14_expiration_word_bockstein_thm2337.py
output: 05-knowledge/results/lrc14_expiration_word_bockstein_thm2337.out
script_sha256: 910a2c4d4bc643d7ea83c5ae135f23ba5221c8c2824d4bed9bb9bba407610c72
output_sha256: 0d66033569526d98635f6fecb01ce55c8acdb75c5ad8d7128534fe7c35d211d3
hash_basis: working-tree bytes (LF)
---

# THM-2337 -- the word survives, but its target response is one layer late

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2333 groups a nonzero zero-word-mode slice by the mod-thirteen target
quotient and proves that some target fibre survives.  That operation
deliberately replaces the terminal word by its positive mass.  Concurrent
THM-2334 recovers the uncoloured full-word target pushforward.  The full
atomic presentation from the audited form of THM-2331 reveals a finer
joint target/word-jet statement and the exact reason the extra colour was
hidden:

```text
base endpoint mode u
 + 13^(lambda_j+1) * word mode beta
 + deepest-comb mode
 - bare endpoint mode v.                             (1)
```

The word mode is present, but the target quotient reduces (1) modulo
thirteen.  It therefore cannot see the word at all.  This is not an
estimate or a possible cancellation.  It is a zero response homomorphism.

The correct next coordinate is a **factor-coloured divided difference**.
It restores the word harmonic one thirteen-adic layer later, has full
two-target response image, and supports a full-semantic-current Abel
landing theorem.  It is not a function of the relation address after the
base/word factor split is forgotten.  This last qualification is
load-bearing.

## 1. The fully atomic marked current

Use one positive strict shallow-owner word stratum from THM-2327.  Write

```text
w=(H,q_1,...,q_5,c_1,c_2,c_3) in Z^9,
E=E_j,
Q=Q_(j,sigma),
sigma in {{a},{b},{a,b}},

k=lambda_j+1,
R=13^k,
W=T^(-k)Q,
E_Q=E intersection W.                              (2)
```

Here `a,b` are the two blocker labels other than `j`.  Up to null
endpoints, THM-2331 gives exact nine-factor presentations

```text
1_E(t)=prod_i chi_i(w_i t),
1_Q(t)=prod_i psi_i(w_i t),

1_(E_Q)(t)
 =prod_i chi_i(w_i t) prod_i psi_i(R w_i t).        (3)
```

Put

```text
a(u)=prod_i (chi_i)_hat(u_i),
b_Q(beta)=prod_i (psi_i)_hat(beta_i).               (4)
```

Let `i_3` be the labelled coordinate index of the deepest speed, so
`w_(i_3)=c_3`, and write `e_(i_3)` for its coordinate vector.

THM-2327 supplies integers `X,Y,m` with

```text
Y=X+m c_3,
gcd(m,91)=1,

(1_(E_Q))_hat(X)
 (1_(D_(c_3)))_hat(m c_3)
 conjugate((1_E)_hat(Y)) !=0.                       (5)
```

Both endpoints have positive owner grade, so

```text
13 divides X,Y,c_3.                                 (6)
```

An atomic term of (5) has indices `u,beta,v` satisfying

```text
(u+R beta).w=X,
v.w=Y.                                              (7)
```

Its exact labelled relation address is

```text
r_full=u+R beta+m e_(i_3)-v,
r_full.w=0.                                         (8)
```

Thus `r_full` lies in the integral relation lattice.  The word factor is
not an aggregate Fourier mode in (7)--(8): `beta` retains all nine
labelled interval-factor indices of the literal terminal word.

## 2. Residue invisibility and the divided-difference jet

Remove the transported word contribution only at the level of the
labelled address:

```text
r_base=u+m e_(i_3)-v.                               (9)
```

Then

```text
r_full-r_base=R beta=13^k beta.                     (10)
```

Consequently

```text
r_full congruent r_base mod 13^k,                  (11)
```

and in particular modulo thirteen.  Let

```text
K_13=(w mod 13)^perp,
G=K_13/L_13,
pi:K_13 -> G                                      (12)
```

be THM-2309's two-target quotient.  Equations (6)--(8) show that
`u,v,m e_(i_3)`, and `r_base` reduce into `K_13`.  Hence

```text
pi(r_full mod 13)=pi(r_base mod 13)
                 =pi(u)+pi(m e_(i_3))-pi(v).        (13)
```

The response of `beta` in `G` is identically zero.  Nonconstant word
modes can reweight a target fibre, but they cannot translate its
mod-thirteen label.

After retaining the decomposition (9)--(10), define the first jet

```text
J_k(r_full;r_base)
 :=(r_full-r_base)/13^k mod 13
 =beta mod 13 in F_13^9.                            (14)
```

This is a Bockstein-style divided difference, not an element of `G`.
In general `beta.w` is nonzero modulo thirteen, so writing `pi(beta)`
would be ill-typed.  The intrinsic two-target restriction is instead

```text
tau(beta)=(beta_a,beta_b) mod 13
          in B:=F_13^{ {a,b} }.                     (15)
```

The depth is sharp: a coordinate with `beta_i!=0 mod 13` is invisible
modulo `13^k` and first appears modulo `13^(k+1)`.  The CRT split also
records

```text
r_full-r_base
 congruent (-1)^k beta mod 7,                       (16)
```

because `13 congruent -1 mod 7`.  Thus, once the factor split is retained,
the septimal difference and (14) recover `beta mod 91`.

## 3. The first-jet response contains every polarizer

The zero response in (13) is not the end of the response-image audit.
The first jet has the largest possible two-target image, even after every
atomic Fourier factor is kept off its septimal zero set.

Assume in this section the primitive and sharp two-coordinate septimal
support hypotheses of THM-2325.  Fix any

```text
z=(z_a,z_b) in B.                                   (17)
```

Choose `beta in Z^9` coordinatewise by CRT so that

```text
beta_a=z_a mod 13,
beta_b=z_b mod 13,
beta_i=1 mod 7 for every i,
||beta||_infinity<=45.                              (18)
```

The other seven mod-thirteen residues may be fixed arbitrarily.  Every
word-factor coefficient in `b_Q(beta)` is nonzero by THM-2331's exact
support law.

Now fix a nonzero target `q in G` and any exact all-`91`-unit relation
`r` in the THM-2325 bank over `q`.  Set

```text
d=R beta+m e_(i_3)-r,
v=u+d.                                              (19)
```

Modulo seven, require

```text
u.w=X-R beta.w,
u_i!=0,
u_i+d_i!=0                  for every i.            (20)
```

Each coordinate forbids at most two residues.  Since the scalar word has
at least two septimally supported coordinates, THM-2331's sharp
two-pivot argument gives at least

```text
3*5^7=234,375                                     (21)
```

solutions for `u mod 7`.  Lift each one exactly with a Bezout vector
`z_w.w=1`:

```text
N=X-R beta.w,
h=(N-u_0.w)/7,
u=u_0+7h z_w,
v=u+d.                                              (22)
```

Then

```text
(u+R beta).w=X,
v.w=Y,
u+R beta+m e_(i_3)-v=r,                            (23)
```

and all `27` endpoint/word indices, as well as `m`, lie in the nonzero
Fourier support required by the atomic current.  The explicit height
invoice is

```text
||u||_infinity
 <=3+B(w)(|X|+R S ||beta||_infinity+3S),

||v||_infinity
 <=||u||_infinity
   +R||beta||_infinity+|m|+||r||_infinity.          (24)
```

THM-2325 supplies `3,134,566,563,840` relations in each nonzero target
fibre.  Therefore every pair

```text
(q,z) in (G minus {0}) x B
```

contains at least

```text
234,375 * 3,134,566,563,840
 =734,664,038,400,000,000                          (25)
```

nonzero atomic term/address pairs with an all-`91`-unit exact address.

In particular the response map `tau` is surjective and contains the
three semantic support polarizers

```text
p_{ {a} }=(1,0),
p_{ {b} }=(0,1),
p_{ {a,b} }=(1,1).                                 (26)
```

This is stronger than the assertion that the word merely moves: the
exact response image contains every prescribed target-side displacement.
It is still termwise, not a noncancellation theorem.

## 4. The full semantic current lands in a joint Abel fibre

The invisibility identity (13) upgrades THM-2333 from its separate
zero-word-mode slice to the actual nonzero current (5).

For `0<rho<1`, Abel-regularize all `28` atomic interval factors by their
base indices.  For `q in G` and `z in B`, let

```text
A_(q,z)(rho)
 =rho^|m| d_hat(m)
  sum a(u)b_Q(beta)conjugate(a(v))
       rho^(||u||_1+||beta||_1+||v||_1),            (27)
```

where the sum is over (7) together with

```text
pi(r_full)=q,
tau(beta)=z.                                        (28)
```

For fixed `rho`, absolute convergence follows from the separate Poisson
regularization of the eighteen marked-source factors, nine bare-source
factors, and the deepest comb.

Let `chi` range over the `169` characters of `G`, and let `eta` range
over the `169` characters of `B`.  Finite orthogonality gives

```text
A_(q,z)(rho)
 =1/169^2 sum_(chi,eta)
    conjugate(chi(q)) conjugate(eta(z))
    chi(pi(m e_(i_3)))
    M_X^(chi,eta)(rho)
    conjugate(F_Y^chi(rho)) rho^|m| d_hat(m),        (29)
```

where

```text
M_X^(chi,eta)(rho)
 =sum_((u+R beta).w=X)
   a(u)b_Q(beta)rho^(||u||_1+||beta||_1)
   chi(pi(u))eta(tau(beta)),

F_Y^chi(rho)
 =sum_(v.w=Y)a(v)rho^||v||_1 chi(pi(v)).            (30)
```

The signs and conjugations in (29) come directly from

```text
1_(pi(r)=q)
 =1/169 sum_chi conjugate(chi(q))chi(pi(r))         (31)
```

and (13).

To prove the limits, extend `chi after pi` from `K_13` to a linear
character of `F_13^9`.  If its coordinate vector is `s`, and the
coordinate vector of `eta after tau` is `t`, then `M_X^(chi,eta)` is the
`X`th coefficient of

```text
prod_i chi_(i,rho)(w_i x+s_i/13)
prod_i psi_(i,rho)(R w_i x+t_i/13),                 (32)
```

where `t` is supported only on `{a,b}`.  Likewise `F_Y^chi` is the `Y`th
coefficient of the first product in (32).  A different extension of
`chi` differs by a multiple of `w`; it is trivial on the indices in
(30) because

```text
u.w=X-R beta.w congruent0 mod 13,
v.w=Y congruent0 mod 13.                            (33)
```

Every shifted Poisson interval factor takes values in `[0,1]`, converges
in `L^1` to its shifted indicator, and composition with an integer circle
cover preserves the `L^1` norm.  The bounded-product telescope gives
`L^1` convergence of (32).  Hence every term on the finite right side of
(29), and therefore every `A_(q,z)(rho)`, has a limit as `rho` tends to
one from below.

Finally,

```text
sum_(q,z) A_(q,z)(rho)
 ->(1_(E_Q))_hat(X)d_hat(m)
   conjugate((1_E)_hat(Y)) !=0.                     (34)
```

Thus at least one joint target/jet fibre has a nonzero Abel limit while
the **literal terminal word `Q` remains in the current**.  The surviving
`q` can still be zero, and the surviving `z` need not be (26).

## 5. Why the jet is a sidecar rather than address data

There is an exact atomic gauge action

```text
(u,beta,v,m)
 ->(u-R gamma,beta+gamma,v,m),
gamma in Z^9.                                      (35)
```

It fixes

```text
u+R beta,
X,
r_full,
pi(r_full),                                        (36)
```

but changes

```text
r_base ->r_base-R gamma,
J_k ->J_k+gamma mod 13.                             (37)
```

Therefore neither the full exact relation address nor any of its residue
classes determines the first jet.  The base/word factor colours must
remain attached.

This loss persists inside the nonzero atomic support.  If `u_i` and
`beta_i` are seven-units, taking `gamma=7t e_i` preserves both support
conditions while cycling `beta_i mod 13` through all thirteen values.
The large bank in Section 3 shows the same phenomenon simultaneously
with every fixed all-`91`-unit address.

The gauge is not coefficient preserving.  For a nonzero danger or
complement interval mode, the one-coordinate `gamma=7e_i` weight ratio
at `rho=1` is

```text
 [(chi_i)_hat(u_i-7R)(psi_i)_hat(beta_i+7)]
 /[(chi_i)_hat(u_i)(psi_i)_hat(beta_i)]

 =u_i beta_i/((u_i-7R)(beta_i+7)).                  (38)
```

At `k=2`, `R=169`, and `u_i=beta_i=1`, this is

```text
-1/9456,                                           (39)
```

not one.  Address-gauge orbits therefore do not pair equal weights and
may reverse signs.  This is the precise obstruction to turning the
surjective response in Section 3 into automatic aggregate survival.

## 6. The ordered word code and the exact support mask

Although `beta` cannot translate `q`, the factor **role** of its two
target coordinates has a gauge-stable code.  Let `d` denote the danger
interval indicator.  At every nonzero supported index `n`,

```text
(1-d)_hat(n)=-d_hat(n),
d_hat(n)!=0 iff 7 does not divide n.                (40)
```

For `t in {a,b}` and `7` not dividing `beta_t`, define

```text
epsilon_t(Q)
 =(psi_t)_hat(beta_t)/d_hat(beta_t).                (41)
```

Then

```text
epsilon_t(Q)=+1 iff t belongs to sigma,
epsilon_t(Q)=-1 iff t does not belong to sigma.     (42)
```

The complete ordered code is

```text
sigma={a}:       (+1,-1),
sigma={b}:       (-1,+1),
sigma={a,b}:     (+1,+1).                           (43)
```

The zero modes give the equally injective ordered code
`(1,6),(6,1),(1,1)` after multiplication by seven.  Multiplying the two
entries would merge the two pure words; the ordered pair is essential.

Write `q=(q_a,q_b)` in THM-2309's canonical target coordinates and put

```text
Z_t(q)=1_(q_t=0)
      =1/13 sum_(lambda in F_13)
         exp(2*pi*i lambda q_t/13).                 (44)
```

The exact word-support masks are

```text
P_{ {a} }(q)=(1-Z_a(q))Z_b(q),
P_{ {b} }(q)=Z_a(q)(1-Z_b(q)),
P_{ {a,b} }(q)=(1-Z_a(q))(1-Z_b(q)).                (45)
```

Their supports have sizes `12,12,144`, exactly the two pure axes and
mixed locus of THM-2309.  Put

```text
A(q)=sum_(z in B)lim_(rho->1-)A_(q,z)(rho).          (46)
```

This is the full-word target aggregate from THM-2334.  The exact
word/support landing statistic is the nonnegative masked energy

```text
E_sigma=sum_(q in G)P_sigma(q)|A(q)|^2.             (47)
```

It satisfies

```text
E_sigma>0
 iff some surviving target q has
     support(q)=sigma.                              (48)
```

There is no cancellation between distinct desired target vectors in
(47).  Moreover it has a finite `169`-twist Gram formula.  Write

```text
A(q)=1/169 sum_(ell in G^)
       exp(-2*pi*i ell.q/13) H(ell)
```

for the uncoloured specialization of (29), and define the unnormalized
mask transform

```text
P_sigma_hat(s)
 =sum_(q in G)P_sigma(q)exp(2*pi*i s.q/13).
```

Then

```text
E_sigma
 =1/169^2 sum_(ell,ell' in G^)
   H(ell)conjugate(H(ell'))
   P_sigma_hat(ell'-ell).                           (49)
```

If

```text
h(x)=12 when x=0,
h(x)=-1 when x!=0,
```

the three exact kernels are

```text
P_{ {a} }_hat(s_a,s_b)=h(s_a),
P_{ {b} }_hat(s_a,s_b)=h(s_b),
P_{ {a,b} }_hat(s_a,s_b)=h(s_a)h(s_b).             (50)
```

Thus their only entries are `-12,-1,1,12,144`, and (49) is an explicit
positive-semidefinite quadratic form in the same shifted currents already
constructed in (32).

Because (44) is a finite character polynomial, the simpler linear test

```text
sum_q P_sigma(q) A(q)
```

is also an explicit finite linear combination of those currents.  Its
nonvanishing is sufficient but not necessary: different desired target
vectors can cancel.  Equations (47)--(50), not that linear test, give the
necessary-and-sufficient remaining word/support landing problem.
Equation (34) proves only that some unmasked fibre survives.  The
first-jet polarizer (26) and target mask (45) are termwise fully occupied
but not yet coupled after summation.

## 7. Response images, the coin checksum, and the zero boundary

THM-2225's fair-coin checksum succeeds because its cyclic response image
contains an additive antipode and the rotation action preserves weight.
The present response audit separates both ingredients:

```text
mod-13 address response of beta:       image {0};
first-jet target response tau(beta):   image F_13^2;
coefficient-preserving gauge action:   absent.       (51)
```

The second image contains every polarizer in (26), and for each fixed
`q` it contains the pointwise reversing displacement `-2q`.  But
`F_13^2` has odd order.  It has no nonzero translation `delta` with
`2delta=0`, so there is no fixed translation antipode.  The involution
`q -> -q` fixes the zero fibre, precisely the fibre left open by
THM-2333 and (34).

Surjective response is therefore necessary but not sufficient.  The
exact finite-group hostile remains sharp even after the target jet is
retained.  On

```text
H=G x B,                    |H|=13^4=28,561,
```

put

```text
U=delta_0+1_H,
V_0=delta_0-(1/(|H|+1))1_H
```

and translate `V_0` as in THM-2333.  Every joint fibre contains `|H|`
nonzero rational term pairs, yet the cross-correlation is `delta_0`.
Only the joint zero label survives.  This is an abstract group-algebra
hostile, not a claim about canonical LRC interval weights.  It proves
that the formal data

```text
full semantic word
 + full target and first-jet response image
 + termwise all-91-unit occupancy
 + nonzero total
```

still need a coefficient-sensitive coupling, such as positivity of
(47), a controlled gauge cocycle, or terminal-component phase.

## 8. Exact remaining boundary

The sharpened ledger is:

```text
every nonzero q and every target first jet
  has all-91-unit atomic addresses                 PROVED;

every joint (q,z) Abel limit for the full word
  exists, and some joint limit survives            PROVED;

word harmonic is invisible in q but appears in
  a surjective first Bockstein response             PROVED;

relation address alone loses that response          PROVED;

survivor lies in the word-support mask P_sigma       OPEN;
survivor has nonzero target q                        OPEN;
bounded visible/Jackson survivor                     OPEN;
terminal-component phase transport                   OPEN;
scalar-row exclusion                                 OPEN.              (52)
```

For THM-2334, the consequence is concrete: pushing the current only by
its full relation residue, even modulo `91`, cannot recover `beta`,
because the exact gauge (35) already fixes `r_full`.  The pushforward
must retain a factor-coloured jet or impose a canonical atomic gauge.

No scalar profile is excluded.  The exact ledger remains `165`, and
LRC(14) remains open.

## 9. Exact companion

The companion exhausts the sharp two-pivot septimal floor; checks word
invisibility, septimal visibility, and first-layer sharpness through
eight expiration depths; constructs all `169` target-jet CRT lifts; and
verifies the three ordered word codes and the `12,12,144` support masks.
It freezes the exact gauge-weight hostile `-1/9456`, proves that the
odd target group has no nonzero involutive translation, and verifies the
symbolic full-support zero-only hostile on the `28,561`-element joint
group.  Every load-bearing check raises explicitly under ordinary and
optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_expiration_word_bockstein_thm2337.py
python3 -O 04-computation/lrc14_expiration_word_bockstein_thm2337.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_expiration_word_bockstein_thm2337.out
```

byte-for-byte after LF normalization.
