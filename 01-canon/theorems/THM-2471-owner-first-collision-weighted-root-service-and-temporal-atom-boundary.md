---
id: THM-2471
title: "Owner first-collision weighted root service and temporal-atom boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. At the
  first positive THM-2306 collision, the penultimate owner-normalized
  source and exact word-current canonically form disjoint rational
  weighted packets on one oriented C_13 root fibre. Their service is
  exactly 169 I_r. Every one of the twelve nonzero root colours is
  nonzero, their sum is -I_r, their squared energy is at least
  I_r^2/12, and some has modulus at least I_r/12. Each colour lands
  syndetically at a common original Fourier frequency of exact
  13-adic valuation lambda_j+r-1. A finite Boolean ancestry stalk
  realizes every weight. Refining after arrival records P_omega(T^k x),
  while refining the source records P_omega(x); the latter collapses
  to the unique owner loop. Thus a nonowner atom can restore aggregate
  endpoint drift only by leaving the semantic owner leg. The remaining
  bridge is an off-diagonal temporal endpoint intertwiner, not a
  positivity, rationality, or common-root existence problem. No scalar
  row is excluded and LRC(14) remains open.
source: codex-2026-07-27-owner-first-collision-weighted-service
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
related:
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
  - THM-2459-four-atom-drift-and-root-service-coarsening
  - THM-2466-delayed-word-simultaneous-drift-service-retention
script: 04-computation/lrc14_owner_first_collision_weighted_service_thm2471.py
output: 05-knowledge/results/lrc14_owner_first_collision_weighted_service_thm2471.out
script_sha256: 7112899f0d84198bca87d72e4f50ce9524382f13e9a72fbbcf5e426aabfca025
output_sha256: 4050c40b1bcf584ee340fc3ca6d6111518d65fd675f7ff710eff4a111f3824a9
independent_script: 04-computation/lrc14_first_collision_full_colour_landing_referee_thm2471.py
independent_output: 05-knowledge/results/lrc14_first_collision_full_colour_landing_referee_thm2471.out
independent_script_sha256: 638004272364368e90fabc66bd65c1b499a8ee4be72a64d8755a7a6350b34a7f
independent_output_sha256: 822bc642a8ef9fe4fb5ee3f4788ed3294fd7042a52d34a0aebe83228c3a3566e
hash_basis: working-tree bytes (LF)
---

# THM-2471 -- the first collision already carries every root colour

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2306 canonically normalizes an exclusive-owner handoff, waits until
the source and exact word-current first overlap, and proves that the
preceding `13`-adic Fourier shell has negative aggregate covariance.  The
open root-service discussion had treated that shell as only a spectral
witness.  It is more geometric: one step before collision, the two
densities are disjoint on a common circle, while their thirteen children
coexist over positive base mass.  Disintegrating that last step gives the
missing oriented root fibre without an external selector.

The packets are rational weights rather than Boolean indicators.  This is
not a defect: THM-2401 explicitly notes that its positive branch uses only
same-gauge disjointness, rationality, and positive service, and the weights
have an exact finite Boolean ancestry expansion.  Prime cyclotomicity then
fires all twelve colours.

The same construction also exposes the remaining semantic obstruction.
There are two inequivalent ways to insert a complete atom:

```text
arrival refinement:  (Q P^K 1_E) P_omega,

source refinement:   Q P^K(1_E P_omega).                    (1)
```

After transfer, the first reads `P_omega(T^K x)` and the second reads
`P_omega(x)`.  Their printed atom label is not an intertwiner.  The
source-refined graph is supported on the unique owner loop; any extra atom
that restores endpoint drift is therefore off the semantic owner leg.

## 1. The weighted positive-root lemma

Let `p=13`, let `m>=1`, and let `f,e:T->R_(>=0)` be nonzero rational-valued
step functions with rational breakpoints.  Put

```text
U=P_m f,
V=P_m e,                                             (2)
```

where

```text
P_a h(y)=(1/a)sum_(b=0)^(a-1) h((y+b)/a).
```

Assume

```text
UV=0                         almost everywhere,

I=integral_T (P_13 U)(P_13 V)>0.                    (3)
```

On the canonical oriented root chart

```text
phi(y,u)=(y+u)/13,                 u in F_13,        (4)
```

define

```text
A(y,u)=U(phi(y,u)),
F(y,u)=V(phi(y,u)),                                (5)

a(y)=sum_u A(y,u)=13 P_13 U(y),
b(y)=sum_u F(y,u)=13 P_13 V(y).                    (6)
```

Then

```text
A(y,u)F(y,u)=0,

M:=integral_T a(y)b(y)dy=169 I.                    (7)
```

Thus `(A,F)` is a rational nonnegative common-root pair with positive
service.  It need not be Boolean.

For `k in F_13`, use normalized root amplitudes

```text
alpha_k(y)=(1/13)sum_u A(y,u) zeta_13^(-ku),

phi_k(y)=(1/13)sum_u F(y,u) zeta_13^(-ku),

J(k)=integral_T alpha_k(y) conjugate(phi_k(y))dy.   (8)
```

Every nonzero colour survives:

```text
J(k)!=0                         for all k!=0.        (9)
```

Moreover the complete signed and energy ledgers are

```text
J(0)=I,

sum_(k!=0) J(k)=-I,                                 (10)

sum_(k!=0)|J(k)|^2>=I^2/12,

max_(k!=0)|J(k)|>=I/12.                             (11)
```

In fact some Hermitian colour pair satisfies

```text
Re J(k)<=-I/12.                                     (12)
```

### Proof

Equation (7) is immediate from (3), (5), and (6).  Define the
unnormalized cyclic cross-correlation

```text
C(s)=integral_T sum_u A(y,u+s)F(y,u)dy.             (13)
```

It is a rational nonnegative vector.  Equations (3) and (7) give

```text
C(0)=0,

sum_s C(s)=M=169I.                                  (14)
```

Hence `C` is nonzero and nonconstant.  If its Fourier transform vanished
at one `k!=0`, the rational polynomial `sum_s C(s)X^s` would vanish at a
primitive thirteenth root.  Divisibility by

```text
Phi_13(X)=1+X+...+X^12
```

would make all thirteen coefficients equal, contradicting (14).
Finite inversion gives

```text
C_hat(k)=13 J(k),                                   (15)
```

so (9) follows.

Root-character orthogonality gives

```text
sum_k J(k)
 =(1/13)integral_T sum_u A(y,u)F(y,u)dy
 =0.                                                (16)
```

Also `J(0)=M/169=I`, proving (10).  Cauchy--Schwarz on the twelve
nonzero colours proves the first inequality in (11), and the second
follows by taking a maximum.  Since real packets give
`J(13-k)=conjugate(J(k))`, taking real parts in (10) proves (12). QED.

This is precisely the weighted extension of the positive branch of
THM-2401.  No historical Boolean filtration is used.  Rationality is
load-bearing for the all-colour assertion (9), while (10)--(12) need only
nonnegativity and same-root disjointness.

## 2. Every colour lands on the original first-collision shell

The physicalizations of (5) through (4) are exactly `U` and `V`.  For a
fixed `k!=0`, equation (9) is therefore

```text
integral_T P_U(y,k) conjugate(P_V(y,k))dy!=0,        (17)
```

in the notation of THM-2323 and THM-2401.  Let `J_U,J_V` be the numbers
of nonzero jumps of `U,V`.  The endpoint-product Vandermonde argument in
those theorems applies to rational-valued step functions, not only
indicators.  For every integer `H`, some

```text
H<=h<=H+J_U J_V-1                                  (18)
```

satisfies

```text
U_hat(k+13h)V_hat(k+13h)!=0.                        (19)
```

Taking `H=0`, put `q=k+13h`.  Then

```text
1<=q<=13 J_U J_V-1,

13 does not divide q.                               (20)
```

The Perron multiplier identity

```text
(P_m h)_hat(q)=h_hat(mq)                            (21)
```

turns (19) into

```text
f_hat(mq)e_hat(mq)!=0.                              (22)
```

Thus *each* nonzero residue `k mod 13`, not merely some unit residue,
lands at a common original frequency.  The allowable `h` is syndetic:
every block of `J_UJ_V` consecutive gauges contains a landing.

## 3. THM-2306 specialization

Use the exact THM-2306 data

```text
E=E_j,
e=1_E,
c=c_j=13^lambda u,             13 does not divide u,

Q=Q_(j,sigma),
K=lambda+1,

f=1_Q P^K e.                                      (23)
```

The same proof works at any fixed positive-return clock allowed by
THM-2306 Section 6.  Let

```text
A=P_c f,
B=P_c e,

I_s=integral_T (P^s A)(P^s B),

r=min{s>=1:I_s>0}.                                 (24)
```

Put

```text
d=c 13^(r-1),

U=P_d f=P^(r-1)A,
V=P_d e=P^(r-1)B.                                  (25)
```

Minimality and nonnegativity give

```text
UV=0,

integral_T(P_13U)(P_13V)=I_r>0.                    (26)
```

Sections 1--2 apply with `(m,I)=(d,I_r)`.  In particular, for every
`k in F_13^*` there is a `q=k+13h` obeying (20) such that

```text
f_hat(dq)e_hat(dq)!=0,

nu_13(dq)=lambda+r-1.                              (27)
```

Equation (27) refines THM-2306 equations (25)--(31): its negative
aggregate shell is the sum of twelve nonzero charged root colours, one in
every unit residue, and every colour has an explicit finite landing
window.

The service retains all of THM-2306's semantic data already present in
`f`: the exclusive owner, prescribed return clock, and exact pure or fork
word.  It does not yet identify this owner-collision root with the
THM-2365 deepest-probe root.

## 4. The weights have an exact Boolean ancestry stalk

The rational packets in Section 1 are marginals of a finite Boolean
object.  This matters because replacing them by arbitrary real weights
would hide the LRC provenance.

Write

```text
R=13^K,

w_u=(y+u)/13,

X_(u,a)=(w_u+a)/d
       =(y+u+13a)/(13d),             a in Z/dZ,

Z_(u,a,b)=(X_(u,a)+b)/R,             b in Z/RZ,

Y_(v,e')=(w_v+e')/d,                 e' in Z/dZ.    (28)
```

For any rational Boolean partition `{P_omega}`, refine the *arrival*
current and source by

```text
f_omega=f P_omega,

e_nu=e P_nu.                                       (29)
```

Let

```text
U_omega=P_d f_omega,
V_nu=P_d e_nu                                      (30)
```

and form their root lifts.  The directed weighted co-support entry is

```text
M_(omega,nu)
 =1/(d^2 R) integral_T
   sum_(u,v,a,b,e')
    1_Q(X_(u,a))
    1_E(Z_(u,a,b))
    P_omega(X_(u,a))
    1_E(Y_(v,e'))
    P_nu(Y_(v,e')) dy.                             (31)
```

Every summand in (31) is Boolean.  The `a,e'` sheets retain the two
owner-normalization branches, `b` retains the prescribed ancestry, and
`u,v` are the collision roots.  Marginalizing these finite sheets produces
the rational weights in (30).  The entire `u=v` slice vanishes by (26),
and

```text
sum_(omega,nu)M_(omega,nu)=169I_r.                 (32)
```

Thus weighted service is not a relaxation of the physical problem.  It is
the exact shadow of a Boolean fibre product.  What is lost is the sheet
coordinate, not integrality.

## 5. Arrival atoms and source atoms are different temporal objects

The refinement (29) is an honest partition of the *current* function
`f`.  It is not automatically the endpoint atom used before the prescribed
clock.  Transfer adjointness makes the distinction exact.

For `P_I=sum_(omega in I)P_omega`, put

```text
W(x)=Q(Rx).
```

Arrival refinement gives

```text
(fP_I)_hat(L)
 =widehat(e (P_I composed_with T^K) W)(RL).         (33)
```

If instead the source is refined first, define

```text
e_omega=eP_omega,

f_omega^src=Q P^K(eP_omega).                       (34)
```

Then `sum_omega f_omega^src=f`, and

```text
(sum_(omega in I)f_omega^src)_hat(L)
 =widehat(e P_I W)(RL).                             (35)
```

Equation (33) reads the arrival-time atom `P_I(T^Kx)`; equation (35)
reads the source/endpoint atom `P_I(x)`.  Equality would require a proved
temporal atom intertwiner.  For a nontrivial Boolean atom, invariance
`P_I=P_I composed_with T^K` for a nontrivial Boolean union is impossible
except in the null/full cases by
THM-2461 Section 4.

For the untwisted complete THM-2452 bank, THM-2461 identifies one unique
source atom

```text
O=(all guard/unit safe, j dangerous, a safe)
```

with

```text
eP_O=e,

eP_omega=0                    for omega!=O.         (36)
```

Consequently the source-refined version of (31) has exactly one possible
entry:

```text
M_(O,O)=169I_r,

M_(omega,nu)=0                otherwise.            (37)
```

So the correctly endpoint-aligned service graph is a loop, not a rich
128-by-128 graph.  A rich graph obtained from (29) is meaningful arrival
data, but it lives at the temporal copy in (33).

## 6. A three-atom coarse coexistence lemma

There is still a useful purely algebraic consequence.  Let `Omega` have
`N>=2` atoms, let

```text
q_omega in H,

q_total=sum_omega q_omega!=0,                       (38)
```

and let `M_(omega,nu)>=0` be any directed matrix with total `M_0>0`.
Choose a largest entry `u->v` and let `S={u,v}`, with repetitions removed.
Write

```text
M(I)=sum_(omega,nu in I)M_(omega,nu).
```

For `J=Omega\S`, `r_*=|J|`, put

```text
z_0=q_S,

z_j=q_(S union {j})=q_S+q_j.                        (39)
```

If `r_*>0`, the exact identity

```text
q_total=sum_(j in J)z_j-(r_*-1)z_0                 (40)
```

shows that some `I=S` or `I=S union {j}` satisfies

```text
|I|<=3,

||q_I||>=||q_total||/(2r_*-1),

M(I)>=M_0/N^2.                                     (41)
```

For `N=128`, uniformly,

```text
||q_I||^2>=||q_total||^2/64009,

M(I)>=M_0/16384.                                   (42)
```

If the largest edge is nonloop, the drift denominator improves from
`64009=253^2` to `63001=251^2`.  Unlike THM-2459, no preselected atom is
retained, so three atoms suffice rather than four.

Applied to the arrival matrix (31), (42) puts nonzero aggregate endpoint
drift and owner-arrival service in one *label union*.  It does not identify
their endpoint legs.  Applied to the source loop (37), forcing `O` and at
most one other atom gives the same drift floor and retains the full service
`169I_r` algebraically.

The last statement must not be overread.  If `P_I` truly acts on the same
semantic owner leg throughout the endpoint table, (36) makes every added
nonowner atom invisible there; then `q_I=q_O`, so it cannot repair a dead
owner drift.  Conversely, if the added atom restores drift, that
contribution is off the semantic owner leg.  This is the exact cospan
dichotomy:

```text
same owner leg  -> added atoms cannot restore q_O;

restored drift  -> the restoring part is not owner-supported.     (43)
```

## 7. The minimal deep-root sidecar

The finite stalk (28) also quantifies when the collision root can descend
to a deepest probe.  Let `C` be the deepest speed and put

```text
delta=gcd(C,d),

C_0=C/delta,

d_0=d/delta.                                       (44)
```

On the current branch,

```text
C X_(u,a)
 = (C_0/d_0)((y+u)/13+a)             modulo 1.     (45)
```

Therefore:

- if `d_0>1`, the sheet residue `a mod d_0` is essential; no root map
  depending only on `u` exists;
- if `d_0=1`, write `h=C/d`.  When `13|h`, the deep phase is independent
  of `u`; when `13` does not divide `h`, the affine invariant is

  ```text
  theta=t-hu mod 13.                                (46)
  ```

For the strict-profile example `c_j=13`, `C=2*13^5`, and `d=13^r`:

```text
r<=4:   the deep root is base-only;

r=5:    theta=t-2u;

r>=6:   a mod 13^(r-5) is indispensable.           (47)
```

The physical endpoint first digit, owner-predecessor digit, collision
root, and deep root are therefore four typed coordinates on the full
stalk.  Divisibility can collapse some of them; equal cardinality cannot.

## 8. Exact stopping boundary

The proved chain is now

```text
exclusive owner + prescribed positive word
  -> first positive owner-normalized collision
  -> canonical weighted common C_13 root service
  -> every nonzero root colour
  -> common original frequency in every root residue
     at exact valuation lambda_j+r-1.               (48)
```

This removes the following former debts:

- existence and orientation of a common owner-collision root;
- positivity of its service;
- restoration of all twelve root colours;
- Boolean provenance of the rational weights; and
- finite landing in each unit residue.

It leaves a smaller, sharply typed obstruction:

```text
source atom P_omega(x)
  versus arrival atom P_omega(T^Kx)
  versus THM-2365 deep-table endpoint/root.         (49)
```

A closure theorem must provide an off-diagonal temporal intertwiner,
retain the full sheet coordinate from (28), or prove a polarized current
coupling the owner loop to the endpoint drift.  Boolean coarsening alone
cannot do this by (43).  Nothing here identifies the common shell
frequency with THM-2365's freshly selected relation address, proves the
owner-loop drift `q_O` nonzero, removes a scalar row, or proves LRC(14).

## 9. Exact companion

Run

```text
python3 04-computation/lrc14_owner_first_collision_weighted_service_thm2471.py
python3 -O 04-computation/lrc14_owner_first_collision_weighted_service_thm2471.py
```

The dependency-free `Fraction` companion:

- checks the weighted root/service identities on a non-Boolean finite
  ancestry stalk;
- verifies all twelve cyclotomic transforms by exact reduction modulo
  `Phi_13`;
- verifies (10)--(12) and the sharp uniform-correlation equality case;
- checks the arrival/source adjoint identities on finite cyclic models;
- exhausts the three-atom identity on hostile rational Hilbert vectors and
  verifies the `253`, `251`, `64009`, and `63001` constants;
- checks the deep-root descent cases in (44)--(47); and
- validates a strict-profile coordinate hostile in which owner digit,
  complete source atom, clock, word, and deep-root set agree while the
  endpoint first digit differs.

Both transcripts must reproduce

```text
05-knowledge/results/lrc14_owner_first_collision_weighted_service_thm2471.out
```

byte-for-byte.

The independently written companion named by `independent_script` tests the
landing statement on a hostile complementary Boolean pair over `65` cells.
All twelve root colours are nonzero, but the ordinary Fourier product is
nonzero precisely when `5|q` and `65` does not divide `q`; consequently ten of
the twelve base representatives `q=k` vanish.  Each progression
`q=k+13h`, `k!=0`, nevertheless lands exactly once in every five consecutive
gauges.  This separates a nonzero root colour from an immediate base-frequency
atom and confirms that the syndetic endpoint-Prony step in Section 2 is
load-bearing.  Normal and optimized runs reproduce

```text
05-knowledge/results/lrc14_first_collision_full_colour_landing_referee_thm2471.out
```

byte-for-byte.

An independent hostile audit rederived the first-collision disintegration,
all-colour cyclotomic argument, signed and energy invoices, finite Boolean
stalk, temporal adjoint correction, unique-owner loop, three-atom lemma,
and deep-root sidecar.  It specifically rejected the false identification
of `P_omega(x)` with `P_omega(T^Kx)` before this theorem was promoted.

QED.
