---
id: THM-2457
title: "Complete-atom root co-support graph and semantic-word hostile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On one
  common oriented physical C_13 chart, the 128 complete THM-2452
  atoms carry an exact directed co-support matrix M. A loop is
  equivalent to same-atom THM-2401 positive service; a directed edge
  gives service after a two-atom Boolean coarsening. Any positive
  semantic parent and any selected nonzero matched table admit a
  coarsening by at most three atoms, with retained co-support at least
  M_0/16384. For service mass M, every nonzero root colour survives,
  the joint charged energy is at least M^2/342732, and one joint
  coefficient has magnitude at least M/2028, sharply. No lower bound
  in terms of THM-2452 drift follows without an explicit root-image
  coupling. A live strict-profile rational hostile has identical
  complete local bits, owner, clock, graft, and deep labels but pure
  and fork THM-2305 terminal words. Coarsening need not preserve
  drift, an old exact address, or a unique literal word. No scalar row
  is excluded.
source: codex-2026-07-26-root-cosupport-graph
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
related:
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2397-clean-root-same-parent-charged-role-partition
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2435-top-blocker-essential-parent-and-punctured-stalk-carrier
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2448-right-endpoint-cospan-transition-atlas
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2456-two-root-replica-uniform-offset-boundary
script: 04-computation/lrc14_complete_atom_root_cosupport_graph_thm2457.py
output: 05-knowledge/results/lrc14_complete_atom_root_cosupport_graph_thm2457.out
script_sha256: ffbc7e9a9a3dc776f4d0a3feea4ac3959476b7761141fe89547b9e94d9fa69d7
output_sha256: 114c3de76f403dcbeb3dede249533036dcf581bda79fe8aad39f81dbc2a137bd
hash_basis: working-tree bytes (LF)
---

# THM-2457 -- complete-atom root co-support graph

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2452 proves that one complete Boolean atom can be copied from the
present endpoint to the bare endpoint without changing its whole
nonnegative table. That is an endpoint identity. It does not identify
the retained data with two semantic root packets in one physical
`C_13` chart.

This theorem isolates the exact missing coordinate. Once one common
oriented root chart and one disjoint semantic pair are supplied, the
complete atoms form a finite directed **co-support graph**. A loop is
exactly the same-atom positive route of THM-2401. An edge is exactly
what is needed after allowing a two-atom Boolean coarsening. The graph
also distinguishes three different levels of retention:

```text
loop
  -> exact literal atom + THM-2401 root service;

edge
  -> two-atom common filter + THM-2401 root service;

arbitrary edge plus a selected positive atom
  -> at most three atoms + positive matched table + root service,
     but no promised drift or old exact address.                 (1)
```

The common chart is an additional hypothesis, not a consequence of
THM-2452.

## 1. The directed 128 by 128 co-support matrix

Put

```text
G=F_13.
```

Fix one lawful target/deep twist of THM-2452. Suppose that a rational
finite-interval base `Y` and one common oriented physical chart have
been identified:

```text
phi:Y x G -> T,

phi(y,u)=(y+u)/13.                                   (2)
```

A common orientation reversal is harmless: it replaces every root
label `u` by `-u`. Separate orientations on the two packets are not
allowed.

Let

```text
A_0,F_0:Y x G->{0,1}
```

be rational Boolean step packets satisfying

```text
A_0(y,u)F_0(y,u)=0                                  (3)
```

almost everywhere. Define

```text
a_0(y)=sum_u A_0(y,u),

f_0(y)=sum_u F_0(y,u),

M_0=integral_Y a_0(y)f_0(y)dy>0.                   (4)
```

Let

```text
{P_omega:omega in Omega},       |Omega|=128,
```

be the complete THM-2452 atom partition. Pull it back through the
same chart:

```text
U_omega(y,u)=P_omega(phi(y,u)).                     (5)
```

The `U_omega` are a rational Boolean partition of `Y x G`. Put

```text
A_omega=A_0 U_omega,

F_omega=F_0 U_omega,

a_omega(y)=sum_u A_omega(y,u),

f_omega(y)=sum_u F_omega(y,u),                      (6)

M_(omega,nu)
 =integral_Y a_omega(y)f_nu(y)dy.                  (7)
```

The `128 by 128` matrix `M` is directed: the first coordinate records
the `A` atom and the second the `F` atom. Every entry is a
nonnegative rational number. Since the atoms partition every root,

```text
a_0=sum_omega a_omega,

f_0=sum_nu f_nu,

M_0=sum_(omega,nu) M_(omega,nu).                    (8)
```

Call `omega -> nu` a directed co-support edge when
`M_(omega,nu)>0`. A loop is allowed.

This graph is not a tournament. Many pairs can be absent or can occur
in both directions, and the direction records semantic packet type,
not an orientation gauge.

## 2. Boolean unions and exact matched-endpoint addition

For any atom set `I subset Omega`, define

```text
U_I=sum_(omega in I)U_omega,

P_I=sum_(omega in I)P_omega.                        (9)
```

Orthogonality of the atoms makes both sums Boolean. Define the common
descendants

```text
A_I=A_0U_I=sum_(omega in I)A_omega,

F_I=F_0U_I=sum_(nu in I)F_nu.                      (10)
```

Their THM-2401 simultaneous mass is exactly

```text
M(I)
 =integral_Y
    (sum_u A_I(y,u))(sum_u F_I(y,u))dy

 =sum_(omega,nu in I)M_(omega,nu).                 (11)
```

Now return to the THM-2452 endpoint table. Write its present packet,
bare envelope, and deepest probe schematically as `F`, `E`, and
`Delta`, with

```text
FE=F.
```

For a complete atom, let

```text
T_omega(z)=integral_T F P_omega E P_omega Delta_z,

z in G^3.                                          (12)
```

For the union `I`, atom orthogonality and idempotence give

```text
T_I(z)
 =integral_T F P_I E P_I Delta_z

 =sum_(omega in I)T_omega(z).                      (13)
```

Consequently:

- `T_I` is again a lawful matched nonnegative table;
- if `I` contains an atom `kappa` with `T_kappa` nonzero, then `T_I`
  is nonzero;
- every factor already fixed outside the complete local atom remains
  fixed.

Equation (13) is a whole-table identity. It does not assert
coefficientwise monotonicity after a Fourier transform.

## 3. Loop, edge, and three-atom service

Fix an atom `kappa`.

### 3.1 Exact-atom service

The complete atom `kappa` itself satisfies the positive branch of
THM-2401 if and only if

```text
M_(kappa,kappa)>0.                                  (14)
```

Indeed, `U_kappa` is one common rational Boolean filter and (3) is
hereditary. Its final simultaneous mass is exactly the left side of
(14). Thus a loop is the sharp atom-preserving service criterion.

### 3.2 Two-atom service

If

```text
M_(omega,nu)>0,                                     (15)
```

then for `I={omega,nu}`, equation (11) gives

```text
M(I)>=M_(omega,nu)>0.                               (16)
```

THM-2401 therefore applies to the two-atom common filter `U_I`.
When `omega!=nu`, this route no longer names one complete literal
truth word.

In particular, if the selected positive atom `kappa` is incident to
an edge in either direction, adjoining the other endpoint gives
two-atom service while (13) retains a nonzero matched table.

### 3.3 Three-atom service around an arbitrary selected table

Equation (8) and (4) imply that some directed edge satisfies

```text
M_(omega,nu)>=M_0/128^2=M_0/16384.                 (17)
```

For an arbitrary atom `kappa` with `T_kappa` nonzero, take

```text
I={kappa,omega,nu}.                                 (18)
```

Repeated vertices are removed, so `|I|<=3`. Equations (11), (13), and
(17) give

```text
T_I nonzero,

M(I)>=M_0/16384.                                   (19)
```

Thus:

> **Three-atom common-root service.** Once the semantic disjoint pair,
> one common oriented physical chart, and its positive initial
> co-support are supplied, every selected positive THM-2452 atom can
> be enlarged by at most two complete atoms to a positive matched
> table satisfying THM-2401.

The factor `16384` is sharp for the matrix statement. Partition a
rational base into `128^2` equal intervals indexed by `(omega,nu)`.
On interval `(omega,nu)`, place the `A` singleton at root zero in
atom `omega` and the `F` singleton at root one in atom `nu`. The two
packets are disjoint, and

```text
M_(omega,nu)=1/16384
```

for every ordered pair.

## 4. A sharp positive-branch energy invoice

Let `I` be any atom set with

```text
M=M(I)>0.
```

Use the THM-2401 cyclic cross-correlation

```text
C(s)
 =integral_Y sum_u A_I(y,u+s)F_I(y,u)dy.            (20)
```

Then

```text
C(s)>=0,

C(0)=0,

sum_s C(s)=M.                                       (21)
```

Normalize the finite transform by

```text
C_hat(k)=(1/13)sum_s C(s)zeta^(-ks).
```

Parseval and Cauchy--Schwarz on the twelve nonzero shifts give

```text
sum_(k!=0)|C_hat(k)|^2
 = (1/13)sum_s C(s)^2-M^2/169

 >=M^2/(13*12)-M^2/169

 =M^2/2028.                                        (22)
```

Write the normalized root coefficients as in THM-2401:

```text
alpha_I(y,k)
 =(1/13)sum_u A_I(y,u)zeta^(-ku),

phi_I(y,k)
 =(1/13)sum_u F_I(y,u)zeta^(-ku),

J_I(k)
 =integral_Y alpha_I(y,k)conjugate(phi_I(y,k))dy.   (23)
```

THM-2401 equation (14) says

```text
C_hat(k)=13J_I(k).                                  (24)
```

Therefore

```text
sum_(k!=0)|J_I(k)|^2>=M^2/342732,                  (25)

max_(k!=0)|J_I(k)|>=M/2028.                        (26)
```

THM-2398, as already used in THM-2401, also implies

```text
J_I(k)!=0                    for every k!=0.        (27)
```

There is no denominator-free lower bound for every individual colour.
The denominator-sensitive norm floor of THM-2401 equation (17)
remains available.

The constants in (22), (25), and (26) are sharp. Take

```text
C(0)=0,

C(s)=M/12                 for s!=0.                (28)
```

This is realized by twelve rational base strata, with one displaced
`A` singleton and one fixed `F` singleton on each stratum. Then every
nonzero root coefficient has

```text
J_I(k)=-M/2028.                                    (29)
```

Because (2) is one common physical chart, THM-2401 Section 3 further
gives a common ordinary Fourier gauge for each root colour, with its
stated jump-count window.

## 5. When drift quantitatively feeds co-support

The graph theorem alone gives no comparison between THM-2452 drift
and `M(I)`. Here is an explicit sufficient sidecar.

For a Boolean matched table `T_kappa:G^3->[0,1]`, put

```text
bar(T_kappa)
 =(1/13^3)sum_z T_kappa(z).                         (30)
```

Let `D(T_kappa)` be its THM-2365 noncirculant drift. Orthogonal
projection and `0<=T_kappa<=1` give

```text
D(T_kappa)
 <=||T_kappa||_2^2
 <=bar(T_kappa).                                   (31)
```

Every normalized finite coefficient `B_kappa` also satisfies

```text
|B_kappa|<=bar(T_kappa).                            (32)
```

Choose a table entry `z_*` with

```text
T_kappa(z_*)>=bar(T_kappa),
```

and let `S_(kappa,z_*)` be its underlying Boolean physical support.
Let

```text
pi_13(x)={13x},

Omega_I={y in Y:
          sum_u A_I(y,u)>0 and sum_u F_I(y,u)>0}.   (33)
```

Assume the explicit root-image condition

```text
pi_13(S_(kappa,z_*)) subset Omega_I
```

up to null sets. Since multiplication by thirteen preserves Haar
measure and

```text
S_(kappa,z_*) subset
  pi_13^(-1)(pi_13(S_(kappa,z_*))),
```

one has

```text
mu(pi_13(S_(kappa,z_*)))>=mu(S_(kappa,z_*)).        (34)
```

The inclusion in (33)--(34) yields

```text
M(I)>=mu(Omega_I)
     >=T_kappa(z_*)
     >=bar(T_kappa)
     >=D(T_kappa).                                 (35)
```

This is the quantitative **root-image coupling**. It is stronger than
separate packet masses or averaged occupancy data.

For the global THM-2452 atom with

```text
D(T_kappa)>=D_0/16384,
```

equations (25)--(26) and (35) give

```text
M(I)>=D_0/16384,

sum_(k!=0)|J_I(k)|^2
 >=D_0^2/92001420705792,

max_(k!=0)|J_I(k)|
 >=D_0/33226752.                                   (36)
```

In the adaptive THM-2452 branch,

```text
|B_kappa|>=|B_alpha|/N_i,

N_i in (1,16,8,4,2,1),
```

the same root-image condition gives

```text
M(I)>=|B_alpha|/N_i,

max_(k!=0)|J_I(k)|
 >=|B_alpha|/(2028N_i).                            (37)
```

If the Boolean delayed-word table satisfies the explicit THM-2452
threshold and the root-image condition, then

```text
D_word>=mu(Q_term)^2D_0/65536
```

and hence

```text
M(I)>=mu(Q_term)^2D_0/65536,

max_(k!=0)|J_I(k)|
 >=mu(Q_term)^2D_0/132907008.                      (38)
```

No inequality of the form

```text
M(I)>=cD_0,          c>0 universal,                (39)
```

follows without a coupling such as (33). Split a rational base into a
large stratum carrying a fixed noncirculant matched table but only one
semantic packet, and a stratum of mass `epsilon` carrying both
semantic packets. The drift stays bounded away from zero while
co-support is `O(epsilon)`. This is the same distinction between
separate and simultaneous mass that is sharp in THM-2401.

## 6. What Boolean coarsening forgets

Equation (13) preserves nonnegative table mass, but the noncirculant
projection is signed:

```text
Q T_I=sum_(omega in I)Q T_omega.
```

The summands can cancel. Likewise, a fixed finite colour or fixed
absolute-address coefficient can cancel after atoms are added.
Therefore the two- and three-atom conclusions do **not** preserve:

- `D(T_kappa)>0`;
- the old finite target/deep colour;
- an old exact `(X,m)`;
- one unique complete literal atom;
- one semantic THM-2305 terminal word.

The result gives root/endpoint service. It is not by itself a
relation-current closure.

This theorem also does not invert the separate averaged
owner-minus-replica observable reserved in THM-2456. Its witness
`Omega_I` is nonlinear and fibrewise: it asks for simultaneous roots
on the same parent. It bypasses an averaged uniform-offset
cancellation because (21) pins `C(0)=0` while `sum C=M>0`; it does not
claim that all owner/replica occupancy densities are pointwise
extreme.

## 7. A live strict-profile semantic-word hostile

The common chart and semantic word cannot be read from the seven local
truth bits alone. Consider the strict blocker profile

```text
(lambda_j,lambda_a,lambda_3)=(1,2,5),

(c_j,c_a,c_3)=(13,338,1113879),

c_3=3*13^5,                                         (40)
```

with scalar speeds

```text
(H;q_1,...,q_5;c_j,c_a,c_3)
 =(9;7,8,10,11,12;13,338,1113879).                 (41)
```

Take

```text
q_*=7,

k_j=lambda_j+1=2,

R=13^2=169.                                         (42)
```

The nine speeds are positive and distinct with gcd one. Moreover

```text
nu_7(q_*)=1,

nu_7(c_3)=0<=nu_7(q_*),                             (43)
```

so this is in the live isolated-graft side of the frontier, not the
closed deep-`c_3` branch.

Define

```text
x_1=1041874/14480427,

x_2=135443621/1882455510.                           (44)
```

At both points:

- the guard and all five ordinary `q` roles are safe;
- `c_j` is dangerous;
- `c_a` and `c_3` are safe;
- the retained `q_*` graft is safe;
- the five split guard/unit roles other than `q_*` are all safe;
- the complete blocker bits are `c_j` dangerous and `c_a` safe.

Thus both points lie in the same source-owner set `E_j` and have the
same seven-bit complete local atom. The deepest labels `t=0,r=2`
also agree, because

```text
||c_3x_1||=2/13,

||c_3x_2||=261/1690=2/13+1/1690,

||c_3x_1-2/13||=0,

||c_3x_2-2/13||=1/1690.                             (45)
```

Now put

```text
y_i=Rx_i.
```

Both terminal points remain guard/unit safe, make `c_j` safe, and
make `c_a` dangerous. The decisive values are

```text
||c_jy_1||=496/6591>1/14,

||c_jy_2||=64481/856830>1/14,

||c_ay_1||=22/507<1/14,

||c_ay_2||=1429/32955<1/14,                         (46)

||c_3y_1||=0<1/14,

||c_3y_2||=1/10>1/14.                               (47)
```

Therefore THM-2305 assigns different intrinsic words:

```text
y_1 in Q_(j,{a,c_3})       -- a genuine fork,

y_2 in Q_(j,{a})           -- a pure transfer.      (48)
```

Every displayed membership is interior. The minimum exact margins
over all source, graft, deep-probe, and terminal inequalities are

```text
353/92274                 at x_1,

11476/2998905             at x_2.                  (49)
```

Small rational intervals around `x_1,x_2` therefore retain the same
two patterns. This is a positive-measure semantic hostile, not a seam
point.

The hostile proves:

```text
complete local atom + source owner + clock + q_* graft
  + deepest (r,t) labels

does not determine the THM-2305 terminal word.       (50)
```

If a word `Q_(j,sigma)` has already been fixed inside the packet
`F_theta`, then it separates the two strata. Thus (50) targets the
inference from a bare/local complete atom; it does not contradict a
branch in which the semantic word was independently retained.

## 8. Quotient-root and mixed-energy boundaries

THM-2435 and THM-2436 do not silently supply the hypothesis of
Section 1:

- THM-2435's translation-equivariant marked root is a derived
  nonlinear selector. At positive depth it also retains an ancestry
  kernel rather than a physical `C_13` section.
- THM-2436's mixed quotient character is an integrated signed
  observable, not two disjoint Boolean packets with positive
  same-parent co-support.

Their quotient colours and energies therefore cannot be substituted
for the matrix `M` without an additional physical ancestry/common-
chart sidecar. The already closed deep-`c_3` shapes are not reopened
by this statement.

## 9. Exact companion

Run

```text
python 04-computation/lrc14_complete_atom_root_cosupport_graph_thm2457.py
python -O 04-computation/lrc14_complete_atom_root_cosupport_graph_thm2457.py
```

The dependency-free exact companion uses only integer and `Fraction`
arithmetic and explicit `require` checks. It:

- verifies the `128 by 128` matrix invoice, loop, edge, and at-most-
  three-atom identities;
- constructs the sharp uniform `16384`-stratum matrix realization;
- checks separate-marginal and selected-atom hostiles;
- verifies the sharp charged-energy constants `1/2028`,
  `1/342732`, and the `M/2028` maximum;
- checks every denominator in (36)--(38);
- verifies positivity, distinctness, gcd and valuation typing of
  (40)--(43);
- checks every guard, ordinary, blocker, graft, deep-probe, and
  terminal inequality in (44)--(49); and
- confirms the pure/fork word split and positive exact margins.

Both transcripts must reproduce

```text
05-knowledge/results/lrc14_complete_atom_root_cosupport_graph_thm2457.out
```

byte-for-byte.

## 10. Scope

This theorem proves an exact conditional service interface:

```text
common oriented physical root chart
  + disjoint semantic packets
  + positive co-support graph edge
  -> THM-2401 root/endpoint service.                (51)
```

It does not construct the chart, prove a loop for the selected
THM-2452 atom, identify a semantic owner or word from local bits,
preserve the selected drift/address under coarsening, or convert the
resulting ordinary endpoint gauge into the canonical relation current.
No scalar profile is excluded, the ledger remains `165`, and LRC(14)
remains open.

An independent hostile audit rederived the directed-matrix
decomposition, sharp edge and root-energy constants, conditional
root-image coupling, split-base boundary, and every exact inequality
in the live pure/fork witness.  It independently reproduced the
normal, optimized, and stored transcripts and both LF hashes.

QED.
