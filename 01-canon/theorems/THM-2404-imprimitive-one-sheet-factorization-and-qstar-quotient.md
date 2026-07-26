---
id: THM-2404
title: "Intermediate cyclic quotient detects primitive one-sheet switches"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. Every
  finite-interval one-sheet section of T_ab factors canonically into
  one-sheet sections through T_a and T_b. An interior T_ab branch
  switch survives as signed boundary of the intermediate section
  exactly when its root displacement is nonzero modulo b; switches
  divisible by b cancel under the first boundary pushforward.
  Equivalently, all switches are imprimitive modulo b iff the
  intermediate inverse-root label is constant on each component of
  the base. On the THM-2396 equality face, Y=T_7(S) lies in D_u^c.
  The 11V components of H_0 have length 6/(91V), so an imprimitive
  section forces u<=91V, equivalently q_*<=c_3. Therefore q_*>c_3
  forces a primitive interior switch and THM-2402 types one endpoint
  by q_*. The threshold is geometrically sharp at u=91V. The signed
  C_7 boundary has no automatic C_13 target colour; a root-constant
  tensor is an exact hostile. No canonical owner/target transport, row
  exclusion, or LRC(14) proof is claimed.
source: codex-2026-07-26-imprimitive-one-sheet-factorization
depends_on:
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
  - THM-2402-orbit-disintegration-equality-and-signed-endpoint-pairing
related:
  - THM-2400-clean-parent-root-gauge-quotient-and-target-slope-boundary
  - THM-2403-clean-toothpick-unequal-slope-target-axis-imbalance
script: 04-computation/lrc14_imprimitive_section_factorization_thm2404.py
output: 05-knowledge/results/lrc14_imprimitive_section_factorization_thm2404.out
script_sha256: 9de876e6c426c996cb9bd365efa6e7aff6bb88742d321d7abcbeba357185d8c9
output_sha256: 8050b83efb63142ea56fd225a45eb4639ade3a9a7c8e41395dc264db03e41ea3
hash_basis: working-tree bytes (LF)
---

# THM-2404 -- the middle cyclic quotient detects primitive switches

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2402 leaves one sharply stated equality debt: force an interior
endpoint switch whose forty-nine-root displacement is primitive modulo
seven. The first guess was that absence of such a switch should force a
two-stage factorization. That guess is too weak: **every** one-sheet
section of a composite cyclic cover factors through the intermediate
cover.

The real invariant is the boundary of the intermediate section.
Imprimitive switches disappear under the first pushforward; primitive
switches survive as an opposite-sign endpoint pair. In the common-core
packet, the quotient also turns `q_*=7u` safety into ordinary `D_u`
safety. Exact component lengths then force a primitive switch throughout
the magnitude branch `q_*>c_3`.

## 1. Canonical factorization of every composite section

Let

```text
q=ab,                         a,b>=2,

T_n(x)=nx                     on R/Z.                  (1)
```

Let `E,H` be finite unions of circular intervals and suppose `E` is a
one-sheet section of `T_q` over `H`:

```text
E subset T_q^(-1)(H),

sum_(j=0)^(q-1)1_E((z+j)/q)=1_H(z)       almost everywhere.  (2)
```

Define the intermediate image

```text
Y=T_a(E).                                             (3)
```

Then

```text
E is a one-sheet T_a section over Y,

Y is a one-sheet T_b section over H.                  (4)
```

To prove the first assertion, two distinct `T_a`-roots in `E` over one
`y` would have the same `T_q=T_b o T_a` image and contradict (2).
Existence holds by the definition of `Y`. For the second, the unique
`x in E` over `z in H` supplies `y=T_a x`; two such `y` values would
again supply two `E` roots over `z`.

Thus the factorization itself carries no primitive-switch information.
It is universal and lossless almost everywhere.

The transfer and Fourier identities factor as well:

```text
L_a 1_E=1_Y,                    L_b 1_Y=1_H,

L_n f(y):=sum_(s=0)^(n-1) f((y+s)/n),

a 1_Ehat(an)=1_Yhat(n),

b 1_Yhat(bn)=1_Hhat(n).                              (5)
```

Composing (5) recovers `q 1_Ehat(qn)=1_Hhat(n)`.

## 2. Sheet indices and the exact primitive boundary

Let `I` be a connected component of the interior of `H`. Away from
finitely many switch points, the unique root of (2) has the form

```text
x(z)=(z+j(z))/(ab),                 j(z) in Z/(ab)Z.  (6)
```

The branch label `j(z)` is locally constant. Write

```text
j(z)=r(z)+b s(z),                   r(z) in Z/bZ.     (7)
```

Then the intermediate root is

```text
y(z)=T_a x(z)=(z+r(z))/b.                              (8)
```

At an interior switch `z_0`, let the left and right labels be
`j_-` and `j_+`. The two exposed endpoints of `E` have opposite signs,
and their `T_q` root displacement is

```text
m=j_--j_+                         modulo ab.           (9)
```

Call the switch **primitive modulo `b`** when `b` does not divide `m`.
By (7)--(8):

- if `b|m`, then `r_-=r_+`; the exit and entry map to the same point
  of `Y` and cancel;
- if `b` does not divide `m`, then `r_-` and `r_+` are distinct; the
  quotient has a genuine exit at `(z_0+r_-)/b` and entry at
  `(z_0+r_+)/b`.

Distributionally this is the first boundary pushforward in

```text
(T_a)_*(D1_E)=D1_Y,

(T_b)_*(D1_Y)=D1_H.                                  (10)
```

The following are therefore equivalent:

```text
(i)   every interior switch of E has b|m;

(ii)  D1_Y has no endpoint above the interior of H;

(iii) r(z) is constant on every connected component I of int(H);

(iv)  for every such I there is r_I in Z/bZ with
      Y intersection T_b^(-1)(I)=(I+r_I)/b
      up to endpoints.                               (11)
```

This is the corrected imprimitive factorization theorem. The quotient
`Y` always exists; condition (11) says exactly when its first-stage
branch is componentwise constant.

## 3. Exact primitive and imprimitive controls

Take `a=b=7` and

```text
H=[0,1/2).
```

An imprimitive switch is obtained from sheet labels

```text
9 -> 23                         at z=1/4.              (12)
```

Their displacement is `14`, while both have residue `2 mod 7`. The
section is

```text
E=[9/49,37/196) union [93/196,47/98),

Y=T_7(E)=[2/7,5/14).                                  (13)
```

The two `E` endpoints above `z=1/4` cancel after applying `T_7`; `Y`
has no interior endpoint.

Replacing the second sheet by `10` gives

```text
9 -> 10,

E=[9/49,37/196) union [41/196,3/14),

Y=[2/7,9/28) union [13/28,1/2).                       (14)
```

Now the displacement is `1`. The quotient has an exit at `9/28` and
entry at `13/28`, both above `z=1/4`. This is the primitive signed
boundary from (10).

The two examples have the same base `H` and exact one-sheet property.
Thus neither total mass nor the abstract two-stage factorization can
distinguish primitive switching.

## 4. The common-core intermediate section

Assume equality in THM-2396's floor:

```text
delta=mu(S)=66/4459.                                  (15)
```

THM-2402 says precisely that `S` is a one-sheet `T_49` section over

```text
H_0=D_V^c intersection D_(13V)^c,

q_*=7u,                    C_3=49V,                   c_3=637V.  (16)
```

Apply Section 1 with `a=b=7` and put

```text
Y=T_7(S).                                             (17)
```

Then

```text
S --T_7--> Y --T_7--> H_0
```

is a pair of one-sheet sections, with exact masses

```text
mu(Y)=mu(H_0)/7=66/637,

mu(S)=mu(Y)/7=66/4459.                                (18)
```

Moreover `S subset D_(q_*)^c`. Since

```text
1_(D_(7u))(x)=1_(D_u)(T_7x),
```

the quotient retains the top-word safety exactly:

```text
Y subset D_u^c.                                      (19)
```

This is the useful data preserved by the intermediate quotient.

## 5. Exact common-core component geometry

Normalize the base high-blocker core:

```text
Hbar=D_1^c intersection D_13^c.                      (20)
```

Its reduced half-open interval decomposition is

```text
Hbar
 =union_(r=1)^11
   [(14r+1)/182,(14r+13)/182).                       (21)
```

Thus it has exactly eleven components, each of length

```text
12/182=6/91,
```

and total mass `66/91`. Since

```text
H_0=T_V^(-1)(Hbar),
```

the actual base has

```text
11V components, each of length 6/(91V).              (22)
```

Suppose the equality section has no primitive interior switch. By
(11), each component `I` of `H_0` is carried by one literal `T_7`
root interval

```text
(I+r_I)/7
```

of length

```text
6/(637V).                                            (23)
```

But every connected component of `D_u^c` has length exactly

```text
6/(7u).                                              (24)
```

Containment (19) and connectedness force (23) not to exceed (24):

```text
6/(637V)<=6/(7u),

u<=91V,

q_*=7u<=637V=c_3.                                    (25)
```

The physical scalar labels are distinct. Hence the imprimitive equality
branch actually satisfies

```text
q_*<c_3.                                             (26)
```

Taking the contrapositive gives the new forced-switch branch:

```text
delta=66/4459                     and q_*>c_3

  => an interior primitive T_49 switch exists

  => at least one endpoint has q_* among its boundary labels.  (27)
```

The last implication is exactly THM-2402's endpoint congruence. It
does not need a new incidence search.

The constant `91` in (25) is geometrically sharp. For `V=1,u=91`,
every root lift of every component in (21) is exactly a safe component
of `D_91^c`:

```text
[(14r+1)/182,(14r+13)/182)+s
---------------------------------------
                   7

 =[(14k+1)/1274,(14k+13)/1274),

k=r+13s.                                             (28)
```

There are `11*7=77` such exact fits. This boundary has
`q_*=c_3` and is therefore an abstract sharpness control rather than a
valid distinct-speed LRC packet.

## 6. What the quotient does not transport to the target axis

The primitive object furnished by (27) is a signed `C_7` boundary pair

```text
beta=delta_(y_+)-delta_(y_-),

(T_7)_*beta=0,                 y_+!=y_-.              (29)
```

Its lift preserves:

- the common `T_49` base point;
- entry/exit orientation;
- a nonzero septimal root displacement; and
- at least one `q_*` boundary label.

The reserved THM-2403 target lane instead needs a nonconstant rational
function on a lawful `C_13` root/target orbit for one retained role and
its gate. There is no canonical map from (29) to that object:

```text
Hom(C_7,C_13)=0.                                     (30)
```

The quotient loses positive mass, the retained role and owner choice,
the thirteen-root coordinate, the unequal-slope target covector, and
terminal complex phase.

This failure is sharp. Tensor either primitive control (14) with a
root-constant rational `C_13` gate. The signed `C_7` boundary survives,
but all twelve nonzero target Fourier colours vanish.

The cheapest positive test is therefore mixed, not scalar. On a
physical primitive packet, form the lawful boundary table

```text
G(s)=sum_e sign(e) 1_(q_*-boundary)(e) W_e(s),

s in F_13,                                            (31)
```

where `W_e(s)` records the retained role/gate at the same endpoint under
the thirteen common-root translates. If `G` is rational and
nonconstant, cyclotomic irreducibility of `Phi_13` forces all twelve
nonzero Fourier colours of `G` to survive. One must then separately
audit the unequal-slope action and terminal endpoint typing.

Thus (27) supplies a signed endpoint anchor, not the missing canonical
common endpoint/root gauge.

## 7. Consequence and boundary

The equality face now has an exact dichotomy:

```text
q_*>c_3:
  forced primitive q_*-typed signed endpoint;

q_*<c_3:
  primitive endpoint may still occur, but component lengths alone
  do not force it.                                    (32)
```

If the THM-2396 floor is strict, THM-2402 instead supplies positive
two-sheet base mass; the present one-sheet branch classification does
not apply directly. No row is excluded, the scalar ledger remains
`165`, and LRC(14) remains open.

## 8. Exact companion

The dependency-free exact companion:

- verifies the two `T_7` one-sheet factors for both controls
  (13)--(14);
- checks that the displacement-`14` switch cancels and the
  displacement-`1` switch survives in `D1_Y`;
- reconstructs all eleven normalized common-core components directly
  from `D_1^c intersection D_13^c`;
- checks their lengths, mass, and all `77` aligned sharp fits (28);
- verifies the strict `u=92,V=1` obstruction; and
- retains the root-constant `C_13` target hostile.

Run

```bash
python3 04-computation/lrc14_imprimitive_section_factorization_thm2404.py
python3 -O 04-computation/lrc14_imprimitive_section_factorization_thm2404.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_imprimitive_section_factorization_thm2404.out
```

after LF normalization. Every executable assertion remains active
under optimized Python.
