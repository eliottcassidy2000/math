---
id: THM-4248
title: "Weight-eleven Z-zero owner-descent planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3992/3997/4012/4045/4218/4222/4232
  + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the inherited b=d=0
  reduced (2,3) seam, the complete exact-M=11 Z=0 wall with
  A*B*(A+B) nonzero contains no nonautomorphic planar Keller pair, for
  arbitrary U and every allowed lower coefficient. With THM-4232 this closes
  the whole exact-M=11 coefficient chart A*B*(A+B) nonzero. The three walls
  A=0, B=0, A+B=0, other cells, seam entry, JC(2), and DC(2) remain OPEN.
source: codex-planar-jacobian-breakthrough-20260826
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
  - THM-4222-dense-weight-eleven-primitive-cm-planar-jacobian-exclusion
  - THM-4232-weight-eleven-u-zero-primitive-cm-planar-jacobian-exclusion
related:
  - THM-4220-weight-ten-zeta-zero-genus-two-planar-jacobian-exclusion
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
mistake_firewall:
  - MISTAKE-487
  - MISTAKE-519
  - MISTAKE-522
scripts:
  - 04-computation/jc23_weight11_z0_owner_descent_exclusion_thm4248.py
  - 04-computation/jc23_weight11_z0_owner_descent_exclusion_thm4248_independent_audit.py
outputs:
  - 05-knowledge/results/jc23_weight11_z0_owner_descent_exclusion_thm4248.out
  - 05-knowledge/results/jc23_weight11_z0_owner_descent_exclusion_thm4248_independent_audit.out
script_sha256:
  - c99e41ad2a814b9a5aa92e3c5d717eb1c0c3a85fe893c2a528626ad51beddbfe
  - 760bd05bc6b5deea5df51f6b97d8bdc59c4710ba861f48c18ceec0e0385a223a
output_sha256:
  - 029bad30ae3c9ebe3cfdd12413ede27001f249f58da9c80c25bd699f7623da9d
  - b6c71aff2def9a0abb44d443b84fd0b1c27bda95be406d823802efcb90f8bc0d
hash_basis: raw LF bytes
audit: >
  PASS/ACCEPT. The primary executes 28,250 exact gates and the standard-
  library-only clean-room audit executes 27,997. Both cover all 26,624
  support/collision hostiles, retain the fixed residual point, and reproduce
  the five owner strata, faces, edges, regular models, CM and genus ledgers,
  degree-zero consequence, and sharp boundaries. Normal, optimized, and
  fixed-hash-seed streams match their frozen outputs.
---

# THM-4248 -- complete exact-weight-eleven `Z=0` owner descent

**PROVED RELATIVE TO THM-3992/3997/4012/4045/4218/4222/4232
+ VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam. Put

```text
A=[p^4y]H, B=[py^3]H, U=[p^5]H, Z=[y^3]H,
K=2848/45-(7/6)Delta.                                  (1)
```

> **Theorem.** The complete exact-weight-eleven wall
>
> ```text
> Z=0,                         A*B*(A+B) != 0           (2)
> ```
>
> contains no nonautomorphic planar Keller pair. The coefficient `U` and
> every allowed lower coefficient are arbitrary.

THM-4232 covers the complementary case `Z!=0`, already for arbitrary `U`.
Consequently:

> **Corollary (complete exact-`M=11` coefficient chart).** No
> nonautomorphic Keller pair exists in this seam when
>
> ```text
> A*B*(A+B) != 0,                U,Z arbitrary.         (3)
> ```

Thus the named coefficient residual is exactly the union of the three walls
`A=0`, `B=0`, and `A+B=0`. This is a statement inside the inherited seam;
it proves neither entry into the seam nor `JC(2)`.

The inheritance quartet is:

- closest proved mechanism: THM-4222/4232's primitive-CM degree-zero
  specialization;
- canonical hostile: THM-4218's hidden elliptic side face;
- corrected near miss: MISTAKE-487/522, requiring the complete weighted
  support and the fixed residual point `(2,0,1)`; and
- least-used sidecar: the affine identity for `K`, which makes the
  first-surviving-owner descent finite and prevents `Delta=K=0`.

## 2. Complete source and the five owner strata

The source is

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y
  +U p^5+Xi p^2y^2+A p^4y+B py^3,                    (4)
```

with no omitted monomial through weight eleven. For

```text
F_Q=(s^2-p)(1-QH)-Q s^2/2,       y=sp,                (5)
```

a row `p^i y^j` contributes the two valued endpoints

```text
(j+2,i+j,1),                     (j,i+j+1,1).          (6)
```

The nonoptional term `-Qs^2/2` contributes `(2,0,1)` and is retained in
both exact implementations.

There are exactly five mutually exclusive first-owner strata:

| stratum | surviving faces | base exponent |
|:---|:---|---:|
| `U!=0, K!=0` | `M,V_5,T_K` | `330` |
| `U!=0, K=0` | `M,V_5` | `330` |
| `U=0, Delta!=0, K!=0` | `M,V_4,T_K` | `132` |
| `U=0, Delta!=0, K=0` | `M,V_4` | `132` |
| `U=0, Delta=0` | `M,V_3,T_K` | `132` |

The last row has `K=2848/45!=0`; when `K=0`, necessarily
`Delta=5696/105!=0`. Hence the table is exhaustive.

The planes are

```text
M:   nu=(r+2l-2)/11,
V_5: nu=(l-1)/5,
V_4: nu=(-r+l-1)/4,
V_3: nu=(-2r+l-1)/3,
T_K: nu=r-l/2-2.                                      (7)
```

The two certificates cross every optional lower-support subset with every
deletion pattern among the five aggregate points

```text
(2,3,1),(2,4,1),(2,5,1),(3,4,1),(3,5,1).             (8)
```

The exact counts by row of the table are

```text
8192, 8192, 4096, 4096, 2048;       total 26624,      (9)
```

with zero hull failure. This is a simultaneous support/collision census,
not five generic spot checks.

## 3. Faces, edge gates, and packets

Up to nonzero torus monomials, the face equations are

```text
g_M=(S^2-P)(1-A S P^5-B S^3P^4)=R*C,
g_T=S^2(1-K(SP)^2-B(SP)^3P),
g_V=P(-1+c_kP^k+A S P^5),                             (10)
```

where

```text
(k,c_k)=(5,U),(4,Delta),(3,-1376/135).                (11)
```

Both `T_K` and every `V_k` are rational: with `T_0=SP`, the tail core is
linear in `P`, while the vertical core is linear in `S`. The only
positive-genus normalization is the unchanged main component `C`.

The face Pick ledgers `(2Area,boundary,interior)` are

```text
M  =(33,5,15),       T_K=(2,4,0),
V_5=(5,7,0),         V_4=(4,6,0),       V_3=(3,5,0). (12)
```

The outer packet is

```text
K!=0: (10,10,5,5,2,2,1),      total 35,
K=0:  (10,10,7,5,1),          total 33,              (13)
```

and both have defect `28=2*15-2`. These packet totals are source-boundary
ledgers, not target-map degrees.

For `K!=0`, all six outer and two internal edge schemes reduce to

```text
X-1, 1-KX^2, -K-BX, (X-1)(AX+B),
A+c_kX, c_k-X^k, AX-1, 1-BX.                          (14)
```

When `K=0`, the tail face and its two schemes are deleted. Their exact
discriminants show that the active stratum hypotheses are sufficient: no
extra unit or coefficient wall is hidden in `(14)`.

The main branches meet at

```text
P=S^2,                         1-(A+B)S^11=0,          (15)
```

and their gradient determinant is `-11(A+B)S^10`. Thus `(2)` gives eleven
distinct transverse torus nodes.

## 4. Exact regular models and genus completeness

For `U!=0`, use `Q=sigma^330`; for `U=0`, use `Q=sigma^132`. The main
charts are respectively

```text
(s,p)=(sigma^-30 S,sigma^-60 P),
(s,p)=(sigma^-12 S,sigma^-24 P).                      (16)
```

After scaling, the local node equations are

```text
R_0 C_0=sigma^330/2       or       R_0 C_0=sigma^132/2. (17)
```

Hence the eleven paths are `A_329` or `A_131` resolutions. Every face
normal is primitive. The outer lifted gcds equal their planar edge lengths;
the internal `M--V_k` and optional `M--T_K` edges are primitive. Their exact
determinant-one fans insert only rational multiplicity-one components:

```text
M--V_5,V_4,V_3: 5,8,19 intermediates;
M--T_K:          74 at base 330, 29 at base 132.       (18)
```

The exact side charts specialize to `(10)`. The edge schemes cover every
compactified boundary point, and the torus derivative determinants on
`C,V_k,T_K` are nonzero under the active hypotheses. The audited toroidal
face/edge input inherited through THM-4045 therefore gives a proper regular
semistable model; no unlisted positive-genus component is discarded.

After contracting rational paths, a `K!=0` core has vertices `R,C,T,V`,
eleven `R--C` paths and one path from `C` to each side vertex. A `K=0`
core deletes `T` and its path. In both cases

```text
b_1=10,                    g_special=5+10=15,          (19)
```

matching every global Pick ledger in the five-stratum atlas.

## 5. Primitive CM and the degree-zero contradiction

The unique nonrational component is

```text
C: 1-A S P^5-B S^3P^4=0.                              (20)
```

With `x=A S P^5`, it has the cyclic presentation

```text
P^11=(B/A^3)x^3/(1-x).                                (21)
```

It is the same genus-five component proved and independently audited in
THM-4222/4232. Its branch residues are `(3,10,9)`, its regular-differential
CM type is

```text
{4,5,8,9,10},                                         (22)
```

and this set has trivial stabilizer in `(Z/11Z)^*`. Together with its
negative, it contains every nontrivial character exactly once. The cited
primitive-CM theorem used in THM-4222 makes `J(C)` simple. Since its
dimension is five,

```text
Hom(J(C),E_0)=0                                       (23)
```

for the good elliptic special target `E_0:Y^2=X^3+1`. Therefore `C` maps
constantly to `E_0`; every other source component is rational and does so by
Riemann--Hurwitz.

The target scalings at bases `330` and `132` are exactly those of THM-4222
and THM-4232 and both give the smooth special fibre `E_0`. For a relative
generic Keller map, resolve the induced rational map from the regular source
surface to the proper target model by point blowups. Every new exceptional
component is rational, so its map to `E_0` is constant by Riemann--Hurwitz,
just like the already listed rational components. Pulling back a relative
degree-one target line bundle to the resolved model, proper-flat degree
conservation now gives

```text
deg(phi_generic)=sum_i m_i deg(phi_i^*L)=0.            (24)
```

A hypothetical Keller pair supplies a finite nonconstant generic morphism,
whose degree is positive. This contradicts `(24)` and proves the theorem.

## 6. Sharp boundary and scope firewall

The remaining three walls change the actual geometry:

- `A=0` creates a genus-five replacement face;
- `B=0` creates a genus-four face with an elliptic-quotient hostile; and
- `A+B=0` changes `(X-1)(AX+B)` to `A(X-1)^2` and coalesces the eleven
  main intersections at the boundary.

By contrast, `K=0` is harmless only after explicitly deleting `T_K`; it
changes the packet total from `35` to `33`. It is not lawful to infer this
from the `K!=0` model by continuity.

The theorem proves the `Z=0` exclusion and the union corollary `(3)`. It
does not prove any of the three remaining walls, another source-fibre cell,
entry into the reduced seam, `JC(2)`, or `DC(2)`.

## 7. Reproduction

From the repository root:

```bash
python3 -B 04-computation/jc23_weight11_z0_owner_descent_exclusion_thm4248.py
python3 -B -O 04-computation/jc23_weight11_z0_owner_descent_exclusion_thm4248.py
PYTHONHASHSEED=4248 python3 -B \
  04-computation/jc23_weight11_z0_owner_descent_exclusion_thm4248.py

python3 -B \
  04-computation/jc23_weight11_z0_owner_descent_exclusion_thm4248_independent_audit.py
python3 -B -O \
  04-computation/jc23_weight11_z0_owner_descent_exclusion_thm4248_independent_audit.py
PYTHONHASHSEED=4248 python3 -B \
  04-computation/jc23_weight11_z0_owner_descent_exclusion_thm4248_independent_audit.py
```

The primary uses exact `Fraction` and SymPy arithmetic. The clean-room audit
uses only the Python standard library, imports neither the primary nor SymPy,
and reconstructs the hull and sparse-polynomial engines independently. Both
print the theorem's source-component and degree-zero consequence, not only
the injected wall assumptions. Their declared hashes are in the frontmatter.
**QED.**
