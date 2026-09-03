---
id: THM-4343
title: "Disjoint A23-cubic contact planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4292/4297/4327/4339 + VERIFIED-EXACT +
  TWICE INDEPENDENTLY AUDITED. In the inherited reduced (2,3), exact-weight-
  twelve seam, Z=beta_11=zeta_3=0, W=-U, and U*K!=0 are extinct. The cubic,
  top A23, and internal root schemes lie in three distinct open toric edge
  orbits; their modifications have disjoint formal support and commute after
  common base change. Deep A23 repetition forces the cubic normal jet B(P)
  to vanish, while the unique terminal A23 packet has squarefree cubic edge.
  All positive-genus components have positive Keller-form order. Only local
  lemmas, not the out-of-scope global conclusions of dependencies, are used.
  Cubic endpoint exits, seam entry, JC(2), and DC(2) are not claimed here.
source: root + A23/cubic and hostile-referee agents / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4292-lambda-zero-repeated-face-keller-extinction
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4339-clean-interior-cubic-edge-planar-jacobian-extinction
related:
  - THM-4337-zeta-zero-exact-weight-twelve-endpoint-wall-extinction
  - THM-4342-clean-cubic-zero-exit-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_thm4343.py
primary_output: 05-knowledge/results/jc2_m12_clean_cubic_a23_contact_extinction_thm4343.out
primary_script_sha256: a7cc281624c784f4dd9d20cf6f7b686be823de2fec370d801d795f14e5db13c6
primary_output_sha256: 11b52dd612da31e94ad2efe2b4c897a976b9558f4356b5f7b546cc77ba1c9e65
independent_audit_script: 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_independent_audit_thm4343.py
independent_audit_output: 05-knowledge/results/jc2_m12_clean_cubic_a23_contact_extinction_independent_audit_thm4343.out
independent_audit_script_sha256: aa5335ef766f0340c8d90ca0188c5dd0c641c3d2f30d9183fb0ca831abf2099a
independent_audit_output_sha256: 3fa6f7627f294a86546a932b8686573fd13d642f8e82a4139ac276218f38c236
hostile_referee_script: 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_hostile_referee_thm4343.py
hostile_referee_output: 05-knowledge/results/jc2_m12_clean_cubic_a23_contact_extinction_hostile_referee_thm4343.out
hostile_referee_script_sha256: 0699b21130d5db2b5731e987ffcc42c619fc5abb942b768d87c1e833ed5f53a0
hostile_referee_output_sha256: 63d3a261359c4007d2dc80a2d1db073a4d0a0698c23553b8c37842be615a46f7
hash_basis: raw LF bytes
audit: >
  PASS AFTER DEPENDENCY-SCOPE HARDENING. The primary reconstructs both local
  charts from all sixteen source rows, verifies disjoint toric owners, solves
  the shortened A23 ladder, and proves terminal cubic squarefreeness. An
  import-free Fraction audit repeats the row transform, edge incidence,
  terminal arithmetic, and genera. A third hostile path independently checks
  Pick data, endpoint units, same-label roots, the delayed correction, a
  mod-11 discriminant witness, all five genus ledgers, and differential-order
  positivity. The proof imports only local cubic and A23 lemmas and supplies
  the disjoint-support proper-flat gluing separately. Normal and optimized
  streams byte-match all three frozen outputs.
---

# THM-4343 -- Disjoint A23-cubic contact planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4292/4297/4327/4339 + VERIFIED-EXACT +
TWICE INDEPENDENTLY AUDITED. THE DISPLAYED `U+W=0` CLEAN-CUBIC GATE IS
EXTINCT. ENDPOINT EXITS, SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OUTSIDE
THIS THEOREM.**

## 1. Statement and constrained inheritance

Work over an algebraically closed field of characteristic zero in the
inherited reduced `(2,3)`, exact-weight-twelve seam. Impose

```text
Z=beta_11=zeta_3=0,        W=-U,        U*K!=0,
K=2848/45-(7/6)Delta.                                      (1)
```

All other allowed lower coefficients are arbitrary. No nonautomorphic
planar Keller pair lies on `(1)`. Together with THM-4339, this closes the
clean cubic gate for arbitrary `U+W` while `U*K*W!=0`.

Dependency use is deliberately local:

- THM-4339 supplies its exact cubic chart, double/triple local
  classification, local form orders, and finite conductor normalization,
  not its global theorem requiring `U+W!=0`;
- THM-4297 supplies Sections 4--5's local `A_23` valuation ladder, not its
  global theorem requiring `U*Z*D!=0`;
- THM-4230/4327 supply the exact source, good target, unchanged charts, and
  proper-flat degree interface;
- no twelve-point simple response packet is transported through the merged
  top root, as forbidden by MISTAKE-531.

The inheritance board was

```text
top edge owner | cubic edge owner | A23 ladder | cubic normal jet
horizontal conductor | delta/graph ledger | good form | proper-flat glue. (2)
```

## 2. Two exact charts from one source

Write the complete sixteen-row source as

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2,              e=-1376/135,

H=-3p+(8/3)p^2+ep^3+Delta p^4+upsilon_5p^5+Up^6
  +Ky^2+Phi p^2y+Theta py^2+eta p^3y+xi_10p^2y^2
  +alpha_11p^4y+beta_11py^3+Wp^3y^2+zeta_3y^3+Zy^4.      (3)
```

### 2.1 Cubic edge

With `Q=sigma^12`, `s=sigma^-6S`, `p=P`, and `delta=sigma^6`, put

```text
A(P)=K+Theta P+xi_10P^2-UP^3,
B(P)=P^3(Phi+eta P+alpha_11P^2),
C_0(P)=-3P+(8/3)P^2+eP^3+Delta P^4+upsilon_5P^5+UP^6.
```

After setting `b_c=1/S`, the exact reciprocal equation is

```text
F_c=(1-delta^2Pb_c^2)
    (b_c^2-P^2A-delta b_cB-delta^2b_c^2C_0)
    -delta^2b_c^2/2.                                      (4)
```

No source term is omitted.

### 2.2 Top contact

On the main face put `b_t=1/S`, `r=P/S^2`, and `t=sigma b_t`. Every source
row `p^i y^j` maps to `t^(12-2i-3j)r^(i+j)`. Direct expansion under `(1)`
gives

```text
Hhat=U(r^6-r^5)
 +t alpha_11r^5+t^2(upsilon_5r^5+xi_10r^4)+t^3eta r^4
 +t^4(Delta r^4+Theta r^3)+t^5Phi r^3
 +t^6(er^3+Kr^2)+(8/3)t^8r^2-3t^10r,

F_t=(1-r)(b_t^12-Hhat)-t^12/2.                            (5)
```

At `q=r-1`, the initial form is `Uq^2(1+q)^5`. Equivalently,

```text
Ctilde=b_t^12-Ur^5(r-1),
Ctilde|_(r=1)=b_t^12,        partial_r Ctilde(0,1)=-U.    (6)
```

Thus the two smooth branches meet to length twelve at a two-branch `A_23`
point. This is one delta-twelve contact, not twelve labelled response points.

## 3. Disjoint toric owners, including the common corner

The global polygon is

```text
Pi=((0,1),(2,0),(4,2),(4,5),(0,7)),
(2Area,boundary,interior)=(42,14,15).                    (7)
```

The three relevant root schemes belong to distinct one-dimensional toric
orbits:

```text
E_c=[(4,2),(4,5)]: A(P)=0,
E_t=[(4,5),(0,7)]: U(X-1)^2=0,
E_i=[(2,0),(4,5)]: 1-WX=0.                              (8)
```

Their endpoint coefficient pairs are `(K,W)`, `(U,U)`, and `(1,-W)`, all
pairs of units under `(1)`. Every root therefore lies in its relative torus,
away from the shared fixed vertex `(4,5)`. Distinct torus orbits are disjoint,
so the root schemes are pairwise disjoint scheme-theoretically, before and
after any base change.

The corner is not omitted: its primitive outer directions have determinant
two, so its fixed toric subdivision is rational and coefficient-independent.
It is neither a mixed centre nor a positive-genus carrier.

Numerical root equality is not incidence. For example,

```text
U=-1, W=1, K=-1, Theta=3, xi_10=-3
    gives A(P)=(P-1)^3.                                  (9)
```

The cubic root, top root, and internal root all have local coordinate `1`,
but the owner pairs `(E_c,1),(E_t,1),(E_i,1)` are three different points.
With `Phi=1,eta=alpha_11=0`, the cubic carrier is the genuine smooth elliptic
curve `Y^2-Y-X^3=0`. Thus `(9)` is a hostile positive-genus separation test.

## 4. Exhaustive cubic strata on the changed top wall

The lower-hull calculation behind THM-4339 does not use `U+W!=0`; that gate
only made the top edge squarefree. Under `(1)`, the main component is

```text
C: 1-UP^6+US^2P^5=0,
```

birational to `y^2=-Ux(x^6-U)`, smooth of genus three. The internal edge is
primitive and simple, and all other outer edges are squarefree.

Because `(4)` is unchanged, the local cubic alternatives remain exhaustive:

| cubic stratum | local object | order / normalization |
|---|---|---|
| squarefree | genus-one `T` | `16` |
| double, `B(a)!=0` | rational `A_11` bridge | `28` |
| double, `B(a)=0` | horizontal node | finite simultaneous normalization |
| triple, `B(a)!=0` | smooth `j=0` elliptic tail | `26` |
| triple, `ord_aB=1` | rational nodal cubic | `46` |
| triple, `ord_aB>=2` | horizontal cusp | finite simultaneous normalization |

This imports only THM-4339's local calculation. Section 3 proves that the
changed top contact cannot alter any entry.

## 5. The complete A23 ladder and delayed correction

Put `A_top=2U+W=U`. The exact identity

```text
Ur^6+Wr^5+Zr^4
 =(A_top/2)(r^6-r^4)-(W/2)r^4(r-1)^2                  (10)
```

shows that the difference from THM-4292's quadratic local model enters
`F_t` as

```text
-(W/2)r^4q^3.                                            (11)
```

After `q=t^6y` and division by `t^12`, it begins at `t^6`. Hence it changes
none of the moving critical coefficients through `t^4`. If those four
coefficients vanish, the terminal gap is `6(nu-s)` while the correction can
first act at `6(s+nu)`; their difference is `12s>0`.

The local THM-4297 ladder therefore applies with `A_top=U!=0`. For
`nu<s=v(sigma)`, cancellation ends in a squarefree quadratic; at `nu=s`,
the only repeated point enters `nu>s`; there the first moving critical value
or terminal `b_t^12q` term splits before `(11)`. The possible good-form
orders are

```text
6s+2nu, (5s+9nu)/2, 2s+4nu,
(3s+7nu)/2, s+3nu,                                      (12)
```

all strictly positive. No simple-root response packet is used.

## 6. Deep compatibility firewall

At `r=1`, the critical-value input is

```text
h(t)=alpha_11t+(upsilon_5+xi_10)t^2+eta t^3
    +(Delta+Theta)t^4+Phi t^5+c t^6+(8/3)t^8-3t^10,
c=e+K=7168/135-(7/6)Delta.                              (13)
```

Entering the repeated range `nu>s` forces the five lower coefficients to
vanish and the balanced relation `c^2+2U=0`. In particular,

```text
alpha_11=eta=Phi=0             ==>             B(P)=0.  (14)
```

Thus any concurrent cubic collision is a horizontal conductor. It can never
be a smoothed bridge, order-46 nodal carrier, or order-26 elliptic tail.

If the last two moving critical coefficients also vanish, the unique
terminal packet is

```text
c=5152/405,              Delta=4672/135,       K=1856/81,
U=-13271552/164025,       W=13271552/164025,
Theta=-4672/135,          xi_10=41216/1215,     alpha_11=eta=Phi=0. (15)
```

Its cubic is

```text
A(P)=64(207368P^3+86940P^2-88695P+58725)/164025,
Disc(A)=-3947729324424178958336/32688603759375!=0.       (16)
```

The integer numerator cubic also has discriminant `2 mod 11`. Therefore the
terminal `A_23` splitter cannot coexist with any cubic collision or endpoint
exit. Before terminality, `(14)` reduces every collision to the already
controlled horizontal rows.

## 7. Commuting modifications and the genus ledger

After one common finite DVR base change, let `Z_c` be the union of cubic
conductor/vertical centres in `E_c^o` and `Z_t` the top centre in `E_t^o`.
Section 3 gives `Z_c cap Z_t=empty` scheme-theoretically. The semilocal formal
completion along their union is therefore the product of the two completed
neighborhoods.

The finite birational normalization algebra at a cubic node or cusp glues to
the identity off `Z_c`, hence is the identity near `Z_t`. Conversely, the
local `A_23` modification is supported over `Z_t` and is the identity near
`Z_c`. The operations commute and are dominated by one proper model. Cubic
normalization is finite, regular, and DVR-torsion-free, hence flat; the local
`A_23` sequence is the proper-flat sequence audited in THM-4297 Sections
4--5. Common ramification multiplies positive orders but cannot mix the two
Newton problems. Actual special-fibre multiplicities are retained.

Replacing twelve simple `R--C` nodes by one `A_23` point preserves delta
twelve. The exact component ledgers are

```text
squarefree:          (0+3+1)+(12+1)-3+1       =15,
singular cubic raw:  (0+3+0)+(12+1+1)-3+1     =15,
horizontal norm:     (0+3+0)+(12+1)-3+1       =14,
rational bridge:     (0+3+0+0)+(12+1+2)-4+1   =15,
elliptic tail:       (0+3+0+1)+(12+1+1)-4+1   =15.       (17)
```

The positive-genus cubic-side orders are `9,16,26`; `(12)` handles every
positive-genus `A_23` tail. All are positive, while every other component is
rational. Thus every component maps constantly to the good elliptic target.
On the common proper-flat model,

```text
deg(phi_generic^*L)=sum_Gamma mu_Gamma deg(phi_Gamma^*L)=0, (18)
```

contradicting the positive generic response degree. This proves `(1)`.
**QED.**

## 8. Exact owner-address lesson

The lossless discrete label here is

```text
(toric edge orbit, local root coordinate, multiplicity, normal Jacobian jet).
                                                                    (19)
```

Mapping `(19)` merely to the natural or Stern--Brocot label `1` identifies
the three distinct events in `(9)`. Equality and multiplicity are preserved
only within one edge orbit; incidence in the toric surface is destroyed. A
tournament on the bare numbers is therefore cosmetic. The native finite
object is the cyclic boundary-incidence graph together with a labelled-root
torsor and normal jet on each edge.

The numerator of `(16)` happens to contain `23^3`; this has no `A_23`
geometric meaning, because coefficient rescaling changes that factorization
without changing the singularity. It is a checksum, not a connection.

The remaining endpoint problems are structurally different: `K=0` moves a
root through the finite endpoint and `W=0` drops the cubic degree at infinity.
THM-4342 closes the former; the latter is not inferred here.

## 9. Reproduction and exact scope

Run all three independent paths in normal and optimized modes and byte-match
the frozen outputs:

```bash
python3 -B 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_thm4343.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_thm4343.py
python3 -B 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_independent_audit_thm4343.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_independent_audit_thm4343.py
python3 -B 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_hostile_referee_thm4343.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_a23_contact_extinction_hostile_referee_thm4343.py
```

The primary reconstructs the full source; the other two paths import none of
its code, and the referee uses only elementary Fraction arithmetic. What is proved is exactly `(1)`,
relative to the inherited seam and the narrowly listed local dependency
interfaces. Endpoint exits and seam entry remain outside the conclusion.
