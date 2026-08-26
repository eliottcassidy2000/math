---
id: THM-4171
title: "Row-A inner-resultant wall planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/
  4155/4157 + VERIFIED-EXACT + INDEPENDENT SOURCE-PAIR/ENDPOINT AUDIT.
  The complete exact-weight-nine repeated-top row-A wall
  zeta=-eta, eta*Delta*(Delta+Theta)!=0,
  4Theta*K^2-27eta^2=0 contains no nonautomorphic planar Keller pair.
  Its four exhaustive source strata have critical lengths 22,21,20,19,
  all with genus 10 and packet (7,7,4,2,2,2,1). Full and finite
  cubic-carrier responses are strictly impossible. Other row-A walls,
  rows B/C, entry, other cells, M>=10, JC(2), and DC(2) remain OPEN.
source: codex-lrc14-jc-sharp-fronts-20260825b
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
  - THM-4157-repeated-top-edge-wall-planar-jacobian-exclusion
related:
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4165-y-only-inner-top-triple-intersection-planar-jacobian-exclusion
script: 04-computation/jc23_row_a_inner_resultant_wall_exclusion_thm4171.py
output: 05-knowledge/results/jc23_row_a_inner_resultant_wall_exclusion_thm4171.out
independent_audit_script: 04-computation/jc23_row_a_inner_resultant_wall_exclusion_thm4171_independent_endpoint_audit.py
independent_audit_output: 05-knowledge/results/jc23_row_a_inner_resultant_wall_exclusion_thm4171_independent_endpoint_audit.out
script_sha256: b55a485d9d393d1c01c55425e8416a50800f85cd94215e36cb0e47ff315ae2b4
output_sha256: cf837f54b2a1d200a7676ef7e5982f59094d6ab6519e9206a4721ad7a0f8fb34
independent_audit_script_sha256: 08671e089f20abe3a60c11ea3cc05f65546e8eaeb08f6d1585b32aa763cdb050
independent_audit_output_sha256: 3a9d86158a92bfec9de67242abb682ca3e3e38c401b6183161d1d172bacda0ba
semantic_sha256: 874df6fe1c3dc79c439ea46ebe70f6d05981f90adaa98064998784c72c747f0a
independent_semantic_sha256: 00d1b3b825527642dc602fb9c79fa116a1a1cbea09fc55fe83b9cf2874710556
hash_basis: raw LF bytes
primary_audit: >
  PASS. Seventy-six exact checks reconstruct the complete wall chart, the
  (A,B) source determinant, four actual endpoint rows, quartic coverage and
  unit firewalls, source and normalized infinity gates, all four universal
  critical points, the carrier discriminant, and both response ledgers.
  Normal, optimized, and hash-seeded outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room (A,C_0) Sylvester determinant independently reads the
  actual endpoint rows, derives the same four residual degrees, and replaces
  the primary gcd gates by terminal resultants. Its scope is an independent
  source-pair/endpoint audit, not a second normalized-resultant audit.
---

# THM-4171 -- row-A inner-resultant wall planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/
4155/4157 + VERIFIED-EXACT + INDEPENDENT SOURCE-PAIR/ENDPOINT AUDIT;
JC(2) REMAINS OPEN.**

Work over `C` in the live `b=d=0` reduced `(2,3)` seam at exact residual
weight nine. Use the repeated-top notation of THM-4157:

```text
zeta=-eta,                  K=2848/45-(7/6)Delta,
C=Delta+Theta,              D_A=4Theta K^2-27eta^2.
```

## 1. The theorem and inheritance pass

> **Theorem.** The complete coefficient locus
>
> ```text
> zeta=-eta,       eta*Delta*C != 0,       D_A=0               (1)
> ```
>
> contains no nonautomorphic planar Keller pair.

There are four exhaustive and disjoint source strata on `(1)`. Their exact
affine critical lengths are:

| stratum | `L` |
|---|---:|
| `J!=0` | `22` |
| `J=0, S!=0` | `21` |
| `J=S=0, P!=0` | `20` |
| `J=S=P=0` | `19` |

Every row retains normalization genus `10` and the complete packet

```text
(7,7,4,2,2,2,1).                                      (2)
```

The closest proved mechanism is THM-4157's exact row-A source, packet, and
labelled cubic carrier. Its canonical hostile is precisely `D_A=0`, where
the normalized leading endpoint and the independent source constant endpoint
both vanish. The corrected near miss is MISTAKE-509: one may not allocate an
arbitrary geometric subset of one irreducible carrier. The least-used
sidecar is the orbit count of the handle subgroup before adjoining all three
carrier meridians. That sidecar makes every finite bound below strict.

## 2. Complete wall chart and exhaustive endpoint algebra

Put `k=K`. On `(1)`, `eta!=0` and `D_A=0` force `k!=0`. There is a unique
parameter

```text
u=3eta/(2k) != 0.                                      (3)
```

The wall and the inherited `K` relation are then equivalent to

```text
K=k,                  Delta=(5696-90k)/105,
Theta=3u^2,           eta=2ku/3,          zeta=-2ku/3,
C=(-90k+315u^2+5696)/105.                              (4)
```

Conversely `(4)` satisfies `D_A=0`. Thus `(k,u,Phi)` with

```text
k*u*(5696-90k)*(-90k+315u^2+5696) != 0                (5)
```

is a bijective chart of `(1)`; there is no square-root sign quotient.

Define

```text
J=3Phi k+8ku+27u^3,

S=-180k^3-135k^2u^2+2752k^2+2160ku^2+4860u^4,

P=16380k^3+12285k^2u^2-115072k^2
  -710640ku^2-9534464k+1999872u^2.                    (6)
```

The first three rows of the table are manifestly disjoint. On `J=0`,

```text
Phi=-(8ku+27u^3)/(3k).                                 (7)
```

To close the last row without a sampled root, put `w=u^2` and write

```text
P=a_P(k)w+b_P(k),
a_P=63(195k^2-11280k+31744),
b_P=4k(4095k^2-28768k-2383616).                       (8)
```

Exact Euclidean division gives `gcd(a_P,b_P)=1`, so `P=0` has no lost
vertical component. Hence `P=0` implies `a_P!=0` and

```text
w=-b_P/a_P.                                            (9)
```

With

```text
H_4(k)=257063625k^4-13336482600k^3+141996773760k^2
       -1407227197440k+47646220845056,                 (10)
```

direct substitution gives the identity

```text
S(k,-b_P/a_P)
 =128k^2 H_4(k)/[49(195k^2-11280k+31744)^2].           (11)
```

The certificate proves that `H_4` is irreducible and squarefree over `Q`,
and that it is coprime to each actual required unit:

```text
k, 5696-90k, a_P, -b_P,
49725k^3-1350640k^2+9717760k-90406912,
the numerator and denominator of the terminal row in (18).              (12)
```

The cubic polynomial in `(12)` is, up to a nonzero scalar, the numerator of
`C` after `(9)`; the denominator is a nonzero multiple of `a_P`. Therefore
`J=S=P=0` is exactly the reduced quartic locus `H_4=0` together with `(7),(9)`,
all hypotheses `(5)` remain units there, and no component has been omitted.
Every coefficient point belongs to exactly one of the four displayed rows.

## 3. Two source projections and the four strict transforms

Use

```text
s=XT,              p=T+s^2,              t=p-s^2,

H=-3p+(8/3)p^2-(1376/135)p^3+k s^2p^2+Phi sp^3
  +Delta p^4+Theta s^2p^3+eta sp^4-eta s^3p^3,

G=-s^2/(2t)+H.                                           (13)
```

Define the polynomial critical rows

```text
A=(-sp+t^2 H_s)/p,
C_0=s^2+2t^2 H_p,
B=(C_0+sA)/t^2.                                         (14)
```

There is no hidden denominator: `A,B,C_0` are polynomials and

```text
t^2G_s=pA,                    2t^2G_p=C_0=t^2B-sA.      (15)
```

Thus either pair `(A,B)` or `(A,C_0)` cuts out the affine critical ideal
where `p*t!=0`. Their `s`-degrees and leading rows are

```text
deg_s(A,B)=(6,3),       LC_s(A,B)=(-2kup^2,-6kup^2)
                                      =(3zeta p^2,9zeta p^2),

deg_s(A,C_0)=(6,7),     LC_s(A,C_0)=(-2kup^2,-4kup^2)
                                      =(3zeta p^2,6zeta p^2).            (16)
```

Since `k*u!=0`, there is no source point at `s=infinity` for `p!=0`.

The primary Sylvester determinant has the exhaustive strict-transform
cascade

```text
J!=0:                 Res_s(A,B)=p^7 R_18(p),
J=0,S!=0:             Res_s(A,B)=p^8 R_17(p),
J=S=0,P!=0:           Res_s(A,B)=p^9 R_16(p),
J=S=P=0:              Res_s(A,B)=p^10 R_15(p).         (17)
```

The independent `(A,C_0)` determinant reconstructs the same residual
degrees with exceptional powers `p^9,p^10,p^11,p^12`. It computes its own
Sylvester determinant and reads the actual first and last rows; it imports no
primary resultant.

Here are the endpoint identities that prove `(17)`. Before `J=0`,

```text
[p^18]R_18
 =524288 k^5u^5(-90k+315u^2+5696)^4/364651875,

R_18(0)=82944k^2u^2 J.                                 (18a)
```

After `(7)`, the leading row is unchanged and

```text
R_17(0)=-(6912/5)ku^3 S.                               (18b)
```

In `Q(k)[u]/(S)`, the actual next row satisfies

```text
[p^1]R_17=-(96/35)k^2u^3 P.                           (18c)
```

Finally the actual following row is the exact identity

```text
[p^2]R_17=-(64/105)(u^3/k) Q(k,u^2),                  (18d)

Q(k,w)=
 -102060k^7-76545k^6w+316224k^6+4286520k^5w
 +12685824k^5+1683990k^4w^2-29961792k^4w
 +909058048k^4-17622360k^3w^2+417392640k^3w
 -51438240k^2w^3+1002637440k^2w^2
 +183708000kw^3+62001450w^4.                          (18e)
```

After substituting `(9)`, both numerator and denominator of `Q` are coprime
to `H_4`; the denominator is `343(195k^2-11280k+31744)^4` up to a nonzero
constant. Therefore `(18d)` never vanishes on `J=S=P=0`. Equations
`(18a)--(18e)`, together with the live top row, prove that the four residual
degrees are exactly `18,17,16,15`. There is no fifth stratum.

## 4. Exact affine critical lengths

The source chart collapses `p=0` and `t=0`. The canonical certificate audits
this loss directly in normalized coordinates. Put

```text
P_N=T+X^2T^2,                  Y_N=XTP_N,

G_N=-X^2T/2-3P_N+(8/3)P_N^2-(1376/135)P_N^3+kY_N^2
    +Phi P_N^2Y_N+Delta P_N^4+Theta P_NY_N^2
    +eta P_N^3Y_N-eta Y_N^3,

f=G_{N,X}/T,                  h=G_{N,T}.                (19)
```

One has

```text
deg_X(f,h)=(7,8),             LC_X(f)=LC_X(h)=8CT^7.   (20)
```

Thus `C!=0` excludes coordinate infinity away from the explicitly audited
collapsed rows. At `T=0`,

```text
f=-X,                         h=-(X^2+6)/2.             (21)
```

At `p=P_N=0` with `T!=0`, so `T=-1/X^2`,

```text
f=-(X^2-6)/X,                h=-(X^2-6)/2.             (22)
```

Consequently the source elimination omits exactly four universal Morse
critical points:

```text
T=0,      X^2=-6,      G_N=0,       det Hess(G_N)=+6,
T=-1/6,   X^2=6,       G_N=1/2,     det Hess(G_N)=-6.  (23)
```

Under a hypothetical Keller realization, the inherited Hessian congruence
from THM-4130/4147 makes every affine critical point Morse. On the open
`p*t!=0` chart, `(15)` is an invertible change of critical generators.
Therefore projected root collisions cannot reduce scheme length, and no
residual discriminant hypothesis is needed. Restoring `(23)` gives

| source residual degree | restored length |
|---:|---:|
| `18` | `L=22` |
| `17` | `L=21` |
| `16` | `L=20` |
| `15` | `L=19` |

This also explains why the saturated three-variable normalized resultant is
not a dependency of the proof: normalized coordinates are used only for the
direct no-infinity and universal-point gates `(20)--(23)`.

## 5. Boundary packet and labelled cubic carrier

Nothing in `D_A=0` contracts the row-A boundary calculation: `(1)` retains
`eta`, `Delta`, and `C`. THM-4157 therefore supplies exact genus `10`, packet
`(2)`, and defect

```text
(7-1)+(7-1)+(4-1)+3(2-1)=18=2*10-2.                   (24)
```

The nonrational face is

```text
q-1/2=kW^2-eta W^3.                                   (25)
```

Because `eta!=0`, this degree-three rational function gives one prime
separable cubic closed point over `C(q)`. Its exact discriminant is

```text
(q-1/2)(4k^3-27eta^2(q-1/2)).                         (26)
```

MISTAKE-509 is firewalled: all three conjugates respond together. If the
cubic maps to the origin, the full response has

```text
n=7+7+4+2+2+2+1=25,             ind(mu_O)=18.         (27)
```

If it maps to a finite horizontal carrier, finite-separable transport splits
it into three distinct conjugate sections, each with a transposition
meridian. Removing the cubic indices from the origin packet gives

```text
finite origin packet=(7,7,4,1),
n=19,                 ind(mu_O)=6+6+3=15,
m=3 carrier transpositions,            total carrier index=3.           (28)
```

There is no quotient loss in the `+3`: after the splitting base change, all
three meridians act on the same `19`-sheet cover, and every one has Cayley
index one. Products or cancellations can only decrease the upper bound used
below.

## 6. Carrier-orbit lemma and both contradictions

The needed refinement is self-contained.

> **Carrier-orbit lemma.** Let `n>m+1`. Suppose
> `X,Y,tau_1,...,tau_m` generate a transitive permutation action on `n`
> letters, every `tau_i` is a transposition, and transported fixed sheets
> give
>
> ```text
> a+b<=2n-L,       a=|supp X|,       b=|supp Y|.       (29)
> ```
>
> If `U=supp X union supp Y` and
> `kappa=|supp X intersect supp Y|`, then
>
> ```text
> |U|>=n-m,                  kappa<=n+m-L.             (30)
> ```

**Proof.** Put `H=<X,Y>`. Adjoining one transposition merges at most two
current orbits, so it lowers the orbit count by at most one. Transitivity
after adjoining all `m` transpositions therefore gives

```text
#Orb(H)<=m+1,                 n-#Orb(H)>=n-m-1.        (31)
```

The hypothesis `n>m+1` forces `U` to be nonempty: otherwise `H` has `n`
singleton orbits and `m` transpositions cannot make the action transitive.
Every point outside `U` is an `H`-fixed singleton, while nonempty `U`
contains at least one further `H`-orbit. Hence

```text
#Orb(H)>=n-|U|+1,             n-#Orb(H)<=|U|-1.        (32)
```

Combining `(31),(32)` yields `|U|>=n-m`. Finally

```text
kappa=a+b-|U|<=(2n-L)-(n-m)=n+m-L.                   (33)
```

Every inequality is in the direction used below. QED.

For a finite response, the punctured elliptic-target relation is

```text
[X,Y] mu_O tau_1 tau_2 tau_3=1.                       (34)
```

The commutator-overlap lemma and the triangle inequality for permutation
index give

```text
ind([X,Y])<=2kappa,
ind(mu_O)<=2kappa+m<=2(n+m-L)+m.                      (35)
```

For a full response, `X,Y` themselves generate transitively, so `|U|=n`,
`kappa<=n-L`, and

```text
ind(mu_O)=ind([X,Y])<=2(n-L).                         (36)
```

The exact ledgers are:

| `L` | full ceiling `2(25-L)` | finite `kappa<=22-L` | finite ceiling `2kappa+3` |
|---:|---:|---:|---:|
| `22` | `6<18` | `0` | `3<15` |
| `21` | `8<18` | `1` | `5<15` |
| `20` | `10<18` | `2` | `7<15` |
| `19` | `12<18` | `3` | `9<15` |

The finite type gate is explicitly satisfied in every row:

```text
n=19>m+1=4.                                            (37)
```

Thus neither response exists in any of the four exhaustive strata. This
proves the theorem.

## 7. Inherited limitation and exact remaining scope

The weakest inherited link is not the new elimination. It is the
finite-separable carrier and fixed-sheet transport of THM-4147, reused in
THM-4155 and THM-4157: after resolution and shrinking, the relative map must
be finite etale off the origin and complete carrier, the three split carrier
sections must remain distinct on the common sheet action, and the two node
fibres must transport all `L` Morse inverse points as fixed sheets. THM-4171
proves the source chart, endpoint coverage, critical lengths, and
carrier-orbit inequalities exactly, but it remains explicitly relative to
that transport package.

This theorem closes the whole `D_A=0` row-A wall under `(1)`, including all
of its endpoint intersections; it imposes no source-resultant discriminant.
It does not close the other row-A coefficient/resultant/discriminant walls
off `D_A`, the surviving walls in rows `B,C`, `eta=0`, `Delta=0`, entry into
the reduced seam, another reduced cell, exact residual weight at least ten,
`JC(2)`, or `DC(2)`.

## 8. Exact artifacts and replay

The primary certificate is the sole normalized-coordinate audit. The
independent certificate is deliberately scoped as a source-pair/endpoint
audit and does not overstate normalized independence.

```bash
python3 -B 04-computation/jc23_row_a_inner_resultant_wall_exclusion_thm4171.py
python3 -B -O 04-computation/jc23_row_a_inner_resultant_wall_exclusion_thm4171.py
PYTHONHASHSEED=271 python3 -B 04-computation/jc23_row_a_inner_resultant_wall_exclusion_thm4171.py

python3 -B 04-computation/jc23_row_a_inner_resultant_wall_exclusion_thm4171_independent_endpoint_audit.py
python3 -B -O 04-computation/jc23_row_a_inner_resultant_wall_exclusion_thm4171_independent_endpoint_audit.py
PYTHONHASHSEED=277 python3 -B 04-computation/jc23_row_a_inner_resultant_wall_exclusion_thm4171_independent_endpoint_audit.py
```

Normal, optimized, and fixed-hash-seed executions byte-match their frozen
outputs. Raw LF-byte and semantic hashes are bound in the front matter.

**QED.**
