---
id: THM-4192
title: "Complete P-only K-zero planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
  4147/4186/4189 + VERIFIED-EXACT + INDEPENDENTLY SOURCE-AUDITED. On the
  live exact-weight-nine b=d=0 reduced (2,3) seam, the complete P-only K-zero
  wall zeta=K=0, eta!=0, Delta=5696/105 contains no nonautomorphic planar
  Keller pair, with Phi and Theta arbitrary. The exhaustive Theta-nonzero,
  Theta-zero Phi-nonzero, and terminal Theta=Phi=0 strata have affine critical
  lengths 22,19,18 and complete rational packets (8,5,5,4,1), (8,4,4,2,1),
  (8,5,4,1). Their response gaps are 1,0,0 and their commutator caps 2,0,0,
  against origin defects 18,14,14. Together with THM-4189 this closes the
  whole P-only exact-weight-nine row zeta=0,eta!=0. Mixed B, other cells,
  entry, M>=10, JC(2), and DC(2) remain OPEN.
source: jc-k-zero-wall-20260826
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
  - THM-4186-complete-p-only-generic-coefficient-critical-wall-planar-jacobian-exclusion
  - THM-4189-complete-p-only-theta-zero-planar-jacobian-exclusion
related:
  - THM-4183-p-only-delta-zero-planar-jacobian-exclusion
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
script: 04-computation/jc23_p_only_k_zero_complete_exclusion_thm4192.py
output: 05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_thm4192.out
independent_audit_script: 04-computation/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.py
independent_audit_output: 05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.out
script_sha256: 45a8e28d859ed1f78743fe3c469de6cba89d3e3cfff6647ec1d90a950cf1aef9
output_sha256: 5921d2c62e42a843e72a0128731c36c5e8d52af96bd6c4fc3117319f89fa2ed1
independent_audit_script_sha256: 8cad1da6f86248c090a1813cab823462d44337d2e65a84b88e0a8b706748b89e
independent_audit_output_sha256: fe237f69cb5ef8c3b162bef2cc59592943b40d017d8027b8c8a34138ece5589c
hash_basis: raw LF bytes
primary_audit: >
  PASS. A normalized (X,T) implementation reconstructs the complete source,
  Hessian bridge, both universal fibres, all three exact residual resultants,
  every Newton hull and face, the rational packet indices, and all three
  strict commutator contradictions. It also verifies two exact projection
  hostiles: a repeated T-value with two Morse points and a residual point on
  the universal T=-1/6 fibre. These are probes, not hypotheses.
independent_audit: >
  ACCEPT. A separate rational (s,p) implementation uses the distinct (A,B)
  critical pair, proves the source Hessian bridge by ideal reduction, excludes
  p=0, t=0, and source-infinity losses, computes p^4R18, p^2R15, and p^3R14,
  enforces direct Sylvester recomputation after both coefficient drops,
  transports both projection hostiles to source coordinates, and independently
  reconstructs the hulls from supporting halfspaces. Normal, optimized, and
  fixed-hash streams byte-match.
---

# THM-4192 -- complete P-only K-zero planar Jacobian exclusion

## Verdict and exact scope

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
4147/4186/4189 + VERIFIED-EXACT + INDEPENDENTLY SOURCE-AUDITED; JC(2) AND
DC(2) REMAIN OPEN.** In the live exact-`M=9`, `b=d=0` reduced `(2,3)` seam,
the complete P-only coefficient wall

```text
zeta=0, eta!=0, K=0,
K=2848/45-(7/6)Delta,
Delta=5696/105
```

contains no nonautomorphic planar Keller pair.  `Phi,Theta` are arbitrary.
The exhaustive strata are

```text
I:   Theta!=0                    (Phi arbitrary),
II:  Theta=0, Phi!=0,
III: Theta=Phi=0.
```

> **Theorem.** Every coefficient tuple in one of strata I--III is excluded
> from being the normalized exact-weight-nine source of a nonautomorphic
> planar Keller pair in the inherited reduced seam.

Together with THM-4189 (complete P-only `K!=0`) this closes the whole
P-only exact-weight-nine row.  It does not close the mixed B chamber, entry,
other reduced cells, `M>=10`, `JC(2)`, or `DC(2)`.

## Inheritance and weakest link

- Source completeness: THM-3992/3997/4007.
- Generic-source connectedness from a finite critical scheme: THM-3827.
- Exact Keller residue identity and toric index rule: THM-4103.
- Target generic-fibre rational points `E_q(C(q))={O}`: THM-4120.
- Target/pencil normalization and polynomial response setup: THM-4122 and
  THM-4130.
- Fixed-sheet transport and the commutator-overlap bound: THM-4147.
- Closest coefficient-wall mechanisms: THM-4186 and THM-4189.

The weakest inherited link remains THM-4147's relative-normalization,
proper-smooth shrinking, and parallel Milnor-core fixed-sheet transport.
The present wall has no finite carrier, so it does not use the prime-carrier
or Mordell--Weil section classification beyond `E_q(C(q))={O}`.

The canonical hostile is exactly the one warned about by MISTAKE-509: the
generic quadratic carrier cannot be specialized as two labelled points.
At `K=0` it is replaced by a **primitive rational toric branch**, so the
correct response has no finite carrier at all.

## Complete source

With

```text
s=XT, p=T+s^2, t=p-s^2=T, y=sp,
```

the source is

```text
G=-s^2/(2t)+H,
H=-3p+(8/3)p^2-(1376/135)p^3
  +Phi*s*p^3+(5696/105)*p^4+Theta*s^2*p^3+eta*s*p^4.
```

Exact weight nine is retained by `eta p^3y`, because `eta!=0`; there is no
filtration exit and no omitted top term.

## Critical divisors and affine lengths

### Normalized path

Put

```text
f=G_X/T, h=G_T.
```

Both are polynomials of `X`-degrees `(8,9)` with common leading row
`9 eta T^8`, so no `T!=0` point is lost at `X=infinity`.  One symbolic
resultant specializes to

```text
Res_X(f,h)=T^56(6T+1)^2 Q_d(T),
```

with

| stratum | `d` | `Q_d(0)` | `[T^d]Q_d` |
|---|---:|---|---|
| I | 18 | `-3^15 eta^7/2^7` | `656100 Theta^6 eta^5` |
| II | 15 | `-3^15 eta^7/2^7` | `6561 Phi^6 eta^5/4` |
| III | 14 | `-3^15 eta^7/2^7` | `59049 eta^9/4` |

Every endpoint is a unit in its stated stratum.  The row-independent bridge

```text
T det D(f,h)=det Hess(G)+f G_XT
```

and the Keller Hessian congruence make every full source critical point
Morse; projected roots may still repeat, and no squarefreeness is assumed.

The coordinate artefact restores exactly

```text
T=0,    X^2=-6, G=0,   det Hess(G)=+6;
T=-1/6, X^2= 6, G=1/2, det Hess(G)=-6.
```

Thus the exact affine critical lengths are

```text
I: L=18+2+2=22;
II: L=15+2+2=19;
III: L=14+2+2=18.
```

### Independent rational-source path

Define

```text
A=(-sp+t^2 H_s)/p,
C0=s^2+2t^2 H_p,
B=(C0+sA)/t^2.
```

These are polynomial and satisfy

```text
t^2 G_s=pA,
2t^2 G_p=t^2B-sA,
p det D(A,B)=2t^2 det Hess_(s,p)(G) mod (A,B).
```

The three direct source resultants are

| stratum | `(deg_s A,deg_s B)` | leading rows | resultant | residual endpoints |
|---|---:|---|---|---|
| I | `(5,2)` | `2Theta p^2,8Theta p^2` | `p^4R_18` | `R_18(0)=-31104Theta^2`, `LC=-65610Theta eta^6` |
| II | `(4,1)` | `p^2(Phi+eta p),p^2(7Phi+9eta p)` | `p^2R_15` | `R_15(0)=1296Phi`, `LC=6561eta^5` |
| III | `(4,1)` | `eta p^3,9eta p^3` | `p^3R_14` | `R_14(0)=1296eta`, `LC=6561eta^5` |

There is a mandatory three-row degree-drop firewall.  If
`Res_s^(5,2)` denotes the row-I Sylvester determinant with its declared
`s`-degrees frozen, then

```text
Res_s^(5,2)(A,B)|_(Theta=0)=0.
```

It is therefore invalid to specialize that determinant into row II.  One
must first specialize the source pair, observe the degree drop `(5,2)->(4,1)`,
and recompute

```text
Res_s(A|_(Theta=0),B|_(Theta=0))=p^2R_15.
```

At `Phi=0` the source pair is specialized and its resultant is recomputed
again, giving `p^3R_14`; the additional `p` is a genuine endpoint loss and
cannot be read from the row-II unit formula.  The independent audit checks
both direct recomputations; it also verifies that the direct row-III
resultant equals the ordinary `Phi=0` specialization of the already-correct
row-II polynomial resultant.

In II the two leading rows have resultant `-2Phi eta`, so they never vanish
together.  All remaining leading rows are visibly units for `p!=0`.  On
`p=0` one has `A=-s,B(0,0)=-6`; on `t=0`, `A=-s,B(0,0)=-6`.  Hence every
residual point lies in `pt!=0`, no source point is lost at `s=infinity`, and
the source calculation independently gives the same `18,15,14` residual
lengths plus the four normalized points.

The different artefact exponents in the two projections are intentional:
they confirm that only the saturated residual degree is geometric.

## Exact projection hostiles (audits, not hypotheses)

The first hostile forbids any squarefree-`T` shortcut in stratum I.  Let
`tau` be the unique real root in `(51/200,32/125)` of

```text
C(tau)=68352tau^3-9632tau^2+1680tau-945,
```

put

```text
B(tau)=68352tau^6+205056tau^5+195424tau^4+49088tau^3
       -7952tau^2+1680tau-315,
eta=-2B(tau)/(315tau^4(tau+1)^2(5tau+2)),
Phi=-eta tau,
Theta=(136704tau^7+410112tau^6+390848tau^5+98176tau^4
       -15904tau^3+3360tau^2+945tau+630)
      /(630tau^4(tau+1)^2(5tau+2)).
```

Exact gcd tests modulo `C` show that all denominators and the required
`Theta,eta` units are nonzero.  Over `Q(tau)`, at `T=tau`,

```text
gcd_X(f,h)=X(X-1).
```

There is no further common root after dividing both equations by `X(X-1)`.
Both points are Morse, and neither has `G`-value `0` or `1/2`.  Their source
coordinates are

```text
(s,p)=(0,tau), (tau,tau+tau^2),
```

so the full source projection separates them although the normalized
`T`-eliminant has a repeated root.

The second hostile forbids deleting the universal fibre from the residual
factor in stratum II.  At

```text
Theta=0, Phi=201808/1575, eta=2766784/875,
```

one has

```text
T=-1/6: gcd_X(f,h)=(X-1)(X^2-6),
Q_15(-1/6)=0, Q_15'(-1/6)!=0,
gcd(Q_15,Q_15')=1.
```

Thus a third, ordinary point `X=1` joins the two universal points while the
residual remains squarefree.  At this point

```text
G=42257/91854,
det Hess_(X,T)(G)=-27296015083/937461924,
(s,p)=(-1/6,-5/36),
det Hess_(s,p)(G)=-27296015083/26040609.
```

The source audit also checks that `p=-5/36` is a simple root of the direct
row-II residual `R_15`.  Neither hostile is an assumption in the exclusion;
they certify that the proof counts the full reduced critical scheme rather
than roots of either projection.

## Newton faces and exact packets

For `Q=q^-1`, put

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2.
```

For a primitive inward normal `(u,v;c)`, THM-4103 gives

```text
e=u+v-c.
```

The vertical face `(1,0;0)` has `s=0,p!=0`, hence `X=0,T=p`; it is affine
and is not part of the infinity packet.

### I: `Theta!=0`

```text
hull=(0,1),(2,0),(4,3),(3,4),(1,5),(0,5),
(2Area,B,g_Pick)=(27,9,10).
```

The five nonaffine faces are

```text
(1,2;2):    s^2(1-Q/2)-p,
(-3,2;-6):  s^2[(1-Q/2)-Q Theta s^2p^3],
(-1,-1;-7): -Qp^3s^3(Theta s+eta p),
(-1,-2;-11):Qeta p^4s(p-s^2),
(0,-1;-5):  Qp^5(5696/105+eta s).
```

Their indices are `1,5,5,8,4`, hence

```text
packet=(8,5,5,4,1), sum=23, defect=18.
```

The old quadratic edge and adjacent index-three edge have become the single
primitive edge `(2,0)->(4,3)`.  Its torus coordinate `W=s^2p^3` satisfies

```text
q-1/2=Theta W,
```

so it is one rational index-five place, not a limiting quadratic pair.

### II: `Theta=0,Phi!=0`

```text
hull=(0,1),(2,0),(3,3),(3,4),(1,5),(0,5),
(2Area,B,g_Pick)=(23,9,8).
```

The changed faces are

```text
(-3,1;-6): s^2[(1-Q/2)-Q Phi sp^3],
(-1,0;-3): -Qs^3p^3(Phi+eta p),
```

with the low, top and horizontal faces unchanged.  The indices are
`1,4,2,8,4`, giving

```text
packet=(8,4,4,2,1), sum=19, defect=14.
```

### III: `Theta=Phi=0`

```text
hull=(0,1),(2,0),(3,4),(1,5),(0,5),
(2Area,B,g_Pick)=(22,8,8).
```

The merged face is

```text
(-4,1;-8): s^2[(1-Q/2)-Q eta sp^4],
```

and the indices are `1,5,8,4`, giving

```text
packet=(8,5,4,1), sum=18, defect=14.
```

Every nonaffine face factor is linear in its primitive edge monomial.  Thus
all displayed places are separable and rational over `C(q)`.  The displayed
defects equal `2g_Pick-2`.  Finiteness of the critical scheme plus THM-3827
gives geometric connectedness; Pick gives the genus upper bound, while
Riemann--Hurwitz over the elliptic target gives the matching lower bound.
Therefore the genera are `10,8,8` and all three packets are complete.

## Response and monodromy contradiction

THM-4120 gives `E_q(C(q))={O}`.  Every boundary place is rational, so every
one maps to `O`; there is no finite carrier response.  Polynomiality keeps
affine source points away from target infinity.  Hence

```text
(n,L,n-L,origin defect)=
I:   (23,22,1,18),
II:  (19,19,0,14),
III: (18,18,0,14).
```

Let `X_0,X_1` be the two transported handle permutations.  The affine Morse
points give `#Fix(X_i)>=r_i`, with `r_0+r_1=L`.  Because this is the full
response, `X_0,X_1` generate transitively and their supports cover all `n`
sheets.  Thus

```text
|supp(X_0) intersection supp(X_1)|<=n-L,
ind([X_0,X_1])<=2(n-L)=2,0,0.
```

The origin meridian is the inverse commutator and has index respectively
`18,14,14`, a strict contradiction in all three strata.

## Exact artifacts, controls, hashes, and replay

Primary normalized certificate:

```text
04-computation/jc23_p_only_k_zero_complete_exclusion_thm4192.py
sha256 45a8e28d859ed1f78743fe3c469de6cba89d3e3cfff6647ec1d90a950cf1aef9

05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_thm4192.out
sha256 5921d2c62e42a843e72a0128731c36c5e8d52af96bd6c4fc3117319f89fa2ed1
```

Independent rational-source certificate:

```text
04-computation/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.py
sha256 8cad1da6f86248c090a1813cab823462d44337d2e65a84b88e0a8b706748b89e

05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.out
sha256 fe237f69cb5ef8c3b162bef2cc59592943b40d017d8027b8c8a34138ece5589c
```

Both paths contain disjoint exact squarefree controls in every stratum.  The
primary controls also avoid `T=-1/6`; the two displayed hostiles separately
show that repeated projection roots and residual/universal-fibre collisions
do occur.  All are audit probes only: squarefreeness and fibre separation are
**not hypotheses**.

Normal, optimized, and fixed-hash executions byte-match the frozen outputs:

```bash
python3 -B 04-computation/jc23_p_only_k_zero_complete_exclusion_thm4192.py \
  | diff -u 05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_thm4192.out -
python3 -B -O 04-computation/jc23_p_only_k_zero_complete_exclusion_thm4192.py \
  | diff -u 05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_thm4192.out -
PYTHONHASHSEED=4192 python3 -B \
  04-computation/jc23_p_only_k_zero_complete_exclusion_thm4192.py \
  | diff -u 05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_thm4192.out -

python3 -B \
  04-computation/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.py \
  | diff -u \
  05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.out -
python3 -B -O \
  04-computation/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.py \
  | diff -u \
  05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.out -
PYTHONHASHSEED=271828 python3 -B \
  04-computation/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.py \
  | diff -u \
  05-knowledge/results/jc23_p_only_k_zero_complete_exclusion_independent_audit_thm4192.out -
```

No unresolved coefficient or projected-critical hostile remains on this
wall.  The only honest residual is the inherited fixed-sheet transport
package named above; globally, the other exact cells and entry remain open.
