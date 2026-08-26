---
id: THM-4189
title: "Complete P-only Theta-zero planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
  4147/4176/4183/4186 + VERIFIED-EXACT + INDEPENDENTLY SOURCE-AUDITED. On
  the live exact-weight-nine b=d=0 reduced (2,3) seam, the complete P-only
  Theta-zero coefficient wall zeta=Theta=0, eta*Delta*K!=0 contains no
  nonautomorphic planar Keller pair, with Phi arbitrary and without any
  critical discriminant, Q_19(-1/6), or projected-coordinate separation
  hypothesis. Exactly one residual critical root escapes in both the
  normalized T and source p gauges; the remaining affine critical length is
  23. The Newton packet is (8,7,4,2,2,1), obtained by blowing the index-three
  and index-five faces down to one primitive index-seven face. Exact hostile
  fibres prove that both repeated projection and collision with the universal
  T=-1/6 fibre really occur on the wall. Together with THM-4183/4186 this
  closes the whole P-only K!=0 locus. The P-only K=0 wall, mixed B walls,
  other cells, entry, M>=10, JC(2), and DC(2) remain OPEN.
source: jc-theta-zero-wall-20260826
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
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
  - THM-4183-p-only-delta-zero-planar-jacobian-exclusion
  - THM-4186-complete-p-only-generic-coefficient-critical-wall-planar-jacobian-exclusion
related:
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
  - THM-4164-y-only-triple-top-root-planar-jacobian-exclusion
  - THM-4165-y-only-inner-top-triple-intersection-planar-jacobian-exclusion
script: 04-computation/jc23_p_only_theta_zero_complete_exclusion_thm4189.py
output: 05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_thm4189.out
independent_audit_script: 04-computation/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.py
independent_audit_output: 05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.out
script_sha256: 2355d423e805201e15526228b1392c7806b8c59dcab69829d4b4ff7985d2c690
output_sha256: 8a61eb9fb270069d826d05cdc23654f6442f8b39fcb12c781f3c269eece80c43
independent_audit_script_sha256: 32b0490929b6fa5fe906d764838dcb2327ba76e0f37a5a99a2ea4a5f8fc5197b
independent_audit_output_sha256: 37644e4fe198bacabb335914b281aab9a9613e9e939530f8df1f260804dc36e2
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone normalized (X,T) calculation reconstructs the complete
  source, row-independent Morse bridge, both universal fibres, and the full
  symbolic three-parameter T-resultant. It independently rebuilds every
  Newton face and packet index, isolates an algebraic repeated-projection
  hostile by Sturm count and quotient-ring gcds, and proves a second exact
  Q_19(-1/6) hostile is squarefree. Normal, optimized, and fixed-hash
  executions byte-match the frozen output.
independent_audit: >
  ACCEPT. A separate rational (s,p) implementation imports no primary
  polynomial. It derives a different critical pair, proves its Hessian bridge
  by ideal reduction, computes both transverse p^2 R_20 and wall p^2 R_19,
  verifies the reciprocal normal escape, reconstructs the polygon from
  supporting inequalities, and separates both normalized hostiles in the p
  coordinate. Normal, optimized, and fixed-hash executions byte-match.
---

# THM-4189 -- complete P-only Theta-zero planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
4147/4176/4183/4186 + VERIFIED-EXACT + INDEPENDENTLY SOURCE-AUDITED;
JC(2) AND DC(2) REMAIN OPEN.**

## 1. Theorem and inheritance pass

Work over `C` in the live `b=d=0` reduced `(2,3)` seam at exact residual
weight nine. Put

```text
P=T+X^2T^2,                    Y=XTP,
K=2848/45-(7/6)Delta.                                    (1)
```

Use the complete P-only source on the `Theta=0` wall,

```text
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K Y^2+Phi P^2Y+Delta P^4+eta P^3Y.                  (2)
```

> **Theorem.** For every coefficient tuple satisfying
>
> ```text
> zeta=Theta=0,                    eta*Delta*K!=0,       (3)
> ```
>
> polynomial `(2)` is not the normalized exact-weight-nine source of a
> nonautomorphic planar Keller pair in the inherited reduced seam.
>
> The coefficient `Phi` is unrestricted. No discriminant of a source or
> projected critical eliminant, no condition at `T=-1/6`, and no projected-
> coordinate separation hypothesis is imposed.

The closest proved mechanisms are THM-4186's full-source Morse/resultant
argument and THM-4183's `Theta=0` Newton blowdown. The canonical hostile is a
reduced source with two distinct Morse critical points sharing one normalized
`T`-coordinate. The corrected near miss is MISTAKE-421: reducedness of the
full source scheme does not make a chosen projection primitive or its
eliminant squarefree. MISTAKE-509 separately requires both geometric branches
of the quadratic carrier to respond together.

The least-used sidecar is the independent source projection

```text
s=XT,                    p=T+s^2,                     (4)
```

which separates the repeated-`T` hostile below. The selected method card is
“divide exceptional multiplicity before judging a wall”: retain `Theta`,
compute the first reciprocal normal coefficient in both projections, and
only then specialize to `Theta=0`.

## 2. Complete source and exact scope

In rational source coordinates put

```text
t=p-s^2=T,              y=sp,
G=-s^2/(2t)+H(p,y),                                      (5)

H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta p y^2+eta p^3y+zeta y^3.            (6)
```

THM-3992/3997 and THM-4007 give `(6)` as the complete normalized source.
For `wt(p)=2,wt(y)=3`, the allowed monomials through weight nine, after the
forced deletions `y,py`, are

```text
p,p^2,p^3,y^2,p^2y,p^4,py^2,p^3y,y^3.                 (7)
```

Only `p^3y,y^3` have weight nine. Consequently

```text
exact M=9  iff  (eta,zeta)!=(0,0).                     (8)
```

On `(3)`, `eta p^3y` retains exact weight nine and `(6)` specializes exactly
to `(2)`. There is no omitted source term or third top row.

## 3. Normalized critical divisor and affine length

Set

```text
f=G_X/T,                         h=G_T.                 (9)
```

Both are polynomials, and exact differentiation gives

```text
(deg_X f,deg_X h)=(8,9),
[X^8]f=[X^9]h=9eta T^8.                                (10)
```

Thus no finite `T!=0` common point is lost at `X=infinity`. The complete
symbolic resultant on `(3)` is

```text
Res_X(f,h)=T^56(6T+1)^2 Q_19(T),                       (11)

Q_19(0)=-(3^15/2^7)eta^7,
[T^19]Q_19=944784 K^6eta^7.                            (12)
```

Both endpoints are units. Hence the residual critical degree is exactly
nineteen for every `Phi` and every coefficient tuple in `(3)`.

The row-independent identity

```text
T det D(f,h)=det Hess(G)+f G_XT                        (13)
```

holds identically. At an affine critical point of a hypothetical Keller
realization `G=E(A_0,C_0)`, the unit source Jacobian and the two Morse target
nodes give

```text
det Hess(G)=det D(A_0,C_0)^2 det Hess(E)!=0.            (14)
```

Equations `(13)--(14)` make every `T!=0` source point reduced, even when
several points share one root of `Q_19`.

The factor `T^56` is a coordinate degree-drop artifact. Direct specialization
restores two reduced points:

```text
T=0,                 X^2=-6,
G=0,                 det Hess(G)=+6.                   (15)
```

The factor `(6T+1)^2` supplies the other universal pair:

```text
T=-1/6,              X^2=6,
G=1/2,               det Hess(G)=-6.                   (16)
```

If `Q_19(-1/6)=0`, its residual valuation counts additional distinct Morse
points in that fibre; it does not erase the universal pair. Therefore the
exact affine critical length is

```text
L=19+2+2=23.                                           (17)
```

Every critical point lies above one of the two target nodes. If `r_i` counts
the points above the node of value `i/2`, then

```text
r_0+r_1=23.                                             (18)
```

## 4. Independent source gauge and the one-root escape

Retain `Theta` transversely in `(6)`. In `(s,p)` coordinates define

```text
A=(-sp+t^2H_s)/p,
C_0=s^2+2t^2H_p,
B=(C_0+sA)/t^2.                                        (19)
```

These are polynomials and satisfy

```text
t^2G_s=pA,                    2t^2G_p=t^2B-sA.         (20)
```

Their `s`-degrees are `(5,2)`, with leading rows

```text
LC_s(A)=2p(K+Theta p),
LC_s(B)=8p(Theta p+3K/4),
LC_s(B)-4LC_s(A)+2Kp=0.                                (21)
```

On `Theta=0,pK!=0`, both leading rows are units, so this source projection
loses no open-chart point at `s=infinity`. Differentiating `(20)` and reducing
the exact numerator modulo `(A,B)` gives

```text
p det D(A,B)=2t^2 det Hess_(s,p)(G)       mod (A,B).   (22)
```

The transverse and wall resultants are

```text
Res_s(A,B)=p^2R_20(p),
[p^20]R_20=-65610eta^6Theta,

Theta=0:
Res_s(A,B)=p^2R_19(p),
R_19(0)=-31104K^2,
[p^19]R_19=-78732Keta^6.                               (23)
```

The exponent `p^2` belongs to this `(A,B)` presentation. THM-4186 used the
different pair `(A,C=2t^2G_p)` and obtained a `p^4` artifact. The differing
artifact exponents do **not** change the affine length; both gauges restore
the same four universal points `(15)--(16)`.

For `v=1/p`, the reciprocal initial form is

```text
v^20R_20(1/v)=-65610eta^6Theta-78732Keta^6v
              +terms in (Theta*v,v^2),

v=-5Theta/(6K)+higher terms.                           (24)
```

Thus exactly one residual root escapes at `Theta=0`. In the normalized
`u=1/T` gauge the inherited top coefficient is
`72900K^4eta^5Theta^4`, while `(12)` is the first live `u` coefficient, so

```text
u=-25Theta^4/(324K^2eta^2)+higher terms.               (25)
```

Equations `(24)--(25)` are projection-specific normal gauges. They certify
the same one-root loss but do not identify an escaping root across
projections.

## 5. Newton blowdown and complete packet

On `Theta=0`, write `(2)` in `(s,p)` coordinates as

```text
H=-3p+(8/3)p^2-(1376/135)p^3+Ks^2p^2+Phi sp^3
  +Delta p^4+eta sp^4.                                 (26)
```

For `Q=q^-1`, put

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2.                              (27)
```

Exact support collection gives

```text
Newton polygon:
(0,1),(2,0),(4,2),(3,4),(1,5),(0,5),

(2Area,B,g_Pick)=(28,10,10).                           (28)
```

The six exact boundary faces, in order, are

```text
s^2(1-Q/2)-p,
s^2((1-Q/2)-KQ(sp)^2),
-Qp^2s^3(Ks+eta p^2),
Qeta p^4s(p-s^2),
Qp^5(Delta+eta s),
p(-1+Q(-3p+(8/3)p^2-(1376/135)p^3+Delta p^4)).        (29)
```

The last face has `s=0,p!=0`, hence consists of affine source points and is
not in the infinity packet. The other primitive inward data and indices are

| inward data `(u,v;c)` | branch data | index |
|---|---|---:|
| `(1,2;2)` | rational | `1` |
| `(-1,1;-2)` | one quadratic closed point | `(2,2)` |
| `(-2,-1;-10)` | rational | `7` |
| `(-1,-2;-11)` | rational | `8` |
| `(0,-1;-5)` | rational | `4` |

The third face is the strict `Theta=0` transform: the former index-three and
index-five faces blow down to

```text
Ks+eta p^2=0.                                           (30)
```

Its edge vector is `(-1,2)` and `(30)` is linear in the primitive torus
coordinate, so it is one smooth index-seven branch. The new horizontal face
is

```text
Delta+eta s=0.                                          (31)
```

Because `Delta*eta!=0`, `(31)` has the unique nonzero torus root
`s=-Delta/eta`, and derivative `eta!=0`; it is simple. The coefficient `Phi`
occurs on no boundary face. Thus there is no hidden `Phi` or coefficient
subwall, and the candidate packet is

```text
packet=(8,7,4,2,2,1),
sum=24,                    defect=18.                  (32)
```

Equation `(11)` makes the affine critical scheme finite. THM-3827's closed-
polynomial factor theorem therefore gives geometric connectedness. Pick gives
normalization genus at most ten, while THM-4103 and Riemann--Hurwitz over the
elliptic target use defect eighteen to give genus at least ten. Hence the
genus is ten and `(32)` is complete: there is no hidden affine ramification
or additional boundary branch.

## 6. Quadratic carrier and both strict contradictions

The only nonrational place is

```text
q-1/2=K W^2.                                            (33)
```

Since `K!=0`, `(33)` is one separable irreducible quadratic closed point over
`C(q)`: `q-1/2` has odd valuation at `q=1/2`. THM-4120 gives
`E_q(C(q))={O}`. If `(33)` has finite horizontal image with residue field
`M`, then

```text
C(q) subset M subset C(q)(W).                          (34)
```

Prime degree and the absence of a finite rational target point force the
full quadratic field. The two geometric branches therefore respond together
as two conjugate carrier points, each with a transposition meridian. The only
exhaustive responses are

```text
full:       n=24;
finite:     n=20,                     beta=2.           (35)
```

THM-4147's proper-smooth/finite-etale construction and parallel Milnor-core
transport apply because Sections 3--5 supply a finite reduced affine critical
scheme under Keller, a connected source, the complete packet, and one
separable prime carrier. If `X_0,X_1` are the two handle permutations, then

```text
#Fix(X_i)>=r_i,
|supp X_0|+|supp X_1|<=2n-23.                          (36)
```

For the full response, the handles generate transitively and their supports
cover all 24 sheets. Their overlap is at most `24-23=1`, so THM-4147's
commutator-overlap lemma gives

```text
ind([X_0,X_1])<=2.                                     (37)
```

The origin meridian has index equal to the complete packet defect eighteen,
a contradiction.

For the finite response, if at least one handle is nonidentity, the handles
and two carrier transpositions have merger capacity at most

```text
2*20-23-1+2=18<19=20-1.                               (38)
```

If both handles are identities, the capacity is only `beta=2<19`.
Transitivity is impossible in either case. This proves the theorem. **QED.**

## 7. Exact hostile projection walls

Both projection conveniences deleted by the theorem are genuinely nonempty
inside `(3)`.

### 7.1 Two reduced points share one normalized coordinate

Let `tau` be the unique real root in

```text
315/1000 < tau < 316/1000                             (39)
```

of

```text
S(t)=16512t^6-174912t^5-235020t^4-35764t^3
     +10607t^2+4530t+1485.                            (40)
```

The exact Sturm count on `(39)` is one, with opposite endpoint signs. Put

```text
N(t)=33024t^6+265344t^5+263496t^4+33064t^3
     +3973t^2-2100t-1485,

Delta=(1376t^2-240t+135)/(180t^3),
Phi=N(t)/(1620t^4(t+1)^2(3t+1)),
eta=-N(t)/(1620t^5(t+1)^2(3t+1)),                     (41)
```

and specialize `t=tau`. Every denominator in `(41)` is nonzero on `(39)`:
`tau>0`, `tau+1>0`, and `3tau+1>0`. Exact polynomial gcd checks against
`S` show that the numerators of `Delta`, `eta`, and

```text
K=(68352t^3-9632t^2+1680t-945)/(1080t^3)              (42)
```

are coprime to `S`; hence `Delta*eta*K!=0` at `tau`.

In `Q(tau)[X]`, the normalized critical fibre at `T=tau` has

```text
gcd_X(f,h)=X(X-1).                                     (43)
```

The audit checks `(43)` by exact remainder reduction modulo `S`, and checks
that the quotient-polynomial resultant is coprime to `S`, so there is no
additional common root. Both Hessian numerators are also coprime to `S`;
therefore the points `X=0,1` are distinct and Morse. Their source coordinates
are

```text
(s,p)=(0,tau),             (tau,tau+tau^2).            (44)
```

The `p`-coordinates differ by `tau^2!=0`. Thus the full source scheme is
reduced and the source projection separates the points, while the normalized
`T`-eliminant has a repeated root. Exact gcd checks also show neither source
value is `0` or `1/2`, so this is a projection-wall hostile, not a Keller
candidate.

### 7.2 An extra point joins the universal fibre

Take

```text
Delta=1,                 Phi=-12416/25,
Theta=0,                 eta=-32856/125,
K=5591/90.                                                (45)
```

Every coefficient required nonzero in `(3)` is a unit. Exact elimination gives

```text
Q_19(-1/6)=0,
gcd_X(f(X,-1/6),h(X,-1/6))=(X-1)(X^2-6).              (46)
```

The residual `Q_19` is squarefree and its root at `-1/6` is simple. The extra
point is Morse, with

```text
det Hess(G)(1,-1/6)=-1375271635/68024448,
G(1,-1/6)=2050529/5038848.                             (47)
```

In the source chart it is

```text
(s,p)=(-1/6,-5/36),                                    (48)
```

separate from the universal `p=0` pair. The complete specialized source
`R_19` is squarefree. Again `(47)` is neither target-node value, so `(45)` is
a wall witness rather than a Keller candidate.

These hostiles show why discriminant and `Q_19(-1/6)` deletions are
projection conveniences, not geometric necessities.

## 8. Combined P-only closure and firewalls

THM-4183 closes `Delta=0` for both `Theta` rows. THM-4186 closes
`Delta*K*Theta!=0`. The present theorem closes the remaining `Theta=0` row
when `Delta*K!=0`. Since `K=2848/45-(7/6)Delta` is nonzero at `Delta=0`, the
three theorems combine to give:

> **P-only scope corollary.** In the live exact-weight-nine reduced `(2,3)`
> seam, the complete P-only locus
>
> ```text
> zeta=0,                    eta*K!=0                  (49)
> ```
>
> contains no nonautomorphic planar Keller pair, with `Delta,Phi,Theta`
> otherwise arbitrary.

The only remaining P-only coefficient wall is

```text
K=0,                       Delta=5696/105.              (50)
```

There the quadratic carrier equation `(33)` collapses, so neither the prime-
carrier response nor the packet ledger may be transported without a new
boundary normalization. This theorem also does not close the mixed `B`
chamber or its coefficient/critical walls, another reduced cell, entry into
the seam, exact residual weight at least ten, a setup outside the inherited
proper-smooth contract, `JC(2)`, or `DC(2)`.

## 9. Exact artifacts and replay

Primary normalized-coordinate certificate:

```text
04-computation/jc23_p_only_theta_zero_complete_exclusion_thm4189.py
sha256 2355d423e805201e15526228b1392c7806b8c59dcab69829d4b4ff7985d2c690

05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_thm4189.out
sha256 8a61eb9fb270069d826d05cdc23654f6442f8b39fcb12c781f3c269eece80c43
```

Independent source-chart certificate:

```text
04-computation/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.py
sha256 32b0490929b6fa5fe906d764838dcb2327ba76e0f37a5a99a2ea4a5f8fc5197b

05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.out
sha256 37644e4fe198bacabb335914b281aab9a9613e9e939530f8df1f260804dc36e2
```

Replay both paths in normal, optimized, and fixed-hash modes:

```bash
python3 -B \
  04-computation/jc23_p_only_theta_zero_complete_exclusion_thm4189.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_thm4189.out -
python3 -B -O \
  04-computation/jc23_p_only_theta_zero_complete_exclusion_thm4189.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_thm4189.out -
PYTHONHASHSEED=4189 python3 -B \
  04-computation/jc23_p_only_theta_zero_complete_exclusion_thm4189.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_thm4189.out -

python3 -B \
  04-computation/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.out -
python3 -B -O \
  04-computation/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.out -
PYTHONHASHSEED=4189 python3 -B \
  04-computation/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_theta_zero_complete_exclusion_independent_audit_thm4189.out -
```

All three primary streams print `checks=70` and
`THM4189_NORMALIZED_PRIMARY_EXACT_ACCEPT`. All three independent streams
print `checks=76` and `THM4189_SOURCE_INDEPENDENT_EXACT_ACCEPT`.
