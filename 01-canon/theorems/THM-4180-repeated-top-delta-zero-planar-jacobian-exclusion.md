---
id: THM-4180
title: "Repeated-top Delta-zero planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
  4147/4155/4157/4171/4173 + VERIFIED-EXACT + INDEPENDENTLY
  NORMALIZED-AUDITED. On the live exact-weight-nine reduced (2,3) seam, the
  whole Delta-zero repeated-top wall zeta=-eta, eta!=0 contains no
  nonautomorphic planar Keller pair. A diagonal W=sp face replaces the lost
  Delta-face and preserves the rational index-four place. Exact source and
  normalized resultants close all three surviving structural rows and the
  complete row-A inner wall. Together with THM-4176, this excludes the whole
  repeated-top locus eta!=0. Eta=0, entry, other cells, M>=10, JC(2), and
  DC(2) remain OPEN.
source: jc-delta-zero-dual-20260826
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
  - THM-4171-row-a-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4173-repeated-top-row-a-complete-planar-jacobian-exclusion
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
related:
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
  - THM-4164-y-only-triple-top-root-planar-jacobian-exclusion
  - THM-4165-y-only-inner-top-triple-intersection-planar-jacobian-exclusion
script: 04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180.py
output: 05-knowledge/results/jc23_repeated_top_delta_zero_complete_exclusion_thm4180.out
independent_audit_script: 04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_repeated_top_delta_zero_complete_exclusion_thm4180_independent_audit.out
script_sha256: d373ef7bb19df199d10818782ce0e8c8e9d396919594b6dd35e423907e4c4094
output_sha256: 33d2c9be7e49db4193a7f085f717b33e163e896462ab0fc1d0597c24c7dd8905
independent_audit_script_sha256: 51eff7406f396bacaf81f69396b8a3c2dac006dba31d45b0bde42e560661e0ff
independent_audit_output_sha256: 5f442520f5a4f9d6cab138900d2a5b5bf069d7da1580054a2b95969d45a26ba3
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone SymPy certificate reconstructs the contracted polygon,
  all five faces, the diagonal replacement, three strict transforms,
  source-chart loss ledger, complete source resultants, row-A inner-wall
  tower, carrier, packets, and response inequalities. Normal, optimized,
  and hash-seeded executions byte-match the frozen output.
independent_audit: >
  ACCEPT. A separate normalized (X,T) projection reconstructs the polygon
  half-spaces and faces, proves the Morse-resultant identity, reads the
  row-A inner wall from the opposite end of its divisor, recomputes B/C
  resultants, and obtains the same lengths. Normal, optimized, and
  hash-seeded executions byte-match the frozen output.
---

# THM-4180 -- repeated-top Delta-zero planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
4147/4155/4157/4171/4173 + VERIFIED-EXACT + INDEPENDENTLY
NORMALIZED-AUDITED; JC(2) AND DC(2) REMAIN OPEN.**

## 1. Theorem, inheritance, and exhaustive rows

Work over `C` in the live `b=d=0` reduced `(2,3)` seam at exact residual
weight nine. Put

```text
P=T+X^2T^2,                         Y=XTP,
K_0=2848/45,                        epsilon=-1376/135,
```

and use the complete normalized source

```text
G=-X^2T/2-3P+(8/3)P^2+epsilon P^3+K_0Y^2+Phi P^2Y
  +Theta P Y^2+eta P^3Y-eta Y^3.                       (1)
```

> **Theorem.** For every `(Theta,Phi,eta) in C^3` with `eta!=0`, `(1)` is
> not the normalized exact-weight-nine source of a nonautomorphic planar
> Keller pair in the inherited reduced seam. Equivalently, the complete wall
>
> ```text
> Delta=0,                    zeta=-eta,             eta!=0           (2)
> ```
>
> is empty of nonautomorphic planar Keller realizations.

The closest packet/carrier mechanism is THM-4157; the closest Delta-zero
face calculation is THM-4155; and the local-multiplicity and carrier-orbit
repairs are THM-4173 and THM-4171. The hostile degeneration is the lost
`Delta p^4` face: reusing the old polygon would be invalid. The least-used
sidecar is the diagonal coordinate `W=sp`, which restores its puncture.

The THM-4176 structural coefficients specialize exactly to

```text
Ctop=Delta+Theta=Theta,
Btop=epsilon+K_0=7168/135!=0.                           (3)
```

Thus its nomenclature gives the exhaustive disjoint partition

```text
A: Theta!=0;
B: Theta=0, Phi!=0;
C: Theta=Phi=0, with Btop=7168/135!=0.                 (4)
```

Row `D` requires `Btop=0`, so it is empty on `Delta=0`. There is no missing
fourth structural row.

## 2. Contracted polygon and replacement face

Use

```text
s=XT,                  p=T+s^2,                  t=p-s^2=T,
```

and write `(1)` as `G=-s^2/(2t)+H`, with

```text
H=-3p+(8/3)p^2+epsilon p^3+K_0s^2p^2+Phi sp^3
  +Theta s^2p^3+eta sp^3(p-s^2).                       (5)
```

For `Q=q^-1`, put `F_Q=(s^2-p)(1-QH)-Qs^2/2`. Exact support collection in
every row of `(4)` gives

```text
Newt(F_Q)=conv{(0,1),(2,0),(5,3),(1,5),(0,4)},
2Area=30,                 boundary=10,       g_arith=11. (6)
```

The polygon contracts, but its arithmetic genus remains eleven. Its edge
ledger and THM-4103 indices are

| edge | length | inward `(u,v;c)` | index `u+v-c` |
|---|---:|---:|---:|
| `(0,1)--(2,0)` | `1` | `(1,2;2)` | `1` |
| `(2,0)--(5,3)` | `3` | `(-1,1;-2)` | `2` |
| `(5,3)--(1,5)` | `2` | `(-1,-2;-11)` | `8` |
| `(1,5)--(0,4)` | `1` | `(1,-1;-4)` | `4` |
| `(0,4)--(0,1)` | `3` | `(1,0;0)` | `1` |

The five exact faces, in order, are

```text
s^2(1-Q/2)-p,
s^2((1-Q/2)-K_0Q(sp)^2+eta Q(sp)^3),
Q eta sp^3(p-s^2)^2,
Q p^4(epsilon+eta sp),
p(-1+Q(-3p+(8/3)p^2+epsilon p^3)).                    (7)
```

The last face is the same `s=0` affine-source edge audited in THM-4155; it
is not an infinity puncture. The first face gives a rational index-one
place. The second gives one separable cubic closed place with three geometric
index-two conjugates. The third is the repeated edge treated below.

The fourth is the new mechanism. In `W=sp`, its root is

```text
W=-epsilon/eta=1376/(135eta)!=0.                       (8)
```

It is one rational index-four place. The old horizontal face
`Delta+eta s=0` is gone, but its puncture and index survive on this diagonal
face. The cubic face gives

```text
q-1/2=K_0W^2-eta W^3,                                 (9)
```

a prime separable degree-three place because `eta!=0`.

## 3. Repeated-edge normalization and packets

The repeated edge has not moved. Set

```text
s=1/z,                    p=(1-a)/z^2,              r=1-a,
L(a,z)=z^11F_Q(1/z,r/z^2).                            (10)
```

Then

```text
L=az^9-Qz^9/2+Qeta a^2r^3-QTheta ar^3z-QPhi ar^3z^2
  -Qa(epsilon r^3+K_0r^2)z^3-Q(8/3)ar^2z^5+3Qarz^7.  (11)
```

The exhaustive local rows are

| row | leading branches | local `delta` | top indices |
|---|---|---:|---:|
| `A` | `a=(Theta/eta)z+...`, `a=-z^8/(2Theta)+...` | `1` | `(7,7)` |
| `B` | `a=(Phi/eta)z^2+...`, `a=-z^7/(2Phi)+...` | `2` | `(6,6)` |
| `C` | `a=(Btop/eta)z^3+...`, `a=-z^6/(2Btop)+...` | `3` | `(5,5)` |

The residue is, up to a unit, `omega=Qz^7 dz/L_a`. The two branches meet to
orders `1,2,3`, and `L_a` has those respective orders. With `(7)--(9)` this
gives

```text
A: (7,7,4,2,2,2,1),        genus <=10, defect=18;
B: (6,6,4,2,2,2,1),        genus <= 9, defect=16;
C: (5,5,4,2,2,2,1),        genus <= 8, defect=14.     (12)
```

Each defect is `2g-2` at the displayed upper genus. Finiteness and
connectedness below make Riemann--Hurwitz force equality, so these are
complete packets.

## 4. Alternate source projection and open-row lengths

Define

```text
A=(-sp+t^2H_s)/p,
C_0=s^2+2t^2H_p,
B=(C_0+sA)/t^2.                                       (13)
```

These are polynomials and

```text
t^2G_s=pA,                    2t^2G_p=t^2B-sA.         (14)
```

On `p*t!=0`, `(A,B)` is the critical ideal. Moreover

```text
deg_s(A,B)=(6,3),           LC_s(A,B)=(-3eta p^2,-9eta p^2). (15)
```

Thus no residual point with `p!=0` is lost at `s=infinity`. On `p=0`, one
has `A=-s,B|_{s=0}=-6`; on `t=0`, again `A=-s,B|_{s=0}=-6`. Hence the source
resultant contains no hidden or spurious point on either collapsed row.

At an `(A,B)` zero, the coordinate change `(X,T)->(s,p)` has Jacobian `t`.
Differentiating `(14)` at that zero gives

```text
p t^2 det D(A,B)=2t^4 det Hess_(s,p)(G).               (16)
```

For a Keller realization the right side is nonzero, so every intersection
is reduced. Resultant multiplicity counts distinct critical points even when
several share one `p`-coordinate.

The normalized chart restores exactly the four omitted points

```text
T=0,      X^2=-6,       G=0,       det Hess(G)=+6;
T=-1/6,   X^2= 6,       G=1/2,     det Hess(G)=-6.     (17)
```

Put `D_A=4Theta K_0^2-27eta^2`. In row `A`,

```text
Res_s(A,B)=p^6R_19(p),
[p^19]R_19=1327104 eta^5Theta^4,
R_19(0)=46656 eta D_A.                                (18)
```

Thus `D_A!=0` gives `L_A=19+2+2=23`. Rows `B,C` are direct specializations:

```text
B: Res_s(A,B)=p^6R_17,
   [p^17]R_17=777924 eta^5Phi^4,       R_17(0)=-1259712eta^3;

C: Res_s(A,B)=p^6R_15,
   [p^15]R_15=(168955354770571264/50625)eta^5,
   R_15(0)=-1259712eta^3.                              (19)
```

All endpoints are units in their rows, so

```text
L_B=17+2+2=21,                    L_C=15+2+2=19.       (20)
```

No critical-resultant discriminant is assumed.

## 5. Complete row-A inner-wall tower

On `D_A=0`, the unit hypotheses give the unique chart

```text
u=3eta/(2K_0)!=0,        Theta=3u^2,       eta=2K_0u/3. (21)
```

Define

```text
J=3K_0Phi+8K_0u+27u^3,
S_0=18225u^4-1515136u^2-129777664,
P_0=26081u^2+12924224.                                (22)
```

The exact source resultant is `p^7R_18`, with nonzero top

```text
[p^18]R_18=(524288/364651875)K_0^5u^5(315u^2)^4,      (23)
```

and bottom cascade

```text
R_18(0)=82944K_0^2u^2J,

J=0:
  R_18=pR_17,
  R_17(0)=-(6912/5)K_0u^3(4S_0/15),

J=S_0=0:
  [p^2]R_18 == -(96/35)K_0^2u^3(3584P_0/15) mod S_0. (24)
```

The terminal gate is

```text
gcd_Q[u](S_0,P_0)=1,
Res_u(S_0,P_0)=3466663036471111680!=0.                 (25)
```

Hence the exhaustive degree/length tower is

| stratum | residual degree | affine length |
|---|---:|---:|
| `J!=0` | `18` | `22` |
| `J=0,S_0!=0` | `17` | `21` |
| `J=S_0=0` | `16` | `20` |

There is no fourth stratum.

The independent audit uses `f=G_X/T,h=G_T`, proves

```text
T det D(f,h)=det Hess(G)+fG_XT,                        (26)
```

and reconstructs

```text
Res_X(f,h)=T^42(6T+1)^2Q_18(T),
Q_18(0)=-8957952u^12,
[T^18]Q_18=unit*u^6(8544Phi+22784u+1215u^3)^2.        (27)
```

On the last factor's zero locus, `[T^17]Q_18` is a nonzero scalar times
`u^6S_0^2`; modulo `S_0`, `[T^16]Q_18` is coprime to `S_0`. Its `X`-leading
rows are `24u^2T^7`. Thus the normalized top-degree cascade independently
gives `18,17,16` and audits every division and infinity unit in `(23)--(25)`.

## 6. Completeness, carrier, and contradictions

Equations `(18)--(25)` make the affine critical scheme finite. THM-3827 then
gives geometric connectedness. Riemann--Hurwitz turns the genus bounds in
`(12)` into equalities, so all three packets are complete.

The cubic place `(9)` has discriminant

```text
(q-1/2)(4K_0^3-27eta^2(q-1/2)).                       (28)
```

Prime residue degree and THM-4120 leave exactly the inherited two responses:

| row | finite origin packet | finite `(n,beta)` | full `n` |
|---|---|---:|---:|
| `A` | `(7,7,4,1)` | `(19,3)` | `25` |
| `B` | `(6,6,4,1)` | `(17,3)` | `23` |
| `C` | `(5,5,4,1)` | `(15,3)` | `21` |

Here `beta=3` is the total index of all three conjugate carrier
transpositions. They respond together. The proper-smooth/finite-etale order,
parallel Milnor cores, and fixed-sheet transport are the inherited
THM-4147/4157 contract.

Let `X,Y` be the handle permutations and `L=r_0+r_1`. Fixed-sheet transport
gives

```text
|supp X|+|supp Y|<=2n-L.                              (29)
```

For full responses, transitivity and the commutator-overlap lemma give
`ind(mu_O)<=2(n-L)`. For open rows `A,B,C`, the ceilings are

```text
4<18,                         4<16,                         4<14. (30)
```

For their finite responses, if a handle is nonidentity, total generator
index is at most

```text
2n-L-1+beta=17<18,        15<16,        13<14.         (31)
```

If both are identities, the index is at most `beta=3`; no action is
transitive.

On the row-A inner wall, `n=19`, `beta=3`, and `19>3+1`. The carrier-orbit
lemma of THM-4171 gives, for `L=22,21,20`,

```text
kappa<=19+3-L=0,1,2,
ind(mu_O)<=2kappa+3=3,5,7<15.                         (32)
```

The full ceilings are `2(25-L)=6,8,10<18`. Every row and coefficient
stratum is impossible. **QED.**

## 7. Combined corollary and firewalls

THM-4176 proves `eta*Delta!=0` on the same normalized repeated-top source.
The present theorem proves the disjoint slice `eta!=0,Delta=0`. Therefore:

> **Corollary.** Inside the live exact-weight-nine reduced `(2,3)` seam, the
> complete repeated-top locus
>
> ```text
> zeta=-eta,                         eta!=0                         (33)
> ```
>
> contains no nonautomorphic planar Keller pair.

This makes no claim on `eta=0` (hence `zeta=0` here), entry into the reduced
seam, other cells, residual weight at least ten, a transport package outside
the inherited proper-smooth setup, `JC(2)`, or `DC(2)`.

## 8. Exact artifacts and replay

```text
04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180.py
sha256 d373ef7bb19df199d10818782ce0e8c8e9d396919594b6dd35e423907e4c4094
05-knowledge/results/jc23_repeated_top_delta_zero_complete_exclusion_thm4180.out
sha256 33d2c9be7e49db4193a7f085f717b33e163e896462ab0fc1d0597c24c7dd8905

04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180_independent_audit.py
sha256 51eff7406f396bacaf81f69396b8a3c2dac006dba31d45b0bde42e560661e0ff
05-knowledge/results/jc23_repeated_top_delta_zero_complete_exclusion_thm4180_independent_audit.out
sha256 5f442520f5a4f9d6cab138900d2a5b5bf069d7da1580054a2b95969d45a26ba3
```

Replay each script against its output in normal, optimized, and fixed-hash
mode:

```bash
python3 -B 04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180.py
python3 -B -O 04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180.py
PYTHONHASHSEED=4180 python3 -B 04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180.py

python3 -B 04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180_independent_audit.py
python3 -B -O 04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180_independent_audit.py
PYTHONHASHSEED=4180 python3 -B 04-computation/jc23_repeated_top_delta_zero_complete_exclusion_thm4180_independent_audit.py
```

The primary and independent runs print `checks=96` and `checks=66` and
byte-match their frozen outputs in all three configurations.
