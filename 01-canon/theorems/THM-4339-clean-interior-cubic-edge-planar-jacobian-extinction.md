---
id: THM-4339
title: "Clean interior cubic-edge planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4220/4230/4327/4337 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. In the inherited reduced (2,3), exact-weight-twelve
  seam, the residual gate Z=beta_11=zeta_3=0 and U*K*W*(U+W)!=0 is extinct.
  The cubic edge is treated through its squarefree, double, and triple-root
  strata. A double root gives an A11 rational bridge or a horizontal node;
  a triple root gives a j=0 elliptic tail of good-form order 26 exactly when
  its first normal coefficient is nonzero, and otherwise a horizontal
  delta-one conductor with rational normalization. Root exits K=0 or W=0,
  the U+W=0 contact, exact-seam entry, JC(2), and DC(2) remain open.
source: root + labelled-cubic agents / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4220-weight-ten-zeta-zero-genus-two-planar-jacobian-exclusion
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4337-zeta-zero-exact-weight-twelve-endpoint-wall-extinction
related:
  - THM-4292-lambda-zero-repeated-face-keller-extinction
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
  - THM-4334-beta-zero-exact-weight-twelve-endpoint-wall-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_clean_cubic_edge_extinction_thm4339.py
primary_output: 05-knowledge/results/jc2_m12_clean_cubic_edge_extinction_thm4339.out
primary_script_sha256: 2adb27b91b9ca712efb8b496806bcc5d92d17bd4d6fd8f971f4f97fb80e8f5b8
primary_output_sha256: b951b0596b5d7ad699b9ad554605bc231b1fba7cd34f1305f584e57fe909434f
independent_audit_script: 04-computation/jc2_m12_clean_cubic_edge_extinction_independent_audit_thm4339.py
independent_audit_output: 05-knowledge/results/jc2_m12_clean_cubic_edge_extinction_independent_audit_thm4339.out
independent_audit_script_sha256: a152a359b4584e88789af73192f794172c1b64595445f4ac14e8ecbee3af99c9
independent_audit_output_sha256: 132060e7cbeb54c3e78a1fd274049274e83eb5bded62aa94e0c884542046c121
hostile_referee_script: 04-computation/jc2_m12_clean_cubic_edge_extinction_hostile_referee_thm4339.py
hostile_referee_output: 05-knowledge/results/jc2_m12_clean_cubic_edge_extinction_hostile_referee_thm4339.out
hostile_referee_script_sha256: 7390db8738b5ae7c08c5c98f414a360b36e2659c9609a22565b2289d6d870553
hostile_referee_output_sha256: f653ab5c2dcc49510081ce90b27330744f96dd3f2fd85205da948975445e2f0e
hash_basis: raw LF bytes
audit: >
  PASS AFTER REPAIR. The primary path pins THM-4327's full-support hull
  engine, exhausts 8,192 optional-row/collision states, and checks the exact
  M,T complex, polygons, edge packet, cubic invariants, reciprocal chart,
  weighted carriers, differential orders, genus ledgers, and root-exit
  hostiles. The import-free audit reconstructs the exact source chart and
  independently checks the labelled-root formula, invariant Hasse gates,
  local normalizations, carrier types, edge schemes, and Pick ledgers. The
  hostile referee caught the false provisional uniform b1=11 claim: the
  smoothed double-root path has b1=12. It then verified the repaired stratum
  table, the general simple-zero remainder, and the finite proper DVR-flat
  horizontal normalization. Normal and optimized streams byte-match all
  three frozen outputs.
---

# THM-4339 -- Clean interior cubic-edge planar Jacobian extinction

**PROVED RELATIVE TO THM-4220/4230/4327/4337 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THE GATE `Z=beta_11=zeta_3=0` AND
`U*K*W*(U+W)!=0` IS EXTINCT. ROOT EXITS, THE `U+W=0` CONTACT, EXACT-SEAM
ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement and inheritance

Work over an algebraically closed field of characteristic zero in the
inherited reduced `(2,3)`, exact-weight-twelve seam. Impose

```text
Z=0,             beta_11=0,             zeta_3=0,                 (1)
U*K*W*Lambda!=0,                    Lambda=U+W.                    (2)
```

Every other allowed lower coefficient is arbitrary, with the inherited
relation `K=2848/45-(7/6)Delta` retained. No nonautomorphic planar Keller
pair lies on `(1)--(2)`. The claim is relative to the inherited exact seam
and proper-flat target interface; it proves neither seam entry nor `JC(2)`
nor `DC(2)`.

The inheritance pass was:

- closest mechanism: THM-4337 Section 8 isolates the cubic edge and gives
  its positive primary differential order;
- canonical hostile: `A(P)=(P-1)^2(P+1)`, a genuine repeated nonzero edge
  root which defeats a squarefree-face argument;
- corrected near miss: MISTAKE-540 forbids computing genus on a ramified
  root cover instead of the invariant exceptional function field;
- least-used sidecar: THM-4220 Section 7.2's finite simultaneous
  normalization of a persistent horizontal conductor.

The live board was

```text
labelled roots | Laurent exits | collision algebra | normal Hasse jet
weighted carrier | conductor normalization | good form | graph genus.    (3)
```

## 2. Exact lower complex, polygons, and edges

Use the sixteen-row support and the inherited lift

```text
(i,j) -> (j+2,i+j,1), (j,i+j+1,1),                      (4)
```

including the fixed residual point `(2,0,1)`. Require the fixed
`p,p^2,p^3` rows and the `U,W,K` rows, and delete `Z,zeta_3,beta_11`.
Seven rows remain optional. Independently delete the six inherited
multiply-owned aggregate points

```text
(2,3,1),(2,4,1),(2,5,1),(2,6,1),(3,4,1),(3,5,1).       (5)
```

Every one of the `2^7 2^6=8,192` conservative hostile states has exactly

```text
M=(1/12,1/6,-1/6),                 T=(1/2,0,-1).         (6)
```

The uniquely `K`-owned point `(4,2,1)` prevents `T` from disappearing.
With `Q=sigma^12`, the face equations, up to torus monomials, are

```text
G_M=(S^2-P)C,             C=1-UP^6-WS^2P^5,             (7)
G_T=S^2(1-S^2P^2A(P)),
A(P)=K+Theta P+xi_10 P^2+W P^3.                         (8)
```

The projected polygons and Pick ledgers are

| object | vertices | `(2Area,boundary,interior)` |
|---|---|---:|
| `M` | `(0,1),(2,0),(4,5),(0,7)` | `(36,10,14)` |
| `C` | `(0,0),(0,6),(2,5)` | `(12,8,3)` |
| `T` | `(2,0),(4,2),(4,5)` | `(6,6,1)` |
| global | `(0,1),(2,0),(4,2),(4,5),(0,7)` | `(42,14,15)` |

The complete nonmonomial boundary list, followed by the internal `M--T`
edge, is

```text
X-1,       1-KX^2,       A(X),       (X-1)(UX+W),
U-X^6;                                      internal: 1-WX.       (9)
```

Every entry except the deliberately variable `A` is active and squarefree
under `(2)`. The top quadratic collides only when `U+W=0`. Since `K*W!=0`,
all roots of `A` are finite and nonzero and cannot hit its toric endpoints.

The birational substitution `x=P^(-1), y=WS/P` turns `C` into

```text
y^2=W x(x^6-U).                                             (10)
```

It is smooth, connected, and genus three under `UW!=0`. The rational factor
`R:S^2=P` meets it at the twelve simple roots of `1-Lambda S^12`; the
primitive shared `C--T` edge contributes one more attachment. Before
refining a cubic collision, the graph has

```text
b1=13-3+1=11.                                               (11)
```

The global edge packet `(11,11,3,3,3,2,2,1)` has ramification sum `28`,
equal to `2*15-2`.

## 3. The three-tax Laurent root lemma

Let `A` be any nonzero polynomial of degree at most three and put

```text
d=deg A,        z=ord_(P=0)A,        Q_A=P^(-z)A,
tau=deg gcd(Q_A,Q_A').                                      (12)
```

If `N_tor` is the number of distinct roots of `Q_A` in `G_m`, then

```text
N_tor=d-z-tau,
3-N_tor=(3-d)+z+tau.                                       (13)
```

The three summands are, respectively, root multiplicity exiting at infinity,
root multiplicity exiting at zero, and collision multiplicity retained in
the Laurent torus. Moreover

```text
tau=length k[P,P^(-1)]/(Q_A,Q_A')
   =ord_(t=0) Disc_P(Q_A+t).                               (14)
```

Thus the collision tax is the localized critical/Jacobian-algebra length.
The polynomials `P(P^2+P+1)` and `(P-1)^2(P+1)` each have two distinct
toric support points, but the first pays a zero-exit tax and is reduced while
the second pays an interior-collision tax and is nonreduced. Support size or
the unsaturated discriminant alone loses this distinction.

On `(2)`, both endpoint taxes vanish. Put

```text
D0=xi_10^2-3WTheta,
D1=2xi_10^3-9Wxi_10Theta+27W^2K.                          (15)
```

Then

```text
4D0^3-D1^2=27W^2 Disc(A),                                (16)

Disc(A)=xi_10^2Theta^2-4WTheta^3-4xi_10^3K-27W^2K^2
        +18Wxi_10Theta K.
```

The exhaustive strata are

```text
Disc(A)!=0:                  three simple roots;
Disc(A)=0, D0!=0:            a double root plus a simple root;
D0=D1=0:                     one triple root.             (17)
```

On the double stratum,

```text
A=W(P-a)^2(P-c),              a*c*(a-c)!=0,               (18)
a=(9WK-xi_10Theta)/(2D0),     c=-xi_10/W-2a.
```

For any perturbation `g`, direct differentiation gives

```text
[epsilon]Disc(A+epsilon g)=4W^3(c-a)^3 g(a).              (19)
```

Evaluation at the labelled repeated root is therefore the intrinsic
conormal coordinate. The double critical algebra has length one; the triple
critical algebra has length two, so the latter cannot be compressed to one
scalar response coordinate.

## 4. Exact reciprocal chart and differential

Write

```text
J(P)=Phi+eta P+alpha_11 P^2,             B(P)=P^3J(P),
C_0(P)=-3P+(8/3)P^2-(1376/135)P^3
       +Delta P^4+upsilon_5 P^5+UP^6.                    (20)
```

In the `T` chart put `Q=sigma^12`, `s=sigma^(-6)S`, `p=P`, `y=sp`, and
`delta=sigma^6`. Multiplying the exact source equation by `sigma^12` gives

```text
G=(S^2-delta^2P)
  (1-S^2P^2A-delta S B-delta^2C_0)-delta^2S^2/2.          (21)
```

At the cubic edge set `b=1/S` and multiply by `b^4`:

```text
F=(1-delta^2Pb^2)
  (b^2-P^2A-delta bB-delta^2b^2C_0)-delta^2b^2/2.         (22)
```

There is no source ellipsis. Since `F_bb(0,a,0)=2`, formal Morse
preparation gives a unique critical section. At a repeated root,

```text
b_crit(a,delta)=delta B(a)/2+O(delta^2),
kappa(a,delta)=-B(a)^2delta^2/4+O(delta^3).               (23)
```

If `B(a)=0`, `(b,P)=(0,a)` is an exact singular section for every `delta`:
`F,F_b,F_P` vanish there identically. There is no later scalar Hasse splitter;
the correct object is a horizontal conductor.

The inherited target identity becomes, up to a unit and sign,

```text
phi^*eta_0=-sigma^16 b^2 db/F_P.                          (24)
```

Indeed `F=b^4G`, so `G_P=b^(-4)F_P` on the curve and
`dS=-b^(-2)db`. If a carrier uses

```text
b=sigma^hY,        P-a=sigma^rX,        F=sigma^D F_0,
```

its good-form order is

```text
16+3h-(D-r).                                              (25)
```

The primary face orders are `C:9` and `T:16`.

## 5. Squarefree and double-root cases

If `A` is squarefree, the `T` function field is `Y^2=A(P)`, a smooth
genus-one curve. The ledger is `3+1+11=15`, and the nonrational components
have orders nine and sixteen.

Now take `(18)` and set `x=P-a`.

### 5.1 `B(a)!=0`

With `b=sigma^6Y`, `x=sigma^6X`, division of `(22)` by `sigma^12` gives

```text
Y^2-B(a)Y-a^2W(a-c)X^2=0.                                (26)
```

This is a smooth projective conic. Completing the square before scaling gives

```text
(b')^2-a^2W(a-c)x^2
 =sigma^12(B(a)^2/4+O(sigma^6,x)).                        (27)
```

Thus the total germ is an `A_11` smoothing whose multiplicity-one rational
path joins the two points of the normalized `T`. It raises the graph Betti
number to twelve:

```text
g=3+12=15.                                                (28)
```

This corrects the provisional uniform `b1=11` assertion rejected by the
hostile referee. Formula `(25)` gives conic order
`16+3*6-(12-6)=28`; rationality already makes its target map constant.

### 5.2 `B(a)=0`

Formal preparation of `(22)` gives

```text
z^2=x^2u(sigma,x),                     u(0,0)!=0.          (29)
```

Extracting the unit square root separates two regular horizontal sheets.
No vertical component appears. The raw toric arithmetic genus is still
fifteen, but the persistent node has delta one on every fibre, so the
normalized family has

```text
g=3+11=14.                                                (30)
```

## 6. Triple-root cases

Let `A=W(P-a)^3`, `a!=0`, and put `x=P-a`.

### 6.1 `B(a)!=0`: the elliptic tail

With `b=sigma^6Y`, `x=sigma^4X`, the leading carrier is

```text
Y^2-B(a)Y-a^2WX^3=0,
(Y-B(a)/2)^2=a^2WX^3+B(a)^2/4.                            (31)
```

Its projective closure is a smooth `j=0` elliptic curve. It attaches once,
so `g=3+1+11=15`. The good target is also `j=0`; the load-bearing fact is

```text
ord_tail=16+3*6-(12-4)=26>0.                              (32)
```

Hence the elliptic-tail map is constant.

### 6.2 `B(a)=0`, `B'(a)!=0`: a rational nodal cubic

Put `b=sigma^18Y`, `x=sigma^12X` and write `B_1=B'(a)`. An arbitrary
quadratic and higher remainder in `B` does not change the face

```text
Y^2-B_1XY-a^2WX^3=0.                                      (33)
```

This cubic is **not elliptic**. It is singular at `(0,0)`, and

```text
(Y-B_1X/2)^2=X^2(B_1^2/4+a^2WX).                         (34)
```

Its normalization is rational; its arithmetic genus and delta are both one.
If retained, its order is `16+3*18-(36-12)=46`.

### 6.3 `ord_a B>=2`: the persistent cusp

On the cusp scale `v(x)=r`, `v(b)=3r/2`, the normal term has excess

```text
6+(ord_a(B)-3/2)r >= 6+r/2>0.                            (35)
```

It never creates a later balanced carrier. The prepared discriminant is
`x^3` times a unit; this is a horizontal `(2,3)` cusp with rational
simultaneous normalization.

## 7. Simultaneous normalization and the corrected ledger

All `B(a)=0` cases normalize over the DVR, not just fibrewise. Formal Morse
preparation writes `z^2=D(x,sigma)` with the exhaustive forms

```text
double:                 D=x^2u,                 u a unit;
triple, ord_a B=1:      D=x^2E,                 E_x(0,0)!=0;
triple, ord_a B>=2:     D=x^3u,                 u a unit. (36)
```

The first splits after a unit square root. In the second, the integral
element `w=z/x` satisfies `w^2=E`, a regular equation because `E_x` is a
unit. In the third, after a unit square root the normalization is

```text
x=t^2,                         z=t^3sqrt(u).              (37)
```

These constructions are finite, birational, regular, and torsion-free over
the DVR. Excellence globalizes finiteness; finite over proper preserves
properness; DVR torsion-freeness gives flatness. Thus no vertical component
or hidden genus is lost at the conductor.

The complete corrected inventory is

| cubic/normal stratum | nonrational carrier | graph `b1` | normalized genus | order |
|---|---|---:|---:|---:|
| squarefree | `T`, genus 1 | 11 | 15 | `16` |
| double, `B(a)!=0` | none; rational bridge | 12 | 15 | `28` |
| double, `B(a)=0` | none; horizontal node | 11 | 14 | -- |
| triple, `B(a)!=0` | elliptic tail | 11 | 15 | `26` |
| triple, `ord_aB=1` | none; rational nodal cubic | 11 | 14 | `46` |
| triple, `ord_aB>=2` | none; horizontal cusp | 11 | 14 | -- |

In every row `C` has order nine. A squarefree `T` has order sixteen and the
only new positive-genus tail has order twenty-six. Their maps to the good
elliptic target are constant; every other component is rational and also
maps constantly.

After finite base change, normalization, and regularization, proper-flat
degree conservation for a positive-degree target line bundle gives

```text
deg(phi_generic^*L)=sum_i m_i deg(phi_i^*L)=0.            (38)
```

This contradicts the positive generic response degree of a nonautomorphic
Keller pair and proves `(1)--(2)`. **QED.**

## 8. Unexpected connections and stopping boundaries

Equation `(13)` is a precise version of the repository's indexed-sequence
lesson. Mapping a labelled packet to the natural number `N_tor` loses why an
address disappeared; the minimal sidecar is

```text
(infinity exit, zero exit, interior Jacobian length).      (39)
```

Two further organizers are exact but do not enter the proof graph:

- modulo polynomial and root scaling, a double-root cubic is
  `A_q=(P-1)^2(P-q)`. Toric inversion sends `q` to `1/q`; `q=1` is triple
  and `q=0,infinity` are endpoint exits. Stern--Brocot words enumerate the
  positive rational part of this slice, but lose the scales, normal Hasse
  coefficient, and surrounding faces;
- signs of real root differences produce only a transitive tournament.
  A collision is a tie, and one Seidel switch creates a directed triangle
  that no one-dimensional root-difference field realizes. The lawful finite
  carrier is the `S_3` labelled-root torsor plus endpoint flags and the local
  Jacobian algebra.

The remaining sharp mechanisms are deliberately separate:

```text
K=0:           a root exits through P=0;
W=0:           degree drops and a root exits through P=infinity;
U+W=0:         the independent top A_23 contact reappears. (40)
```

For `W!=0`, the exact zero-exit filtration `m=ord_0A=1,2,3` has global Pick
ledgers

```text
m=1: (40,12,15),       m=2: (38,12,14),
m=3: (36,10,14).                                            (41)
```

These are next-task signals, not an extinction extension. The next audits
are: combine `(36)` with the `K=0` fan; build the `W=0` degree-drop atlas;
merge cubic collision with the `U+W=0` ladder without importing simple-edge
responses; and independently referee the new `U=0` `(2,5)/(2,3)`
conductor-Morse candidate.

After this theorem the exact `Z=0,U!=0` residual remains

```text
beta_11=0 and W*zeta_3=0,                                (42)
```

but its subgate `zeta_3=0,W*K*(U+W)!=0` is closed.

## 9. Reproduction and honest scope

Run all three certificates in normal and optimized modes and byte-match the
frozen outputs:

```bash
python3 -B 04-computation/jc2_m12_clean_cubic_edge_extinction_thm4339.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_edge_extinction_thm4339.py
python3 -B 04-computation/jc2_m12_clean_cubic_edge_extinction_independent_audit_thm4339.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_edge_extinction_independent_audit_thm4339.py
python3 -B 04-computation/jc2_m12_clean_cubic_edge_extinction_hostile_referee_thm4339.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_edge_extinction_hostile_referee_thm4339.py
```

The primary path reuses only the pinned THM-4327 hull engine. The independent
path reconstructs `(21)--(22)` directly from the displayed source. The
referee targets the false-genus hazards in the general simple-zero triple
cubic and the double horizontal tangent discriminant.

What is proved is exactly `(1)--(2)`, relative to the inherited seam and
proper-flat interfaces. Equations `(39)--(41)` are organizers and next-task
signals, not extensions of the theorem gate.
