---
id: THM-4342
title: "Clean cubic zero-exit planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327/4339 + VERIFIED-EXACT + TWICE
  HOSTILE-AUDITED. In the inherited reduced (2,3), exact-weight-twelve seam,
  the gate Z=beta_11=zeta_3=K=0 and U*W*(U+W)!=0 is extinct; K=0 forces
  Delta=5696/105. Laurent saturation separates zero-exit depth from the sole
  possible interior collision. The three exit depths have normalized genera
  15,14,14; only depth one retains a possible elliptic T face, of Keller-form
  order 16. Its double-root degeneration is a rational A11 bridge or a
  horizontal node. All complementary blowup charts are explicit. W=0,
  U+W=0, seam entry, JC(2), and DC(2) are not claimed here.
source: root + cubic-root/second-referee agents / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4339-clean-interior-cubic-edge-planar-jacobian-extinction
related:
  - THM-4340-u-zero-repeated-cusp-planar-jacobian-extinction
  - THM-4341-odd-self-similar-cusp-reciprocal-tail-duality
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_thm4342.py
primary_output: 05-knowledge/results/jc2_m12_clean_cubic_zero_exit_extinction_thm4342.out
primary_script_sha256: 55dd4f19c405ab656a803ada884ed4a97c4fbb6ef13cb63b489f1b6343807214
primary_output_sha256: 7ced31516748487151e3ec472c5429190e5de9f95b632a4aad30e90073dbc5dc
hostile_referee_script: 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_hostile_referee_thm4342.py
hostile_referee_output: 05-knowledge/results/jc2_m12_clean_cubic_zero_exit_extinction_hostile_referee_thm4342.out
hostile_referee_script_sha256: afca0243d8f9299cc907dd31eb685ba878a68d5f20274bc4fdfe2f1c5671b460
hostile_referee_output_sha256: 524bc6beb3aaeb6dfc0d8202383a3adf434a95a1de31f47d7639016999da487d
second_referee_script: 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_second_referee_thm4342.py
second_referee_output: 05-knowledge/results/jc2_m12_clean_cubic_zero_exit_extinction_second_referee_thm4342.out
second_referee_script_sha256: 182af72ac5066bc944bcf687b288deb0d29be3f8bad570753f5f344999c65249
second_referee_output_sha256: 8d6d8f0bb96e9aa8f4d5e050dda40862582812a5dd8518abeb4526f0a9718da0
hash_basis: raw LF bytes
audit: >
  PASS AFTER SCOPE REPAIR. The pinned primary hull path exhausts 4,096,
  2,048, and 2,048 conservative support masks at exit depths one, two, and
  three; checks the saturated fan, Pick/Riemann-Hurwitz ledgers, exact strict
  transforms, normal Hasse gate, form transport, and graph genera. The first
  import-free referee reconstructs the local chart and both double-root
  alternatives. The second import-free referee catches the forced Delta
  value and verifies every complementary ordinary-blowup chart and glue.
  Pick and Riemann-Hurwitz are corroborating checks; completeness comes from
  the exhaustive fan and chart coverage. Normal and optimized streams
  byte-match all three frozen outputs.
---

# THM-4342 -- Clean cubic zero-exit planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327/4339 + VERIFIED-EXACT + TWICE
HOSTILE-AUDITED. THE DISPLAYED `K=0` GATE IS EXTINCT. `W=0`, `U+W=0`,
SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OUTSIDE THIS THEOREM.**

## 1. Statement and inheritance

Work over an algebraically closed field of characteristic zero in the
inherited reduced `(2,3)`, exact-weight-twelve seam. Impose

```text
Z=beta_11=zeta_3=K=0,                  U*W*(U+W)!=0.     (1)
```

The seam relation is load-bearing:

```text
K=2848/45-(7/6)Delta,
K=0  ==>  Delta=5696/105!=0.                              (2)
```

Every other allowed lower coefficient is arbitrary. No nonautomorphic
planar Keller pair lies on `(1)--(2)`. The claim is relative to the inherited
seam and proper-flat target interface; it proves neither seam entry nor
`JC(2)`.

The inheritance pass was:

- closest mechanism: THM-4339's reciprocal cubic chart, local collision
  forms, good-differential orders, and horizontal normalization;
- canonical hostile: `A=P^2(P+1)`, whose full discriminant vanishes although
  its Laurent-saturated factor is squarefree;
- corrected near miss: THM-4339's theorem cannot be imported across its
  failed `K!=0` hypothesis; the endpoint fan and all blowup charts must be
  proved anew;
- least-used sidecar: the toric strict transform at the exited root.

The live board was

```text
exit depth | saturated polynomial | endpoint fan | complementary chart
interior normal jet | conductor | good form | proper-flat degree.        (3)
```

## 2. Honest zero-exit fan

On the cubic edge put

```text
A(P)=Theta P+xi_10 P^2+WP^3=P^m Q_m(P),                 (4)
```

where, because `W!=0`, exactly one case holds:

```text
m=1: Theta!=0,              Q_1=Theta+xi_10P+WP^2;
m=2: Theta=0, xi_10!=0,     Q_2=xi_10+WP;
m=3: Theta=xi_10=0,         Q_3=W.                       (5)
```

Using the full sixteen-row source and independently deleting the inherited
multiply-owned aggregate points gives the conservative exact atlas:

| `m` | required exit owner | hostile masks | lower faces |
|---:|---|---:|---|
| 1 | `Theta` | `4096/4096` | `M,T` |
| 2 | `xi_10` | `2048/2048` | `M,T` |
| 3 | none | `2048/2048` | `M` only |

The enumeration deliberately also allows deletion of the `Delta` owner.
Those masks are impossible under `(2)`, but the overcount is conservative:
every mask, including the actual `Delta`-present subfamily, has the displayed
complex.

The honest polygons are

| `m` | `T` vertices; `(2Area,B,I)` | global vertices; ledger |
|---:|---|---|
| 1 | `(2,0),(4,3),(4,5)`; `(4,4,1)` | `(0,1),(2,0),(4,3),(4,5),(0,7)`; `(40,12,15)` |
| 2 | `(2,0),(4,4),(4,5)`; `(2,4,0)` | `(0,1),(2,0),(4,4),(4,5),(0,7)`; `(38,12,14)` |
| 3 | no two-dimensional `T` | `(0,1),(2,0),(4,5),(0,7)`; `(36,10,14)` |

Their edge packets are

```text
m=1: (11,11,5,3,3,1),
m=2: (11,11,3,3,3,1),
m=3: (11,11,7,1),                                      (6)
```

with ramification sums `28,26,26=2g-2`. These Pick and
Riemann--Hurwitz identities are exact checks, not substitutes for the fan and
chart completeness below.

## 3. Laurent saturation and boundary schemes

The complete outer edge lists, followed by the internal edge when present,
are

```text
m=1:
X-1, 1-Theta X, Theta+xi_10X+WX^2,
(X-1)(UX+W), U-X^6; internal 1-WX;

m=2:
X-1, 1-xi_10X^2, xi_10+WX,
(X-1)(UX+W), U-X^6; internal 1-WX;

m=3:
X-1, 1-WX, (X-1)(UX+W), U-X^6.                         (7)
```

Every endpoint coefficient is nonzero under `(1)` and `(5)`. The top
quadratic is squarefree because `U+W!=0`, and `U-X^6` is squarefree because
`U!=0`.

Only `m=1` can have an interior collision:

```text
Disc(A)=Theta^2(xi_10^2-4WTheta),
Disc(Q_1)=xi_10^2-4WTheta.                              (8)
```

The factor `Theta^2` is the zero-exit tax. Treating it as an interior
collision is exactly the error exposed by `P^2(P+1)`. The saturated `Q_2`
and `Q_3` are linear and constant.

## 4. Exact strict transforms and complete blowup coverage

Let `delta=sigma^6` and put

```text
J(P)=Phi+eta P+alpha_11P^2,
D(P)=-3+(8/3)P-(1376/135)P^2+Delta P^3+upsilon_5P^4+UP^5.
```

The exact reciprocal equation is

```text
F=(1-delta^2Pb^2)
  (b^2-P^(m+2)Q_m-delta bP^3J-delta^2b^2PD)
  -delta^2b^2/2.                                         (9)
```

Set `e_1=1`, `e_2=e_3=2`, and take the `P`-chart `b=P^e y`. Dividing by
`P^(2e)` gives

```text
E_m=(1-delta^2P^(2e+1)y^2)
    (y^2-P^(m+2-2e)Q_m-delta P^(3-e)yJ-delta^2Py^2D)
    -delta^2y^2/2.                                      (10)
```

Its special models are

```text
m=1: y^2=P(Theta+xi_10P+WP^2),
m=2: y^2=xi_10+WP,
m=3: y^2=WP.                                             (11)
```

For `m=1`, `(E_1)_P(0,0)=-Theta`, so one ordinary blowup resolves the exit.
For `m=2`, the second `P`-chart has exceptional equation

```text
(1-delta^2/2)y^2=xi_10,                                 (12)
```

giving two simple points. For `m=3`, `(E_3)_P(0,0)=-W`; `(11)` is an
optional-refinement rational binomial collar.

The complementary charts close the coverage gap. On the first blowup put
`P=bz` and divide by `b^2`. For all three depths its exceptional equation is

```text
1-delta^2/2=0,                                           (13)
```

which is empty near the DVR origin because the left side is a unit. At the
second blowup for `m=2,3`, the complement of `y_1=Py` is `P=y_1z`. For
`m=2` its exceptional equation is

```text
(1-delta^2/2)-xi_10z^2=0,                               (14)
```

and `z=1/y` glues its two points to `(12)`. For `m=3` the complementary
equation is again the unit `(13)` and has no point. Thus no endpoint or
exceptional component is omitted. The remaining toric subdivision pieces
are rational.

There is no delayed exit carrier. On the balance
`b^2~P^(m+2)Q_m(0)`, the first normal term `delta bP^3J` has excess

```text
m=1: 6+3v(P)/2,
m=2: 6+v(P),
m=3: 6+v(P)/2,                                          (15)
```

in sigma units. Every value is positive; extra vanishing of `J` only delays
the term.

## 5. The sole interior collision

For `m=1`, if `Disc(Q_1)!=0`, `(11)` is a smooth genus-one double cover. If
`Disc(Q_1)=0`, write

```text
Q_1=W(P-a)^2,              a=-xi_10/(2W)!=0.             (16)
```

Its normalization `v=y/(P-a)`, `v^2=WP`, is rational. The invariant first
normal gate is

```text
N=4W^2Phi-2Wxi_10 eta+xi_10^2 alpha_11=4W^2J(a).        (17)
```

If `N!=0`, the scale `P=a+sigma^6X`, `b=sigma^6Y` gives

```text
Y^2-B(a)Y-a^3WX^2=0,                  B=P^3J.            (18)
```

This is a smooth rational conic and the `A_11` bridge of THM-4339. Its
good-form order is `28` and it restores one graph cycle.

If `N=0`, `(b,P)=(0,a)` is an exact horizontal singular section. With
`x=P-a`, the tangent discriminant is

```text
x^2[delta^2B'(a)^2+4q_0a^3W],
q_0=1-delta^2(C_0(a)+1/2).                              (19)
```

The bracket is a DVR unit. Formal preparation gives `z^2=x^2u`, and a unit
square root separates two regular horizontal sheets. The normalization is
finite, proper after globalization, and DVR-flat; it produces no vertical
tail. This imports THM-4339's local conductor mechanism, not its failed
`K!=0` global theorem.

## 6. Differential transport and complete genus ledger

The inherited exact form in the reciprocal chart is, up to a unit and sign,

```text
phi^*eta_0=-sigma^16 b^2 db/F_P.                         (20)
```

For `b=P^e y`, write `F=P^(2e)E(P,y)`. Differentiation at fixed original
`b`, followed by `dE=0`, gives exactly

```text
b^2db/F_P=P^e y^2dy/E_P.                                (21)
```

Thus no new sigma power is introduced. The only positive-genus endpoint
model is the squarefree `m=1` face, whose order remains `16`. The unchanged
central component is

```text
Y^2=WX(X^6-U),                  genus 3, order 9,         (22)
```

and is smooth under `UW!=0`.

The twelve `R--C` nodes remain simple because their scheme is
`1-(U+W)S^12`. The `C--T` edge is primitive for `m=1,2`; at `m=3` its
scheme belongs to the outer boundary. Hence `b1=11`, except that the
smoothed double bridge raises it to twelve. The exhaustive normalized ledger
is

| stratum | nonrational pieces | graph `b1` | genus |
|---|---|---:|---:|
| `m=1`, squarefree | `C:g3`, `T:g1` | 11 | 15 |
| `m=1`, double, `N!=0` | `C:g3`; rational bridge | 12 | 15 |
| `m=1`, double, `N=0` | `C:g3`; horizontal node | 11 | 14 |
| `m=2` | `C:g3`; rational `T` | 11 | 14 |
| `m=3` | `C:g3`; optional rational collar | 11 | 14 |

Fan exhaustion plus the complete charts in Section 4 prove component
completeness; `(6)` is an independent genus checksum. Every positive-genus
component has positive form order, and every other component is rational.
After finite base change, normalization, and regularization, the unchanged
proper-flat degree interface gives, with actual multiplicities,

```text
deg(phi_generic^*L)=sum_Gamma mu_Gamma deg(phi_Gamma^*L)=0. (23)
```

This contradicts the positive generic response degree of a nonautomorphic
Keller pair and proves `(1)`. **QED.**

## 7. Exact indexing lesson and next boundaries

At fixed cubic degree, `m=ord_0A` is a lawful natural-number index for the
zero-exit depth. It is not a complete degeneration label. The minimum
sidecar is

```text
(infinity exit, zero exit m, interior Jacobian length,
 saturated polynomial Q_m, toric edge owner).            (24)
```

Equation `(8)` explains why: the full discriminant multiplies the zero-exit
tax and the interior collision tax, while Laurent saturation separates them.
A Stern--Brocot label for a root of `Q_m` similarly loses the exit depth and
edge owner. No tournament on bare numerical root labels restores that data;
the native object is the edge-addressed root scheme with its local Jacobian
algebra.

The next sharp gates are genuinely different:

```text
W=0:       degree drops and a root exits at infinity;
U+W=0:     twelve R--C nodes merge into the A_23 contact;
U=0:       the sixfold edge degenerates.                                (25)
```

The first two have separate candidates under audit; `(25)` is not used to
enlarge this theorem.

## 8. Reproduction and scope

Run all three certificates in normal and optimized modes and byte-match the
frozen outputs:

```bash
python3 -B 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_thm4342.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_thm4342.py
python3 -B 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_hostile_referee_thm4342.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_hostile_referee_thm4342.py
python3 -B 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_second_referee_thm4342.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_zero_exit_extinction_second_referee_thm4342.py
```

The primary path pins only THM-4327's hull engine. Both referee paths are
import-free. What is imported from THM-4339 is limited to unchanged open
charts, the interior-double and differential mechanisms, horizontal
normalization, and the final degree interface. Endpoint completeness is the
new content of Sections 2--4.
