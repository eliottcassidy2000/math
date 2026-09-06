# The valuation-eleven response is constant across its varying background

**Status: PROVED FINITE-ROW + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The complete source response with valuation at least eleven has a sharp
weight-22 replacement threshold through row fifteen. Its weight-22 source
family over the actual boundary is `G_m x A^9`, with affine ten-dimensional
terminal fibres. The sixth background row genuinely varies along `G_m`,
but the entire residual response matrix is constant in that varying
coefficient. A polynomial inverse certificate proves this at every
specialization, including values outside the actual boundary chart.
Source valuation below eleven, later equations, termination, chart entry,
and JC(2) remain open.

The [independent audit](planar_jc_long_20260906_memory_earlier_audit.md)
recovers all seven equations in 81 tangent coordinates, proves a separate
polynomial minor certificate, and verifies every affine family direction
and the ten-dimensional coordinate kernel. Its 5,868 always-active checks
pass in byte-identical normal/optimized runs.

## 1. Inheritance and the first nonconstant background row

The closest proved mechanism is the
[valuation-twelve classification](planar_jc_long_20260906_memory_valuation12.md),
which enlarged the source family to `G_m x A^5` without lowering its sharp
replacement weight. Its ancestor is the
[complete valuation-thirteen classification](planar_jc_long_20260906_memory.md).
The nonzero weight-22 transport itself was already proved in the incoming
[compensated transport](planar_jc48_sep06_weight14.md); it is a positive
control here, not a new existence claim.

Targeted exact-constant and statement searches recovered the earlier
row-eleven continuation and the valuation-twelve result; neither classifies
the complete moving-prefix universe below. No external priority claim is made.

The canonical hostile is failure to replace the high channel at weight 21.
The corrected near miss would be to infer a statement on every boundary
point from a rank over the rational function field. The least-used sidecar
is the actual sixth background row, including its source parameter. The
five-concept board is source valuation, residual weight, raw depth
constraints, varying background, and nilpotent correction of an inverse.
The map changes unknown rows starting at ten and source rows starting at
eleven; it retains actual monomial coefficients and every finite equation.
The sidecar is the full low background, rather than an unlabelled sample.

Use the coordinates and depth modules of
[THM-4308](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md):

```
p=t(1+x^2t), y=xtp, u=x^2t,
P(A,C)=C^2-A^3+(3/4)A+1/4,
P_d=span{x^a u^b p^c y^e:a+b<=d}.
```

Rows zero through four are precisely the constant rows displayed in the
valuation-twelve report. The next actual rows, on
`Phi=eta=0,alpha11=1,Delta=896/15`, are

```
A5=128x^6/9-8576x^4/75-(6/11)Xi*x^2-203776x^2/7425
   -x/2+3528704/6075,
C5=-32x^7/3-128x^5/5+(9/22)Xi*x^3+130816x^3/275
   +3x^2/8+(9/11)Xi*x-4024832x/4455+3/4,              (1)
```

where `Xi=xi10`. This coefficient is not constant on the actual boundary.
Reconstructing the source graphs through row twelve of
[THM-4426](../../01-canon/theorems/THM-4426-source-normal-row-fourteen-weight-eighteen-memory-repair.md)
at `c51=1087/135` gives

```
Xi=(108391820625 h+3765431711250 z+324974300165767168)
     /24542012448000.
```

In that theorem's actual parameter
`s=801h+27826z+85855050266495746048/37533020625`, this simplifies to

```
Xi(s)=(1485/269322496)s+34969998848/55604475.           (2)
```

Thus `Xi` is itself an affine coordinate on the boundary `G_m`, omitting
only `34969998848/55604475`. Rows zero through four cannot distinguish its
points, while the `x^2` coefficient of row `A5` already recovers `Xi` and
hence `s`. The constant response quotient below therefore loses a background
coordinate which is visible in the actual rows; the polynomial inverse
retains that coordinate when constructing lifts.

In particular, the new row really sees a varying background. The proof
below retains `Xi` as an arbitrary field element; it does not specialize it
to zero, to a rational test point, or to a generic point without auditing
the excluded specializations.

The [background reconstruction](../../04-computation/planar_jc_long_20260906_memory_earlier_background.py)
pins the audited THM-4426 script. Specialization changes six *primitive
integer normalizations*. Each affected gate is replaced by an explicit
nonzero rational pivot gate, rather than suppressed. The reconstruction
checks all raw row-nine through row-twelve compatibilities after substitution,
reconstructs (1), and proves (2). No later equation is needed for this input.

## 2. Complete response classification

Work over any characteristic-zero field, with arbitrary `Xi` in (1).
Other background rows may be arbitrary. Allow perturbations supported in
rows ten through fifteen, with
`deg dA_n<=n+1`, `deg dC_n<=n+2`, and impose

```
J(A+dA,C+dC)-J(A,C)=0 mod t^15,
P(A+dA,C+dC)-P(A,C)=R(p,y) mod t^16,
pi_15(dA) in pi_15(P_2), pi_15(dC) in pi_15(P_3).      (3)
```

The complete monomial universe for `R_low`, of valuation at least eleven
and weight at most 23, is

```
valuation11: p^11, p^9y, p^7y^2, p^5y^3, p^3y^4, py^5;
valuation12: p^10y, p^8y^2, p^6y^3, p^4y^4, p^2y^5, y^6;
valuation13: p^7y^3, p^5y^4, p^3y^5, py^6;
valuation14: p^4y^5, p^2y^6, y^7;
valuation15: py^7.                                   (4)
```

Adjoin only the designated high control `k p^3y^6`, of weight 24.
Coefficients in the following equations are coefficients of the
*deformation source* `R`, not the total inherited Hamiltonian. Put

```
a=[p^11]R, b=[p^7y^2]R, c=[p^3y^4]R,
d=[p^8y^2]R, e=[p^4y^4]R,
f=[y^6]R, o=[py^5]R,
v=[p^5y^4]R, w=[py^6]R, r=[p^2y^6]R,
D=119712905708881603.
```

**Classification.** System (3) is solvable if and only if `f=o=0` and

```
D a= 736701988512897 v-156822363286128 w-128022192499536 r
     -(315619079794353/2) k,
D b=-(17667890668710639/2) v+8012657615180013 w
     -1940411140875144 r+(1696910148536463/2) k,
D c= 6603015506480190 v-(15756867503321415/2) w
     +2521182637658160 r-(2500583114637405/8) k,
D d=-40574567706597180 v-17040074980207200 w
     -5161415972895171 r+1095964361654112 k,
D e= 67697466072152616 v-13525083050179824 w
     +(37920301110828423/2) r-691333448130801 k.         (5)
```

The three even coefficients `v,w,r`, the high coefficient `k`, and all ten
odd coefficients other than `o` are free. Every admissible source has an
affine ten-dimensional perturbation fibre. All equations and dimensions
hold for **every** `Xi` in every characteristic-zero field.

At weight at most 22 the six free odd channels are exactly

```
p^9y, p^5y^3, p^6y^3, p^2y^5, p^3y^5, y^7.          (6)
```

Together with `v,w,r`, these give nine free low coefficients for each
prescribed high coefficient. At weight at most 23 the extra four odd
channels are `p^10y,p^7y^3,p^4y^5,py^7`, so the corresponding dimension
is thirteen. The conditions `f=o=0` remain essential. In particular the
previous coefficient of `y^6` is preserved; the theorem does not assert
that the total inherited `y^6` coefficient vanishes.

## 3. The sharp threshold and actual boundary family

The exact ranks of the low source image, and after adding the high column,
are

| Weight bound | Low response rank | Rank with high column |
|---:|---:|---:|
| 17 | 1 | 2 |
| 18 | 3 | 4 |
| 19 | 3 | 4 |
| 20 | 6 | 7 |
| 21 | 6 | 7 |
| 22 | 7 | 7 |
| 23 | 7 | 7 |

Thus weight 21 cannot replace a nonzero multiple of `p^3y^6` in the declared
valuation filter, while weight 22 can. The earlier normalized positive
control is retained exactly:

```
L=-(27945/235202)p^5y^4
  +(39123/470404)py^6+(52578/117601)p^2y^6;
R=L-p^3y^6 satisfies (5).                             (7)
```

The change from valuation twelve to eleven therefore enlarges the complete
family, but does not lower the sharp threshold.

Let `H_pre(s)+j(s)p^3y^6` be the actual positive source of
[THM-4438](../../01-canon/theorems/THM-4438-jc2-row-fifteen-relative-response-on-boundary-gm.md).
Take every low source `R_low` obtained from (5) with `k=-j(s)`, using only
weight-22 monomials. Then precisely the replacement sources in this filter
are

```
H_pre(s)+j(s)p^3y^6+R_low-j(s)p^3y^6
 =H_pre(s)+R_low.                                    (8)
```

The nine free coefficients described in (6) and after (5) identify their
source space with `G_m x A^9`; unchanged `z,h` recover `s`. All unknown rows
zero through nine are fixed to the inherited boundary. For each source,
the complete row-fifteen fibre relative to that prefix is affine
ten-dimensional. Formula (2) introduces no rank divisor: the universal
certificate below applies to its entire image and even to its missing
value at `s=0`. At a zero of `j(s)`, no nonzero payer is needed and the same
affine classification remains valid. At every rational `s!=0`, THM-4438's
irreducible quartic gives `j(s)!=0`, so the sharp replacement statement
applies.

## 4. Why no background specialization is lost

The producer begins with all 180 coordinate coefficients of `dA,dC`, not a
chosen tangent ansatz. Expanding the literal variations in (3) and the
complete 91-row depth bank gives a homogeneous system

```
M(Xi) q+S r_source=0,
```

with 301 raw equations and 21 source coordinates. Exact elimination over
`Q(Xi)` first finds tangent rank 170 and source rank seven. This computation
alone would only justify a generic claim. The following additional
polynomial identities remove that gap.

Choose the displayed 170 tangent pivot columns `M_c` and a fixed set `I`
of 170 rows selected at `Xi=0`. Its square minor is

```
B(Xi)=B0+Xi B1=B0(I+Xi K), K=B0^(-1)B1.
```

The certificate constructs `B0^(-1)` over `Q` and verifies

```
K^2=0,
B(Xi)^(-1)=(I-Xi K)B0^(-1),
M=M_c B(Xi)^(-1) M_I.                                (9)
```

Thus the full tangent rank is 170 at every specialization, with ten free
unknown coordinates and no nonconstant denominator. After removing these
columns, the *entire raw residual response matrix*

```
F=S-M_c B(Xi)^(-1)S_I                                (10)
```

is constant in `Xi`. Its seven-row reduced matrix `E` is exactly (5),
including `f=o=0`. The certificate checks `F=F[:,P] E`, where `P` is the
seven source pivot columns. A fixed seven-by-seven minor has determinant

```
122586015445894761472/7942904748804375 != 0.           (11)
```

Consequently these seven rows generate the full source condition space
for every `Xi`. Equations (9)–(11) prove both necessity and sufficiency over
`Q[Xi]` and after every characteristic-zero specialization. They exhibit the
mechanism: variation of this background coefficient changes the choice of
lift by a nilpotent correction, while leaving all residual measurements
unchanged.

The earlier filtered superposition principle is exact here as well.
Corrections start in row ten, so nonlinear terms in `P` begin at row twenty,
and in the Jacobian at row nineteen. Through the orders in (3), the literal
variation is therefore linear in the correction. Only the background
through row five can enter. The linear classification is consequently an
exact finite-parameter theorem, not merely a tangent-space statement.

## 5. Hostiles, reproduction, and the next boundary

Parity alone fails as a response criterion. The odd channel `py^5` is
jointly neutral through row thirteen, but the complete row-fourteen system
forces its deformation coefficient to zero. This delayed failure is checked
by two fresh prefix systems. Later odd channels in (6) remain neutral.
The raw joint constraints identify this distinction; a scalar parity test
would miss it. The old `y^6` deformation condition is also retained, and
all weight-21 replacement attempts fail by the exact rank increase above.

The [complete response source](../../04-computation/planar_jc_long_20260906_memory_earlier.py)
imports no repository mathematical implementation. Its
[output](planar_jc_long_20260906_memory_earlier.out) retains the complete
monomial list and relation matrix, exact polynomial inverse tests, source
minor, weight ranks, all 301 original equations after the positive lift,
and the odd-channel onset. The separate background reconstruction has its
[own output](planar_jc_long_20260906_memory_earlier_background.out) and pinned
inherited implementation. Reproduce both with

```
python3 -B 04-computation/planar_jc_long_20260906_memory_earlier.py
python3 -B -O 04-computation/planar_jc_long_20260906_memory_earlier.py
python3 -B 04-computation/planar_jc_long_20260906_memory_earlier_background.py
python3 -B -O 04-computation/planar_jc_long_20260906_memory_earlier_background.py
```

The precise next frontier is valuation ten, where corrections start at row
nine and see the seventh actual background row. Its source-field dependence
and complete response have not been assumed or tested here. Nor do these
finite additive translations assert a polynomial flow, a full depth lift,
termination, or a Keller pair. The incoming Hamiltonian realization and
non-local-nilpotence boundary remain separate, inherited results.

Both normal/optimized pairs byte-match. The response certificate has **329**
counted always-active gates; the background reconstruction has **276**,
including its 210 raw compatibility checks. SHA256 pins are

```
response source fe00f61262d3ca006e2582fbe49aef07a62b5843fc0089030abccb3b2618aa5d
response output e992cf8e2be79f60fd10ff8fe126b9b6d994862ec86a37c19b5e587fc5fa43f4
background source 3bfa772f2116e6c7d9cf4f43e10d865bfe99f24e0276e706c807b63f4b57b379
background output 628bf9e473efc2b06b95976c0029815ab20409d06b677b34bf159cb2a6edf136
response semantic 892042bb9a634c4201dacbb935df0d331f2c9624a7fc23c443ca3becc8c47f5d
```
