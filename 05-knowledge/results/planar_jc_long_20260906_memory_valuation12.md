# Complete valuation-twelve source response and its affine five-parameter replacement

**Status: PROVED FINITE-ROW + FINITE-EXACT + INDEPENDENTLY AUDITED.**
Allowing one earlier source valuation enlarges the weight-22 replacement
family from two to five parameters. The sharp minimum remains 22 in the
new declared universe. The classification holds on the entire inherited
boundary `G_m`, with an affine ten-dimensional terminal fibre. Source
valuation below twelve, later rows, polynomial termination, chart entry,
and `JC(2)` remain open.

## 1. Inheritance and the enlarged actual operator

The closest proved mechanism is the complete
[valuation-thirteen response classification](planar_jc_long_20260906_memory.md),
whose nonzero weight-22 transport was already recovered from the incoming
[compensated transport](planar_jc48_sep06_weight14.md). The canonical hostile
is its weight-21 failure. The corrected near miss is to extrapolate that
minimum or its two-parameter family after changing the fixed prefix. The
new sidecar is the actual fifth background row, which enters precisely
when the perturbation begins in row eleven. The five live concepts are
valuation, residual weight, the complete tangent space, full projected
depth, and nonlinear onset. This report recomputes their joint system.

Work over a characteristic-zero field `k`, with the source coordinates and
depth modules of
[THM-4308 / source-normal-bracket-hasse-truncation-through-row-eight](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md):

```text
p=t(1+x²t), y=xtp, u=x²t,
P(A,C)=C²-A³+(3/4)A+1/4,
J(A,C)=A_x C_t-A_t C_x,
P_d=span{x^a u^b p^c y^e:a,b,c,e>=0,a+b<=d}.
```

Fix rows zero through four of the background to

```text
A0=1+x²/4,                     C0=-3x/4-x³/8,
A1=4/3+2x²,                    C1=-4x-3x³/2,
A2=-32/9-4x²/5,                C2=88x/15-12x³/5,
A3=2176/135+64x²/9-32x⁴/15,
C3=-1408x/45+64x³/15+8x⁵/5,
A4=224x⁴/9-256x²/75-37376/405,
C4=-184x⁵/15-2048x³/25+98944x/675.                    (1)
```

These are the actual first five rows on the boundary of
[THM-4426 / source-normal-row-fourteen-weight-eighteen-memory-repair](../../01-canon/theorems/THM-4426-source-normal-row-fourteen-weight-eighteen-memory-repair.md)
and
[THM-4438 / jc2-row-fifteen-relative-response-on-boundary-gm](../../01-canon/theorems/THM-4438-jc2-row-fifteen-relative-response-on-boundary-gm.md).
They are independent of its `G_m` parameter. A separate inheritance check
for the newly needed row is given below. Background rows from five onward
can otherwise be arbitrary for the response theorem.

Consider all perturbations supported in rows eleven through fifteen,
with `deg_x delta A_n<=n+1`, `deg_x delta C_n<=n+2`. Impose the entire system

```text
J(A+delta A,C+delta C)-J(A,C)=0 mod t^15,
P(A+delta A,C+delta C)-P(A,C)=R(p,y) mod t^16,
pi_15(delta A) in pi_15(P_2),
pi_15(delta C) in pi_15(P_3).                         (2)
```

The allowed low source has valuation at least twelve and weight at most
23, where `val(p^a y^b)=a+2b` and `wt(p^a y^b)=2a+3b`. Its complete
monomial list is

```text
p^10y, p^8y², p^6y³, p^4y⁴, p²y⁵, y⁶,
p^7y³, p^5y⁴, p³y⁵, py⁶,
p^4y⁵, p²y⁶, y⁷, py⁷.                              (3)
```

Adjoin only the designated high channel `p³y⁶`, of valuation 15 and
weight 24. This is not a classification of every weight-24 monomial.

## 2. Five exact necessary and sufficient relations

Let `c_(a,b)` denote the coefficient of `p^a y^b` in the allowed source.
For readability put

```text
M=383445497,
a=c_(8,2), b=c_(4,4), c=c_(0,6),
d=c_(5,4), e=c_(1,6), v=c_(2,6), h=c_(3,6).
```

**Theorem.** The full system (2) is solvable if and only if

```text
M a+85731129 v+38329362 h=0,
2M b-200039301 v-89435178 h=0,
c=0,
2M d-225793440 v-192065985 h=0,
4M e-869307264 v-261093861 h=0.                       (4)
```

All eight odd-`y` coefficients in (3) are free. The two even coefficients
`v,h` are also free, and (4) determines the other five even coefficients.
For each admissible source, the complete perturbation fibre in (2) is
affine ten-dimensional. These statements hold over every characteristic-zero
field, including all specializations of the actual boundary parameter.

The primary reconstruction starts with all 155 coordinate coefficients
of `delta A,delta C`. Its tangent matrix has rank 145 and its source
quotient has rank five. The independent reconstruction first solves each
bracket row by its complete tangent parameterization. It retains 70
tangent parameters, all 58 remaining source coefficients and the 91
depth constraints obtained from literal source-generator matrices. Its
tangent matrix has rank 60 and its source quotient has rank five. Both
produce exactly (4), not merely the same ranks. This proves necessity,
sufficiency, and the fibre dimensions by rational linear elimination.

## 3. Complete signed replacements and sharp minimum

Define

```text
L0=(38329362/M)p^8y²-(44717589/M)p^4y⁴
   -(192065985/(2M))p^5y⁴-(261093861/(4M))py⁶,

K12=p²y⁶-(85731129/M)p^8y²+(200039301/(2M))p^4y⁴
    +(112896720/M)p^5y⁴+(217326816/M)py⁶.             (5)
```

The sign convention is explicit: to replace a positive source term
`p³y⁶` by a lower source, use the response difference `R=L-p³y⁶`,
so that `h=-1` in (4). Every weight-at-most-22 replacement is exactly

```text
L=L0+v K12+alpha p^6y³+beta p²y⁵+gamma p³y⁵+delta y⁷,
v,alpha,beta,gamma,delta in k.                        (6)
```

Thus the complete replacement family is `A^5`, with four odd directions
and one even direction. At weight at most 23, the four additional neutral
directions are `p^10y,p^7y³,p^4y⁵,py⁷`, giving `A^9` for one fixed high
coefficient. Setting `v=52578/117601` and the new odd directions zero
recovers the earlier valuation-thirteen replacement. The new family
therefore extends the inherited one rather than contradicting it.

No nonzero high response can be replaced by weight at most 21 in this
universe. Such a source has `a=d=v=0` because those three monomials have
weight 22; the first relation in (4) then forces `h=0`. Equivalently the
source ranks before and after adjoining the high column are

| Weight cap | Low source rank | With the high column |
| --- | --- | --- |
| 18 or 19 | 1 | 2 |
| 20 or 21 | 3 | 4 |
| 22 or 23 | 5 | 5 |

The minimum 22 is sharp when the high coefficient is nonzero. A zero
coefficient requires no replacement and creates no rank singularity.

## 4. Completeness, inheritance, and the finite nonlinear range

At a fixed valuation, distinct `p^a y^b` have distinct leading `x^b`
coefficients. Hence terms of lower valuation cannot cancel to create a
missing allowed source. Enumerating `12<=a+2b<=15` and `2a+3b<=23`
gives exactly (3). Terms of valuation above 15 do not enter (2), and none
of weight at most 23 have that valuation. The source enumeration is
therefore exhaustive in the stated filtered polynomial space.

For a perturbation beginning in row eleven, quadratic source terms
start in row 22 and quadratic Jacobian terms in row 21. Thus (2) is
exactly linear in the perturbation, including for arbitrary scalar
parameters, and depends only on background rows zero through four.
This justifies the affine conclusion for actual finite jets. It does not
justify an affine interpretation at those later nonlinear rows.

The row-`n` bracket operator is

```text
n[(x/2) delta C_n+(3/8)(x²+2) delta A_n].            (7)
```

It is surjective on the capped polynomials, with kernel
`theta_n(A0',C0')`, `deg theta_n<=n`. For example if its required value
divided by `n` is `F(x)`, choose

```text
delta A_part=4F(0)/3,
delta C_part=2[F-(3/8)(x²+2)delta A_part]/x.
```

The numerator vanishes at zero, and the degree caps hold. Adding the
kernel gives every solution. This gives `12+13+14+15+16=70` tangent
coordinates. The independent depth matrices enumerate all
`x^a u^b p^c y^e` with `a+b<=d` and initial row at most 15, and compute
their actual left-nullspaces. No proposed incomplete bank or deduplicated
response replaces those constraints.

For the newly needed fifth background row, THM-4308 on the inherited
boundary has

```text
Phi=eta=0, Delta=896/15, Theta=512/75, K=-32/5,
G5=-731648/2025+(6144/25)x²-(1952/45)x⁴.              (8)
```

The independent audit substitutes `A4,C4` from (1) into the entire row-four
bracket and the predicted row-five source and obtains (8). The homogeneous
change between any two such solutions is `theta_4(A0',C0')`, and its next
source response is the injective `m=5` Student operator from THM-4308.
Thus (1) is uniquely determined. Every additional source face in the
THM-4426/4438 construction begins after row five, so the row is constant
on the entire actual boundary. This supplies the required inheritance
step rather than inferring it from a single rational boundary point.

## 5. Consumer on the actual boundary `G_m`

Use THM-4438's entire inherited partial source
`H_pre(s)+j(s)p³y⁶`, for `s in G_m`; `j(s)` is its normalized high
coefficient. Fix its actual unknown rows through ten. Applying the
complete response classification gives exactly the sources

```text
H_pre(s)+j(s)L0+v K12
  +alpha p^6y³+beta p²y⁵+gamma p³y⁵+delta y⁷,         (9)
```

with all five displayed new parameters arbitrary. Every such source has
an affine ten-dimensional row-fifteen perturbation fibre relative to that
fixed prefix. Conversely every replacement within the declared source
valuation and weight bounds and with that fixed unknown prefix has form
(9). This is an exact finite lift on `G_m x A^5`.

The old `z,h` coefficients of valuations nine and ten are unchanged by
these responses and recover the actual boundary parameter. The four odd
coordinates and the `p²y⁶` coordinate then recover the five new parameters,
so the product describes actual source parameters. All linear pivots are
nonzero rational constants; no extra boundary denominator or exceptional
rank locus occurs. At a zero of `j(s)`, (9) still describes the same
five-dimensional neutral family. For every rational `s!=0`, THM-4438's
irreducible-quartic hostile gives `j(s)!=0`, so the sharp weight threshold
applies at every such rational point.

The first failed shortcut beyond this result would be to keep the same
operator while allowing source valuation eleven: a new unknown row and
a new background row would enter. Likewise a later-row conclusion would
require the nonlinear terms and all subsequent depth compatibilities.
Neither extension is asserted here.

## 6. Reproduction and independent verification

The [primary source](../../04-computation/planar_jc_long_20260906_memory_valuation12.py)
and [output](planar_jc_long_20260906_memory_valuation12.out) use the literal
full coefficient system. The independent
[audit](planar_jc_long_20260906_memory_valuation12_audit.md),
[source](../../04-computation/planar_jc_long_20260906_memory_valuation12_audit.py),
and [output](planar_jc_long_20260906_memory_valuation12_audit.out) use the
different tangent-row and literal depth-generator route described above.
Neither engine imports the other or a later boundary producer.

```bash
python3 -B 04-computation/planar_jc_long_20260906_memory_valuation12.py
python3 -B 04-computation/planar_jc_long_20260906_memory_valuation12_audit.py
python3 -B -O 04-computation/planar_jc_long_20260906_memory_valuation12_audit.py
```

The audit checks the entire five-parameter affine family by lifting six
source columns (one affine base and five directions) and substituting
them into every original equation. It separately constructs all ten
terminal kernel directions. The primary has1,646 live checks and the audit
5,146; both normal/optimized output pairs agree. Final byte hashes are in
the session manifest. The proved statement remains explicitly finite-row and
uses projected, rather than full polynomial, depth.
