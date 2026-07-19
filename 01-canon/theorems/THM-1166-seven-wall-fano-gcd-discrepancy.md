---
id: THM-1166
title: Seven-wall quadratic/Fano gcd discrepancy at the actual one-fourteenth radius
status: PROVED analytic cover implications and FINITE-EXACT three-speed overlap floor.  Every seven-comb packet has global uncovered mass at least 1/12; a covered lower-LRC needle forces common dilation G/m<=77/12, every labelled Fano plane obeys sum m/G_line>=32/231, and every forest obeys an exact adaptive edge-gcd budget.  This does not prove a uniform local overlap constant, crown collapse, or LRC(14)
source: codex-2026-07-18-S75
depends_on: [LEM-042, LEM-043, THM-1153]
related: [THM-856, THM-1156, HYP-7678]
script: 04-computation/lrc14_seven_wall_fano_gcd_codex_20260718.py
output: 05-knowledge/results/lrc14_seven_wall_fano_gcd_codex_20260718.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCSevenWallFanoGCD.lean
referee:
  - 04-computation/lrc14_fano_gcd_discrepancy_referee_codex_20260718.py
  - 05-knowledge/results/lrc14_fano_gcd_discrepancy_referee_codex_20260718.out
---

# THM-1166 -- seven-wall Fano/gcd discrepancy

Work at the actual LRC(14) radius.  For a positive integer `s`, put

```text
D_s={t in R/Z: ||st||<1/14}.
```

Let `s1,...,s7` be distinct positive speeds, let `I` be an interval of length
`L`, and suppose their union covers `I`.  In the LRC application, the
protected-interval construction and the already-known lower-case LRC(7)
input supply the separate hypothesis

```text
L>=1/(7m),                 m=max(retained core).         (1)
```

THM-1153's first-order coefficient is zero here.  The point of this theorem is
that the exact arithmetic pair mass is nevertheless large enough to impose
scale constraints once periods and Fano incidence are retained.

## 1. Exact pair and three-speed floors

Write

```text
rho(a,b)=measure(D_a intersect D_b).
```

After dividing by `gcd(a,b)`, take coprime `a<=b`, and put
`F(r)=r(14-r)` for the representative `0<=r<14`.  The folded form of
LEM-042's trapezoid identity is

```text
rho(a,b)=
 [4ab+F((a+b) mod 14)-F((b-a) mod 14)]/(196ab).         (2)
```

Since the folded difference is at least `-49`,

```text
rho(a,b)>=1/49-1/(4ab).                                 (3)
```

For `ab>=27`, the right side is greater than `1/91`.  The complete 36-row
coprime bank `ab<=26` has unique minimum

```text
rho(1,13)=1/91.                                         (4)
```

Thus (4) is the universal pair floor, with equality only at reduced ratio
`1:13`.

There is a stronger collective fact.

> **Three-speed lemma.**  For any three distinct speeds,
>
> ```text
> rho12+rho13+rho23>=1/24.                               (5)
> ```

**Proof.**  Suppose the sum were below `1/24`.  By (4), every individual edge
would then be below

```text
tau=1/24-2/91=43/2184.                                  (6)
```

LEM-043's defect bound

```text
rho(a,b)>=1/49-1/(7b')
```

for reduced `a'<b'` forces `b'<=198`, because

```text
1/49-tau=11/15288,       b'<2184/11<199.                (7)
```

Exact evaluation leaves 106 oriented reduced ratio types below `tau`.  For
ordered speeds `x<y<z`, all three rational ratios

```text
r=y/x,       s=z/x,       s/r=z/y
```

would have to lie in that list.  Its exact quotient closure has 173 triples.
Their minimum pair sum is

```text
11/252=1/24+1/504,                                      (8)
```

attained only by the two ratio words `(9/4,27,12)` and `(12,27,9/4)`.
This contradicts the assumption and proves (5).  The 106-type and 173-triple
enumerations are frozen, rational, and independently recompute (2). ∎

Summing (5) over the 35 triples of seven vertices counts every edge five
times.  Therefore

```text
R:=sum_(i<j)rho(si,sj)>=7/24.                            (9)
```

## 2. The optimal quadratic certificate and the j=4 chamber

At a point `t`, let `C(t)` be the number of active danger combs and put

```text
Q(C)=C-(2/7)binom(C,2)=C(8-C)/7.                        (10)
```

For `C=0,...,7`, its values are

```text
0, 1, 12/7, 15/7, 16/7, 15/7, 12/7, 1.                (11)
```

Thus `1_{C>=1}<=Q(C)`.  Integration gives

```text
measure(union_i D_si)<=1-(2/7)R<=11/12,                 (12)
measure(uncovered)>=2R/7>=1/12.                         (13)
```

The coefficient `2/7` is optimal within the normalized family
`C-alpha*binom(C,2)`: imposing the majorant at `C=7` forces
`7-21alpha>=1`, hence `alpha<=2/7`.

This also has a multiplicity-debt reading.  If `mu_k=measure{C=k}`, the mean
`integral C=1` gives

```text
mu_0=sum_(k>=2)(k-1)mu_k,
R=sum_(k>=2)binom(k,2)mu_k<=(7/2)mu_0.                  (14)
```

Equality in the last comparison can put nonprivate debt only at depth seven.
On a covered chamber,

```text
Q(C)-1=(C-1)(7-C)/7.                                   (15)
```

The slack vanishes at private depth one and full depth seven and peaks exactly
at `C=4`.  Hence the open `j=4` flood tail is not an unrelated side problem:
it is the central chamber of the unique symmetric quadratic certificate whose
two boundary depths are tight.

## 3. Common-dilate branch: a genuine positive seventh crown

Suppose a common integer `G` divides all seven deleted speeds.  Their union is
`1/G`-periodic.  By (13), every period contains uncovered mass at least
`1/(12G)`.  Hence a connected interval contained in the union has length at
most the total covered mass of one period:

```text
L<=(1-2R/7)/G<=11/(12G).                                (16)
```

Combining (16) with the lower-LRC needle (1) yields the sharper exact form

```text
G/m<=7-2R<=77/12.                                       (17)
```

This crosses the zero-density wall on the common-dilate branch.  It is a
local continuous Mirsky-Newman statement: global pair mass creates a hole in
every common period, and a covered protected needle cannot cross that hole.

## 4. Fano decomposition for mixed triple periods

Fix any labelled Fano plane on the seven comb vertices.  Its seven lines are
triples, every vertex lies on three lines, and every pair lies on exactly one
line.  For a line `ell`, define

```text
c_ell(t)=sum_(i in ell)1_{D_si}(t),
f_ell(t)=c_ell/3-(2/7)binom(c_ell,2),
R_ell=sum_({i,j} subset ell)rho(si,sj),
G_ell=gcd(si:i in ell).                                 (18)
```

The four values of `f_ell` at `c_ell=0,1,2,3` are

```text
0, 1/3, 8/21, 1/7.                                     (19)
```

Fano incidence gives the exact pointwise identity

```text
Q(C)=sum_(ell) f_ell.                                   (20)
```

The function `f_ell` has period `1/G_ell`, maximum `8/21`, and mean

```text
nu_ell=(1-2R_ell)/7.                                    (21)
```

Use the sharp periodic-remainder lemma: if `0<=f<=M` is `P`-periodic with
mean `nu`, then every interval `J` satisfies

```text
integral_J f<=nu|J|+nu(1-nu/M)P.                        (22)
```

Indeed, remove all full periods.  On a remainder of length `r`, its mass is
at most `min(Mr,nu P)`; relative to `nu r`, the maximum excess occurs at
`r=nu P/M` and equals `nu(1-nu/M)P`.

Put

```text
e_ell=nu_ell(1-21nu_ell/8).                             (23)
```

Since the seven combs cover `I`, (10), (20), and (22) give the exact necessary
metric-Fano inequality

```text
(2LR)/7<=sum_(ell) e_ell/G_ell.                         (24)
```

Every line is a three-speed packet, so (5) gives `R_ell>=1/24`.  Consequently
`0<=nu_ell<=11/84<4/21`, and the quadratic in (23) is increasing on this
range.  Hence

```text
e_ell<=11/128.                                          (25)
```

Insert (1), (9), and (25) into (24): for **every** labelled Fano plane,

```text
sum_(ell) m/G_ell>=32/231.                              (26)
```

In particular every labelled plane has a line satisfying

```text
G_ell/m<=1617/32.                                       (27)
```

Equivalently, if all seven triple gcds of some plane are at least `G`, then
`G/m<=1617/32`.  This gives the Fano/`chi_7` probe a concrete metric consumer:
`chi_7` may address or compare the 30 labelled planes, but it cannot discard
their numerical line periods and pair masses.

## 5. Adaptive forest certificate

The same periodic-remainder argument yields, on any interval of length `L`,

```text
measure(I intersect D_s)<=L/7+6/(49s),
measure(I intersect D_si intersect D_sj)
 >=L rho_ij-rho_ij(1-rho_ij)/gcd(si,sj).                (28)
```

For every forest `F` on the seven comb vertices, the pointwise forest-Hunter
inequality and the cover of `I` therefore force

```text
sum_({i,j} in F)
 [L rho_ij-rho_ij(1-rho_ij)/gcd(si,sj)]
 <=(6/49)sum_i 1/si.                                    (29)
```

This is an adaptive exact edge-gcd certificate.  Negative edge lower bounds
may be omitted; a proof search should maximize the positive forest weight,
not preselect a speed-order path.

## 6. Why global pair credits still do not localize for free

There is an exact protected-needle guardrail.  Let the retained core contain
six distinct odd speeds at most an odd `m>=11`, and put

```text
I=[1/2-1/(14m),1/2+1/(14m)].                            (30)
```

This interval is safe for every core speed.  Choose an arbitrarily large `a`
and the six further deleted speeds

```text
a+1, 2a+1, 2a+3, 3a+1, 4a+1, 5a+1.                    (31)
```

If `b=qa+c` and both `a,b` were dangerous at `t`, then

```text
||ct||<(q+1)/14.                                        (32)
```

On (30), odd `c` instead gives

```text
||ct||>=1/2-c/(14m).
```

For `(q,c)=(1,1),(2,1),(2,3),(3,1),(4,1),(5,1)`, this contradicts (32).
Thus the six pair intersections with the central speed `a` are all empty on
`I`; they form an empty spanning star while all seven deleted speeds remain in
a bounded ratio cone and `a/m` is arbitrarily large.

This construction does **not** say that the seven combs cover `I`, and so it
does not refute (24), (26), or a cover-conditioned adaptive Fano theorem.  It
does refute any cover-free rule that inherits a fixed global Hunter tree on an
arbitrary protected needle.

## 7. Carrier and verification audit

Runner order, Fano line-mass `R_ell` order, and line-period order all give
transitive tournaments after ties are lexicographically resolved: score
histogram `0,...,6`, no directed cycles, seven singleton SCCs, and one
Hamiltonian path.
They lose which edges share a Fano line, the values `R_ell,G_ell`, and the
endpoint discrepancy in (24).  Candidate vertices were challenged as runners,
pairs, Fano lines, slow gaps, wall events, and proof obligations.  The faithful
current carrier is

```text
labelled Fano incidence hypergraph
  + (R_ell,G_ell)
  + protected-needle endpoint/owner data.               (33)
```

The exact script checks (2) against an independent trapezoid evaluator on
6,360 pairs, the 36-row pair tail, all 106 strict-low ratios and 173 compatible
quotient triples, all 30 labelled Fano planes and all 128 activity patterns,
every displayed rational constant, and the empty-star inequalities.  Normal
and optimized runs are byte-identical.  Frozen hashes are

```text
source  34ad825245f7288bd2253db1472b875f6303c2e604d62473be05312586c0ee0b
output  b2781fdf746bf3a9584508b9f8b0ab4f73369c5da5138cfa77323bf4fa37c802
```

The separate referee script independently checks the sharp period-cell
discrepancy, the full pair-density range `1/91<=rho<=1/14`, all 128 Fano
activation masks, and the weaker density-only selected-line consequence

```text
m[H_ell+(13/28)sum_(pairs in ell)1/gcd]>=3/91.
```

Its source/output hashes are respectively
`51516e7156ada06424af7cf4b594d70a73fb3e13d87db5db1f8170063615a340`
and `da0cea50f3eb8998c36c0e6a1d280262d6a50ca61b23b0c8be84848d2504d3a4`.
The resulting common-dilate constant `3523/36` is a cross-check only and is
strictly dominated by (17).

What has closed is the common-dilate seventh crown and a mixed-period Fano/gcd
necessary law.  What remains open is a uniform way to turn (24) or (29) into a
contradiction for arbitrary seven deleted speeds, together with the broader
`r=5`, `r=6`, crown-collapse, n=12 equality, and LRC(14) obligations.
