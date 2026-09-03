---
id: THM-4380
title: "Source-normal weight-twelve row-twelve extinction"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4308/4315 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; FULL SECOND IMPLEMENTATION PENDING. In the complete fixed source-normal
  residual-weight-at-most-twelve family, the row-nine bracket hypersurface
  has automatic joint projected P_2/P_3 depth. Its Phi=0 stratum has two
  reduced row-ten points and dies at row eleven. Its Phi!=0 stratum has a
  row-ten ratio-presented carrier, exactly fourteen reduced row-eleven points with
  affine-eight terminal fibres, and dies at row twelve by a literal mod-29
  Bezout certificate. This is a finite-chart extinction theorem only. Chart
  or seam entry, all weights, an all-row lift, a Keller pair, JC(2), and
  DC(2) remain OPEN.
source: root + source_seam / JC2 continuation session, 2026-09-03
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
related:
  - THM-4357-source-normal-row-eight-endpoint-pullback-stratification
  - THM-4358-source-normal-s4339-row-ten-delayed-depth-extinction
  - THM-4360-source-normal-zeta-zero-row-ten-delayed-depth-extinction
  - THM-4361-source-normal-beta-zero-row-ten-joint-depth-extinction
  - THM-4364-sharp-binomial-diagonal-annihilator-hierarchy
  - THM-4366-source-normal-u-zero-row-eleven-hierarchy-selected-extinction
  - THM-4376-source-normal-u-zero-row-eleven-depth-hierarchy-completeness-and-bracket-blindness
mistake_firewall:
  - MISTAKE-354
  - MISTAKE-522
  - MISTAKE-540
primary_script: 04-computation/jc2_source_normal_weight12_row12_extinction_thm4380.py
primary_output: 05-knowledge/results/jc2_source_normal_weight12_row12_extinction_thm4380.out
primary_script_sha256: 15cd129452c3da033fa59985dda435077fe4526dbdf9d9df6feae32a8cb0ac6a
primary_output_sha256: 5b37b424137dba93fe7e7c4c621cf96d0068a9247a8bc9a3ff962bbf74279429
hash_basis: raw LF bytes
audit: >
  PASS WITH RATIONALITY-WORDING AND RELATED-SLUG REPAIRS; FULL SECOND
  IMPLEMENTATION PENDING. The primary has 314 exact checks. An independent
  theorem/code/branch-algebra review and triple replay found no proof failure;
  it is an audit, not an independent full implementation. The
  primary rebuilds THM-4308's capped source rows, bracket selectors, and
  P_2/P_3 projection matrices; checks at every new row that degrees 0 through
  the row exhaust the bracket polynomial; audits both Phi strata, terminal
  ranks and pivots, all denominator localizations, the fourteen-point radical
  carrier, all 45 row-eleven depth residuals, the row-twelve obstruction, and
  finite-field positive controls. Normal, optimized, and hash-seeded replays
  agree byte-for-byte with the frozen LF output. Scratch broadcast commit
  84590342a2 independently corroborates the strict Phi*U!=0 branch only; it
  is evidence, not a checked-in proof dependency or a full-family audit.
---

# THM-4380 -- Source-normal weight-twelve row-twelve extinction

**PROVED FINITE-ROW RELATIVE TO THM-4308/4315 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; FULL SECOND IMPLEMENTATION PENDING. THIS CLOSES THE FIXED THM-4308
SOURCE-NORMAL RESIDUAL-WEIGHT-AT-MOST-TWELVE FAMILY THROUGH ROW TWELVE.
CHART OR SEAM ENTRY, ALL WEIGHTS, AN ALL-ROW LIFT, POLYNOMIAL TERMINATION,
A KELLER PAIR, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Exact universe and theorem

Work over an algebraically closed field `k` of characteristic zero in the
fixed source-normal chart and residual-weight-at-most-twelve universe of
THM-4308. Retain its notation and abbreviate

```text
P=Phi,       e=eta,       X=xi_10,       a=alpha_11,
b=beta_11.                                                    (1)
```

THM-4315 gives the complete row-nine bracket equation

```text
E_9=613527750P^2-511211250Pa-3154140000Pe
    -255605625e^2+6736896000X-46483785515008=0.                (2)
```

The conclusion is

```text
there is no point of the fixed THM-4308 source family that satisfies
the bracket equations and the required projected P_2/P_3 depth equations
through row twelve.                                               (3)
```

This is a direct full-family calculation inside the declared finite chart.
It neither assumes `U=0` nor uses the proved counterfactual restricted-depth
completeness theorem THM-4376, which studies an already bracket-dead `U=0`
row-ten carrier. THM-4376 is related provenance only and is not a dependency
of `(3)`.

The exact bracket audit is coefficientwise and exhaustive. After every prior
row has been selected, the inherited caps give

```text
deg_x(G_n-predicted_G_n)<=n.                                  (4)
```

The primary certificate recomputes `(4)` before each row-ten, row-eleven, and
row-twelve selection in both strata. It then compares every coefficient in
degrees `0,...,n`. Hence the reported bracket residuals are the whole capped
row obstruction, not a chosen projection of a larger polynomial.

## 2. The entire row-nine hypersurface has automatic joint depth

On `(2)`, including `P=0`, the exact row-nine projection data are

```text
pi_9(P_2): 75 rows, 160 columns, rank 59, left nullity 16;
pi_9(P_3): 85 rows, 251 columns, rank 73, left nullity 12.      (5)
```

The terminal coefficient ranks for `P_2`, `P_3`, and the combined system are

```text
3, 2, 3,                                                     (6)
```

with joint pivots in coordinates `7,8,9`. All 28 projected-depth residuals
vanish after the selected terminal solution. Thus depth imposes no source
equation beyond `(2)` and leaves an affine-seven terminal fibre.

This automaticity is important: neither later extinction is a disguised
row-nine emptiness statement.

## 3. The `P=0` stratum: two row-ten points die at row eleven

Set `P=0`. Equation `(2)` gives the complete row-nine graph

```text
X=(255605625e^2+46483785515008)/6736896000.                    (7)
```

The row-ten bracket selector has rank seven. Its only possible residual
positions are 6 and 8, and they are equivalent to

```text
F(e)=18612736875e^2-4820239249178624=0,

a=-(854861882349375e^2+27743789318253707264)
    /(278844730740000e).                                      (8)
```

The quadratic `F` is squarefree and has nonzero constant term. Therefore it
has exactly two geometric roots and neither has `e=0`; the division by `e`
in `(8)` is lawful at both points.

The row-ten depth universes and ranks are

```text
pi_10(P_2): 88 rows, 193 columns, rank 68, left nullity 20;
pi_10(P_3): 99 rows, 304 columns, rank 83, left nullity 16,

terminal ranks P_2/P_3/joint = 3/3/3,
terminal pivots             = (8,9,10).                     (9)
```

The two standalone selected solutions agree in coordinates 9 and 10. Their
coordinate-8 compatibility is exactly

```text
b=-6e/5.                                                       (10)
```

Consequently `(7)--(10)` are exactly two reduced row-ten
bracket-and-joint-depth source points, each with an affine-eight terminal
fibre. Their inherited response coordinates are the fixed nonzero values

```text
U=-54752794624/8110125,
W= 6445195264/180225,
Z=-307951616/11125.                                           (11)
```

The row-eleven bracket selector has rank eight. Modulo `F`, every selected
residual vanishes except position 8, where the remainder is

```text
12353759608682044232849562614853632
------------------------------------------------ !=0.          (12)
137875982427855174564984375
```

Thus both genuine `P=0` row-ten points die at row eleven.

## 4. The `P!=0` stratum: the row-ten ratio-presented carrier

Assume `P!=0`. Equation `(2)` first solves

```text
a=(613527750P^2-3154140000Pe-255605625e^2
   +6736896000X-46483785515008)/(511211250P).                 (13)
```

The rank-seven row-ten bracket selector then gives

```text
X=-(13365000P^2+15035625Pe-964604821504)/57672000,

a=-(69009144750P^2+357574500000Pe+18612736875e^2
    -4820239249178624)/(37225473750P),

b=(11296175760000P^3+49737870463500P^2e+6793915500000Pe^2
   -631918028977864704P+353642000625e^3
   -91584545734393856e)/(670058527500P^2).                   (14)
```

The row-ten projection universes have the same sizes and ambient ranks as
`(9)`. Both terminal systems have rank three with pivots `(8,9,10)`, and the
combined coefficient rank is three. Coordinates 9 and 10 agree identically;
the coordinate-8 difference, `P_2` minus `P_3`, is

```text
-D(P,e)/(502543895625P^2),                                   (15)
```

where

```text
D(P,e)=7231154026500P^3+50541940696500P^2e
      +6793915500000Pe^2-631918028977864704P
      +353642000625e^3-91584545734393856e.                   (16)
```

Hence row-ten joint depth is exactly the ratio-presented affine curve `D=0` in the
open chart `P!=0`, with an affine-eight terminal fibre over each source
point.

## 5. Row eleven leaves exactly fourteen reduced source points

The row-eleven bracket selector has rank eight. Its selected residual support
is `{8,9}` and the two residuals are

```text
R_9=-2D/(1842660950625P^2),
R_8= E/(1971983693758919400000000P^2),                       (17)
```

where

```text
E(P,e)=
 35898624648545625000000P^6
+80771905459227656250000P^5e
+45434196820815556640625P^4e^2
-5071037898252514202802000000P^4
-5410571375319483376728000000P^3e
+106595566955131983930000000P^2e^2
+179564695363704166433073974476800P^2
+16899675818606082900000000Pe^3
-2249959222183587936424427520000Pe
+851303752055234564062500e^4
-268196544113786432750223360000e^2
+12360863815079359770761606339756032.                         (18)
```

Introduce the invariant coordinates

```text
r=e/P,                     Y=P^2.                            (19)
```

Then `(16)` divided by `P` is

```text
A(r)Y-B(r)=0,                                                  (20)

A(r)=353642000625r^3+6793915500000r^2
    +50541940696500r+7231154026500,

B(r)=91584545734393856r+631918028977864704.                  (21)
```

Here `gcd(A,B)=1`. Thus a solution of `(20)` can never have `A=0`; every
solution has the unique affine value `Y=B/A`. Conversely this substitution
loses no point of the `P!=0` carrier.

Substitute `Y=B/A` into the even invariant form of `(18)`, multiply by
`A^3`, and take the primitive numerator. The result is the septimic

```text
K(r)=
 21252176198679866250754006276839556755825r^7
+799311827675117522149997435401077574131600r^6
+14384863896403857958176347858723924433398460r^5
+142730433788981669223548142320603830110956220r^4
+786944231209420107856657052701244375708027892r^3
+1852564642916723756803328543705267257790149632r^2
-147733098192443646925107791876239203619548432r
+3714269896529642422852685702214695613036368.                (22)
```

Exact Euclidean calculations give

```text
deg K=7,       gcd(K,K')=1,       gcd(K,A)=gcd(K,B)=1.        (23)
```

There is also no hidden `U=0` point on this carrier. The inherited response
on `(14)` is

```text
U=(32805000P^2+36905625Pe-1752089427968)/259524000.
```

After `(19)--(21)`, clearing the single factor `A` gives, up to the nonzero
content `531505152000`,

```text
Ucar(r)=-1165769340615r^3-16036650988830r^2
        -117079277403336r+15165313804484,                    (23a)

gcd(K,Ucar)=1.                                                (23b)
```

Thus every surviving point has `P*U!=0`; this conclusion is derived here and
is not assumed as an entry gate.

The leading and constant coefficients of `K` are nonzero. Thus all seven
roots of `K` are distinct and affine; there is no lost projective root at
infinity. At each root, `A` and `B` are nonzero, so `Y=B/A!=0`. There are
exactly two distinct choices of `P` with `P^2=Y`, and then `e=rP` is unique.
It follows that `(D,E)` has exactly fourteen geometric points on `P!=0`.

This also proves radicality rather than merely a length count: after
localizing at `P`, the changes `(P,e)<->(P,r)` and `(20)` identify the carrier
with the squarefree `K` algebra followed by `P^2-Y`; the latter is separable
because `Y!=0` in characteristic zero. Hence the fourteen points are reduced.

## 6. Row-eleven depth is automatic on the fourteen points

After selecting row eleven, the exact depth data are

```text
pi_11(P_2): 102 rows, 228 columns, rank 77, left nullity 25;
pi_11(P_3): 114 rows, 361 columns, rank 94, left nullity 20.   (24)
```

The terminal coefficient ranks are

```text
P_2/P_3/joint = 4/3/4,
P_2 pivots     = (8,9,10,11),
P_3 pivots     = (9,10,11).                                  (25)
```

A Groebner reduction modulo `(D,E)` proves that the three shared selected
coordinates agree and that all `25+20=45` projected-depth residuals vanish.
Every rational terminal formula used in this reduction is audited before
numerator reduction: its denominator is a single monomial in `P` and is
independent of `e`. Thus the computation localizes only at the already
declared `P!=0`; it does not silently invert `A`, `B`, `K'`, or another
possible carrier divisor. Each of the fourteen reduced source points
therefore has an affine-eight row-eleven terminal fibre.

## 7. The row-twelve bracket kills all fourteen points

On the row-ten graph `(14)`, the literal weight-twelve source row is

```text
G_12=-13x^12(3918375P^2+1366875Pe-3318235136)/32440500.       (26)
```

The row-twelve selector has rank eight. Reduce the numerators of all selected
residuals modulo `(D,E)`. Only positions 8 and 9 can remain, and position 9
is exactly

```text
-1650P C(P,e),                                                (27)

C(P,e)=1245684290400P^2-592839494100Pe
      -18612736875e^2+4820239249178624.                       (28)
```

On the ratio graph `Y=B/A`, equation `(28)` satisfies

```text
A(r)C(P,e)=7228470067200J(r),                                (29)

J(r)=-4607940726340893r^2-2340230825466590r
     +113720641620096958.                                    (30)
```

The extinction certificate is literal. Reduction modulo 29 preserves the
degrees and leading coefficients of both `K` and `J` and gives

```text
Kbar=8r^7+9r^6-10r^5+12r^4+2r^3+8r^2-12r+3,
Jbar=-11r^2+8r+5,                                            (31)

(11r+11)Kbar
 +(8r^6+7r^5+13r^4+12r^3-3r^2+8r+11)Jbar = 1    in F_29[r]. (32)
```

Because the primitive integer polynomials retain their degrees at 29,
`(32)` certifies `gcd_Q(K,J)=1`. At every one of the fourteen points,
`P!=0` and `A!=0`; if `(27)` vanished then `(29)` would force `J=0`,
contradicting `K=0` and coprimality. Thus every row-eleven point dies at row
twelve.

Combining the strata proves `(3)`: points off `E_9=0` die at row nine, the
two `P=0` points die at row eleven, and all fourteen `P!=0` points die at row
twelve.

## 8. Nonvacuity controls, replay, and scope firewall

The exact finite-field control over `F_31` is

```text
r=5, Y=8, P=15, e=13, X=20, a=16, b=1,
U=26, W=22, Z=5,
D=E=K=0,                    J=8, C=14.                       (33)
```

It is a positive control for the complete staged carrier through row eleven
and a hostile to premature extinction; its nonzero `J,C` detect the intended
row-twelve cut.

The primary script performs 314 exact assertions. Its normal, `python -O`,
and `PYTHONHASHSEED=17` streams are byte-identical to the frozen output. Raw
LF hashes are

```text
primary script:
15cd129452c3da033fa59985dda435077fe4526dbdf9d9df6feae32a8cb0ac6a

primary output:
5b37b424137dba93fe7e7c4c621cf96d0068a9247a8bc9a3ff962bbf74279429
```

These canonical hashes supersede the pre-check-in scratch pair recorded in
the first promoted body; see MISTAKE-541. The theorem, output, and
frontmatter values were unchanged by that documentation-only repair.

Replay from the repository root with

```text
python3 -B 04-computation/jc2_source_normal_weight12_row12_extinction_thm4380.py
python3 -B -O 04-computation/jc2_source_normal_weight12_row12_extinction_thm4380.py
PYTHONHASHSEED=17 python3 -B 04-computation/jc2_source_normal_weight12_row12_extinction_thm4380.py
```

Scratch broadcast commit `84590342a2` independently corroborates the
`P*U!=0` portion: the nonempty row-ten curve, squarefree degree-seven
eliminant, avoidance of `A*B*U`, fourteen row-eleven points, all 45 depth
residuals with affine-eight fibres, and the mod-29 row-twelve coprimality
cut. No checked-in artifact accompanied that broadcast, and it did not
rebuild the `P=0` stratum or the full inherited source construction. It is
therefore corroborating evidence, not a proof dependency or a full
clean-room audit. An independent full-family implementation remains pending;
this is stated in the status rather than hidden.

Finally, `(3)` is strictly conditional on entry into the fixed THM-4308
source-normal residual-weight-at-most-twelve finite chart. It proves no chart
or seam entry, no exclusion of higher residual weights, no all-row statement,
no polynomial termination theorem, no Keller-pair theorem, and no conclusion
for `JC(2)` or `DC(2)`. All of those scopes remain **OPEN**.
