# Product-Gamma Young-cover hazard order: an all-degree closure target

**Status: HYPOTHESIS / VERIFIED EXACT FINITE EVIDENCE, not a theorem and not
a proved dependency.**  The exact census below exposes a plausible all-degree
closure of the arbitrary-anchored width-three product-Gamma problem.  It also
rules out three simpler proofs.  Nothing in this note promotes THM-3110 or
settles the Gaussian Moment Conjecture by itself.

## Inheritance pass

The closest maintained result is the **PROVED** theorem
`THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction.md`.
Its exact 24/25-atom banks, chamber common-root deletion, dominant residual
alphabets, dual-Cauchy expansion, and degree-5 through degree-12 Schur/Newton
front face are the only inputs used here.  Its canonical near miss is the
first-order row-majorization transport: it leaves coefficient mass `9` and
`11`, so termwise Ferrers dominance cannot close the signed bank.

The least-used relevant sidecar is Schur-ratio monotonicity.  Belton--Guillot--
Khare--Putinar, *Matrix positivity preservers in fixed dimension II*
([arXiv:2310.18020](https://arxiv.org/abs/2310.18020)), Theorems 1.4 and 2.11,
prove coordinatewise monotonicity, including the required boundary form, for
the Schur ratio attached to coordinatewise ordered exponent tuples.  In the
partition notation needed here, the same conclusion says that
`s_mu(x)/s_nu(x)` is coordinatewise decreasing when `mu` is strictly contained
in `nu`.  Herbig--Herden--Kolehmainen--Seaton give an independent tableau-
injection proof of this formulation in Theorem A of *The partial derivative of
ratios of Schur polynomials*
([arXiv:2504.19466](https://arxiv.org/abs/2504.19466)).  These results turn
THM-3110's coordinatewise root domination `S<=Q_i` for every negative atom
into the exact one-atom hazard inequality used below.  They do **not** supply
the signed grouped inequality.

## The observed cover law

Let `Phi_i(f)=sum_R c_R f(S_R)` be the signed residual-alphabet functional.
The distinguished positive atoms are

```text
Q_1 = S_(3b,2a,2a,a),             c_1=1,
Q_2 = S_(3b,a+b,2a,a,a),          c_2=2.              (1)
```

Let

```text
B_1=a^4 b^2 (b-a)^3,              B_2=a^3 b^2 (b-a)^4. (2)
```

For every partition `lambda` with `s_lambda(Q_i)>0`, define

```text
R_i(lambda)=Phi_i(s_lambda)/(B_i c_i s_lambda(Q_i)).  (3)
```

The exact finite signal is the strict Young-cover law

```text
R_i(lambda + one box) > R_i(lambda).                  (4)
```

The first three cover ranks have a stronger **global exact certificate**.
For a cover `mu<nu` with `|mu|=n`, cross-multiply `(4)` and put

```text
D_(i;mu,nu)=Phi_i(s_nu)s_mu(Q_i)-Phi_i(s_mu)s_nu(Q_i). (4a)
```

The collision divisor and degree bounds in THM-3110 imply that `D/B_i` is
a chamber polynomial of degree at most `4n-7`.  Exact Gregory--Newton
interpolation, an entire excess-degree shell, and three off-grid points prove

```text
cover       degree     Newton slots     positive     structural zero
5 -> 6        13           7,980          7,980              0
6 -> 7        17          20,520         20,518              2
7 -> 8        21          45,540         45,534              6.       (4b)
```

There is no negative coefficient.  Thus `(4)` is already proved for every
integer support `0<a<b` through target degree eight; the larger census below
is evidence for the unproved continuation.  The coefficient degrees follow
the stable linear law `4n-7`, but no stable factor or recurrence beyond the
forced `B_i` has yet been extracted.

The eight displayed zero slots are harmless Newton constant terms, not zero
cover polynomials: all eight occur in the wide chamber on column or
near-column covers, and every one of those polynomials has other strictly
positive Newton coefficients.  They therefore record a boundary cancellation
pattern worth explaining in a future recurrence, while preserving strict
positivity of every certified cover.

The scout tested `(4)` on all `115` supports

```text
1<=a<=10,                  a<b<=min(3a+4,21),         (5)
```

both invariants, and every admissible cover from degrees `5..15` to degrees
`6..16`.  All `555,662` comparisons were strict.  The smallest additive slack
was

```text
125872476615673643 / 522634777460648197207357773696
```

at `I_2,(a,b)=(10,21),(5)->(6)`.  The smallest multiplicative ratio was

```text
53329143091 / 53305526863
```

at `I_1,(a,b)=(1,6)`, from `(2,1^13)` to `(3,1^13)`.

There are `49` added-box `(coarm,coleg,content)` profiles in the tested
universe, with coarm and coleg in `0..15` and content in `-15..15`.  The
smallest additive increment in every bank/chamber group is the first-row box
`(coarm,coleg,content)=(5,0,5)` on `(5)->(6)`.  The smallest multiplicative
ratio in the two wide groups is the first-row hook extension `(2,0,2)`.
The tight groups are different: their minima add a box near the bottom, at
`(1,6,-5)` for `I_1` and `(1,7,-6)` for `I_2`.  Thus a row/hook extremal
reduction is plausible for additive slack and wide ratios, but is false as a
literal all-chamber ratio statement.

The normalization by `B_i` is needed only to state the additive extremum; it
is independent of shape, so it changes neither the cover direction nor the
multiplicative ratio.

## A sharper survivor: minimize the response, not every cover increment

A second exact scan points to a narrower and potentially much cheaper
all-degree statement.  At each fixed degree `n`, compare every admissible
shape directly to the row shape:

```text
R_i(lambda) >= R_i((n)),             |lambda|=n.       (5a)
```

Across the same `115` supports, both banks, and degrees `5..16`, `(5a)` passes
all `206,454` admissible shape comparisons.  Exactly `2,760` are equalities,
one row shape for each support, bank, and degree; all `203,694` nonrow
comparisons are strict.  The smallest strict margin is

```text
87165569201173 / 17844075025881441684247232640
```

at `I_2,(a,b)=(10,21),n=5,lambda=(4,1)`.  The smallest row response in the
finite universe is `R_2((5))=17/259124718450048`, again at `(10,21)`.  If
`(5a)` holds generally, then
all-shape positivity reduces to positivity of the single complete-homogeneous
ray `Phi_i(h_n)`; no chain of Young covers is needed.

Three nearby strengthenings are false and therefore help type the target.

* The row cover `(n)->(n+1)` does **not** minimize every normalized cover
  increment.  There are `35` failures in the census.  The earliest is at
  source degree `14`, `I_1,(a,b)=(9,10)`, on the column cover
  `(1^14)->(1^15)`; the most negative tested gap is again a column cover, at
  degree `15`, `I_1,(a,b)=(5,6)`.
* Full dominance monotonicity is false earlier.  At `(a,b)=(1,2)`, degree
  eight, `(3,1,1,1,1,1)` dominates `(2,2,1,1,1,1)`, but `R_mu-R_lambda` is
  `-7/975` in `I_1` and `-461/100050` in `I_2`.
* Added-box content supplies no prefix proof.  Universal content intervals and
  even adjacent content-class minima cross already at `I_1,(a,b)=(1,2)`,
  degree five; the exact covers and fractions are in the transcript.

Nor is `(5a)` termwise in the residual atoms.  Already for `(a,b)=(1,2)`,
`lambda=(4,1)`, a positive non-dominant atom has

```text
s_lambda(S)/s_lambda(Q)-h_5(S)/h_5(Q)
 = -2532173/22477467                  in I_1,
 =   -83023/9923040                   in I_2.           (5b)
```

Thus the live statement is a grouped signed-circuit inequality

```text
Phi_i(s_lambda) h_n(Q_i)-Phi_i(h_n)s_lambda(Q_i) >= 0, (5c)
```

not Schur majorization, atomwise ratio order, content order, or a cover-ray
bound.  This makes a bank-preserving Jacobi--Trudi or tableau injection the
right proof object.

## Exact row generating function and positive limiting residue

The row ray has a compact exact rational representation.  Put

```text
W_L(t)=product_(r=0)^(L-1) (1-r t)^(-1),
C_1(t)=W_a(t)^3 W_max(2a,b)(t),
C_2(t)=W_a(t)^3 W_max(2a,b)(t)^2.                    (5d)
```

For signed directions `D_1={(a,+1),(0,-1)}` and
`D_2={(b,+1),(a,-1)}`, define

```text
G_pq = sum_(x in D_p,y in D_q) eps_x eps_y
       W_(x+y)/(W_x W_y),
T_pqr = sum_(x in D_p,y in D_q,z in D_r) eps_x eps_y eps_z
        W_(x+y+z)/(W_x W_y W_z),
I_1(t)=3 T_112 G_11 G_22-T_222 G_11^2-2 T_111 G_12 G_22,
I_2(t)=3 T_122 G_11 G_22-2 T_222 G_12 G_11-T_111 G_22^2. (5e)
```

Then the exact row numerator and dominant denominator series are

```text
F_1(t)=-W_a^5 W_b^3 I_1(t)/C_1(t),
F_2(t)=-W_a^5 W_b^4 I_2(t)/C_2(t),
H_Q1(t)=W_(3b) W_(2a)^2 W_a/C_1(t),
H_Q2(t)=W_(3b) W_(a+b) W_(2a) W_a^2/C_2(t),
R_i((n))=[t^n]F_i(t)/(B_i c_i [t^n]H_Qi(t)).        (5f)
```

This turns row positivity into a concrete rational-recurrence problem.
Moreover its leading pole factors completely.  Let `M=3b-1` and

```text
P(L)=product_(r=0)^(L-1)(1-r/M)^(-1)
    =M^L (M-L)!/M!,
x=P(a)^2/P(2a),
y=P(2a)P(b)/(P(a+b)P(a)).                            (5g)
```

The length-`a` quotient defining `P(L+a)/P(L)` increases with `L`, so
`0<x,y<1`.  Extracting the unique top pole from the three `I_1` atoms and four
`I_2` atoms which retain the row `3b` gives the exact limits

```text
lim_(n->infinity) R_1((n))=(1-x)^2/B_1,
lim_(n->infinity) R_2((n))=(1-x)(1-y)/B_2.           (5h)
```

Both are strictly positive.  This does not yet rule out a finite dip of the
row sequence, but it supplies an exact positive endpoint and suggests attacking
`[t^n]F_i` by total positivity or a recurrence barrier rather than another
Young-cover census.

The cheap recurrence scout finds no such dip: on all `230` support/bank pairs
and every degree `5..128`, all `28,520` row responses are positive, lie
strictly below `(5h)`, and all `28,290` consecutive steps are strictly
increasing.  The smallest unnormalized cross-product is `29,433,312`, at
`I_1,(a,b)=(1,2)`, from degree five to six.  This is exact finite evidence,
not an all-degree recurrence proof.

## Why a symbolic cover theorem would close every degree

For every partition `mu` of `5`, exact character interpolation gives

```text
Phi_1(s_mu)/B_1 = f^mu ((3a+5b)/2 + chi_mu/4),
Phi_2(s_mu)/B_2 = f^mu ((3a+4b)/2 + chi_mu/4),         (6)
```

where `f^mu` is the standard-tableau dimension and `chi_mu` is the sum of
cell contents.  Since `chi_mu>=-10` in degree five, both right sides are
strictly positive for `1<=a<b`.  This notation deliberately separates two
conventions.  THM-3110 defines
`kappa_mu=sum_i mu_i(mu_i-2i+1)=2 chi_mu` and labels the dual-Cauchy
coefficient `A_(j,lambda)=Phi_j(s_(lambda'))`; its formula is
`base_j f^lambda/2 (L_j-kappa_lambda/4)`.  Setting `lambda=mu'` and using
`kappa_(mu')=-2 chi_mu` gives exactly `(6)`.  Thus the apparent factor and
sign changes are entirely convention plus conjugation.

Every partition `lambda` of size at least five contains a size-five partition
and can be reached from it by a chain of Young covers.  Therefore `(4)` and
`(6)` imply `Phi_i(s_lambda)>0` whenever `s_lambda(Q_i)>0`.  If the length of
`lambda` exceeds the common residual-alphabet degree, every atom vanishes.
Thus a proof of `(4)` would upgrade the finite Schur face of THM-3110 to
all-degree Schur positivity and remove its residual histogram bank.

There is a representation-theoretic reason to take the content signal
seriously.  At degree five, `(6)` says that the central symmetric-group
operator detected by `Phi_i/B_i` is exactly

```text
(L_i/2) * identity + (1/4) * Omega,                  (7)
```

More precisely, with `Omega` the sum of all transpositions, `(6)` is
`(L_i/2) identity+(1/4) Omega`, because `Omega` acts on shape `mu` by
`chi_mu`.  A useful next construction is therefore the degree-`n` central
operator associated to `Phi_i`, followed by the `Q_i` Doob/Schur transform.
The target is positivity of each cover increment, not positivity of the raw
central operator alone.

## Exact proof-mechanism audit

For a cover `mu<nu` and atom `S`, put

```text
d_S(mu,nu)=s_mu(S)/s_mu(Q_i)-s_nu(S)/s_nu(Q_i).      (8)
```

THM-3110 gives `S<=Q_i` coordinatewise for every negative atom.  The cited
Schur-ratio theorem proves `d_S>=0`; the scout independently checked all
`6,945,664` negative-atom instances.  The desired cover law is exactly the
grouped inequality

```text
sum_(c_S<0) |c_S| d_S
    >= sum_(c_S>0, S!=Q_i) c_S d_S.                  (9)
```

Three cheaper mechanisms fail:

1. **Fixed pair transport fails.**  Intersecting all pairwise orders
   `d_N>=d_P` over the test universe gives maximum flows

   ```text
                  tight             wide
   I_1           17/36             23/36
   I_2           18/37             22/37.             (10)
   ```

   Thus `(9)` requires grouped cancellation; no shape-independent matching of
   negative atoms to positive atoms can prove it.

   Requiring in addition the theorem-friendly root order
   `S_negative<=S_positive` leaves **zero** transport edges in every
   bank/chamber intersection.  The one-row Schur-ratio theorem is therefore a
   genuine ingredient but cannot be iterated as a negative-to-positive atom
   matching.

2. **Positive quotient-Cauchy factorization fails.**  If
   `F=sum_R c_R C_(S_R)` and `C_Q` is the dominant dual-Cauchy kernel, then
   `F/C_Q` already has negative Schur coefficients at `(a,b)=(1,2)`:

   ```text
   [s_(6)] F_1/C_(Q_1)=-216,       [s_(6)] F_2/C_(Q_2)=-184. (11)
   ```

3. **Young-diamond supermodularity fails.**  For `I_1,(a,b)=(1,2)`, the
   diamond with bottom `(3,1,1,1,1)`, sides `(4,1,1,1,1)` and
   `(3,2,1,1,1)`, and top `(4,2,1,1,1)` has mixed `R`-difference

   ```text
   -75157440544709 / 243760997247100830.              (12)
   ```

   Hence cover increments do not propagate by ordinary lattice convexity.

These hostiles locate the missing object precisely: a bank-level
Jucys--Murphy or tableau injection proving `(9)` after the `Q_i` normalization.
It must couple several atoms at once and need not have a fixed atom transport
or a fixed diamond sign.

The incoming exact scout
`arbitrary-anchored-product-gamma-order5-transport-circuit-codex-20260802.md`
now explains why this grouping is forced.  After adjoining the reverse
dominance core to THM-3110's forward transport, each chamber has one
full-support minimal circuit: it annihilates every Schur flag through degree
four and exits positively in degree five, while deleting any circuit column
kills the kernel.  Its raw `2 x 2` coefficient minors already have both signs.
This is not a proved dependency, but it rules out both a proper subcircuit
decomposition and ordinary coefficientwise total positivity.  The row-ray
inequality `(5c)` should therefore be sought on that full rooted rank-four
circuit, with the chamber Newton cone retained.

There is also a concrete repo-local alternant route.  Each numerator in `(8)`
is the `2 x 2` evaluation determinant

```text
s_mu(S)s_nu(Q_i)-s_nu(S)s_mu(Q_i),                    (13)
```

an adjacent generalized-alternant minor.  Section 3 of proved THM-3089
(`THM-3089-logarithmic-moving-gap-cluster-cone-and-condition-number-boundary`)
shows how Schur/Vandermonde factorization and total positivity control a whole
generalized alternant uniformly.  It does not orient the signed sum `(9)`, but
it suggests the cheapest next exact test: expand the *grouped* bank of minors
`(13)` and ask whether it has a new dominant alternant row after the rank-four
edge-slide cancellations of THM-3110, rather than matching the old atoms.

A deliberately remote repo analogue reinforces the grouped formulation.
Proved THM-1143's `A_12` chipwalk has local edge slides for which arbitrary
tie ordering fails, while simultaneous-wall grouping preserves total mass and
a prefix invariant.  It is not a GMC dependency, but it suggests a precise
experiment here: label the rank-four Ewens insertion edges by the added box's
`(coarm,coleg,content)`, group equal-content/tableau slides, and search for a
prefix or central-content invariant.  The zero theorem-friendly pair-edge
count above is a counterindication to any ungrouped version of that analogy.

A related holotopy probe must be performed one layer earlier than the present
bank.  The exact constructor's `response_moment`, polynomial multiplication,
and final collection successively merge equal row indices and then sort the
row multiset.  Consequently, attaching formal `epsilon_j` labels to the rows
of `BANKS` would be noncanonical: the macro-row incidence has already been
forgotten.  The lawful next experiment is to retain a label on every numerator
and uncancelled denominator row through the three invariant summands, replace
its interval alphabet by `{r+epsilon_j}`, and only then compute the degree-five
augmented Schur minors.  Coefficientwise positivity before label collapse would
be meaningful; positivity after an arbitrary relabelling of the collected bank
would not.  No epsilon-lift claim is made here.

## Reproduction

```powershell
python3 04-computation/gmc_product_gamma_young_cover_scout.py
python3 -O 04-computation/gmc_product_gamma_young_cover_scout.py
```

Both runs must match
`05-knowledge/results/gmc_product_gamma_young_cover_scout.out` exactly.
LF-normalized SHA-256 at this checkpoint:

```text
script  920b1f789f59df7cbb4b0c525ffbe5a7d9dee8107c70d951f60f66da10821cb6
output  2150757883fefbf61003b27b5a434de0ca10f5ec14847a4734a99d9b0442eb11
```
