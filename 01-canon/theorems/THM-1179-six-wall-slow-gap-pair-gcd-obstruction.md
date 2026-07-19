---
id: THM-1179
title: Six-wall slow-gap graph debt, rooted Fano periods, and the quantitative chi7 pair floor
status: PROVED analytic necessary inequalities and FINITE-EXACT same-chi7 pair floor.  A six-comb cover of an a-slow gap obeys an adaptive graph/arboricity pair-gcd law; the complete graph closes every common-killer period with gcd(b_1,...,b_6)/a<397/432; every carrier-rooted Fano plane obeys a mixed pair/triple-period inequality; and same-chi7 pairs have the sharp global overlap floor 1/77.  These close large-period and low-harmonic-excess branches, not the full harmonic-crowded slow-gap cone or LRC(14)
source: codex-2026-07-18 post-THM-1176 frontier session
depends_on: [THM-1156, THM-1166, THM-1176]
related: [THM-856, THM-1153, HYP-7678]
script: 04-computation/lrc14_six_wall_slow_gap_pair_gcd_codex_20260718.py
output: 05-knowledge/results/lrc14_six_wall_slow_gap_pair_gcd_codex_20260718.out
---

# THM-1179 -- six-wall slow-gap graph debt

At the actual LRC(14) radius put

```text
D_s={t in R/Z: ||st||<1/14}.
```

Let `a<b_1<...<b_6` be distinct positive integers.  On a real lift, an
`a`-slow gap is

```text
G_k=[(14k+1)/(14a),(14k+13)/(14a)],    L=|G_k|=6/(7a).  (1)
```

The endpoints in (1) are not in `D_a`, because the danger teeth are open.
Assume that the other six combs cover `G_k`.  Write

```text
H=sum_i 1/b_i,                 delta=aH-1,
rho_ij=measure(D_bi intersect D_bj),
g_ij=gcd(b_i,b_j),             eta_ij=rho_ij(1-rho_ij). (2)
```

THM-1176 proves `delta>0` and supplies the slow-gap/toothpick recursion.  The
new point here is that the excess harmonic pressure must also pay for local
pair multiplicity, pair periods, and carrier-rooted Fano line periods.

## 1. The periodic endpoint ledger

We use both sides of one elementary periodic-remainder lemma.

> **Periodic remainder lemma.**  Let `f` have period `P` and mean `nu`.
> If `f` is an indicator, then every interval `J` satisfies
>
> ```text
> integral_J f >= nu|J|-nu(1-nu)P.                     (3)
> ```
>
> More generally, if `0<=f<=M`, then
>
> ```text
> integral_J f <= nu|J|+nu(1-nu/M)P.                   (4)
> ```

For (3), remove the full periods and let the remaining interval have length
`r<=P`.  It contains at least

```text
max(0,r-(1-nu)P)
```

of the indicator.  This is at least `nu r-nu(1-nu)P` on both sides of
`r=(1-nu)P`.  For (4), the remainder mass is at most
`min(Mr,nu P)`; its maximum excess above `nu r` is attained at
`r=nu P/M`.  Both constants are sharp for a single block in one period.

Apply (4) to one danger comb, whose mean is `1/7`, maximum is one, and period
is `1/b`.  Apply (3) to a pair indicator, whose mean is `rho_ij` and whose
period is `1/g_ij`.  For every interval `J` of length `x`,

```text
measure(J intersect D_b) <= x/7+6/(49b),
measure(J intersect D_bi intersect D_bj)
 >=x rho_ij-eta_ij/g_ij.                               (5)
```

This derives the pair discrepancy used below rather than treating it as an
asymptotic `O(1/g)` term.  Open, closed, or half-open endpoint conventions
give the same measures in (3)--(5).

We will also use the exact reduced pair formula from THM-1166.  For coprime
`x<y`, with `F(r)=r(14-r)` for `0<=r<14`, it is

```text
rho(x,y)=[4xy+F((x+y) mod 14)-F((y-x) mod 14)]/(196xy). (6)
```

In addition to the lower floor `rho>=1/91`, distinct speeds satisfy

```text
rho<=1/14,                 eta<=13/196.                 (7)
```

Indeed the folded correction in (6) is at most `49`.  If `xy>=5`, then
`49<10xy`, which makes (6) strictly less than `1/14`.  The only coprime rows
with `xy<5` are `1:2`, `1:3`, and `1:4`, with densities `1/14`, `1/21`, and
`1/28`.  Since `rho(1-rho)` is increasing on `[0,1/14]`, (7) follows.

## 2. The adaptive graph theorem

Let `Gamma` be any nonempty simple graph on the six killer-comb labels.  Put

```text
d(Gamma)=max_(A subset V, |A|>=2) |E(Gamma[A])|/(|A|-1),
R_Gamma=sum_(ij in E(Gamma)) rho_ij.                    (8)
```

Thus `d(Gamma)` is the maximum induced edge density with the tree rank in the
denominator.  It is one for a nonempty forest, `3` for `K_6`, and is the
quantity that a pair-only majorant actually sees.

> **Theorem A (slow-gap graph debt).**  Every six-comb cover of (1) obeys,
> for every nonempty graph `Gamma`,
>
> ```text
> a sum_(ij in E(Gamma)) eta_ij/g_ij
>   +(6 d(Gamma)/49) delta >= (6/7)R_Gamma.             (9)
> ```

**Proof.**  At `t in G_k`, let `A(t)` be the active killer labels and
`C(t)=|A(t)|`.  Coverage says `C(t)>=1`, and (8) gives the pointwise graph
majorant

```text
|E(Gamma[A(t)])| <= d(Gamma)(C(t)-1).                  (10)
```

Let `S=sum_i measure(G_k intersect D_bi)` and let `P_Gamma` be the sum of the
local pair intersections over the edges of `Gamma`.  Integrating (10), then
using the one-comb half of (5), gives

```text
P_Gamma <=d(Gamma)(S-L)
 <=d(Gamma)[-L/7+(6/49)H]
 =(6d(Gamma)/49)(H-1/a).                               (11)
```

The pair half of (5) gives

```text
P_Gamma >= L R_Gamma-sum_(ij in E(Gamma)) eta_ij/g_ij. (12)
```

Combine (11)--(12), multiply by `a`, and use `aL=6/7`.  This is (9). ∎

This is the slow-gap analogue of THM-1166's adaptive forest certificate, but
its endpoint budget is the **excess** `H-1/a`, not the uncentered harmonic
sum `H`.  Forests remain available with `d=1`; denser graphs are charged
exactly by their induced density rather than discarded.

For `Gamma=K_6`, the pointwise majorant behind (10) is

```text
1_(C>=1) <= Q_6(C):=C-(1/3)binom(C,2),                 (13)
```

whose values at `C=0,...,6` are

```text
0, 1, 5/3, 2, 2, 5/3, 1.                              (14)
```

The coefficient `1/3` is optimal in this quadratic family: the value at
`C=6` forces `6-15 alpha>=1`.

## 3. Complete-pair consequences and the common-period closure

THM-1166 proves that every triple of distinct speeds has pair-density sum at
least `1/24`.  Summing over the twenty triples of six vertices counts each
of the fifteen pairs four times, so

```text
R:=sum_(i<j)rho_ij >=5/24.                              (15)
```

The complete-graph instance of (9) is therefore

```text
a sum_(i<j) eta_ij/g_ij +(18/49)delta >=5/28.           (16)
```

Using (7) gives the simpler, weaker gcd-only consequence

```text
a sum_(i<j) 1/g_ij >=(35-72delta)/13.                  (17)
```

Thus (17) is nonvacuous throughout the low-excess cone `delta<35/72`.  If
all pair gcds are at least `G_pair`, then

```text
G_pair/a <=195/(35-72delta).                            (18)
```

THM-1191 subsequently refuted the proposed black-box improvement obtained by
averaging a uniform four-comb floor: the self-similar `13`-adic quartet makes
that route too weak.  Thus `5/24` is the current proved input to (16)--(18).
A stronger **direct** six-speed floor is not refuted, but would have to retain
six-way compatibility, the exponent stalk, or the nonuniform graph/Fano data
rather than averaging one quartet constant.

There is a stronger conclusion when the six killers themselves share one
period.  Put

```text
G_0=gcd(b_1,...,b_6).
```

Integrating (13) over the whole circle and using (15) yields

```text
measure(union_i D_bi)
 <=6/7-(1/3)R <=397/504.                               (19)
```

The union is `1/G_0`-periodic.  A connected component in the periodic lift
has length at most the total covered mass in one period, hence at most
`397/(504G_0)`.  Because the open union contains the closed interval `G_k`,
the containing component is strictly longer than `L`.  Consequently

```text
G_0/a <397/432.                                        (20)
```

If coverage is weakened to almost-everywhere coverage or endpoints are
closed by fiat, (20) becomes non-strict `<=`.  In the actual open LRC
convention it is strict.  In particular, the entire branch `G_0>=a` is
empty.  Unlike THM-1166's seven-wall common-dilate result, (20) requires a
common divisor only among the six killers; the slow carrier `a` need not
share it.

## 4. The quantitative chi7 face

Strip the full `7`-adic factor and let

```text
epsilon(v)=Legendre_7(v/7^nu7(v)) in {+1,-1}.           (21)
```

THM-1156 shows that exact tooth seams are cross-colour.  The same colour also
removes the exceptional global pair-overlap well.

> **Lemma (sharp same-colour pair floor).**  If `epsilon(b)=epsilon(c)` and
> `b!=c`, then
>
> ```text
> rho(b,c)>=1/77.                                       (22)
> ```
>
> Equality occurs exactly when the reduced unordered ratio is `1:11` or
> `2:11`.

**Proof.**  The character in (21) is multiplicative, so a common gcd cancels
and the reduced coprime pair still has equal colour.  THM-1166's defect bound
for a reduced `x<y` is

```text
rho(x,y)>=1/49-1/(7y).                                 (23)
```

For `y>=20`, the right side is at least `13/980>1/77`.  The remaining finite
bank `y<=19` contains 58 equal-colour coprime rows.  Exact substitution in
(6) has minimum `1/77`, only at `(x,y)=(1,11),(2,11)`.  This finite bank is
replayed independently from (6) by the companion script. ∎

There is also a sharp floor for the full complement of the exact-seam graph.
For a reduced pair `x<y`, THM-1156's seam condition is `14|(x+y)`.  If a
pair is not seam-compatible, then

```text
rho(b,c)>=1/84.                                        (22a)
```

Same-colour pairs satisfy the stronger (22).  For cross-colour nonseam
pairs, (23) is greater than `1/84` once `y>=17`; at the first tail value it
is `10/833`.  Exact evaluation of the remaining 33 coprime cross-colour
nonseam rows has unique minimum `1/84` at reduced ratio `1:12`.

Consequently, if `Gamma_ns` is the graph of non-seam-compatible killer pairs,
with `N_ns` edges and induced density `d_ns`, Theorem A also gives

```text
a sum_(ij in E(Gamma_ns)) eta_ij/g_ij
  +(6d_ns/49)delta >=N_ns/98,

a sum_(ij in E(Gamma_ns)) 1/g_ij
  >=(2N_ns-24d_ns delta)/13.                           (22b)
```

Thus the exact seam graph and its complement both have quantitative
consumers: seam-compatible edges carry the endpoint/third-support alternative
of THM-1156, while every complementary edge carries the stronger global-pair
floor (22a).

Let the two colour classes among the six killers have sizes `p<=q`, and put

```text
E_chi=binom(p,2)+binom(q,2).                            (24)
```

Take `Gamma` in Theorem A to be the disjoint union `K_p union K_q`.  Its
induced density is `q/2`: concentrating an active set in the larger clique
attains this value, while splitting it between the two cliques cannot exceed
it.  Equations (9) and (22) give

```text
a sum_(epsilon(b_i)=epsilon(b_j)) eta_ij/g_ij
  +(3q/49)delta >=6E_chi/539.                           (25)
```

Using (7) once more,

```text
a sum_(same colour) 1/g_ij
 >=24E_chi/143-(12q/13)delta.                           (26)
```

This is a genuine metric consumer of the seam colour: same-colour pairs not
only fail to abut exactly, but must spend a larger global pair mass and hence
more local endpoint/gcd budget.  It does not assert a uniform seam
penetration; THM-1156's crossing/containment caveat remains essential.

## 5. A carrier-rooted Fano period inequality

Label `{a,b_1,...,b_6}` by the seven points of an arbitrary Fano plane, with
`a` distinguished.  Three lines pass through `a` and contain two killers;
four lines avoid `a` and contain three killers.  For a line `ell`, set

```text
B_ell=ell minus {a},      r_ell=|B_ell| in {2,3},
G_ell=gcd(b_i:i in B_ell),
R_ell=sum_({i,j} subset B_ell) rho_ij,
nu_ell=r_ell/21-R_ell/3,
e_ell=nu_ell(1-3nu_ell).                               (27)
```

> **Theorem B (rooted Fano clock).**  Every labelled Fano plane obeys
>
> ```text
> (6/(7a))(1/7+R/3) <=sum_ell e_ell/G_ell.              (28)
> ```

**Proof.**  On `G_k`, the carrier comb is inactive.  If `c_ell(t)` is the
number of active killers on line `ell`, define

```text
f_ell(t)=c_ell(t)/3-(1/3)binom(c_ell(t),2).             (29)
```

The values at `c_ell=0,1,2,3` are `0,1/3,1/3,0`, so
`0<=f_ell<=1/3`.  Every killer lies on three Fano lines and every killer pair
lies on exactly one line.  Therefore the exact pointwise identity is

```text
sum_ell f_ell=C-(1/3)binom(C,2)=Q_6(C).                (30)
```

The mean and period of `f_ell` are `nu_ell` and `1/G_ell`.  By (4),

```text
integral_(G_k) f_ell <=L nu_ell+e_ell/G_ell.           (31)
```

Coverage and (13) give `L<=integral Q_6`.  Incidence in (27) also gives

```text
sum_ell nu_ell=6/7-R/3.                                (32)
```

Sum (31), use (30)--(32), and rearrange.  This is (28). ∎

The exact density floors make (28) explicit.  On a three-killer line,
`R_ell>=1/24`, so

```text
nu_ell<=65/504,        e_ell<=6695/84672.               (33)
```

On a two-killer line, `rho>=1/91`, so

```text
nu_ell<=25/273,        e_ell<=550/8281.                 (34)
```

The function `nu(1-3nu)` is increasing throughout these ranges.  Combining
(15), (28), (33), and (34) gives, for every rooted Fano plane,

```text
107/(588a)
 <=(6695/84672) sum_(4 lines off a) 1/G_ell
   +(550/8281) sum_(3 lines through a) 1/G_ell.         (35)
```

If all seven active line gcds are at least `G`, then

```text
G/a <=1844255/650988 <2.834.                            (36)
```

The `chi_7` colour slightly strengthens this genuinely mixed-period probe.
Any perfect matching of the six killers can be realized as the three Fano
lines through `a`.  Pair within each colour class as far as possible.  Since
the two class sizes have the same parity, this produces at least two
same-colour matched pairs.  For such a pair, (22) replaces (34) by

```text
nu_ell<=1/11,          e_ell<=8/121.                    (37)
```

Choosing that labelled Fano plane and again assuming all seven active line
gcds are at least `G` sharpens (36) to

```text
G/a <=222893927/78769548 <2.830.                        (38)
```

The numerical improvement is small; the structural gain is that (35)--(38)
constrain pair and triple periods even when the six killers have no useful
common gcd.

## 6. Endpoint guardrail: the discrepancy tax is indispensable

There is no positive local pair floor before the endpoint terms in (5) are
paid.  Take

```text
a=14,       k=2,       (b_1,...,b_6)=(18,19,20,21,22,23),
G_k=[29/196,41/196].                                   (39)
```

Direct endpoint sorting gives total single mass

```text
sum_i measure(G_k intersect D_bi)=850237/16959096,
```

but all fifteen local pair intersections vanish.  The occupied pieces occur
in the strict order

```text
b20, b19, b18, b23, b22, b21, b20, b19,                (40)
```

with positive gaps between consecutive pieces.  Thus replacing (12) by the
false strengthening `P_Gamma>=L R_Gamma`, or dropping the `eta/g` terms, is
impossible even on a slow gap.  This packet does not cover (39); it is a
guardrail on proof technology, not a counterexample to the slow-gap
conjecture.

For an actual cover, the strict inequality `delta>0` also has a topological
reading.  Equality in the one-comb sum would force zero pair mass.  A finite
family of relatively open, pairwise-disjoint nonempty pieces cannot cover a
connected interval unless one member covers it alone, which the one-comb
bound forbids.  The parity clock of THM-1176 gives additional arithmetic
rigidity, but strictness itself is already visible in the owner-event
hypergraph.

## 7. Tournament, Paley, and carrier audit

On runner vertices, the pairwise observable

```text
O(i,j)=a/b_i-a/b_j
```

with positive-sign gauge gives the speed-order tournament.  It is transitive:
score histogram `0,1,2,3,4,5`, no directed cycles, six singleton SCCs, and
one Hamiltonian path, with no ties.  It preserves the ordering used by the
H-drift of THM-1176, but destroys every `rho_ij`, `g_ij`, endpoint position,
multiple owner, and Fano line period.

The Paley tournament on seven residue labels is not a faithful replacement.
The stripped residues of speeds may repeat, its orientation does not preserve
exact-seam indices or the local-cover predicate, and quotienting to the
binary character retains only the two predicates actually used here:

```text
same colour forbids an exact seam;
same colour forces rho>=1/77.                           (41)
```

Likewise, ordering Fano lines by `G_ell` or `e_ell/G_ell` produces a
transitive tie-resolved tournament but loses line incidence and the sum in
(28).  No new LRC predicate is preserved by either tournament quotient, so
neither is promoted to the proof carrier.

The challenged vertex sets are runners, danger arcs, slow-gap clocks,
labelled tooth endpoints, wall-crossing events, pair obligations, and Fano
lines.  The smallest object retaining all proved predicates is

```text
slow-gap clock k
 + interval hypergraph of endpoint/owner events
 + selectable graph Gamma with (rho_ij,g_ij)
 + rooted Fano incidence with (nu_ell,G_ell)
 + chi7 colour sidecar.                                (42)
```

This is a one-dimensional Kakeya-needle view: a protected slow gap is the
needle, tooth intervals are its labelled cuts, and pair/Fano periods measure
how much apparent global density can drift off the needle through its two
endpoints.  The owner chronology, not a bare runner tournament, decides
coverage.

## 8. Exact replay and honest frontier

The dependency-free `Fraction` replay checks:

```text
7       Q_6 activation values,
4       chi7 colour splits and their exact graph densities,
58+33   reduced same-colour and cross-colour-nonseam low-ratio rows,
4,004   sharp periodic-envelope samples,
3,612   exact one-comb interval bounds,
378     exact pair means from independently intersected teeth,
48,762  exact pair interval-discrepancy bounds,
64      carrier-rooted Fano activation masks,
1       zero-local-pair endpoint guardrail.
```

Normal and optimized outputs are byte-identical.  Frozen hashes are

```text
source  8215d074541a8110bd97b45e9d8cba73ea55497cf647b677b89ea35b074e56cc
output  25ad2f2a03b72e3f02f7ca8162eeb618fad82a1dafc79ebd4efffcec797403eb
```

What is closed is the common-killer-period branch (20), the low-excess
pair-gcd branches (16)--(18), the same-colour branch (25)--(26), and the
large rooted-Fano-line-period branch (35)--(38).  What remains open is the
mixed-period harmonic-crowded cone in which enough pair discrepancies are
small, the Fano line gcds are sufficiently mixed, and the endpoint owner word
still covers every slow gap.  No universal six-comb slow-gap noncoverage,
CrownCollapse14, or LRC(14) conclusion is claimed.
