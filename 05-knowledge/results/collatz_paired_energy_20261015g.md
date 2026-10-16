# Paired reset energy and the cost of coupling rows

2026-10-15. **PROVED:** the global conformal allocation formula, its angular
debt, the joint-bin discounted bound, and the conditional energy implication.
**FINITE-EXACT:** 105 paired transitions in the inherited native census
`k=8..23`. **OPEN:** the uniform paired inequality, independent weighted
overlap, unconditional amplitude nonvanishing, and universal ROOT entry.
Self-audited; no independent referee is claimed.

## Inheritance and the changed target

The nearest proved mechanism is P1--P2 of
[native survivor packets](collatz_survivor_packets_20261015e.md). The hostile
is the actual one-reset failure at `k=23` in
[reset debt](collatz_reset_debt_20261015f.md); the corrected near miss is
discarding the original row weighting when counting coalescences. The
underused coordinate is the *joint vector of debts on the same physical bins*.
The old one-reset `7X` hypothesis restricted to `k>=38` remains OPEN: the
smaller counterexample does not refute it.

Anchor: weighted cancellation of paired native stages. Niche: six-colour
boundary extension after a safe singleton palette. Wildcard: the difference
between a small Hamiltonian cyclic space and access to its marked readout.
Board: original coefficients; common terminal fibres; literal clocks;
row directions; boundary lists; operation closure. The latter lanes are in
[compression and marked contexts](compression_contexts_20261015g.md).

The live remote was fetched at startup and during this session. It remained
at `8c630640e0`; there was no new incoming commit to merge. The inheritance
pass recovered THM-2856's earlier factorial/Laguerre rank-one boundary
defect and THM-3305's observer-visible one-edge update. We use their
mechanisms, rather than claim new versions of those results.

## 1. An exact optimum when rows share the contraction budget

Let one checked signed native map act on all rows. Put
`u_i=||B_i||_1`, `v_i=||(BA)_i||_1`, and
`||x||^2=sum_i x_i^2/b_i`, with `b_i>0`. Thus `0<=v_i<=u_i`.
Keep every source coefficient and literal path clock. Fractions of a source
weight may be assigned to a core and residual, summing to the original.

Fix `0<rho<1`, and require only the *global* core condition

    Q(CA) <= rho Q(C),                 B=C+R conformally.       (1)

Set

    eta=1-rho^2,
    F=||v||^2-rho^2||u||^2,
    beta_i=(v_i-rho^2 u_i)_+,
    BETA=||beta||.

**G1 (PROVED).** If `F<=0`, the minimum `Q(R)` is zero. If `F>0`, it is

\[
 r_* = \frac{\mathrm{BETA}-\sqrt{\mathrm{BETA}^2-\eta F}}{\eta}. \tag{2}
\]

It is attained by taking row residual masses
`r_i=(r_*/BETA) beta_i` from the authenticated unmatched survivors.
This is an optimum for (1) and conformal splits of this fixed map, not for
changing the native map or arbitrary signed decompositions.

**Proof.** Opposite-sign matching in each actual target fibre leaves total
mass `v_i`. Any `0<=r_i<=v_i` can be selected from it, giving exact core
norms `u_i-r_i`, image norms `v_i-r_i`, and residual image norms `r_i`.
For any other conformal split of row residual mass `r_i`, the core image
is at least `(v_i-r_i)_+`. Replacing `r_i>v_i` by `v_i` improves both the
residual cost and the core constraint, so no optimizer needs `r_i>v_i`.

The necessary and attainable constraint is therefore

\[
 F-2\langle v-\rho^2u,r\rangle+\eta\|r\|^2\le0.        \tag{3}
\]

For `r>=0`, Cauchy--Schwarz bounds the middle inner product by
`BETA ||r||`. The smaller root of the resulting quadratic is (2).
Equality is attained in the direction `beta`. Coordinatewise
`beta_i/eta<=v_i`, so the smaller root stays in the allowed box.
The discriminant is nonnegative: the convex quadratic is minimized at
`r=beta/eta`, which lies in the box, while the feasible choice `r=v`
has value `-rho^2||u-v||^2`. This also proves attainment. QED.

The radical is avoidable in budget checks. For a proposed `K>=0`, if
`F>0` and `K<BETA/eta`, then

\[
 r_*\le K\quad\Longleftrightarrow\quad
       (F+\eta K^2)^2\le4K^2\mathrm{BETA}^2.             \tag{4}
\]

If `K>=BETA/eta`, feasibility is automatic. These are rational/integer
comparisons on the frozen rows.

### The discount coefficient must be proved again

Global contraction does **not** imply P2's recurrence with residual
coefficient `1-rho`. For survivor allocations the valid statement is

\[
 Q(BA)\le\rho Q(B)+(1-\rho)Q(R)+\rho\mathcal A,
 \quad \mathcal A=\|u-r\|+\|r\|-\|u\|\ge0.             \tag{5}
\]

This is the triangle inequality followed by (1). The additional term is
the triangle defect of the row-norm vectors. It vanishes when the core
and residual norm vectors have the same direction. P1 avoids it by proving
the stronger coordinatewise inequality before taking the Euclidean norm.

An exact hostile is `u=(5,3)`, `v=(5,0)`, `r=(3,0)`, `rho=7/10`.
The core contracts because `4 <= (49/100)*13`, yet the incorrectly
discounted conclusion would be

    5 <= (7/10)sqrt(34)+9/10,

which is false (`49*34<41^2`). The signed rows `(5,0)` and `(3/2,-3/2)`
on columns `3,8` realize these norms under the actual path `3->5->8`,
with terminal clocks `2,0`. This recovers the native two-edge hostile from
[Gram slack](collatz_gram_slack_20261015b.md). These are test rows, not the
original Collatz annular rows.

G1 also does not repair the actual one-reset failure: the exact tests still
fail `7X` at `k=23`, stages 7, 8 and 9. Their optimum residual ratios are
approximately `7.60845`, `7.89396`, `7.17716`. This locates the remaining
problem in the observed action, beyond merely enforcing separate row budgets.

## 2. A cleaner paired stopping rule

Use exactly the original physical bins, source coefficients and native
graph-closure policy of the preceding notes. Write `E_j=Q(B_j)^2` for
the row-envelope energy; this `E_j` is **not** the dyadic `E_X=sum H(N)^2`.
Each ordinary stage exposes one first reset from each current terminal,
then closes the full observed graph. Pair two already closed stages.
After the checked `Q<7X` stopping event, pad with identities.

The new uniform candidate is

\[
          E_{j+2}\le\max\{E_j/2,49X^2\}.                \tag{6}
\]

This is a direct sufficient target, not a proved bound. It says: either
halve the current energy, or enter the stated energy ball. It permits the
failing single stages and avoids borrowing the rowwise discount formula
for a global core.

**FINITE-EXACT.** Equation (6) holds for all 105 consecutive two-stage
windows in the frozen `k=8..23` census, including both starting parities.
A final singleton is followed by one declared identity. Windows already
ending below the floor are automatic; the nontrivial test is on windows
ending above it: 74 windows. This candidate was formulated using the frozen
census, so this is not an out-of-sample prediction test. The largest ratio is

\[
 \frac{E_{14}}{E_{12}}=
 \frac{444344772471356261}{895543619307167669}
 =0.4961732325\ldots<\frac12,
 \qquad k=22.                                           \tag{7}
\]

Thus the margin is small. This is evidence for a precise hypothesis, not
an extrapolation to untested scales. The stronger affine inequality
`3E_(j+2)<=E_j+98X^2` also passes this finite census; no special arithmetic
significance is assigned to that empirical constant.

**G2 (PROVED implication; premise OPEN).** If (6) holds for the first
`k+1` paired stages at scale `X=2^k`, then

\[
 E_{2(k+1)}\le49X^2.                                    \tag{8}
\]

Indeed iteration gives `E_(2m)<=max(2^-m E_0,49X^2)`, while
`E_0<=64X^3`, so `2^-(k+1)E_0<=32X^2`.
This uses at most `2(k+1)` reset/closure stages. A reset contains a
variable number of literal shortcut steps. Their integer clocks are retained;
neither two resets nor one pair is renamed as one shortcut step.

The inherited fine-bin price then gives

\[
 D_X\le(49+4/3)X^2=(151/3)X^2.                          \tag{9}
\]

If the premise is established for every `k>=38`, statement 3 of
[THM-4658, two-octave memory energy descent](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md)
and its two exact seeds give `E_X<=52X^2` for all those scales.
The induction seam is valid for `X>=1024`: use
`sqrt(52)-sqrt(151/3)>5/(6sqrt(52))` and
`3*144*52/25<1024`. The seeds at `2^38,2^39` are respectively
`70317261195513340128655`, `23439502132608158412453` and are below `52X^2`.
This remains an energy/amplitude implication, **not universal ROOT entry**.
The tested observed graphs contain no nonROOT cycles. A uniform controller
must either retain that scope or authenticate cycle deletion with the
ROOT-supported zero label; it must not assume the Collatz conclusion.

### The weighted overlap that would prove the rule

For one paired signed map without deletion let `M_i=(u_i-v_i)/2`, exactly the mass matched against
opposite sign in its genuine terminal fibres. Then

\[
 E_j-E_{j+2}=2\sum_i M_i(u_i+v_i)/b_i.                   \tag{10}
\]

For `E_j>49X^2`, (6) is equivalent to an independent lower bound

\[
 \sum_i M_i(u_i+v_i)/b_i
       \ge\tfrac12\min\{E_j/2,E_j-49X^2\}.              \tag{11}
\]

The identity is proved; the needed uniform arithmetic estimate is not.
When `E_j>=98X^2`, the right side is simply `E_j/4`. It is this
source-weighted opposite-sign mass, not the number of mergers or the
marginal reset-length distribution, that an independent proof must deliver.

## 3. Discounted debt can retain its bin direction

For the original P1 allocations set `r_(i,j)=d_rho(u_(i,j-1),u_(i,j))`
and retain the nonnegative vector

\[
 d_m=(1-\rho)\sum_{j=1}^m\rho^{m-j}r_j.
\]

Coordinatewise P1 gives the proved refinement

\[
 Q(B_m)\le\|\rho^m u_0+d_m\|,
 \qquad \|d_m\|\le(1-\rho)\sum_j\rho^{m-j}\|r_j\|.       \tag{12}
\]

The first bound keeps which physical bins carry debt at different times.
Its square contains the actual inner products `〈r_i,r_j〉`; replacing them
by products of norms loses that information. Identity padding still has
residual equal to the actual final row norms. The exact rational computation
verifies `||d_m||<=7X(1-rho^m)` and a final bound below `7.5X` for every
`k=8..23`, `m=k+8`. No random signs or independence are invoked.

There is also a weaker quadratic discount route. For paired energies put
`a_j=(E_(2j)-E_(2j-2)/2)_+`. Then

    E_(2m) <= 2^-m E_0 + sum_j 2^-(m-j) a_j.

An independently proved bound of `49X^2(1-2^-m)` on this sum would give
`E_(2m)<=49.25X^2` at `m=k+8`; the same dyadic induction yields `E_X<=52X^2`
for the required tail. This alternative can tolerate individual bad pairs.
Defining the debt is not a proof that it is bounded.

## 4. ROOT separation does not explain away the residual

Fresh native runs at `k=12,14,16` retain signed ROOT coefficients `a_i`
and unresolved row norms `z_i`, giving the exact invoice

    Q(B)^2 = sum_i (a_i^2+z_i^2+2|a_i|z_i)/b.

At their stopping stages, the three terms divided by `X^2` are:

| k | ROOT square | unresolved square | cross term |
|---:|---:|---:|---:|
| 12 | 2.649610 | 16.651931 | 12.745461 |
| 14 | 1.001328 | 30.319069 | 8.498498 |
| 16 | 0.639436 | 27.155852 | 5.581785 |

The unresolved term dominates these cases; separating the known ROOT label
does not make the debt disappear. These are three finite controls, not a
claim about ROOT's contribution at `k=23` or in the limit.

## 5. Refining a physical bin has a second price

Let `v_N` be the full terminal vector of the original annular increment,
after one common native graph refinement. Retaining the initial halving
edges gives the exact three-label increment

\[
 v_{N+1}-v_N=w(4N+1)\phi(4N+1)-w(4N+2)\phi(2N+1)
                         +w(4N+3)\phi(4N+3),           \tag{13}
\]

where `phi(n)` is its signed terminal basis vector with the actual clock.
The terms at `4N+4` and `N+1` cancel because their weights agree modulo
three and their clocks differ by two. The middle term uses one halving.
The three absolute weights sum to four, so `||v_(N+1)-v_N||_1<=4`,
uniformly over these common native refinements.

For a bin `I` of length `b`, put `u_N=||v_N||_1`,
`q_I=sum_(N in I)u_N`, and `s_I=||sum_(N in I)v_N||_1`.
The exact row-envelope refinement invoice is

\[
 \sum_Nu_N^2-\sum_I s_I^2/b
 =\sum_I\sum_{N\in I}(u_N-q_I/b)^2
                         +\sum_I(q_I^2-s_I^2)/b.        \tag{14}
\]

Both terms are nonnegative. The first is at most
`(4/3)X(b^2-1)`, by the four-Lipschitz bound and the pairwise formula for
variance. The second is a distinct *within-bin sign-mixture price*.
It records cancellation between different physical rows, not a new native
merger. It cannot be discarded by importing the existing scalar fine-bin
variance theorem, which uses one shared terminal label vector.

For example, take the two abstract terminal rows `(L,1)` and `(L,-1)`,
`L>=1`, with the first coordinate marked ROOT. Their increment norm is two,
their norm variance is zero, but the envelope refinement gap is `4L+2`.
For either one **shared** second-terminal sign `z=+/-1`, the actual scalar
energy gap is only

    (L+z)^2+(L-z)^2-2L^2=2.

The large gap comes from separately maximizing the two rows with
incompatible labels. This is an abstract hostile to a generic envelope
inference, not a claimed family of actual annular source rows.

**FINITE-EXACT refinement test.** On the same native stage maps at
`k=8,10,12,14`, every dyadic bin width `1,2,...,X` was tested. All 130
paired comparisons satisfy (6), including the unchanged policy's final
identity padding. Every full vector increment satisfies (13)'s bound,
and (14) is verified exactly. At `k=14`, stage one and the default bin
width, the variance and sign-mixture terms are respectively
`0.01566519X^2` and `0.44979983X^2`. The latter dominates this example.
There is no all-scale bound on that second price here. Thus the finite
refinement success does not justify dropping bin scale or shared labels
from a proposed uniform theorem.

## Reproduction and next decisive work

    python3 04-computation/experiments/collatz_paired_energy_20261015g.py
    python3 -O 04-computation/experiments/collatz_paired_energy_20261015g.py
    python3 04-computation/experiments/paired_root_profiles_20261015g.py
    python3 -O 04-computation/experiments/paired_root_profiles_20261015g.py
    python3 04-computation/experiments/paired_bin_refinement_20261015g.py
    python3 -O 04-computation/experiments/paired_bin_refinement_20261015g.py

[Paired output](collatz_paired_energy_20261015g.json) includes hashes of every
inherited input, all 105 transitions, exact global-budget failures and joint
debt squares. Its 980 rational controls independently compare the optimum
test with feasible quarter-grid allocations. [ROOT profiles](paired_root_profiles_20261015g.json)
reproduce the frozen rows from actual native graphs, audit every internal
source edge and replay literal clocks on a fixed grid.
[Bin-refinement output](paired_bin_refinement_20261015g.json) retains all 130
comparisons and both exact prices in (14), with the same graphs and clocks.

The next independent target is a family of explicitly matchable opposite-sign
source packets proving (11), or a discounted version of it. The useful
recovered lesson from THM-3305 is to seek a small *actual* extra edge that
changes an observable fibre, rather than a smaller statistic of an unchanged
graph. The inverse-clock and palette lanes provide checks on the preservation
contract, not missing Collatz inequalities.
