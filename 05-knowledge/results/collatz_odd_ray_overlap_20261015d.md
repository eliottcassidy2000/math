# Odd-ray telescoping gives a uniform first overlap payment

2026-10-15. **PROVED:** the exact odd-ray formula, a uniform half-payment
for the initial annular cone, its failure to iterate automatically, and a
fixed-depth lower bound. **FINITE-EXACT:** the declared controller experiments.
**CONDITIONAL:** the cone-with-defect and updated energy implications.
**OPEN:** iterated uniform weighted overlap, unconditional amplitude
nonvanishing, and universal ROOT entry. The first-block result does not close
any of those three targets. Self-audited; no independent referee is claimed.

[Program](../../04-computation/experiments/collatz_odd_ray_overlap_20261015d.py),
[exact output](collatz_odd_ray_overlap_20261015d.json).
The [terminal-cube companion](collatz_cube_context_20261015d.md) certifies
the exact optimum of an older failed controller, then repairs its budget
with one additional native edge into a retained 112-step suffix.

## 1. Inheritance, incoming work, and the board

The closest mechanism is the exact signed pair-fibre invoice in
[Gram/slack §2](collatz_gram_slack_20261015b.md), applied to the literal rows in
[annular overlap §2](collatz_annular_overlap_20261015c.md). The canonical hostile
is the halving-only label chi_3 in
[annular overlap §5](collatz_annular_overlap_20261015c.md): it
satisfies every halving identity while its annular energy is cubic. The corrected
near miss is MISTAKE-618 and the false branch-square-to-joint-contraction
inference in [the incoming audit](collatz_incoming_audit_20261015b.md).
The least-used sidecar is the **position along an inverse odd ray**, including
its residue phase and the two physical annular endpoints.

The ray S(u)=4u+1 and its period-six valuation law are inherited, not new:
see [the October 5 canonical section](collatz_fusion_helpers_20261005.md#3-q2-the-canonical-3-parent-section-and-its-induction-prices).
The new use below is a quantitative payment for the actual annular weights.
The October 13 [wavelet packet §3](collatz_wavelet_energy_20261013d.md)
already reduces sibling triples by their median; its coefficient norm and the
present native-terminal Gram norm are distinct objects and neither is silently
substituted for the other. The fixed-source finite-attainment theorem F1 in the
address/lattice note already removes the need for another compactness argument.

Incoming commits adcd7eed1f and 75c7c49746 were read. The latter independently
audits THM-4658, records exact character inheritance, retracts the stationary
sqrt-scale bias interpretation, and enlarges the *conditional* future-energy
margin after the 2^37 census. Its larger census was not rerun here. Its referee
report explicitly excludes that new census from its audited population. No
stationarity, white sampling, or empirical sign persistence enters our proof.

Anchor: arithmetic odd-step overlap. Niche: the cost of a ternary cyclic
potential. Wildcard: a Pell constant arising from an honest defect invoice.
Board: source multiplicity; ray residue phase; actual native clock; terminal
sign across physical bins; observed-edge glue; conductor cost. Each connection
below says which of these coordinates it retains.

## 2. Halving, then one odd block: an exact telescope

Use the shortcut map T, ROOT-supported t, and weight w(n)=3[n=2 mod3]-1
from [THM-4658, two-octave-memory-energy-descent-critical-lines](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md).
Put chi_3=(0,1,-1) on residues 0,1,2. For every N>=1,

    Delta(N)=sum_(N<n<=4N) w(n)t(n).

First send each source n=2^j u to its odd part, with sign (-1)^j and the
literal j-step receipt. The resulting coefficient of an odd u is exactly

    a_N(u) = -3 chi_3(u),        u<=2N,
             w(u),              2N<u<=4N,
             0,                 otherwise.                     (1)

For u<=2N exactly two consecutive powers of two lie in the annulus. Their
sum is (-1)^j[w(2^j u)-w(2^(j+1)u)]=-3 chi_3(u). In the outer shell only u
itself occurs. Thus a_N(u) has the sign of w(u) whenever it is nonzero.
Its exact l1 norm is

    L_0(N)=4N-3 floor((2N+3)/6)
                  +floor((4N+1)/6)-floor((2N+1)/6)
          =(10/3)N + epsilon_N,   |epsilon_N|<=2/3.              (2)

Next apply one accelerated odd block U(u)=oddpart(3u+1), whose native
clock is a(u)=v_2(3u+1), with sign (-1)^a(u). Stop at ROOT; skipping the
formal two-step return of u=1 does not change its sign.

Odd integers sharing U are exactly the rays

    u_j=S^j(b)=4^j b+(4^j-1)/3,
    b odd, b mod8 in {1,3,7}.                                  (3)

Indeed inverse valuations differ by an even integer. The unique primitive
b has valuation 1 or 2; adding S increases the valuation by exactly 2.
Hence all points of one ray have the same native sign to U(b).
Different primitive b have different U(b).

The useful identity is

    w(Su)-w(u)=3 chi_3(u),       S(u)=u+1 mod3.                  (4)

For a primitive b<=4N, let j be maximal with u_j<=4N. Summing (1) on its
ray therefore leaves the coefficient

    q_N(b)= w(b),                       u_j>2N,
            w(b)-w(S^(j+1)b),            u_j<=2N.               (5)

This is the coefficient *before* the final sign (-1)^a(b) to U(b).
For b>4N set q_N(b)=0. Equation (5) proves the full row identity; it does
not guess a terminal label. Since w takes only -1 and 2, q_N(b) always has
the sign of w(b), or is zero. Consequently the final terminal sign
(-1)^a(b) sign w(b) is independent of N.

**Retained information.** The map from original sources to ray heads loses the
valuation index if used alone. The controller retains n, its initial halving
count, the actual valuation at its odd part, and U(b). The path length is at
most floor(log_2 n)+2 for this first block. A ray head is a bookkeeping port;
the actual forward terminal is U(b). ROOT remains marked.

## 3. A logarithmic-time count and an all-scale half-payment

For p=4^j define

    h_j(Y)=floor((Y-(4^j-1)/3)/4^j),
    P_r(l,h)=#{b: l<b<=h, b mod8 in {1,3,7}, b=r mod3}.

Each P_r is a sum of three arithmetic-progression counts modulo 24. Formula
(5) gives the following exact, O(log N)-term evaluation without enumerating
sources or testing their ROOT membership:

    L_1(N)=sum_(j>=0) sum_(r=0..2) [
       |w(r)-w(r+j+1)| P_r(h_(j+1)(4N),h_j(2N))
       +|w(r)| P_r(h_j(2N),h_j(4N)) ].                          (6)

Stop when h_j(4N)<1. All endpoints occurring before that stop are nonnegative.
The first term is the inner band in (5); the second is its outer band.

**Proposition O1 (PROVED).** For every N>=1,

    |L_1(N)-(16/7)N| <= 32(1+floor(log_4(4N)))+1.              (7)

Proof with explicit rounding costs: each P_r differs from its interval length
/8 by less than 3. The outer interval length differs from 2N/4^j by less
than 1, so its weighted error is at most 12+1/2. The inner interval length
differs from N/4^j by at most 5/4; its weighted error is at most 18+15/16.
Their sum is less than 32 per level. The omitted geometric tail is less
than 1: the first omitted p exceeds 3N+1/4. The density contribution is

    N sum_(j>=0) 4^(-j)
      +(3N/4) sum_(j>=0, j !=2 mod3) 4^(-j)
      = (4/3+20/21)N = (16/7)N.                               (8)

Here 1-4^(-3)=63/64 occurs for a precise reason: the ray dilates by 4 and
its residue phase repeats after three letters. This is a checked geometric
series, not a numerical association with 63. The remainder depends on the
integer digits; it was not discarded.

**Theorem O2 (PROVED, uniform first overlap).** For every N>=50,

    2 L_1(N)^2 <= L_0(N)^2.                                  (9)

More generally, take any finite collection of nonnegative physical weights
rho_i(N), all supported on N>=50, and any positive row denominators d_i.
Make row i the exact sum sum_N rho_i(N) Delta(N). Let G_0 be its Gram matrix
after all initial halvings and G_1 after the one odd block above. Then

    ||G_1||_entry1 <= (1/2) ||G_0||_entry1.                    (10)

In particular this holds for every dyadic annular mesh X=2^k, k>=6,
with its original multiplicities and integer bins. Thus the ordered-pair
opposite-sign overlap invoice satisfies sum min(P,N)>=||G_0||_1/4: at
least half the old Gram mass is paid. No ROOT convergence premise is used.

Proof: for 2^k<=N<2^(k+1), use E_k=32(1+floor((k+3)/2))+1 in (7).
The rational inequality

    2((16/7)2^k+E_k)^2 <= ((10/3)2^k-2/3)^2

holds at k=12; E_k/2^k and (2/3)/2^k both decrease thereafter, so it
holds for all k>=12 and all N in that octave. This proves (9) analytically
for N>=4096. Below that threshold the independent literal-row and floor-count
checks agree at every N. They identify exactly the twelve exceptions
{1,7,11,13,23,25,31,35,36,37,43,49}. This finite exact component extends (9)
to N>=50; the program records the complete failure set.
Both stages have one fixed sign per terminal across *all* physical N.
For a row family B with this sign property,

    ||sum_i B_i B_i^T/d_i||_1
        =sum_i (sum_e |B_i(e)|)^2/d_i.                        (11)

Nonnegative source weights give the row norms sum_N rho_i(N)L_j(N).
Apply (9) before squaring to get (10). Finally use the inherited exact
signed pair-fibre identity; this controller deletes no labels. Any additional
authenticated common extension can only improve (10).

For the standard square-root mesh, (2), (7), and the exact sum of squared
bin means also show

    ||G_0||_1 = (700/27)X^3 + O(X^2),
    ||G_1||_1 = (256/21)X^3 + O(X^2 log X),
    ||G_1||_1 / ||G_0||_1 -> 576/1225.                       (12)

Thus the asymptotic first payment is 649/1225 of the halving Gram mass.
The coarse-bin subtraction is X(b^2-1)/12 in sum_N N^2; because b^2<=X
it fits the displayed error. This is still **cubic** scale, not the desired
quadratic or logarithmic budget.

## 4. What breaks on repetition, and what information repairs the statement

The sign cone in (11) is not preserved by another odd block. After two
blocks, the coefficient at terminal 11 is

| physical N | original contributing sources (coefficient; native clock) | total |
|---:|---|---:|
| 4 | 9 (+1;3) | +1 |
| 5 | 9 (+1;3), 18 (-1;4), 19 (-1;4) | -1 |

This is an actual annular-row witness, not an arbitrary vector. At X=256,
the exact first and second Grams are respectively 3263584788/16 and
1856150820/16, and the second exceeds half the first; therefore
blind iteration of (10) is false even in the declared native family.
The script records the integers. The strongest survivor is native norm
nonincrease and the exact signed-overlap invoice.

This failure can produce useful *physical-clock* cancellation. For arbitrary
rows B, let R=sum_i ||B_i||_1^2/d_i and G=sum_i B_i B_i^T/d_i. For each ordered
terminal pair split the products B_i(e)B_i(f)/d_i into positive mass P_ef and
negative absolute mass N_ef. Then, independently of terminal transport,

    R-||G||_1 = 2 sum_(e,f) min(P_ef,N_ef).                    (13)

Thus a terminal whose sign changes across bins must retain its entire column
profile. Replacing that column by occupancy or total magnitude destroys a
possible payment. Equation (13) is the same elementary fibre identity on a
different index: physical bins instead of incoming terminal pairs. It does
not assume that every sign-changing column gives a strictly positive invoice.

There is a second exact boundary. Replace w by w_0(n)=3[3 divides n]-1.
Initial halving cancels the whole inner part and leaves only

    sum_(2N<u<=4N, u odd) w_0(u)t(u).                        (14)

Across N in [X,2X), this uses odd u in (2X,8X). U is injective on that entire
set: two distinct points in one U-fibre differ by at least S(u)=4u+1>8X.
Therefore one odd block pays **zero** Gram mass in this orthogonal mod-3
sector. It is incorrect to extend O2 uniformly to all mean-zero mod-3 weights.
The different channel budgets in the older wavelet note are essential.

## 5. Why a fixed number of odd blocks cannot settle every scale

**Proposition O3 (PROVED, representation-relative).** Use the controller
which strips the initial powers of two and then takes exactly r odd blocks,
stopping at ROOT. Do not add graph-union closure or other terminal relations.
For X=2^k, k>=2, and 0<=r<=k-2, its standard annular coarse Gram satisfies

    ||G_(k,r)||_1 >= X^3/(2*4^r).                            (15)

Choose original sources

    4X<n<=6X,   n=-1 mod 2^(r+1).

There are exactly X/2^r of them. They are odd and their next r valuations
are all 1, so their terminal is e=3^r(n+1)/2^r-1 and native clock r.
Each such terminal is unique among **all** original sources m<=8X.
To see this, in v=u+1 coordinates a valuation-one step multiplies v by 3/2.
A valuation >=2 instead gives at most (3/4)v+1/2. A path with any such
step has endpoint+1 at most (3/2)^r((8X+1)/2+1/3), while an initially even
source has endpoint+1 at most (3/2)^r(4X+1). Both are strictly below our
selected endpoint. Additional non-one valuations only lower the bound;
moving the first extra division earlier maximizes its additive term.
The r=0 case follows directly from initial halving. Other all-one paths are
injective. A ROOT-stopped path cannot have a selected large terminal.

Every bin lying in [3X/2,2X) includes all these sources with multiplicity b.
Their isolated Gram submatrix has no source cancellation and each |w(n)|>=1.
Its entrywise norm is at least (X/2)(X/2^r)^2, proving (15).
Consequently a budget B_k X^2 in this restricted representation requires
r >= (k-1-log_2 B_k)/2 whenever the bound applies. For fixed r the cubic
obstruction persists as k grows. This does **not** lower-bound the cost of
arbitrary controllers: retained interior edges and unequal-time joins can
remove the isolation, precisely the reason to keep graph glue.

## 6. A constructive reroute: overlap with a priced defect

Define the row norm Q(B)=(sum_i ||B_i||_1^2/d_i)^(1/2). A native signed
identification A is nonexpanding in Q. Suppose at a stage one supplies a
literal decomposition B=C+R on the same terminal coordinates, with a checked
native controller A satisfying Q(CA)<=rho Q(C), rho=1/sqrt(2). Then

    Q(BA) <= rho Q(B)+(1+rho)Q(R).                            (16)

This follows by triangle inequalities and native nonexpansion, with no
assumption on the signs of R. Theorem O2 supplies the initial positive
annular cone C and its A. A **later** cone or reencoding needs its own actual
native map and proof of that contraction; arbitrary changes of generators
or labels do not inherit it.

If a changing-controller construction can supply Q(R_j)<=C_* X at every
stage while preserving the original source rows, iteration gives

    Q(B_m) <= 2^(-m/2)Q(B_0)+(3+2sqrt(2)) C_* X.             (17)

After m>=log_2 X+O(1), the first term is O(X), since Q(B_0)=O(X^(3/2)).
Thus this **sufficient but unproved decomposition premise** would give an
O(X^2) Gram budget. Equivalently, a proportional defect epsilon contracts
when epsilon<3-2sqrt(2). The reciprocal Pell numbers occur because
(1+1/sqrt(2))/(1-1/sqrt(2))=3+2sqrt(2), an exact error-amplification factor.
This does not import the Q covariance hypothesis or a changing-basis norm.

There is a concrete numerical sufficient version. If the required native
decompositions exist for m=k+8 stages with Q(R_j)<=7X, then
Q(B_0)<=8X^(3/2) gives

    Q(B_m)<=[1/2+7(3+2sqrt(2))]X < (83/2)X,
    ||G_m||_1 <= Q(B_m)^2 < (6889/4)X^2 < 1798X^2.          (17a)

The rational comparison uses sqrt(2)<10/7. Thus a residual allowance of
7X would suffice for the incoming constant-budget route in §7. Producing
those decompositions and contracting native maps remains the open part;
(17a) is a fully priced target, not their construction.

There is a useful hostile for proposed phase-only models. S acts as one
cycle of length m=3^q on residues modulo 3^q, by
v_3(S^j(x)-x)=v_3(j). Define a bounded zero-mean forcing h along that cycle
as (m-1)/2 plus ones, one zero, then (m-1)/2 minus ones. Every solution of
p(Sx)-p(x)=h(x) has range (m-1)/2, and its minimum possible sup norm after
centering is (m-1)/4. Thus a mean-zero digit pattern has a cyclic potential,
but the potential cost can grow exponentially with ternary depth. Merely
recovering phase or showing a finite mean zero does not pay (16).
This is a test of a proposed model class, not a claim about the actual
Collatz residual forcing. It identifies a measurable sidecar: the potential
oscillation or Q(R_j), not just its conductor or mean.

The next concrete helper question is therefore: can actual row profiles be
split into certified contracting packets and a residual of Q-size O(X),
using (13) to retain their cancellation across physical clocks? It is weaker
than insisting on the original cone at every stage and stronger than a
formal reencoding that merely renames the unpaid residual.

## 7. Incoming margin: a less expensive sufficient constant target

The incoming interval program
[collatz_pole_certificate_margins_20261015.py](../../04-computation/experiments/collatz_pole_certificate_margins_20261015.py)
reports a rigorously tested constant-energy threshold above 1806.59 at
N=2^37. We checked its tail formula, absolute boundary term, and scope; we
did not rerun the large moments computation or its Arb certificate. Use the
strictly smaller budget 1800, not the rounded threshold 1806.6.

A new **PROVED implication**, conditional on the inherited exact seeds, is:
if for every k>=35 a finite controller gives the standard annular Gram

    ||G_D,k||_1 <= 1798 X^2,  X=2^k,                         (18)

then E_X<=1800 X^2 for every k>=35. The fine envelope contributes <=4/3 X^2,
so delta_k<=1800-2/3. THM-4658 gives

    a_(k+2)<=(a_k+sqrt(delta_k))/2+sqrt(3/X).

At a_k<=sqrt(1800), the reserved gap in the average is greater than
1/(6sqrt(1800)); sqrt(3/X) is smaller once X>=194400. The inherited seeds
E_(2^35)=93862178945193578345 and E_(2^36)=172802101944936995215 are both
below 1800 X^2. Induct on the two parity chains. The incoming conditional
certificate then yields a nonzero Mercator residue. Neither (18) nor the
iterated defect premise (17) has been proved. This route concerns amplitude,
not universal ROOT entry.

## 8. Computation, scope, and reproduction

The exact controls independently compare (6) with a literal annular row at
every N<=4096, compare full ray coefficients and strict integer-clock paths
at every N<=256, and evaluate (6) on three much larger integers without
materializing their source windows. They also check the all-one high-terminal
obstruction for k=4..9, every r<=k-2, and ternary potentials for q=1..7.

The fixed-depth census uses every source 1..8X for k=8,10,12,14 and odd-block
counts 0,1,2,4,8,16. Exact full Grams are computed at k=8 and at k=10 for
blocks 1 and 2; other entries explicitly remain row-envelope upper bounds.
A second census retains the whole observed native graph and compiles its
closure for k=8,10,12,14,16 and blocks 1,2,4,8. Full Grams are computed only
when at most 3500 columns remain. All source weights and integer clocks are
retained; no unknown terminal label is sampled. Native closure is checked
edge by edge by the inherited compiler. Integer product bounds guard NumPy's
int64 operations, with independently summed sample entries.

Observed first-block row constants approach the proved 256/21. Later fixed
blocks have mixed-sign column profiles. Graph glue improves the finite
bounds substantially, but no fitted rate is promoted to a theorem. For
example the glued k=16, eight-block row envelope is 423.4071 X^2, with
25,067 original sources still unresolved; that is an upper bound for a
finite Gram, not ROOT coverage and not an all-scale estimate.

```sh
OPENBLAS_NUM_THREADS=1 python3 04-computation/experiments/collatz_odd_ray_overlap_20261015d.py
OPENBLAS_NUM_THREADS=1 python3 -O 04-computation/experiments/collatz_odd_ray_overlap_20261015d.py
```

The remaining target has narrowed: establish a scale-uniform *iterated*
invoice, either directly on pair fibres or through a source-preserving
packet decomposition with the defect in (16) paid. The first arithmetic
block is now quantitatively paid at every scale; the subsequent phase and
clock changes remain explicit obligations.

Normal and optimized Python outputs match byte for byte. SHA256:
`2f92c7fb168c996881a024e8821e1af4b6a2dbcd2d0db80de26e9af52670bc35`.
