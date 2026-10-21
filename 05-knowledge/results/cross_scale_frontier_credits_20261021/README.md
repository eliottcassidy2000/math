# Cross-scale frontier credits and optimal permanent capacity

2026-10-21. **PROVED, scoped:** complete fixed-frontier inverse compiler,
three-phase scale law after the explicit frontier cutoff, mod105 period
stratification, and the weighted permanent-capacity criterion below.
**FINITE-EXACT:** optimal extracted budgets273 and474 complete the frozen
N4096 and N8192 sibling banks with79 and154 distinct interaction owners.
The universal native153N inequality for every N>=4096, uniform control of
the remaining R', and arbitrary-source Collatz entry remain **OPEN**.

## 1. Inheritance and what changes

The inherited native source is a positive odd integer s not5 modulo8.
P uses branch a=2 at s=1 modulo8, otherwise a=1, then strips B(x)=4x+1
from (3s+1)/2^a. Its first arrival in the fixed627-point grounded frontier V
has an actual P-edge depth and an integer braid clock. These are not the
fully odd Collatz U clock. Sources already in V stop at depth and clock0.

The closest mechanisms are the complete inverse classification and the
at-most-one-member-per-branch annulus lemma in
[sibling guards](../sibling_annulus_guards_20261020i/README.md), the affine
phase law in [native phase lifts](../session18_native_phase_lifts/README.md),
and source-owned temporal handoffs in
[residual flow](../annulus_residual_flow_20261020k/README.md).
The [previous external extension](../residual_external_credit_20261020k/README.md)
used94 distinct weight-3 negative pairs, totaling282, to pay94 individually
omitted intervals at N4096. Its finite certificate remains valid.
The present improvement reoptimizes the enlarged bank jointly, rather than
requiring one permanent owner for each previously omitted whole interval.

Canonical hostile: native N7 at actual depth8 exceeds153N by13. Corrected
near miss: the94 minimum whole-job omissions inside the sibling bank is not
the minimum permanent capacity needed after enlarging that bank. Least-used
sidecars: the exact negative weight, its unique pair address, and the
mod105 register period, which need not equal the charge period.

Board: inverse branch address / annular bands / temporal capacity /
transferred remainder / joint prime registers / full braid clock.
The fixed-frontier depth-one supply is bounded; this does not prove that
the sibling deficit grows without bound or that all-scale payment fails.

## 2. Complete all-N compiler, with literal first-arrival guards

Let A_N(s) be the inherited Eisenstein charge, Q(x,y)=x²-xy+y² and
beta(v,w)=Q(v+w)-Q(v)-Q(w). The native palettes are

    L=((1,2),(-2,-1),(1,-1)), U=((0,1),(-1,-1),(1,0)).

For every integer N>=1, old parent u in V with u<4N, and a in{1,2}, put

    j0 = the least residue of 3-a-u modulo3,
    s_j = [2^a(3u+1)4^j-(2^a+3)]/9,  j=j0+3r, r>=0.       (1)

Take the least r for which s_j>=4N. This gives exactly one candidate in
[4N,256N). Then exclude it if s_j lies in V. Keep the pair (u,s_j) only
when its literal signed interaction

    b_N(u,a)=beta(A_N(u),(-1)^a A_(64N)(s_j))               (2)

is negative. Its pair address is (u,s_j), not merely (u,a). It contributes
b_N(u,a) at every actual P-depth d>=1 and zero at d=0. Its old endpoint
has route (u,0,0); its new endpoint has route (u,a+2j,1).

Proof of completeness: divisibility of the numerator in(1) by9 gives
exactly j=j0 modulo3. The formulas have the required primitive first
branch. Consecutive values satisfy

    s_(j+3)=64s_j+c_a, c_1=35, c_2=49.                    (3)

Thus two candidates cannot both lie in[4N,256N). The initial candidate is
positive and below256N whenever u<4N: j0<=2 gives
s_j<=(64(3u+1)-7)/9<256N. If a previous candidate is below4N, it is an odd
integer at most4N-1; its next value is at most256N-64+c_a<256N.
This proves existence as well as uniqueness. Inverting the literal
native edge gives the stated clock; excluding V prevents a false depth1
assignment to a source already stopped at depth0.

Different pairs may share an old source. They remain distinct quadratic
interaction terms, so this does not reuse source norm. The pairs are
old/new and therefore disjoint from every all-new sibling pair address.
The compiler returns at most2|V|=1254 pair addresses. This cap is a direct
specialization of the inherited same-parent annulus lemma, not a new
inverse mechanism. Every charge has norm0,1 or3: complete lower-palette
triples vanish, and the remaining palette sums give those three norms.
Consequently |beta|<=6 and the total negative mass is at most7524. The
finite exact ray values below are substantially sharper.

This bounded supply is an obstruction to treating a fixed finite set of
lifted credits as a growing supply proportional to N. It is not a proof
that the actual required supply is unbounded; that comparison remains open.

## 3. Exact scale phases and the separate register clock

Assume N>=15917, equivalently4N>max V=63665. All627 parents are old and
all annular candidates lie outside V, so the complete catalogue has1254
addresses before testing negativity. The bijection between catalogues at
N and64N is

    (u,s,a,j) -> (u,64s+c_a,a,j+3).                     (4)

The six native annular bands have endpoints4N,8N,16N,32N,64N,128N,256N.
Every source is odd and every endpoint even. Since0<c_a<64, (4) preserves
its band exactly. Also64=1 modulo3, so s modulo3 rotates by c_a modulo3.
For R(x,y)=(-y,x-y), the band table therefore gives

    A_(4096N)(64s+c_a)=R^(c_a) A_(64N)(s).               (5)

The fixed old charge obeys A_(64N)(u)=A_N(u). The braid clock increases
by6 and its sign is unchanged; the first-arrival depth remains1.
After three lifts, R^(3c_a)=I, so all corresponding scalar weights,
negative-owner counts and total negative mass repeat. The scalar period
therefore divides3. Actual source labels grow and actual clocks do not
repeat. The complete-bank theorem is not asserted below N15917: for
example, N4096 has1206 candidates but its three-lift scale has1254.

For the exact ray N=4096*64^k, k>=1, the three repeating phases are:

| k modulo3 | Negative owners | Available negative mass |
|---|---:|---:|
|1|325|708|
|2|375|946|
|0|383|858|

These three distinct masses give exact aggregate period3 on this ray.
The finite computation authenticates three base phases; equations(4)-(5)
prove their continuation to every k>=1. It does not enumerate the other
annulus interactions or assert that these masses cover their deficits.

The joint guards carry more information. For T_c(s)=64s+c, modulo105,

    T_c²(s)=s+65c,
    T_c³(s)-s=63s+66c,
    T_c⁶(s)=s.                                         (6)

The exact register period is3 when a=1 and s=0 modulo5, or a=2 and
s=2 modulo5; it is6 otherwise. Indeed (6) gives the period3 criterion;
period1 is impossible because21 does not divide c, and period2 is
impossible because65c is70 or35 modulo105. These exclude every proper
divisor of6 except3. Compatible mod15 and mod35 registers jointly encode
this one mod105 state, not525 independent states.

The full-denominator reader computes the numerator in(1), with j+3k,
modulo9m and only then divides by9. It is checked at m=15,35,105.
The braid clock is a+2j+6k, so a six-step register return adds36 rather
than resetting the clock. The inherited
[finite-observer clock theorem](../clock_channels_20261021/FINITE_OBSERVERS_INTEGER_CLOCKS.md)
already proves the broad period-six/mod840 versus unbounded-clock
phenomenon. Equation(6) refines it to the exact period3/6 strata needed
by this compiler; it does not claim a new general clock-loss principle.

## 4. Sharp permanent capacity, not one credit per omitted job

Let the fixed sibling bank have p_d positive and b_d negative active
intervals. Its signed debit is J_d=3(p_d-b_d). Every sibling positive
interval starts at depth>=1. Let E be any set of distinct negative
old/new pairs from section2 and let

    K=sum_(e in E) -beta_e, C_d=-K for d>=1, C_0=0.      (7)

Then the entire sibling bank, with every whole positive interval retained,
is source-owned payable exactly when

    J_0<=0 and K>=max(0,max_(d>=1) J_d).                 (8)

Necessity is the scalar sum at a maximizing cut. For sufficiency subdivide
each positive or negative sibling interval into three labeled unit
weights. A weight-w negative external pair has exactly w labeled units,
not extra copies of that pair. At each integer depth, (8) permits an
injective matching of active positive units to active negative units.
Keep assignments until either owner ends, then reassign released units.
The resulting half-open segments partition every positive unit's entire
interval, and no negative unit is used twice simultaneously. Finite
horizon endpoints suffice; unused permanent negative tail is retained
in the ledger. This proves whole-interval completion, although different
units or temporal pieces of a single interval can use different owners.

The deterministic receipt supplies the actual segmented assignments,
not just the capacity inequalities. It reconstructs the original sibling
bank and full source rows and matches their authenticated hashes.

| N | Complete negative catalogue | Available mass | Minimum extracted mass | Minimum owners achieving it | Whole positive intervals paid | Unit segments |
|---:|---:|---:|---:|---:|---:|---:|
|4096|376|849|273|79|8852|50670|
|8192|294|628|474|154|17142|98289|

The maximum original J is273 at depth9 and474 at depth16, proving the
minimum extracted masses. At N4096 the chosen owners comprise twelve
weight6 and67 weight3 pairs. The78 largest available weights sum to270,
so79 owners are necessary as well as sufficient. At N8192 the chosen
owners comprise ten weight6,127 weight3,16 weight2 and one weight1 pair.
The153 largest weights sum to473, so154 owners are necessary and
sufficient. These owner minima are inside the complete negative depth-one
frontier catalogue, not among all possible external interaction families.

The previous94 omitted jobs are a bank-relative whole-job statistic.
They do not force282 units of permanent capacity: the new matching may
reassign the old negative owners too. The273 optimum does not contradict
the previous94 minimum omissions without external owners. Likewise the
N8192 sibling-bank optimum159 omissions does not force477 permanent
units;474 suffice after the enlarged joint matching.

## 5. Exact replacement ledger and the untouched R'

Write the original physical identity as

    Delta E_d=A_d+I_d, A_d+M_d=153N,
    I_d=R_d+J_d.

Here R contains every old/new and new/new interaction outside the original
sibling-pair addresses. Extract the selected external terms exactly once:

    R'_d=R_d-C_d, J'_d=J_d+C_d,
    I_d=R'_d+J'_d, J'_d<=0.                            (9)

Thus I_d<=R'_d, and R'_d<=M_d remains a sufficient obligation.
At d>=1 the new remainder is R+273 or R+474, not R-273 or R-474.
Retaining J' gives exactly the original energy, with no manufactured
improvement. Every source diagonal and every other interaction remains.

This is a replacement certificate, not an additional credit layered onto
the previous94 pairs. Those94 pairs are a subset of the N4096 catalogue;
their282 mass is included in849. Relative to that already-used set, only
567 of the full catalogue's mass is fresh. Either use the old282 bank,
replace it by the new optimal273 bank, or explicitly subtract overlaps
before augmentation. Never add282+849 as if they were disjoint.

All257 physical cuts at N4096 give max(R'-M)=-201401; the final physical
arrival cutoff authenticates every later finite horizon. This is nine
units sharper than the previous282-extraction bound, because nine fewer
negative units were moved out of R, not because physical energy changed.
The N8192 result here is full sibling-bank payment, not a new full R'
or153N census. Any independently supplied full-coherence calculation
must identify the same selected pair addresses before using(9).

Native N7 remains hostile. At actual depth8 the replay gives old energy19,
new-only energy1086, mixed term-2, arrived-new norm1054 and moving norm17.
The unchanged increment1084 exceeds153*7=1071. Repartitioning its pair
terms cannot cure that fact. N7 lies outside the intended N>=4096 target.
No bounded control or scale-periodic credit catalogue proves that target.

## 6. Portable verification and declared scope

Files: `check.py`, `receipt.json`, normal and optimized verification summaries.
The checker uses only the Python standard library and an explicit repo root:

    python -B check.py --repo-root /path/to/math --verify
    python -B -O check.py --repo-root /path/to/math --verify

Its inputs are pinned in `PINS` and repeated in the receipt. The complete
N4096/N8192 sibling banks are rebuilt, including actual terminal, P-depth,
braid clock and joint registers for every retained source. Each external
owner is checked directly from its native edge. Every handoff segment is
checked for positive-owner coverage, both birth/end guards and absence of
negative-unit overlap. The receipt contains all selected source rows.

Independent direct scans for N=1,...,32 check inverse completeness; these
are finite controls, while section2 proves the all-N compiler. Large
symbolic controls include10^40+7 without enumerating that annulus. Source
scale and inverse arguments reject booleans, floats and negative inputs.
Hostiles include duplicate segments, premature credit, small-scale failure
of complete-bank period3, the nonmatching register/charge periods and
native N7. All checks are explicit exceptions and remain active under-O.
The source hashes identify evidence, not substitutes for these replays.

No new theorem ID is assigned. No literature priority claim is made.
The honest next comparison is the native sibling deficit against available
owned capacity at unbounded scales, together with the separate bound on
R'. Bounded fixed-frontier capacity by itself neither proves nor refutes
that comparison.
