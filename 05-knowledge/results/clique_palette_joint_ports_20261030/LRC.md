# Complete shared-bank phase classification

**Status: PROVED for arbitrary finite positive-integer cores and banks at rational radius; the listed 78 all-pair families and 2,184 dilation families use exact finite certificates. LRC(14) for arbitrary inputs remains OPEN.**

For \(0<r\le1/2\), let
\[
G_V(r)=\{t\in[0,1]:\|vt\|\ge r\text{ for all }v\in V\}.
\]
Retain its complete union of closed interval components, including isolated points. The earlier terminal-fees packet gives the one-speed first-arrival test, and the odd-eighths packet gives joint-bank arrival in special top windows. We apply the same mechanism to the complete phase cover of an arbitrary finite core and an arbitrary jointly dilated bank.

## An exact formula for every multiplier

Let \(G_C(r)=\bigcup_i[L_i,R_i]\). Let \(B=G_A(r)\), regarded periodically modulo one. For nonempty B define \(g_B(\phi)\) as the forward circular distance from phase \(\phi\) to B, with distance zero inside B.

**L1.** For every positive integer d,
\[
\boxed{G_{C\cup dA}(r)\ne\varnothing
\iff \exists i:\quad d(R_i-L_i)\ge g_B(\{dL_i\}).}
\]
Whenever the inequality holds, an actual common time is
\[
t=L_i+\frac{g_B(\{dL_i\})}{d}.
\]

Indeed the bank condition is exactly \(dt\bmod1\in B\). As t traverses a core interval, its unwrapped bank phase starts at dL_i and moves forward a distance d(R_i−L_i). It hits B iff it reaches its first periodic arrival. This proves both directions and preserves the actual time, including equality at endpoints. If either complete safe set is empty, every d fails.

Let Q be the least common multiple of the denominators of the L_i. In each residue class d≡s mod Q, all forward gaps in L1 are constant. For each positive-width component, its condition is the exact threshold
\[
d\ge \frac{g_B(\{sL_i\})}{R_i-L_i}.
\]
For a point component it is either always true or always false in that residue, according as the gap is zero or positive. Taking the minimum threshold gives a complete residue-wise classification, with the residue restriction still enforced.

**L2 (finite head or purely periodic obstruction).** Suppose both safe sets are nonempty. If the largest core width W is positive, let h be the largest complementary circular gap of B. Then every
\[
d\ge D=\max(1,\lceil h/W\rceil)
\]
is accepted. Thus testing 1,…,D−1 gives the *entire* rejection set. The proof is \(g_B(\phi)\le h\) for every phase and the use of a width-W core interval. If the core consists only of points, the rejection set is purely periodic modulo Q instead.

This is a finite-to-infinite certificate. It does not infer an all-size conclusion from a sampled prefix. Once the complete core and bank covers have been compiled, a query costs at most their product in rational comparisons, independently of d's numerical value apart from integer bit complexity. Compiling the covers still depends on the original speeds; no strongly polynomial general LRC algorithm follows.

## Thirteen-speed families at radius 1/14

The checker certifies, for every choice
\[
C\subseteq\{1,\ldots,13\},\ |C|=11,\qquad
A\subseteq\{1,\ldots,8\},\ |A|=2,
\]
and every positive integer d,
\[
\exists t:\quad \|vt\|\ge1/14\quad(v\in C\cup dA).
\]
There are \(\binom{13}{11}\binom82=2,184\) declared families. Every finite rejection list is empty. The largest proved tail threshold is 23. Individual C, A, and tail certificates are in lrc_families.json. Some small d give overlaps, hence fewer than thirteen distinct speeds; all are covered.

The native comparison universe independently constructs the full safe set for all two-element C,A⊂{1,…,6}, all r∈{1/5,1/4,1/3}, and d=1,…,16: 10,800 comparisons. Every returned first-arrival time is checked against the original modular inequalities. A hundred-digit multiplier is checked by its literal rational witness without enumerating its teeth.

### Live integration: all positive pairs, not only fixed ratios

Incoming commit fd5a8148a4 supplies the stronger [two-speed compiler and two-tail conductor](../precision_lifts_ramsey_deficits_20261030/LRC.md). Its theorem and original four certified core families receive precedence. Applying its conductor to our full core bank gives:

**L4 (78 all-pair families).** For every eleven-element C⊂{1,…,13} and all positive integers u,v, there is a time t with
\[
\|st\|\ge1/14\qquad(s\in C\cup\{u,v\}).
\]
This subsumes the 2,184 fixed-ratio family assertions above. The core {1,…,11} was already covered by the incoming theorem; the present extension checks all other 77 cores too. It is not an unrestricted thirteen-speed assertion.

For completeness, the incoming stopping argument at r≤1/6 is short. If a core has a width-W interval, put B=ceil(1/W). For B≤u≤v, that interval contains a full u-period. Its u-safe part has at most two pieces and total length (1−2r)/u. One piece therefore has length at least (1−2r)/(2u)≥2r/u≥2r/v, enough to reach the closed v-safe set from any starting phase. Thus both-large speeds always work.

For each u<B, compile the complete joint cover for C∪{u}. If its largest width is W_u>0, every v≥ceil(2r/W_u) works by the one-tail gap bound. It remains only to check u≤v<ceil(2r/W_u). Every branch in our 78-core certificate has W_u>0. All 6,397 remaining pairs pass; the largest both-large cutoff is160. all_pair_cores.py checks every finite pair both by a closed-tooth inequality and by a literal comb intersection, and validates its rational witness against all original speeds. Exact individual branch cutoffs are saved in all_pair_receipt.json; normal and -O receipts agree.

We do not reimplement or claim the incoming Euclidean floor-sum algorithm. The finite bank check here uses a simpler independent interval path; the incoming compiler supplies an efficient way to answer new large pair queries without enumerating either attached comb.

Hostile cases delimit the theorem:

- For C={1,2}, A={3,12}, r=1/4, the complete rejected set is {1}; every d≥2 succeeds. At d=1 each individual bank runner can be attached, but they have no common feasible time.
- For C=A={1,2}, r=1/3, the core consists of isolated points and rejection is exactly d≡0 mod3. Deleting point components would destroy the true periodic classification.

## The correct clique vertices

For each additional speed v, take the closed component intervals of \(G_C(r)\cap G_{\{v\}}(r)\), retaining color v and their coordinates. Merge touching intervals of the same color. Join two graph vertices iff their intervals intersect.

**L3.** A clique containing one component of every speed exists iff the full runner system has a common feasible time. For intervals, pairwise intersection implies \(\max L\le\min R\), proving the one-dimensional Helly property directly. Distinct components of one color are disjoint, so every clique uses at most one from that color. Therefore the clique number equals the number of added speeds iff the system is feasible.

Closed endpoints and point intervals are essential. An event implementation must count starts before deleting ends at an equal coordinate.

For speeds {1,2,3} at r=1/3, the graph connecting whole runners whose safe sets intersect is K3, yet all three safe sets have empty intersection. The literal component graph has six vertices and clique number two. The checker also tests all 35 three-speed subsets of {1,…,7} at four radii.

This is a faithful exact reformulation, not a claim that pairwise runner information proves LRC. THM-658's Goddyn–Wong example also warns that several graph-coloring surrogates can disagree with the rotation threshold. Its displayed chain has three distinct levels among four quantities; this session corrects the old wording “four distinct,” retaining its mathematical conclusion.

**Next decisive test.** Find a structural class of thirteen-speed inputs whose complete phase components can be generated compactly and whose native bank gap/width ratio is uniformly controlled. L1 is exact once those covers are available; it does not prove that arbitrary inputs have nonempty core or bank covers.
