# Good-cut recovery: SCC run coefficients and a generating-function repair

**CORRECTION + PROVED RECOVERY:** the run-composition formula is already in canon; the printed bivariate generating function in THM-341 is missing a factor1/(1−x). **PROVED COROLLARY:** fixed-good-cut conditioned morphology and positive finite differences. **FINITE-EXACT:** every tiling through n=8, independent SCC checks through n=7, formal coefficient checks through n=20. No new theorem ID or Git mutation.

## 1. Inheritance and exact attribution

Recovered sources, all cited by ID plus path:

- THM-336, `01-canon/theorems/THM-336-good-cuts-structure.md`: small good-cut buckets.
- THM-340, `01-canon/theorems/THM-340-SC-composition-formula.md`: ordered block convolution Q(k,r)=[u^k]A(u)^r, including low-excess diagonals.
- THM-341, `01-canon/theorems/THM-341-bivariate-GF-good-cuts.md`: its scalar formula has the error repaired below; its displayed intermediate sum retains the missing factor and identifies the failure exactly.
- THM-349, `01-canon/theorems/THM-349-good-cut-interval-union-recurrence.md`: actual good-run factorization and interval-union recurrence.
- THM-354, `01-canon/theorems/THM-354-good-cut-scc-count.md`: g(T)=n−#SCC(T), independently of Hamiltonian path.

The closest mechanism is the SCC/run factorization, not an analogy to geometric dimensions. The canonical hostile is n=4, g=2: its two possible adjacent good-cut pairs require two presentations, while the printed THM-341 GF counts only one. The underused sidecar is the number r of nontrivial SCCs, which determines the polynomial degree. The live board is good-support geometry / SCC sequence / fixed-path presentation weights / generating functions / conditioning / black path-orbit deficit.

There is also an indexing/numeric correction to the incoming proposed sequence: the strong fixed-path count at six vertices is903, not874. No inherited canon source inspected here supports874; THM-340 and THM-349 both correctly use903.

## 2. The exact formula, recovered rather than newly claimed

Let s_j count strong tournaments with the chosen directed Hamiltonian path fixed, on j vertices. Then

s_1=1, s_2=0, s_3=1, s_4=5, s_5=50, s_6=903, s_7=30773.

Define A(u)=sum_{ell>=2} s_(ell+1)u^ell. A good run of ell cuts is a nontrivial SCC of ell+1 vertices. No upward tile may cross a bad cut; choices inside different runs are independent. Thus a support with ordered run lengths ell_1,...,ell_r has product s_(ell_i+1) realizations, and no run may have length1.

For n vertices and k good cuts, fix those ordered lengths summing to k. There are n−1−k bad cuts. The r−1 internal gaps must each contain at least one bad cut, while the two endpoint gaps may be empty. Subtracting the mandatory gaps and distributing the remainder among r+1 slots gives binom(n−k,r) placements. Therefore

F_n(k)=sum_{r=0}^{floor(k/2)} Q(k,r)binom(n−k,r),
Q(k,r)=[u^k]A(u)^r,

with Q(0,0)=1 and Q(k,0)=0 for k>0. Here n>=k+1; ordinary binomial coefficients vanish when their upper index is smaller than r. F_n(1)=0. This formula is the explicit placement synthesis of THM-340 and THM-349, not a newly discovered counting law.

For example:

F_n(5)=903(n−5)+10 binom(n−5,2),
F_n(6)=30773(n−6)+125 binom(n−6,2)+binom(n−6,3),
F_n(7)=2032504(n−7)+2306 binom(n−7,2)+15 binom(n−7,3).

At n=10, F_10(6)=123846. These are raw fixed-path presentations, not tournament isomorphism classes.

## 3. Exact repair to THM-341

Let T(x,u)=sum_{N>=0}sum_k F_(N+1)(k)x^N u^k count by the number N=n−1 of cuts. The correct formal identity is

T(x,u)=(1+A(xu))/(1−x(1+A(xu))).

One proof is the run recurrence of THM-349. Another uses the ordered SCC sequence: a singleton block contributes x, and a nontrivial SCC contributes x*A(xu), so the vertex-count GF including the empty tournament is

W(x,u)=1/(1−x−xA(xu)).

Retaining the nontrivial-component count r gives the stronger formal marker W(x,u,t)=1/[1−x−t*x*A(xu)]. Its coefficient [x^n u^k t^r] is Q(k,r)binom(n−k,r): there are exactly n−k total SCCs, and one chooses which r of their ordered positions are nontrivial. This supplies a second direct explanation of the placement factor.

Since W(x,u,1)=1+xT, the displayed cut-count identity follows. The vertex-count GF restricted to positive good-cut count is

Fpositive(x,u)=W(x,u)−1/(1−x)
             =x*A(xu)/[(1−x)(1−x−x*A(xu))].

THM-341 prints the same expression WITHOUT the first denominator factor1−x. Its proof goes from

(1/(1−x))*[x*A(xu)/(1−x)]/[1−x*A(xu)/(1−x)]

to x*A(xu)/(1−x−x*A(xu)), dropping the prefactor1/(1−x). The correct simplification retains it.

**Minimal witness:** [x^4 u^2] of the printed expression is1, whereas both THM-336 and direct enumeration give2. At n=3,k=2 both are1, so n=4 is the first mismatch. This refutes the displayed GF and its associated claimed brute-force verification, not THM-340's composition formula, THM-349's recurrence, or THM-354's SCC identity. THM-341's C1 specialization must use the repaired denominator too. Its C4/C5 bucket statements survive but do not follow from the erroneous printed formula.

Strongest survivor of the erroneous expression: x*A(xu)/(1−x−x*A(xu)) correctly counts tilings whose FIRST SCC is nontrivial. It forces that first block and then appends an arbitrary SCC sequence. The missing geometric factor restores the arbitrary leading singleton prefix, so the error has a precise combinatorial meaning rather than an indexing repair.

The repaired scalar formula, the SCC-sequence formula and the cut-word formula are equivalent. This is a formal-power-series identity: A has superexponentially growing coefficients and radius of convergence zero. Explicitly, for j>=3, forcing the closing arc of the fixed Hamiltonian path creates a directed Hamiltonian cycle and leaves all other nonconsecutive arcs free, giving s_j>=2^(binom(j−1,2)−1). There is no claim of an ordinary analytic generating function near positive x. For fixed k, its individual x-column is rational and has denominator (1−x)^(floor(k/2)+1), after cancellation only when the column vanishes.

## 4. Polynomial degree and conditioned morphology

For fixed k>=2, F_n(k) is a polynomial in n of degree floor(k/2), on its natural range n>=k+1. There can be at most floor(k/2) nontrivial SCCs, and their weights are positive.

For k=2r, the top term Q(2r,r)=1 comes from r runs of length2, so the leading coefficient is1/r!. For k=2r+1, the top term Q(2r+1,r)=5r comes from one run of length3 and r−1 runs of length2, so the leading coefficient is5/(r−1)!. This recovers the diagonal mechanism already present in THM-340, rather than inferring the degree from a fit.

All forward finite differences in n of orders j<=floor(k/2) are positive for sufficiently large n; those above this degree vanish. Indeed Δ binom(n−k,r)=binom(n−k,r−1), and all Q(k,r) are nonnegative. This supplies an integer-valued polynomial certificate and a finite-difference diagnostic.

**Conditioned morphology corollary.** Choose uniformly from all fixed-path tilings with exactly k good cuts, with k fixed and n->infinity. The probability of the maximal possible number floor(k/2) of nontrivial SCCs tends to1, with error O_k(1/n). For k=2r the dominant structure is r three-vertex strong components, with all remaining SCCs singletons. For k=2r+1 it is one four-vertex strong component and r−1 three-vertex strong components, with the rest singleton. Among maximal-run odd-k patterns, the position of the four-vertex block and its five fixed-path strong presentations are exactly uniform over5r choices: every such ordered pattern has the same placement count.

Proof: all smaller-r terms have polynomial degree at most r−1, while the top coefficient is positive. The description follows from the minimal run length2. This is a conditional raw-presentation statement; it is not a distribution on unlabelled tournament classes or automorphism-weighted classes.

## 5. Why deficit-six is not this coefficient

`05-knowledge/results/deficit_six_centers_20261021b/README.md`, §1, defines black deficit Δ=H/|Aut(T)|−b, where b counts blue automorphism-orbits of Hamiltonian words, and its R_6 counts strong self-converse centers with Δ=6. Here k=g=n−#SCC is an entirely different invariant. In particular, on EVERY strong n-vertex source, g=n−1 regardless of its black deficit. R_6 is not the coefficient [u^6] of this run series, and903 or30773 must not be substituted into it.

The shared mechanism is ordered SCC decomposition, with different block alphabets and different weights. A lawful bridge must retain H/Aut, blue orbit count, converse constraints, and boundary attachments in addition to SCC sizes. Without those sidecars there is no exact deficit-six transfer. We do not use the older packet's superseded degree bounds as current claims; only its definition is needed here.

## 6. Reproduction and scope

`python3 session35-arithmetic/check_goodcuts.py` uses only Python standard libraries and explicit exception checks. It enumerates every subset of nonconsecutive path chords through n=8 (2^21 presentations at n=8), computes interval unions, and compares every bucket with the placement formula. Through n=7 it also constructs the actual directed graph and computes transitive closure to count SCCs independently of the cut representation. The strong counts are separately generated by first-SCC decomposition of all presentations, and total counts/formal-GF repair are checked through n=20.

Outputs: `goodcuts_receipt.json`, `goodcuts_check.out`. The receipt preserves the exact smallest GF counterexample, the corrected903, and all small buckets. Run with `python3 -O` as a control; checks do not disappear.

Research outcome: recovery prevented duplicate canon; a hostile exact coefficient exposed a real old error; the surviving SCC sequence gave the repaired GF and a precise conditioned asymptotic. No major open problem is claimed closed.

Independent audit: `session35-algebra-audit/GOODCUT-SYNTHESIS-AUDIT.md` reports PASS, independently enumerating interval subsets through n=7 and confirming the GF hostile, polynomial coefficients, conditioned morphology, and formal convergence boundary.
