# Refinement debt under paired native resets

2026-10-16. **PROVED:** a sharp general bound for the row-envelope
refinement cost, a uniformly bounded Gram remainder that contracts under
native extensions, and the conditional transport rules below.
**FINITE-EXACT:** a new individual-row test at `X=2^22` passes all fifteen
paired windows of the previously fixed controller.
**OPEN:** uniform paired half-energy, uniform paired `7X` residual,
independent all-scale weighted overlap, unconditional amplitude
nonvanishing, and universal ROOT entry. These are distinct obligations.
Self-audited; no independent referee is claimed.

## Inheritance and the active board

The nearest proved mechanism is the original rowwise survivor allocation
in [native packets](collatz_survivor_packets_20261015e.md). The actual hostile
is the one-reset `7X` failure at `k=23` in
[reset debt](collatz_reset_debt_20261015f.md). The corrected near miss is
the angular term lost when rowwise contraction was replaced by a global
core bound. The underused sidecar is the signed terminal vector of each
physical source row before taking its norm.

The [preceding paired-energy note](collatz_paired_energy_20261015g.md)
separates scalar norm variance from within-bin sign mixture. The
[earlier annular theorem A1--A4](collatz_annular_overlap_20261015c.md)
already pays the **scalar** refinement cost by `4X(b^2-1)/3` and warns that
an entrywise-absolute Gram score need not increase under a bin split.
The result here preserves those facts: it does not claim a new scalar
energy constant. It lifts the refinement account to a matrix that can
be transported through every subsequent native map.

Anchor: paired weighted cancellation. Niche: component-wise six-colour
palette compatibility. Wildcard: inverse clocks and projected operator
defects. Board: original source coefficients; literal integer clocks;
shared terminal labels; physical-bin refinement; matched cancellation;
boundary component orientations. The connection contract is in
[component palette matching](six_colour_component_pairs_20261016a.md).
The META-PATTERNS cards used are “Retain the product and shared labels
before taking scalar totals” and “Audit and close sections under their
next native operation.” No new general meta-pattern is claimed.

The live remote was fetched at startup and after the first computation.
It had no incoming changes beyond this session's own checkpoint. The
inheritance pass recovered
[THM-2238, terminal-cover-moment-hinge-gram-and-fibrewise-convex-hierarchy](../../01-canon/theorems/THM-2238-terminal-cover-moment-hinge-gram-and-fibrewise-convex-hierarchy.md)'s
distinction between a global Gram and fibrewise capacity data, and the stronger prior scalar refinement
theorem, rather than relabel their conclusions as new.

## 1. The exact objects and preserved clock

Use the shortcut map `T(n)=n/2` for even n and `(3n+1)/2` for odd n,
stopped at ROOT 1. Let `w(n)=3[n=2 mod3]-1`. Fix `X` and every original
source through `8X`. One authenticated common observed graph includes all
initial halving edges. Its signed terminal column is
`phi(n)=(-1)^d(n) e_(endpoint(n))`, retaining the full integer `d(n)`.
No convergence oracle or independent choice of terminal sign is used.

For `X<=N<2X`, define

\[
 v_N=\sum_{N<n\le4N}w(n)\phi(n).
\]

The initial halving edges imply

\[
 v_{N+1}-v_N=\sum_{a=1}^3 w(4N+a)\phi(4N+a),\qquad
 \|v_{N+1}-v_N\|_1\le4.                              \tag{1}
\]

The fourth entering term cancels the departing one: their weights agree,
their terminal is the same, and their literal clocks differ by two.
Every common signed native extension is an l1 contraction on these vectors,
so (1) survives all such extensions. Witnessed zero-deletion can also be
allowed; the computations below use no nonROOT cycle or deletion.

For a physical bin `I` of length `b`, write
`s_I=||sum_I v_N||_1`, `q_I=sum_I ||v_N||_1`. Let

\[
 E_f=\sum_N\|v_N\|_1^2,\qquad E_c=\sum_I s_I^2/b,
 \qquad Q_c=\sqrt{E_c}.
\]

The original source weights are unchanged. The factor `1/b` is the
original physical averaging measure, not a reweighted terminal population.

## 2. An independent bound for the entire envelope refinement cost

**R1 (PROVED).** More generally let consecutive vectors in any normed
space differ in norm by at most `L`. For uniform bins of length `b`, put

\[
 m_b=\frac{b^2-1}{3b},\qquad
 t_b=\frac{7b^4-15b^2+8}{60b^2}.
\]

Then

\[
 E_f\le E_c+2L m_b\sum_I s_I+L^2X t_b
     \le E_c+2L m_b\sqrt X Q_c+L^2X t_b.              \tag{2}
\]

For the native rows, `L=4` and `b<=sqrt(X)` give the scale-uniform bound

\[
 \boxed{E_f-E_c\le\frac83 XQ_c+\frac{28}{15}X^2.}     \tag{3}
\]

The exact finite-b formula is sharper than (3). The same relaxed bound
holds for unequal bins of maximum length at most `sqrt(X)`.

**Proof.** Number a bin's rows `i=1,...,b`, and let `bar v` be their mean.
With

\[
 a_i=b^{-1}\sum_{j=1}^b|i-j|,
\]

the triangle inequality gives `||v_i||<=s_I/b+L a_i`. Direct summation
gives `sum_i a_i=b m_b` and `sum_i a_i^2=b t_b`. Squaring and summing
proves the first bound in (2). Cauchy--Schwarz gives
`sum_I s_I<=sqrt(X) Q_c`. For unequal bins, use
`m_b<=b_max/3`, `t_b<=7b_max^2/60` and the same weighted Cauchy inequality.
This proves (3). QED.

The first inequality in (2) is sharp for unrestricted l1-valued paths,
even with a prescribed nonzero mean. On one bin take a constant coordinate
`a>=0` and one additional coordinate for each cut `j=1,...,b-1`:

\[
 v_i=\left(a,\ (L(1_{i>j}-(b-j)/b))_{j=1}^{b-1}\right).
\]

Each step changes exactly one coordinate by L. All cut coordinates have
mean zero and `||v_i||=a+L a_i`, giving equality. These rational test paths
are not asserted to be actual annular rows. In particular, sparse increments
alone do not remove the refinement price.

This closes a deterministic pricing question, not the reset problem.
For example, **if** a native controller establishes `Q_c<=7X`, (3) gives

\[
 E_f\le\frac{1043}{15}X^2,
 \qquad E_f-E_c\le\frac{308}{15}X^2.                  \tag{4}
\]

The estimate needs the coarse energy input. A large common constant
coordinate has zero increment and arbitrarily large mean; no increment
bound by itself can prove `Q_c<=7X`.

## 3. A better carrier for refinement is a Gram remainder

Define matrices on the **same** terminal coordinates,

\[
 G_f=\sum_N v_Nv_N^T,\quad
 G_c=\sum_I\frac{(\sum_Iv_N)(\sum_Iv_N)^T}{b},\quad K=G_f-G_c.
\]

**R2 (PROVED).** The matrix K is positive semidefinite and has the exact
pairwise representation

\[
 K=\sum_I\frac1b\sum_{i<j\in I}(v_i-v_j)(v_i-v_j)^T.
                                                               \tag{5}
\]

Consequently, with entrywise absolute norm,

\[
 \boxed{\|K\|_{\rm ent}\le\frac43X(b^2-1).}           \tag{6}
\]

Under a further common signed native map A, the new remainder is
`K'=A^T K A`, and

\[
 \|K'\|_{\rm ent}\le\|K\|_{\rm ent}.                 \tag{7}
\]

**Proof.** Expand (5). For each outer product its entrywise norm is
`||v_i-v_j||_1^2<=16(j-i)^2`. Summing
`sum_(i<j)(j-i)^2=b^2(b^2-1)/12` proves (6). Each old terminal has one
signed target column (or an authenticated zero), so triangle inequality
in every ordered target-pair fibre proves (7). Averaging the physical
rows commutes with the common linear map, giving the stated congruence.
QED.

In particular

\[
 |\|G_f\|_{\rm ent}-\|G_c\|_{\rm ent}|\le\|K\|_{\rm ent}.
                                                               \tag{8}
\]

No monotonicity of the entrywise Gram norm under physical refinement is
claimed. For each one shared label vector z, however,
`z^T G_f z-z^T G_c z=z^T Kz>=0`; its maximum over `|z_e|<=1` is at most
the right side of (6). This recovers the earlier scalar allowance while
also retaining a matrix that can be carried through later controllers.

### Why the envelope invoice cannot be treated the same way

Use the actual native test columns `3,5`, which coalesce at 5 via clocks
`1,0`. Take physical rows `(1,0)` and `(0,1)`. Before merging, the fine
and coarse envelope energies both equal 2. Afterwards they equal 2 and 0.
Thus the envelope refinement gap grows from 0 to 2. Its norm-variance
part is zero throughout, so this is entirely sign mixture.

In contrast K has entrywise norm 2 both before and after the merger.
This is a test of a genuine native map on abstract rows, not a counterexample
inside the annular source family. The missing coordinate is the shared
terminal sign: taking separate row maxima destroys the linear commuting
square used in (7).

## 4. A source-preserving bank of refinement rows

For dyadic bins, let `h_J=sum_(left half J)v_N-sum_(right half J)v_N`
for every internal dyadic subinterval J. Then

\[
 K=\sum_J h_Jh_J^T/|J|,
 \qquad R=\sum_J\|h_J\|_1^2/|J|
       \le\frac43X(b^2-1).                             \tag{9}
\]

The identity is ordinary orthogonal mean/detail decomposition, with all
source labels retained. Its l1 budget follows either by bounding each
half-mean difference by the four-Lipschitz ramp or by telescoping the
scalar variance allowance down the dyadic tree. Every individual h_J
contracts under native extension; hence R is nonincreasing as well.

This is a useful reroute: keep the coarse mean rows as the unpaid part,
and keep the detail rows as a uniformly priced bank. If separately chosen
finite native controllers pay different rows, their observed-edge union,
followed by common closure and repushing every row, preserves each payment.
That last fact is inherited from the
[finite-address controller lattice](collatz_finite_address_lattice_20261015c.md);
it provides no uniform bound on how large the union must be.

There is no cumulative refinement surcharge at every reset. For a paired
recurrence `C_j<=C_(j-1)/2+a_j` on the coarse score, the certificate
`C_j+R_j` obeys

\[
 C_m+R_m\le2^{-m}C_0+\sum_{j=1}^m2^{-(m-j)}a_j+R_0. \tag{10}
\]

Thus the **entire** continuing refinement account is at most the original
`4X(b^2-1)/3`, even with changing native generators. The independent
bound on the coarse discounted debt remains open. Equation (10) does not
turn a small refinement term into a mean-current bound.

## 5. The tight paired window survives individual-row refinement

The frozen controller is unchanged: expose one first reset from each
terminal, close the observed graph, and pair two already closed stages.
Use exactly the stopping stage selected previously by the default bins;
the final singleton is followed by identity. A reset remains a variable
number of literal T steps.

At `k=22`, all fifteen individual-row comparisons pass

\[
 E_{j+2}\le\max(E_j/2,49X^2).
\]

Thirteen end above the floor and are nontrivial. The previously tight
window `12 -> 14` has the new exact point-energy ratio

\[
 \frac{867888579386948}{1749166918803912}
 =\frac{216972144846737}{437291729700978}
 =0.4961725322\ldots<\frac12.                         \tag{11}
\]

At those two stages the sign-mixture terms, divided by X squared, are
approximately `0.00307666759` and `0.00133982137`. The final point energy
is approximately `36.74542648 X^2`. This is a held-out **refinement** test
at a previously tested source scale, not an untested-scale prediction or
a uniform theorem. It is evidence that the old tight pair was not created
by coarse-bin sign mixing in this example.

The streaming event sweep uses all `8X=33,554,432` original sources.
Every default row is reproduced from the frozen run; every physical
increment checks the exact halving cancellation and four-unit bound.
The native graph audits all internal edges and replays a
fixed source grid with full clocks. Small controls compare the event sweep
against the earlier full point-matrix implementation. Normal and optimized
Python executions agree on the small controls; the large k22 file records
one full exact run, not a second independent large-scale implementation.

## 6. Reproduction and the remaining arithmetic target

    python3 04-computation/experiments/collatz_point_energy_20261016a.py --k 8 14 22
    python3 -O 04-computation/experiments/collatz_point_energy_20261016a.py --k 8 14
    python3 04-computation/experiments/refinement_gram_controls_20261016a.py
    python3 -O 04-computation/experiments/refinement_gram_controls_20261016a.py

The point sweep compiles its small C++ event kernel with `clang++` into
a temporary directory. It uses no ROOT oracle. Outputs:
[k8](collatz_point_energy_20261016a_k8.json),
[k14](collatz_point_energy_20261016a_k14.json),
[k22](collatz_point_energy_20261016a_k22.json), and
[algebra/Gram controls](refinement_gram_controls_20261016a.json).
The last file checks 11,718 scalar paths, 72 sharp l1 paths, exact Haar
outer products and native congruence at `k=8,10`. The k22 SHA256 is
`7b81608ef3c240e221524ec2ae3ea85c7979cd5260e3846ee20dc353cbd2dd9d`.

For a paired map let `u_i,v_i` be the original and terminal row norms,
and `M_i=(u_i-v_i)/2` the actually matchable opposite-sign mass. The
unpaid source-specific statement is still the earlier bound

\[
 \sum_i\frac{M_i(u_i+v_i)}{b_i}
 \ge\frac12\min\{E/2,E-49X^2\},\qquad E>49X^2.        \tag{12}
\]

The strongest new separation is operational: refinement now has an
independent, nonaccumulating matrix budget, while the arithmetic mean
still needs source-weighted matching. A proof of (12), or a sufficient
discounted substitute, must use more than the four-Lipschitz fact.
The next useful helper is a family of disjoint, authenticated opposite-sign
source pairs whose weighted capacity gives that lower bound; a count of
terminal collisions alone cannot do so.
