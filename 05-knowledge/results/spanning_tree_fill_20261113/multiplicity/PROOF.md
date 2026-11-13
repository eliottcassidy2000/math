# Exact growth forests, permanent ports, and the support obligation

Status: **PROVED** collision-controlled growth and marked-count consequences
below; **INHERITED** SCC factorization, fixed-H polynomiality, the ready inverse,
and the prime/7p reduction; **FINITE-EXACT** the declared inherited-bank replay.
Chan–Pak remains **OPEN**. No literature-priority claim is made.

## 1. What is inherited, and what is being added

Let H(T) count actual directed Hamilton paths. Throughout, ordinary
isomorphism classes are counted, without identifying converse. The known
fixed-H result already answers the basic polynomial-growth question:

- [`actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md),
  equations(1),(4),(5), proves the rational generating function and eventual
  polynomial degree, with the 7-adic deficit.
- [`spectrum_multiplicity_signed_tables_20261111/multiplicity/PROOF.md`](../../spectrum_multiplicity_signed_tables_20261111/multiplicity/PROOF.md)
  adds the exact labelled/automorphism-weighted and cumulative conventions.
- The less-used
  [`session153_gadgets_20261103/species/SPECIES_TROPICAL.md`](../../session153_gadgets_20261103/species/SPECIES_TROPICAL.md)
  already gives the formal Dirichlet identity and minimum-order convolution.
- [`THM-4094, hamiltonian-matching-deficit-and-two-prime-lane-completeness`](../../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md),
  Sections5–6, already isolates ordinary primes and7p, with the finite carry
  atoms49,63,343. This arithmetic reduction is not rediscovered here.

Our additions make the inherited count formula into an actual one-vertex
forest, retain permanent named ports, and quantify how total rooted
multiplicity can diverge from useful ready-source multiplicity.

Concept board: actual-H fibres / ordered SCC words / gap-lattice growth /
permanent marks / target support versus repeated occurrences. The hostile
is H9: many ordinary and rooted objects coexist with only one unit-ready
rooted class at every sufficiently large order. A second hostile is H39:
its3·13 component words must be counted in addition to its strong atoms.

## 2. The inherited arithmetic compiler, with its exact scope

Write W=(C₁,…,Cᵣ) for an ordered word of non-singleton strong classes. Set
`m(W)=Σ|C_i|`, `h(W)=∏H(C_i)`, `a(W)=∏|Aut(C_i)|`.
For fixed h, the set W_h of such words with product h is finite: Camion
implies `3≤|C_i|≤H(C_i)`, so both factors and their possible graph orders are
bounded. The empty word occurs only at h=1.

Each actual class has a unique description

\[
TT_{g_0}\to C_1\to TT_{g_1}\to\cdots\to C_r\to TT_{g_r},
\qquad g_i\ge0.
\tag{1}
\]

All arcs between successive SCCs point forward. Every Hamilton path traverses
these SCCs in that order, so H is exactly h(W), and the automorphism group is
the direct product of those of the C_i. Identical components do not admit an
extra permutation factor: their positions in the condensation order are fixed.

With t=n−m(W), the number of gap vectors is `binom(t+r,r)` for t≥0. Thus the
inherited exact formulas are

\[
N_n(h)=\sum_{W\in W_h}\mathbf1_{t\ge0}\binom{t+r}{r},\qquad
\frac{L_n(h)}{n!}=\sum_{W\in W_h}\frac{\mathbf1_{t\ge0}}{a(W)}\binom{t+r}{r}.
\tag{2}
\]

Here L counts tournaments on the fixed labelled set[n]. If h is attained
and `r_h=max |W|`, N_n(h) eventually has degree r_h. Its leading coefficient
is `b_h/r_h!`, where b_h counts top-length words. This is inherited, not a new
asymptotic assertion. Absent h has zero polynomial, not degree0.

Equivalently, for `A_d(z)=Σ_C z^{|C|}` over nontrivial strong classes of
count d, the ordinary formal Dirichlet series is

\[
\sum_{h\ge1}\left(\sum_nN_n(h)z^n\right)h^{-s}
=\frac1{1-z-\sum_{d>1}A_d(z)d^{-s}}.
\tag{3}
\]

This is a coefficientwise formal identity, not analytic convergence.
A finite strong bank computes exactly its own SCC-closed language; it
computes the unrestricted h-fibre only when all strong divisors of h are
represented completely. Missing strong composite atoms can change both
minimum size and multiplicity.

## 3. A collision-free physical one-vertex forest

Fix W and orient every new singleton according to(1). For a nonzero gap
vector g, define its parent by subtracting1 from its **first nonempty gap**.
This deletes an actual singleton SCC, preserves H, and keeps every core
vertex and every surviving arc. The only root is g=0.

Conversely, if j is the first nonempty gap of g, the canonical children are
`g+e_i` for0≤i≤j. At g=0, permit every0≤i≤r. A child has exactly the declared
parent, because its first nonempty gap is precisely the selected insertion
gap. Therefore this is an actual spanning tree of all padding classes of W,
with one new vertex per edge and no isomorphism collisions.

Unrestricted insertion into all r+1 gaps is a different graph. A child g has
exactly `#{i:g_i>0}` parents there. Consequently at depth t+1,

\[
(r+1)\binom{t+r}{r}
=\sum_{|g|=t+1}\#\{i:g_i>0\}.
\tag{4}
\]

Counting every insertion history as a fresh class would overcount. The
canonical parent removes this collision exactly, not by an approximate
multiplicity correction.

For t≥1, a node has more than one canonical child exactly when g₀=0. Hence

\[
\#\{\text{branching nodes at depth }t\}=\binom{t+r-1}{r-1},\qquad
\frac{\#\text{branching}}{\#\text{nodes}}=\frac r{t+r}.
\tag{5}
\]

For r=0 the count is0. The mean outdegree is
`(t+r+1)/(t+1)=1+r/(t+1)`. Thus this exact fixed-core growth tree becomes
asymptotically unary; polynomial growth comes from sparse branching. Summing
over W_h gives branching fraction `r_h/n+O_h(n^{-2})` for h>1. This is a
statement about this specified canonical forest, not every possible
one-vertex augmentation graph.

## 4. Permanent ordered ports: exact positive formula and lawful growth

Fix k distinct named ports, counted up to isomorphisms preserving each name.
Let `f_W(g)` be the number of fixed core vertices of an automorphism g of W,
and put

\[
O_j(W)=\frac1{a(W)}\sum_{g\in\operatorname{Aut}(W)}(f_W(g))_j,
\]

where `(x)_j=x(x−1)…(x−j+1)`. This counts ordered distinct j-tuples of core
vertices modulo automorphisms. The fixed-point distribution, not merely the
group order, is the required sidecar.

For one padding vector with t singleton vertices, every singleton is fixed
by every automorphism. Burnside therefore gives exactly

\[
\frac1{a(W)}\sum_g(t+f_W(g))_k
\tag{6}
\]

port choices. Multiplying by the number of gap vectors yields the k-marked
ordinary count. There is a more constructive positive expansion:

\[
\boxed{
P_{n,k}(h)=\sum_{W\in W_h}\sum_{j=0}^k
\binom{k}{j}\frac{(r+j)!}{r!}\,O_{k-j}(W)
\mathbf1_{n\ge m+j}\binom{n-m+r}{r+j}.}
\tag{7}
\]

**Proof and physical interface.** Choose j names to live on singleton SCCs.
Place the remaining names on distinct core vertices, with O_{k-j} orbits.
Interleave the j distinctly named singleton anchors with the r ordered core
components: there are `(r+j)!/r!` choices. These r+j immutable anchors leave
r+j+1 gaps for **unmarked** singleton SCCs. Distribute n−m−j vertices among
them. This proves(7). On each anchor pattern, the parent rule in Section3
deletes only an unmarked singleton. Every port is retained as the same
physical vertex at every edge. A canonical representative for each finite
core-mark orbit removes the remaining finite symmetry.

Equivalence to(6) follows from the falling-factorial Vandermonde identity and

\[
\binom{t+r}{r}(t)_j=\frac{(r+j)!}{r!}\binom{t+r}{r+j}.
\]

For any attained h and fixed k, the eventual degree is exactly r_h+k and the
leading coefficient is `b_h/r_h!`. The leading terms have all k ports on
singleton SCCs. The fraction having any port inside a nontrivial SCC is
`O_h,k(1/n)`. Thus the large population of port choices is mostly padding
positions, not new strong-core access.

The formula concerns ordered distinct marks. Unordered marked subsets need
cycle-index data and cannot be substituted into(6).

## 5. Unit-ready occurrences stabilize instead of growing

For an indegree-one root r with sole predecessor p, put
`δ=H(T−r)−s_T(r)`. The inherited birth identities identify this with `F−H`,
where F counts ordered covers `(P,Q)` with P nonempty ending at r and Q possibly empty.
Here “unit-ready” means δ=1. Let U_n(h) count ordinary rooted classes of this
type, and R_h(z) count **strong** unit-ready rooted classes by order.

For h>1 the inherited ready inverse gives

\[
\boxed{\sum_n U_n(h)z^n=\frac{R_h(z)}{1-z}.}
\tag{8}
\]

For completeness, the SCC reason is short. If r is a singleton SCC, its
indegree-one condition leaves exactly one vertex before it. It has no full
starting path, and `δ=H(T)=h`, so it cannot be unit-ready at h>1. Otherwise
its nontrivial strong component already supplies an in-neighbor, forcing it
to be the first SCC, with no preceding singleton. The exact product formula
then gives

\[
\delta_T=\delta_C\,H(\text{remainder}).
\]

Equality1 forces δ_C=1 and a transitive remainder. Conversely append any
transitive tail to such a strong ready source. This preserves the root,
predecessor, exact H and defect. The unique SCC decomposition prevents
collisions. This is the same inherited classification used in
[`two_entry_funnels_20261103/READY.md`](../../two_entry_funnels_20261103/READY.md)
and later ready-inverse packets; (8) records its counting consequence.

R_h is a finite polynomial, so U_n(h) becomes an exact constant. In contrast,
all ordinary rooted occurrences P_{n,1}(h) grow with degree r_h+1. If any
unit-ready source exists, its proportion is `Θ_h(n^{−r_h−1})`; otherwise it
is zero at every order. For h=1, separately, the second vertex of the
transitive tournament is its unique indegree-one unit-ready root for n≥2,
so the series is `z²/(1−z)`.

A concrete all-order example is h=9. The inherited global strong inventory
has two rigid order5 atoms plus the word(C₃,C₃), giving

\[
N_n(9)=\frac{(n-4)(n-1)}2,\qquad
P_{n,1}(9)=\frac{(n-4)(n^2-5n+20)}2\quad(n\ge5),
\]

while `U_n(9)=1` for every n≥5. At n=14 this is65 ordinary,730 rooted,
7870 ordered-two-port classes, but one unit-ready rooted class.
The all-order values U_n(15), U_n(25), U_n(27) are identically zero despite
nonempty and growing ordinary fibres. This is a sharp obstruction to
interpreting multiplicity growth as a ready-source supply theorem.

## 6. Fair scheduling and what a coverage certificate still owes

Given a finite complete strong divisor inventory for h, the finite forest
roots are exactly W_h. Breadth-first enumeration in total vertex order,
using the canonical rule, lists every ordinary class exactly once and never
starves a padding branch. Equation(7) gives the same schedule while retaining
any fixed number of named ports. With an enumerated infinite strong bank,
diagonal scheduling over discovery stage, total core size and padding depth
visits every finite constructed object. Fairness says nothing about labels
not present in the bank's support.

For any finite target set E, one actual witness per h∈E suffices to put all
of E at a common vertex order: pad each to at least the largest witness
order. Conversely any row with every h∈E supplies those witnesses. Thus
perfect fill on a finite window has the exact missing obligation

\[
\text{one actual witness for every allowed target in the window}.
\]

After that obligation is discharged, the forest supplies perpetual growth,
exact multiplicities for a complete bank, and common-order scheduling. It
does not discharge the obligation itself. In particular, closing a spanning
tree on object addresses is different from covering every odd integer label.

For a bank used only by ordinal concatenation and singleton insertion,
its target support is its multiplicative semigroup. Any finite bank has only
finitely many prime divisors among its counts and therefore misses infinitely
many odd prime targets. This obstruction does not apply to changing-core or
selective-port operations, which can create new count primes.

The prior prime/7p reduction pinpoints the infinite witness obligations; it
is not a finite counting identity. At a **fixed vertex budget**, missing
composite values also matter: e.g. at order4,3 and5 are realized but9 is not.
A factorization3·3 costs six vertices, while the strong9 atom first appears
at order5. Thus global multiplicative support and minimum-order fill must
not be conflated.

One further inherited hostile makes the conventions concrete. Fixing
`H/|Aut|=1` allows both K₁ and C₃ as unit letters, with generating function
`1/(1−z−z³)` for their ordinal-word sublanguage. Fixing actual H=1 leaves only
K₁ and one class per order. Exponential growth in the normalized fibre is
not growth of one actual-H fibre.

## 7. Reproduction and finite universe

```
python3 05-knowledge/results/spanning_tree_fill_20261113/multiplicity/check.py
```

The script reuses the inherited ordinary/converse-separated strong atlas for
H≤27. The cited Busch bound gives strong order≤7 in this range, so this
finite bank is globally complete for the declared targets. No new tournament
census is run. It checks:

- 35 gap-tree layers, including unrestricted collision multiplicities;
- 963 exact anchor/Burnside identities against the inherited fixed-H counts;
- 30 direct ordered-mark orbit computations on inherited classes through5.

`growth_table_h_le27.csv` gives exact ordinary, labelled, rooted, two-ordered-
port and unit-ready-root counts at every n=1,…,14 and every odd h≤27, including
zero cells. `receipt.json` preserves the strong atoms and fixed-point
sidecars. The proof of the growth and interface laws is independent of this
small replay.
