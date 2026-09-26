# Gilbreath, universal edge traces, and the metric sidecar of a tournament

**PROVED elementary identities / FINITE-EXACT controls / OPEN arithmetic
selection problem, 2026-09-26.** No Collatz or Gilbreath proof, no claim of
literature priority. The exact statements below have self-contained proofs.
The scope of each finite computation is stated separately.

## 1. Inheritance and the live board

Anchor: the [excursion forest](forest_20260926_excursions.md), which keeps one
actual Collatz source and its chronological affine/carry data. Niche:
Gilbreath's absolute-difference triangle. Wildcard: the inverse binomial
transform, Fermat's sequence, and the metric meaning of a tournament path.

Closest inherited mechanism: the additive `0/2` sublattice and the absorbing
edge `(1,0/2,0/2,...)` in the
[incoming CA note](collatz_oscillation_gilbreath_20260926.md). Canonical
hostile: an actual parity cylinder has arbitrarily large sources; a finite
pattern is not one bounded integer. Corrected near miss: additivity does
**not** forbid encoding a nonadditive orbit into an additive automaton's
edge. Least-used sidecar: the complete triangular sign schedule, together
with metric gaps and the original index labels.

The board has six objects: absolute differences, binomial transforms,
binary parity, scalar tournaments, ordered metric gaps, and one-source
excursion forests. Each connection below specifies its information loss.

Classical context is **CITED**, not rediscovered: Odlyzko's 1993
[*Iterated absolute values of differences of consecutive primes*](https://www-users.cse.umn.edu/~odlyzko/doc/arch/gilbreath.conj.pdf)
studies the prime difference triangle and its heuristics. Rédei's odd
Hamiltonian-path count is treated in
[Schweser--Stiebitz--Toft, *The Tournament Theorem of Rédei revisited*](https://arxiv.org/html/2510.10659v1).
The [Grinberg--Stanley Rédei--Berge paper](https://arxiv.org/abs/2307.05569)
also derives that theorem and a modulo-four refinement; its posted draft
explicitly distinguishes outline portions. None of those deeper results
is needed for the elementary proofs below. The repo's classical route is
[THM-001-redei](../../01-canon/theorems/THM-001-redei.md).

## 2. Every binary edge trace is realizable

Work with one-sided infinite binary sequences. Put `D=1+S` over F2, where
`(Sx)_i=x_(i+1)`. If `b_r=(D^r x)_0`, then

    b_r = sum_(j=0)^r binom(r,j) x_j mod 2.

Let B denote this lower-triangular Pascal transform. It is an involution:

    (B^2)_(n,j)
      = sum_(k=j)^n binom(n,k) binom(k,j)
      = binom(n,j) 2^(n-j) = delta_(n,j) mod 2.

Thus **every binary sequence b is the edge trace of exactly one initial
binary row x=Bb**. The statement holds directly for infinite sequences
because each coordinate uses a finite sum; no limit interchange is needed.
The same proof applies to every finite prefix.

The dyadic recursion is exact:

    D^(2^m) = 1 + S^(2^m).

It gives the familiar Pascal pattern modulo two. Universality concerns an
edge observable with an unconstrained infinite seed; it is not a claim that
this additive CA simulates arbitrary local dynamics on finite inputs.

For shortcut Collatz `T(n)=(3n+1)/2` on odd n and `n/2` on even n, take
`b_r=T^r(n) mod 2`. The inverse transform gives a **source-derived causal
code**: its j-th seed coordinate is computed from the first j iterates of
that very integer. This needs no assumption of termination. It preserves
the complete parity trace, but pays for it with an infinite seed whose
construction already follows the orbit. It does not produce the primes.

A precise resource boundary is available. If a binary seed is supported
in indices `<N`, its edge is globally periodic with any power-of-two period
`2^m>=N`. For `j<2^m`, the identity
`(1+z)^(r+2^m)=(1+z)^r(1+z^(2^m))` over F2 proves that the coefficient of
`z^j` is unchanged. If a positive Collatz orbit converges, its parity trace
is eventually alternating. A globally periodic eventually alternating
trace alternates from the start. For an odd state u in such a trace,

    u_(j+1)=(3u_j+1)/4,
    u_j=1+(3/4)^j(u_0-1).

Integrality for every j forces `u_0=1`; an initial even state must be 2.
Consequently the code of a **convergent source n>2 necessarily has infinite
support**. No assertion of convergence for other sources is used here.
Conversely, an eventually alternating positive-integer parity trace does
force eventual arrival at the 1--2 cycle by the same calculation.

## 3. Stronger: every nonnegative integer trace has an absolute-difference seed

For any prescribed sequence `b_0,b_1,...` of nonnegative integers, define

    x_i = sum_(j=0)^i binom(i,j) b_j.

Let `Delta x_i=x_(i+1)-x_i`. Pascal's identity gives, by induction,

    Delta^r x_i = sum_(j=0)^i binom(i,j) b_(r+j) >= 0.       (1)

Every intermediate forward difference is nonnegative, so taking an
absolute value changes nothing. Therefore the iterated absolute-difference
triangle has left edge **exactly b_r**, including the magnitudes.

This is an exact conjugacy between the shift on nonnegative sequences and
the absolute-difference map restricted to the cone of sequences with all
forward differences nonnegative. In that cone, the seed is unique: ordinary
binomial inversion recovers `b_r=Delta^r x_0`. Outside that cone, other sign
chambers can give the same edge (§4).

Taking `b_r=T^r(n)` embeds every actual positive Collatz orbit into this
cone. If `C(n)` is the resulting infinite row, then

    abs-Delta C(n) = C(T(n)),       C(n)_0=n.

This is injective, computable coordinate by coordinate, and retains the
entire original value sequence at the edge. It still gives no simplification
of the stopping problem: the j-th coordinate requires the first j iterates,
and `C(n)_j>=2^j`. It supplies neither a finite digit input nor the initial
row of consecutive primes.

The clean Fermat connection is obtained by prescribing `b_0=2` and
`b_r=1` for every `r>=1`:

    x_i=2^i+1,            Delta^r x_i=2^i  for r>=1.

Thus `(2,3,5,9,17,33,65,...)` has the desired unit edge. Its subsequence at
indices `i=2^j` consists of the Fermat numbers, while composites already
occur at index 3. This proves that edge 1 does not characterize primes. It
also gives the unique such seed in the nonnegative-forward-difference cone.
For formal exponential generating functions the same transform is
`X(z)=exp(z)B(z)`; this algebraic identity alone says nothing about prime
values or Bernoulli denominators.

## 4. Keeping signs makes the full edge transform unimodular

For a finite real or integer row of length N, write

    x_i^(r+1)=|x_(i+1)^r-x_i^r|,
    sigma_(r,i)=+1 if x_(i+1)^r>=x_i^r, and -1 otherwise.

The choice at equality is merely an algebraic convention; it does not
turn a tie into a tournament edge. With the full triangular sign schedule
fixed, each update is an integer linear map
`diag(sigma_r) Delta`. The full left edge
`b=(x_0^0,x_0^1,...,x_0^(N-1))` has the form

    b=E_sigma x,
    (E_sigma)_(r,r)=product_(s=0)^(r-1) sigma_(s,r-1-s).

E is lower triangular. Every diagonal entry is +1 or -1, and each triangular
sign occurs in exactly one diagonal product. Hence

    det E_sigma = product_(r,i) sigma_(r,i) in {+1,-1}.      (2)

The initial row is recovered exactly over the integers by forward
substitution. Modulo two every E_sigma is the same Pascal matrix B, since
both subtraction and signs disappear.

For an arbitrary formal sign schedule equation(2) still holds, but the
reconstructed input need not satisfy the inequalities defining that
schedule. **Sign feasibility is a separate sidecar.** Each feasible chamber
is cut out by the corresponding linear inequalities. No single chamber is
being asserted to contain the prime sequence.

Discarding signs really loses information: the increasing integer rows
`(2,5,7)` and `(2,5,9)` both have edge `(2,3,1)`. The second-row comparison
is reversed. Keeping only the initial row's order is therefore insufficient
even for invertibility. The exact carrier is edge + triangular signs, not
the edge or a tournament alone.

## 5. Gilbreath's parity is automatic; the conjecture concerns magnitude

If the first initial entry is even and all subsequent entries are odd,
the left entry at every positive depth is odd, because

    sum_(j=1)^r binom(r,j) = 2^r-1 = 1 mod 2.

All other entries are even after the first row. This uses only the initial
parities and holds for vastly more rows than the consecutive primes.
Gilbreath asks for **exactly 1**, not merely an odd number. For example,
both `(2,3,5)` and `(2,3,7)` consist of primes and have identical strict
orders in each nonsingleton difference row; their bottom values are 1 and 3.
Consecutiveness and actual prime gaps cannot be discarded.

The unsigned rule has an unconditional residue quotient only modulo 1 or 2.
For any `m>2`, the input pairs `(0,1)` and `(m,1)` are coordinatewise
congruent modulo m, but their absolute differences are 1 and m-1, which
are not congruent. A deeper congruence rule must retain signs or bounds
that fix those signs. On a common `2^v` sublattice one can divide by `2^v`
and retain one further parity bit; that is not a sign-free map at arbitrary
2-adic precision. In particular, ordinary absolute value is not a
continuous 2-adic operation: `2^k-1 -> -1`, but the absolute values tend
to -1, not to `|-1|=1`.

Rédei's theorem gives an odd number of Hamiltonian paths in any tournament.
It does not give exactly one. The analogy between two parity residues cannot
supply the missing size bound in Gilbreath. Conversely, the elementary
parity calculation above does not prove Rédei's counting theorem.

## 6. A genuine Collatz hostile preserving every comparison and parity

For every integer `s>=0`, put `n=23+32s`. Its first five shortcut states are

    23+32s, 35+48s, 53+72s, 80+108s, 40+54s.

Their parity word is 11100. Put `t=1+2s`; their complete difference triangle
is

    7+16t, 11+24t, 17+36t, 26+54t, 13+27t
       4+8t,   6+12t,   9+18t,   13+27t
          2+4t,    3+6t,    4+9t
             1+2t,    1+3t
                  t.

For every `t>=1`, all row entries are distinct, the initial order of
indices is `0<1<4<2<3`, and all later rows increase from left to right.
Thus the **complete comparison tournament of every row is unchanged**.
For odd t, the whole parity triangle is also unchanged. Nonetheless the
final magnitude `t=1+2s` is unbounded. The concrete starts 23 and 55 yield
bottom values 1 and 3. These are actual orbit prefixes of actual integers,
not formal parity words or resampled paths. No difference in eventual
convergence behavior is asserted.

A general reason underlies this witness. In any fixed length-L parity
cylinder `n=r+2^L t`, each state through time L is affine in t. Any fixed
finite circuit formed from these states using affine combinations,
absolute differences, and comparisons is piecewise affine on the real
ray. Induction on the finite circuit shows finitely many breakpoints:
each absolute value introduces at most one zero per existing affine piece.
Every comparison therefore stabilizes for sufficiently large t, unless
the compared affine functions coincide, in which case it is a persistent
tie. This forbids a bound on source height from such an ordinal output
alone. It does not forbid a rule whose depth grows with the source or an
estimate retaining the affine coefficients and metric magnitudes.

## 7. What an honest tournament map preserves

For distinct scalar entries, use vertices carrying their original indices
and orient `i -> j` exactly when `x_i<x_j`. This is a **transitive**
tournament. Its unique directed Hamiltonian path is the sorted order;
Rédei adds no additional existence result for this particular carrier.
Equal entries give a comparison preorder, with genuine ties retained.

The metric sidecar is `g_ij=x_j-x_i`, including its sign. Heights can be
recovered from one anchor and spanning-tree gaps precisely when the gaps
are cycle-consistent. On a complete graph the triangle conditions suffice:

    g_ik=g_ij+g_jk.

To prove sufficiency, fix `x_0`, set `x_i=x_0+g_0i`, and use the triangle
identity to recover every gap. Necessity is telescoping. In the transitive
case it is enough to retain the sorted adjacent gaps, the sorted list of
original indices, and one anchor. Positive gap labels along a directed
cycle cannot come from scalar heights: their sum would have to be both
positive and zero. A general tournament Hamiltonian path does not establish
this missing cycle consistency.

Source: a labelled metric row. Target: scalar tournament plus path gaps.
Map: compare values, sort, and record adjacent differences.
Preserved with the sidecar: the entire row, hence every later difference.
Lost without it: scale, additive origin, and the actual differences; without
original index labels, the original neighboring pairs are also lost.
The actual family in§6 is a decisive test of that loss.

The link to the excursion forest is structural and exact at this level:
composition must retain chronological order and the quantities needed to
undo a quotient. Root's affine-block comparator is also transitive when
its curvature values are distinct; sorting changes source residues. Here
sorting changes adjacency unless the original indices are retained. Neither
order theorem produces an allowed move of one fixed Collatz integer.

For a binary absolute-difference tree with disjoint labelled leaves, the
root is a signed sum of the leaves, and modulo two it is their XOR,
independently of the tree shape. Induction at each internal node proves
both statements. Gilbreath's triangle is a DAG with overlapping leaf
occurrences; those multiplicities give the Pascal coefficients. Thus a
variable-depth forest needs leaf/source incidence as well as signs and
metrics. Replacing it by the parity of an unlabelled tree loses more than
just an orientation. No descent-preserving map from Collatz affine
composition to this tree operation has been established.

## 8. Corrected statement and next targets

The incoming note's additivity obstruction should be replaced by:

> Additivity alone does not obstruct an edge encoding. Every binary trace
> is the edge of a unique additive difference seed by the self-inverse
> Pascal transform modulo two; indeed every nonnegative integer trace is
> the edge of an absolute-difference seed by the positive binomial
> transform. Applied causally to a Collatz orbit these constructions
> preserve its trace while placing the orbit computation in an infinite
> seed. They do not produce consecutive primes or a finite-input reduction.
> A useful reduction must constrain seed arithmetic and resources and
> preserve the target magnitude or stopping predicate.

The Rédei sentence should distinguish Gilbreath's exact unit value from its
automatic odd parity. Shared binomial-tail shapes in random models are a
heuristic comparison, not an exact identification of the two conjectures.

Three concrete next questions survive the hostile tests:

1. For consecutive primes, can a feasible triangular sign chamber be
   controlled together with a magnitude bound? Parity is already settled;
   an estimate must bound the odd edge by 1, not merely keep it odd.
2. For an actual Collatz excursion forest, can a boundary quantity use an
   unbounded ordered carry path while admitting a composition inequality?
   Any proposal depending only on the fixed-depth comparisons/parities in
   §6 has already failed the source-height test.
3. Can a constrained seed class, specified independently of the unknown
   future orbit, retain the universal trace while providing a simpler
   verifiable arithmetic invariant? The exact encodings in§§2--3 supply a
   benchmark and an explicit stopping reason, not such an invariant.

## 9. Exact reproduction

Run from the repository root:

    python -X utf8 04-computation/experiments/forest_20260926_gilbreath.py
    python -X utf8 -O 04-computation/experiments/forest_20260926_gilbreath.py

The [script](../../04-computation/experiments/forest_20260926_gilbreath.py)
uses only the standard library and explicit exceptions. Its independent
paths are direct absolute-difference evolution versus Pascal/Newton sums,
and direct initial rows versus signed-matrix inversion. The exhaustive
universes are 8190 binary traces through length 12; 3279 nonnegative traces
over 0..2 through length 7; 19530 actual rows over 0..4 through length 6; 1099
formal sign schedules through length 5; 1099 labelled tournaments through
size 5; and 873 scalar permutations through size 6. Formal schedules are
not counted as feasible chambers.

It also checks the unit edge for the first 1000 primes only, 1092 synthetic
prime-parity rows, 1001 members of the actual family in§6, 128-step causal
codes for 11 specified Collatz starts, dyadic recursions, finite-support
periods, and the explicit hostile controls. These checks do not establish
either open conjecture. Exact normal and optimized stdout is retained in
[forest_20260926_gilbreath.out](forest_20260926_gilbreath.out).

Independent theorem audit: the forest tilings lane accepted §§2--4 and §6,
including the infinite-support boundary, sign determinant, exact orbit family,
and eventual comparison stabilization. The audit found no theorem-level gap.
