# The bivariate packing polynomial and Hadamard-style flatness

**Proved.** The regular-tournament autocorrelation defect is
read exactly by two coefficients of the bivariate odd-cycle packing
polynomial, but those contributions cancel in every diagonal cover
read at that degree. A first vertex-deletion read recovers them.
Uniform module substitutions preserve the diagonal, and each cyclic
wrap consumes one matched cover coordinate of the compensated ×7
construction. This does not prove all-odd coverage.

## 1. Same method, opposite existence statements

For the user's normalization let G_T(x,y) count packs of disjoint
directed odd cycles of length at least three, with x per
cycle and y per uncovered vertex (singleton). Every odd
cycle or singleton is an independently colorable block in
the [Rédei–Berge specialization](https://arxiv.org/abs/2307.05569).
A nontrivial odd cycle also has weight 2. Therefore P_T(q)=G_T(2q,q):
q optional ordered path slots give P(1)=H, P(2)=R and P(3)
the exact third cover count.

In a real **circulant Hadamard** matrix the periodic autocorrelation
of its first sign row must be n at shift zero and zero at every
other shift. Writing its row as a group-ring element h, this is
h h*=n times the identity; in the Fourier domain every squared
modulus is n. The [recent proof](https://github.com/openai/math/tree/main/preprints/The-circulant-Hadamard-conjecture-September-23-2026)
uses cyclotomic and coefficient constraints to eliminate orders
other than 1 and 4, and the [Lean scope note](https://github.com/mathcommunity/oai-math-revised/blob/main/lean/docs/179.md)
records the exact classification and witnesses. The analogy with
Chan–Pak is **not** that two sets of size two coincide: 1 and 4
are the only allowed *orders* there, whereas 7 and 21 are
excluded *counts* here. Seven-vertex tournaments are plentiful.

The useful transfer is to retain the observer that detects
the obstruction before projecting down to one scalar.

## 2. An exact Hadamard defect inside G

For a regular tournament of odd order n≥5 let A be the skew
sign matrix (zero diagonal, +1 for an outgoing arc, −1
otherwise). Put J=11^T and

    E(T) = ||AA^T − (nI−J)||_F².

The bordered sign matrix

    M = I + [[0, 1^T], [−1, A]]

has ||MM^T−(n+1)I||_F²=E(T), because A1=0. Its defect vanishes
exactly for doubly regular tournaments. **This is a bordered
skew-Hadamard matrix, generally not circulant.** For a cyclic
tournament with skew row epsilon, the same E equals

    n times sum over s≠0 of (C_epsilon(s)+1)²,

so its target correlation is −1 rather than zero. With the
unnormalized Fourier transform it also equals the sum of
(|hat epsilon(j)|²−n)² at nonzero frequencies. This is an exact
centered fourth-moment read.

Write C5(T) for directed 5-cycles, D2(T) for unordered pairs
of disjoint directed triangles. Under the unweighted
cycle convention of G these are respectively its coefficients
of x y^(n−5) and x² y^(n−6) (the latter is zero at n=5).
For two regular tournaments T and U of the same order n≥7,
with Delta denoting T minus U, we obtain

    Delta C5 = −Delta E/32,    Delta D2 = Delta E/64.

**5-cycle proof.** The adjacency matrix is (J−I+A)/2.
Its constant eigenvalue is (n−1)/2; on the sum-zero subspace
its eigenvalues are (−1+i t_j)/2 with the t_j real and paired
by sign. The trace of its fifth power counts 5 times C5:
a repeated vertex would split a closed walk of length five
into shorter closed walks, each of length at least three.
The linear and cubic imaginary terms cancel in pairs, and
the second moment sum t_j²=n(n−1) is fixed. The fourth
moment sum t_j^4 is n²(n−1)+E. The real fifth-power
coefficient at t_j^4 is −5/32. Dividing by 5 gives the law.
The flat value is n(n²−1)(n−2)(n−3)/160, consistent with
[the known extremal theorem of Savchenko](https://arxiv.org/abs/2403.19555)
and the repo's [Paley count](../../../01-canon/theorems/THM-309-five-cycle-design-paley.md).

**Disjoint-triangles proof.** A regular tournament has fixed
triangle totals T3 and fixed per-vertex incidence t3:
for outdegree d=(n−1)/2 the latter is d(d+1)/2.
The triangle incidence of an arc i→j is
t_ij=(n−(AA^T)_ij)/4. Indeed regularity makes the common
out- and in-neighbor counts equal, and t_ij is d minus the
common out-neighbor count. Inclusion–exclusion gives

    D2 = binom(T3,2) − n binom(t3,2) + sum_arcs binom(t_ij,2).

The mean t_ij is (n+1)/4, so only the squared deviations vary.
They contribute E/64. The factor 2 comes from the two
matrix entries for each unordered pair.

Consequently the entire total-degree n−4 part of G_T−G_U is

    (Delta E/64) x(x−2y)y^(n−6).

On x=2y it disappears. This is why the coefficient of q^(n−4)
in the **full** slot polynomial cannot see this Hadamard defect:
2 Delta C5+4 Delta D2=0. Higher or lower terms need not agree.
At n=7, interval7 and Paley7 have the same 14 triangles;
their (C5,D2) are (28,14) and (42,7). Interval7 has E=448,
Paley7 has E=0. Their q³ coefficients are identical, while
their q coefficients are 34 and 48. Both H=175 and H=189
are realized, not examples of 7 or 21.

For cyclic n≡3 mod4, any nonflat row has E≥64n:
C(s)+1 is in 4Z, is symmetric under s→−s, and sums to zero
off the origin. It must contain a pair at least +4 and a
pair at most −4. Thus D2 rises by at least n above the flat
value and C5 drops by at least 2n; interval7 attains both.

## 3. What a vertex-deletion read sees

Differentiate while holding x fixed. Marking a singleton
is exactly deleting a vertex:

    partial_y G_T(x,y) = sum_v G_(T−v)(x,y).

More generally (1/k!) partial_y^k marks a k-set of singletons,
so Taylor expansion around y=x/2 gives the exact recovery

    G_T(x,y) =
      sum_k (y−x/2)^k [sum_{|S|=k} P_(T−S)(x/2)].

This is one way to recover the transverse packing information:
the first k aggregate deletion polynomials agree iff the
difference in G is a multiple of (x−2y)^k (levels 0 to k−1).
For equal-order regular T,U the coefficient of q^(n−5) in
the first deletion-sum difference is −Delta E/16. Thus
deletion recovers the defect the diagonal coefficient lost.

A stronger physical control is on six vertices. With outgoing
masks A=(38,8,2,5,15,30) and B=(58,4,1,6,14,30), the
P's are identically q^6+8q^4+8q², but G_A−G_B=x(x−2y).
The first aggregate deletion polynomial separates them by
−4q. The direct check does not assign the two tournaments
the same internal cycle structure just because they have 17 paths.

## 4. Why the missing coordinate stays missing in uniform composition

For m labelled nonempty runs selected inside one module,
let epsilon_i be the sign attached per vertex in its ith run.
Expanding the Hamilton count on each chosen run's vertex set
as an odd-pack sum yields a marked transfer:

    sum over J subset of [m] of (−1)^(m−|J|)
        G_B(x+2 sum_{i in J} epsilon_i, y+sum_{i in J} epsilon_i).

Indeed each odd atom chooses the leftover pack or one
run. A selected nontrivial cycle has weight 2 epsilon_i
and a singleton weight epsilon_i, since both sizes are odd.
Inclusion–exclusion makes the selected runs nonempty.

In the full exterior, list each directed pack cycle that
crosses the module by its ordered exterior vertices and
module-run placeholders. Anchor at its least exterior vertex
to avoid rotational duplication, and order mixed cycles
by those anchors. Arcs at the module boundary are uniform,
so the run words are free once their lengths and sets are
chosen. The total length of each mixed cycle must be odd:
average the two signs on its group of placeholder lengths
with the required parity; multiply by x per mixed cycle
and by G of the unused exterior. This enumerates each actual
pack exactly once and proves that uniform substitution uses
only shifts (x,y)→(x+2s,y+s) of G_B in a linear operator.

All these shifts preserve x−2y. **A difference divisible by
(x−2y)^k stays hidden to order k under every common uniform
exterior.** In particular the six-vertex pair is
indistinguishable in *every* uniform scalar Hamilton read,
although a deletion separates it. Conversely equality of
all strong uniform Hamilton reads forces equality of P:
the [matching-exterior tests](../atomic_symmetry_descent_20261117/README.md)
recover every ordered nonempty cover count by triangular
inversion. So the full P, not just H or P(1),P(2), is exactly
the uniform scalar interface.

## 5. A concrete cost for the ×7 proof route

Let C=C3[A,B,1] be the cyclic gluing with a singleton. A mixed
cycle must use that singleton and one run in each arm, whose
two lengths sum to an even number. The bivariate transfer
above gives, with Delta_+ P(q)=P(q+1)−P(q) and
Delta_− P(q)=P(q−1)−P(q),

    P_C(q)=q [P_A(q)P_B(q)
               +Delta_+ P_A(q) Delta_+ P_B(q)
               +Delta_− P_A(q) Delta_− P_B(q)].

If B' scales B by 7 through color r and first fails at r+1,
then one common cyclic wrap with **any nonempty A** scales
through r−1 and first fails at r. The discrepancy at r is

    r [P_A(r+1)−P_A(r)] [P_B'(r+1)−7P_B(r+1)],

which is nonzero: extra colors allow new assignments for
any nonempty A. P_B(0)=P_B'(0)=0. Thus each wrap consumes
one matched cover coordinate exactly, not merely an upper bound.

The [previous compensated ×7 switch](../fermat_flip_stern_20261117/FACTOR_SEVEN.md)
matches at 1 and 2 but its exact residual at 3 is
23,651,115,982,200. So after one wrap its H scales exactly,
but its new two-cover count can **never** scale by 7 for any
common nonempty A. To make it repeat once more, a physical
third-coordinate compensation is necessary; to permit k common
wraps by this route one needs matching through k+1. Infinite
reuse of a fixed nontrivial ×7 pair is impossible in this
chain because finite monic P's cannot satisfy P'=7P identically.
Changing blocks or reading nonuniform contacts remains possible.

This is the useful analogy to Hadamard rigidity: matching a
small set of observations can leave a hidden defect; a
second operation may expose it. Neither the circulant
order theorem nor these identities supply the absolute
atom selector needed for Chan–Pak.

See [ORBIT_PHASES.md](ORBIT_PHASES.md) for an exact cyclic
symmetry law and a hostile autocorrelation control. Run
`python replay.py` for the literal polynomial checks
(normal and optimized Python give identical outcomes).
