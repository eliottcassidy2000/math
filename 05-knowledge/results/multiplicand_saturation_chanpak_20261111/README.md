# Product saturation classification, a completed grid, and a simultaneous injection

2026-11-11. **PROVED IN NOTES:** the finite-seed saturation classification,
fixed-grade compilers, the complete d=35 tournament construction, physical
offset transport, and the all-order quantitative domination inequality.
**FINITE-EXACT:** the explicitly scoped controls. **OPEN:** the full
Chan–Pak conjecture and uniform coverage by the constructed interfaces. Existing
repository theorems are inherited inputs; no literature-priority claim or
numbered canon promotion is made.

## 1. Arbitrary finite multiplicand seeds are completely classified

Let S be a nonempty finite set of positive integers. At each simultaneous
generation retain old values and adjoin products of two distinct currently
available values. Let C be the full attainable closure and M_t the current
maximum. Then

    S_t = C intersection [1,M_t] for every sufficiently large t

holds **if and only if the prime-exponent vectors of the nonunit seeds
have rational rank at most one**. Equivalently, all nonunit seeds are
powers of one common integer base. The unit1 produces no new birth.
For a singleton seed the conclusion is immediate; empty input is excluded
because its maximum is undefined.

The higher-rank obstruction is stronger than occasional failure. There
are attainable missing values X_t for every sufficiently large t with
`X_t/M_t -> 0`. The proof splits into two exhaustive cases on the extreme
rays of the seed valuation cone:

* A ray with one seed a has an exposed first positive layer. If B is the
  largest seed in that layer, its exact maximum at time t is B*a^t.
  Therefore B*a^(t+1) is delayed while the whole construction grows faster.
* If every extreme ray has at least two seeds, each ray has its own exact
  eventual dyadic clock. At a common stage these clocks are
  `c_R^d_R * Q_R^(2^j)` with distinct integer Q_R. A missing next power
  on a slower ray eventually lies far below a faster ray's maximum.

These arguments allow arbitrary relations among the original seeds.
In particular, `{2,3,12}` is not treated as a free seed alphabet.
The complete proof and exact factorization/rank classifier are in
[classification/PROOF.md](classification/PROOF.md).

## 2. Every fixed boundary strip has finite memory

For two independent product seeds a,b, write a^u b^k as (u,k). Define

    c_1=0,
    c_k=k+min_(i+j=k, i<=j)(c_i+c_j+[i=j]).

For every t>=k, the full positive row is exactly

    1 <= u <= k*t-c_k.

The first c-values are `0,3,6,10,14,19,23,28`. The equal-split correction
is the actual forbidden repeated parent, not a fitted loss term. This
gives exact later arrival times and a linear count for any fixed strip.

The [general boundary compiler](history/BOUNDARY_COMPILER.md) permits
relations: near any face with one zero-grade seed A, every bounded set
of positive grades has finitely many cosets modulo ZA. Each coset
eventually becomes an interval `[m_q,s_q*t+b_q]` in its longitudinal
coordinate, with a finite exact recursive algorithm. Representative
carries and equality of actual parent vectors remain in the state.

The [historical recovery](history/HISTORY.md) identifies the older support,
protected-face, context-ideal, carry, and capacity mechanisms used to
derive this compiler. The transfer is proved on arithmetic constructions;
ordinal tournament joins allow equal-count factors and therefore do not
inherit the deleted numerical diagonal automatically.

## 3. The whole d=35 tournament grid is now constructed

For **every r,s>=1**, there is an explicit tournament with

    H = 1+2^s+35*2^(r+s).

The incoming seed-recursion work had s=1 and s>=6; this session closes all
four remaining uniform levels s=2,3,4,5.

At s=2 the exact target factors as `5*(1+7*2^(r+2))`. Actual ordinal
joining of an H5 graph with the inherited coefficient-seven constructor
realizes the product, using two explicit finite heads and a repeatable
tail. At s=3,4,5, three literal six/seven-vertex selected sources have
offset/multiplier pairs `(D,U)=(9,140),(17,70),(33,280)`. The physical
tail identity `H=D+2^t U`, with t=r+1,r+3,r+2 respectively, realizes every
remaining target. Requiring D=1 would exclude exactly these useful sources.

[chanpak/PROOF.md](chanpak/PROOF.md) gives all rows, selected sets, proofs,
and construction orders. `chanpak/realize35.py --r 5 --s 4` emits a literal
example; its default controls directly count all 21 instances with
r=1,2,3 and s=1,...,7 through order17. This is a complete infinite grid,
not a claim that arbitrary odd d or every odd Hamilton count is covered.

The continuation [OFFSET_TRANSPORT.md](chanpak/OFFSET_TRANSPORT.md) supplies
three reusable carriers on which adding a transitive twin implements
`(D,U) -> (2D-1,2U)`. Internal cloning and external tail growth commute as
actual graph operations. Their complete run coefficients also permit an
arbitrary nonempty tournament M at the port, constructing

    H(M)+(a+b*2^t)*Q_2(M),  t>=1,

for `(a,b)=(16,24),(4,42),(32,35)`, where Q_2 counts ordered two-path
covers allowing an empty slot. Only H and Q_2 are needed for replacement
in these carriers. This uses the incoming Session261 module-replacement
mechanism, with a smaller interface than the full cover profile. It gives
new ways to continue a construction, not independent choices of its two
statistics or a universal source bank.

## 4. The quantitative domination gap is closed

For every tournament U and every newborn vertex z, the formerly open
Session260 inequality is now proved:

    H(U+z) >= H(U)+d(G_z)-1.

The [exchange graph](exchange/PROOF.md) retains both old-path deletion
fibres and dominating-cut fibres of every newborn Hamilton path. Its
exact defect is

    unused_paths + 1 - closed_tree_components + nonnegative_surplus.

Here G_z is the bipartite graph on A=N^+(z), B=N^-(z), retaining old arcs
a->b, and d counts dominating sets. The [star theorem](exchange/STARS.md)
proves that every closed tree component is a star and at most one exists.
A branching upper cut would force a strictly smaller lower branching cut,
then a strictly larger upper branching cut. A maximal cut rules this out.
The argument retains actual paths, domination, and deletion seams at every
step; those data are lost by scalar fibre counts.

Matching each component to its edges and stubs now assigns distinct newborn
Hamilton paths to **all old Hamilton paths and all but at most one dominating
cut**. Consequently `H(U+z)-d(G_z)>=H(U)-1>=2*t_3(U)`, where t_3 counts
directed triangles and the second inequality is the inherited odd-cycle
formula. Both bounds are sharp on the inherited 2K2 example with
`H(U+z)=13, H(U)=5, d(G_z)=9, t_3(U)=2`.

The [endpoint lemma](history/ENDPOINT_LEMMA.md), independent review, and
explicit failed stronger shortcuts are retained. This is a universal
lower bound and simultaneous injection, not an arbitrary selectable
Hamilton-count increment. Full odd-spectrum existence remains open.

**Concurrent independent proof integrated:** Session261 commits40eba89fd5
and d42d89d3f9 prove the same inequality by grouping inclusion-exclusion
into even induced components, with a stronger single-cycle bound. Our
[two-decomposition bridge](exchange/DUAL_DEFECT.md) identifies the exact
common slack. Their expansion gives local obstruction certificates; our
star theorem gives a literal simultaneous injection. The resulting
equality class also inherits their effective finite forbidden-basis
theorem, whose complete list has not been computed. Their two-port
bilinear kernel retains the same H,Q2 module data used by our carriers;
it does not supply a freely selectable pair of statistics.

## Reproduction and scope

From the repository root:

```sh
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/classification/classify.py
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/history/check.py
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/chanpak/realize35.py
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/chanpak/offset_ladders.py
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/exchange/probe.py
```

Each receipt names its finite universe. Unbounded conclusions come from
the displayed proofs and actual constructor identities. The concept board
and inheritance pass are in [BOARD.md](BOARD.md). The known exclusions of
Hamilton totals7 and21 are inherited; the remaining global existence
problem is still open.
