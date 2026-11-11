# Complete product saturation classification and a completed Chan–Pak grid

2026-11-11. **PROVED IN NOTES:** the finite-seed saturation classification,
fixed-grade compilers, and the complete d=35 tournament construction below.
**FINITE-EXACT:** the explicitly scoped controls. **OPEN:** the full
Chan–Pak conjecture and the remaining exchange-component lemma. Existing
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

## 4. A new exact reduction for the quantitative domination gap

The remaining Session260 inequality is

    H(U+z) >= H(U)+d(G_z)-1.

The [exchange graph](exchange/PROOF.md) retains both old-path deletion
fibres and dominating-cut fibres of every newborn Hamilton path. Its
exact defect is

    unused_paths + 1 - closed_tree_components + nonnegative_surplus.

At most one closed tree component would give a literal simultaneous
injection and prove the inequality. New local lemmas force branching cuts
to lie below or above the distinguished cut B; closed stars are unique.
Thus it is enough to prove that every closed tree is a star. That final
step is **OPEN**, with exact small controls and explicit hostile examples
against two stronger proposed shortcuts. The
[branching analysis](classification/EXCHANGE_IDEAS.md) and
[two-sided endpoint lemma](history/ENDPOINT_LEMMA.md) isolate the remaining
geometry. These partial theorems do not claim Chan–Pak completeness.

## Reproduction and scope

From the repository root:

```sh
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/classification/classify.py
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/history/check.py
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/chanpak/realize35.py
python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/exchange/probe.py
```

Each receipt names its finite universe. Unbounded conclusions come from
the displayed proofs and actual constructor identities. The concept board
and inheritance pass are in [BOARD.md](BOARD.md). The known exclusions of
Hamilton totals7 and21 are inherited; the remaining global existence
problem is still open.
