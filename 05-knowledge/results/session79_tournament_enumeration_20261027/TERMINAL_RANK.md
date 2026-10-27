# Unit-tail linear rank: strict compression and a sharp scoped obstruction

2026-10-27. **PROVED below:** exact minimal linear dimension of the degree-respecting unit-tail output interface, its core-substitution analogue, and unboundedness. **FINITE-EXACT:** checks and ordinary-count benchmarks. This is not a new Burnside formula, a new future-word equivalence, or an all-algorithm lower bound. No priority claim.

## Inheritance and scope

The closest mechanisms are the [October 20 future-response DP](../tournament_future_boundary_20261020/README.md), Proposition 1, and its unit-tail kernel. The [October 21 integer/packed implementation](../tournament_enumeration_20261021a/README.md) already counts through 200; its packing and odd-factorial denominator are inherited, not new results here.

The inherited quotient is minimal for equality of individual continuation weights. It did not claim minimality under **linear combinations of histories**. Our source is the weighted vector at a fixed committed degree; target is its vector of all final coefficients through the horizon. The map below preserves that output vector, destroys labelled historical weights and arbitrary future marking, and requires stage, degree, horizon and the fixed-point-only future alphabet as sidecars. No individual tournament, core, automorphism, or metagraph is constructed. The board is wordwise equivalence / linear observable rank / finite horizon / unit-tail interpolation / actual computational cost.

## 1. Exact rank theorem

Fix degree `d` and residual bound `r=N-d`, with only cycle length 1 unprocessed. Let `t_1,...,t_h` be distinct historical actual cycle counts. Allow arbitrary weights `w_i` on these existing states. (This is the ambient linear interface, not the single weight vector produced by one fixed run.) Appending `k` fixed points contributes

`C[k,i] = 2^(k(k-1)/2) * (2^t_i)^k / k!`, for `0<=k<=r`.

**Theorem.** Over the rationals, `rank(C)=min(h,r+1)`. Therefore this is exactly the least dimension of any linear carrier that preserves **all residual coefficients** of this fixed-degree block.

**Proof.** Nonzero row scalings turn `C` into the rectangular Vandermonde matrix on the distinct nodes `2^t_i`. Every square initial-row minor on distinct nodes is nonzero. A factorization through a q-dimensional linear carrier requires rank at most q, and storing a basis of the row space achieves equality. QED.

An explicit online carrier uses the original h weights if `h<=r+1`, and otherwise the `r+1` moments `M_k=sum_i w_i 2^(k t_i)`. The latter outputs `2^(k(k-1)/2)M_k/k!`. This is a strictly smaller mathematical linear representation where `h>r+1`; its bit sizes and aggregation costs are not assumed equal to those of original weights.

For the **core** fixed-point substitution `p_1=x/(1+x)`, the row at offset `s>0` is

`Ccore[s,i] = sum_{k=1}^s (-1)^(s-k) binom(s-1,k-1) C[k,i]`,

and row zero is unchanged. This is an invertible lower-triangular row transformation with diagonal one, so the identical rank theorem holds. This substitutes only the unprocessed coordinate p_1; the contributions and signed substitution excess from processed nontrivial cycle lengths have already been incorporated in the committed degree d. It does not apply a scalar substitution to the full ordinary generating function, and no extra factor (1+x)^(-d) is introduced. Historical actual count t must still exclude substitution excess. This rank statement does not assert that the ordinary-state reachability count below also holds for signed core states.

### Strict witness and hostile

At degree 18 the histories `(9,9)`, `(3,3,3,9)`, and `(3,3,3,3,3,3)` have counts `t=(2,4,6)`. For bound 19, their three wordwise-distinct states admit a two-dimensional linear carrier: `(4,-5,1)` is in the kernel since `4-5+1=0` and `4*4-5*16+64=0`. At residual offset 2 its image is 2880, so increasing the horizon invalidates this relation. These are linear weight perturbations, not an assertion that the actual three Burnside weights satisfy this relation. A sufficiently small rational perturbation of any strictly positive weight vector in this direction also remains positive.

## 2. Reachable ordinary blocks and a uniform obstruction

Before processing fixed points, a nonempty history uses only odd lengths at least 3. At degree d, every positive count t with `t≡d mod 2` and `3t<=d` is attainable: take `t-1` parts equal to 3 and one part `d-3(t-1)`. Thus

`h(0)=1; h(d)=floor((floor(d/3)+(d mod 2))/2)` for `d>0`.

When `d=N`, the inherited finite-horizon quotient already merges all reachable t because no fixed point can be added. Apart from this last block its wordwise size is h(d). The explicitly **degree-respecting block-diagonal** carrier therefore has dimensions:

| N | Existing unit-tail states | Minimal sum of block ranks |
|---:|---:|---:|
| 19 | 29 | 28 |
| 31 | 77 | 72 |
| 100 | 819 | 722 |
| 200 | 3302 | 2872 |

The sums are `N²/12+O(N)` versus `N²/14+O(N)`. Indeed `h(d)=d/6+O(1)` and `min(d/6,N-d)` switches at `d=6N/7`; integration gives `N²/14`. This is an asymptotic dimension calculation for this terminal degree-respecting interface, **not** a complexity improvement for the full DP (whose peak occurs earlier).

**No bounded-dimensional universal linear carrier for this interface.** Take `d=6m`, `N=7m-1`. There are exactly m reachable counts `2,4,...,2m`, and `r=m-1`, so the rank is m. This explicitly supplies arbitrary finite lower bounds; no complexity-theoretic assumption is used.

Scope is essential: a request for only one scalar `A_N` has a one-row terminal map and rank at most one. Moreover a direct projection of all degree blocks onto the final coefficient vector has dimension at most `N+1`; the block-diagonal sum above does not prohibit such cross-degree mixing. Precomputing that projection can cost the original problem. Our result rules out a uniform bounded linear interface preserving the declared full residual output vector, not every symbolic summation algorithm and not a polynomial-time recurrence.

## 3. What this does and does not improve

The theorem improves the dimension of a declared terminal linear interface. It
does not reduce the earlier peak of the counting DP, and storing large moments
can increase bit sizes and arithmetic cost. An exploratory implementation that
summed short residual tails directly did not give a meaningful reduction in
work: it repeated continuations that the inherited elimination had shared.
That experiment is not the advertised algorithmic result. No practical speedup
or polynomial all-order running time follows from this paper.

The important distinction is between equality of each future word's response
and equality after linear combinations and the chosen terminal observation.
This is the same distinction that prevents a metagraph cardinality from serving
as a reconstructible tournament carrier. The allowed future operations, the
observation, and the field must all be named before claiming minimality.

## 4. Reproduction and independent evidence

Run from this directory:

    python3 check_terminal_rank.py
    python3 -O check_terminal_rank.py

This standalone checker uses rational elimination, explicit nonremovable
validation, and no external packages or production imports. It independently
builds reachable cycle-count supports through degree 60, tests 2,379 exact
ordinary and triangularly transformed matrices for residual bounds 0 through
12, checks the displayed witness, verifies lower-bound examples through m=12,
and evaluates the dimension formulas. It also tests the different full-series
triangular transform as a comparison, without confusing that transform with
the p1-only core substitution used by the theorem. Normal and optimized outputs
agree. The all-order assertions follow from the proofs, not these finite checks.

Separate development checks reconstructed complete output vectors from the
hybrid weight/moment carrier through order 30 and compared them with direct
rational Burnside sums. The published compact checker focuses on the rank and
support theorem rather than reproducing the pre-existing order-200 benchmark.
