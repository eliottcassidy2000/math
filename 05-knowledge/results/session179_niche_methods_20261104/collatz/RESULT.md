# A periodic-template rank atlas and an exact mismatch compiler

**Status: proved elementary construction; finite exact tests; Collatz remains open.**
This is a synthesis extending an older, useful rank construction, not a claim of historical novelty or a proof of universal convergence. No new source-coverage theorem is claimed.

## 1. Recovered older techniques and choice

The repository snapshot inspected was `1b91742ce55a202b73c80d56561ba56d00b503a6`. The following deliberately precede the recent honeycomb / detector work. All source links below refer to that snapshot.

- [Branch-toll rank, October 4](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_branch_toll_rank_20261004.md): proper energy `3^K t^2` at `n-1=2^K t`, invariant along valuation-two steps; precision supplies strict rank decrease. It permits some larger children. This is the main construction recovered.
- [Run-block family, October 5](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_run_block_family_20261005.md): maximal valuation-one runs, exact ternary inverse phases, and a source-bounded expense. The proposed parser below generalizes maximal-run parsing to arbitrary periodic words.
- [Weak-reset family, October 5](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_weak_reset_family_20261005.md): positive rational cycle centers encode weak first descents. This motivates centering the rank at a rational cycle rather than only at 1.
- [Bounded inverse cover, October 7](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_bounded_inverse_cover_20261007e.md): exact native guards and positive smaller-child safety. Its finite guard mass does not prove child grounding. This route was surveyed but not selected for another bank census.
- [Guarded affine lifts, October 4](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_affine_guarded_lifts_20261004.md): legal inverse realization of every finite affine group element, with explicit endpoint/height limitations.
- [Atomic rank interface, October 11](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_atomic_rank_interface_20261011.md): unbounded ranks can be integrated against source atoms, but require authenticated transitions and a genuine terminal.
- [Backward measurement compiler, October 5](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_backward_measurement_compiler_20261005.md): independent deadline to positive measurement floor, with the hidden valuation budget retained.
- [Polynomial atom dual, October 5](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_atom_polynomial_dual_20261005.md): explicit signed polynomial atom minorants; the independent moment oracle remains missing.
- [Chronological gcd rank, September 26](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4493-collatz-gcd-height-rank-certificate.md): high no-dip sources have long full-rank multiplicative blocks, without excluding divergence.
- [Value-time poset and tight spine, September 27](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4514-collatz-value-time-poset-spine-blocks-and-two-place-descent-tree.md): Beatty-length spine blocks and two-place descent-tree structure.
- [Backward-tree rigidity, October 1](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4523-collatz-functional-graph-is-rigid-backward-trees-identify-vertices.md): quantitative 3-adic identification, not connectivity.
- [Finite valuation potential obstruction, September 26](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md): a finite feature correction to logarithmic height cannot decrease globally. A template-specific rank is not a global way around this theorem.

The restricted search also checked related rank/credit/receipt notes for neutral-tower duplication. The forbidden MISTAKES file was neither read nor scanned. This is a focused older-lane survey, not an exhaustive literature search.

## 2. Exact periodic-template parser

Let `U(n)=(3n+1)/2^{v2(3n+1)}` on positive odd integers. Fix a nonempty word `w=(a_0,...,a_(L-1))`, all letters positive. Put `A=sum a_i`, `D=2^A-3^L`, and let `B_i` be the positive carry of the cyclic rotation of w starting at i. Thus the full rotated affine composition is `(3^L n+B_i)/2^A`.

Every `B_i` and D is odd, D is nonzero, and

    2^{a_i} B_(i+1) = 3 B_i + D                  (indices cyclic).

For the supplied source n at phase 0 define `z=D n-B_0`.

**Theorem 1 (exact mismatch without trajectory search).** If z is nonzero and `K=v2(z)`, the number J of initial letters of the actual valuation stream agreeing with `w w w ...` is exactly the largest j for which

    S_j < K,

where `S_j` is the cumulative valuation cost of the first j template letters, allowing repetitions. Equivalently, write `K-1=qA+r`, `0<=r<A`; then J is qL plus the longest prefix of w of cost at most r.

At the mismatch, let k be the remaining precision `K-S_J`, and let a be the next template letter. If k<a, the actual valuation is exactly k. If k=a, the actual valuation is strictly greater than a.

**Proof.** At phase i set `z_i=D n_i-B_i`. The carry identity gives

    D(3n_i+1)=3 z_i+2^{a_i} B_(i+1).

If `v2(z_i)>a_i`, the right side has valuation exactly a_i. If it is smaller than a_i, the right side has that smaller valuation. If it equals a_i, cancellation raises the valuation strictly above a_i. D is odd. Therefore the actual letter equals a_i if and only if `v2(z_i)>a_i`. On a matching step,

    z_(i+1)=3 z_i/2^{a_i},

so its valuation falls by exactly a_i. Induction proves every assertion. Positivity ensures the actual next valuation is finite. QED.

If z=0, n equals the rational cycle center `B_0/D`; if this is a positive odd integer, the same identities prove that it follows the entire word and returns to itself. Such a center is an exceptional cycle case, **not automatically ROOT**. The theorem does not establish that 1 is the only positive cycle.

Example: w=(1,1,3), D=5, B_0=19, n=55. Then z=256, K=8. Template prefix costs are 1,2,5,6,7,10, so exactly five letters match. The actual states are `55,83,125,47,71,107`, and at 107 the actual valuation is 1 rather than the expected 3. No completed ROOT route is required to make this prediction.

## 3. Phase-corrected proper integer rank

Write `s_i=sum_(h<i) a_h`, `e_i=L s_i-A i`, and `c_i=3^{e_i-min_h e_h}`. All c_i are positive integers. For `z_i=D n-B_i !=0`, write `z_i=2^K t` with t signed odd and define

    E_i(n)=c_i 3^{L K} |t|^A,
    R_i(n)=(E_i(n),K), ordered lexicographically.

**Theorem 2.** Every matching step preserves E exactly and decreases K by its positive valuation letter. Hence R strictly decreases on the phase-labelled template path. Each phase has finite energy sublevel sets.

**Proof.** A matching step sends `(K,t)` to `(K-a_i,3t)`, multiplying the unweighted energy by `3^{A-L a_i}`. The phase weight changes by `3^{L a_i-A}`. At the wrap from L-1 to 0 this remains true because `s_L=A`. The product is one. Moreover `c_i>=1`, `3^L>=2`, and `A>=1`, whence

    E_i(n) >= 2^K |t| = |D n-B_i|.

Since D is fixed and nonzero, this bounds n on every energy sublevel. QED.

For w=(2), D=B=1, the construction is exactly the recovered `3^K t^2` rank centered at n=1. For w=(1), it is the maximal-rise rank centered at -1. Arbitrary ordered words give a rational-cycle-centered atlas, with the phase correction retaining the order information lost by a scalar word multiplier.

This is a usable certified run-compression component: compute the mismatch location from one integer valuation, carry the proper rank along the compressed run, and expose the actual deviation at the boundary. It is not a global Collatz rank. Changing templates changes both the energy and its center; no inter-template decrease has been proved.

## 4. Canonicalization and a rejected false improvement

Replacing w by its r-th power does not supply a new cycle center. The full carry and D are both multiplied by

    H=sum_(h=0)^(r-1) 2^{A(r-1-h)} 3^{Lh},

which is odd. Consequently `v2(D_r n-B_r)=v2(D n-B)`, and the parser is unchanged. On phase 0 the rank energy is a positive constant times the r-th power of the original energy. Repeated templates therefore add no new ordering within that phase. Primitive cyclic words are the appropriate atlas entries.

A tempting larger-child extension was checked and deliberately not promoted as coverage progress. The old common-future family `n=223+62208t`, `m=233+65536t` has strictly smaller energy at m. Whenever `3^r | m-1`, the inverse valuation-two tower

    m_r=1+4^r(m-1)/3^r

is integral, positive, has the same energy as m, and maps back to m by r genuine valuation-two steps. Thus it gives lower-energy dependencies with unbounded child/source ratio as r varies, on exact parameter classes of density `3^{-r}`. But every such dependency factors through the original m. It adds no source coverage or new terminal grounding. Neutral towers should be collapsed, not counted as independent proof directions.

## 5. Tests, next obligation, and stopping boundary

`periodic_template_rank.py` uses only standard-library exact integers and explicit exception-raising checks (not optimization-removable asserts). It independently replays actual U trajectories to test the predicted mismatch, phase energy transport, cyclic carry identities, mismatch direction, and repeated-word identities.

The test universe is every word of length 1 through 5 with letters 1 through 4: 1,364 words, 51,702 nonfixed source controls, 145,920 checks, five fixed-point controls. Additional exact long-shadow tests use five templates repeated 10, 50, and 200 times, including at least 800 matched letters. Normal and `python -O` outputs are byte-identical. See `test_results.json` and `test_results.optimized.json`.

The precise next research obligation is a rule for choosing primitive templates at mismatch boundaries together with an independently proper common rank that decreases across the changes. The atlas provides exact local certificates and a simulation-free boundary parser, not that rule. Proving a global selector cannot be replaced by observing long finite runs or by adjoining infinitely many redundant template powers.
