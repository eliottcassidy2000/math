# Independent audit: synchronized source pencils

**PASS, without mathematical correction.** The general rank-one composite obstruction, four-probe equivalence, quantitative slope gap, and rank-pruned prime inverse are valid within the stated physical-source grammar. No claim that every prime or every admissible odd target occurs is established. The preprint's major conclusions are not dependencies.

## Frozen inputs

The source packet is `session201-rank-pencils/`. `INPUTS.json` records every frozen source and dependency hash. Principal hashes are:

- `PROOF.md`: `549ec35531687e38c9f63758eb36ba5603c79fbd6e5036ca300ab4f9f3998986`.
- `verify.py`: `5959d4dc5ab8cf0be857c8102fc128501ef871e44969545a20b83785a8c050b9`.
- `SOURCE-CONTEXT.md`: `ccff9229453bf2366e168e35e6f250d7a8c87011e1a701c11690a9167d753786`.
- Both author receipts: `0de68d3f40247f9da6be50f8a7a2ba6cf6fffe4fb7b20c12f7908d7fe7f7b90c`.

The three frozen dependencies match the separately audited physical-adjoint, synchronized-joint-inverse, and two-return release proofs byte-for-byte. Their scopes are inherited as written. A later release status promotion needs its own exact change/hash verification.

## Mathematical review

### 1. Primitive factorization does not require primitive physicality

For a synchronized rank-one pair, w is an actual source. Its coordinate gcd g is positive and odd because its path count D is odd. The primitive integral quotient p=w/g has nonnegative integral positive coordinates: the coordinate change uses integer linear combinations of the raw coordinates. Its deletion coordinate is a positive integer, since the physical deletion d_w is positive and divisible by g.

Rank one implies q=c p over the rationals. Primitivity of p and Bezout imply c is an integer a, without asserting that p is an actual source. The positive bilinear matrix has all entries nonnegative and its deletion-support two-by-two block has entries at least three. Consequently `K(v,p)>=3 d_v d_p>=3` for every actual source v. The statement is a cone inequality and is valid even if p has no physical recipe.

Apply the exact joint cone to v=t. Since `m=g K(t,p)` and `h=a K(t,p)` with positive common factor, `h<m<2h` gives `g/2<a<g`. In particular g cannot equal one, and odd g is at least three. Therefore both the balanced readout m and every displayed Mersenne target have two integer factors exceeding one. No false divisibility-to-physicality inference appears.

The same common generator acts linearly on both columns, so it cannot increase their rank. Rank zero cannot occur on the synchronized orbit, because w is physical and the auxiliary readout h is positive. Rank-one collapse is therefore absorbing, proving the prefix exclusion. The non-return example OJOOJ is explicitly rank one but not on the triangle ray, so the new statement extends the earlier repeated-return setting.

### 2. Four genuine sources detect rank over the rationals

The four displayed source words are empty, O, Z, and ZJ. Their raw row matrix P has determinant -4. Q has determinant -1, and the displayed response matrix PQ has determinant 4. This is an invertible rational measurement map; the proof does not claim that it is unimodular over the integers.

For every synchronized pair, m_0 is positive. Vanishing of the three determinants Delta_i makes the four-component response vectors proportional. Invertibility of PQ then makes q and w proportional. The converse and the interpretation as independence of h/m from every physical source follow immediately. Thus the four-source bank proves separation, not target occupancy.

In rank two at least one determinant is nonzero. Every m_i and h_i is odd because it is a common-word physical joint readout, so this determinant is a nonzero even integer and its absolute value is at least two. The physical product upper bound gives `m_0<=2U_w` and `m_i<=6U_w`. Dividing yields exactly the stated gap `1/(6U_w^2)`. This bounds a rational slope difference; it is not a height-independent gap and does not manufacture two estimates on a target fiber.

### 3. Completeness and costs of the pruned inverse

Any word reaching a rank-two final pair has rank two at every prefix. Thus pruning each rank-one node removes no rank-two state, no prime witness, and no minimum-cost prime witness. This also explains why it is legitimate to stop extending a collapsed node, rather than merely filtering final outputs.

The complete decorated pair, not w alone, is the deduplication key. The example OJOJ versus ZZOJ has the same w but different auxiliary q and different rank; forgetting q would invalidate pruning. Root-change cost belongs to the complete pair's common word. The inherited physical inverse and target decomposition retain the sum of the two source costs, so the minimum-J conclusion is exact in its specified precursor convention.

The height bound is the inherited numerical cutoff on both physical factors, not just on the smaller factor or on the auxiliary vector separately. The target n=3 is handled outside the nontrivial balanced inverse. Rank two is only a necessary condition for prime output: the displayed target 51 is a literal composite counterexample to sufficiency.

## Full numerical comparison and independent graph controls

The author verifier was replayed serially at cap 5,001 under ordinary Python and Python -O. Both outputs match both frozen author receipts byte-for-byte.

The separate `independent.py` imports no author modules. It uses raw matrix operations, exact rational rank elimination, increasing-height bucket inventories rather than the author's heap, trial-division primality, and direct depth-first traversal of literal directed paths rather than subset dynamic programming.

Both independent runs are byte-identical. Their complete cap 5,001 comparison gives:

- 126,910 actual source profiles through U=5,000.
- Square-root factor height 141: 641 source profiles and 499 synchronized states.
- 42 rank-one states removed; all 457 rank-two states and their minimum J costs preserved.
- 136,144 in-cap source-pair presentations before pruning, 126,785 after pruning.
- Exactly 649 prime C1 targets and 639 balanced prime targets at least 5, with the same minimum precursor J cost in direct, unpruned, and pruned constructions.

The forward cutoff is complete because F=D-E+U>=U+1 and G=D+2U>=2U+1. Each transition strictly increases cover height, so the finite inventory is not a guessed word-length sample. The square-root cutoffs rely on the separately reviewed inverse theorem. The frozen author checker uses the common maximum cutoff to compare all targets through 5,001; it does not claim a separate all-target existence result.

Independent checks additionally reconstruct 1,996 common-word probe readouts, test 26,922 physical-source/primitive-cone factor pairs, verify the cone inequality on generating rays, verify all four-probe determinants and rational gaps, and retain the equal-w/different-rank hostile.

Five independently built graphs were enumerated by direct directed-path DFS:

- Rank-one Mersenne composite 171, 10 vertices, one Hamilton cycle.
- Rank-two Mersenne prime 139, 9 vertices, one Hamilton cycle.
- Rank-two composite 51, 7 vertices, one Hamilton cycle.
- Rank-two prime 11, 5 vertices, one Hamilton cycle.
- Rank-one balanced composite 45, 8 vertices, one Hamilton cycle.

The Hamilton-path counts agree with the literal compiler in every case. Hamilton cycles are counted once by anchoring the first vertex, not once per rotation.

## Attached-paper boundary

The supplied PDF hash is `f85bf570334601930ae862615c4f13365d2cd45370bb0b853252bf15445d9abf`, matching the source-context note; the PDF has 39 pages. Its first page states the 12+epsilon genus-two bound and the 24+epsilon abc-type exponent, not the standard 1+epsilon abc bound. The cited passages do use two independent sections, a dimension-two space away from an exceptional set, and a separate basepoint-free argument before concluding no common zero.

The two inequalities quoted from Sections 6.3.6-6.3.7 concern the same base-changed fiber and field. The difference between coefficients 5/3 and 3/2 is 1/6; the source-context note correctly presents this only as the formal coefficient mechanism, conditional on the manuscript's uniform geometric, metric, and base-change inputs. This audit does not verify the preprint's proof of those inputs or its abc/Mordell/Siegel-zero conclusions.

The new tournament statements follow from their own physical grammar, integral cone, common-word pairing, and finite-dimensional linear algebra. No line-bundle or adelic-lattice identification, finite-exception transfer, field-descent theorem, or preprint headline is used. The connection is a motivated independence test with an exact target-side rank map, not a theorem transporting the manuscript's geometry to tournaments.

## Reproduction

From this audit directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

Finite checks certify the implementation and declared examples. The all-order results rest on the reviewed proof and frozen physical dependencies. Universal prime coverage and the complete {7,21} spectrum remain OPEN.
