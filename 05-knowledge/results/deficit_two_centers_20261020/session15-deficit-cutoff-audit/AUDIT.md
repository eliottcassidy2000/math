# Independent audit: effective deficit-two cutoff

Date: 2026-10-20. Verdict: **APPROVED** for the revised cutoff `deg D <= 28` and the residual support stated below. The originally requested cutoff 44 is valid but superseded by the candidate's independently checked stronger orbit-connectivity lemma.

Audited source: `../session15-deficit-centers/PREFIX_CONNECTIVITY_CUTOFF.md`.
SHA-256: `f4b36f9ce5a7549311a1d2a1d54cf140f59395e6f32177c8fb274d0b5980cdf4`.

No repository or remote edits were made in this audit. All new files are in the separate `session15-deficit-cutoff-audit` stage.

## 1. Exact approved conclusion and dependencies

Let `D(z)` count isomorphism classes of strong self-converse tournaments with black Hamiltonian-path-orbit deficit two, as in the inherited leaf-reduction packet. The revised argument proves

    D(z) = z^4 + d_10 z^10 + d_12 z^12 + d_14 z^14
                 + d_15 z^15 + d_20 z^20 + d_21 z^21 + d_28 z^28,

where the seven uncomputed coefficients are nonnegative integers. In particular `deg D <= 28`. This does **not** determine those coefficients or prove `D=z^4`.

Dependencies are:

1. The explicitly checked local six-vertex finite lemma.
2. Moon's vertex-pancyclic theorem, checked in the primary paper; its Hamiltonian case includes the needed Camion input.
3. The inherited cycle-cut quotient and odd/even restrictions in `session14-black-leaf-reduction/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md`.
4. The previously audited complete small-order exclusions through nine and positive `S4` control. These census computations were inspected as inherited results, not rerun wholesale in this audit.

The generic bounded-defect finiteness theorem is not needed to obtain this explicit cutoff. No factorial asymptotic estimate, tournament-group order upper bound, or assumed global regularity enters the new proof.

## 2. Converse action and the half-supported permutation lemma

Fix an anti-automorphism `sigma`. The operation sending an ordered directed path `P` to `sigma(reverse(P))` induces an involution on automorphism-orbits, because its square is the automorphism `sigma^2`. The fixed orbits are exactly the blue ones: an automorphism identifying `P` with its image gives precisely the positional reversal of `P` as an anti-automorphism. Thus two black orbits at deficit two really form one nonfixed paired orbit pair.

For `Q=qP`, the relative permutation is uniquely determined by the full ordered words. If its support `S` misses `r_P(S)`, the relevant conjugate has disjoint support. Consequently:

- If both words are blue, `r_Q r_P = q (r_P q^{-1} r_P)` is an automorphism of the same order as `q`, because the two disjoint factors each have that order.
- If the words are in the same ordinary path orbit, the unique relative permutation `q` is itself an automorphism.
- If their orbits are paired, `Q=a sigma r_P P`, hence `q r_P=a sigma` is anti and its square is `q (r_P q r_P)`. Again the disjoint factors each have the same order as `q`.

An even order in any case contradicts the elementary fact that a tournament has no nonidentity involutory automorphism. This needs no odd-order group theorem. Even-order permutations need not have negative sign; the proof correctly uses order, not sign.

The clique bound follows: at most one word can be blue; all black words represent distinct unpaired black orbits; each needs its own distinct mate. Thus a clique of size `t` requires deficit at least `2(t-1)`.

The half-support condition is essential. The independent boundary script finds in regular cyclic `C5` the blue paths `(0,1,2,3,4)` and `(0,2,3,4,1)`, with relative order four but overlapping support and reflected support. It also verifies the two disjoint-support composition identities for all 15 even-order permutations of the first four positions in eight positions.

## 3. Four-cycle prefix and local finite lemma

For a directed four-cycle on `S`, cutting after each vertex that beats `v` gives distinct Hamiltonian words on `S` that can precede one fixed Hamiltonian word of `T-S` starting at `v`. Relative maps are nonidentity powers of a four-cycle, so have order two or four. They are supported entirely on the first four vertices. At `n>=8`, their positional reflections lie among the last four vertices, disjointly. Three available exits therefore force deficit at least four.

The independent verifier imports no producer code. It represents tournaments by adjacency matrices and tests strongness both by exhaustive dominating cuts and by Floyd transitive closure. It finds witnesses by searching for actual directed four-cycles, avoiding dependence on a separate claim that an induced four-tournament is strong. Each retained witness is further checked by constructing the exit words and verifying their pairwise relative permutation orders.

Complete independent counts:

- All `2^15 = 32,768` labelled six-vertex tournaments: `22,320` strong and `10,448` nonstrong. Every strong input has a witness.
- All `2^9 = 512` ways to adjoin two vertices to a differently labelled fixed `S4`: `450` strong and `62` nonstrong. Every strong input has a witness.

Both normal Python and optimized `python -O` executions passed. Required checks use explicit exceptions, not removable assert statements. Exact summaries and deterministic certificate hashes are in `independent_local_results.json`. The full scan alone proves the lemma, so coverage does not depend on uniqueness of the strong four-vertex type; the 512-case scan is an additional independent control.

The primary source was opened directly: J. W. Moon, *On subtournaments of a tournament*, Canadian Mathematical Bulletin 9 (1966), Theorem 1, page 297, DOI [10.4153/CMB-1966-038-7](https://doi.org/10.4153/CMB-1966-038-7). [Primary PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf). Its statement covers cycles of every length from three through the full order at every vertex. A directed six-cycle makes its induced six-vertex tournament strong; the full-order case supplies Hamilton cycles.

If `kappa(T)>=5`, deleting the four vertices of the local witness leaves a strong tournament. A Hamilton cycle of that complement can be cut at the specified vertex `v`. For `n>=8` the prefix certificate applies, contradicting deficit two. Hence `kappa(T)<=4`. No claim is made at orders six or seven by this support argument.

## 4. Rotation blocks and the stronger connectivity estimate

For any Hamilton cycle successor `R`, let `m=|Aut(T) intersect <R>|` and `q=n/m`. The subgroup is `<R^q>` and `m` is odd. Its `q` vertex orbits have common size `m`; each induced tournament is regular because the subgroup acts transitively on it. This is internal regularity, not regularity of the whole tournament. The cycle supplies a perfect directed matching from each orbit to its cyclic successor, including the last-to-first transition.

For a regular orbit with deletion set of size `s_i` and a surviving split `X_i -> Y_i`, averaging internal degrees gives

    |X_i| <= s_i,    |Y_i| <= s_i,    m <= 3 s_i.

Both inequalities have the correct direction. They imply `kappa(T)>=ceil(m/3)` even for one orbit.

When `q>=2`, suppose a total of `s<m/2` vertices is deleted and the survivors split `X->Y`. Every orbit is nonempty. Two split orbits would already require at least `2m/3` deletions, so there is at most one.

- With no split orbit, a cyclic `Y`-label to `X`-label transition occurs. Its perfect matching cannot be destroyed with fewer than `m` deletions, contradicting `X->Y`.
- With exactly one split orbit `O0`, the other labels cannot transition `Y` to `X`. They form an `X` segment followed by a `Y` segment.
- If all other labels are `X`, the matching from `O0` forces every surviving `Y` vertex there to have a deleted mate in `O1`. Thus `b0<=s1`, while regularity gives `a0<=s0`; hence `m<=2s0+s1<=2s`.
- If all other labels are `Y`, the incoming matching gives the dual contradiction.
- If both segments occur, `q>=3`, so `O1` and `O(q-1)` are distinct. The outgoing and incoming matching constraints give `b0<=s1` and `a0<=s(q-1)`, whence `m<=s0+s1+s(q-1)<=s`.

All cases are impossible, giving `kappa(T)>=ceil(m/2)`. The `q=2` case is covered by the two uniform-label cases; it never enters the mixed-segment case. For odd `m`, an integer `s<ceil(m/2)` indeed satisfies `s<m/2`.

This estimate is sharp in the broader self-converse cyclic-block class. Take two identical regular cyclic tournaments of odd order `m`, orient `O0_i->O1_i`, and orient every off-diagonal cross pair from `O1` to `O0`. Alternating the blocks gives a Hamilton cycle with rotation quotient `q=2`; swapping blocks while negating indices is anti-automorphic. A vertex of `O0` has outdegree `(m+1)/2`, so removing its outneighbors disconnects the tournament. Thus `kappa=(m+1)/2`. The independent boundary script verifies the construction, anti-map, rotation map, every smaller deletion set, and a minimum separator at `m=3,5,7,9`. These are not asserted to have deficit two.

## 5. Cutoff and parity bookkeeping

The inherited cycle-cut lemma gives `q<=4`. For `q>=2`, the two connectivity bounds imply `ceil(m/2)<=4`, so odd `m<=7` and `n=qm<=28`.

If `q=1`, the weaker estimate gives `n=m<=11`. Unknown centers have order at least ten, and `n` is odd here, leaving only eleven; the inherited odd-order divisibility-by-three obstruction excludes it. Therefore every unknown center has `q` equal to two, three, or four.

Using odd `m`, the earlier congruence restrictions, and the complete exclusions through nine leaves exactly:

    q=2: n=10,14;
    q=3: n=15,21;
    q=4: n=12,20,28.

This establishes the approved support. No cyclic presentation is treated as an isomorphism-unique parameterization, and no uncomputed coefficient is silently set to zero.

## 6. Separately checked m=7 refinement

The producer proposed a further statement after freezing the audited cutoff manuscript. It is also valid, independently of any unfinished census:

> If a deficit-two tournament of order at least eight has a Hamilton cycle with `q>=2` and rotation-block size `m=7`, its minimum outdegree and minimum indegree are both four.

Proof: the preceding bounds force `kappa=4`; take a separator of size four. Every split orbit needs at least three deletions, so at most one splits. No split is impossible by cyclic matching. A unique split orbit with both types among the other orbits would force `7<=4`, so all other orbit survivors lie on one side.

In the all-`X` case, `s0>=3` and `b0<=s1`. If `s0=4`, then `s1=0`, contradicting `b0>=1`. Hence `s0=3`, `s1=1`, `b0=1`, `a0=3`, and all other deletion counts vanish. Thus `Y` is a single global vertex whose outneighbors are among the four deleted vertices. Every vertex has three outgoing edges within its regular seven-orbit and a further outgoing Hamilton-cycle edge to a different orbit, so its outdegree is exactly four. The all-`Y` case is dual and gives indegree four. The same internal/cycle-edge argument gives minimum semidegree at least four globally, and self-converseness interchanges indegrees with outdegrees, yielding both minima equal to four.

This refinement can safely filter residual cyclic presentations. It does not itself eliminate any of the seven remaining orders.

## 7. Reproduction

Run from the workspace root:

    python3 session15-deficit-cutoff-audit/independent_local_audit.py
    python3 -O session15-deficit-cutoff-audit/independent_local_audit.py
    python3 session15-deficit-cutoff-audit/independent_boundary_controls.py

Normal/optimized local logs, exact JSON summaries, and boundary-control outputs are included. `AUDITED_SOURCE_SHA256SUMS` pins the revised manuscript and inspected inherited dependencies; `SHA256SUMS` pins this separate audit packet.
