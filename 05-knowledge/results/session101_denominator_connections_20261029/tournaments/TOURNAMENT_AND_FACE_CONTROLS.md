# Session101: exact 175/35 tournament and addressed-face controls

**Status: FINITE-EXACT; independently reproduced, 2026-10-29.**
Historical source pin: `e871dd4513c7dcaa75cb29bb9847aafac35a795c`.
This is a bounded verification of two explicitly defined tournaments and one
35-by-35 rational matrix, not a new full order-seven census.

## 1. Actual tournaments, not number matching

Vertices are `Z/7`; orient `i → j` iff `j−i` lies in the displayed difference
set. There are no ties. Write `A` for the 0/1 adjacency matrix and `S=A−Aᵀ`.
`H` counts directed Hamiltonian paths as ordered vertex sequences, not cycles.

| Tournament | Differences | H | Paths at each start | Aut order | det(xI−S) |
|---|---|---:|---:|---:|---|
| rotational R7 | 1,2,3 | 175 | 25 | 7 | `x(x⁶+21x⁴+35x²+7)` |
| Paley7 | 1,2,4 | 189 | 27 | 21 | `x(x²+7)³` |

All `7!` permutations were checked for paths and automorphisms for each
construction. Exact integer matrix powers followed by Newton identities give
traces `tr(S^k)`, `k=1,…,7`:

- R7: `(0,−42,0,742,0,−14154,0)`.
- Paley7: `(0,−42,0,294,0,−2058,0)`.

Translation acts transitively, explaining the equal per-start counts and
`175=7·25`, `189=7·27`. R7's automorphisms are exactly translations; Paley7's
are exactly `i↦ai+b`, `a∈{1,2,4}`. Reflection `i↦−i` reverses every arrow in
both constructions and is not an oriented automorphism. In particular the
geometric dihedral symmetries of R7 must not be called its oriented
order-14 automorphism group.

**Positive/hostile comparison.** Both objects are regular of order seven,
but their H and skew polynomials differ: the construction and orientation
matter. The 35 coefficient and 175 path count coexist in R7; their ratio is
not itself a five-to-one map or a theorem about Collatz denominators.

Historical result independently recovered here:
`05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md`,
S15 and S17 (read at the pin via `git show`). Its larger regular-tournament
classification is not re-enumerated in this check. See also
`history/PROVENANCE_MAP.md`, B3, B5, B6, for source roles and symmetry cautions.

## 2. The two 35-dimensional middle face layers

Let `V={0,…,6}` and index rows and columns by actual triples `S,R⊂V`.
Define `K[S,R]=1` iff `S∩R=∅`. This is the disjointness matrix of `KG(7,3)`.
The exact rational candidate

`B[S,R] = (−1)^|S∩R| / (4·binom(3,|S∩R|))`

has entries `1/4,−1/12,1/12,−1/4` for intersections `0,1,2,3`.
The script checks **all 1225 entries of both `KB=I` and `BK=I`**, and an
additional nonconstant rational-field round trip. Every row of K has sum 4.

For an addressed additive triple field `f`, set
`g(Q)=Σ_{R⊂Q, |R|=3} f(R)` on four-sets Q. Complementing the four-set row
address, `h(S)=g(V\S)`, gives exactly `h=Kf`, so `f=Bh` over Q (and any
characteristic-zero field). This is a lossless linear transform of
**addressed additive face fields**, not a claim that an unlabelled histogram
of tournament types reconstructs its tournament. Forgetting addresses loses
overlap information. Nor does this show independent higher orientations arise
from edge tournaments.

The elementary identity can also be checked using only four intersection
cases: if `t=|S∩R|`, then `(KB)[S,R]` equals

`Σ_u binom(3−t,u) binom(1+t,3−u) (−1)^u / (4 binom(3,u))`,

which is 0 for `t=0,1,2` and 1 for `t=3`. This explains why the full matrix
product succeeds; it is not merely a numerical inversion heuristic.

Historical reflection (provenance, with the specific identity independently
verified here):
`07-reflections/tournament-subset-decks-higher-order-faces-and-kneser-middle-duality-codex-20260821.md`,
§1 and equation (1.1); routed by `history/PROVENANCE_MAP.md`, C3.

## 3. Scope and reproduction

The common higher-order lesson with the Collatz lanes is an information-loss
diagnosis: an operation may require a marked seam/carry/phase or an addressed
overlap coordinate. No map equating Collatz denominators 35/175 with path
counts or face dimensions is constructed or asserted. No Collatz convergence
consequence, full seven-vertex census, unlabelled reconstruction, or general
higher-tournament theorem follows from this bounded computation.

From this directory:

```sh
python3 verify_tournament_controls.py > verification.json
python3 -O verify_tournament_controls.py > /tmp/tournament_controls_optimized.json
cmp verification.json /tmp/tournament_controls_optimized.json
```

Both modes passed with byte-identical JSON. The checker uses only the Python
standard library and explicit exceptions, so optimization cannot disable its
verification. No publication was performed by this lane.
