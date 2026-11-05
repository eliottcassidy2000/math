# The reachable set of the span-2 ribbon transfer: its congruence lattice, an exact tripling block, sparse response sets, and the local-limit form of the collective invariant

**Session** mac-mini-2026-11-04-invariant (worktree `math-wt-invariant-20261104`, a sibling of the repository).
**Directive** (owner, 2026-11-04): clear every worktree whose contents are already on origin/main and stop the
practice; let their contents inspire the session; find the collective invariant of the reachable set in Z^13.
**Status in one line.** Process: 26 worktrees removed (all clean, all HEADs in origin/main), one kept for the owner's
decision, the hygiene rule and tooling committed (17d7038509). Mathematics: THM-4881 (A), (B) PROVED (the exact
congruence coset of the reachable set; eventual exact tripling of the block `1-.-2`), (C), (D), (E) FINITE-EXACT
(response sets are sparse and saturated, so the near-join density route is closed inside the family; residue bias and
exponential-sum decay; proximality), and HYP-9474, which states the invariant in the only form the data supports: an
integer-scale local limit law for the linear random walk on Z^13. Chan-Pak remains OPEN; the ribbon tiling law
(HYP-9473) remains open. Canon: THM-4881, HYP-9474. Scripts: `04-computation/experiments/ribbon_spectrum_20261104/`
(reach.py, reach2.py, bsets.py, collective.py, fourier.py, with .log receipts).

## 1. Worktree hygiene (process)

At the start of the session `git worktree list` showed 27 worktrees besides the main checkout, 24 of them under
`/private/tmp` created by the old default of `agents/new_session_worktree.sh`; the data volume had hit 100 % during
the previous session. For every worktree the audit recorded the branch, `git status --porcelain --untracked-files=all`,
whether HEAD is an ancestor of `origin/main`, and `git cherry origin/main HEAD` (receipt:
`04-computation/experiments/ribbon_spectrum_20261104/worktree_inventory_20261104.txt`). All 27 HEADs were contained in
`origin/main` with no unmerged patch. 24 were clean and were removed; two more had only a whitespace-only line-ending
phantom (`lrc14_removal_measure_kps_S38.out`) or a log that grew after close-out, and were removed after that was
verified. One is kept: `math-wt-collatz-procgen-20260922` holds 699,167 untracked files (8.4 GB, `scratch/lrc_certs/`)
that were never pushed; deleting never-pushed work is the owner's call. Free space went from 2 GB to 28 GB.

Tooling committed: `agents/worktree_hygiene.py` (reports every worktree; `--remove [--path <wt>]` removes only clean
worktrees whose HEAD is in `origin/main`, tolerating the line-ending phantom; dirty or unmerged ones are kept and
reported), `new_session_worktree.sh` now creates worktrees as siblings of the repository, refuses `/tmp`, and warns at
six or more worktrees or under 20 GB free; `cleanup_session_worktrees.sh` wraps the Python script;
`finish_session.py` prints the removal command; AGENTS.md and agents/README.md state the rule.

## 2. Inspiration from the removed lanes

The removed worktrees were the klein/codex spectrum lanes of 2026-11-02/03 (boundary carry returns, clock boundary,
digit memory in all bases, multiport grammar, nonzero reset, prime inverse ports, six roots and hidden bits, triple
re-entry interfaces, MRB square spectrum, tournament orders 1499). Their common object is the ternary return map
`h -> 3h + {0, 2, 4}` (the C2 births with digits 0, 1, 2): every odd number other than 7, 21 descends in that forest
from one of the six roots `1, 23, 25, 63, 65, 67` (triple_reentry_interfaces_20261103/RETURN_FOREST.md), so Chan-Pak
would follow from "every attained h has a realiser ready for each digit". The ribbon family has no such digit system
(THM-4880 (H)), but it does have an exact *multiplicative* block: the two letters `1-.-2` triple the count after at most
thirteen repetitions, for every prefix (THM-4881 (B)); it is Moon's chain of 2-blocks seen from the automaton, and its
matrix has characteristic polynomial `x^13 (x - 3)`. The lanes' "readiness" question becomes, in the automaton, the
question of which root responses occur at a given count, answered in section 5: sparsely.

## 3. The reachable set and its congruence coset (THM-4881 (A), PROVED)

Let `y(u) in Z^14` be the reduced state of THM-4880 (`y(empty) = e_1`, `H = y_2`, `y_1 = y_2` on nonempty words,
letters acting by the integer matrices `A_m` of span2_intrep.txt), and `Y_n = {y(u) : |u| = n}`:
`|Y_n| = 1, 2, 8, 24, 68, 194, 561, 1709, 5309, 16775, 52741, 165816, 525630` for `n = 1..13`. Let `y_ref = (1, 1, 0, 2,
0, ..., 0)` (the least state of level 7) and let `L` be the lattice with echelon basis `2e_1 + 2e_2`, `2e_3`, ...,
`2e_12`, `4e_13 + 4e_14`, `8e_14` (index 2^16 in the hyperplane lattice `{d_1 = d_2}`).

**Theorem (A).** For every word of length `n >= 7`, `y(u) in y_ref + L`; that is, `H` is odd, `y_3, ..., y_12` are
even, `y_13 = 0 (mod 4)` and `y_14 = y_13 (mod 8)`. *Proof.* `L` is the Z-span of the differences of the 561 states
of level 7 (flint HNF, collective.py); `A_m L subset L` for all four letters (each basis row's image reduces to zero)
and `A_m y_ref - y_ref in L`, so `A_m (y_ref + L) subset y_ref + L`; level 7 lies in the coset, hence every later
level. The levels 3-6 lie in proper sublattice cosets (index 128 at level 3, 16384 at level 4; the lattice of
differences stabilises at level 5). The coset is the whole congruence structure: with 400 random differences from each
level the HNF is identical from level 5 on, and no coordinate is constant modulo 16.

## 4. The block `1-.-2` triples exactly (THM-4881 (B), PROVED)

The matrix `A_B = A_{-2} A_{1-}` of the two-letter block has characteristic polynomial `x^13 (x - 3)` (exact, flint),
rank 4, so `A_B = 3P + N` with `P` a rank-one idempotent, `PN = NP = 0`, `N^13 = 0`, and `A_B^k = 3^k P` for `k >= 13`.
Hence for every word `u` and `k >= 13`, `H(u B^k) = 3^{k-13} H(u B^13)`; on 40 random prefixes the exact tripling
starts already at `k <= 3` (`H(u B^k)`: `1, 11, 33, 95, 285, 855, ...`; `23, 63, 185, 555, 1665, ...`;
`H(B^k) = 1, 1, 5, 15, 45, 135, ...`). The same holds for `A_{1-} A_{-2}`. Reading: the first vertex of the block
(`q+1 -> q`) is beaten by everything before `q` and beats everything after, so it is a singleton strong component; the
second (`q+2 -> q`) continues a chain of 2-blocks in the sense of Moon (THM-4744: each block of size 2 contributes the
factor `2^(2-1) + 1 = 3`). No additive structure: the residual `H(uB) - 3H(u)` is unbounded (THM-4880 (H1)).

## 5. Root-response sets are sparse and saturated (THM-4881 (C), FINITE-EXACT)

For a word `u` with last vertex `t`, the post-sink state `A_{--} y(u)` reads off the rooted triple `(h, e, R)` with
`h = H(u)`, `e = e_t(u)` and `R = R_t(u) = sum_X h(u[X]) s_t(u - X)` (the lane's `(h, e, B)`; `(A_0 y)_6 = R - e`,
`(A_0 y)_4 = 2R`). The near-join `u . -- . -2` has count `h + 2R` (HYP-9453 at `k = 1`), so HYP-9453 (2)'s mechanism
- "the set `B(h) = {R}` is dense in `[h, ch]`" - would give the tiling law by induction. bsets.py computes, for every
odd `h <= 99` and over all span-2 words of length `<= 20` (BFS on reduced states pruned by `H <= 99`), the set of
pairs `(e, R)` with `R <= 4h`:

| h | pairs (e, R) | distinct R | parity-admissible R in [h, 4h] | density | distinct e | unchanged from length 16 |
|---|---|---|---|---|---|---|
| 9 | 11 | 8 | 28 | 0.29 | 6 | yes |
| 15 | 13 | 10 | 46 | 0.22 | 5 | yes |
| 33 | 20 | 17 | 100 | 0.17 | 11 | yes |
| 45 | 28 | 21 | 136 | 0.15 | 11 | yes |
| 63 | 6 | 6 | 190 | 0.03 | 5 | yes |
| 65 | 29 | 22 | 196 | 0.11 | 16 | yes |
| 99 | 44 | 34 | 298 | 0.11 | 22 | yes |

For all 47 attained `h <= 99` the sets are identical at lengths 16 and 20 (saturated; a proof of completeness would
bound the length of a word with given `(h, R)`, not attempted). The all-tournament sets from every root of every
tournament of orders `<= 8` are denser (17-56 % for `h <= 45`) but far from full; the ribbon sets are proper subsets
(e.g. `h = 9`: all tournaments reach 14 of 28 values, ribbons 8). So inside the span-2 family the covering of every odd
number by `h + 2R` cannot come from dense `B(h)`: the tiling law (THM-4880, every odd `<= 100,000` except 7, 21, 47)
holds while every per-`h` response set is sparse. The same holds for the three letters directly: at level 10 the
extension sets `{H(ua) : H(u) = h}` have 3-40 members with gaps (extdensity.py, THM-4880 note section 9).

## 6. Residues, multiplicities, exponential sums (THM-4881 (D), FINITE-EXACT) and proximality ((E))

Let `r_n(N) = #{u : |u| = n, H(u) = N}` (four letters, `4^(n-1)` words). On the solid window `[49, d_n]` of level
`n`: median `r_n` = 79, 108, 100, 304, 232 and maximum 666, 2,226, 9,635, 40,156, 158,580 for `n = 9..13`; the
minimum is 2-8 (attained near the window's upper end, where the solid prefix stops). Residues of `H(u)` on the window
(collective.py): modulo 7, 8, 16 the class frequencies are uniform to within 0.03 at `n = 13`, but modulo 3, 5, 9 the
deviation grows with `n` (0.07, 0.02, 0.08 at `n = 9`; 0.21, 0.16, 0.13 at `n = 13`): ordinal sums of 3-cycles and of
the strong 4-tournament over-represent multiples of 3 and 5. The normalised exponential sums
`|S_n(theta)| / 4^(n-1)`, `S_n(theta) = sum_{|u| = n} e(theta H(u))` (fourier.py):

| n | 1/3 | 1/5 | 1/9 | 1/7 | 1/8 | 1/16 | 3/11 | 5/13 | sqrt2-1 | golden | pi-3 | 0.1234... |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 9 | 0.190 | 0.165 | 0.081 | 0.106 | 0.063 | 0.057 | 0.120 | 0.054 | 0.109 | 0.052 | 0.109 | 0.062 |
| 11 | 0.265 | 0.186 | 0.073 | 0.081 | 0.039 | 0.050 | 0.067 | 0.024 | 0.054 | 0.052 | 0.085 | 0.031 |
| 13 | 0.324 | 0.227 | 0.112 | 0.043 | 0.020 | 0.021 | 0.044 | 0.022 | 0.032 | 0.025 | 0.049 | 0.003 |

At `1/3` and `1/5` (and slowly at `1/9`) the sums grow: these are genuine major arcs carrying the multiplicative bias;
everywhere else they decay geometrically. `(1/n) E log H` over the `4^(n-1)` words is 0.386, 0.403, 0.416, 0.427,
0.436, 0.443 for `n = 8..13` (still rising; the Lyapunov exponent of the random walk is above 0.44). Proximality
(E): the two largest eigenvalue moduli are 1.839/0.737 (`1-`), 2.288/1.255 (`-2`), 2.563/1.346 (`12`), and 3.000/0.000
(`1-.-2`), 4.745/0.604 (`1-.12`), 5.236/1.576 (`-2.12`): every letter and two-letter product is proximal. With the
irreducibility of THM-4880 (G) these are the standing hypotheses of the limit theory of products of random matrices
(Guivarc'h-Le Page; Benoist-Quint, *Random walks on reductive groups*, Chapters 10-17: CITED as the frame; strong
irreducibility of the semigroup is not verified here and no theorem is applied).

## 7. The collective invariant (HYP-9474)

Every candidate "local" invariant is now excluded inside the family: a quotient numeration or a bounded-carry digit
block (THM-4880 (G), (H)), an arithmetic progression (THM-4880 (H2)), a dense response set at each count (section 5),
a dense extension set at each count (THM-4880 note, section 9). What remains, and what the data describe, is a
*distributional* invariant: `r_n(N)` behaves like `4^(n-1)` times a density that is log-normal in `N` around
`exp(n lambda)` (the random-walk picture) times an arithmetic factor in the residues of `N` modulo powers of 3 and 5
(the major arcs). The tiling law is the statement that this density never vanishes at integer scale on the bulk of
the window. HYP-9474 states it, lists the minor-arc decay that a proof would need (a spectral gap for the twisted
transfer operators `f -> sum_a e(theta phi(A_a y)) f(A_a y)` away from the major arcs), and records the numerical
programme that would test it (decay rates in `n` on a frequency grid; the Lyapunov exponent and variance; the major-arc
measure at `a/3^i 5^j`).

## 8. Audit

Independent adversarial audit: requested at the pre-audit commit; outcome and corrections recorded in the close-out
commit (see the theorem's status block).

## 9. Scripts and receipts (`04-computation/experiments/ribbon_spectrum_20261104/`)

reach.py (+ reach_n11.log: projections, cone, post-sink triples, next-level cover), reach2.py (+ reach2_n12.log:
lattice stabilisation, invariance, cumulative multiplicities, rooted triples), bsets.py (+ bsets.log), collective.py
(+ collective_n13.log: exact basis, invariance, multiplicities, residues, proximality), fourier.py (+ fourier_n13.log),
worktree_inventory_20261104.txt (the audit of the 27 worktrees before removal). Theorem: THM-4881; hypothesis: HYP-9474.
