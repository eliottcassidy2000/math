# Five October-2026 papers read against the Collatz work: a dimension-free spectral gap for the loosened Collatz walk modulo `3^k` (PROVED, the KLS theme realised exactly), the stopping-time route to the last-dip residual ("wanted entries only"), Mumford–Shah blow-ups as degree-`≤3` trees, Golod–Shafarevich towers as cylinder towers, and Miller's non-diagonalizable sequences as measurement independence

**Session:** mac-mini (claude) 2026-10-06, worktree `codex/session-trees-tournaments-20261005`, fourth wave.
**Owner's seed:** merge ideas related to arXiv:2609.26732v2, 2610.06461v1, 2610.01447v2, 2404.01256v3, 2610.06783v1 into
the exploration, leveraging analogous ideas and themes toward Collatz.

**The five papers (all read 2026-10-06 through their HTML abstracts/theorem statements; none mentions Collatz).**

| paper | result | mechanism borrowed here |
|---|---|---|
| Deangelis, *Solution of the Mumford–Shah conjecture* (2609.26732v2) | blow-up limits of planar minimizers: the discontinuity set is locally `m ∈ {1,2,3}` `C^1` arcs — an endpoint, a smooth arc, or a 120° triple junction; endpoints and triple junctions locally finite; no bounded component (compact translation rigidity) | local models are **trees of maximum degree 3**; a degree-4 crossing is never minimal |
| Keller–Pagano, *Golod–Shafarevich for arithmetic surfaces* (2610.06461v1) | towers `S_{n+1} -> S_n` of étale double covers with a section splitting completely at every level; `g(D_n) - 1 = 2^n (g-1)`; infinitude from the Golod–Shafarevich inequality `dim H^2 > g^2/4` with obstruction codimension `≤ 5000 g^2/log g` | a **pro-2 tower with a completely splitting section** = a 2-adic point lying in every cylinder |
| Song–Zhang, *An `O(1)` bound for the KLS constant* (2610.01447v2) | universal Poincaré/Cheeger constant for isotropic log-concave measures via stochastic localisation, Appell polynomials, iterated logarithmic refinement, repeated height reduction | a **dimension-free spectral gap** |
| Bauer–Hanson, *The countable reals* (2404.01256v3) | a topos where the Dedekind reals are countable, built from Miller's non-diagonalizable sequences (every real computable from the sequence already appears in it); Cauchy reals stay uncountable; countable choice fails | **a bank that already contains everything derivable from it**; two notions of "all" that split |
| Alman–Vassilevska Williams, *Truly subquadratic 3SUM and truly subcubic APSP via triangles in sparse lopsided graphs* (2610.06783v1) | thin matrix products via Schönhage's identity, visiting only the recursion leaves that feed the **wanted entries**; lopsided tripartite triangle problems | **compute only the entries that can matter** |

**Status.** PROVED (elementary): Proposition 3.1 (exact spectrum of the loosened Collatz walk mod `3^k`), Proposition 5.1
(stopping-time route). FINITE-EXACT: spectral data `k ≤ 9`; stopping-time records to `2^33` (§5.1: the last-dip lemma holds for all
`n < 2^33` without the CST input); the `Δ ≤ 3` tree census. CITED: the five papers as listed; THM-4520 (level operator with spectrum on
the half-circle); Chung–Diaconis–Graham 1987. Everything else is typed ANALOGY or DICTIONARY in §7. Collatz OPEN.
Scripts: `04-computation/experiments/collatz_loosened_graph_mod3k_spectral_gap_20261006.py`,
`collatz_stopping_time_records_20261006.c`; outputs in this directory. Audit OWED on Proposition 3.1's proof.

## 0. Inheritance

- E-SCC / HYP-9120: the Le–Smith loosened Collatz graph `E` (arcs `n -> 3n+1` for every `n`, `n -> n/2` for even `n`);
  Q2 = `1` reaches every `m` prime to 3; reduced to the classes `1, 14 (mod 27)` (Lean), verified to `10^18`
  (`collatz_procgen_20260922_q2_endgame.md`). The mod-6 synthesis has the greedy 3-adic map as an exact Markov chain mod
  `3^(J+1)` with drift `log(2/3)`.
- THM-4520: the level-`n` Syracuse frequency operator is a Gauss-twisted circulant on `Z/L_n`, `L_n = 2·3^(n-1)`, with
  characteristic polynomial `(2^L - 1)λ^L + λ^(L/2) - 1`, every eigenvalue of modulus `1/2`; "Parseval rate = non-normal
  transient".
- The last-dip lemma and its residual (`collatz_last_dip_cst_alphabet_exotic_20261005.md` §1); THM-4512 (coefficient
  descent); the delay/path records to `10^12` used as test inputs in `procgen_landing_20260926` (A006877, A006884).
- The first note's trees `A_5, D_5, D̃_4` and the cubic hosts; the Sumner degree wall.

## 1. Mumford–Shah: minimal singular sets are trees of degree `≤ 3`

Deangelis's classification says a minimizer's discontinuity set is locally one of three models — a crack tip (one arc),
a smooth arc (two opposite arcs), or three arcs at 120° — and never a crossing. In graph language the local vertex types
of a minimal singular set are exactly the vertex types of the **fork** `S(2,1,1)` (degrees 1, 2, 3), while the **star**
`K_{1,4}` is the forbidden crossing (four arcs at a point). This is the Steiner-tree fact (a crossing is beaten by two
120° junctions) and it is the same dichotomy that organised the earlier notes: at `t = 5` the star is the one tree that
sees the Erdős–Sós half (the degree wall `3n/2`) and the Sumner wall (regular `7`-tournaments), the cubic hosts `K_{3,3}`,
Petersen and the Collatz functional graph are star-free, and the affine `D̃_4` is the star. Census (FINITE-EXACT,
A000672): trees with maximum degree `≤ 3` on `n = 4..12` vertices number `2, 2, 4, 6, 11, 18, 37, 66, 135` of
`2, 3, 6, 11, 23, 47, 106, 235, 551`; the Mumford–Shah-admissible topologies are a vanishing fraction, and the three
5-vertex trees split `2 + 1`. Compact translation rigidity ("no bounded component") has the trivial Collatz twin that no
component of the functional graph is finite (every `n` has the preimages `2^j n`). Typed: ANALOGY with an exact
combinatorial core (vertex types `{1,2,3}`); no transfer of the variational mechanism.

## 2. Golod–Shafarevich: towers with a completely splitting section

Keller–Pagano build `S_{n+1} -> S_n` with `[S_n : S] = 2^n` and a section meeting every level, and prove infinitude by a
counting inequality (`dim H^2 > g^2/4` against bounded obstructions). The Collatz twin is the cylinder tower
`Z/2^(k+1) -> Z/2^k` of parity words: the `T`-word of length `k` is determined by `x mod 2^k`, each level is a double
cover of the previous one, and a word `w` with its 2-adic realizer `x_w ∈ Z_2` is a **section that splits completely**
— `x_w` lies in every cylinder of its prefixes. The Golod–Shafarevich move "few relations relative to generators ⟹ the
tower never stops" is the entropy count of alive cylinders: the no-descent classes at level `k` number `Θ(2^{hk} k^{-3/2})`
with `h = h(log_3 2) = 0.95` (THM-4495, THM-4476), so the tower of undecided cylinders is infinite for a counting reason,
exactly as a GS tower is. What the analogy makes vivid is the obstruction the repo keeps meeting: the arithmetic-surface
tower has a *rational* section by construction, while the Collatz tower always has a *2-adic* section (every word is
realised) and the question is whether an *integer* lies on it — a property no level of the tower sees (S7's measurement
independence; the shadow theorem's "integrality is the class condition"). The genus-doubling identity
`g(D_n) - 1 = 2^n (g - 1)` has as its Collatz twin Terras's exact cylinder law `#{odd n < 2^K : word w} = 2^{K-A-1}`.
Typed: ANALOGY; the two towers share the pro-2 structure and the counting criterion, not the splitting mechanism.

## 3. KLS: a dimension-free spectral gap, realised exactly on the loosened Collatz walk mod `3^k`

Song–Zhang prove a universal constant independent of the dimension. Asked of the Collatz objects: does the loosened
Collatz graph have a dimension-free gap on its 3-adic quotients? Modulo `3^k` both generators of `E` descend —
`f(x) = 3x+1` and `h(x) = x·2^{-1}` (2 is a primitive root mod `3^k`, so `h` is a single cycle of length
`L_k = 2·3^{k-1}` on the units `U_k`; `f` maps `U_k` 3-to-1 onto the class `1 (mod 3)`). Let `P_k` be the random walk
on `U_k` that applies `f` or `h` with probability `1/2`.

**Proposition 3.1 (PROVED).** The spectrum of `P_k` is
`{1, -1/2} ∪ { e^{2πi j/L_m}/2 : 2 ≤ m ≤ k, 0 ≤ j < L_m, 3 ∤ j }`
(with multiplicity; `1 + 1 + Σ_{m=2}^{k} 4·3^{m-2} = L_k`). In particular every non-trivial eigenvalue has modulus
exactly `1/2`: the directed walk has a **dimension-free spectral gap `1/2`**. Its characteristic polynomial is
`(λ - 1)(λ + 1/2) Π_{m=2}^{k} Φ_{3^{m-1}}((2λ)^2)` (since `Π_{3∤j}(x - ζ_{L_m}^j) = (x^{L_m}-1)/(x^{L_m/3}-1) = Φ_{3^{m-1}}(x^2)`).

*Proof.* Reduction mod `3^{k-1}` commutes with `f` and `h`, so `P_k` projects onto `P_{k-1}`; the lifts of `P_{k-1}`'s
eigenmeasures are eigenmeasures of `P_k`. The complementary invariant subspace `W_k` is the space of measures with zero
sum on every fibre `{x, x + 3^{k-1}, x + 2·3^{k-1}}` of `U_k -> U_{k-1}` (dimension `L_k - L_{k-1} = 4·3^{k-2}`). The
`f`-part of `μ P_k` at `y` is `(1/2) Σ_{f(x)=y} μ(x)`, a sum over a whole fibre (all three fibre elements have the same
image `3x + 1 mod 3^k`), hence zero for `μ ∈ W_k`; the `h`-part is `(1/2) μ(h^{-1} y)`, and `h` permutes fibres, so `W_k`
is `h`-invariant. Thus `P_k|_{W_k} = (1/2) H|_{W_k}` with `H` the permutation matrix of the single `L_k`-cycle `h`. The
eigenvectors of `H` are the characters `χ_j(x) = e^{2πi j·pos(x)/L_k}`; `H^{L_{k-1}}` acts within fibres as a 3-cycle
(because `2^{L_{k-1}} ≡ 1 (mod 3^{k-1})` but not mod `3^k`), so `χ_j` has zero fibre sums iff `χ_j^{L_{k-1}} ≠ 1` iff
`3 ∤ j`. Hence the new eigenvalues at level `k` are `e^{2πi j/L_k}/2`, `3 ∤ j`. The base `U_1 = {1, 2}` has
`P_1 = [[1/2, 1/2],[1, 0]]` with eigenvalues `1, -1/2`. ∎ (FINITE-EXACT: the predicted spectrum matches the computed
one to `10^{-11}` for `k = 2..8`; `|λ_2| = 0.500000` throughout; the stationary distribution is far from uniform,
`max π / min π = 8, 28.8, 78.2, 138, 303, 659, 1210` for `k = 2..8`; total variation from a point mass after `4k` steps
`5.8·10^{-3}, …, 1.2·10^{-6}` (`k = 2..9`) and after `8k` steps `≤ 10^{-14}`: mixing in `Θ(k) = Θ(log N)` steps.)

**The undirected contrast.** The normalised-Laplacian gap of the *undirected* Schreier graph decays:
`0.667, 0.333, 0.173, 0.109, 0.080, 0.065, 0.056, 0.048` for `k = 2..9` (ratios `2.00, 1.93, 1.59, 1.36, 1.23, 1.17, 1.16`;
log-log slope `-0.22` against `|U_k|`), so `E_k` is **not** an expander family in the Cheeger sense, while the directed
walk has the exact gap `1/2`. The same split holds for the classical Chung–Diaconis–Graham walk `x -> 2x, 2x+1` on
`Z/3^k`: directed `|λ_2| = 1/2` exactly (`k ≤ 7`; the mechanism is `Π_{a ∈ U_k}(1 + ζ^a) = Φ_{3^k}(-1) = 1`), undirected
gap `0.35 -> 0.032` over `k = 2..8`. This is the KLS theme in its Collatz form: the universal constant lives in the
non-normal directed operator (THM-4520's "Parseval rate = non-normal transient", now with the Gauss twist removed), not
in the geometry of the underlying graph. Typed: PROVED proposition + FINITE-EXACT; the relation to Song–Zhang is ANALOGY
(dimension-freeness), the relation to THM-4520 is EXACT (same cyclotomic mechanism, `Φ_{3^{m-1}}((2λ)^2)` here against
`(2^L-1)λ^L + λ^{L/2} - 1` there).

**What it does and does not give for Q2.** Strong connectivity of `E` on the units mod `3^k` is immediate (one `h`-cycle),
and the walk forgets its start in `O(k)` steps; but Q2 asks about the integer graph, where halving is allowed only at even
numbers, and the mod-`3^k` shadow erases parity entirely. The proposition is a statement about the 3-adic shadow, consistent
with the repo's standing verdict that residue arguments are side-blind. The stationary measure `π_k` (a 3-adic measure
concentrating on `1 (mod 3)` by the folding) is the natural next object: its limit on `Z_3` and its relation to S10's visit
measure `π` of the minimum map are not computed here.

## 4. The countable reals: a bank that already contains everything derivable from it

Miller's non-diagonalizable sequence `μ` has the property that any real computable from `μ` already appears in `μ`;
inside the topos this makes the Dedekind reals countable while the Cauchy reals stay uncountable, and countable choice
fails. Three Collatz twins, all typed ANALOGY:
- **Measurement independence (S7, Codex-corrected form).** For a certified bank `C` of orbits and a source `m` whose orbit
  avoids `C`, every quantity determined by `C` takes the same value whether `m` is rooted or not: "a readout certifies
  exactly its bank". This is the negative form of Miller's property — the bank already contains every integer its
  identities can certify. The repo's refuel-bill note says the same about weights ("the weights cannot be more generous
  than the chain's own code probability").
- **Two notions of "all".** Dedekind-countable vs Cauchy-uncountable is the split between *density-one coverage* (Terras:
  almost all integers have finite stopping time; the renewal notes: every proved family has density zero or one in the
  wrong sense) and *universal coverage* (every integer rooted). The repo already records that under the atomic prior
  "mass-one root coverage is exactly universal coverage" (MISTAKES.md:36) — the one prior where the two notions agree,
  as the two kinds of real agree classically.
- **No countable choice.** The refuel-bill note's "Collatz = totality of a prefix-free machine" and S7's "every
  certificate encodes the orbit" are the statement that no uniform choice of certificates exists short of the orbits
  themselves; the topos where the reals are enumerable but choice fails is the logical shadow of that.

## 5. 3SUM / APSP: compute only the wanted entries

Alman–Vassilevska Williams win by never computing the full product: only the recursion leaves feeding the wanted entries
`W` are visited. The last-dip residual has the same shape. A violation in the window `(l, A, [X, N_l])` needs a source
`m* < N_l` with `σ_T(m*) > A`; everything else about the orbit is irrelevant. So the wanted entry is a single number:

**Proposition 5.1 (PROVED).** If `max { σ_T(m) : m odd, m < N }` is `< A` for every residual length `l` with `N_l ≤ N`,
then the last-dip lemma holds for every `n ≤ N`. *Proof.* A violation `(m*, n)` with `n ≤ N` has `m* < n ≤ N_l` for its
word's length `l ∈ B(X)` and stays above `n > m*` for `A` `T`-steps, so `σ_T(m*) > A`. ∎

FINITE-EXACT: the record sweep (`collatz_stopping_time_records_20261006.c`) gives the stopping-time records
`σ_T(362343) = 165, σ_T(381727) = 173, σ_T(626331) = 176, σ_T(1027431) = 183, σ_T(1126015) = 224, σ_T(8088063) = 246,`
`σ_T(13421671) = 287, σ_T(20638335) = 292, σ_T(26716671) = 298, σ_T(56924955) = 308, σ_T(63728127) = 376,`
`σ_T(217740015) = 395, σ_T(1200991791) = 398, σ_T(1827397567) = 433, σ_T(2788008987) = 447`, so
`max σ_T < 2^32` is `447`, against `A ≥ 4701` (`X = 2^20` windows), `24727` (`2^28`), `125743` (`2^32`): all residual windows
with `N_l ≤ 2^32` are empty by the stopping-time route alone, independently of the direct sweep and of the CST input. The
`2^33` run (window of `l = 79335`, `N_l = 7.2·10^9`) is reported in §5.1. This is exactly how Rozier–Terracol push the CST
verification to `2.8·10^19` (excursion and delay records), and it is the lopsided structure of the problem: two short
lists (sources, targets) against one astronomically long list (hovering words), with the wanted entries living on the
short side.

### 5.1 The `2^33` run

`max { σ_T(m) : m odd, m < 2^33 } = 447` (no new record between `2^32` and `2^33`; `collatz_stopping_time_records_2e33_20261006.out`).
The only residual window at `X = 2^32` with `N_l ≤ 2^33` is `l = 79335` (`A = 125743`, `n ∈ [2^32, 7.2·10^9]`), and
`447 < 125743`, so it is empty. **Consequence (FINITE-EXACT, no CST input):** the last-dip lemma, hence
`D = ⋃ rising cones`, holds for every `n < 2^33` by the direct sweep below `2^32` plus Proposition 5.1 on `[2^32, 2^33)`.
The next window on this route is `l = 190537` (`A = 301994`, `N_l = 9.8·10^11`), which a stopping-time sweep to `10^12` would
close; the published delay records to `10^12` (A006877, maximum delay `1348`-ish Collatz steps) already imply `σ_T < 10^4`
there, but they were not re-verified here and are left as CITED-not-used.

## 6. The two mechanisms that transferred, and the three that did not

Transferred (with proofs here): **dimension-free directed gap** (§3; the KLS theme), **wanted-entries-only verification**
(§5; the fine-grained-complexity theme). Did not transfer beyond typed analogy: the variational blow-up (§1: only its
combinatorial shadow, degree `≤ 3`), the GS counting criterion (§2: the Collatz tower is infinite for the same counting
reason, but infinitude is not the question), the realizability topos (§4: a logical mirror of measurement independence).
Nothing here changes any Collatz status; the honest gain is one exact theorem about the 3-adic shadow of `E` and one
clean verification principle.

## 7. Verdicts

| theme | paper | Collatz counterpart | type |
|---|---|---|---|
| local models are degree-`≤3` trees; crossings excluded | Mumford–Shah | fork vs star; cubic hosts; Erdős–Sós / Sumner degree walls | ANALOGY with exact combinatorial core |
| pro-2 tower with completely splitting section | Golod–Shafarevich | cylinder tower, 2-adic realizers `x_w`; alive-cylinder entropy `h = 0.95` | ANALOGY |
| dimension-free constant | KLS | Proposition 3.1: `|λ| = 1/2` exactly for the loosened walk mod `3^k`; THM-4520 | PROVED (prop.), ANALOGY (theme) |
| undirected geometry vs directed operator | KLS/Cheeger | Cheeger gap decays `~N^{-0.22}`, directed gap exact | FINITE-EXACT |
| a sequence containing all it computes | Bauer–Hanson / Miller | measurement independence; density-one vs universal coverage; no uniform certificate | ANALOGY |
| wanted entries only | Alman–Vassilevska Williams | Proposition 5.1; stopping-time records exclude all residual windows below `2^32` | PROVED + FINITE-EXACT |

## 8. Reproduction

```
python3 04-computation/experiments/collatz_loosened_graph_mod3k_spectral_gap_20261006.py 9 > 05-knowledge/results/collatz_loosened_graph_mod3k_spectral_gap_20261006.out   # ~5 min (dense eigen k <= 8)
cc -O2 -o stoprec 04-computation/experiments/collatz_stopping_time_records_20261006.c && ./stoprec 4294967296 > 05-knowledge/results/collatz_stopping_time_records_2e32_20261006.out   # ~10 min
./stoprec 8589934592 > 05-knowledge/results/collatz_stopping_time_records_2e33_20261006.out   # ~20 min
```

References (CITED). F. Deangelis, arXiv:2609.26732v2 (2026). T. Keller, C. Pagano, arXiv:2610.06461v1 (2026). Z. Song,
X. Zhang, arXiv:2610.01447v2 (2026). A. Bauer, J. E. Hanson, arXiv:2404.01256v3. J. Alman, V. Vassilevska Williams,
arXiv:2610.06783v1 (2026). F. R. K. Chung, P. Diaconis, R. L. Graham, *Random walks arising in random number generation*,
Ann. Probab. 15 (1987). O. Rozier, C. Terracol, arXiv:2502.00948 (2025). OEIS A000672, A006877, A006884.
