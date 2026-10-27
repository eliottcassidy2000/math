# Order-11 prime cycle-index types by symmetry-first enumeration, the order-12 ledger, and a census of the dipole law's potentials (λ, δ⁰, L, the fractional dipole LP) on all classes of orders 4–8

**Session:** klein-2026-10-27 (worktree `codex/session-overnight-20261027`), 2026-10-27, continuing
[order10_prime_types_n11_ledger_spectral_blackcore_20261027/](../order10_prime_types_n11_ledger_spectral_blackcore_20261027/README.md)
(THM-4814, HYP-9422). Canon: [THM-4818](../../../01-canon/theorems/THM-4818-anti-automorphism-cycles-are-2-mod-4-order-11-prime-types-by-symmetry-order-12-ledger-dipole-census-condensation-criterion-and-p7-class-distance-six.md);
hypotheses: [HYP-9424](../../hypotheses/HYP-9424-singleton-distance-is-delta0-iff-an-eligible-endpoint-exists-linear-potentials-never-exact-for-the-transitive-target-and-the-order-12-types.md);
companion note on the three external preprints: [external_preprints_mub6_local_smoothing_eth_merge_20261027.md](../external_preprints_mub6_local_smoothing_eth_merge_20261027.md).

**Status.** FINITE-EXACT: the order-11 type census (two exact global validations), the order-12 ledger given four OEIS
inputs, the class move graphs and all distances of orders 4–8, the labelled BFS at order 7, the all-regular order-11
black witnesses. PROVED (elementary): the 2 mod 4 law for anti-automorphism cycles; the condensation criterion
(one direction); the root-word reading of the score bound. Not independently audited.

**Concurrency.** mac-mini-2026-10-27-dipole reserved THM-4815 (dipole principle), THM-4816 (order-11 types) and
THM-4817 (order-12 ledger) before this packet was started. The ledger here is filed as an **independent cross-check**
(different method: invariant enumeration of symmetric primes rather than a full census) and takes no theorem id for
it; the numbers were sent to mac-mini (agents/mac-mini/inbox, 2026-10-27). THM-4818 covers only what is not in those
reservations: the anti-cycle law, the census of potentials, the criterion, and two corrections.

## 1. The order-11 prime cycle-index types, symmetry first (7 seconds)

**Why no full census is needed.** A prime tournament of order 11 with a nontrivial automorphism has one of odd prime
order `p`; a lone `p`-cycle with the rest fixed is a nontrivial module (THM-4814's lemma), so the cycle type is
`3^2 1^5`, `3^3 1^2`, `5^2 1` or the single `11`-cycle. Every self-converse tournament has an anti-automorphism `ψ`,
and the **2 mod 4 law** (new here, two lines): a pair `{u, v}` swapped by `ψ^k` forces `u → v ⟹ v → u` when `k` is
even, so every cycle of `ψ` has length `2m` with `m` odd (plus at most one fixed point). THM-4814 had only "all cycles
even"; the types `4 2^3 1`, `4^2 2 1`, `6 4 1`, `8 2 1` never occur, and `scripts/anti_types_check.py` confirms the
law for every partition of every `n ≤ 12` (`data/anti_types_check.out`). So the symmetric primes of order 11 are
among the invariant tournaments of seven permutations, one orientation bit per orbit on pairs, with two reductions
(the free five-point block taken from the 12 classes of 5-tournaments; 2-cycles of `ψ` oriented `p → q`):

| type | invariant tournaments (candidates) |
|---|---:|
| `3^2 1^5` | 12 · 2^15 = 393,216 |
| `3^3 1^2` | 2^19 = 524,288 |
| `5^2 1` | 2^11 = 2,048 |
| `11` | 2^5 = 32 |
| anti `2^5 1` | 12 · 2^15 = 393,216 |
| anti `6 2^2 1` | 2^12 = 4,096 (of 2^14) |
| anti `10 1` | 2^6 = 64 |

`scripts/tourn_types_dedup.c` (THM-4814's census program plus a canonical-form hash set) keeps the strong, prime
candidates, deduplicates them by nauty canonical form and records |Aut|, the automorphism cycle types and the
anti-automorphism cycle types: 1,316,960 candidates, 1,111,416 strong, **270,999 distinct symmetric prime classes**,
of which 259,307 self-converse (`data/s11_sym.txt`, `data/t11_sym.txt`). With the rigid non-self-converse remainder
`P(11) − 270,999 = 796,605,430` the twelve types are (`data/prime_types_n11.txt`):

| classes | |Aut| | automorphism types | SC | anti-automorphism types |
|---:|---:|---|:-:|---|
| 796,605,430 | 1 | `1^11` | no | — |
| 258,844 | 1 | `1^11` | yes | `2^5 1` ×1 |
| 10,456 | 3 | `3^2 1^5` ×2 | no | — |
| 176 | 3 | `3^2 1^5` ×2 | yes | `2^5 1` ×3 |
| 176 | 3 | `3^2 1^5` ×2 | yes | `6 2^2 1` ×2, `2^5 1` ×1 |
| 1,224 | 3 | `3^3 1^2` ×2 | no | — |
| 96 | 3 | `3^3 1^2` ×2 | yes | `2^5 1` ×3 |
| 12 | 5 | `5^2 1` ×4 | no | — |
| 6 | 5 | `5^2 1` ×4 | yes | `2^5 1` ×5 |
| 5 | 5 | `5^2 1` ×4 | yes | `10 1` ×4, `2^5 1` ×1 |
| 3 | 11 | `11` ×10 | yes | `2^5 1` ×11 |
| 1 (P_11) | 55 | `11` ×10, `5^2 1` ×44 | yes | `10 1` ×44, `2^5 1` ×11 |

|Aut| histogram of the 796,876,429 primes: 1: 796,864,274; 3: 12,128; 5: 23; 11: 3; 55: 1. The four classes with an
11-cycle are the circulant tournaments on `Z/11` (Burnside over the multiplier group `(Z/11)^×` of order 10:
`(32 + 4·2)/10 = 4`), all self-converse (`x ↦ −x`).

**Two exact validations.** (i) The self-converse primes number 259,307, exactly the value forced by OEIS A002785(11)
= 279,968 through THM-4814's self-converse decomposition identity (SC non-strong 2,768 + SC strong non-prime 17,893,
`data/species12b.out`, which now also reproduces the census at `n = 11`: `total=279968 census=279968 OK`).
(ii) `Σ_{non-rigid primes} (1 − 1/|Aut|) = 1337728/165 = P(11) − P_lab(11)/11!` exactly, with `P_lab(11) =
31,808,433,417,949,440` from the labelled species identity. Any miscount of a non-rigid class or of its group order
would break (ii); any miscount of a self-converse class would break (i).

## 2. The order-12 ledger

The species identity (THM-4803/THM-4814) with the types through order 11 determines the entire non-prime side at
order 12 (`data/species11_fast_12.out`, `data/species12b.out`; `species11_fast.py` replaces THM-4814's `euler_mset`,
which looped `c_s(11) ≈ 4.4·10^8` times, by the Euler-transform recurrence). Inputs taken from OEIS on 2026-10-27:
A000568(12) = 154,108,311,168; A002785(12) = 1,492,288; A051337(12) = 152,310,149,735; A054946(12) =
72,926,863,133,112,198,144.

| quantity at order 12 | value | depends on |
|---|---:|---|
| all classes without the order-12 primes | 11,127,875,377 | types ≤ 11 |
| strong classes without the order-12 primes (= strong classes with a nontrivial module) | 9,329,713,944 | types ≤ 11 |
| non-strong classes | 1,798,161,433 | types ≤ 11; **= A000568(12) − A051337(12)**, an independent check |
| prime tournaments P(12) (A280772, OEIS stops at 8) | **142,980,435,791** | A000568(12) (and A051337(12) agrees) |
| self-converse non-strong / strong non-prime | 8,839 / 312,500 | types ≤ 11 |
| self-converse primes | **1,170,949** | A002785(12) |
| labelled strong s_lab(12) | 72,926,863,133,112,198,144 | none; = A054946(12) |
| labelled strong non-prime | 4,439,096,863,965,333,504 | types ≤ 11 |
| labelled primes P_lab(12) | **68,487,766,269,146,865,664** | none |
| strong intransitivity classes c_s(12) | 4,663,947,439 + (P(12) + SCP(12))/2 = **76,154,750,809** | A000568, A002785 |
| intransitivity classes c(12) | 5,111,156,796 + 71,490,803,370 = **76,601,960,166** | A000568, A002785 |
| labelled intransitivity structures ℓ(12) | 1,281,079,575,450,125,699 + P_lab(12)/2 = **35,524,962,710,023,558,531** | none |

(The order-12 primes enter `c_s` and `c` only through their own `x^12` term, one per self-converse prime and one per
converse pair, hence `(P + SCP)/2`; they enter `ℓ` through `P_lab/2`.) At order 11 the same run returns the values of
THM-4814 (`c(11) = 447,114,461`, `c_s(11) = 442,336,617`, `ℓ(11) = 16,901,802,400,214,078`), now from the complete
type files rather than from the OEIS back-substitution.

## 3. The dipole law's potentials on every class of orders 4–8

`scripts/movegraph.c` builds the native-move graph on isomorphism classes: for every class `T` and every directed
Hamiltonian path `Q` of `T`, `τ_Q T` (keep `Q`, reverse every other pair) canonicalised with nauty, with `Q` marked blue
when the positional reversal is an anti-automorphism (`data/mg4.txt`–`mg7.txt`; `mg8.txt`, 1.8 M edge lines, is
regenerated in seconds). The move is an involution, so the class graph is undirected and BFS from a locus gives exact
distances. `scripts/analyze_mg.py` compares them with the potentials of THM-4811 (A): `λ` (arc connectivity, target
the non-strong locus), `δ⁰ = min_v min(s(v), n−1−s(v))` (target: a vertex of score 0 or n−1), and THM-4810's
`L = ½ Σ |s_(i) − i|` (target: the transitive class). (`data/an5.txt`–`an8.txt`.)

| n | classes | `λ`, both colours | `λ`, black only | `δ⁰` singleton, both | `δ⁰`, black | `L` to TT exact | L-local minima | max dist to TT |
|---:|---:|---|---|---|---|---|---|---:|
| 5 | 12 | exact 12/12 | exact 11/12 (R5 unreachable) | exact | exact | 6/12 | — | 3 |
| 6 | 56 | exact | exact | 55/56 | 55/56 | 23/56 | 22 | 4 |
| 7 | 456 | exact | exact | 453/456 | 453/456 | 274/456 | 110 | 6 |
| 8 | 6880 | exact | exact | 6864/6880 | 6864/6880 | 5449/6880 | 904 | 7 |

1. **`λ` is exact everywhere**, both colours and black-only (THM-4811 (B), (C) reproduced on classes).
2. **Singleton distance and the condensation criterion.** The 1, 3, 16 failures (distance 2, `δ⁰ = 1`) are all
   non-strong, and **every one is explained by the condensation** (`scripts/singleton_failures.py`,
   `data/singleton_failures_n6_8.out`): a legal word starts in the initial strong component and ends in the terminal
   one; by the dipole law a vertex `w` gains `|d(w)|` only as the start when `s(w) ≤ (n−1)/2` or as the end when
   `s(w) ≥ (n−1)/2`; so if no vertex attaining `δ⁰` is *eligible* (small score in the initial component, or large
   score in the terminal one) the first move cannot progress and the distance exceeds `δ⁰` (PROVED). Conversely, at
   `n ≤ 8` every class with an eligible vertex is exactly `δ⁰` from a singleton (FINITE-EXACT; the strong case is
   THM-4811 (D)). One failure at `n = 8` (scores `6,6,5,5,1,2,2,1`, components 4+4) is rigid and not self-converse:
   the mechanism is forced structure, not forced symmetry.
3. **The transitive target is different in kind.** `L` is exact on a shrinking majority, its local minima are
   generic (at `n = 8`: 748 of the 904 are rigid non-self-converse non-strong classes, 42 rigid strong), and the
   distance exceeds `L` by up to 3 already at `n = 7`. No symmetry explains this: the transitive class is a single
   point, and the unique Hamiltonian path of `TT_n` makes its move blue, so `TT_n` has exactly one neighbour and the
   black-only distance to it is infinite for every other class.
4. **Linear potentials cannot be exact for the transitive target.** `scripts/lp_dipole.py` computes the best bound
   obtainable from *any* linear negation-invariant potential: a word of length `ℓ` gives
   `(K' − (−1)^ℓ K)/2 = Σ_i (−1)^{ℓ−i} A_{Q_i}` with `A_Q` the signed matrix of a directed Hamiltonian path, and
   relaxing both the integrality and the requirement that `Q_i` be a path *of the current tournament* gives the LP
   `min Σ c_Q` over all `n!` directed paths (its dual variables are exactly the linear potentials; the LP dominates
   every score, cut and sorted-score bound, and the trace-norm tariff is dominated by the same relaxation). At
   `n = 5` the LP is exact except on the three classes with `|Aut| = 3`; at `n = 6` it has gaps on rigid
   non-self-converse classes (`data/lp5.txt`, `lp6.txt`). The binding constraint is that the kept path must be
   directed in the current tournament, an orientation condition no linear functional sees.
5. **Blue moves are not loops on classes** (32 of 32 at `n = 5`, 88 of 96 at 6, 780 of 780 at 7, 4,500 of 4,512 at
   8 change the class); the few blue loops are THM-4729's.
6. **P_7 is six moves from the transitive class, and `L` is exact there.** The class BFS gives 6 = `L(P_7)` for all
   three regular 7-tournaments; `scripts/lbfs7.c` repeats this on all 2^21 labelled tournaments: from the standard
   `P_7` the 5,040 transitive labellings split 2,520 at distance 6 and 2,520 at distance 7, and the standard
   `i < j` labelling is at distance 7 (`data/lbfs7_paley.out`; eccentricity 7, the graph is connected, the
   eccentricity of `TT_7` is 13 as in THM-4810). THM-4810's "P_7 at distance 7 from the transitive tournament" is
   the distance to one fixed labelling; the reflection table's "P_7 no (7 > 6)" is therefore a labelled/class
   conflation (MISTAKE-682). The root-word reading of `L` and the two-sided score bound for arbitrary targets are in
   section 4 of the external note.
7. **All 1,223 regular 11-tournaments** have an all-black word of exactly `m = 5` moves to a singleton plus a strong
   core (`data/bc11all.out`, THM-4814's `blackcore.py` on `gentourng -d5D5 11`, 0.4 s per class): HYP-9422 (1) is
   FINITE-EXACT through order 11 and is the regular case of THM-4811 (D2) at the first order beyond its census.

## 4. Corrections and sharpenings recorded

- THM-4814's anti-automorphism lemma is sharpened to the 2 mod 4 law (update line in THM-4814; no error).
- THM-4810: class distance from `P_7` to the transitive class is 6, not 7 (update line; the labelled statement stands).
- The reflection `why-the-advances-work-ii` row "distance to the transitive class": `L` is exact at `P_7`
  (MISTAKE-682, filed against the inference, not against THM-4810's labelled number).
- HYP-9422 (1) at order 11: closed; HYP-9422 (4): resolved by section 1.

## 5. Reproduction

`scripts/run.sh` (timings there). Receipts: every Python script was run under `python3` and `python3 -O`; outputs are
byte-identical (`receipts/`). `SHA256SUMS` covers every file of the packet. Omitted from `data/` for size: the seven
candidate files (75 MB) and `mg8.txt` (40 MB), both regenerated by `run.sh` in under a minute.
