# Classifying the repository's open problems by logical form and by machine size

2026-10-25, klein. **Status:** a classification ledger (CITED external data; the logical-form
assignments are elementary and proved inline where nontrivial), plus a repaired version of the
owner's proposal to classify problems "the way NP-complete problems are" through Turing machines.
External data are from the bbchallenge Cryptids page and the tinycryptids repository read on
2026-10-25 (sources at the end); counts labelled "not independently verified" there are repeated with
that label.

## 0. Inheritance (the repository had most of this already)

The logical-form column below extends, and in several rows merely repeats, three existing
documents: the implication atlas `05-knowledge/results/procgen_atlas_20260924_collatz_implication_atlas.md`
(Collatz, twin primes `Pi^0_2`; no-cycle, Goldbach, RH, LRC `Pi^0_1`; Kurtz–Simon
`Pi^0_2`-completeness and Tao's finite check CITED there), THM-4622
`logical-status-dictionary-absoluteness-choice-and-automorphisms-for-jc2-collatz-hadwiger-nelson`
(`chi_alg <= k` is `Pi^0_1`; `chi(R^2) = chi_alg` in ZFC), the camion–busch note
`camion_busch_gaps_polyhedra_collatz_20261002.md` (`{7,21}`-completeness is `Pi_1`; Collatz is `Pi_2`
and its generalization `Pi_2`-complete, so no uniform reduction), and the reflection
`07-reflections/lrc-collatz-famous-problems-not-pvsnp-but-the-cancellation-family-s599.md`, which
already rejected the P-versus-NP framing. New here: the busy-beaver size column with its October
2026 values, the three-band reading, and the identification of the only tournament–halting contact
in canon: the CITED undecidability of tournament density inequalities (TANGENTS T068,
Chen–Lin–Ma–Wei), which is where "machines built from tournament substructures" is an actual theorem
rather than a metaphor — tournament densities can encode halting; the converse composition is not a
construction anyone has recorded.

## 1. The repair: machines give two coordinates, not a completeness class

NP-completeness is mutual polynomial-time reducibility inside one class. Encoding a conjecture as
"machine `M` halts iff the conjecture is false" places the conjecture in `Pi^0_1` (its negation is
the halting of `M`, a `Sigma^0_1` event). Every `Pi^0_1` sentence is many-one reducible to every
`Pi^0_1`-complete one, so inside `Pi^0_1` all such conjectures are "equivalent" in the trivial sense
and in no useful sense: no reduction between Goldbach and the Riemann hypothesis is known, and none
is implied by the two machines. What the machine encoding *does* give is two exact coordinates:

1. the **level** in the arithmetical hierarchy (`Pi^0_1`: a counterexample is a finite certificate;
   `Pi^0_2`: "for all n there exists k", which halting cannot encode and *quasihalting* (beeping busy
   beaver) can; `Sigma^0_2`, …), and
2. the **descriptive size** (states and symbols of the smallest known machine), a Kolmogorov-type
   measure that the Cryptids page organizes as `BB(n,m)`.

A third, repository-specific coordinate is whether the problem is already **decidable in principle**
(a finite but infeasible check), which several of the repository's anchors are. "Fractal composition
from tournament substructures" is not a mechanism found in the machines; what the small cryptids
share is the **parity/valuation shift** (Collatz-like affine maps selected by residue), which the
repository studies in its Collatz lanes and connects to tournaments only through THM-4754's
3-adic Sierpiński tournament and the score lattices of THM-4752. That is the honest substructure.

## 2. External ledger (bbchallenge Cryptids, Oct 2026)

| statement | level | smallest known machine | note |
|---|---|---|---|
| Goldbach | `Pi^0_1` | BB(25) (2016, Lean-verified); BB(21) (DeepMind, Oct 2026, no explicit statement) | the "25-state machine isomorphic to Goldbach" |
| Riemann hypothesis | `Pi^0_1` (via a decidable inequality) | BB(744) (Matiyasevich–O'Rear 2016); BB(120), BB(95) entries (2026) | also MBB(129) |
| Sloane ternary (`2^n` has a `0` in base 3 for `n > 15`) | `Pi^0_1` | BB(6,4) (Oct 2026) | the fifteen of the exception set `{0,1,2,3,4,15}`; checked to `5.9·10^21` (Saye 2022) |
| Erdős #406 ternary (`2^n` has a `2` in base 3 for `n > 8`) | `Pi^0_1` | BB(15) (Stérin–Woods 2021) | fifteen states |
| weak Collatz (no nontrivial cycle) | `Pi^0_1` | BB(14) (Reams, Oct 2026; tinycryptids); earlier BB(124) | repo: THM-4591 (no cycle of period ≤ 301,993), THM-4633 |
| full Collatz | `Pi^0_2` | BBB(17) quasihalting (not independently verified) | halting cannot encode it |
| twin primes | `Pi^0_2` | BBB(18) quasihalting (not independently verified) | |
| odd perfect number exists | `Sigma^0_1` | BB(37) | non-existence is `Pi^0_1` |
| odd multiply perfect `> 1` exists | `Sigma^0_1` | BB(23,3) (kg583, 2024) | see the multiply-perfect ledger |
| Antihydra / Hydra / Bigfoot never halt | `Pi^0_1` | BB(6) / BB(2,5) / BB(3,3) | "probviously" non-halting; this session's note |
| smallest Sierpiński / Riesel number | `Sigma^0_1` (finitely many candidates remain) | BB(60) / BB(68) | |
| Fermat's last theorem | `Pi^0_1` | BB(400) claimed (1965), BB(90) compiled | a former cryptid, now a theorem |
| PA / ZF / ZFC consistency | `Pi^0_1` | BB(372) / BB(432) / BB(215) | ZFC entry not independently verified |

## 3. The repository's problems

| problem (PROBLEM-LEDGER / frontier) | logical form | level | remark |
|---|---|---|---|
| LRC(14) | for all 13 integer speeds, some time is `1/14`-lonely | **decidable in principle** (`Delta^0_1`): for integer speeds the good-time set is a finite union of intervals, and Tao's reduction (CITED, "Some remarks on the lonely runner conjecture") bounds the speeds that must be checked by `n^{O(n^2)}` | the repository's finite-height/finite-box theorems (THM-763, THM-1290, THM-2100's row box) are this reduction made partly explicit; a machine for it halts |
| GMC(2) / NC2 | proved (THM-2022) | — | was `Pi^0_1`-type: a nonzero constant term for every member |
| JC(2), DC(2) | every Keller map of `C^2` is invertible | `Pi^0_1`: for each degree the statement is first-order in `ACF_0`, decidable by Tarski; the quantifier over degrees is the only unbounded one | the repository's degree-by-degree closures (18, 22, …) are the inner decidable checks |
| Mahler 3/2 (no `Z`-number) | `not exists xi` with all `{xi (3/2)^n} < 1/2` | `Pi^0_2` after Mahler's reduction to integer orbits of `ceil(3x/2)` (no `g_0` has an all-small orbit: for all `g_0` exists `n`) | THM-4469 ties it to Collatz block pairs, also `Pi^0_2` |
| Collatz | for all `n` exists `k` | `Pi^0_2` | THM-4499/4581 lanes |
| Hadwiger–Nelson `chi = 7` | exists a finite unit-distance graph with `chi >= 7` (`Sigma^0_1`, de Bruijn–Erdős; unit-distance realizability is decidable by Tarski) and every finite one is 7-colourable (`Pi^0_1`) | `Sigma^0_1 ∧ Pi^0_1` | currently `5 <= chi <= 7` |
| integer `3x3` magic square of distinct squares | exists | `Sigma^0_1`; non-existence `Pi^0_1` | Sessions 57–61 |
| Sun 2-4-6-8 | for all `n` a representation exists (bounded search) | was `Pi^0_1`; **refuted** by a `Sigma^0_1` witness (THM-4026) | the one row where the certificate shape did its job |
| `H >= disc` for tournaments | for all tournaments | `Pi^0_1` | THM-1950 |
| black-core connectivity of the metagraph | for all `n` | `Pi^0_1` | ports packets |
| Rule 30 centre column not eventually periodic | not `(exists p, N for all n > N ...)` | `Pi^0_2` | THM-3778 |
| ABC (every `epsilon`) | for all `epsilon` exists `C` for all triples | `Pi^0_3` on its face; effective forms `Pi^0_1` | contested IUT |
| Euclidean Kakeya in `R^4` | all Besicovitch sets have full dimension | not arithmetical as stated (quantifies over subsets of `R^4`) | finite-field forms are theorems |
| AMM 12592 biased coin `C* < 3/2` | an infimum over strategies is below a constant | real-analytic; `Pi^0_2`-type after discretization | THM-4467/4494 bounds |

The chart has three bands: the decidable anchor (LRC(14)), the `Pi^0_1` conjectures (JC(2), `H >= disc`,
black-core, Goldbach, RH, weak Collatz, Sloane/Erdős ternary, odd-perfect non-existence), and the
`Pi^0_2` dynamical problems (Collatz, Mahler, Rule 30). The second band is where a small machine is a
faithful encoding; in the third only quasihalting encodes the statement, and in the first a machine
would halt — the "cryptid" status there measures feasibility, not logic.

## 4. Where tournaments actually enter the machines

- The small cryptids are affine maps selected by parity with a counter (Antihydra: `floor(3a/2)`
  and `O > 2E`; the BB(7) Antihydra-like machine: `2n -> 3n`, `2n+1 -> 3n+5`, same counter rule).
  Their state is a valuation sequence; the repository's Collatz, Mahler and Antihydra notes read
  that sequence as a full shift with 2-adic-ball cylinders.
- THM-4754 realizes a Collatz arc-flip law on the Cayley tournament of `Z/3^k`; THM-4752 lifts
  tournament scores into `II_{25,1}`. These are the only tournament-to-dynamics contacts in canon
  that pass the connection contract; neither composes machines from tournaments.
- The repository's own counting machines (the ordinary A000568 recurrence, the clone-excess
  recurrence) are terminating computations, hence level zero; their "hardness" is cost, measured by
  state graphs (13,445 states at order 100), the same kind of coordinate as a busy-beaver size.

## Sources

- bbchallenge wiki, Cryptids: https://wiki.bbchallenge.org/wiki/Cryptids (read 2026-10-25).
- tinycryptids, weak Collatz machines: https://github.com/calr3/tinycryptids/tree/main/weakcollatz .
- Goldbach 25-state machine: Wikipedia "Busy beaver" and the 2025 BB(5) determination paper
  (https://arxiv.org/pdf/2509.12337); the 2016 Yedidia–Aaronson paper (https://arxiv.org/pdf/1605.04343)
  for the 4,888-state original and the 27/31/47-state code-golf lineage.
- Sloane's ternary conjecture: OEIS A004642, A102483, A346497; Saye, J. Integer Seq. 25 (2022)
  (https://cs.uwaterloo.ca/journals/JIS/VOL25/Saye/saye3.html; https://arxiv.org/pdf/2202.13256).
- Lonely runner decidability for fixed `n`: T. Tao, "Some remarks on the lonely runner conjecture"
  (Contrib. Discrete Math. 2018), CITED from memory of the statement; verify before relying on the bound.
