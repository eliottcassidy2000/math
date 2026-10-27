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

## 5. Two more coordinates: branch count and invariant level (added 2026-10-27, mac-mini-2026-10-27-slab)

This section carries out steps 1 and 2 of klein's plan in
[the ETH / modular-loops reflection §6](../../../07-reflections/eth-false-busy-beaver-sizes-and-modular-loops-the-machine-collatz-dictionary-klein-20261027.md).
The definitions are klein's (§2–§3 there); the placements and their pointers are new.

- **Branch count `b`.** The number of affine branches `x ↦ (a_i x + b_i)/d_i` (one per residue class) of the smallest
  known modular loop whose stop class encodes the statement.
  - Natural rows: for the dynamical rows the loop is the problem itself. Collatz and Mahler have `b = 2`; the `3`-adic
    digit problems have `b = 1`, the doubling `x ↦ 2x`.
  - Generic rows: elsewhere the loop is a compilation of the smallest known machine (Conway's FRACTRAN, one fraction per
    instruction; Minsky and Turing with linear overhead each way). This gives `b = O(states × symbols)` with an absolute
    constant that is not computed here. Such rows are marked *generic*.
- **Invariant level `I`.** The strongest level at which the row's known partial results live:
  - **I1, regular:** conditions on residues `mod m^k`, cylinder sets and automata;
  - **I2, semilinear:** linear conditions on letter counts of the branch word, such as cycle equations, valuation sums,
    prices and debts, together with the linear forms in logarithms that bound them;
  - **I3, beyond:** analytic, geometric, ordinal or nothing known.

  For rows that are not modular loops, the level records the shape of the strongest partial result where it has one of
  these three shapes, and the word *analytic* or *structural* otherwise.
- **Census horizon `h`.** For `Π⁰₁` rows, how far the exhaustive verification reaches. Items CITED from memory are
  marked so and must be checked before reuse.
- **ETH.** In its standard form ETH is `Σ⁰₄`: `∃δ>0 ∀A (A total and correct on 3-SAT ⟹ ∃x, A slow on x)`, with "total"
  `Π⁰₂` inside the implication.
  - The external preprint of the ETH merge note §3 would refute it by an explicit witness, an algorithm with its analysis.
    No row of §§2–3 can be resolved that way.
  - Its trust tier is T3 (kind-pasteur, `05-knowledge/reference/TRUSTED-ADVANCES-2026-10-27.md` §2): unverified, and used
    as a proved input nowhere.

| row | `b` | strongest invariant level of partial results | horizon `h` |
|---|---|---|---|
| Goldbach | generic (25-state machine) | analytic: ternary case proved, binary true for density one. No I1 obstruction: every even number is a sum of two units modulo every `m` | `4·10^18` (CITED from memory) |
| Riemann hypothesis | generic | analytic (zero-free regions; zeros on the line to a large height) | — |
| Sloane ternary; Erdős #406 ternary | 1 (`x ↦ 2x` on `Z_3`) | **I1**: `2` is a primitive root mod `3^k`, so the `k` low digits of `2^n` equidistribute over units. The exceptions therefore have density zero. Statements about all digits need I3 (low-digit 3-adic structure against the Archimedean top digits) | Sloane: `5.9·10^21` (Saye 2022, as in §2) |
| weak Collatz (no nontrivial cycle) | 2 | **I2**: the cycle equation forces `|K log 2 − d log 3|` small. Steiner's 1-cycles and the `m`-cycle bounds use linear forms in logarithms | period `≤ 301,993` (THM-4591); convergent budget (THM-4633) |
| full Collatz | 2 | **I1** (Terras: parity vectors form a full shift; Tao's 3-adic mixing; the finite-section transfer operators mod `2^K` and mod `3^J` of THM-4821/4822) and **I2** (the cycle row above) | verified to at least `2^68` (CITED from memory) |
| twin primes | generic | analytic (bounded gaps) | — |
| odd perfect / odd multiply perfect exist | generic | structural: Euler's form, prime-factor counts and size bounds (CITED from memory) | odd perfect `> 10^1500` (CITED from memory) |
| Antihydra / Hydra / Bigfoot | 2, plus a counter | **I2**: halting iff `Σ(2v_i − 3) = −6`, a Parikh condition; the drift heuristic is I3 | simulation (bbchallenge) |
| Sierpiński / Riesel | per candidate, primality of `k·2^n ± 1` | **I1**: known members are certified by covering systems, which are cylinder conditions mod `lcm`. The open candidates are exactly those without a covering, so they need I3 | PrimeGrid ranges (CITED from memory) |
| Fermat's last theorem | — | resolved at I3 (modularity); its classical partial results (Sophie Germain, Kummer) were I1 | — |
| PA / ZF / ZFC consistency | generic | none inside the theory (Gödel) | — |
| LRC(14) | not a loop (a linear flow on `T^13`) | **I1** (finite box and height theorems: good-time sets are finite unions of intervals, speeds mod `q`) and **I2** (Bonferroni certificates, exact iff maximum multiplicity `≤ D`, LEM-020) | the repository's finite boxes |
| JC(2) / DC(2) | not a loop | structural (degree-by-degree Tarski-decidable closures) | degrees 18, 22 (§3) |
| Mahler 3/2 | 2 (`⌈3x/2⌉`) | **I1** (THM-4654's slack automaton) and **I2** (its block law: `≤ ⌊ℓ/3⌋+1` ones) | — |
| Hadwiger–Nelson `χ = 7` | not a loop | `Σ⁰₁` witnesses (finite unit-distance graphs). The residue colourings of field planes are I1 objects (THM-4558) | the repository's HN censuses |
| `3×3` magic square of squares | not a loop | I1 (no congruence obstruction) and I3 (the geometry of the parametrizing surfaces) | Sessions 57–61 searches |
| `H ≥ disc` | not a loop | structural (reduced to the strong base, THM-1950) | exhaustive `n ≤ 6`, sampled `7, 8` (THM-1950) |
| black-core connectivity | not a loop | structural. **Regular starts are PROVED at every order `≥ 7`**: exactly `m` black moves, with an unprescribed singleton (Session75) and with a prescribed singleton in either role (THM-4819 (D)) | every strong class through order 10 (THM-4811 (D2)) |
| Rule 30 aperiodicity | 8 local rules | **I1** (periodic points by transfer matrices; the dyadic atlases); THM-3778 is a finite-scale no-go | — |
| ABC | not a loop | I3 | — |
| Kakeya in `R^4`; AMM 12592 | not arithmetical / real-analytic | — | — |

**One new tournament row.** "Every regular tournament of odd order `n` is exactly `(n²−1)/8` native moves from the
transitive class" (HYP-9423 (7)).
- **Level.** It is `Π⁰₁`, since each order is a finite search.
- **Invariant.** Its strongest invariant is **I2**: THM-4819 (A)'s slab `L + R = (n²−1)/8` is a linear functional of the
  sorted score vector, the tournament's Parikh vector of out-degrees, and it confines every geodesic.
- **Horizon.** Order 11 complete. Order 13: all circulants and 1,000 random tournaments. Order 15: all circulants, three
  special classes and a random sample. Order 19: `P_19`.
- **What remains.** The open part is the I3 step of realizing each descending pair by a Hamilton path while avoiding
  sources and sinks.

**Reading.** The table sorts the rows two ways:
- **Two-branch rows.** The dynamical rows with `b = 2` (Collatz, weak Collatz, Mahler, Antihydra) all have partial
  results at I1 and I2 and none at I3.
- **Census rows.** The number-theoretic `Π⁰₁` rows with generic `b` have analytic partial results and huge horizons.

The ETH algorithm, a level-2 device over `F₂` (klein §3), meets neither group: the two-branch rows are bounded-treewidth
carry chains, and the census rows are not sparse `F₂` systems.

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
