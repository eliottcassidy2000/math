# Why the multicolour Ramsey targets R(3,4,4) ≤ 73, R(3,4,5) ≤ 137, R(3,3,6) ≤ 86, R(4,4,4) ≤ 198 are out of reach of every recursion-type argument: the saturated-port cascade, its exhaustion on the 2026 table, the (3,4,4;76) ledger, and the (3,3,4;29) ports

**Session:** klein-2026-10-30 (worktree `codex/session-ramsey-20261030`). Canon:
[THM-4859](../../../01-canon/theorems/THM-4859-saturated-port-cascade-for-multicolour-ramsey-upper-bounds-exhausted-on-the-2026-table-at-the-unsaturated-334-port-the-344-76-ledger-and-the-unique-circulant-334-29-colouring.md);
open questions: [HYP-9450](../../hypotheses/HYP-9450-rigidity-of-334-29-ports-exclusion-of-344-76-and-335-56-and-what-the-ramsey-targets-73-137-86-198-would-require.md).

**Status.** No new Ramsey bound is proved or claimed. PROVED (elementary, generalising a CITED argument): the depth-`t`
saturated-port cascade of §2. FINITE-EXACT: its mechanical evaluation on the April 2026 table (§3), the `(3,4,4;76)`
ledger identities (§4), the circulant `(3,3,4;29)` census (§5). CITED: Radziszowski, *Small Ramsey Numbers*, DS1
revision 18 (24 April 2026), sections 2.1 and 6.1; L. Boza, *New upper bounds for the classical Ramsey numbers
R(4,4,4), R(3,4,5) and R(3,3,6)*, arXiv:2603.10851 (March 2026); Codish–Frank–Itzhakov–Miller 2016 (`R(3,3,4) = 30`);
McKay's Ramsey database (uniqueness of the `(4,4;17)` graph Paley17 and of the `(3,5;13)` graph). Same-agent review;
not independently audited.

**Sources read in the repository.** [ramsey_saturated_ports_20261030](../ramsey_saturated_ports_20261030/PROOFS.md)
(saturation and exact edge counts inside the 29-vertex blue ports of a hypothetical `(3,3,5;56)` colouring; P7 on
`R(3,3,6)`), [terminal_fees_phase_covers_20261030](../terminal_fees_phase_covers_20261030/COLLATZ.md) (the nine closed
E11 price families and the seven open E12 families; LRC.md on arbitrary attached speeds),
[return_kernel_core_transplants_20261030](../return_kernel_core_transplants_20261030/LRC.md),
[entry_cosets_resonant_echoes_20261030](../entry_cosets_resonant_echoes_20261030/README.md). The repository's other
"Ramsey" packets of this week ([ramsey_circuit_interfaces_20261030](../ramsey_circuit_interfaces_20261030/README.md),
[session117_ramsey_circuit_interfaces_20261030](../session117_ramsey_circuit_interfaces_20261030/SYNTHESIS.md)) concern
transitive-subtournament capacities and affine Ramsey geometry over `F_2`, not graph Ramsey numbers; they are not inputs.

**Owner's prompt (verbatim).** "the seven twelfth-layer price families and arbitrary-input lonely runner coverage. also
attempt to decrease th eupper bound for ramsey numbers R(3,4,4) to 73 R(3,4,5) to 137 R(3,6,6) to 86 and R(4,4,4) to
198 etc pull in recent repo work and reroute creatively".

Vocabulary map. "Seven twelfth-layer price families" = the seven native E12 price families left open in
terminal_fees COLLATZ.md §5 (the nine E11 families are closed there). "Arbitrary-input lonely runner coverage" = the
open universal case of the attached-speed coverage algorithm (terminal_fees LRC.md: every thirteenth speed certified
for 28 declared twelve-speed cores at radius `1/14`; arbitrary cores open). "R(3,6,6) to 86" is read as **R(3,3,6)**:
the survey gives `338 ≤ R(3,6,6)`, so 86 cannot be an upper bound there, while `86 = 2·18 + 51 − 1` is exactly the
recursive bound for `R(3,3,6)` that would follow from `R(3,3,5) ≤ 51`.

## 1. The state of the art (DS1.18, April 2026, with Boza 2026)

| number | lower | upper | where the upper bound comes from |
|---|---:|---:|---|
| `R(3,3,4)` | 30 | 30 | Codish–Frank–Itzhakov–Miller 2016 (exact; the recursion gives 34) |
| `R(3,3,5)` | 45 | 57 | recursion `14 + 14 + 30 − 1` |
| `R(3,4,4)` | 55 | 77 | recursion `18 + 30 + 30 − 1` |
| `R(3,3,6)` | 61 | 91 | recursion `18 + 18 + 57 − 1 = 92`, strict (Boza 2026; also parity, see §3) |
| `R(3,4,5)` | 89 | 157 | recursion `25 + 57 + 77 − 1 = 158`, strict (Boza 2026) |
| `R(4,4,4)` | 128 | 229 | recursion `3·77 − 1 = 230`, strict (Boza 2026) |
| `R(3,3,3,3)` | 51 | 62 | Fettes–Kramer–Radziszowski (recursion gives 66) |
| `R(3,3,3,4)` | 97 | 149 | recursion `3·30 + 62 − 2 = 150`, parity strict |

Exact two-colour inputs: `R(3,3) = 6, R(3,4) = 9, R(3,5) = 14, R(3,6) = 18, R(4,4) = 18, R(4,5) = 25`; the unique
`(4,4;17)` graph is Paley17 and the unique `(3,5;13)` graph is `Cay(Z/13, {1,5,8,12})`.

## 2. The saturated-port cascade (PROVED; depth 1 is Boza's Theorem 2.1)

Write `P(k) = 2 − r + Σ_i U(k^{(i)})`, where `k^{(i)}` lowers `k_i` by one and `U` is any valid upper bound (an exact
value where known). Then `R(k) ≤ P(k)` (Greenwood–Gleason). Suppose `R(k) = P(k)`. Then every term attains its bound,
`R(k^{(i)}) = U(k^{(i)})`, and in a colouring on `n_0 = P − 1` vertices every vertex has colour-`i` degree exactly
`n_1^{(i)} = U(k^{(i)}) − 1` (saturation). Fix a colour `i` with `k_i ≥ 4` and write `k^{(i,j)}` for `k` with `k_i`
lowered by `j`. If moreover `U(k^{(i,j)}) = P(k^{(i,j)})` for `j = 1, …, t−2` (the nested colour-`i` ports are
themselves saturated), then every colour-`i` port of level `j` has forced size `n_j = U(k^{(i,j)}) − 1` and forced
colour-`i` degree `n_{j+1}`, so

> **the number of colour-`i` cliques `K_t` is exactly `n_0 n_1 ⋯ n_{t−1} / t!`, which must be a non-negative integer.**

If it is not, `R(k) ≤ P(k) − 1`. The case `t = 3` is Boza's Theorem 2.1 (`3 ∤ n_0 · n_1 n_2 / 2`); the parity
strictness of the classical statement is the case `t = 2` (`2 ∤ n_0 n_1`). Each forced-count condition is one unit of
improvement and no more: once `R(k) ≤ P(k) − 1` the colouring on `P − 2` vertices is no longer saturated, and the
method has nothing further to say about it.

## 3. The scan (FINITE-EXACT; `scripts/ramsey_cascade.py`, `data/ramsey_cascade.out`)

With the table of §1 as `U` and recursion-plus-parity as the fallback for the rest, every three-colour tuple with
entries `3..7` and every four-colour tuple with entries `3..5` was tested. Results:

- Boza's three bounds are reproduced exactly: `R(3,3,6)`: `91 · 56 · 29 / 6 = 73892/3`; `R(3,4,5)`: `157 · 56 · 13 / 6 =
  57148/3` (colour 4, port `(3,3,5;56)`) and `157 · 76 · 29 / 6 = 173014/3` (colour 5, port `(3,4,4;76)`); `R(4,4,4)`:
  `229 · 76 · 17 / 6 = 147934/3`. The parity rule alone already gives `R(3,3,6) ≤ 91` (`P = 92` even, term 18 even).
- **The cascade is dead for all four targets at the same place.** For `R(3,4,4) = 77` the colour-4 ports are `(3,3,4;29)`
  colourings with `U = 30 < P = 34`: unsaturated, degrees not forced, no forced count. For `R(3,3,5) = 57` the colour-5
  port is the same `(3,3,4;29)`. `R(3,3,6)` and `R(3,4,5)` and `R(4,4,4)` now sit one below their recursive bounds, so
  their extremal colourings are unsaturated and the method cannot touch them again until `R(3,3,5)` or `R(3,4,4)`
  moves. Every other flagged case in the scan (`R(3,3,4,5) ≤ 1079`, `R(3,6,6) ≤ 1623`, `R(4,5,5) ≤ 1441`, …) lives on
  tuples whose bounds are themselves recursive; these are implicit one-unit improvements of the kind Boza's Theorem
  2.5 records, not new information about small numbers.
- **What each target would require** (recursion-level arithmetic, exact):

| target | current | what would deliver it by recursion | gap |
|---|---|---|---|
| `R(3,4,4) ≤ 73` | 77 | nothing: `77 = 18 + 30 + 30 − 1` with all three terms exact | direct exclusion of `(3,4,4;n)`, `n = 76, 75, 74, 73` |
| `R(3,3,6) ≤ 86` | 91 | `R(3,3,5) ≤ 51` (`2·18 + 51 − 1 = 86`; `≤ 52` gives 87) | six units on `R(3,3,5)`, whose lower bound is 45 |
| `R(3,4,5) ≤ 137` | 157 | `R(3,3,5) + R(3,4,4) ≤ 113` (now 134; e.g. `45 + 68` or `51 + 62`) | twenty-one units across the two |
| `R(4,4,4) ≤ 198` | 229 | `R(3,4,4) ≤ 66` (`3·66 − 1 = 197`; `≤ 67` gives 200) | eleven units on `R(3,4,4)`, whose lower bound is 55 |

The intermediate steps are also exact: `R(3,4,4) ≤ 76` gives `R(4,4,4) ≤ 227` and `R(3,4,5) ≤ 157` (no change);
`R(3,3,5) ≤ 56` gives `R(3,3,6) ≤ 91` (no change) and `R(3,4,5) ≤ 157` (no change); `R(3,3,5) ≤ 55` gives `R(3,3,6) ≤ 89`.

## 4. The `(3,4,4;76)` saturation ledger (PROVED, elementary)

If `R(3,4,4) = 77`, a `(3,4,4;76)` colouring exists, colour 1 triangle-free, colours 2 and 3 `K_4`-free, and every
vertex has degrees `(17, 29, 29)`. Then: the colour-1 graph is 17-regular triangle-free (646 edges); colours 2 and 3 are
29-regular `K_4`-free (1102 edges each); every `N_1(v)` carries the unique `(4,4;17)` colouring, Paley17, with 68
colour-2 and 68 colour-3 edges, 8-regular in each; for every `r ∈ N_1(v)`, `|N_1(r) ∩ N_2(v)| = |N_1(r) ∩ N_3(v)| = 8`
(v has eight colour-2 and eight colour-3 neighbours in the Paley port `N_1(r)`), hence `e_1(N_1(v), N_2(v)) =
e_1(N_1(v), N_3(v)) = 136`; summing colour-1 degrees over each 29-port, `e_1(N_2(v)) = e_1(N_3(v)) =: e^*(v)` with
`2e^*(v) + e_1(N_2(v), N_3(v)) = 357`, so the colour-1 cross count between the two 29-ports is **odd**, and since colour
1 inside a 29-port has degree at most 8 (a colour-1 neighbourhood inside the port is a `(3,4)`-set), `e^* ≤ 116` and
`e_1(N_2, N_3) ≥ 125`; the mean colour-1 codegree over colour-2 edges is `136/29`, non-integral, so it is not
constant. Each 29-port is a `(3,3,4;29)` colouring with in-port degree caps `(8, 8, 16)` summing to at most 32 against
the actual 28: a slack of four per vertex. A Boza-type contradiction would need `e_2(N_2(v))` forced to a constant `c`
with `3 ∤ 76c`; nothing forces it, exactly because the port is unsaturated. (The same ledger for `(3,3,5;56)` is the
saturated-ports packet's P1–P2, where the 13-ports are unique and the 29-port counts `(104, 104, 198)` are forced, and
still no contradiction appears.)

## 5. The `(3,3,4;29)` ports (FINITE-EXACT census; SAT probes)

Every route in §3 and §4 bottoms out in the same objects: the `(3,3,4;29)` colourings, the ports of both
`R(3,4,4) = 77` and `R(3,3,5) = 57`. Two finite facts about them:

- **Circulant census.** Colouring the 14 differences of `Z/29` with two triangle-free colours and one `K_4`-free colour
  (`scripts/circ_probe_334.py`, exhaustive over `3^14`) gives exactly 14 labelled cyclic `(3,3,4;29)` colourings, all
  with difference-class split `(4, 4, 6)`, i.e. edge counts `(116, 116, 174)` and degrees `(8, 8, 12)`; up to the
  multiplier group and the swap of the two triangle-free colours there is **one** cyclic colouring (classes
  `D_0 = {1,3,7,12}`, `D_1 = {2,5,6,14}`, `D_2 = {4,8,9,10,11,13}`), presumably Kalbfleisch's 1966 construction (DS1.18
  credits the lower bound to [Ka2]). In it both triangle-free colours are at their degree cap 8 and the `K_4`-free
  colour has degree 12.
- **SAT probes (inconclusive).** `scripts/port_probe_334.py` encodes a `(3,3,4;29)` colouring as a 1218-variable CNF (32,685
  clauses: two triangle-free colours, one `K_4`-free colour, two symmetry-breaking units) and was run in two modes, plain
  with Glucose4 and with the cardinality constraint "vertex 0 has colour-0 degree ≤ 7" with CaDiCaL (the question whose
  positive answer shows the ports are not degree-rigid). Neither run returned a model within about twenty minutes
  before the session was interrupted (`data/probe_deg7.log`, `data/probe_glucose.log`). The existence of `(3,3,4;29)`
  colourings is not in doubt (the cyclic one above satisfies the plain instance); the SAT difficulty is the solver's,
  not the problem's. Rigidity is therefore OPEN (HYP-9450 (1)); the instance is in the packet.

## 6. Verdict on the four targets (typed)

OPEN, and not approachable by any recursion, parity or forced-count argument: all such arguments are one-unit tools
and all four cascades terminate at the unsaturated `(3,3,4;29)` port. The only levers are (i) direct exclusion of
`(3,4,4;76)` and `(3,3,5;56)` colourings, which with the ledgers of §4 and the saturated-ports packet is a structured
but very large SAT problem (2850 and 1540 pair variables with `K_4`/`K_5` clauses), and (ii) a census of `(3,3,4;29)`
colourings, which would turn the ledgers' free quantities (`e^*(v)`, `e_2(N_2(v))`) into finitely many cases. Each
unit of `R(3,3,5)` or `R(3,4,4)` propagates as in §3; the targets ask for six, eleven and twenty-one such units.

## 7. The other two anchors (typed reading, no closure)

- **Seven E12 price families** (terminal_fees COLLATZ.md §5). Their reduction (C3) asks, for `(b, a) ∈ {(1,0), (1,1),
  (1,2), (1,3), (2,0), (2,1), (3,0)}`, that `Σ_{i≤b} g_11(X_i) ≥ 1 − a/5` whenever `Σ X_i = −a(4,4)`, where
  `g_11 > γ = 5878847/19531250 ≈ 0.301` is the proved uniform floor. Exact thresholds: a uniform floor `g_11 ≥ 1/3`
  closes `(3,0)`, `≥ 2/5` closes `(2,1)`, `≥ 1/2` closes `(2,0)`; the four `b = 1` families need tangent-specific prices
  `g_11(−a(4,4)) ≥ 1 − a/5`, i.e. `1, 4/5, 3/5, 2/5` at `a = 0,1,2,3`, which no uniform floor can supply. So the
  closure is a terminal-conditioned refinement at the four axis tangents plus a uniform lift from 0.301 to 0.5, as the
  packet itself says; nothing here changes that.
- **Arbitrary-input lonely runner coverage.** The attached-speed algorithm certifies every thirteenth speed for a
  declared core; "arbitrary cores" is LRC(14) itself. In the information vocabulary of THM-4840's packet the
  complete interval cover is the annihilator of all future attached speeds; its compression (LRC.md next question 1)
  is a question of which dual words have bounded support, the same shape as the cryptid criterion. No progress.

## 8. Reproduction

`scripts/ramsey_cascade.py` (standard library; seconds), `scripts/circ_probe_334.py` (seconds), both run under
`python3` and `python3 -O` with byte-identical outputs in `receipts/`. `scripts/port_probe_334.py` needs `pysat`
(CaDiCaL/Glucose); its outputs are solver-dependent samples saved in `data/`, not receipts.
