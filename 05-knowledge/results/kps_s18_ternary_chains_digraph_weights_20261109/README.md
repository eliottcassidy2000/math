# Ternary Chan-Pak by port chains, and the weighted cycle formula for digraph path totals

**Session:** kind-pasteur-2026-11-09-S18. **Canon:** THM-5070, THM-5071, HYP-9640.

**Status.**
- PROVED: the transfer and the apex bridge (THM-5070 (1)-(3)); the weighted formula, parity law and one-defect
  exclusions (THM-5071 (1)-(4b)).
- FINITE-EXACT: coverage, the tower obstruction, the censuses and the core decompositions.
- REFUTED: the owner's matrix conjecture; Session234 refuted it first.
- OPEN: ternary Chan-Pak (HYP-9530 (A)) and binary Chan-Pak.

Proofs are in [PROOF.md](PROOF.md); the predictions are in [PREREGISTRATION.md](PREREGISTRATION.md) (P1-P5, all
CONFIRMED).

**Independent audit:** PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md)). 2 MAJOR and 15 MINOR findings, all
applied; MISTAKE-748. No wrong number or proof was found, and every count reproduced. The majors:
- (4b) "one reciprocal pair never gives 1, 4, 7": the cases 1 and 7 and the method are Session234's; 4 is new.
- The tower-direction check was glossed beyond its finite range (length <= 3; the audit adds length 4).

The audit extended the ternary 15-generator coverage to 3*10^6.

## Owner prompt (2026-11-09)

> keep going and prove the ternary Chan-Pak conjecture and beyond. spend the next session considering hamiltonian
> path count totals on digraphs, those graphs represented as square binary matrices, and consider the conjecture that
> all totals of 7 or 21 contain at least a single 1 on the main diagonal (one reflexive node) or one pair of zeros
> (coordinates mirrored across the main diagonal that are both 0), and study precisely how these families are
> generated and explain tournament structure. [...] consider increasing understanding of digraphs and their core
> subset tournaments at the same time [...] a tournament can be represented by a triangular tiling where every cell
> is binary, and a digraph can be represented by one where the hypotenuse is binary and the remaining triangle is
> ternary [...]

## Concurrent work: read first and credited

- **Session234** (864f45f40e, cdb77c02e9) was first on the digraph prompt.
  - It refuted the conjecture: semicomplete matrices with 7 and 21 paths exist at orders 4 and 5.
  - It classified the loopless semicomplete 7/21 families at every order, by strong cores.
  - It proved the sharp necessary conditions "7: a missing pair or exactly two reciprocal pairs" and
    "21: a missing pair or 1 <= b <= 3".
  - None of that is redone here. Its packet says it did not use the Grinberg-Stanley expansion.
- **opus-2026-11-09-S89** holds THM-5040-5042 and HYP-9610: even graphs, matrices as pairs of (n+1)-tilings, and
  ternary triangles with the binary/ternary dictionary.
- **mac-mini-2026-11-09-digraph721** holds THM-5055/5056 and HYP-9625.
- **The ternary inheritance:**
  - opus-S86 (THM-4940, HYP-9530);
  - opus-S87 (THM-4950: lane reduction, order-10 census, coverage to 29,035);
  - mac-mini-2026-11-08-ternary (THM-4952: ribbon automaton, base-3 ribbon census to 10^5);
  - the seam transfer (SEAM.md, THM-4875 (B), THM-4866 (B));
  - the one-apex recurrence (session25, Session191).

## 1. Ternary Chan-Pak (THM-5070, HYP-9640 (A)-(B))

**Port chains.** A port chain is built from one vertex u, the port, by repeatedly gluing a small block (B, w, u'):
- every old vertex beats every block vertex, except that w -> u;
- the port moves to u'.

This is SEAM.md's seam operation, up to converse.

**(1) The transfer, PROVED.** The state (I(Omega(T), x), E, O) evolves by a nonnegative integer 3 x 3 matrix that
depends only on the block, for every base x. Here E and O sum the packing counts of T - V(P) over paths P from the
port of even and odd order. The base-2 case is SEAM.md's transfer, read through OCF. The packing form for every x is
what is used here.

**(2) Apex bridges, PROVED.** If the odd cycles of each side all pass through its glued vertex, then

    I = I_1 I_2 + x c.

For one-apex sides this is explicit in the binary words W, W'.

**(3) The tower, PROVED; recurrence KNOWN from THM-4952 (C).** The one-vertex block is the tower, with
characteristic polynomial t^3 - t^2 - t - (x - 1). In base 3 this is (t - 2)(t^2 + t + 1).

**(4) Coverage, FINITE-EXACT; pre-registered P1-P3, CONFIRMED.**

| family | base | misses (all of them, up to the bound) | bound |
|---|---|---|---|
| blocks of order <= 5 (271 matrices) | 3 | **10, 31** | **10^6** |
| blocks of order <= 4 (57 matrices) | 3 | 10, 31, 61, 106, 217, 274, 337, 346, 607, 865 | 2*10^5 |
| 15 generators (a greedy subset of the 57) | 3 | the same ten | 10^6 (64,420,757 states) |
| blocks of order <= 5 (264 matrices) | 2 | **7, 21** | 2*10^5 |
| blocks of order <= 4 (56 matrices) | 2 | 7, 21, 47, 67, 127, 173, 491 | 2*10^5 |
| 13 generators | 2 | 7, 21, 47, 67, 127, 163, 173, 331, 491, 511, 799, 959 | 2*10^5 |

**Reading the table.**
- The first row is the new ternary record: every v = 1 (mod 3) up to 10^6 other than 10 and 31 is H_3 of a port
  chain. The previous range was 10^5 (THM-4952 (D)) and 29,035 (THM-4950 (3)).
- One explicit family realises both Chan-Pak spectra, and its only misses are the permanent gaps.
- The binary row is a control: THM-4875 (D) already reaches 2^32 with one-arc joins.

**(6) Beyond base 3: the gap sets of bases 4-8, FINITE-EXACT, not pre-registered.** For each base b = 4..8:
- m <= 80 is THM-4940 (6)'s exact table.
- Every m in (80, 30000] is settled by port chains with blocks of order <= 5, together with 61 hill-climbed
  tournaments of orders 8-11 for the values the chains miss (`hunt.c`).
- The witnesses were checked by a second script of this session (`verify_witnesses.py`), and independently by the
  audit (all 61 pass).
- Order-join closure was tried and settled nothing in this range.

So Gap_b below 30,000 is exactly the known small gap set:

| b | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|
| gaps (alpha-sums) | 3, 14, 27 | 3, 10, 27 | 3, 10, 14, 21, 46 | 3, 10, 14, 29 | 3, 35, 52 |

This is strong evidence for HYP-9530 (B) and its sharpening HYP-9640 (E).

**(5) Why there is no easy digit proof, FINITE-EXACT.** A digit (numeration) proof along the tower would need block
products that preserve the tower direction (3, 1, 2), the eigenvalue-2 eigenvector.
- Among the products of length <= 3 of all 57 order-<= 4 blocks, only the tower powers M_T^1, ..., M_T^12 do so.
- The audit adds length exactly 4: only M_T^4, ..., M_T^16.
- Longer products are not checked.

Separately, exploratory runs (not claimed) show that the 2-dimensional "sink-normalised" subsystem and the apex-bridge
family alone leave gaps that keep growing in the computed range.

**HYP-9640 (A), OPEN.** The 15 explicit base-3 generators reach every v = 1 (mod 3) with v >= 866. Checked to 10^6
here and to 3*10^6 by the audit.

**HYP-9640 (B), OPEN.** The 56 order-<= 4 base-2 matrices reach every odd v >= 493.

**Consequence.** With the finite censuses below those thresholds, (A) implies ternary Chan-Pak and (B) implies
binary Chan-Pak. Both conjectures become cofiniteness statements about the orbit of (1, 0, 1) under an explicit
finite set of nonnegative 3 x 3 integer matrices.

## 2. Digraph path totals (THM-5071)

### The formula

By Grinberg-Stanley's expansion of the Rédei-Berge function (arXiv:2307.05569), evaluated at (1, 0, 0, ...), every
loopless digraph D satisfies

    H(D) = sum over packings of vertex-disjoint admissible cycles of prod w(C).

The weights depend on the type of each unordered pair:

| pair type | matrix cells | cycles carrying it | weight |
|---|---|---|---|
| T-pair (tournament pair) | exactly one 1 | odd cycle of D using only T-pairs | 2 |
| | | even cycle using only T-pairs | 0 |
| B-pair (reciprocal) | 1 1 | any cycle of D through a B-pair, length >= 2 | 1 |
| Z-pair (missing) | 0 0 | any cycle of the complement through a Z-pair | (-1)^(length-1) |
| loop | diagonal 1 | none | invisible |

The formula was checked against direct counts on 232,263 digraphs (P5 gate, 0 failures).

### How this explains tournament structure

- **The factor 2 of OCF.** It is the pairing of each odd cycle with its reverse in the complement. A tournament is
  exactly a matrix in which every cycle has such a partner, so H = 1 + 2(...) is odd (Rédei).
- **What a reciprocal 11 pair does.** It removes the partner of every cycle through it, so those cycles weigh 1
  instead of 2.
- **What a missing 00 pair does.** It creates unpartnered cycles in the complement, with signed weight.
- **Loops.** The diagonal never enters (Session234 DG1). The "reflexive node" clause of the conjecture is therefore
  vacuous: a loop can be added to any matrix without changing its path total.

**Parity law (PROVED).** H(D) is congruent, mod 2, to the number of packings of defect cycles. This refines Berge's
theorem H(D) = H(complement) (mod 2); Rédei's theorem is the tournament case. The formula itself is Grinberg-Stanley's
Theorem 1.30 (arXiv:2307.05569v1), specialised at (1, 0, 0, ...).

### The 7/21 families, decomposed (FINITE-EXACT; Session234's cores)

| H = 7 core (order 4, b = 2) | decomposition |
|---|---|
| two classes | 7 = 1 + 6 x 1: six pairwise intersecting reciprocal cycles, no odd tournament cycle |
| one class | 7 = 1 + 5 x 1 + 1 x 1: one disjoint pair of reciprocal cycles |

- **Why 7 is out of reach for tournaments.** In tournaments 7 would need alpha = (3, 0): three pairwise
  intersecting odd cycles of weight 2, which is impossible (THM-4940 (3), THM-343).
- **How the cores escape.** Each core uses six weight-1 cycles instead.
- **The 21-cores.** The sixteen order-5 cores (b = 2 or 3) decompose as 1 + 18 + 2, 1 + 17 + 3, 1 + 12 + 6 + 2,
  and so on (`out/core_weights.out`). Each mixes weight-1 reciprocal cycles with at most three odd tournament cycles.

### One defect at a time (PROVED / FINITE-EXACT; pre-registered P4, CONFIRMED)

**Exactly one reciprocal pair, no missing pair.** The total is never 1, 4 or 7, at any order.
- **Credit.** The cases 1 and 7, and the method, are Session234's: component factorisation DG3, strong spanning
  tournament DG4, and DG5's "7 needs exactly two reciprocal pairs". New here: the exclusion of 4.
- **Proof.** Factor over strong components. A strong semicomplete digraph with one reciprocal pair has a strong
  tournament resolution (DG4; a self-contained four-cell check is in PROOF.md). The Busch floor f(5) = 9 then leaves
  orders <= 4, and those are checked exhaustively.
- **What is attained.** Every other value 2..197 is attained at orders <= 7. In particular 21 is attained with one
  reciprocal pair, already at order 6 (first found by Session234).

**Exactly one missing pair, no reciprocal pair.** Every value 0..129 is attained at orders <= 7. The H = 7 case at
order 5 was already in Session234's census (600 labelled digraphs). The examples are
- 7 at order 5: rows 00001 10000 01000 11100 01110, with pair {0, 2} missing;
- 21 at order 6.

**Every positive integer is a semicomplete total.** Use TT_n plus back arcs from the sink to the positions of the
binary digits of h - 1; then H = h.

**What survives of the conjecture.** The diagonal clause is vacuous, and an 11 pair can do what the 00 pair does.
The sharp form is:
- 7 and 21 force a defect pair;
- a single 00 pair suffices for both;
- a single 11 pair suffices for 21 but never for 7 (both Session234's);
- for 7 one needs a 00 pair or exactly two 11 pairs (Session234 DG5).

1, 4 and 7 are excluded by one 11 pair at all orders (proved; 4 is the new case). That nothing else is excluded is
HYP-9640 (C), verified below 198.

## 3. Binary and ternary (typed ANALOGY; S89 holds the dictionary lane)

- **Binary triangle (tournament).** Cycle weights are {0, 2}; H = I(Omega, 2) is odd; the gaps are 7 and 21.
- **Ternary triangle (semicomplete digraph; cell states ->, <-, both).** Cycle weights are {0, 1, 2}, and every
  positive integer occurs. The binary hypotenuse is invisible. The binary digits of h - 1 can be written onto the
  sink's back arcs.
- **The other "ternary".** The base-3 packing count H_3 = I(Omega, 3) of a tournament is a different object. It is
  the x = 3 specialisation of the same port-chain transfer, and it has the same gap alpha-sums {3, 10} as base 2.
- **"2 and 3" in the tower.** In base 3 the tower grows as exactly 2^n with a period-3 correction (eigenvalues 2, ω,
  ω^2). The exact count is H_3 = (3 * 2^n + c_n)/7 (THM-4952 (C)).

Whether some digraph construction D(T) has H(D(T)) = H_3(T) is left open. The weighted formula suggests it would
need weight 3 on every odd cycle. This is a heuristic, since nothing forces a construction to act cycle by cycle. No
such construction is claimed.

## 4. Files

| file | purpose |
|---|---|
| `chain_blocks.py` | block matrices N_x(B, w, u') for all blocks of a given order; `.reps` lists a representative block for each |
| `chain_cover.c` | saturated breadth-first coverage of the port-chain semigroup |
| `chain_check.py`, `packing_count.c` | independent brute-force check of the transfer on random chains |
| `tower_direction.py` | the tower-direction obstruction (5) |
| `digraph_weights.c` | check of the weighted formula; one-defect censuses |
| `core_weights.py` | weighted decompositions of the 7/21 cores |
| `one_defect_witness.py` | explicit one-defect witnesses; check of the sink construction |
| `bases_gaps.py` | chain coverage in bases 4-8 and the open alpha-sums after order joins |
| `hunt.c`, `verify_witnesses.py` | hill-climb for the open alpha-sums; independent from-scratch check of the 61 witnesses |
| `out/` | every output quoted above, plus the 15 and 13 generator lists with block representatives |

**Reproduce:**
- Build the block matrices: `python3 chain_blocks.py 3 4 m3k4.txt`.
- Compile: `gcc -O2 -o chain_cover chain_cover.c`.
- Run the coverage: `./chain_cover m3k4.txt 3 200000`.
- The other checks run the same way.

The run at 10^6 with the 15 generators needs about 1.5 GB of memory.
