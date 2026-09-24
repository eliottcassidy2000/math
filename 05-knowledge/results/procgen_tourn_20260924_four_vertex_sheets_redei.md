# The four-vertex reading of `3n±1`: the two diamonds are time reversal, not the sheet swap, and Rédei's parity has no Collatz counterpart

**Status: the owner's claim is REFUTED as stated and REAL in a corrected form. The two converse 4-tournaments `(1,1,1,3)` and `(0,2,2,2)` do arise from the Collatz pair structure, but they are the forward map and the inverse tree of one AM-fair pair. Reversing all 6 arcs is time reversal (up to the pair-sum reflection), and exactly negation ∘ time reversal. The sheet swap `x -> -x` fixes the class. "3 = H(C3)" is an ANALOGY. Among the candidates tested, no Collatz count has Rédei's or Sperner's odd-count property, and "fixed point chain growth" is an ANALOGY with precise differences.**

* **PROVED** (hand proofs below; the scripts re-check every finite ingredient):
  1. **Theorem A (the pair quadruple).** Take the AM-fair pair `{s_o, s_e = s_o + b}` of the sheet `3n+b`. The tournament "map arcs + numerical order" on `{s_o, s_e, T(s_o), T(s_e)}` has class:
     * `(0,2,2,2)` (3-cycle over a sink) iff `b s_o > 0`;
     * transitive iff `b s_o < 0`.

     The inverse-tree tournament is its converse, via `rho(x) = s_o + s_e - x`. Negation fixes the class. Exactly, `(F+)^op = F-∘nu`. `H = 2 + sgn(b s_o)`.
  2. **Proposition B.** No local construction (map arcs + any mix of order and fixed-label orientations) gives the two sheets the converse pair on a common side of 0. Their two tournaments differ in at most one arc.
  3. **The involution dictionary.** For every construction built from dynamical arcs and order arcs, converse = (negation) ∘ (time reversal), exactly. On the AM-fair crossing pattern, negation ≅ id and time reversal ≅ converse. On Mahler's pattern, negation ≅ time reversal ≅ (TT ↔ STRONG) and converse ≅ id.
  4. **Proposition D and Theorem C (windows).**
     * On a fixed vertex set of integers `>= 3`, the map+order tournaments of `3n+1` and `3n-1` coincide (Proposition D).
     * For generic orbit windows, reversing the parity word gives the converse tournament (Theorem C).
  5. **Rédei as one fixed point.** The OCF expansion `H = sum_sigma 2^psi(sigma)` carries an involution with exactly one fixed point, the identity, so `H` is odd. By contrast, the 2-adic periodic points of `T` number `2^p` and are paired by a **free** involution `iota` that commutes with `T`.
  6. **Type graph.** The forward type graph of `3n+1`, `3n-1` and `5n+1` mod `2^k` is the de Bruijn graph `B(2,k)`. Its Hamiltonian cycles number `2^(2^(k-1)-k)` (classical count, re-verified here).
* **FINITE-EXACT:**
  * the 512-rule census: 0 rules separate the sheets on a common side; 24 gauge rules give the converse pair on ν-related quadruples, and only 8 (Type II) make negation the converse on every pattern;
  * window crossing census, both sheets, `m <= 64`;
  * identical statistic distributions near `10^6`;
  * `dH` takes both signs at crossings;
  * `H` separates the two sides of 0 for every generic word with `m <= 11` except `101` and `01110`;
  * preimage-count moments for `d <= 15`, and `R_p(n)` parities for `p <= 16`.
* **CITED:** Rédei 1934; Sperner 1928; Knaster–Kuratowski–Mazurkiewicz 1929; Cohen 1967; Kuhn 1968; Papadimitriou 1994; Flye Sainte-Marie 1894 / de Bruijn 1946 / van Aardenne-Ehrenfest–de Bruijn 1951; Lagarias 1985. §9 lists exactly what was read. Repository: THM-001, 002, 070, 466, 510, 511, 582, 589, 647, 2228, 4057, 4470, and the inverse-tree and order-laws notes.
* **REFUTED:**
  * "3n+1 and 3n-1 are the two diamonds, swapped by `x -> -x`", for every intrinsic construction;
  * "3 = H(C3)" as an identity: every odd multiplier `q >= 3` gives `H = 3`.
* **ANALOGY:**
  * Rédei / Sperner parity versus Collatz fixed points;
  * `H`-growth versus Collatz counts, where the second moments differ.
* **OPEN:** the limit of the Haar ratio `E[N_d^2]/E[N_d]^2`, about 1.09 on units through `d = 15`.

Session `collatz-procgen-20260922`, lane "four-vertex tournaments / Rédei", 2026-09-24. No HYP or THM file was created. No independent-agent audit yet.
* **Scripts:** `04-computation/experiments/procgen_tourn_20260924_{four_vertex,windows,redei}.py`, runner `procgen_tourn_20260924_run.py`.
* **Output:** [`procgen_tourn_20260924.out`](procgen_tourn_20260924.out), sections A1–A4, B0–B5, C0–C4. Every check raises on failure.

## 0. The claim, the probe, the answer

The owner's words (verbatim): "the '3' and '±1' is an underlying tournament on four vertices, specifically the 2 isomorphism classes which swap between each other when all 6 arcs are inverted", and "connect Rédei, Hamiltonian path, fixed point chain growth".

**The coordinator's probe is correct in every detail (A2, A3).** Let `i >= 2`, and use map arcs plus the order "smaller → larger" on the other pairs:
* On `{i, 2i-1, 2i, 3i-1}` (`3n+1`) this gives `(0,2,2,2)`: the 3-cycle `i -> 2i-1 -> 2i -> i` over the sink `3i-1`.
* On `{i, 2i, 2i+1, 3i+1}` (`3n-1`) it gives the transitive tournament. So does `3n+1` on negative pairs.
* The order "larger → smaller" gives `(0,2,2,2)` for `3n+1` again.

| owner's statement | verdict | where |
|---|---|---|
| the diamonds `(1,1,1,3)`, `(0,2,2,2)` arise from the Collatz pair structure | REAL (Theorem A) | §2 |
| they swap when all arcs are reversed | REAL: converse = time reversal (up to `rho`), and exactly = negation ∘ time reversal | §2, §3 |
| `3n+1` is one diamond, `3n-1` the other, swapped by `x -> -x` | REFUTED for every intrinsic construction. It is realisable only with a gauge arc, and then the gauge decides which sheet gets which diamond | §3 |
| "3" is the 3-cycle, `H(C3) = 3` | ANALOGY: every odd `q >= 3` gives the same class | §4.7, §5.1 |
| "±1" is the source or sink | REAL once read as the **time direction**: forward = sink, inverse tree = source. The offset really does flip under inversion, since `3n+b` inverts to `(y-b)/3` | §2.2 |
| Rédei ↔ Collatz fixed points / chain growth | ANALOGY: no odd-count law; the parity mechanisms are opposite | §6, §7 |

## 1. The four tournaments on four vertices (A1)

| class | score | labelled | `H` | `c3` | converse |
|---|---|---|---|---|---|
| TT (transitive) | `(0,1,2,3)` | 24 | 1 | 0 | TT |
| SRC+C3 (source over a 3-cycle) | `(1,1,1,3)` | 8 | 3 | 1 | C3+SNK |
| C3+SNK (3-cycle over a sink) | `(0,2,2,2)` | 8 | 3 | 1 | SRC+C3 |
| STRONG | `(1,1,2,2)` | 24 | 5 | 2 | STRONG |

* `H = 1 + 2 c3` on all 64 tournaments (the OCF at `n = 4`; THM-466(ii)).
* **Distance lemma.** One arc reversal changes exactly two scores, by `+1` and `-1`. From `(0,2,2,2)` it reaches only `(1,1,2,2)` or `(0,1,2,3)`. So the labelled arc distance between the two diamonds is at least 2, and it equals 2 (A1).

## 2. The AM-fair pair quadruple

### 2.1 Theorem A

Notation:
* `T_b(x) = x/2` for `x` even and `(3x+b)/2` for `x` odd, with `b = ±1`.
* The AM-fair pair of sheet `b` is `{s_o, s_e}` with `s_o` odd and `s_e = s_o + b`: `{2i-1, 2i}` for `b = +1` and `{2i, 2i+1}` for `b = -1` (THM-4470 §1).
* Images: `t_o = T_b(s_o) = (3s_o+b)/2` and `t_e = s_e/2 = (s_o+b)/2`.

Three identities:

* **(I1)** `t_o + t_e = s_o + s_e` (AM-fairness).
* **(I2)** `t_o - t_e = s_o`.
* **(I3)** `(s_o - s_e)(t_o - t_e) = -b s_o`.

The four numbers are distinct iff `|s_o| >= 3`. At `|s_o| = 1` the quadruple collapses onto the cycles of length at most 2:
* for `+`: `{1,2}`, and `{-1,0}` (both fixed);
* for `-`: `{0,1}` (both fixed), and `{-1,-2}`.

On the labels `s_o, s_e, t_o, t_e` put:
* the map arcs `M = {s_o -> t_o, s_e -> t_e}`;
* the inverse-tree arcs `M^op = {t_o -E-> s_o, t_e -D-> s_e}`, with `E(x) = (2x-b)/3` and `D(x) = 2x`;
* on the other four pairs, the order `O` (smaller → larger) or `O^op`.

Then `F+ = M ∪ O` is the forward tournament and `F- = M^op ∪ O` the inverse-tree tournament.

**Theorem A (PROVED; A2 checks all 4000 non-degenerate quadruples with `|s_o| <= 2001`).** Let `|s_o| >= 3` and `u = b s_o`.
* **(a)** `F+` is C3+SNK if `u > 0` and transitive if `u < 0`. The same holds for `M ∪ O^op`.
* **(b)** `rho(x) = s_o + s_e - x` maps the quadruple onto itself. It swaps `s_o <-> s_e` and `t_o <-> t_e`, reverses the order and preserves the set `M`. Hence `rho: M∪O ≅ M∪O^op` and `rho: (F+)^op ≅ F-`. So `F-` is SRC+C3 if `u > 0` and transitive if `u < 0`: **time reversal = converse**.
* **(c)** `nu(x) = -x` carries `Q_b(s_o)` onto `Q_{-b}(-s_o)` label by label and reverses the order. So `F+(nu Q) = M ∪ O^op ≅ F+(Q)`: **negation fixes the class.**
* **(d)** `(F+(Q))^op = F-(nu Q)` as labelled tournaments: **converse = negation ∘ time reversal, exactly.**
* **(e)** `H(F+) = H(F-) = 2 + sgn(u)` and `c3 = [u > 0]`. On positive pairs, `H = 2 + b`.

*Proof.*
1. Multiply every value by `b` (for `b = -1` this is `nu`). The quadruple becomes `(u+1)/2, u, u+1, (3u+1)/2` for `t_e, s_o, s_e, t_o`.
2. For `u >= 3` the increasing order is `t_e < s_o < s_e < t_o`. Each member of the pair jumps over its partner: the pair is reversed by `T`, which is (I3).
   * In positions 1–4, the map arcs are `2 -> 4`, which agrees with the order, and `3 -> 1`, which goes against it with range 2.
   * This is the transitive tournament with the arc between positions 1 and 3 reversed: the 3-cycle `(t_e, s_o, s_e)` over the sink `t_o`.
3. For `u <= -3` the order is `t_o < s_o < s_e < t_e`: each member moves away from its partner.
   * The map arcs are `2 -> 1` and `3 -> 4`. The only arc against the order joins adjacent positions.
   * Reversing an adjacent arc of a transitive tournament leaves it transitive.
4. The true order is the scaled one (`b = +1`) or its reverse (`b = -1`). By (b) the class does not see this reversal.
5. For (b): `rho(s_o) = s_e` by definition and `rho(t_o) = t_e` by (I1), so `rho` maps `s_o -> t_o` onto `s_e -> t_e`.
6. For (c): `nu T_b = T_{-b} nu` (transport theorem, [inverse-tree note](collatz_procgen_20260924_inverse_tree_mod192.md) Theorem 6), and `-s_o + (-b) = -s_e` is the AM-fair partner of `-s_o` on sheet `-b`.
7. For (d): both sides equal `M^op ∪ O^op` on the labels.
8. For (e): at `n = 4`, `H = 1 + 2c3` (THM-466(ii)). ∎

### 2.2 What Theorem A says about "3" and "±1"

* **The class is the side bit `sgn(b·x)`.**
  * Diamonds sit exactly on the side of the contracting 2-cycle: `{1,2}` for `3n+1`, `{-1,-2}` for `3n-1`. The transitive tournament sits on the other side.
  * This is the inverse-tree note's `T`-definable half-line (§2.2 there), in tournament costume. It is the sign law, and it is consistent with Corollary 7 there: the statistic uses the order, so it may see the side.
  * On the positive integers it separates `3n+1` (diamonds) from `3n-1` (transitive), because the positive half-line is the contracting side of one sheet and the expanding side of the other.
* **Source versus sink is the time direction, not the sheet.** Under either order convention the forward map gives C3+SNK and the inverse tree gives SRC+C3.
* **The "±1" does flip under inversion.**
  * The inverse of `x -> 3x+b` is `y -> (y-b)/3`, and the odd inverse branch of `T_b` is `E_b(x) = (2x-b)/3`.
  * `nu` (sheet swap) and `tau` (inversion) both flip the sign of the offset. The tournament converse follows `tau`, and exactly `nu∘tau`, never `nu` alone.
* **In THM-4470's pairing family** (bit `eps_i` per pair; every member AM-fair), Theorem A applies pair by pair. On positives, the class of pair `i` is C3+SNK iff `eps_i = 0`.
  * Collatz is all diamonds; `3n-1` (the all-flipped pairing, shifted) is all transitive.
  * The class sequence *is* the pairing word. So THM-4470 §4 applies verbatim: flipping a density-zero set of pairs changes the class on a density-zero set and creates a divergent orbit. The class is DEFECT-blind.

## 3. Every local construction (A3)

A **local rule** chooses `M` or `M^op`, and orients each non-map pair by one of four options:
* the pairs are `A = {s_o,s_e}`, `B = {t_o,t_e}`, `C = {s_o,t_e}` and `D = {s_e,t_o}`;
* the options are ascending order, descending order, or one of the two fixed label directions (a gauge).

There are `2·4^4 = 512` rules.

**Lemma (comparison vector).** For non-degenerate pairs:
* the comparison on `A` depends only on the sheet, since `s_e - s_o = b`;
* those on `B`, `C`, `D` depend only on the side, since `t_o - t_e = s_o` and `s_o - t_e = t_o - s_e = (s_o - b)/2` all have the sign of `s_o`. ∎

**Proposition B (PROVED).** Under any local rule, the tournaments of the two sheets on a common side of 0 differ at most in the arc on `A`. By the distance lemma, **no local rule assigns the converse pair to `3n+1` and `3n-1` on a common side.** ∎

**Census (FINITE-EXACT, A3).** Rules giving exactly `{SRC+C3, C3+SNK}`:

| quadruples compared | rules | of which pure-order |
|---|---|---|
| the two sheets, positive side | 0 | 0 |
| the two sheets, negative side | 0 | 0 |
| ν-related: `3n+1` at `x > 0` vs `3n-1` at `-x < 0` (both contracting sides) | 16 | 0 |
| ν-related: `3n-1` at `x > 0` vs `3n+1` at `-x < 0` (both expanding sides) | 16 | 0 |
| one sheet, both sides | 16 | 0 |

The ν-related hits come in three types of 8 rules each, 24 distinct rules in all:
* the contracting-side comparison has Type I + Type II;
* the expanding-side comparison has Type I′ + Type II.

| type | level pairs `A`, `B` | cross pairs `C`, `D` | classes (+,pos), (-,pos), (+,neg), (-,neg) |
|---|---|---|---|
| I | one order convention | **lineage gauge**: both arcs from `{s_o,t_o}` to `{s_e,t_e}`, or both back | e.g. SRC, STR, TT, SNK |
| I′ | opposite order conventions | lineage gauge | e.g. TT, SNK, SRC, STR |
| II | **label gauge** (`s_o -> s_e`, `t_e -> t_o`, or both reversed) | one order convention | e.g. SNK, SNK, SRC, SRC |

With diamonds present, `nu` acts as the converse on all four patterns for exactly the 8 Type II rules. Another 104 rules satisfy this vacuously, because all their classes are the self-converse TT and STRONG.

* **Type I realises the owner's picture.**
  * `3n+1` on its contracting side (positives) and `3n-1` on its contracting side (negatives) get opposite diamonds, swapped by `nu`.
  * But reversing the lineage gauge swaps which sheet gets which diamond, and the two expanding sides get TT and STRONG.
* **Type I′ does the same on the expanding sides** (`3n-1` on positives, `3n+1` on negatives), leaving TT and STRONG on the contracting sides.
* **Type II makes the class a function of the side alone.** The two sheets are identical on each side.
* **By THM-4057's guardrail**, an orientation is a gauge unless an intrinsic pairwise observable exists.
  * The order is intrinsic.
  * A lineage orientation of the cross pair `{s_o, t_e}` is not: both directions are equally admissible and give mirror assignments.
  * The owner's version therefore exists only in a gauge: **ANALOGY**.

**Pure-order rules (32; no gauge arc).**
* `(F_R)^op = F_{tau R}∘nu` exactly for all 32. `nu` reverses exactly the order arcs and `tau` exactly the dynamical arcs. In THM-511's language, on the arc cube `nu` flips the order coordinates, `tau` the dynamical coordinates, and the converse flips all of them.
* The 16 `rho`-symmetric rules orient `C` and `D` by the same convention. The uniform rule of Theorem A is one of them. For all 16, `nu ≅ id` and `tau ≅ converse` on classes.
* The 16 other rules produce no diamond at all.
* `nu ≅ converse` holds for **none** of the 32.
* On positive pairs the `(3n+1, 3n-1)` class pairs are `SNK/TT, SNK/STR, SRC/TT, SRC/STR, TT/SNK, STR/SNK, TT/SRC, STR/SRC` (2 rules each) and `TT/TT, STR/STR` (8 each). A diamond can sit on either sheet, but the two sheets are never the converse pair.

**Involution dictionary (PROVED for the families; FINITE-EXACT tables in A3–A4).**

| construction | negation `nu` | time reversal `tau` | converse `kappa = nu∘tau` |
|---|---|---|---|
| AM-fair crossing pattern (`3n±1`, contracting side) | fixes the class (via `rho`) | = converse (via `rho`) | C3+SNK ↔ SRC+C3 |
| AM-fair outward pattern (expanding side) | TT | TT | TT |
| parallel pattern (Mahler `ceil(3n/2)`, `floor(n/2)`) | TT ↔ STRONG | TT ↔ STRONG | fixes both self-converse classes |
| step quadruple `{n, 3n, 3n+b, (3n+b)/2}` (the two sheets agree on each side) | SRC+C3 → STRONG | SRC+C3 → STRONG | SRC+C3 ↔ C3+SNK |
| generic orbit windows | complements the order | reverses the path | = word reversal (Theorem C) |

## 4. Other natural quadruples (A4)

**4.1 The four up/down choices on a pair `{2i-1, 2i}`.** Each member moves by its length `i`.

| choice | AM-fair | positive fwd / bwd | negative fwd / bwd |
|---|---|---|---|
| odd up, even down = `3n+1` | yes | C3+SNK / SRC+C3 | TT / TT |
| odd down, even up = all-flipped = `3n-1` shifted (THM-4470 §3) | yes | TT / TT | C3+SNK / SRC+C3 |
| both up = Mahler `ceil(3n/2)` (THM-2228) | no | TT / STRONG | STRONG / TT |
| both down = `floor(n/2)` | no | STRONG / TT | TT / STRONG |

**4.2 All 12 placements of two disjoint directed arcs on four ordered points** (other pairs ascending). The diamonds occur only for the crossing matching `{13, 24}` with antiparallel arcs:
* `3 -> 1, 2 -> 4` gives C3+SNK (the Collatz pattern);
* `1 -> 3, 4 -> 2` gives SRC+C3.

**4.3 Step quadruple `{n, 3n, 3n+b, (3n+b)/2}`** (arcs `x3`, `+b`, `/2`; odd `n`). Both sheets give the same class:

| side | fwd/asc | fwd/desc | bwd/asc | bwd/desc |
|---|---|---|---|---|
| positive | SRC+C3 | STRONG | STRONG | C3+SNK |
| negative | STRONG | SRC+C3 | C3+SNK | STRONG |

Here the 3-cycle is literally `3n -> 3n+b -> (3n+b)/2 -> 3n`: the `+b` arc, the halving arc and one order arc, with `n` as the source. This is the closest reading of "3 and ±1". It is still sheet-blind, and its converse is `nu∘tau`.

**4.4 Inverse-tree fork `{T(x), x, 2x, E(x)}`.**
* There are three dynamical arcs.
* The class is SRC+C3 for odd `x`. For even `x` it is STRONG on positives and TT on negatives. This holds identically on both sheets (checked for `|x| < 3000`).
* `nu` is exactly order reversal, and converse = `nu∘tau` holds exactly.

**4.5 Length-4 orbit windows `{x, Tx, T^2x, T^3x}`.** The class is a function of the word, identical on the two sheets:

`000:S 001:S 010:S 011:C3+SNK 100:S 101:S 110:SRC+C3 111:TT`

The diamonds are the words `011` and `110`, reverses of each other (Theorem C at `m = 4`).

**4.6 `b = 0` and Mahler.**
* `T_0(x) = 3x/2` on odd `x`, over `Q`, has no AM-fair pairing, so the partner offset `beta` is a free gauge. The class is `beta·sgn(s_o)`: `beta = +1` puts the diamonds on the positive side, `beta = -1` on the negative side.
* AM-fairness is exactly what forces `beta = b` for `3n±1`.
* Mahler's `ceil(3n/2)` and its `nu`-conjugate `floor(3n/2)` give TT / STRONG on the positive side and STRONG / TT on the negative side, for either pairing.

**4.7 Multiplier test** (`(qn+r)/2` on odds, `10 <= i <= 300`):
* Every odd `q >= 3`, with `r = ±1`, gives C3+SNK / SRC+C3 on the pairing whose odd member is the lower one, and TT / TT on the other.
* `q = 1` gives STRONG / TT or degenerates.
* So **`H = 3` does not see the multiplier 3.** `q = 3` is special only through AM-fairness, which makes the order-reversing isomorphism the affine reflection `rho` and forces the offset.

**4.8 The literal kernel of "negation = converse".**
* On the cyclic triangle `Z/3` with arcs `x -> x+b`, `nu` maps `C_b` onto `C_{-b} = C_b^op`, label by label. So negation is the converse there. The proof is one line: `-(x+b) = (-x) - b`.
* But `C_b ≅ C_b^op` as classes. To turn this into the diamond swap, a fourth vertex must be a source for one sign and a sink for the other: a sign-dependent gauge. This is the same obstruction as in §3.

## 5. `H` and the Collatz graph

**5.1 Exact identities on the quadruple.**
* `H(F±) = 2 + sgn(b s_o)`, so on positive pairs `H = 2 + b`.
* The second Rédei digit `H mod 4 = 1 + 2 alpha_1` (THM-466) is the side bit.
* The "3" is `2 + 1` and counts one 3-cycle, which every `q >= 3` produces (4.7).

**5.2 Proposition D (PROVED; B0 checks 3000 random sets).** Let `V` be a finite set of integers `>= 3`. Put the map arcs `{x, T_b(x)}` inside `V` and the order on the other pairs. The result is the same for `b = +1` and `b = -1`: the ascending tournament with the halving pairs `{x/2, x}` inside `V` reversed.
* *Proof.* For odd `x >= 3`, `T_b(x) = (3x+b)/2 > x`, so that arc agrees with the order. For even `x`, `T_b(x) = x/2` on both sheets. ∎
* A map+order tournament therefore sees the sheet only through the choice of `V`. In the quadruple, `V` is the AM-fair pair, whose offset is `b`.

**5.3 Theorem C (PROVED; B1 checks every word, `m <= 12`).** Let `G_w` be the generic window tournament of a parity word `w` of length `m-1`:
* the path arcs are `j -> j+1`;
* for non-adjacent `j < k`, the arc is `j -> k` iff `y_k > y_j`, where `y_j = a_j log 3 - j log 2` is the multiplicative skeleton ("growth window goes up").

Then `G_{w^R} = r(G_w)^op` with `r(j) = m-1-j`: **reversing the word gives the converse tournament.** Hence `H(G_{w^R}) = H(G_w)` and `c3` agrees.
* *Proof.* `a^R_j = a_{m-1} - a_{m-1-j}`, so `y^R_j = y_{m-1} - y_{m-1-j}`. The reversed word's skeleton is the original read backwards in time and upside down. The path arcs match after `r` and `^op`. Ties never occur, because `log_2 3` is irrational. ∎
* The same proof works for Syracuse windows, with `y_j = j log 3 - K_j log 2`.

Generic-window `H` (B1):

| `m` | 4 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|
| max over words | 5 | 33 | 293 | 3377 | 53729 |
| mean over words | 4.0 | 17.5 | 82.7 | 477.8 | 3178.9 |

The minimum is always 1 (the word `1^(m-1)`, transitive). For comparison, random tournaments have `E[H] = m!/2^(m-1) = 233887` at `m = 12`. Windows are near-transitive.

**5.4 Crossings (B2–B5).**
* **Direction lemma** ([order-laws note](collatz_procgen_20260922_order_laws.md) §3.1, here in `T`-map form). A window's tournament equals its word's `G_w` except on crossing pairs.
  * Plus crossings are decay pairs that go **up**, i.e. forward-in-time arcs.
  * Minus crossings are growth pairs that go **down**.
  * 0 violations over 92,103 plus and 6,677 minus (window, crossing-pair) incidences, all lengths `m <= 64`, starts in `[3, 3+2^16)` and `[10^6, 10^6+2^16)`.
* **Where crossings sit.** Every crossing pair starts at a value `<= 4611`. Windows from near `10^6` cross only after descending there: none crosses for `m <= 24`, and 1 plus and 1 minus crossing appear at `m = 32`. Near the root the plus sheet has far more crossing windows (13,460 against 520 at `m = 64`): the 27-highway of the order-laws note §5.
* **Distributions (B3).** On complete residue systems `[10^6, 10^6 + 4·2^(m-1))`, `m = 6..14`, the distributions of sorted scores, `c3`, `phi` and `H` are **identical** on the two sheets (total variation 0). Near the root (`x0 = 3`) they differ, by terminal windows and crossings (TV 0.002–0.04).
* **What a crossing does (B4).**
  * `phi` (forward order arcs) moves by `+1` per plus crossing and `-1` per minus crossing.
  * `dH` is even (Rédei) and takes **both signs on both sheets**: T-window plus crossings mostly raise `H` (60 of 61 at `m = 16`); the 12 Syracuse lag-17 plus windows lower it or leave it (e.g. `171`: `222539` vs `343241`); minus crossings go both ways.
  * So `H` and `c3` do not carry the direction of the sign law.
* **Side separation (B5).** For `m <= 11`, `H(nu G_w) != H(G_w)` for every word except `101` and `01110`. So `H` of a generic window tells the two sides of 0 apart word by word. It is not a mod-4 law: all four pairs of residues mod 4 occur.

**5.5 Typing (PROVED from `kappa = nu∘tau`; B5 checks all words at `m = 7, 9`).**

| statistic | under `nu` (sheet swap; order reversed) | under `tau` (path reversed) | under converse | sheet-separating? |
|---|---|---|---|---|
| `H`, `c3`, OCF `alpha_k` | not invariant (side-aware) | `f∘tau = f∘nu` | invariant (THM-511 even) | only through crossings; direction not carried |
| scores | changes | changes | `s -> m-1-s` (odd) | only through crossings |
| `phi` (forward order arcs) | `#pairs - phi` | invariant | `#pairs - phi` | yes at crossings, with the sign of `b`: the sign law itself |
| class of a quadruple | invariant (Theorem A(c)) | converse | converse | yes, as the side bit `sgn(b s_o)` |

No window or quadruple statistic separates the sheets on a fixed side without the order. With the order, every separation is the sign law, as Corollary 7 of the inverse-tree note and the order-laws transfer theorem demand.

## 6. Rédei as a fixed-point parity, and Sperner

**6.1 Rédei through the OCF (PROVED in the repository; C0 re-checks the involution on 1147 tournaments).**
* THM-002 in the Grinberg–Stanley form: `H(T) = sum_sigma 2^psi(sigma)`. The sum runs over permutations `sigma` whose nontrivial cycles are directed odd cycles of `T`, and `psi` counts the nontrivial cycles.
* Let `X` be the set of pairs `(sigma, bits)` with one bit per nontrivial cycle. Then `|X| = H`.
* Flipping the bit of the cycle through the least moved point is an involution on `X`. Its only fixed point is `(id, ())`. Hence **`H ≡ 1 (mod 2)`: Rédei's theorem says that the identity is the only permutation without a nontrivial cycle.**
* The higher digits form THM-466's tower, `H ≡ sum_{k<m} alpha_k 2^k (mod 2^m)`.
* The twisted version, THM-582 and THM-647 ("anti-Rédei"): for an involutory anti-automorphism `phi`, the map `rev∘phi` is an involution on Hamiltonian paths, and `H ≡ #Fix(rev∘phi) (mod 2)`.
* Rédei's own proof (1934) is by induction. The classical arc-reversal route is summarised in THM-001. **Berge's proof: UNVERIFIED (not read, not used).**

**6.2 Sperner's lemma, the combinatorial core of Brouwer.**
* **Statement.** Every Sperner labelling of a triangulated `n`-simplex has an **odd** number of fully labelled cells.
* **Proof in dimension 2 (door-in / door-out).**
  * A door is an edge labelled `{1,2}`. Every cell has at most two doors, and a fully labelled cell has exactly one.
  * So the cells form paths and cycles, and their ends are the fully labelled cells together with the boundary doors.
  * The boundary side `1-2` carries an odd number of doors, and the handshake lemma finishes the proof.
* **Consequences.** Sperner implies the KKM lemma, which implies Brouwer. Following the doors is PPAD (Papadimitriou 1994).
* C1 checks 400 random Sperner labellings of a 14-subdivided triangle: rainbow counts odd, from 19 to 65. Every cell has at most 2 doors.

**6.3 Comparison.**

| | Rédei | Sperner | Collatz periodic points (`T` on `Z_2`) |
|---|---|---|---|
| counted objects | Hamiltonian paths | fully labelled cells | fixed points of `T^p` |
| pairing | involution on OCF configurations | door paths (degree <= 2) | `iota = Phi^-1∘(bit complement)∘Phi` |
| unpaired | exactly `(id, ())` | ends of paths: cells plus odd boundary | **none** (`iota` is free) |
| parity | odd | odd | **even** (`2^p`) |
| existence mechanism | parity | parity | Banach: the inverse branches `D`, `E` contract `|.|_2` by `1/2` |

The parity mechanisms are opposite. Rédei and Sperner leave one unpaired configuration. Collatz's 2-adic periodic points are paired off completely by a free involution that commutes with `T`, and their existence comes from contraction, not parity.

## 7. Collatz odd-count candidates and "chain growth" (C2–C4)

**(a) Preimage counts** `N_d(a) = |T^-d(a)|`, a function on `Z/3^d` (inverse-tree note, Proposition 9).
* `N_d` is odd on about 2/3 of the classes.
* The recursion increment `[a ≡ 2 mod 3] N_{d-1}(E(a))` is odd on about 0.22 of them.
* So there is **no parity law and no even-increment law**. Contrast Claim A for tournaments: `H(T) - H(T-v) = 2 sum_{C ∋ v} mu(C)` is always even (THM-002, THM-070), which is how Rédei parity survives vertex addition.

**(b) Cycle candidates.** Let `R_p(n) = #{w in {0,1}^p : x_w in [1, n]}`. Both parities occur for every `n` tested (`1, 2, 10, 100`; `p <= 16`), so there is no Rédei-type law.

**(c) `#Fix(T^p)` on `Z_2` is exactly `2^p` (C3, `p <= 16`, exact rationals).**
* Every periodic point is rational, `x_w = c_w/(2^p - 3^a)`.
* `iota(x_w) = x_{w-bar}` is a free involution commuting with `T`. Proof: the parity-vector map `Phi` conjugates `T` to the shift (Lagarias 1985; [sibling-dimension note](collatz_procgen_20260922_sibling_dimension_ladder.md) Lemma 1), and the complement commutes with the shift.
* So the count is even.
* The integer periodic points for `p <= 16` are exactly the five known cycles.

**(d) Hamiltonian "chains" of the type graph (C4).**
* The forward type graph mod `2^k` (`r -> T(r')` over both lifts) is the de Bruijn graph `B(2,k)` for `3n+1`, `3n-1` and `5n+1`, via Terras's word bijection. This is checked for `k <= 8`.
* Its Hamiltonian cycles are the binary de Bruijn sequences, `2^(2^(k-1)-k)` of them: `1, 1, 2, 16, 2048, ...`. This is checked by dynamic programming for `k <= 4` and by BEST / matrix-tree for `k <= 8`.
* The count grows doubly exponentially and is an even number for every `k >= 3`. It is **sheet-blind and drift-blind**: it belongs to the 2-shift, not to Collatz.

**(e) Growth and fluctuation.**

| | mean | second-moment ratio |
|---|---|---|
| tournaments (THM-589) | `E[H] = n!/2^(n-1)` (factorial) | `W(n)/n! -> 1` (`1.234` at `n = 8`, `1.050` at `n = 40`): self-averaging |
| Collatz preimages, Haar root | `E[N_d] = (4/3)^d` (exponential) | `1.622` at `d = 15` (all roots); `1.07` at `d = 4` rising to about `1.09` from `d = 9` on (units; reproduces the tree note's `1.070`/`1.089`) |
| independent-subtree 3-type model | `(4/3)^d` | `1.87` (all) / `1.26` (units) at `d = 15` |

* The D- and E-subtrees share 3-adic digits, which makes the true counts more concentrated than independent branching. The ratio does not decrease toward 1.
* **"`N_d` behaves like `H`" is false in growth and in fluctuation.**

**Verdict on "fixed point chain growth": ANALOGY.**
* The fixed-point lane's Banach fixed points `x_w` of inverse-branch words are exactly the `2^p` periodic points of (c). Their count is even, with a free involution: the opposite of Rédei.
* On the tournament side, the only exact statements attaching Rédei to Collatz data are Rédei itself, applied to tournaments built from Collatz data (quadruples, windows). There `H` is odd whatever the data are, and the second digit `H mod 4` is the side bit on quadruples (5.1).

## 8. Dictionary: what each phrase maps to

| owner's phrase | exact object | status |
|---|---|---|
| "tournament on four vertices" | the AM-fair pair quadruple with map arcs + order | REAL (Theorem A) |
| "the 2 isomorphism classes" | forward map → C3+SNK and inverse tree → SRC+C3, on each sheet's contracting side | REAL |
| "swap when all 6 arcs are inverted" | converse = time reversal (via `rho`), and exactly = negation ∘ time reversal | PROVED |
| "3" | `H(C3) = 3 = 2 + 1`; given by every odd `q >= 3` | ANALOGY |
| "±1" | the offset `b`: it fixes the AM-fair partner (`s_e = s_o + b`), hence the side carrying the diamonds, and it flips under inversion `(y-b)/3` | REAL as side bit / time direction; REFUTED as "`3n+1` = one diamond, `3n-1` = the other" |
| "Rédei" | `H` odd = one fixed configuration of the OCF involution | PROVED (repository) |
| "Hamiltonian path" | window `H`: a word function, sheet-blind on a side, side-aware; word reversal = converse | PROVED / FINITE-EXACT |
| "fixed point chain growth" | `2^p` periodic points (free involution); `(4/3)^d` preimages (not self-averaging); `2^(2^(k-1)-k)` de Bruijn cycles; versus `n!/2^(n-1)`, self-averaging | ANALOGY |

Related repository readings:
* THM-510's "B₂ atom" maps the four 4-tournaments to subsets of `{x+1, x/2}` for triangular numbers. It is a different correspondence, and nothing here depends on it.
* The mod-6 lane's "non-tournament verdict" concerned score profiles of cell counts. The quadruples here are genuine tournaments.

## 9. Sources (exactly what was read)

* **Wikipedia, "Tournament (graph theory)"** (raw wikitext read, 2026-09-24). Read: Rédei's existence theorem, cited as L. Rédei, "Ein kombinatorischer Satz", *Acta Litteraria Szeged* 7 (1934) 39–43. The odd-count form is taken from repository THM-001. **The original paper is UNVERIFIED (not read).**
* **Wikipedia, "Sperner's lemma"** (raw wikitext read). Read:
  * the statement (odd number of rainbow simplices) and the handshake proof;
  * the equivalence with Brouwer;
  * PPAD-completeness attributed to Papadimitriou.
* **Wikipedia, "Knaster–Kuratowski–Mazurkiewicz lemma"** (raw wikitext read): B. Knaster, C. Kuratowski, S. Mazurkiewicz, "Ein Beweis des Fixpunktsatzes für n-dimensionale Simplexe", *Fund. Math.* 14 (1929) 132–137, doi:10.4064/fm-14-1-132-137. Content beyond the citation is UNVERIFIED.
* **Crossref metadata** (bibliographic data only; the papers themselves were not read, so their content is UNVERIFIED beyond what the Wikipedia pages state):
  * E. Sperner, "Neuer Beweis für die Invarianz der Dimensionszahl und des Gebietes", *Abh. Math. Sem. Univ. Hamburg* 6 (1928) 265–272, doi:10.1007/BF02940617;
  * C. H. Papadimitriou, "On the complexity of the parity argument and other inefficient proofs of existence", *JCSS* 48(3) (1994) 498–532;
  * H. W. Kuhn, "Simplicial approximation of fixed points", *PNAS* 61(4) (1968) 1238–1242;
  * D. I. A. Cohen, "On the Sperner lemma", *J. Combin. Theory* 2(4) (1967) 585–587;
  * J. C. Lagarias, "The 3x+1 problem and its generalizations", *Amer. Math. Monthly* 92(1) (1985) 3–23. The conjugacy theorem is inherited from the sibling-dimension note, not re-read.
* **Wikipedia, "De Bruijn sequence"** (raw wikitext read). Read:
  * the count `2^(2^(n-1)-n)`;
  * attributions to C. Flye Sainte-Marie (*L'Intermédiaire des Mathématiciens* 1 (1894) 107–110) and to N. G. de Bruijn (*Proc. KNAW* 49 (1946) 758–764);
  * T. van Aardenne-Ehrenfest and N. G. de Bruijn, "Circuits and trees in oriented linear graphs", *Simon Stevin* 28 (1951) 203–217 (BEST).

  The primary sources are UNVERIFIED (not read); the count is re-verified computationally in C4.
* **Repository (read):**
  * THM-001, THM-002, THM-466, THM-510, THM-511, THM-582, THM-589, THM-647, THM-2228, THM-4057, THM-4469, THM-4470;
  * the inverse-tree note (Theorem 6, Corollary 7, Propositions 9–11);
  * the order-laws note (§§1, 3, 4, 5);
  * the sibling-dimension note (Lemma 1, cited only);
  * LEM-020 and the mod-6 cell-ordering note, for context.

  Two items were used only as stated in other files: THM-070 (Claim A), as stated in THM-002 Proof 2; and the Grinberg–Stanley / Irving–Omar derivation, as recorded in THM-002.

## 10. Reproduction

```bash
cd <worktree>
python3 04-computation/experiments/procgen_tourn_20260924_run.py   # writes 05-knowledge/results/procgen_tourn_20260924.out
# or one part at a time (stdout = sections, stderr = timing):
python3 04-computation/experiments/procgen_tourn_20260924_four_vertex.py   # A1-A4, ~1 s
python3 04-computation/experiments/procgen_tourn_20260924_windows.py       # B0-B5, ~55 s
python3 04-computation/experiments/procgen_tourn_20260924_redei.py         # C0-C4, ~8 s
```

* The full run takes about 65 s, one process at a time. Peak RSS is about 530 MB (the `N_15` arrays and the `m = 64` window batches).
* All randomness is seeded, and two consecutive runs produced byte-identical output.
* Each script ends with `ALL CHECKS PASSED`, and every check is an explicit `raise`.
* Scratch explorations (rule inspection, `dH` signs, side-equality words) live in `scratch/procgen_tourn/`. They are not deliverables; their conclusions are re-checked inside the scripts.
* SHA-256 (raw LF bytes):
  * `four_vertex.py` `bb042ed26b3c3fc1d567395ba9532423c4c2cc87a662a5ac676e32fafa8186c9`
  * `windows.py` `636bd54101940ae476adb30c74a0ab1d3c59af079b3456cbe52e7cc5d07c0f39`
  * `redei.py` `29fb358a92ed6784ff3b4ab505c70ea6652ad13ada821a6ae752637b5467d513`
  * `run.py` `f4846704d351cb4c8275f0e031a2f6a25b8a6ab153a5fcd266cb4a5a8f8f224e`
  * `procgen_tourn_20260924.out` `ebe55081f3f18391b34e089387d3fc84b794a0d6f8abff50bc5774a3d851987d`
