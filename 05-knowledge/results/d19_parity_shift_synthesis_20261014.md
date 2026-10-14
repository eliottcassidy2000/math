# Guy D19, Giuga and Korselt, and the parity / shift-by-one theme (THM-4650, THM-4651)

Session `mac-mini-2026-10-14-wilson`, second half, 2026-10-14. Not yet independently audited.

**Owner's brief (verbatim):**

> "keep relevant jobs running and ideas in mind as you make this session double length by now considering the following: richard guys unsolved number theory problems D19, and additionally think how guiga numbers and korselts criterion relate with guigas conjecture, which revolves around a single value which is both, see how this duality relates deeply with twin primes and even vs odd goldbach. again think in terms of linked lists and a number sometimes embodying the qualities of one greater or less than it. think about prior work on the multiplication/division even/odd half/whole and rational/irrational dualities (as in dedekind cuts) and think about the real/transcendental duality and prior work regarding irrationality exponents and consider elegant solutions to the continuum hypothesis, seeing how syncing together with prior work and bott periodicity and maybe an assumed partial or extended riemann hypothesis we can prove collatz or hadwiger nelson or the planar jacobian by seeing how number theory permeates through the structure of everything and discrete graphs and geometry/topology/analysis/abstract algebras just embody those structures through their lenses in a predictable way, think prior equidistribution, equidecomposability etc related work connects get creative, be free in your reasoning, pursue most relevant mathematical proof progress goals identified as they continue to appear"

**Concurrent work.** codex-tiling is executing the same brief:
* `giuga_korselt_defects_20261014f.md` and `giuga_block_lattice_20261014f.md`: signed Giuga defects, primary-pseudoperfect completions, a mod-8 lattice;
* `collatz_successor_guards_20261014f.md`;
* `continuum_finite_witness_boundary_20261014f.md`: forcing cannot create finite witnesses, and a mod-24 local guard for D19.

This note avoids repeating them and records what is new here.

Typing: PROVED, KNOWN, CITED, FINITE-EXACT, HEURISTIC and COINCIDENCE as usual. Analogies are marked as analogies.

## 1. Guy D19 (THM-4650)

D19 asks for a point at rational distance from all four corners of the unit square. It is open.

* **Rationality.** A D19 point has rational coordinates: `x = (a² − b² + 1)/2` and `y = (a² − c² + 1)/2`. So D19 is about rational points on a `(Z/2)⁴`-cover `S` of the plane.
* **Special loci (PROVED, THM-4650).** No D19 point lies on:
  * a diagonal (`√2` is irrational);
  * the circumcircle (Ptolemy gives `PC = PA + √2·PB`);
  * a side line: the conductor-40 curve `y² = x³ − 2x + 1`, of rank 0;
  * a midline: the conductor-32 curve 2-isogenous to `y² = x³ − x`.
* **The midline is Fermat's right-triangle theorem.** A midline point would make 1 a congruent number (Guy D27). This is Fermat's first infinite descent, now sitting inside a rational-distance question about the unit square.
* **Rectangles (PROVED reduction).** On the midline of a `1 × a` rectangle, D19 is exactly a rational triangle with base `2a` and height 1. That is the base/height-ratio problem of Guy's D19 remarks (MacLeod), and the square is ratio 2.
  * The Jacobians have `j = 16(16a⁴ + 16a² + 1)³/(a²(a² + 1))`.
  * Rank-1 ratios include `1/3, 5, 6, 8`; ratios 1, 2, 3, 4 have rank 0 (FINITE-EXACT, 46 ratios).
* **Where the open problem sits (HEURISTIC).**
  * Dropping one corner gives Berry's K3 surface, on which three-distance points are dense.
  * With four corners, `S` is a complete intersection of four quadrics in `P^6`, with `K = O(1)` before resolving the two singular points over the circular points `I`, `J`. Over each of them, four branch lines meet. The local model `w_k² = ℓ_k(u, v)` is a cone over a quartic elliptic curve, i.e. simple elliptic singularities of degree 4. A rough count gives `K² = 16 − 4 − 4 = 8 > 0`, so general type is expected.
  * Bombieri–Lang would then confine all but finitely many D19 points to finitely many curves of genus 0 or 1, as in the conditional Erdős–Ulam theorems (Shaffaf; Tao).
  * THM-4650 eliminates the natural genus-1 curves. A complete list of genus-0/1 curves on `S` would turn D19 into a conditional finite search.
* **Rational versus irrational (PROVED in THM-4650; reading).** The special loci are guarded by the two oldest proofs in number theory.
  * The irrationality of `√2` guards the diagonals and, through Ptolemy, the circumcircle.
  * Fermat's infinite descent, that 1 is not a congruent number, guards the midlines.
  * D19 asks for a rational shadow of a square whose own diagonal is irrational.
* **Irrationality exponents and rational points (KNOWN: Vojta's dictionary, 1987).**
  * Vojta's conjectures put two kinds of statement on one footing:
    * Diophantine approximation: Roth's theorem says every algebraic irrational has irrationality exponent 2;
    * finiteness of rational points: Faltings for curves of genus at least 2, and Bombieri–Lang for varieties of general type.
  * D19's open part is a Bombieri–Lang question for `S`.
  * The Collatz cycle question is an approximation question for `log₂3` (THM-4633/4634; HYP-9263 asks whether `μ(log₂3) = 2`). That number is transcendental (Gelfond–Schneider), so Roth gives nothing and only Baker-type bounds are available.
  * Both are Diophantine approximation in Vojta's sense, on opposite sides of the algebraic/transcendental divide.
* **Link to Hadwiger–Nelson.** In the language of `Q(i)`, `N(z) = □` means `z = s·w²` (Hilbert 90). So D19 asks for `z ∈ Q(i)` such that `z`, `z − 1`, `z − i` and `z − 1 − i` are each a rational multiple of a Gaussian square.
  * The 5-chromatic unit-distance graphs of this session (§5) are built from rotations `e^(iθ)` with rational `cos θ` (`5/6, 7/8`), i.e. norm-1 elements of quadratic fields.
  * Both problems are statements about multiplicative groups of norm-1 or norm-square elements meeting additive translates.

## 2. Giuga, Korselt, Wilson (THM-4651)

* **"A single value which is both" (KNOWN).**
  * Korselt at `p` is `n ≡ 1 (mod p − 1)`. Giuga at `p` is `n/p ≡ 1 (mod p)`.
  * Both hold at every `p | n` iff `n ≡ p (mod p²(p−1))` for every `p | n`. That is the counterexample condition for Giuga's conjecture.
  * Primes satisfy it trivially.
* **The parity divide (PROVED).**
  * An even squarefree `n` has Korselt score exactly 1, so the Korselt world is odd. Every known Giuga number is even.
  * A counterexample would have to be an odd Giuga number. Then `ω(n) ≡ k + 1 (mod 2)`, so with index 1 it has at least 10 prime factors (elementary).
  * Census to `10^8`: no Carmichael number has more than half its primes Giuga (THM-4651 (4)).
* **The p²-lift (KNOWN + PROVED consequences).**
  * `S(p) ≡ p + (p−1)! ≡ p − 1 + p w_p (mod p²)`.
  * Wilson primes are where `S(p) ≡ p − 1`. Giuga's congruence holds at `p²` only for the defect-1 near-Wilson primes: only 1277 below `2·10^5`.
  * This connects THM-4641's Wilson dictionary to Giuga's conjecture.
* **Linked lists under the arithmetic derivative (KNOWN).**
  * `p ↦ 1`; primary pseudoperfect `N ↦ N − 1`; Giuga `g ↦ g + 1`; `p^p ↦ p^p`.
  * A primary pseudoperfect number times a prime neighbour is a Giuga number: `6·5 = 30`, `42·41 = 1722`, `47058·47057` (Sondow–Tsukerman; codex).
  * This, not "Giuga numbers are twin-prime centres", is the twin-prime mechanism. The latter holds only for 30, 858, 1722 (COINCIDENCE).

## 3. Twin primes, even versus odd Goldbach, and the parity barrier

* **KNOWN.**
  * Selberg's parity problem: sieve weights cannot distinguish integers with an even number of prime factors from those with an odd number.
  * This is why sieve methods alone cannot produce twin primes or the binary (even) Goldbach theorem. Bounded gaps (Zhang; Maynard) stop at 246, at 12 under Elliott–Halberstam and at 6 under GEH; the parity barrier blocks gap 2.
  * The ternary (odd) Goldbach problem escapes: with three primes, the circle method's main term dominates. Vinogradov proved it for large odd `n` and Helfgott completely (2013).
  * The odd problem is an analysis problem; the even one is a parity problem.
* **The divide in this session (PROVED vs analogy).**
  * Korselt lives on odd `n` and Giuga (as far as known) on even `n`, and a Giuga counterexample must cross the divide. That is PROVED (THM-4651).
  * Reading it as "the same" parity barrier as Goldbach's is an analogy only.
* **Collatz.**
  * The Collatz map is itself a parity machine; the parity vector is the orbit (THM-4634).
  * The corpus has density and almost-sure theorems for shifts:
    * `n` and `n + 1` merge almost surely (THM-4581);
    * `n` and `n + h` merge for a density-one set (THM-4627);
    * −1 switches cover almost every trailing-ones source (THM-4640 (8)).
  * In each case the exceptional set is null or of density zero. Covering every integer would prove Collatz.
  * **Primes and twin primes inherit the merging (PROVED corollary; FINITE-EXACT data).**
    * By THM-4627 (2)(a), the set of `n` for which `n` and `n + h` merge by Terras time `K` (equal values, equal odd counts) is a union of residue classes mod `2^K`. Its density `μ_h(K)` tends to 1.
    * Primes are equidistributed over the odd classes mod `2^K` (prime number theorem in progressions), and the odd classes have density at least `1 − 2(1 − μ_h(K))` in it. Hence, for every `h ≥ 1`, **almost all primes `p` (relative density 1) have Collatz orbits of `p` and `p + h` that merge**.
    * The same holds for twin-prime pairs `(p, p + 2)` **conditionally** on the Hardy–Littlewood conjecture for twin primes in progressions mod `2^K`. All odd classes are admissible and have equal local factors at 2.
    * Data (`prime_shift_merging_1e7.out`). The fractions merging before either orbit reaches 1, below `10^7`:
      * `h = 2`: 0.31485 (all odd `n`), 0.31378 (primes), 0.31540 (58,980 twin pairs);
      * `h = 1`: 0.37007 (all odd `n`), 0.36922 (primes).
    * Primes and twin primes are statistically indistinguishable from generic odd numbers here. The approach to 1 is very slow, like THM-4581's `T^(−1/2)` law.
    * In the brief's language: almost every prime "embodies" the Collatz future of its neighbour.
  * This is the corpus's version of the parity barrier: measure and density methods are blind to the null set where a single counterexample could live (HEURISTIC parallel, not a theorem).

## 4. Shift by one: the common thread

* **Consecutive smooth numbers (KNOWN).** Störmer: for a finite prime set `S`, only finitely many consecutive pairs `(n, n+1)` are both `S`-smooth. For `S = {2, 3}` they are `(1,2), (2,3), (3,4), (8,9)`. These are exactly Gersonides' pairs, the `{2, 4, 8}` of THM-4640/4641: the repair relations at −1 and the unit-denominator Collatz cycles.
* **Gaps (KNOWN).**
  * Catalan–Mihăilescu (`8, 9`) and Tijdeman (gaps between `{2,3}`-units tend to infinity, effectively) control how close `2^A` can come to `3^k`.
  * That closeness is the cycle question (THM-4633: `k ≥ 72,057,431,991` odd steps from Eliahou's budget and Barina's verification).
* **abc (KNOWN conjecture).** `n` and `n + 1` cannot both be highly composite in the radical sense.
* **The "−1" characterizations of primes.**
  * Wilson: `(n−1)! ≡ −1`. Giuga: `S(n) ≡ −1`. Agoh: `n B_(n−1) ≡ −1`.
  * All hold for primes. At `p²` their defects are governed by the Wilson quotient (THM-4651 (1)).
* **Hadwiger–Nelson.** The unit distance is "one greater" in the metric. Rotations with rational cosine `(2N−1)/(2N)` are the hinges, with rotation field `Q(√−(4N−1))` (corpus rotation-field law). `N = 3` gives the Moser spindle and `N = 4` gives Parts' `ρ`.

## 5. Hadwiger–Nelson progress this session (SOLVER verdicts)

Details are in `hn_type_m_hinge_rotations_and_forced_radii_20261014.md`.

* **Observed rule (SOLVER; no proof).**
  * At the hinge radius 2 (`t = 4`), among Parts' 25 base pairs, the type-M graph is 5-chromatic exactly when one piece is rigid (`V49 ⊕ V37` or `V49 ⊕ V49`, with very few centre-coloured patterns) and the other forces radius 2.
  * The six spindle-free graphs tested behave the same way.
  * At the nine other radii tested (225 graphs), no type-M graph was 5-chromatic, even where both pieces force the radius (`t = 4/9`). Forcing alone is not enough, and rigidity at those radii was not measured.
* **Spindle-free.**
  * Rebuilt exactly: the Voronov–Neopryatnaya–Dergachev graph (64,513 vertices, 542,472 edges, not 4-colourable).
  * Its hinge `ψ = φ_0^(−1)ρ` is Parts' `ρ`.
  * Replacing one side by the clock piece `M_2 ⊕ M_1^clock` gives the **48,121-vertex** spindle-free graph `M_3 ∪ ψ(M_2 ⊕ M_1^clock)`, which is still not 4-colourable (CaDiCaL, 1,514 s).
  * This is the first reduction step towards the 1,441-vertex record (Heule 2021).
* **Type M with Parts' pieces.** Vertex-critical subgraphs of 611–618 vertices, against the record 509.

## 6. Continuum hypothesis, Bott periodicity, Riemann hypotheses: what can and cannot be leveraged

* **CH (KNOWN, plus codex).**
  * CH is independent of ZFC (Gödel; Cohen). codex's finite-witness note proves that forcing cannot create or destroy finite witnesses.
  * Collatz counterexamples, finite 5- or 6-chromatic unit-distance graphs (de Bruijn–Erdős), JC2 counterexamples of fixed degree and D19 points are all finite witnesses. So no choice of CH, ¬CH or a forcing axiom can decide them.
  * Axiom-sensitivity is real only for non-finite variants. Examples: the measurable chromatic number of the plane (Falconer: at least 5), and Shelah–Soifer graphs whose chromatic number depends on choice.
  * "Elegant solutions" of CH (Freiling's symmetry axiom, Woodin's programmes) therefore cannot be the missing ingredient for these four problems.
* **Bott periodicity (KNOWN, as analogy).**
  * The period-8 phenomena that genuinely appear are arithmetic: `χ_8(p) = (2/p)` in the tournament count (THM-4641 row 5), Milgram's Gauss-sum formula (signature mod 8), and codex's mod-8 lattice.
  * Real Bott periodicity (Clifford modules, `KO`) shares the 8 through the Witt group and signature mod 8. No mechanism transfers a topological periodicity to Collatz, Hadwiger–Nelson or JC2 here.
* **Riemann-type hypotheses (KNOWN conditional results).**
  * GRH gives Artin's conjecture (Hooley 1967), hence infinitely many primes with 2 as a primitive root (THM-4641 row 7's pattern recurs). It also gives ternary Goldbach for large odd `n` (Hardy–Littlewood; now unconditional).
  * Bombieri–Lang (not RH) gives the conditional Erdős–Ulam and D19-type finiteness of §1.
  * abc gives many finiteness statements for `{2,3}`-type equations.
  * **No known form of RH implies Collatz, Hadwiger–Nelson or JC2.** The obstructions are respectively the Diophantine/measure-blindness of §3, the finite constructions of §5, and the Newton-polygon combinatorics.
  * **Where Riemann-type statements already act in this corpus (PROVED there).**
    * THM-4628: Weil's Riemann hypothesis for curves over finite fields bounds the Kloosterman sums. These are the eigenvalues of the residue graphs of the field plane `Q(ζ_55)`. So residue colourings at a fixed prime with residue field `F_(q²)` need at least `1 + (q+1)/(2√q)` colours (Hoffman), and none uses six or fewer.
      * This is a genuine RH-to-Hadwiger–Nelson implication. It applies to the restricted class of residue colourings, not to all colourings of the plane.
    * THM-4520: the Collatz level operator on `Z/(2·3^(n−1))` has characteristic polynomial `(2^L − 1)λ^L + λ^(L/2) − 1`, with every eigenvalue on two circles of modulus about `1/2`. This is a circle law of Ramanujan/RH type for a finite Collatz model.
    * Such spectral statements give mixing and equidistribution, i.e. "almost all" results. Like §3, they do not reach every integer or every colouring.
  * **JC2 and arithmetic (CITED from memory; checked by audit B).**
    * Tsuchimoto (2005) and Belov-Kanel–Kontsevich (2007) prove `JC_(2n) ⟹ D_n`: the Jacobian conjecture in dimension `2n` implies the Dixmier conjecture for the `n`-th Weyl algebra. The proof reduces modulo primes and uses the Frobenius `p`-centre.
    * With the classical `D_n ⟹ JC_n`, this makes the two conjectures stably equivalent. For the plane it gives `JC_2 ⟹ D_1` and `D_2 ⟹ JC_2`, **not** an equivalence of `JC_2` with `D_1`.
    * Reduction mod `p` and the Frobenius `p`-centre are the clearest case of number theory permeating JC: they turn an analytic statement into a Frobenius-twisted algebraic one. They give implications, not a proof.
* **Equidistribution and equidecomposability (KNOWN).**
  * Laczkovich's circle squaring, and Marks–Unger's Borel version, need discrepancy bounds for translation actions: Weyl equidistribution with an irrationality-measure input.
  * This is the same input as the corpus's irrationality-measure dictionary (THM-4633/4634: `μ(log₂3)` controls Collatz cycles).
  * Equidistribution gives "almost all" theorems. As in §3, it does not by itself reach every integer.

## 7. What would be proof progress next

1. **Hadwiger–Nelson, spindle-free record.** Core-guided minimisation of the 48,121-vertex graph and further rigid/forcing reductions, aiming below 1,441. This is the most realistic notable result in view.
2. **D19, conditional.** Resolve the singularities of `S` rigorously (general type), and classify its genus-0/1 curves.
3. **Giuga.** The parity divide suggests targeting "odd Giuga numbers do not exist" with index constraints. The index-1 parity bound is elementary; literature bounds are far stronger. No realistic proof in view.
4. **Collatz.** The rate in THM-4640 (8) and the structure of the exceptional null set (see the first note of this session).
