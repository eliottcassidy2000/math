---
id: THM-4804
title: "Complement lines shatter every module, and the carousel fibre is the Frobenius action on a normal basis. (1) For every tournament T on n >= 3 vertices, every Hamiltonian path P and every module M of T with 2 <= |M| <= n-1, M is not a module of the line partner K_P(T) = T^op xor E(P) (proof: the path leaves M at p_j -> x = p_(j+1); M beats x in T, so x beats all of M in T^op; the partner restores only the path pair, and by the no-immediate-re-entry lemma (p_(j+2) is outside M) every other pair {y, x}, y in M, keeps x -> y, so x sees M non-uniformly). Hence the two ends of a line share no proper non-singleton module, their strong-module trees share no internal node, the LCA structure of session56's arc interface is rebuilt by every line, and no substitution fibre Q[T_1..T_k] with a block of size >= 2, no stack with a non-singleton block and no cyclic-triple fibre of a non-prime tournament is closed under complement lines (the owner's caution, in its strongest form). (2) Normal-basis dictionary for the carousel fibre T_w = Q_m[A_(w_i)] of cyclictriple_composition_20261026: with a normal element alpha of F_(2^m), the Frobenius x -> x^2 rotates the coordinate word, so the necklace classes are the Frobenius orbits, the least period of w is the degree of x_w = sum w_i alpha^(2^i) over F_2, the surviving symmetry Aut(T_w)/C_3^m = C_(m/d) is Gal(F_(2^m)/F_2(x_w)), the period-d class counts are Gauss's counts of irreducible polynomials of degree d, and the presentation count is proportional to the degree: tc(T_w) = H_m d/(m 3^m) (F4 path-cover invariance), every presentation black; the same holds over F_q for q child types, q a prime power, and the one-type case is the F_(1^m) structure of Kapranov-Smirnov with the full rotation group as Galois group. (3) FINITE-EXACT at n = 12 (all 1,851,660 presentations of the four carousel classes): no old module survives in any partner, no partner is in the fibre, 98.6% of partners are prime and the rest have a module of size 2, every presentation is black, and the tiling-to-partner-class map is injective (5715, 17145, 17145, 5715 distinct partner classes = tc), with 45,424 partner classes in all, 296 of them shared between fibre nodes (multiplicities 27, 54, 81, 108)"
status: >
  PROVED: (1) (Theorem 1 of the note with its lemma and corollaries); (2) with the normal basis theorem, Gauss's
  irreducible-polynomial count, Theorem C / formula (2) of cyclictriple_composition_20261026 and F4 of
  tournament_composition_classification_20261026 as CITED inputs.
  FINITE-EXACT: (3), carousel_line_partners_20261026.cpp (+ _m3.out, 14 s): literal check of (1) on every
  presentation, primality by pair-closure, canonical forms by colour refinement with 128-bit words.
  CITED: THM-643 (tilings), THM-4729 (lines), THM-4731 (black multiplicities), THM-4732 (K_P conventions), THM-4741
  (tc = H/|Aut| counts tilings); session56_rmq_lca_connections_20261025 (strong-module tree, LCA arc interface);
  tournament_composition_classification_20261026 (decorated-tree lemma, stabilizer formula (2.1), F4);
  cyclictriple_composition_20261026 (Theorems C-G; its block-cut metagraph is a different operation, as it states);
  Boussaïri-Chergui-Ille-Zaidi and Boussaïri-Ille-Lopez-Thomassé through those packets; the normal basis theorem;
  Gauss's formula; Kapranov-Smirnov's F_(1^m) = {0} u mu_m convention.
  CONCURRENT: the two packets above are same-week Codex/klein lanes; this entry cites and does not modify them; the
  owner's finite-field sentence is stated and proved here, and was not found in the repository text.
  OPEN: HYP-9414 (injectivity of the tiling-to-partner map for all odd m and in general; the pair-module partners;
  connectivity of the fibre's line neighbourhood; the orbit-size law for every decorated fibre).
session: opus-2026-10-26-S59
source: 05-knowledge/results/module_shattering_lines_carousel_finite_field_20261026.md
depends_on:
  - 01-canon/definitions.md
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
  - 01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md
  - 05-knowledge/results/session56_rmq_lca_connections_20261025/tournament-hierarchy/REPORT.md
  - 05-knowledge/results/cyclictriple_composition_20261026/README.md
  - 05-knowledge/results/tournament_composition_classification_20261026/CLASSIFICATION.md
related:
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
  - 01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md
  - 05-knowledge/results/tournament_composition_classification_20261026/fibre/RESULT.md
  - 05-knowledge/hypotheses/HYP-9414-module-shattering-consequences-injective-fibres-pair-module-partners-and-the-orbit-size-law.md
scripts:
  - 04-computation/experiments/carousel_line_partners_20261026.cpp (+ _m3.out)
---

# THM-4804: lines shatter modules; the carousel fibre is a Frobenius action

Conventions: THM-4732's `K_P(T) = T^op xor E(P)` (the packets' `tau`), blue/black presentations, `tc = H/|Aut|`; modules, strong modules and the tree `D(T)` as in session56; the carousel `Q_m`, the blocks `A+ = (C3 -> 1)`, `A- = (1 -> C3)` and `T_w = Q_m[A_(w_i)]` as in `cyclictriple_composition_20261026`. Proofs in the [session note](../../05-knowledge/results/module_shattering_lines_carousel_finite_field_20261026.md), sections 2-4.

**(1) Module shattering.** *Lemma.* If `P` leaves a module `M` at `p_j -> p_(j+1)`, then `p_(j+2) notin M`: every vertex of `M` beats `p_(j+1)`, so `p_(j+2) in M` would contradict `p_(j+1) -> p_(j+2)`. *Theorem.* With `x = p_(j+1)`: in `T` the module `M` beats `x`; in `T^op` `x` beats `M`; `K` restores `p_j -> x` and, since `x`'s other path neighbour `p_(j+2)` is outside `M`, leaves `x -> y` for every other `y in M`; so `x` sees `M` non-uniformly and `M` is not a module of `K` (the entering case is symmetric). *Corollaries.* `K_P` is an involution, so the module lattices of the two ends of a line meet only in singletons and `V`; `D(T)` and `D(K)` share no internal node; the quotient at `LCA(u, v)` that explains an arc in `T` is never transported; any family with a common proper non-singleton module pattern is closed under no line.

**(2) Dictionary.** `x_w = sum w_i alpha^(2^i)` for a normal element `alpha`; `x_w^2` has the rotated coordinates; `x_w in F_(2^d)` iff `w` is `d`-periodic; so necklaces = Frobenius orbits, least period = degree, `C_(m/d) = Gal(F_(2^m)/F_2(x_w))`, `P_2(d)` = irreducible polynomials of degree `d`, and `tc(T_w) = H_m/(3^m . m/d) = (H_m/(m 3^m)) d` because `H` is constant on the fibre (F4); `bm(T_w) = 0`. Over `F_q` (`q` a prime power, `q` child types) the same; for one child type the class is the `F_(1^m)` structure with Galois group `C_m` and `tc = H/(|Aut(A)|^m m)`.

**(3) Order 12.**

| class | presentations (black) | partners prime | pair-module partners | old modules surviving | in the fibre | distinct partner classes (= tc) |
|---|---|---|---|---|---|---|
| `T_000` | 462915 | 456435 | 6480 | 0 | 0 | 5715 |
| `T_100` | 462915 | 456462 | 6453 | 0 | 0 | 17145 |
| `T_110` | 462915 | 456462 | 6453 | 0 | 0 | 17145 |
| `T_111` | 462915 | 456435 | 6480 | 0 | 0 | 5715 |

`45,424` partner classes in all (`44,790` prime), reached by `27, 54, 81, 108` presentations (`33846, 148, 11282, 148` classes): `148` partner classes are common to a converse pair and `148` to a `|Aut| = 81` class and a `|Aut| = 27` class.
