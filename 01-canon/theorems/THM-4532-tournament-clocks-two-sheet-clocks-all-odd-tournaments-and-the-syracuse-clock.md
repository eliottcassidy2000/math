---
id: THM-4532
title: "Tournament Clock Prime Collatz (the owner's program), made concrete. Clock digraphs Clk(G, C) on finite abelian groups (loops, doubled arcs and missing pairs allowed) obey a reflection law: every Cayley digraph of an odd-order abelian group is arc-even, and H = nF (mod 2n). Two-sheet clocks D_r (two rotational r-clocks, the second turned a quarter and run backwards) have every arc on an odd number of Hamiltonian paths for every odd r <= 19. This gives all-odd tournaments at N = 30, 34, 38; N = 34 and 38 are not reachable by Paley, since 35 and 39 are not prime powers. D_r is anti-circulant, D_3 = QR_7 - v, D_5 = QR_11 - v, D_7 = QR_127[mu_14]. The Syracuse map is a clock: log_2 S(A) = 2(A mod 3) - v (mod 6) in (Z/9)^x"
status: >
  PROVED + INDEPENDENTLY AUDITED:
  - the clock definitions and Proposition 1.1 (Paley minus a vertex is a two-sheet clock in discrete-log coordinates);
  - L1, the reflection law (THM-4529 Theorem 5.1 reached independently: Cayley digraphs of odd-order abelian groups are
    arc-even, loops and doubled arcs allowed), plus its corollary H = nF (mod 2n);
  - T1, the vertical-arc theorem for reversed-sheet two-sheet clocks with symmetric cross set;
  - the quadrant and circle forms of D_r; D_r's symmetries (Z/r automorphisms, self-converse, and a single-2r-cycle
    anti-automorphism, so D_r is anti-circulant);
  - the Syracuse clock step law mod 9.
  PROVED (lane; read): L2 (lexicographic decomposition), L3 (power-residue clocks mod 3^k are blown-up cycles), L4
  (multiplicative combination), L5 (complement duality; mod 2 = Berge's Stronger Theorem, CITED from Schweser-Stiebitz-Toft
  arXiv:2510.10659); the doubling law and the square-clock IFS for Tao's Syracuse measure (exact reformulations, no new
  bound).
  FINITE-EXACT:
  - D_r is all-odd for every odd r <= 19: independently audited for r <= 13 (N <= 26); r = 15 (N = 30) and r = 17
    (N = 34) by two independent inclusion-exclusion engines (lane, reproduced on re-run); r = 19 (N = 38) by one engine run
    (karp4, 14.5e9 subset representatives; a second-engine confirmation is still owed);
  - the N = 14 two-sheet census (exactly one all-odd class, D_7).
  REFUTED: the lane's own draft one-term parity formula for lexicographic products (corrected in the note).
  OPEN:
  - C1, D_r is all-odd for every odd r, which would settle the existence half of HYP-9167(b) (folded into HYP-9170);
  - the non-existence side at N = 13;
  - why the apex completions Q_r are all-even for r <= 9.
  ANALOGY: the owner's triplet principle (p^2qr = "p and p^3 combined"). It is exact only where it reduces to Z/3 with
  2 = -1 or to the parity split of 2^(-v) weights.
source: collatz-procgen-20260922 session, tcpc lane (2026-10-01), the owner's "Tournament Clock Prime Collatz" proposal; audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md
  - 01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md
related:
  - 05-knowledge/hypotheses/HYP-9167-paley-minus-a-vertex-is-redei-rigid.md (existence at N = 30, 34, 38)
  - 05-knowledge/hypotheses/HYP-9170-all-odd-anti-circulant-tournaments-exist-for-every-n-2-mod-4.md (D_r as the candidate family)
note: 05-knowledge/results/procgen_tcpc_20261001_tournament_clock_prime_collatz.md
scripts: 04-computation/experiments/procgen_tcpc_20261001_{run,lib}.py, procgen_tcpc_20261001_{hp,par26,karp3,karp4,karp2}.c
script_audit: 04-computation/experiments/procgen_tcpc_20261001_orchestrator_check.py
output: 05-knowledge/results/procgen_tcpc_20261001.out (--long --n26 --n34), procgen_tcpc_20261001_n38.out
output_sha256: 5a76471fbb6d432e6103e5194812f33007f083d70435d1f45d42ad47869d4d02 (n38 file: 983db6ecd24ca20fc3b996ae03bfd479b0642e9fd89ae6abf0c9d41cc30788f4)
output_audit: 05-knowledge/results/procgen_tcpc_20261001_orchestrator_check.out
output_audit_sha256: 30335e9126bcee5a2315051d9e7e0474ee88ba93d59379166b4eca01ce49ee8e
hash_basis: raw bytes
audit: >
  Independent code (procgen_tcpc_20261001_orchestrator_check.py; the orchestrator's own exact and mod-2 engines; the lane's
  code was not read), 16 checks:
  - D_r for r = 3..13 is a tournament with translation automorphisms and a single-2r-cycle anti-automorphism;
  - all arcs are odd (exact for r <= 9, with H(D_7) = 24540117 and H(D_9) = 116670839805; mod 2 for r = 11, 13);
  - D_3 = QR_7 - v, D_5 = QR_11 - v, D_7 = QR_127[mu_14], and D_9 is not QR_19 - v (networkx isomorphism);
  - T1 on 375 vertical arcs of 75 random symmetric two-sheet clocks;
  - the interval-control pattern;
  - the mod-9 step law on all 133333 odd A < 4e5 prime to 3.
  D_9's H equals the unexplained second all-odd anti-circulant class at N = 18 in THM-4529's census, so that class is D_9.
  The lane's runner was re-run with --long --n34 (1657 s, 671799 checks, ALL CHECKS PASSED). D_15 (N = 30) by karp3 and karp4,
  and D_17 (N = 34) by karp3 on T and karp4 on T^op, were reproduced. The output is identical to the lane's apart from the --n26 lines
  not requested. The N = 38 run (karp4 only, 72 minutes) was not repeated.
---

# THM-4532 — tournament clocks, two-sheet clocks and all-odd tournaments

**PROVED / FINITE-EXACT + INDEPENDENTLY AUDITED.** Full note: [procgen_tcpc_20261001_tournament_clock_prime_collatz](../../05-knowledge/results/procgen_tcpc_20261001_tournament_clock_prime_collatz.md).

## 1. What "Tournament Clock Prime Collatz" became

A *clock* is a circulant-type digraph on a cyclic group or another finite abelian group. Its arcs come from a residue class (squares, cubes, quadratic residues, a discrete-log class), and loops, doubled arcs and missing pairs are allowed. These are the owner's "tournaments with missing, doubled, or self-looped edges on partitions of a modulus".

**Combination laws.**
- **Reflection.** Odd clocks are arc-even.
- **Lexicographic products.** These have an explicit Hamiltonian-path formula.
- **Power-residue clocks mod `3^k`.** These are blown-up cycles.
- **Products.** Clocks multiply by unions of scaled copies.
- **Complements.** Complements are dual, and mod 2 this is Berge's theorem.

## 2. The new existence result

Take two rotational clocks on `Z/r`. Turn the second a quarter and run it backwards, then join them by the quadrant rule. The result `D_r` has every arc on an odd number of Hamiltonian paths, for every odd `r` up to 19.
- **What `D_r` is.** It is Paley minus a vertex with the finite field replaced by the circle.
  - It coincides with Paley minus a vertex for `r = 3, 5`.
  - It equals THM-4529's `QR_127[mu_14]` for `r = 7`.
- **Where it is new.** At `N = 34` and `38`, where no Paley tournament exists.

All-odd tournaments now exist for every `N = 2 (mod 4)` up to 38. The conjecture that `D_r` works for every odd `r` would settle the existence half of HYP-9167.

## 3. Collatz

The Syracuse map is a clock mod 9:
- the step `3A + 1` always lands in the square coset `{1, 4, 7}`;
- the halvings turn the hands by `2^(-v)`;
- `log_2 S(A) = 2(A mod 3) - v (mod 6)`.

This is the mod-9 shadow of the 3-adic structure behind Tao's 2019 theorem. It is an exact reformulation and gives no new bound on Collatz.
