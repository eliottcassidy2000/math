---
id: THM-4852
title: "PRIMEGAME over the first ten primes: the owner's list is Conway's 1987 PRIMEGAME (Theorem 1), Guy's 1983 list is the earlier 15/14 form; at each entry into the state-free set both pass through 2^n 7^(D(n)-1) (Conway 1987, section 5, Figure 2) and reset to 3^n 5^(n+1) 11, so they emit the same primes; the owner's pure fractions 1/7, 1/17, 55/1 = {7, 17, 55}; every FRACTRAN program in lowest terms over primes <= P is a selector (quasi-)cellular automaton on base-P# digits with the rule chosen by the units digit (new to the repository only in this generality)"
status: CITED (Conway 1987 Theorem 1 and section 5 Figure 2 - the phase law and the roles of 1/7, 55/1, 1/17; Guy 1983; OEIS A007546/A007547; Kari 2012, Kari-Kopra 2017, Cloney-Goles-Vichniac 1987; repository collatz_oscillation_gilbreath_20260926 section 3) + PROVED (elementary: the reduction of the two lists to the law; the general selector-automaton statement) + FINITE-EXACT; metagraph comparison a REMARK, readings NUMEROLOGY; independently audited SOUND WITH CORRECTIONS (MISTAKE-697)
session: opus-2026-10-31-S70
date: 2026-10-31
note: 05-knowledge/results/primegame_variant_selector_automaton_metagraph_20261031.md
script: 04-computation/experiments/primegame_variant_phase_law_20261031.c (+ .out), primegame_selector_automaton_metagraph_20261031.py (+ .out)
audit: 04-computation/experiments/primegame_selector_automaton_20261031_audit_S70/AUDIT.md
related: [HYP-9448, THM-4729, THM-4730, THM-4838, THM-4759]
external:
  - "J. H. Conway, FRACTRAN: a simple universal programming language for arithmetic (1987): Theorem 1 (p. 4, the owner's list), section 5 Figure 2 (pp. 9-10, the proof), section 11 (residue-class structure, modulus)"
  - "R. K. Guy, Conway's prime producing machine, Math. Mag. 56 (1983) (the 15/14 list); OEIS A007546 (fast version; A. P. Heinz 2009 step formula), A007547 (slow version)"
  - "J. Kari (DLT 2012) and J. Kari, J. Kopra (RAIRO-ITA 51(4), 2017): cellular automata multiplying by p/q in base pq; Cloney, Goles, Vichniac (1987): quasi cellular automata"
---

# THM-4852 - PRIMEGAME over the first ten primes

**Setting.**
- States are exponent vectors over 2..29: data registers 2, 3, 5, 7 and state registers 11, 13, 17, 19, 23, 29.
- The owner's list is Conway's 1987 PRIMEGAME (Theorem 1), ending 15/2, 1/7, 55/1.
- Guy's 1983 list ends 15/14, 15/2, 55/1. It is the "fast version" (OEIS A007546); Conway's is the "slow version"
  (A007547).

**(1) Phase law and the two lists** (CITED + elementary). Each entry into the state-free set (no state prime present)
is 2^n 7^(D(n)-1) for n = 1, 2, 3, ..., with D(n) = n/spf(n). This is Conway's proof of Theorem 1, section 5,
Figure 2, which also gives the reset word (15/2)^n (1/7)^(D(n)-1) (55/1).

The two lists differ only in their last three rules:
- From the entry 2^a 7^b (b = D(a)-1 < a), Guy's list resets by (15/14)^b (15/2)^(a-b) (55/1) and the owner's by
  (15/2)^a (1/7)^b (55/1). Both reach 3^a 5^(a+1) 11.
- Conversely, in any configuration that contains a state prime, the first applicable rule is among A-K, the rules the
  two lists share. The denominators of F-K are exactly 29, 23, 19, 17, 13, 11, and every state-bearing configuration
  has exactly one state prime.

So the orbits coincide at every phase boundary for all n. As computational confirmation, both lists were re-run to
n = 2500 (HYP-9448), emitting exactly the 367 primes below 2500.

**(2) Anatomy over 2..29.**
- The lcm of the denominators is 29# (Conway section 11; repository procgen_atlas .out).
- Every part is 1, a prime, a semiprime, or 78 = 2*3*13. Under the repository's divisor-count convention F = U + S
  (`arithmetic_braids_20260917_divisors.md`; `aliquot_66_five_fifteen_20261028/PROOFS.md`) the balanced parts are
  exactly the primes.
- The owner's list has exactly three pure fractions, 1/7, 1/17 and 55/1. Their roles follow Conway's semantics:
  - 1/7 erases the counter 7^(D(n)-1);
  - 55/1 starts the next division;
  - 1/17 leaves state 17.
- The product is 55^2/(2*17^2) for the owner's list; Guy's is 15/2 times that.

**(3) Selector automaton.** Let B = P# and let f be a FRACTRAN program whose fractions are in lowest terms with all
primes <= P.
- For each fraction p/q, the map x -> x p/q on the multiples of q is digit-local on base-B digits:
  - divide by q, reading the right neighbour;
  - multiply by p (p | B^k), reading the left neighbour, with bounded carry.
- The fraction is chosen by x mod q, i.e. by the units digit.
- So f is a **selector (quasi-)cellular automaton**: each step applies one of finitely many local rules to all digits,
  the rule being selected by the units digit. It does not commute with the shift, so it is not a cellular automaton in
  the strict sense: F(29# * 2) = 2417248020, but 29# * F(2) = 97045398450.
- Prior art: Kari 2012 (3x/2 in base 6), Kari-Kopra 2017 (p/q in base pq), Cloney-Goles-Vichniac 1987 (quasi cellular
  automata), Conway 1987 section 11 (residue selection), and the repository's Collatz version
  (`collatz_oscillation_gilbreath_20260926.md` section 3).
- New to the repository is only this general FRACTRAN statement. It was checked on PRIMEGAME (20,000 steps of each
  list) and on 400 random programs by the audit.

**(4) Remark (the n = 5 merged metagraph).** The register-flow graph (an arc p -> q when p divides a denominator and q
the numerator of the same fraction: 24 arcs, 19 pairs) and the merged tiling line metagraph of THM-4729 (18 pairs, two
black loops) differ in size, so no adjacency-preserving bijection exists.
- That comparison is nearly vacuous. For the prime co-occurrence graph there are 70 adjacency-preserving bijections
  (audit), and none sends the owner's {7, 17, 55} registers to the pure-blue nodes TT5, 1->C3->1, T5.
- The count matches (10 = 4 + 6; 3 pure fractions = 3 pure-blue nodes) are NUMEROLOGY.

**Novelty.** Almost everything here is classical or prior repository work. New to the repository are:
- the general selector-automaton statement for all FRACTRAN programs;
- the explicit two-list reduction;
- the anatomy observations, the roles of the owner's three pure fractions, and the remark.
