---
id: THM-4858
title: "Exact translations of the Goldbach machines between Turing machines, cellular automata and FRACTRAN. (a) A 53-symbol radius-1 cellular automaton that writes HALT from the blank-tape configuration iff strong Goldbach fails: Leng's 25-state Turing machine, one CA step per TM step. (b) LENGGAME, a 568-fraction FRACTRAN program over 349 primes that halts iff strong Goldbach fails: the same machine in the two-stack encoding, with 4 fractions per state and 9 or 10 per instruction. (c) The square of Goldbach machines: GOLDBACHGAME (37 fractions, THM-4856), LENGGAME (568), the TM (25 states, 49 instructions), the CA (53 symbols). Their sizes measure encodings, not the predicate."
status: >
  PRE-AUDIT (mac-mini-2026-11-01-lenses).
  PROVED (elementary): the step-for-step simulation lemmas (a) and (b), using standard constructions (TM -> CA:
  folklore, e.g. A. R. Smith III 1971; two-stack counter encoding: Minsky 1967; FRACTRAN as a Minsky machine:
  Conway 1987).
  CITED: Leng's Lean theorem that the 25-state machine halts from the blank tape iff some even n >= 4 is not a sum of
  two primes (session65 audit, `goldbach_tm/GoldbachTm/Tm25/Content.lean`).
  FINITE-EXACT:
  - (a) the CA equals the TM at every step for 200,000 steps (pre-registered T1, PASSED);
  - (b) 568 fractions, PASSED (pre-registered T2a); literal integer FRACTRAN equals the dispatch simulator on every
    configuration of the first 300 TM steps (T2c, PASSED).
  FAILED AS STATED: pre-registered T2b and T3. The checker counted returns from the head-bit loop's helper state as
  TM steps. Post-hoc with the alignment fixed (`machine_translations_posthoc.out`):
  - LENGGAME's stacks equal the TM's at each of the first 3,000 TM steps (24,405,662 FRACTRAN steps);
  - compiled BB(2) halts with four 1s.
source: mac-mini-2026-11-01-lenses (owner directive of 2026-11-01: analogies and seamless translation between Turing machines, cellular automata, FRACTRAN and the tournament research)
depends_on:
  - 05-knowledge/results/session65_cryptid_tournament_interfaces_20261026/session65-twentyfive/REPORT.md
related:
  - 01-canon/theorems/THM-4856-goldbachgame-fermatgame-legendregame-halts-iff-counterexample-fractran-programs-built-on-primegames-divisor-loop.md
  - 01-canon/theorems/THM-4852-primegame-phase-law-both-published-lists-agree-selector-cellular-automaton-in-base-29-primorial.md
  - 05-knowledge/results/session135_fractran_tm_structure_20261101/SYNTHESIS.md
  - 05-knowledge/results/lenses_tm_ca_fractran_tournaments_20261101.md
scripts:
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/machine_translations.py
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/machine_translations_posthoc.py
outputs:
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/machine_translations.out
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/machine_translations_posthoc.out
---

# THM-4858 - exact translations of the Goldbach machines

## (a) Turing machine -> cellular automaton

Let M be a Turing machine with states Q, tape symbols Γ, and a partial transition function δ (write, then move).
The radius-1 cellular automaton on Z has alphabet Γ ∪ (Q × Γ) ∪ {HALT}. Its local rule is F(a, b, c):

- If b = (q, g): output HALT if δ(q, g) is undefined, else output the written symbol w.
- Else if a = (q, g) with δ(q, g) = (q', R, w): output (q', b).
- Else if c = (q, g) with δ(q, g) = (q', L, w): output (q', b).
- HALT persists; otherwise output b.

Start the CA from the configuration with one head cell (q_0, 0) and blank tape. Then the CA configuration after t
steps encodes M's configuration after t steps: the head cell carries the state, and the other cells carry the tape.
Proof: a configuration of M has exactly one head, and F changes only the head cell and the neighbour it moves into.
The CA writes HALT at step t + 1 iff M has no transition at its configuration after t steps.

For Leng's machine, |Q| = 25 and Γ = {0, 1}, so the alphabet has 2 + 50 + 1 = 53 symbols. By Leng's theorem,
**strong Goldbach fails iff this 53-symbol radius-1 automaton, started from one head cell (0, 0) on a blank line,
ever writes HALT.**

## (b) Turing machine -> FRACTRAN (LENGGAME)

**Encoding.** For a configuration with head position h:
- R = sum_{i >= 0} tape(h+i) · 2^i, so the head cell is bit 0;
- L = sum_{i >= 1} tape(h-i) · 2^(i-1).

These are stored as the exponents of 3 (R) and 2 (L), with temporaries 5 (A) and 7 (C). The control is one-hot: one
state prime per control state. FRACTRAN cannot express self-loops, so every loop s -> s goes through a helper state
with the single fraction s/s'.

**Fractions per TM state q** (4).
- P_q: [R >= 2: R -= 2, A += 1 -> Ph_q]; [R >= 1: R -= 1 -> B1_q]; [else -> B0_q].
- Ph_q -> P_q.

Afterwards A = R div 2 and the read bit is the index of B.

**Fractions per instruction δ(q, g) = (q', d, w)**, at state Bg_q.
- **d = R** (9 fractions): move A into R; double L through A; move A back into L, adding w; go to P_q'. That is,
  new R = R div 2 and new L = 2L + w.
- **d = L** (10 fractions): split L into its low bit, which is added to R, and L div 2, kept in C; move C into L;
  move A into R four times over, adding 2w; go to P_q'. That is, new L = L div 2 and
  new R = 4(R div 2) + 2w + (L mod 2).
- An undefined δ(q, g) is a state Bg_q with no fractions, so FRACTRAN halts there.

**Step-for-step simulation.** Every reduced fraction keeps exactly one state prime in its numerator and one in its
denominator, because register and state primes are distinct. So from the start integer P_{q_0} (blank tape) the
FRACTRAN run is the Minsky run of this program. Its k-th arrival at a state prime P_q from a non-helper state carries
(L, R) = the TM's two stacks after k steps. The arithmetic of each block is the displayed update of the two stacks.

**Halting.** FRACTRAN halts iff M halts, so LENGGAME halts iff strong Goldbach fails.

**Size.** Leng's 49 instructions are 22 right moves and 27 left moves, so the program has
25·4 + 22·9 + 27·10 = 568 fractions, over 4 register primes and 345 control primes, the largest being 2351.

**Cost.** One TM step costs O(L + R) = O(2^(tape length)) FRACTRAN steps. The first 3,000 TM steps take 24,405,662
FRACTRAN steps.

## (c) The square of Goldbach machines

| machine | size | data representation | cost growth |
|---|---|---|---|
| Leng TM (2024) | 25 states, 2 symbols, 49 instructions | unary blocks on one tape | about n^3.26 steps per even n (n <= 74) |
| its CA (a) | 53 symbols, radius 1 | the same tape, the state in one cell | identical to the TM |
| LENGGAME (b) | 568 fractions, 349 primes | the tape as two binary numbers in exponents | exponential per TM step |
| GOLDBACHGAME (THM-4856) | 37 fractions, 30 primes | counters as exponents (unary) | about 9N^2 per target |

All four halt, or write HALT, iff strong Goldbach fails.

**Translation map.**
- (a) and (b) are exact simulations of the TM.
- GOLDBACHGAME is a different algorithm. It shares only the trial-division design, with divisors tested in the
  opposite order (THM-4856, audit).
- The sizes differ by an order of magnitude. Each measures one encoding of one algorithm, not the Goldbach predicate.
  Session65 showed that reachable padding changes them, and a 21-state machine exists (DeepMind cryptids, 2026).
