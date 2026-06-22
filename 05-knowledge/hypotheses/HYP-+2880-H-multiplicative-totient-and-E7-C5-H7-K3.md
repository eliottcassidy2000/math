---
id: HYP-+2880
title: H multiplicative over strong-component ATOMS (the totient/multiplicative-function analogy; {7,21} the finite forbidden gaps, NOT a 7-ramification -- 35=5*7,49=7^2 ARE atoms); + the E_7 C_5-hole <-> H=7 K_3-obstruction question (apex-7 odd-cycle, perfect-graph odd-hole<->odd-clique)
status: Thread B VERIFIED (atom H-values by n; {7,21} forbidden gaps; ramified-prime instinct CORRECTED). Thread A grounded (H=7=K_3 of 3 conflicting cycles, THM-029; E_7 has C_5 holes, kps S31e); bijection map is the open literal-test (shared with kps).
source: mac-mini-2026-06-22-S37 (user: E_7 C_5<->H=7 K_3 bijection + H-multiplicative<->totient)
related:
  - HYP-2873   # forbidden-H finite {7,21}, achievable cofinite
  - h21-moon-reduction-s617  # H = prod H(strong components)
  - THM-029    # H=7 impossible; THM-079: H=21 impossible
---

# HYP-+2880 -- H multiplicative (totient analogy) + the E_7 C_5 <-> H=7 K_3 question

## Thread B: H multiplicative over atoms; the totient analogy (VERIFIED)
H(T) = prod H(strong components) (Moon). The achievable-H = the multiplicative semigroup generated
by the ATOM (irreducible single-strong-component) H-values. ATOM H-values by n (verified): n=3:{3},
n=4:{5}, n=5:{9,11,13,15}, n=6:{15,17,...,45}, n=7:{25,...,35,...,49,51}. 
- 7 is NEVER an atom (H=7 forbidden, THM-029); 21 never (THM-079). So {7,21} = the ONLY forbidden
  H-values = the finite Kuratowski gaps (HYP-2873).
- **CORRECTION (discipline):** 7 is NOT a 'ramified prime' -- 35=5*7 and 49=7^2 ARE atoms (n=7), so
  7 DOES divide atoms. {7,21} are forbidden simply as non-atom-non-product GAPS, not a 7-power rule.
- **TOTIENT ANALOGY:** H is a multiplicative function (atoms = its 'irreducibles'), like phi(n)=
  n*prod(1-1/p) (primes = irreducibles). The previous euler-totient work (LRC Mertens Sum phi(b)/b ~
  6q/pi^2, ze(2) floor HYP-2856) is the SAME multiplicative-function machinery; the apex prime 7 is
  special in BOTH (H: {7,21} gaps; LRC: D=14=2*7 obstruction). Both governed by multiplicativity +
  the apex-7. (Open: a sharp Dirichlet-series / Euler-product for the H-distribution via the atoms.)

## Thread A: E_7 C_5 holes <-> H=7 K_3 obstruction (the literal pentagon test, for kps)
- **H=7 obstruction = K_3 (PROVED, THM-029):** H=7 requires 3 independent odd cycles all pairwise
  CONFLICTING = a triangle (K_3) in the conflict graph -- which cannot exist. The K_3 (3-cycle) of
  conflicting cycles is the obstruction.
- **E_7 has C_5 holes (kps S31e):** the even-graph metagraph at n=7 loses chordality, gaining
  chordless 5-cycles (C_5 odd holes) -- the apex-prime non-chordality.
- **THE QUESTION (literal match):** do E_7's C_5 holes map, via the even-graph<->tournament cycle-space
  bijection (cycle_space_bijection_s20ge), onto the H=7 K_3? The C_5 (odd hole) vs K_3 (odd clique/3-
  cycle) is the PERFECT-GRAPH odd-hole<->odd-clique duality (C_5 self-complementary; both the minimal
  'odd' obstruction at apex 7). If the bijection sends C_5 -> K_3, the thematic pentagon match becomes
  LITERAL: one apex-7 odd-cycle phenomenon underlying {7,21}, E_7 non-chordality, AND LRC(14). 
  @kps (you have E_7 + 7-adic atoms) -- this is the test to run with your E_7 hole data + the bijection.
-> HYP-2873, h21-s617, THM-029, kps E_7 (S31e), kps 7-adic atoms (S31f).


## Thread A VERIFIED (S37): the C_5 <-> H=7 K_3 match IS literal via the cycle space
COMPUTED in K_5/K_7 cycle space (GF(2)):
- **C_5 = XOR of 3 pairwise-VERTEX-CONFLICTING triangles** (the fan {123}+{134}+{145} = the 5-cycle
  1-2-3-4-5; all share vertex 1). ROBUST: ALL 5 triangle-triples in K_5 XORing to C_5 are pairwise
  vertex-conflicting (5/5).
- The **H=7 obstruction = 3 pairwise-vertex-conflicting cycles = a K_3** in the conflict graph (THM-029).
- The bijection is even-graph = XOR of fundamental cycles (cycle_space_bijection_s20ge), so a C_5
  even-graph vertex of E_7 corresponds, in the tournament cycle space, to exactly the 3-conflicting-
  triangle (H=7 K_3) configuration. **=> the thematic pentagon match is LITERAL:** C_5 (even-graph /
  E_7 side) = the H=7 K_3 (tournament / H side) = the apex-7 odd-cycle obstruction. One phenomenon
  underlies {7,21}, E_7 non-chordality, and (via the apex-7) LRC(14).
- **HONEST distinction:** this is the C_5 even GRAPH (a 5-cycle, one vertex of E_7). Whether kps's E_7
  METAGRAPH HOLES (chordless 5-cycles in the iso-class ADJACENCY, S31e) are these same C_5 (or reduce
  to them) is the connecting step -- @kps to check with your E_7 hole data + the bijection. The
  cycle-space identity (C_5 = 3 vertex-conflicting triangles = H=7 K_3) is the literal core, verified.