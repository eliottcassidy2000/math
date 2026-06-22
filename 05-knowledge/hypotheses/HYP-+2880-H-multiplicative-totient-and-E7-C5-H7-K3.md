---
id: HYP-+2880
title: H multiplicative over strong-component ATOMS (the totient/multiplicative-function analogy; {7,21} the finite forbidden gaps, NOT a 7-ramification -- 35=5*7,49=7^2 ARE atoms); + the E_7 C_5-hole <-> H=7 K_3-obstruction question (apex-7 odd-cycle, perfect-graph odd-hole<->odd-clique)
status: Thread B VERIFIED (atom H-values by n; {7,21} forbidden gaps; ramified-prime instinct CORRECTED). Thread A RESOLVED WITH A TWO-LEVEL DISTINCTION: S37 verifies directed C5 support = three vertex-conflicting triangles = H=7 K3 in cycle space; S100/HYP-2881 shows an E7 metagraph C5 hole is not that single support object but a five-class quotient cycle.
source: mac-mini-2026-06-22-S37 (user: E_7 C_5<->H=7 K_3 bijection + H-multiplicative<->totient)
related:
  - HYP-2873   # forbidden-H finite {7,21}, achievable cofinite
  - HYP-2881   # S100 fixed-path audit of E7 C5 vs H=7 K3
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

## Thread A: E_7 C_5 holes <-> H=7 K_3 obstruction (literal pentagon test resolved)
- **H=7 obstruction = K_3 (PROVED, THM-029):** H=7 requires 3 independent odd cycles all pairwise
  CONFLICTING = a triangle (K_3) in the conflict graph -- which cannot exist. The K_3 (3-cycle) of
  conflicting cycles is the obstruction.
- **E_7 has C_5 holes (kps S31e):** the even-graph metagraph at n=7 loses chordality, gaining
  chordless 5-cycles (C_5 odd holes) -- the apex-prime non-chordality.
- **TWO-LEVEL ANSWER:** the directed pentagon support really is the H=7/K3 cycle-space object,
  but an E7 metagraph C5 hole is one quotient level higher.

## Thread A VERIFIED (S37): the directed C_5 <-> H=7 K_3 support match is literal
COMPUTED in K_5/K_7 cycle space (GF(2)):
- **C_5 = XOR of 3 pairwise-VERTEX-CONFLICTING triangles** (the fan {123}+{134}+{145} = the 5-cycle
  1-2-3-4-5; all share vertex 1). ROBUST: ALL 5 triangle-triples in K_5 XORing to C_5 are pairwise
  vertex-conflicting (5/5).
- The **H=7 obstruction = 3 pairwise-vertex-conflicting cycles = a K_3** in the conflict graph (THM-029).
- The bijection is even-graph = XOR of fundamental cycles (cycle_space_bijection_s20ge), so a C_5
  even-graph vertex of E_7 corresponds, in the tournament cycle space, to exactly the 3-conflicting-
  triangle (H=7 K_3) configuration. **=> the directed support match is literal:** C_5 as a cycle-space
  vector equals the H=7 K3 obstruction support.
- **S100/HYP-2881 quotient audit:** HYP-2881 then tests the E7 metagraph-hole reading with the same
  fixed-path fundamental-cycle map: fix
  `0->1->...->6`; a reversed free arc `(i,j)` maps to the path-fundamental even cycle
  `i,i+1,...,j,i`.  Under this map, a directed pentagon support is one E7 vertex class
  (class 3), whereas an E7 C5 hole is a five-class quotient cycle.  The H=7 point
  `alpha=(3,0)` has zero masks/classes; `k3_forces_pentagon` occupies five classes and
  hits `835/1496` E7 C5 holes, but no E7 C5 hole is entirely made from those classes.
  Therefore "C5 = K3" is literal for the directed cycle-space support, while "E7 metagraph C5 hole
  = H=7 K3" is false as object equality.  The next test is incidence-profile comparison against
  the E7 C5/C7 hole classes, not rerunning either literal support or metagraph-hole question.
-> HYP-2873, HYP-2881, h21-s617, THM-029, kps E_7 (S31e), kps 7-adic atoms (S31f).
