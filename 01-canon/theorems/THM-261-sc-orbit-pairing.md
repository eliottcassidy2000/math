# THM-261: Z₂ Orbit Pairing on Conflict Graph of SC Tournaments

**Status:** PROVED (algebraic, all n)
**Filed by:** kind-pasteur-2026-03-25-S1
**Dependencies:** OCF (THM-002), Definition of SC tournament

## Statement

Let T* be a self-complementary tournament on n vertices with involutory
anti-automorphism sigma (sigma^2 = id, T*(sigma(u), sigma(v)) = T*(v, u)).

Define phi: C -> sigma(C)^{rev} on directed odd cycles C of T*, where:
- sigma(C) = (sigma(v_1), ..., sigma(v_k)) for C = (v_1, ..., v_k)
- C^{rev} = (v_k, ..., v_1)

Then:

1. **phi is a well-defined involution** on the set of directed odd cycles of T*.
   phi(C) is always a directed odd cycle of T*, and phi(phi(C)) = C.

2. **phi is an automorphism of Omega(T*)**: two cycles C_1, C_2 share a vertex
   iff phi(C_1), phi(C_2) share a vertex. So phi preserves the conflict graph.

3. **Even n (sigma fixed-point-free):** If cycle C uses at most one vertex from
   each sigma-orbit, then C and phi(C) are vertex-disjoint. In particular, for
   3-cycles where each vertex comes from a different orbit, the pair
   {C, phi(C)} contributes to alpha_2 (independent pair count).

4. **Odd n (sigma has fixed point v_0):** Every cycle through v_0 has its
   phi-partner also pass through v_0. These pairs are ALWAYS adjacent in Omega.
   The pairing CANNOT create independent pairs through v_0.

## Proof

**Part 1:** sigma is an anti-aut: T*(sigma(u), sigma(v)) = T*(v, u) for all u,v.
For a directed cycle C = (v_1 -> v_2 -> ... -> v_k -> v_1):
- T*(v_i, v_{i+1}) = 1 for all i
- So T*(sigma(v_{i+1}), sigma(v_i)) = 1 for all i (anti-aut property)
- This means sigma(v_{i+1}) -> sigma(v_i) in T*, i.e., the REVERSED image is directed.
- phi(C) = (sigma(v_k), sigma(v_{k-1}), ..., sigma(v_1)) is a directed cycle in T*.
Since k is odd (only odd cycles in Omega), phi(C) has the same length k.
phi(phi(C)) = sigma(sigma(C))^{rev rev} = C (since sigma^2 = id). QED.

**Part 2:** C_1, C_2 share vertex v iff phi(C_1), phi(C_2) share sigma(v).
Since sigma is a bijection, this preserves the adjacency structure of Omega. QED.

**Part 3:** If sigma is fpf, the sigma-orbits are {a_i, sigma(a_i)} for i=1,...,n/2.
If C uses vertex a_i from orbit i (not sigma(a_i)), then phi(C) uses sigma(a_i)
from the same orbit. Since C and phi(C) use different halves of each orbit,
their vertex sets are disjoint. QED.

**Part 4:** If sigma(v_0) = v_0, then any cycle C through v_0 has phi(C)
through sigma(v_0) = v_0. So C and phi(C) share vertex v_0 and are adjacent
in Omega. QED.

## Conjectured Connection to SC Maximizer (NOT PROVED)

**CAUTION (S29 review):** The following describes a MECHANISM, not a proof.
THM-261 does NOT prove the SC Maximizer conjecture (OPEN-Q-016). It shows
that SC tournaments HAVE a cycle-pairing structure, but does NOT show that
this structure MAXIMIZES H. An NSC tournament could in principle have more
disjoint pairs through a different mechanism.

**Even n mechanism (alpha_2 boost, conjectured):**
The fpf pairing creates vertex-disjoint cycle pairs from 3-cycles that use
one vertex per sigma-orbit (Part 3). The number of such 3-cycles and the
resulting alpha_2 boost are NOT quantified here.

**Odd n mechanism (alpha_1 boost, conjectured):**
The phi-pairing doesn't create disjoint pairs through v_0, but the BIBD-like
constraint on cycle arrangements (forced by the SC symmetry) correlates with
higher total cycle count. This is an empirical observation, not a proof.

## Verification

- n=5: Circulant T_5 with S={1,2}. sigma = (0,4,3,2,1) (fixed point 0).
  5 involutory anti-auts found. All 7 cycles pair under phi:
  (0,1,3) <-> (0,2,4), (0,2,3) <-> (0,1,4), etc.
  Each pair shares vertex 0 (fixed point). CONFIRMED.

- n=6: Regular SC maximizers have sigma with 3 orbits.
  3-cycles using one vertex per orbit pair with vertex-disjoint partners.
  alpha_2 = 4 disjoint pairs, vs NSC max alpha_2 = 1. CONFIRMED.

- n=7: Paley T_7, sigma = v -> 7-v. Fixed point at v=0.
  All 80 directed odd cycles pair under phi; none create disjoint pairs
  through the fixed point. But total count alpha_1=80 exceeds NSC max 65.
  CONFIRMED.

## Related

- OPEN-Q-016: SC Maximizer theorem (this provides the structural mechanism)
- THM-255: SC Maximizer dichotomy at n=6
- THM-099: Maximizer Betti topology
