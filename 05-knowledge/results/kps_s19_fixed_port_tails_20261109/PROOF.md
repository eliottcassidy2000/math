# THM-5085: fixed-port (O/Z) tails in every base (proofs)

kind-pasteur-2026-11-09-S19.

**Priority (audit M1).** At x = 2, Lemma 1 and Theorem 2 below are KNOWN:
- session195/combined-growth section 2 (the K/D laws);
- Session197 ONE_CORE_BINARY_INVERSE (1), (6) and ROOTED_ORDERED_AMALGAMATION (7)-(8);
- Session198 DIGIT_BLOCK (1), whose pair sum C(A) is Q(I);
- Session205 finite-oz-tail-bank ARITHMETIC_REDUCTION;
- THM-4888.

The dictionary is e = -Delta, d = h_u. The proofs below are those proofs with remainders counted by I(Omega(.), x)
instead of H(.). Only the base-x statement is new.

**Notation.**
- I(Omega(T), x) is the base-x odd-cycle packing count (THM-4940). At x = 2 it is H(T), by OCF.
- Paths are vertex sequences, and one-vertex paths are allowed.
- For a tournament T with a vertex u (the **port**):
  - E(T, u) is the sum of I(Omega(T - V(P)), x) over the paths P of T starting at u with |V(P)| even;
  - O(T, u) is the same sum over paths with |V(P)| odd, including P = (u).
  - R = E + O, Delta = E - O, h_u = I(Omega(T - u), x).

At x = 2, (H, S, R) of SEAM.md is (h, O - E, E + O) (observed by the S18 audit).

## 1. The two letters

Let (T, u) be a ported tournament. Define:
- **S (sink).** Append a vertex s that every existing vertex beats.
- **B (bridge).** Append a vertex w that every existing vertex beats, except that w -> u.

The port stays u for ever. A **fixed-port tail** over a core (T0, u) is T0 followed by a word in {S, B}. Bridges are
the letters I_0 of the pre-registration. A letter I_t (t transitive vertices, then a bridge) is literally the word
S^t B: the t vertices are sinks appended before w. So {S, B} generates the whole pre-registered family.

**Lemma 1 (fixed-port transfer).** Let the state of (T, u) be (h, R, Delta) and let h_u be as above. Then:
- (S): (h, R, Delta) -> (h, 2R, 0);
- (B): (h, R, Delta) -> (h + x (R + Delta)/2, 2R - h_u, -h_u);
- h_u is unchanged by both letters.

*Proof.*
- **h_u.** T' - u is (T - u) followed by one vertex that all of T - u beats. A bridge vertex w beats only u, which is
  deleted. So no new odd cycle appears, and I(Omega(T' - u), x) = I(Omega(T - u), x).
- **S, h.** The new sink lies on no cycle.
- **B, h.** Every new odd cycle uses w -> u, the only arc into u from the new vertex. Such a cycle consists of a path
  P of T from u, then w. It is odd iff |V(P)| is even. All new cycles contain u, so a packing holds at most one. With
  one new cycle, the rest of the packing lies in T - V(P), and w's only other arc points into the deleted vertex u.
  This is THM-5070 Lemma B with B = {w}, giving x Σ_{|P| even} I(Omega(T - V(P)), x) = x E = x (R + Delta)/2.
- **Paths from u in T' (both letters).** Such a path is a path P of T from u, optionally followed by the new vertex
  v_new. The only exception is that for a bridge, (u) cannot be followed by w, since w -> u.
  - The remainder is (T - V(P)), plus v_new if unused. If v_new = w is unused, w beats only the deleted u, so it is a
    sink of the remainder and the count is unchanged.
  - Hence the path generating function f(z) = Σ_P z^|P| I(Omega(T - V(P)), x) becomes:
    - for S: f'(z) = (1 + z) f(z);
    - for B: f'(z) = (1 + z) f(z) - z^2 h_u. The excluded path is (u, w), whose remainder is T - u.
  - Evaluating at z = 1 and z = -1 gives the stated (R, Delta) updates. QED

**Checks.**
- `fp_check.py` builds 800 random tails, with cores of order <= 6, every port, and words of total order <= 15.
- It compares the predicted h with an independent count, in bases 2 and 3: Hamiltonian-cycle DP per odd vertex set,
  then a set-partition DP, plus a direct H DP for OCF.
- Result: 0 mismatches.

## 2. Closed form

**Theorem 2.** Let (T0, u) have state (h0, R0, Delta0) and h_u. Apply a word ℓ_1 ... ℓ_L over {S, B} with bridge set
I = {i : ℓ_i = B}. Put:
- N(I) = Σ_{i in I} 2^(i-1);
- Q(I) = Σ_{j < i, both in I} c(i - j), where c(1) = 2 and c(d) = 2^(d-1) for d >= 2.

Then the final I(Omega, x) is

    h = h0 + (x/2) [ R0 N(I) + Delta0 [1 in I] - h_u Q(I) ],

and the final port state is R = 2^L R0 - h_u Σ_{i in I} 2^(L-i), Delta = -h_u [L in I] (0 if ℓ_L = S, or Delta0 if
L = 0).

*Proof.* By Lemma 1, R_i = 2 R_{i-1} - h_u [i in I], so R_{i-1} = 2^(i-1) R0 - h_u Σ_{j in I, j < i} 2^(i-1-j).
- Delta_{i-1} is Delta0 if i = 1.
- Otherwise Delta_{i-1} is -h_u if i - 1 is in I, and 0 if ℓ_{i-1} = S.
- The bridge at i adds (x/2)(R_{i-1} + Delta_{i-1}).

Summing over i in I, the term -h_u 2^(i-1-j) for j = i - 1 combines with Delta_{i-1} = -h_u into c(1) = 2. QED

**Corollaries.**
- **(a) One-apex words.** T0 = ({u}) has (h0, R0, Delta0, h_u) = (1, 1, -1, 1). Every tail is a one-apex tournament
  with apex u: S is an out-neighbour and B an in-neighbour of the apex. Theorem 2 then reproduces the inherited
  one-apex recurrence (session25; Session191 one_apex/PROOF.md (1); written as Session236 offsets/PROOF.md (3)).
- **(b) Single bridge.** A single bridge after s sinks, I = {s + 1}, gives h = h0 + (x/2) 2^s R0 for s >= 1. This is
  Session236's ray C + 2^m A with B a singleton (Session236 offsets/PROOF.md (2)).
- **(c) Signed offset transport.** Prepending s sinks replaces R0 by 2^s R0 and Delta0 by 0, without changing I's
  relative pattern. So for every bridge pattern I,

      h_s = o(I) + 2^s (x/2) R0 N(I),     o(I) = h0 - (x/2) h_u Q(I),     s >= 1.

  The offset o(I) is invariant along the dyadic ray, as in Lemma 6 of THM-5055. It can be negative, and it can be 7
  or 21.
  - Example: T0 = ({u}) with I = {1, 2} gives o = -1 and N = 3, so h_s = 3 * 2^s - 1. The pre-audit text said
    3 * 2^(s+1) - 1, which was wrong (audit M3).
  - Session236 already exhibits the same values as one-apex negative offsets (offsets/PROOF.md (3)).
- **(d) Bridge rays (the audit's formula).** From a state (h, R, -h_u), the word B^k gives
  [h - (x/2)(R - h_u)] + 2^k (x/2)(R - h_u). Each further bridge doubles R - h_u.

## 3. Pure dyadic rays: sink rays versus bridge rays (corrected after audit M4)

- **Sink rays.** For even x, if the leading coefficient 2^s (x/2) R0 N(I) is a power of 2, then N(I) is one. So
  |I| = 1, Q(I) = 0, and the offset is h0, the value of T0. For odd x the leading coefficient is never a power of 2.
- **Bridge rays.** These differ.
  - **Example.** The order-7 core with rows [94, 84, 32, 6, 12, 27, 60] and port 0 has base-2 state (71, 287, -5, 31).
    After one B the state is (353, 543, -31), and R - h_u = 512 = 2^9.
  - **The ray.** The word B^(k+1) has value exactly 2^(k+8) - 159. This was brute-forced for k = 2..8 at orders 9-15
    (audit `out/dyadic_ray_bridges_bruteforce.out`).
  - **Consequence.** Pure-dyadic rays inside the fixed-port family can carry negative offsets.
- **Withdrawn.** The pre-audit sentences "extends the barrier to multi-bridge tails" and "values just below powers of 2
  need non-dyadic leading parts" were wrong and are withdrawn.
- **Session236's barrier** is scoped to its single-final-bridge mechanism and remains correct there.

## 4. Coverage (finite computations)

- **Program.** `fp_cover.c` takes the closure, under {S, B}, of every distinct core state (h0, R0, Delta0, h_u) over
  all tournament classes of order <= n0 and all ports, from classgen + port_stats.
  - It keeps states with h <= LIM, since h is nondecreasing.
  - A state is a leaf when its cheapest continuation B already overflows: (x/2)(R + Delta) > LIM - h, with Delta <= R.
  - The program can also use I_1, I_2 directly, which are the words S B and S S B.
- **Pre-registered outcomes.** See PREREGISTRATION.md.
