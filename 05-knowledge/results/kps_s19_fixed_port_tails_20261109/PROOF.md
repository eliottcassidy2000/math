# THM-5085: multi-bridge fixed-port tails (proofs)

kind-pasteur-2026-11-09-S19.

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
  - Example: T0 = ({u}) with I = {1, 2} gives o = -1, so h_s = 3 * 2^(s+1) - 1.
  - Session236's barrier (offsets are positive H-values) concerns rays whose leading coefficient A is a single
    product R(q_e + q_o). Here the leading coefficient (x/2) R0 N(I) has N(I) with at least two bits whenever
    Q(I) > 0.

## 3. Pure-dyadic rays keep Session236's barrier

If the leading coefficient (x/2) R0 N(I) is a power of 2, then N(I) is a power of 2. So |I| = 1 and Q(I) = 0, and the
offset is h0, an H-value of T0 (at x = 2). The family therefore realises 2^K + c along a pure dyadic ray only for
offsets c that are themselves values. Values just below powers of 2 need non-dyadic leading parts. This agrees with
Session236's barrier and extends it to multi-bridge tails.

## 4. Coverage (finite computations)

- **Program.** `fp_cover.c` takes the closure, under {S, B}, of every distinct core state (h0, R0, Delta0, h_u) over
  all tournament classes of order <= n0 and all ports, from classgen + port_stats.
  - It keeps states with h <= LIM, since h is nondecreasing.
  - A state is a leaf when its cheapest continuation B already overflows: (x/2)(R + Delta) > LIM - h, with Delta <= R.
  - The program can also use I_1, I_2 directly, which are the words S B and S S B.
- **Pre-registered outcomes.** See PREREGISTRATION.md.
