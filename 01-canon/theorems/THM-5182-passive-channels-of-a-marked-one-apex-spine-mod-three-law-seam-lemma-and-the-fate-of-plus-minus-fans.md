---
id: THM-5182
title: "Passive channels of a marked one-apex spine, and the fate of plus/minus fans. (A) Pair expansion: the apex count of a marked ordinal spine is N = A + sum_i (n_i - a_i) A/a_i + sum_{i<j} s_i t_j A_{<i} Q_{ij} A_{>j} + sum_i e_i f_{i+1} A_{<i} A_{>i+1}; blocks compose by (a,q,s,e)(YZ) = (a_Y a_Z, q_Y q_Z, q_Z s_Y + a_Y s_Z, a_Y e_Z) and by the mirror rule (t,f)(YZ) = (q_Y t_Z + a_Z t_Y, a_Z f_Y); the channels (n,t,f) of a block are passive, with weights (h, S, d) of the incoming state. (B) Mod-3 law: for every tournament of order m and every vertex v, T_v = (-1)^(m+1) S_v (mod 3); hence t = (-1)^m s + [m odd] a (mod 3) for every mark set. (C) Seam lemma S_v >= H and T_v >= H; a slide that matches the right-facing data moves t geometrically. (D) FINITE-EXACT: no tournament of order <= 9 has both an armed +1 child and an armed -1 child under one-vertex extensions with arbitrary out-sets. (E) FINITE-EXACT: marked two-block one-apex tournaments with blocks of order <= 7 attain every odd N in [23, 1,542,343), with block orders (8,7) every odd N in [23, 8,667,349)"
status: >
  PROVED (A)-(C) (elementary). FINITE-EXACT (D), (E) and the slide census.
  (A) unrolls the exact cut transfer of chanpak_selective_target_20261114/CUT_STATE.md (KNOWN, credited). The
  mirror rule follows from converse duality. (B) combines the squarefree path identity
  H(X) P_v(-X) = -E_v(X) (chanpak_leaf_fans_20261115/FAN_CHANNELS.md and chanpak_heisenberg_target_20261115,
  credited) with the congruence P_T(2) = (-1)^n H(T) (mod 3) of THM-5172 (B), applied to every induced
  subtournament.
  Checks:
  - composition rules and the per-vertex mod-3 law on 400 random two-block spines (scripts/compcheck.c);
  - the cut transfer itself on 300 random two-block spines against direct DP (scripts/palette.c check);
  - the mod-3 law on every marked tournament of order <= 8 (1,761,280 + smaller; scripts/palette.c run).
  (D): all classes of orders 3-9 as seeds (children to order 10); adaptive trees from orders 3-8 to order 10
  (scripts/armed.c, tree2.c, childpair.c).
  (E): 2,736,649,969 profile pairs at (7,7) and 88,574,069,828 at (8,7) (scripts/twoblock.c); 160 sampled
  witnesses re-counted by direct Hamiltonian-path DP on the assembled tournaments (scripts/twowit.c,
  verify_two.py).
  NOT CLAIMED: Chan-Pak (OPEN); the completeness of the two-block family (HYP-9702); that unit extensions never
  branch twice at higher order (HYP-9712).
session: mac-mini-2026-11-15-slides
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-26-fans-slides-and-reservoirs--unit-extensions-never-branch-twice-a-mod-3-law-for-passive-channels-and-marked-two-block-completeness
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/armed.c, tree2.c, childpair.c, xsize.c, base4.c
  - 05-knowledge/results/chanpak_waves_20261113/scripts/palette.c, compcheck.c, palstats.py
  - 05-knowledge/results/chanpak_waves_20261113/scripts/slides.c, slides2.c
  - 05-knowledge/results/chanpak_waves_20261113/scripts/twoblock.c, twowit.c, verify_two.py, twomult.c
related:
  - THM-5172 (profile congruences), HYP-9692 (band fullness), HYP-9702, HYP-9712,
    chanpak_selective_target_20261114/CUT_STATE.md (the cut transfer),
    chanpak_leaf_fans_20261115 (fan law, cross-port law, movable-port census),
    chanpak_window_descent_20261115 (matched slides, CRT amplifier, exposure theorem)
---

# Passive channels of a marked one-apex spine, and the fate of plus/minus fans

## Setting and notation

A **marked block** is a tournament R of order m together with a mark set V ⊆ V(R): the vertices beaten by a
common apex. Following CUT_STATE (credited), its data are:
- a = H(R) and q = P_R(2), the number of ordered two-slot path covers, empty slots allowed;
- for a nonempty path P, its weight is H(R − P);
- S_v and e_v are the total and the odd-minus-even weights of paths starting at v;
- T_v and f_v are the same for paths ending at v;
- s = Σ_{v∈V} S_v and e = Σ_{v∈V} e_v (starts from marked vertices);
- t = Σ_{v∉V} T_v and f = Σ_{v∉V} f_v (ends at unmarked vertices);
- n = H(R + apex).

A **marked spine** is an ordinal sum B_1 ⇒ B_2 ⇒ ... ⇒ B_r of marked blocks with one common apex. Its count N is
the number of Hamiltonian paths of the whole tournament. The exact cut transfer is

    (h, N, S, d) -> (a h, a N + (n − a) h + t S + f d, q S + h s, h e),

started from (1, 1, 0, 0).

Put A = Π a_i, A_{<i} = Π_{k<i} a_k, A_{>j} = Π_{k>j} a_k and Q_{ij} = Π_{i<k<j} q_k.

## (A) Pair expansion, composition and passivity

**Pair expansion.**

    N = A + Σ_i (n_i − a_i) A/a_i + Σ_{i<j} s_i t_j A_{<i} Q_{ij} A_{>j} + Σ_i e_i f_{i+1} A_{<i} A_{>i+1}.

*Proof.* Unroll the transfer. After block j, S = Σ_{i≤j} s_i A_{<i} Π_{i<k≤j} q_k and d = A_{<j} e_j. Every later
block multiplies N by its a. ∎

Each term is one class of odd apex cycles: a cycle inside block i, or a cycle that leaves the apex into a marked
start in B_i and returns from an unmarked end in B_j. The intermediate blocks contribute q_k each: their vertices
split freely between the cycle's path and the complement. Their signed parts cancel, which is why the signed term
survives only between adjacent blocks.

**Composition.** Write YZ for the spine Y ⇒ Z. Then

    (a, q, s, e)(YZ) = (a_Y a_Z, q_Y q_Z, q_Z s_Y + a_Y s_Z, a_Y e_Z)            (right-facing data)
    (a, q, t, f)(YZ) = (a_Y a_Z, q_Y q_Z, q_Y t_Z + a_Z t_Y, a_Z f_Y)            (left-facing data, mirror rule)

*Proof.* The right-facing rule is the transfer's (S, d) update for a prefix started at (1, ·, 0, 0). For the mirror
rule, reverse every arc, including the apex arcs. This maps the marked spine Y ⇒ Z to Z^op ⇒ Y^op with
complemented marks. It preserves every Hamiltonian count and every apex count. Reversing a path swaps its start
and end, and complementing marks swaps "starts at a vertex the apex beats" with "ends at a vertex that beats the
apex". So converse duality exchanges (s, e) with (t, f). ∎

**Passivity.** A block's right neighbour reads only (a, q, s, e). The channels (n, t, f) enter N linearly, with
coefficients (h, S, d) of the incoming state, times the later ordinary counts. They never enter the outgoing
state. So any family of blocks with equal right-facing data (a palette) can be switched freely. Switching block i
changes N by

    (Δn_i h_{<i} + Δt_i S_{<i} + Δf_i d_{<i}) A_{>i}.

- The n-channel gives the flat weights A/a_i. These are the CRT amplifier of chanpak_window_descent_20261115.
- The t-channel gives positional weights S_{<i} A_{>i}. They grow by about q_i/a_{i+1} per block.
- The last block of a spine is passive in all of (a, n, t, f).

## (B) The mod-3 law

**Theorem.** For every tournament R of order m and every vertex v,

    T_v ≡ (−1)^(m+1) S_v (mod 3).

Consequently, for every mark set,

    t ≡ (−1)^m s + [m odd]·a (mod 3),

and dually s ≡ (−1)^m t + [m odd]·a.

*Proof.* Work in the squarefree algebra on the vertices (x_v² = 0). Let H(X) = Σ_U H(R[U]) x^U, and let P_v(X) and
E_v(X) be the generating series of paths starting, respectively ending, at v. The identity H(X)P_v(−X) = −E_v(X)
holds (credited). Then

    S_v = [x^V] H(X)P_v(X),     T_v = [x^V] H(X)E_v(X) = −[x^V] H(X)² P_v(−X).

Now [x^U] H(X)² = Σ_{A⊆U} H(A)H(U∖A) = P_{R[U]}(2). By the odd-cycle form P(q) = q^k + Σ_j w_j q^(k−2j) (THM-4900,
THM-5172 (B)),

    P_{R[U]}(2) ≡ P_{R[U]}(−1) = (−1)^|U| H(U) = [x^U] H(−X)   (mod 3).

So H(X)² ≡ H(−X) coefficientwise mod 3, and

    T_v ≡ −[x^V] H(−X)P_v(−X) = −(−1)^m [x^V] H(X)P_v(X) = (−1)^(m+1) S_v.

For the consequence, sum over the unmarked vertices. Use Σ_v S_v = q − a and q ≡ (−1)^m a (mod 3). ∎

**Consequences.**
1. In any palette with fixed (a, q, s, e), all the t-values lie in one residue class mod 3. In the census of
   orders ≤ 8 the t-values of a fiber are never consecutive; their gcd of differences is 3.
2. The law is mark-blind. It constrains every passive t-reading by the right-facing data alone.

## (C) Seam lemma and geometric t-slides

**Seam lemma.** For every vertex v, S_v ≥ H(R) and T_v ≥ H(R).

*Proof.* Cut every Hamiltonian path at v. The suffix from v is a path starting at v, and the prefix is a
Hamiltonian path of its complement. This injects the Hamiltonian paths into the pairs counted by S_v. For T_v,
cut after v instead. ∎

Hence in a two-block layout L ⇒ B the seam term satisfies s_L t_B ≥ |V_L| · H(L) · t_B. If B has an unmarked
vertex preceded by a transitive run of ℓ vertices, then t_B ≥ 2^ℓ H(B).

**t-slides are geometric.** Suppose a core C and transitive words z, w of width d satisfy Cz ≡ wC in the
right-facing data. This is the first matching condition of the matched slides. Place the slide after the prefix
w^j. Since q(w^j) = 2^(dj), the mirror rule gives

    t(w^j C z^(k−j)) − t(w^(j+1) C z^(k−j−1)) = 2^(dj) (t(Cz) − t(wC)).

So any t-palette built by sliding is a geometric progression, never an arithmetic one. The matched slides of
chanpak_window_descent_20261115 avoid this by also matching t. They then move only the n- and signed channels, at
flat weights.

## (D) Plus/minus fans with arbitrary out-sets (FINITE-EXACT)

A **unit extension** of T adds one vertex x with an arbitrary out-set X ⊆ V(T) and has H(T + x) = 2H(T) ± 1. By
the odd-cycle expansion, H(T + x) = H(T) + 2W(X), where

    W(X) = Σ_{u∈X, v∉X} Π(u,v),   Π(u,v) = Σ over even-order paths u → v of H(T − path).

So a unit extension is a cut of value (H ± 1)/2. T is **armed** if it has unit extensions of both signs. The leaf
model of chanpak_leaf_fans_20261115 is the special case |X| = 1. Here every out-set is allowed, and about half of
all unit out-sets are not leaves or reverse leaves.

| seed order | classes | eligible (H ≠ 3, 11) | armed | armed fraction |
|---:|---:|---:|---:|---:|
| 5 | 12 | 8 | 5 | 0.625 |
| 6 | 56 | 50 | 18 | 0.360 |
| 7 | 456 | 448 | 50 | 0.112 |
| 8 | 6,880 | — | 134 | 0.019 of all classes |
| 9 | 191,536 | — | 308 | 0.0016 of all classes |

For all seeds of orders 3–8, an exhaustive search allowed the out-set to be re-chosen at every node:
- **No tournament of order ≤ 9 has both an armed +1 child and an armed −1 child.** At order 8, 8 armed seeds
  have an armed +1 child and 18 an armed −1 child; at order 9, 10 and 27 do. None has both.
- Adaptive binary trees of unit extensions have depth at most 1 from every seed with H ∉ {3, 5, 11}, which is all
  6,868 order-8 classes with H ≥ 9 except the H = 11 ones.
- When the branches whose target is 7 or 21 are skipped, depth 2 occurs only at H ∈ {5, 11} and depth 3 only at
  H = 3. From the 3-cycle the reachable values are 3 → 5 → {9, 11} → {17, 19, 23}.
- The children show **sign freezing**. A +1 child usually has only +1 unit extensions, and a −1 child only −1
  ones. Branching survives on one side at most.

Two-vertex steps (H → 4H + ε, ε ∈ {±1, ±3}, any two out-sets) from seeds of orders 3–6 also stop at depth ≤ 1,
except H = 3 (depth 2).

## (E) The marked two-block family (FINITE-EXACT)

Take L ⇒ B with one apex. L is the reservoir: its right-facing data are (g, n_L, s_L, d_L) = (a, n, s, e) of L.
B is the repair: its data are (a_B, n_B, t_B, f_B). By (A),

    N = a_B (n_L − g) + g n_B + s_L t_B + d_L f_B,

a bilinear pairing of the reservoir's right-profile with the repair's left-profile. Over all marked tournaments of
the given orders:

| block orders | profiles per side | every odd N attained in |
|---|---:|---|
| (≤ 6, ≤ 6) | 2,831 | [23, 67,967) |
| (≤ 7, ≤ 7) | 52,313 | [23, 1,542,343) |
| (≤ 8, ≤ 7) | 1,693,156 × 52,313 | [23, 8,667,349) |

Together with 1, 3, ..., 19, this means every odd N < 8,667,349 except 7 and 21 is the count of a tournament of
order ≤ 16 that has a vertex whose removal leaves an ordinal sum of two blocks. (8,7) and (7,8) give the same
spectrum, by converse duality.

Two further measurements:
- **Palettes are trivial.** For a last block, a palette is a fiber over (a, t, f). These fibers have mean size
  1.05–1.09 at orders 5–8, and their longest run of consecutive odd apex counts is 4. So the coverage comes from
  the variety of interfaces, not from repair windows.
- **Multiplicities.** On odd N in [10^6, 1.4·10^6] the (7,7) representation count averages 215.9 for N ≡ 0
  (mod 3) and 171 for N ≡ ±1 (mod 3). Primes average 171.1, exactly their residue classes. Residues mod 5 and 7
  move the average by less than 2%. The minimum is 16 representations.

## Slide census (FINITE-EXACT)

Matched slides need Cz ≡ wC in (a, q, s, t). Every strong marked core of orders 3–7 was searched, with transitive
tiles z, w, the guards of chanpak_window_descent_20261115, and the step Δ = N(pCzq) − N(pwCq):
- **Orders 3–4:** no slide with Δ ≠ 0.
- **Order 5, width ≤ 14:** only the count-13 core has slides with Δ ≠ 0. All 48 have width 12, with
  Δ ∈ {±10, ±20, ±30}, and each step was confirmed by direct DP.
- **Order 6, width ≤ 30:** no slide with Δ ≠ 0.
- **Order 7, width ≤ 28:** 16 slides with Δ ≠ 0, from three cores:
  - two count-29 classes (Q = 580) at width 28, with Δ ∈ {±6, ±32};
  - a count-115 core (Q = 1460) at width 24, with Δ = ±46.

  These steps come from the verified transfer, since direct DP was feasible only up to 20 vertices (width ≤ 13
  here). The Δ = 6 slide of a count-29 core agrees with the literal replay of chanpak_window_descent_20261115.

The half-step 23 of the count-115 core divides 115, so that slide alone meets only 5 of the 115 residues.

No unit pump (|Δ| = 2) was found.
