# What the carrier checksum is: the odd-cycle formula of a one-apex tape

Session klein-2026-11-15-ocf. **Chan–Pak remains OPEN.**

## The unifying frame (KNOWN)

The repo's foundational THM-002 (OCF) states H(T) = I(Ω(T), 2). Here Ω(T) is the conflict graph whose vertices are the
directed odd cycles of T, adjacent when they share a vertex. This is the Grinberg–Stanley (Rédei–Berge) expansion.

Its path-cover form is

    P_T(q) = Σ_K 2^|K| q^(n − |V(K)| + |K|),

summed over packings K of vertex-disjoint directed odd cycles. This is THM-4900's odd-cycle form. It was re-checked
here on 80 random tournaments of order ≤ 7 at q = 1, 2, 3, with 0 mismatches (`oddcycle.py`).

## 1. The apex identity (NEW; verified exhaustively for tapes of length ≤ 9)

Let T be the transitive tournament on positions 0 < 1 < … < m−1 (i → j for i < j), plus one apex v. Let u be the
m-bit word with u_i = 1 iff position i beats v. Then

    H(T) = 1 + p(u),       p(u) = 2 · #{directed odd cycles of T}.

Here p is the carrier checksum p(0) = 0, p(2w+1) = p(w), p(2w) = p(w) + w + (w mod 2).

*Proof.*
- T − v is acyclic, so every directed cycle passes through v and Ω(T) is a complete graph. Hence
  H = I(K_c, 2) = 1 + 2c, where c is the number of odd cycles.
- A cycle through v has the form v → a₁ → ⋯ → a_k → v with a₁ < ⋯ < a_k. Here v beats a₁, and a_k beats v. The cycle
  is odd iff k is even.
- For fixed ends i = a₁ < j = a_k, the number of choices is the number of even subsets of the j − i − 1 positions in
  between: 2^(j−i−2) if j > i + 1, and 1 if j = i + 1.
- So c = Σ_(i<j, u_i=0, u_j=1, j≥i+2) 2^(j−i−2) + #{i : u_i = 0, u_(i+1) = 1}.
- The checksum's own expansion is p(u) = Σ_(i<j, u_i=0, u_j=1) 2^(j−i−1) + #{i : u_i = 0, u_(i+1) = 1}.
- Separating the adjacent pairs in p gives p = 2c exactly. ∎

All 1,022 tapes with m ≤ 9 agree with brute force (`oddcycle.py`). For an interval out-set it reproduces
H = 2^m − 2^t − 2^(m−t) + 3.

**What it explains.**
- *Parity.* p is always even: this is Rédei's theorem for the apex tournament.
- *THM-5232's Lemma 1, p(rc u) = p(u), is the converse symmetry H(T) = H(T^op).* Reversing every arc reverses the tape
  order and complements the set of positions beating v. That turns u into rc(u) = 2^m − 1 − rev_m(u). So the
  reverse-complement collision family is a family of *partial converses* of a tape block.
- *The split identity* p(Au + b) = p(u) + p(b) + (A − 1 − rev(b))u + … sorts odd cycles by whether they stay in the
  high block, stay in the low block, or cross the cut.
- *Carrier rows.* The surcharge a(p(X) + p(Y)) in every carrier form (THM-5192 (A)) is therefore an odd-cycle count:
  the cycles that run through the apex inside one tape. The cut-transfer coordinate p(X) + 1 of that proof is H of the
  X-tape with its apex.

## 2. How far one apex reaches

A tournament whose odd cycles all share a vertex v has T − v free of odd cycles. A tournament with no odd cycle is
transitive, so T is exactly an apex over a transitive tape. Complete conflict graphs from a common vertex therefore
give exactly the values {1 + p(u)}.

Over tapes of length ≤ 24 (`apexspec.py`), this family misses a stable half of the odd numbers:

| range | odd N missed |
|---|---:|
| N < 2^10 | 44.5% |
| N < 2^14 | 48.7% |
| N < 2^18 | 50.0% |
| N < 2^21 | 50.2% |

The first misses are 7, 13, 21, 29, 35, 41, 43, 49, ….

## 3. Chan–Pak as conflict-graph engineering

In the OCF frame every construction designs a conflict graph whose independence polynomial at 2 is computable:
- *Ordinal sums* give disjoint unions, so H multiplies; this is also why bottlenecks factor.
- *One apex* gives a clique: 1 + 2c.
- *Two apexes* give two cliques joined according to which cycles are disjoint: N = 1 + 2(a + b + c) + 4d.
  - a and b count the cycles through each apex alone, and c those through both.
  - d counts the disjoint pairs, a bilinear digital function of the tape.
  - This is the swarm's two-apex hosts.
- *Mod 8.* N ≡ 1 + 2m + 4d (mod 8), where m is the number of odd cycles and d the number of disjoint pairs. So the
  diagonal joins' class 3 (mod 8) of the synthesis packet reads m + 2d ≡ 1 (mod 4).

**The exceptions.** The odd numbers are exactly the values I(G, 2) over all graphs G, since cliques give every 2c + 1.
- 7 = I(K₃, 2) requires exactly three pairwise-intersecting odd cycles and no others. That conflict graph cannot occur
  in a tournament (odd-cycle crowding, THM-5192 (E)).
- 21 likewise has only non-realizable conflict graphs.

**The non-branching program.** Chan–Pak is the statement that I(·, 2) is onto the odd numbers other than 7 and 21 on
the class of *realizable* conflict graphs. A proof should come from a closure property of that class, not from more
families.
- One apex realizes cliques K_c for half of all c.
- The missing half must come from adding a second apex, or a disjoint cycle, at controlled overlap. That turns
  1 + 2c into 1 + 2c + 2(a + b) + 4d.
- The key lemma to seek is that the overlap counts (a, b, d) of two-apex tapes realize every residue class needed
  to fill the half left by one apex.
