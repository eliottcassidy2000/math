# The Skew-Sylvester Doubling: One Cayley-Dickson Step, Seen Five Ways

**Source:** kind-pasteur-2026-06-09-S1 (THM-447..448, 450..453, T767, HYP-2332..2337 + 2344..2349, MISTAKE-065)
**Seed:** the human's directive — "three copies of the original and one negated copy; skew-Hadamard
matrices normalized with first row all 1s are analogous to the tiling model with the Hamiltonian
path fixed."

## The directive was exactly right

Normalizing a skew-Hadamard matrix (first row all +1) IS fixing the frame: border row = the
dominating source, fixed superdiagonal = the base Hamiltonian path, free C(n−1,2) core bits =
the tiles. The tiling model and the normalized skew-Hadamard core are the same moduli space.
And the doubling that respects this frame is the skew analog of Sylvester's [[H,H],[H,−H]]:
D(T) = [[M, M+I],[M−I, −M]] — three copies, one negated. Everything below unfolded from taking
that analogy literally.

## One operation, five shadows

The same doubling, looked at through five lenses, produced five structures that all turned out
to be the SAME fact wearing different clothes:

1. **Algebra (THM-450):** M' = P⊗M + Q⊗I with {P,Q} = 0, Q² = −I — a Cl(1,1) Clifford pair.
   The doubling is literally one Cayley-Dickson step applied to a tournament. The complete
   2×2-block family has exactly three orbits: two rank-1 projectors (T[K₂] on the
   swap-symmetric line, SCblow on the swap-antisymmetric line) and one mixer (D, the 2×2
   Hadamard). Skew-Hadamard preservation ⟺ Chebyshev spectral law ⟺ tr P = 0 ⟺ det P = −2.
   **The 2×2 Hadamard matrix is the doubling's DNA.**

2. **Spectra (THM-447):** eigenvalues obey μ'² = T₂(μ) — Chebyshev doubling; the conserved
   charge λ²+1 (= Hadamard order) exactly doubles; asymptotic stretch √2 = the triangle's
   hypotenuse ratio. The rank-1 orbits instead give unimodular lift pairs (product exactly 1):
   boost pairs, echoing the THM-251/252 rapidity lattice. Projectors give boosts; the mixer
   gives Chebyshev. One pencil det(λP + Q) controls both.

3. **Tilings (THM-452):** in the canonical frame (path, twin, reversed path — which for the
   iterated tower is the BINARY REFLECTED GRAY CODE), the doubled tiling carries the Sylvester
   sign pattern at tile level: two σ-exchanged positive copies on the diagonal blocks, the
   ordered arc matrix interleaved as complementary σ-pairs in the cross block (the negated
   copy), and an all-ones twin anti-diagonal — the hypotenuse, fractally repeated at every
   scale of the tower. The blue conjecture failed but failed PERFECTLY: the grid-symmetry
   defect is the CONSTANT vector c_n, independent of T. D embeds tiling space as an affine
   section of one fixed coset of the blue subspace, and grid transpose implements the
   op-functor as translation by c_n. The doubling doesn't make things blue — it makes the
   FAILURE of blueness rigid.

4. **Number theory (THM-448):** the tower visits the Mersenne numbers 2^k − 1 and stays DRT
   forever (proved: it is BORN normalized — closed form S[i,j] = (−1)^(popcount((i&j)≫(v+1)) +
   bit_v(i)), the "skew-Walsh function": Sylvester character truncated at the twin scale).
   It greets Paley at 7 and then leaves: T31 is a non-Paley, non-circulant DRT(31). At 15,
   where NO circulant DRT exists at all, the tower still delivers one — with
   H(T15) = 198,335,025 ≈ 2.485·E[H], the natural n=15 H-maximizer candidate.

5. **Hadamard classification (THM-451):** the tower tracks Sylvester through order 8, then at
   16 lands precisely in the UNIQUE transpose-split pair {had.16.3, had.16.4} of Hall's five
   classes — it is not even equivalent to its own transpose.

## The two big morals

**Chirality is conserved across all five lenses.** D is the unique non-op-equivariant doubling
(H(D(T)) ≠ H(D(T^op)) in 50/74 classes); the tower cores have ZERO anti-automorphisms (pure
black tilings, never self-converse, unlike Paley); and the Hadamard equivalence classification
detects exactly this: the tower lands in the only chiral (transpose-split) class pair. The
op-asymmetry introduced by choosing WHICH copy gets negated is not a blemish — it is the
signature the construction carries through every representation. Where Paley DRTs are
maximally symmetric (self-converse, vertex-transitive, |Aut| growing), the tower DRTs are
maximally rigid (Aut frozen at F_21 forever, the Frobenius group inherited from Paley T_7 at
the moment the tower last touched Paley). The doubling preserves the group but never feeds it.

**Recursion beats basis.** The tower matrix is provably MAXIMALLY DENSE in the Walsh basis
(every Walsh-domain entry ≡ 2 mod 4 — no zero can exist), yet admits an exact O(m log m)
butterfly because S = W + Corr with a self-similar correction tower. Sparsity lives in the
recursion, not in any fixed basis. This is the same lesson the metagraph keeps teaching
(Mode A/Mode B descend recursively where no single invariant suffices) — and it is now an
ENGINEERING deliverable: the skew-Walsh transform (27× vs BLAS at m=4096, 82.5 ms at m=2^20),
a fast transform whose kernel is a DRT fingerprint at Mersenne sizes.

## The parity mixer

The deepest single finding may be the smallest: twin insertion converts EVEN cycles of T into
ODD cycles of T[K₂] (c5' = 32c5 + **32c4** + 6c3). The OCF sees only odd cycles; the doubling
pumps the invisible even-cycle sector into the visible odd sector. That is why
I(Ω(T),x)-determination of H(T[K₂]) fails (the one refuted hypothesis of the session that
stays refuted) and why the cycle SPECTRUM works. Doubling mixes parity the way the Walsh
recursion mixes frequencies — the d=1 wiggly layer and the d=m blue/black layer of the waggly
spectrum are the two ends of exactly this mixing. A future session should chase: what does
the doubling do to the waggly LAYER decomposition as a whole?

## Resonances to carry forward

- T15's vertex links have H ∈ {189, 171} — two of the three rigid regular-7 classes of the
  S18h BIBD trichotomy. The tower's interior remembers the trichotomy.
- Link-H (H of the out-neighborhood) is a strictly stronger vertex invariant than every cheap
  per-vertex statistic (splits T31's vertices 5 ways where scores/c3/cycle-counts are constant
  and equal to Paley's). New tool for the metagraph: fiber H along links.
- B_0(T_{2m−1}) = T_{m−1} as a LABELED submatrix: the tower is its own link recursion — a
  Mode-A-like descent (vertex deletion) that exactly inverts the doubling. Mode A removes the
  hypotenuse; the doubling ADDS a hypotenuse of twins. The triangle picture closes.
- The Mersenne tower 2^k − 1 and the Cayley-Dickson tower 2^k + 1 are the two sides of the
  same doubling, one THM-067 vanishing apart. Both lose a property per level (CD: algebraic;
  here: Paley-ness at 15, transpose-equivalence at 16 — the SAME level).

## What to do next (ranked)

1. **HYP-2346:** prove Aut ≡ F_21 and the labeled link recursion from the closed-form arc law
   (both look inductive; the closed form makes them finite-check-per-level).
2. **HYP-2347:** n=7 stress test — separate cycle-spectrum from pair statistics in THM-453(5);
   also first nonzero-c7 test of the 128 coefficient.
3. **HYP-2345:** is H(T15) the n=15 maximum? Needs the score-class reduction (SC maximizer
   theorem machinery) — brute force is out of reach.
4. **HYP-2344:** the consecutive-circulant K₂≅SCblow coincidence — find the bijection, test n=9.
5. Waggly-layer transport under doubling (the parity-mixer question above).
