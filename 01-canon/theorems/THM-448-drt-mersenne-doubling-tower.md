# THM-448 — The DRT/Mersenne Doubling Tower (STUB — claimed by kind-pasteur-2026-06-09-S1)

**Status:** CLAIMED / CONJECTURE — computation in progress this session. Honest stub: statement
below is the proof target; only the skew-Hadamard-order algebra (THM-447 claim 2) is proved.
**Source:** kind-pasteur-2026-06-09-S1.
**Related:** THM-447, HYP-2333, THM-067 (Mersenne vanishing), Cayley-Dickson tower
(`07-reflections/lrc-cayley-dickson-tower-s387.md`).

## Statement (target)

Let T be a doubly regular tournament (DRT) on n vertices, so S(T) bordered with a +1 row /
−1 column is skew-Hadamard of order n+1. Applying the skew-Sylvester doubling (THM-447) to the
bordered matrix and re-normalizing yields a skew-Hadamard matrix of order 2(n+1) whose core is a
DRT on **2n+1** vertices. Hence:

1. DRTs propagate along the tower n → 2n+1; from the empty/1-vertex seed the tower visits
   exactly the **Mersenne numbers 2^k − 1**: 1 → 3 → 7 → 15 → 31 → ⋯
2. **C_3 → Paley T_7** (conjecture: the doubled 3-cycle is isomorphic to the unique DRT on 7,
   the Paley tournament).
3. The 15-vertex member: identify it among the known DRTs of order 15; compute H, |Aut|,
   spectrum; compare with Paley-type H-maximization (memory: Paley = H-maximizer at p=3,7,11).
4. The Mersenne tower 2^k − 1 is the complement of the Cayley-Dickson tower 2^k + 1; both are
   driven by doubling. THM-067's Mersenne vanishing (c_1^{(f,d)} = 0 iff n = 2^{f+1} − 1) should
   be re-read as a statement about the doubling tower's fixed frame.

## Evidence so far

- THM-447(2): skew-Hadamard order doubles under D (proved).
- Spectral charge λ²+1 doubles; tower spectral radii √(2^k − 1) match DRT on 2^k − 1 (consistent).
- Computational verification: PENDING (this session, skew_doubling_tower_kps1.py).
