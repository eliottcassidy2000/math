# Trierments: the char-3 edge algebra, the RXTX/Gram split, and a released Redei-extension

**opus-2026-07-28. Provenance note / wildcard record, not a truth source.**
Owner supplied: 3-state tournaments (missing and/or both-way edges), XOR,
and arXiv:2505.09814 (RXTX: a faster bilinear algorithm for `X X^t`,
found by ML + combinatorial search). Per the analogy-bridge pattern:
structural shape, not subject overlap.

## The structural analogy (typed)

RXTX beats generic matrix multiplication on `X X^t` by exploiting the
SYMMETRY of the output — decomposing a structured tensor more cheaply
than the unstructured one. The repo's tournament frontier has spent a
week discovering the same split from the physics side: THM-2501 (the
first even moment is signed-C4/GRAM energy while the skew tournament
quadratic form vanishes identically), THM-2346 (pair couplings are
symmetric cohabitation energies, not orientations), THM-2504 (the
intrinsic endpoint object is a typed directed weighted co-support
matrix — symmetric Gram part + antisymmetric current part + missing
entries — and forgetting types forces exactly 1,652 ties). A
**trierment** (3-state complete graph: forward / backward / tie-or-
both) is this decomposition made into a first-class object: the
symmetric support pattern (ties) is the Gram side; the orientation on
the rest is the current side.

**The char split.** Over F_2 (XOR) symmetric and antisymmetric
coincide — which is exactly WHY mac-mini-S143's Smith-theory no-go
killed mod-2 forcing of the 91-line (Z/2-equivariance with odd
7 (x) 13 coefficients vanishes): char 2 degenerates the Gram/current
split. Over F_3 = {-1, 0, +1} the split is canonical, and edge values
land naturally in F_3. The H-spectrum's 3-power rigidity (THM-2450's
C3-towers, THM-2444's pure-blue alphabet {1, C3, T5}, HYP-9029's
semigroup) then reads as the char-3 shadow of trierment structure, in
the same way Redei parity is the char-2 shadow of tournaments. The
correctly-primed forcing instruments (S144's rotation localization at
7 and 13) fit the same frame.

## The mini-census (exact, n <= 5) and a released law

Exhaustive Hamiltonian-path counts over all 3-state complete graphs,
both tie semantics (missing vs both-way), mod 2 and mod 3:

- Redei reproduced at t = 0 (odd count) in both semantics.
- n = 3 exact congruences (FINITE-EXACT, trivially provable):
  one both-way edge forces #HP in {2,3} mod 6; TWO both-way edges
  force #HP = 4 mod 6; all-both-way gives 0 mod 6. Missing semantics:
  one tie forces #HP in {0,1} mod 6.
- **RELEASED:** the n <= 4 pattern "odd #HP in a proper partial
  tournament forces #HP = 1 mod 3" is REFUTED at n = 5
  (5,424 violating configurations; all (t, mod-3) pairs realized
  through t = 3). A small-n artifact, recorded per kps-S134 method.

## What survives as questions (for a future session)

1. Is there a genuine F_3-valued path invariant of trierments
   (a "mod-3 Redei") — presumably not a naive count congruence, but a
   SIGNED count (orientation-weighted, the F_3 analogue of the skew
   determinant/Pfaffian pairing THM-2290 uses at char 0)?
2. The RXTX-shape question proper: the bilinear rank of the
   Gram-part tensor of THM-2457/2504's co-support matrix — does the
   structured decomposition (the analogue of RXTX's 5%) have a
   physical meaning as a smaller sidecar basis?
3. THM-2504's 1,652 forced ties: reread as a trierment invariant —
   is 1,652 = the Gram-part support count of a canonical trierment,
   and does its mod-3 class mean anything?

No claims; a frame, three exact facts, one honest refutation.
