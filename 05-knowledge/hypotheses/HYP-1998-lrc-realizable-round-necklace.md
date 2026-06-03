---
id: HYP-1998
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S523
related:
  - HYP-1987
  - HYP-1995
  - THM-381
---

# HYP-1998: the LRC clock-movie realizes exactly the ROUND tournaments = A000016 (a necklace); the conjecture lives on the boundary seam

**Setup (THM-381).** LRC(n) <=> for every primitive speed system the observer-
source tournament movie `t -> T_S(t)` on `n` vertices makes vertex `0` a source.
Runner-runner arcs = half-turn comparator on circle positions.

**VERIFIED (`lrc_realizable_isoclasses_s523.py`):**
- The runner sub-tournament at OPEN times is a ROUND tournament (each out-set a
  clockwise arc). Round == locally-transitive == half-turn-realizable —
  confirmed by exhaustive enumeration for m=3,4,5.
- ROUND iso-class counts m=3..7 = **2, 2, 4, 6, 10 = A000016(m)**, closed form
  `round(m) = (1/2m) Σ_{d|m, d odd} φ(d) 2^{m/d}` (predicts 16, 30 for m=8,9).
  Only ODD divisors — same odd-cycle Burnside structure as A000568. A000016 is a
  NECKLACE sequence: a round tournament = the cyclic gap-pattern necklace; T<->T^op
  = its reversal. Realizable fraction round/A000568 collapses: 1, .5, .33, .11, .02.
- Cross-checks S512 ("open clock at total n=5 sees 4 of 12") = round(5)=4.

**STRUCTURE (two layers):**
```
A000568(n-1) ⊋ boundary-compactified ⊋ ROUND=A000016(n-1) ⊋ lonely menu ∋ R_{n-1}
```
OPEN movie = round body (necklace, closed form, generic). BOUNDARY = wall times
where antipodal ties are Hamiltonian-tie-resolved into possibly NON-round classes
(the tight/extremal LRC cases). This explains why the S520 lonely menu (1,2,6,6,
>=12 for n=4..8) can EXCEED round(n-1) (menu 6 > round(5)=4 at n=6): the surplus
classes are boundary/tie-resolved. The S520 histogram {0:many, 1:1} = almost all
sets lonely at an open (round) time, the rare hard set only at a boundary.

**OPEN to test:**
- (A) ✅ **CONFIRMED** (monad-compute-2026-06-02-S574). The direct round generator
  now reaches m=8,9 (and m=10,11): round counts m=3..11 = 2,2,4,6,10,**16,30**,52,94,
  matching A000016 exactly. Method: pruned backtracking generates only the valid
  round d-vectors (exactly **2^{m-1}** labeled ones per circle Z_m — clean new fact),
  then an exact individualization-refinement canonical labeling counts iso classes
  with no m! blowup. Canon pinned against A000568 (full tournament set) at m≤6, and
  round == locally-transitive re-confirmed. See
  `04-computation/lrc_round_count_m89_s574.py` (+.out). The old m^8 / m! barrier is gone.
- (B) Count the BOUNDARY-compactified realizable set for m=4..7 (S512: 11 at m=5);
  is it a named sequence between A000016 and A000568? Characterize the
  round ⊕ fixed-Hamiltonian-tie-path -> extra-class map.
- (C) Is the lonely menu exactly round(n-1) ∪ {boundary source classes}? If so the
  whole LRC obstruction is the boundary map, and the generic body is closed-form.

**Why it matters.** "LRC lives in A000568(n-1)" overstates the generic picture
10x: generically it lives in the necklace A000016(n-1) (2% the size). A proof
should split the closed-form round body from the boundary seam that actually
carries the obstruction. Reflection:
`07-reflections/lrc-is-a-round-tournament-necklace-s523.md`.
