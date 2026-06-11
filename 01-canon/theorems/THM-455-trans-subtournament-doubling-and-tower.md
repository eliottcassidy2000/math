# THM-455 — Transitive Subtournaments under Skew-Sylvester Doubling: the +1 law (with one exception), and the Mersenne tower's extremality window 7–31

**Status:** (1) PROVED. (2) PROVED for transitive T; the universal +1 conjecture REFUTED with
exact counterexample AND shown Ramsey-FORCED to fail (2b). (3)–(4) COMPUTED and INDEPENDENTLY
VERIFIED (source ≡ sink recursion; explicit witnesses; brute force 18/18 at n≤5; my
trans(Paley_31) = 7 matches the published Momihara–Suda Table 2 value — external validation
of the solver). Literature fully verified (agent report with sources; Erdős problem #1216).
**Source:** kind-pasteur-2026-06-09-S2 (T769, HYP-2356/2357; computed by session owner after
two branch-agent API failures).
**Related:** THM-447/448 (doubling, tower, labeled link recursion), THM-454, Erdős–Moser
(1964), mac-mini THM-453 Part A (Erdős 592's red target = order-faithful transitive
subtournament — this theorem is its finite shadow).

## (1) Lower bound: trans(D(T)) ≥ trans(T) + 1 (PROVED)

If S = {v₁ → ⋯ → v_t} is transitive with source v₁, then S ∪ {v₁'} is transitive in D(T):
v₁ → v₁' (twin), v₁' → v_j for j ≥ 2 (cross arcs follow T). ∎

## (2) The +1 law holds in 17/18 classes at n ≤ 5 — but is NOT universal

trans(D(T)) = trans(T) + 1 for 17 of the 18 iso classes n=3..5 (verified by fast recursion AND
brute force). For transitive T it is exact (the interleaving obstruction: primed elements must
appear in reversed T-order while cross arcs pin each v_j' to the slot after v_j — two primes
cannot straddle an unprimed element, so a chain-comparable configuration gains at most 1).

**EXCEPTION (refutes universality):** n=5 idx10, scores (1,2,2,3,2), arcs
{(0,4),(1,0),(1,3),(2,0),(2,1),(3,0),(3,2),(3,4),(4,1),(4,2)}: trans(T) = 3 but
trans(D(T)) = 5 (+2). Witness TT5 in D(T): the chain 1', 3, 2', 0, 4' — primed and unprimed
ALTERNATE, using cross arcs in both directions. The interleaving obstruction binds only
chain-comparable configurations; alternating mixed chains evade it. **n=6 census (ALL 32768 labeled tournaments): delta = trans(D)−trans(T) distribution
{+1: 27248, +2: 5520} — delta NEVER exceeds 2.** The +2 rate grows with n (1/18 classes at
n≤5 → 16.8% labeled at n=6), concentrated on (t,t+2) = (4,6) and (3,5). This upgrades the
trivial bound trans(D(T)) ≤ 2·trans(T) to the SANDWICH CONJECTURE
trans(T)+1 ≤ trans(D(T)) ≤ trans(T)+2 (HYP-2360).

## (2b) The exceptions are partly FORCED by directed Ramsey numbers (new, verified)

The literature notes no doubling lower-bound construction can exist (it would contradict
Reid–Parker). Sharply: if T is TT_{t+1}-free on n vertices and 2n ≥ R(t+2) (the directed
Ramsey number forcing TT_{t+2}), then EVERY 2n-vertex tournament — including D(T) — contains
TT_{t+2}, so the +1 law MUST fail. First instance: Paley T₇ (t=3, 2n=14 = R(5), Reid–Parker
1970). VERIFIED: trans(D(Paley₇)) = 5 = 3+2 exactly; also trans(SCblow(Paley₇)) = 5,
trans(Paley₇[K₂]) = 6 — the skew double is the SLOWEST-growing doubling on Paley₇. The n=5
and n=6 exceptions are NOT Ramsey-forced (2n < 14): two species — Ramsey-forced and
structural (alternating-chain). Conversely, the sandwich (HYP-2360), if proved, gives
R(k+2) ≥ 2R(k) − 1 (weaker than known records; structural value only).

## (3) The tower values (all independently verified, witnesses explicit)

```
           T7    T15    T31    T63          controls (10 random each)
trans:      3      5      7     11          n=15: 7,7,7,7,7,7,7,7,8,8
                                            n=31: 9,10,10,10,10,10,10,10,10,11
Paley_31:                7                  n=63: 11,12,12,12,12,12,12,12,13,13
```
Tower-step growth: +2, +2, +4. TT7 witness in T31: (8,25,22,4,11,29,18); TT11 witness in
T63: (23,56,32,2,26,61,37,6,30,57,33).

## (4) Reading — literature-confirmed (Erdős problem #1216; R(2..6) = 2,4,8,14,28;
## 34 ≤ R(7) ≤ 47, NMH 2022 + McKay; trans(QR_q) = 3,4,5,5,5,7 for q = 7,11,19,23,27,31)

- **T7 = Paley T₇ = ST₇, THE unique extremal TT₄-free tournament** (R(4) = 8; Erdős–Moser 1964).
- **T15 achieves trans = 5 = f(15) exactly** (R(5) = 14 forces TT₅ at 15; TT₆-free exists at 15
  inside ST₂₇). Pointwise f-extremal — and T15 matches at order 15 the value trans = 5 that QR
  tournaments only achieve at orders 19–27 (Momihara–Suda Table 2).
- **T31 = 7 = trans(QR₃₁) exactly** (published value confirmed): tower and Paley tie at 31
  despite non-isomorphism (THM-448c). Both are beaten by the Sanchez-Flores 1998 record: a
  NON-QR circulant on 31 vertices (increments {2,3,6,10,11,13,15,17,19,22,23,24,26,27,30})
  that is TT₇-free (trans = 6) — so f(31) = 6 and both T31 and QR₃₁ are +1 above pointwise
  minimal. Notable: the extremal objects at 13 (ST₁₃, fourth-power-residue circulant) and at
  31 are NOT QR tournaments either — extremality keeps passing to sporadic circulants.
- **T63 = 11 decouples**: R(7) ≤ 47 forces TT₇ at 63, and T63 carries TT₁₁ — far above the
  floor, at the low end of the random range (11–13) but no longer extremal.
  The tower's Erdős–Moser extremality window CLOSES between 31 and 63 — mirroring how the
  tower leaves Paley after 7 (THM-448) and leaves Sylvester-equivalence at order 16 (THM-451).
  Every lens (isomorphism, Hadamard class, now extremal combinatorics) sees the tower stay
  classical for a few levels and then go its own way; remarkably its |Aut| = F₂₁ stays frozen
  through all of it.
- The tower remains BELOW the random median at every computed level (3<7-8, 5<7-8, 7<10,
  11<12) — the doubling structure suppresses transitive subtournaments at every scale, just
  not extremally at 63.

## (5) Bridge to Erdős 592 (mac-mini's thread — coordination, not duplication)

THM-453 Part A reframes ω^β → (ω^β, 3)² via order-faithful transitive subtournaments. This
theorem is the finite-quantitative shadow: explicit tournaments with controlled trans give the
finite obstruction landscape that 592's witness colorings must navigate. The doubling +1 law
(and its alternating-chain exception) is exactly the kind of "pattern-evading mixed chain"
mac-mini's E.3 found at ω³ (witnesses must use gap magnitudes, not pure patterns).

## Literature (all web-verified with sources; full report in session record)

Erdős–Moser 1964 (Magyar Tud. Akad. 9, 125–132); Stearns 1959 (f(n) ≥ ⌊log₂n⌋+1);
Reid–Parker JCT 9 (1970) 225–238 (R(5)=14, ST₁₃ unique); Sanchez-Flores GC 10 (1994)
(R(6)=28, ST₂₇ = QR/GF(27) unique) and GC 14 (1998) (31-vertex record circulant);
Neiman–Mackey–Heule, Graphs Combin. 38 (2022), arXiv:2011.00683 (SAT: 34 ≤ R(7) ≤ 47;
TT₆-free catalogs; 23-vertex count 25 = 22 ⊂ ST₂₅ + 3 DRTs);
Lidický–Pfender SIDMA 35 (2021) (flag algebras, R(TT_k) ≤ 53·2^{k−7});
Momihara–Suda arXiv:1605.02469 (trans(QR_q) table; DRT eigenvalue bounds ≈ √(3v) — Tabib 1986);
asymptotic gap [√2, 2] in the growth base of R(k) remains OPEN (no constant improvement since
1964/1970; the conjecture f(n) = ⌊log₂n⌋+1 fails first at n=14 yet holds on 16 ≤ n ≤ 27).
NOTE: the lexicographic product is multiplicative (trans(G[H]) = trans(G)·trans(H)) — our D is
an ADDITIVE-increment doubling, a different species from anything in this literature.

## Scripts

`trans_tower_erdos_moser_kps2.py`, `verify_trans_tower_kps2.py`,
`trans_doubling_n6_kps2.out` (n=6 census), `trans_paley7_doubles_kps2.out` (forced exception).

---

## ADDENDUM (kind-pasteur-2026-06-11-S2): HYP-2360 REFUTED; the sandwich is replaced by the zigzag law

THM-482 refutes the +2 sandwich with the unbounded alternating family A_l
(trans(A_l) = l+1, trans(D(A_l)) ≥ 2l+1; A_2 = this file's n=5 idx10 exception;
delta = 3 already at n = 7, one vertex past the n=6 census). The exact law is
trans(D(T)) = z(T), the zigzag number (THM-482 B). Consequences for this file:
(2b)'s conditional remark "the sandwich would give R(k+2) ≥ 2R(k)−1" is moot;
the tower's step pattern +2, +2, +4 is now read as the zigzag numbers of the
levels (z(T₃₁) = 11). The trivial bound trans(D(T)) ≤ 2·trans(T) is
asymptotically tight. See THM-482.
