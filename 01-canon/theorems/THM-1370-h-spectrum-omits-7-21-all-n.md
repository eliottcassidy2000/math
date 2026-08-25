# THM-1370: The Hamiltonian-path spectrum of tournaments omits 7 and 21 — for ALL n

**Status:** VERIFIED (complete proof; ingredients: two classical theorems + one-line
lemma + exhaustive small spectra, all frozen in-repo)
**Author:** boxeph-2026-07-20-S152 (HYP-8220)
**Predecessors:** S146 (n ≤ 6 + n=7 sample), S151 (n=7 symmetric strata),
S152 (n = 7 and n = 8 EXHAUSTIVE via augmentation covers)

## Statement

No tournament, on any number of vertices, has exactly 7 or exactly 21
Hamiltonian paths. Moreover the attained spectrum is monotone in n (padding),
and by n = 8 it covers every odd value in [1, 609] except 7 and 21.

## Proof

**(1) Multiplicativity.** h(T) = Π h(T_i) over the strong components of the
condensation (paths concatenate in the transitive order; S146, 300/300 checks
+ standard argument). Since 7 is prime and h = 1 only for transitive
components, h(T) = 7 requires a STRONG subtournament (a component) with h = 7;
h(T) = 21 requires a strong component with h ∈ {7, 21} (21 = 3·7 or 21·1).

**(2) Insertion lemma (one line, new).** For any tournament T and vertex v:
h(T) ≥ h(T − v). *Proof:* every Hamiltonian path of T − v extends to one of T
by inserting v at a feasible position (the classical Rédei insertion argument
guarantees at least one), and removing v recovers the original path, so the
map is injective. ∎

**(3) Strong floor is nondecreasing.** Let `f(n)=min{h(S):S strong on n
vertices}`. Moon's theorem says every strong tournament on `n>=4` vertices
contains at least two induced strong subtournaments on `n-1` vertices. With
(2), this gives `f(n)>=f(n-1)` for `n>=5`. [CITATION: J. W. Moon, "On
Subtournaments of a Tournament," Canadian Mathematical Bulletin 9(3) (1966),
297--301, Theorem 2, DOI 10.4153/CMB-1966-038-7; its formula
`s(n,k)=n-k+1` gives `s(n,n-1)=2`.]
This bibliographic correction is recorded as MISTAKE-505.

**(4) Exhaustive small spectra (frozen outputs).** strong(3) = {3},
strong(4) = {5}, f(5) = 9, f(6) = 15, f(7) = 25, f(8) = 45 — exhaustive at
n ≤ 6 (S146), n = 7, 8 by augmentation covers (S152: every (n+1)-tournament
minus a vertex is an n-tournament, so augmenting all canonical n-reps reaches
every (n+1)-class; A000568 counts verified: 56 → 456 classes; n=8 via the 456).
Also 21 ∉ strong(6) (else 21 were attained at n = 6 — it is not).

**(5) Conclusion.** By (3)+(4): f(n) ≥ 9 for all n ≥ 5, so no strong
tournament of any size has h = 7 (spectra {3}, {5} at n = 3, 4 miss it too);
by (1), h = 7 is never attained. For 21: 21 ∉ strong(n) for n ≤ 6 (4), and
f(7) = 25 > 21 with monotonicity kills n ≥ 7; so a strong-21 component never
exists, and 21 = 3·7 needs the impossible 7. Hence h = 21 is never attained. ∎

## The why, structurally

7 and 21 die in the CREVASSES of the strong floor: 7 ∈ (f(4), f(5)) = (5, 9),
21 ∈ (f(6), f(7)) = (15, 25), and both have no rescue by factorization
(7 prime; 21 = 3·7 recurses into the first gap). The floor data
3, 5, 9, 15, 25, 45 fits f(n) = min{3^a·5^b : 2a + 3b = n − 1} (Busch-type
3-cycle/5-block chains; prediction f(9) = 75, f(10) = 125). NOTE HONESTLY:
through n = 7 the floor also matches 2F(n)−1 (Leonardo); the two laws split at
n = 8 (41 vs 45) and the EXHAUSTIVE value is 45 — the Leonardo reading is
REFUTED, another instance of the repo's fifth-term-break motif (two laws
agreeing on exactly five terms). Literature (pinned S152 via web check): the
floor law IS the Moon–Busch theorem — Moon 1972 gave the chain construction
and bounds; A. Busch, "A note on the number of Hamiltonian paths in strong
tournaments," Electron. J. Combin. (2006), proved Moon's construction optimal.
Our exhaustive n ≤ 8 floors independently re-derive it. The spectrum-gap
statement ({7, 21} never attained, any n) we have NOT found stated anywhere —
it appears to be new as a stated theorem, though its proof is Moon–Busch plus
elementary gluing.

## Corollaries and data

- h = 63 IS attained at n = 8 (its absence at n = 7 was a strong-spectrum
  hole, not structure): the "7·3^k gap tower" reading is REFUTED at k = 2.
- The spectrum is monotone in n (pad with a dominated vertex: h unchanged), so
  spectrum(∞) ⊇ all odds in [1, 609] ∖ {7, 21} already.
- CONJECTURE (spectrum completeness): spectrum(∞) = ℕ_odd ∖ {7, 21}.
- n = 7 top structure: attained {159, 171, 175, 189} with holes between —
  h(QR₇) = 189 = max; n = 8 max = 661 with near-top holes {611,...,655} —
  the "ceiling gaps" move with n while the floor gaps {7, 21} are permanent.
- h(R_n) = 2^{n−2} + 1 for the almost-transitive family (transitive with the
  extreme arc reversed) — exactly the repo's hypotenuse constant H = 1 + 2^d:
  R_n is the transitive class's big SC neighbor on the principal line.
- With S146's dictionary (h-monoid ↔ Keller-degree monoid): the permanence of
  {7, 21} sharpens the tournament side of the odd-degree analogy — the h-monoid
  has EXACTLY two non-monoid odd holes, both explained by one floor law.

Scripts: `h_spectrum_n7_exhaustive_boxeph_S152.py`, `..._n8_...S152b.py`,
`almost_transitive_leonardo_boxeph_S152c.out`; reps file
`n7_class_reps_boxeph_S152.txt`. All outputs frozen.
