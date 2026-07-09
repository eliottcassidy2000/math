# THM-674 — The domination theorem for descent blocking: blocked mod prime k ⟺ the inverse-class set T-dominates the ±-class group (T2 is the j=1 case; k=29 is cycle domination)

**Status:** PROVED (Theorem + Corollaries 1–3 + the general ledger, elementary, proofs below) +
VERIFIED with ZERO violations: the T·D = G characterization at k = 29, 31, 37, 41 (j=2) AND
k = 43, 53 (j=3), 170k sampled subsets + 20k multisets each (blocked counts 571–52,761 per
modulus); Corollary 3's per-cycle bound 0 violations and TIGHT at k = 31, 37, 41 (min observed
#classes = bound: 9, 9, 10); the general ledger formula exact for every k ∈ [29,42], every r.
**Dodger data:** cap 150 could not fully dodge [15,42] (min 3 open descents); the cap-250 full
[15,42]-dodger was caught at k = 49 = 7² (q = 98) — one band deeper again. Blocking probability
for random covering sets collapses with k: 13.8% (k=29) → 6.2% (31) → 0.75% (37) → 0.25% (41) —
the SPREAD demand is rapidly unaffordable without AP-like tuning. The SIMULTANEITY conjecture
(Part 4) remains OPEN — it is the covering branch.
**Source:** mac-mini-2026-07-09-S65 (cont. 3). Answers the owner's ask: the k≥29
cycle-domination statement is Corollary 2, proved.
**Depends on:** THM-668 (C2 descent), THM-672 (window anatomy = the j = 1 case).

**Setup.** Modulus `k > 14`, danger radius `j := ⌈k/14⌉ − 1 ≥ 1` (the closed band
`[⌈k/14⌉, k − ⌈k/14⌉]` has danger complement `{0, ±1, …, ±j}`). Residue multiset
`R = {v_l mod k}`, 13 entries, no zeros (a zero residue makes the descent dead). "Blocked" = no
`s ∈ [1, k−1]` puts all of `R·s` in the band; C2 (THM-668) fires through `k | q` iff unblocked.

## Theorem (prime k: the domination form)

> Let `k` be prime, `G = Z_k^*/{±1}` (cyclic of order `m = (k−1)/2`), `C ⊆ G` the occupied
> classes of `R`, `D = C^{−1}`, and `T = {classes of 1, 2, …, j}` (these are `j` distinct
> classes since `j < k/2`). Then
> **blocked ⟺ T·D = G**, i.e. in log coordinates `⋃_{t ∈ ind(T)} (ind(D) + t) = Z/m`.

*Proof.* For `s ∈ [1, k−1]` and unit `r`: `rs ≡ 0` is impossible, so `s` is bad for `r` ⟺
`rs ≡ ±d` for some `d ∈ [1, j]` ⟺ `s ≡ ±d·r^{−1}` ⟺ `class(s) ∈ T·{class(r^{−1})}`.
Hence `s` is bad for SOME `r ∈ R` ⟺ `class(s) ∈ T·D`, and blocked ⟺ every class is bad
⟺ `T·D = G`. ∎

## Corollary 1 (j = 1 ⟹ THM-672's T2)

`k ∈ [15,28]` prime: `T = {1̄}`, so blocked ⟺ `D = G` ⟺ `C = G` (inversion is a bijection):
every ±-class occupied — exactly T2. The window theorem is the degenerate domination.

## Corollary 2 (k = 29: the cycle-domination statement — PROVED)

`j = 2`, `T = {1̄, 2̄}`. Since `2` is a primitive root mod 29 (`2^14 ≡ −1`), `2̄` generates
`G ≅ Z/14` and we may take `ind(2̄) = 1`. Then
> **blocked mod 29 ⟺ `ind(D) ∪ (ind(D) + 1) = Z/14` ⟺ the complement of `ind(D)` contains no
> two consecutive elements of the 14-cycle** (a dominating set of the directed 14-cycle).
∎ (Sanity: the AP `{1..13}` occupies 13 of 14 inverse classes, missing only one — dominated —
blocked ✓, matching its known everywhere-blocked behavior off its own modulus.)

## Corollary 3 (per-cycle counting bound, j = 2 primes)

Let `o = ord(2̄)` in `G` (additive order of `t₂ = ind 2̄`). Translation by `t₂` decomposes `Z/m`
into `m/o` cycles of length `o`; dominating a length-`o` cycle with `{0, +1}` needs `⌈o/2⌉`
elements. Hence
> **blocked ⟹ #occupied ±-classes ≥ (m/o)·⌈o/2⌉.**
Table: k=29 (o=14): ≥ 7 of 14; k=31 (2⁵ ≡ 1, o=5): ≥ 9 of 15; k=37 (o=18): ≥ 9 of 18;
k=41 (2¹⁰ ≡ −1, o=10): ≥ 10 of 20. With only 13 residues available, blocking k = 41 commits
10 of 13 class-slots to a rigid pattern — the SPREAD demand (Part 4).

## General ledger (all k, all j — the j ≥ 2 correction to THM-672's Lemma 1)

> For `g = gcd(r, k)`: `|A_r \ {0}| = (g − 1) + 2g·⌊j/g⌋`.

*Proof.* `rs ≡ d` is solvable iff `g | d`, with exactly `g` solutions. `d = 0` gives `g − 1`
nonzero solutions; each `d = ±e`, `e ∈ [1, j]`, `g | e` (there are `⌊j/g⌋` such `e`) gives `g`
solutions, all distinct across distinct `d`. ∎
NOTE: for `j ≥ 2` non-units reach the danger elements divisible by `g` — the [15,28] nesting
(`A_r = (k/g)Z`) breaks; the necessity ledger `Σ_classes |A \ 0| ≥ k − 1` is the correct
composite generalization, with per-stratum domination (units mod `k/h` for strata `h ≤ j`) as
the exact composite characterization — instance-checked by script rather than tabulated.

## Part 4 — the spread-vs-concentration tension (the remaining core, exactly)

The certificate hierarchy now has a proved shape at every modulus:
- `j = 1` (k ∈ [15,28]): blocking demands **torsion CONCENTRATION** (THM-672 T3) — which the
  covering condition SUPPLIES (odd multiples of 13 sit at 13 mod 26, etc.).
- `j ≥ 2` (k ≥ 29): blocking demands **class SPREAD in a dominating pattern** — ≥ (m/o)⌈o/2⌉
  well-placed classes of at most 13; nothing in the covering condition supplies this. The tight
  AP achieves it because the inverse set of an interval is a maximally-spread Farey fan — one
  more face of "the extremal is the maximally coherent object".

**The remaining open statement (the covering branch of LRC(14), sharpest known form):** no
primitive covering 13-set can simultaneously (a) torsion-occupy every window divisor `k ∈
[15,28]` of its pair sums, and (b) T-dominate every prime `p ≥ 29` dividing a pair sum, for all
the ~91 sums at once. The extended dodger search (results file) measures how deep adversaries
get; every dodger of [15,28] found so far fails (b) at the first prime tried.

**Verification:** `04-computation/lrc14_domination_theorem_macmini_S65cont3.py` (+ .out).
**Related:** THM-668, THM-672, klein THM-671 (quintic Bonferroni — the counting complement on
these moduli), kps-S116 LRCPairSumDispatch (+ ratioBand = C0 in Lean; the (q,p) consumer),
boxeph-S3 (mid-band Erdős–Turán form), opus LEM-015 / kps LRCSchurRigidity (global coherence
rigidity).
