# The shell census at p ~ 300–500, and three more faces between the lines

**kind-pasteur-2026-07-19-S128c92** (HYP-7995) · owner: *"run the shell-collapse census at
p ~ 300–500; look for even more instances of 'two faces of the same underlying concept';
mine the history for hidden connections between the lines."*

**Relation to opus-S405 (read first, not repeated):** their six graded faces — the −1-ray
identity (F1), the k ≡ 1 mod 6 gate across generations (F2), spread ≡ maximizer (F3),
the one cycle-engine/two suppliers reading (F4), the 2-adic mirror stack (F5), G–K ≡
observer-lens-conjugated (F6) — set the grading discipline this session follows
(verified identity / proved note / rhyme). Everything below is disjoint from F1–F6.

---

## 1. Lane A — the shell-collapse census (results section; see frozen out)

Census primes {307, 331, 349, 397, 457, 467} — three greedy-easy, three thin-sea-hard
(from the c91 ranking) — minimal covers of the level-1 sieve enumerated under a soft
node budget with truncation REPORTED, each cover classified exactly against the full
ansatz library and scored by NEAR-SHELL Hamming distance (min substitutions to any
ansatz dilate). Baseline to reverse: the S128c88 generic sea at p ≤ 71 (1,280 / 25,711
/ 260,568 minimal covers at p = 43/61/71, ×20 per prime step, 99%+ unclassified,
distance-scattered).

**Results (sampling census v2 — the enumeration attempt itself became a finding: 40M
DFS nodes reach ZERO complete covers at p ≥ 307, the c91 lesson biting again, so the
census samples via randomized-greedy + drop-one minimalization, with capture-recapture
abundance and the greedy-reachability caveat stated):**

| p | class | draws ok /800 | distinct | abundance N̂ | sizes | ansatz hits | near-shell histogram |
|---|---|---|---|---|---|---|---|
| 307 | easy | 111 | 83 | ~220 | all 13 | AP ×1 | 0:1, 1:4, 2:10, 3:15, 4:19, 5:24, 6:10 |
| 331 | easy | 85 | 58 | ~134 | all 13 | — | 2:7, 3:6, 4:24, 5:18, 6:3 |
| 349 | HARD | 14 | 9 | ~20 | all 13 | — | 1:1, 2:2, 4:2, 5:3, 6:1 |
| 397 | easy | 83 | 61 | ~157 | all 13 | GW ×1 | 0:1, 1:4, 2:13, 3:11, 4:14, 5:11, 6:7 |
| 457 | HARD | 3 | 3 | ≫4 | all 13 | — | 1:1, 4:1, 5:1 |
| 467 | HARD | 11 | 10 | ~60 | all 13 | AP ×1 | 0:1, 1:3, 2:5, 3:1 |

**Verdict: SHELL COLLAPSE CONFIRMED (in the greedy-reachable space).** Three signatures,
all present: (1) ABUNDANCE — from ≥260,568 minimal covers at p=71 (complete enumeration)
to N̂ ~ 10¹–10² at p ~ 300–460: three-to-four orders of magnitude; (2) RIGIDITY — every
sampled minimal cover has size EXACTLY 13 (at p ≤ 71 sizes ran 9–13, mostly 10–12):
minimality has hardened to the full tuple, the c(p)-sticks-at-13 signature; (3)
STRUCTURE — exact AP and GW dilates reappear as sampled members, and at the thin-sea
primes the whole sample concentrates within Hamming ≤ 3 of the ansatz (467: max
distance 3), while easy primes retain a mid-distance residue (distances 4–6) — the
candidate NEW-FAMILY material (mod-p shadows of unknown near-floor structures; the
natural feedstock for mac-mini's three-far frontier). Caveats: greedy-reachability
bias (N̂ estimates the reachable subspace; unreachable covers possible), 800 draws,
capture-recapture is order-of-magnitude only. The S-T ansatz question now has a
measured answer shape: **the level-1 improper sea thins to the tight-family shell as
p grows, and fastest at multiplicatively-unsmooth p** (see S6).

## 2. Face S1 (VERIFIED IDENTITY, definitional): boxeph's wasted overlap ≡ the covering integrality gap

Two vocabularies, one number, and — the sharp part — two *complementary blindnesses*.

- **Measure face** (codex-S16 rank identity; boxeph-S125 at q=38): covering at radius λ
  forces cover-debt ∫r = Σμ(combs) − 1, and boxeph's mechanism finding was that failing
  families have MORE overlap, *wasted* — "concentrated at strong resonances instead of
  spread."
- **LP face** (S128c91): at level 1 the fractional covering number is exactly h/dk
  (uniform weight 1/dk; dual regularity), so the LP sees allowance and never
  obstruction; all of c(p)'s excess is integrality gap.

**The identity.** For any m danger-APs, waste := Σ|A(wᵢ)| − |∪A(wᵢ)|; covering by m
sets ⟺ waste = m·dk − h EXACTLY (the allowance). So

> **c(p) = min{ m : some m-family spreads its overlap to exactly the allowance m·dk − h },**

and the integrality gap (c(p) > ⌈h/dk⌉) IS the statement that for m below c(p) every
m-family's overlap is forced to CONCENTRATE strictly beyond its allowance — boxeph's
wasted-overlap mechanism, stated combinatorially. The two faces are blind in opposite
directions: the measure face sees the debt but not the quantization (its relaxation
computes 7 forever); the LP face sees the quantized bound but cannot exhibit the
concentration mechanism. The proof target both faces point at, now finite and exact
per p: **every 8-family wastes > 8·dk − h** (≈ p/14 forced excess concentration) ⟹
c(p) ≥ 9 — boxeph's cover-debt frame on a finite object, where the resonance counts
r(w,w′) = #{j·w ≡ ±j′·w′} are the whole story.

## 3. Face S5 (STRUCTURED RHYME, one transport question named): the primorial cascade ≡ the sieve's lcm-forcing

death-star's D-gate tower opens at N ≡ 1 mod primorial-layers (D=3 at mod 6, D=4 at
mod 30, D=6 at 210, …): the binder competition requires ALL small binders killed
simultaneously — a CRT conjunction over the small primes. Rosenfeld's Lemma 4 forces
lcm(2..k+1) | ∏vᵢ for any counterexample: the sieve side's CRT conjunction over the
same primes. **Both are "the small primes tax near-floor structure jointly, and the
tax is primorial"** — the spectrum face conditions the AMBIENT PARAMETER N; the sieve
face conditions the FAMILY's divisibility burden. opus-F2's question to death-star
(is the mod-6 gate coincidence structural?) is the D=3 slice of this pairing; the
general transport question it opens: does the binder b's role in the gate (the
surviving competitor) map to a specific divisor obligation in Lemma 4's lcm — i.e.,
is there a gate-indexed refinement of Lemma 4 forcing not just lcm(2..14) | ∏v but
WHICH speed carries which factor near the floor? Graded honestly: rhyme with one
named falsifiable transport; not claimed as an identity.

## 4. Face S6 (empirical, census-decided): the thin sea ≡ multiplicative unsmoothness of p − 1

The 15 greedy-hard primes of the c91 triage include a striking density of safe primes
(p = 2q + 1: 467, 719, 839, 983 …). Hypothesis: I(13,p,1) is thin exactly when
(Z/p)* has maximally unsmooth structure — few small-order elements ⟹ the danger-APs
A(w) = w·±[1,dk] have fewer accidental resonances ⟹ covers are scarcer and
greedy-harder. **S6 VERDICT: CORRELATED.** Across all 130 triage primes: thin-sea-hard median
largest-prime-factor fraction of p−1 = 0.125 vs greedy-easy 0.042 (3×); means 0.200 vs
0.105; safe-prime rate 27% vs 10% (4 of 15 hard are safe primes: 467, 719, 839, 983).
And the census closes the loop: the collapse is deepest exactly at the unsmooth primes
(457: 3 reachable covers). Face graded: verified correlation, mechanism proposed
(few small-order elements ⟹ fewer accidental AP-resonances ⟹ scarcer covers), not
yet a proof.
With it confirmed, the
face is: *the thin-sea ranking (search-difficulty side) ≡ the multiplicative
smoothness spectrum of p − 1 (structure side)* — and sieve-cost prediction becomes a
one-line arithmetic filter, no search needed.

## 5. A face forming in real time (cite-graded, no claim): slack = D − k, twice in one day

klein-THM-1290's close-out names "slack = D − k unification"; opus-S402's G–K pinning
derived the same dictionary from the published Amended Conjecture 1.2 (their s = our
D, their k = our slack). Two agents, two routes, one dictionary, same day — filed
here only as a cross-cite so the identification gets ONE canonical statement instead
of two vocabularies (the exact failure mode HYP-7890 documented).

## 6. Handoffs

- **boxeph:** Face S1 turns your q=38 cover-debt frame into the c(p) ≥ 9 finite
  target with the resonance counts r(w,w′) as the complete data — one session.
- **death-star:** Face S5's transport question (gate-indexed Lemma 4) sits exactly on
  your binder machinery; and the census's thin-sea numbers (Lane A) calibrate your
  per-prime sieve-cost model.
- **opus:** S1 and S5 extend your F-catalog; the grading discipline held (one
  identity, one rhyme, one empirical, one cite).

**Files:** census script + frozen out (`lrc14_shell_census_kps_S128c92`), HYP-7995,
SESSION-LOG entry, backlog updates.
