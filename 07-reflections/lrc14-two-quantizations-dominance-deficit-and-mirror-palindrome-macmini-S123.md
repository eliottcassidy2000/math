# The two quantizations: how the LRC(14) near-misses converge, and what the tournament/metagraph frame adds

**Instance:** mac-mini-2026-07-19-S123 (owner-directed: "work to finish LRC(14); see the times we
have come close; synthesize multiple valid perspectives; focus on tournament and metagraph
analogies").

**Standing on, not repeating:** opus-S399 (era map, five-axis triage, 2.5-wall census),
boxeph-S130 (fourteen near-completion episodes, certificate-rung ladder), death-star-S58 (five
failure genera), kind-pasteur-S128c86 (three-kernel cross-walk). Read those first; this document
adds (i) the same-day convergence picture, (ii) two new proved/verified instruments from the
tournament lens, (iii) the typed next moves they shape.

---

## 1. Where the campaign stands tonight (2026-07-19, after a single extraordinary day)

Four independent lines advanced Wall A *today*:

- **THM-1289 (opus-S402, citation-grade):** the floor is ISOLATED from above — some δ>0 has no
  13-speed M-value in (1/14, 1/14+δ) — via Giri–Kravitz Thm 1.4 (published). Ineffective δ.
  The accumulation horn of THM-1268's dichotomy is dead.
- **HYP-7920 (opus-S400, the S-T tightness cage):** below 113/1466 every primitive 12-family is
  dilated-AP-congruent mod every sieve prime; rigidity + micro-gap emptiness to height 258,276.
  Effective but bounded height — the exact complement of THM-1289's all-heights-but-ineffective.
- **kps-S128c87/c88 (feasibility ledger + acceptance test):** the S-T architecture priced for
  n=14 at ~50 core-years; level-1 pruning REFUTED (B/C≈1.06 at p=71) ⟹ the missing ingredient is
  lift-level kill theorems — pruning theorems of exactly the repo's kind (mod-p spread, cages).
- **boxeph-S132 (in-band census):** the {3,35}-straddled band slice is EMPTY of gap families
  outright (min over all five-rung survivors = 3/29, above the whole gap) — the blocking-budget
  argument measured end-to-end on a complete slice.

Read against the wall census: **Wall A is being squeezed from four sides simultaneously** —
isolation (ineffective, all heights), congruence cage (effective, bounded height), finite-check
pricing (needs pruning theorems), and complete slice censuses (bounded certificates saturate at
rung values). The four have non-overlapping weaknesses. That is what convergence looks like; none
of them alone finishes, and each names the same residual: **effectivize the isolation radius δ /
kill the finitely many rungs in (1/14, 3/41]**.

Wall B (codex's phase-transport program) moved concurrently: THM-1283 gained the Hunter-star
joint tail tax (in-flight in this worktree during this session), THM-1287-selected-prefix is
reserved. Wall B remains genuinely distinct and protected.

## 2. The convergence, read as the tournament corpus reads Rédei

The repo's founding theorem (Rédei) isolates the transitive tournament two ways at once:

- **globally**, by a parity quantization — H is always odd, so H=2 is impossible and the H=1
  pole has a spectral gap below H=3, proved by an involution pairing Hamiltonian paths with a
  unique fixed point;
- **locally**, by an exact flip formula — one-arc flips of the transitive class land at
  H = 1+2^d: the neighborhood spectrum is computed, not just bounded.

LRC(14)'s Wall A is now visibly the SAME two-part statement about the M-spectrum at the AP pole:

| | metagraph G_n at transitive pole | M-spectrum at the AP pole (n=14) |
|---|---|---|
| global quantization | H odd (Rédei involution) | D = M·s ∈ ℤ, the rung lattice (THM-1269); G-K: no value approached from above (THM-1289) |
| isolation | H=1 isolated, no H=2 | (1/14, 1/14+δ) empty — but δ INEFFECTIVE |
| local flip spectrum | H = 1+2^d, exact | rung realizability: ladder c/(cN+1) (THM-633), gate tower D/Q (THM-1286), slack-1 realized at D=1,2,3 — INCOMPLETE map |
| the open content | (none — both parts proved) | make δ effective = complete the local flip-spectrum map |

boxeph-S110 proposed transporting "transitive-class isolation" to M; today the transported
statement EXISTS (THM-1289) — by citation, not by transport. What the metagraph frame still owns
is the shape of the missing half: **Rédei's isolation became effective only when the local flip
formula was computed. The M-side analog is the rung-realizability map near the AP (which
(D,s)-rungs are attained, at which heights) — death-star's cross-N gate census and the slack-1
realization work ARE the H=1+2^d computation of this problem.** That is the precise, typed
next target for Wall A, and it is already an assembly line, not a new idea.

## 3. The two new instruments from this session (both exact, both 2-adic)

### 3a. The dominance deficit law (HYP-7960, resolved)

Define the pair-dominance measure P(a,b) = meas{t : ‖at‖ > ‖bt‖}. Then for coprime a<b:

    P(a,b) = 1/2 − 1/(2(b²−a²))   if a+b is odd (mixed parity),
    P(a,b) = 1/2                  exactly, if a,b both odd.

Proof: (both-odd) for odd a, ‖a(1/2−t)‖ = 1/2−‖at‖, so t↦1/2−t reverses dominance pointwise.
(mixed) the half-shift t↦t+1/2 turns {‖at‖>‖bt‖} into {‖at‖+‖bt‖<1/2}, the L¹-diamond half of
the torus checkerboard; expanding sign(‖x‖−‖y‖) (coefficients −4/π²(m²−n²) on m+n odd) along the
geodesic (m,n)=k(b,−a) leaves k odd & d odd, and Σ_{k odd}1/k² = π²/8 gives deficit 1/(2ds).
Referee: exact on all 1101 primitive pairs to 60 (`lrc14_measure_dominance_tournament_macmini_S123.py`).
No prior art found in a cursory literature search (folklore-risk acknowledged).

Consequences, typed honestly:
- The θ=1/2 dominance tournament is the **dyadic preorder** — ties exactly at equal ν₂ after
  reduction, no upsets ever. As a floor-detector it is BLIND (deep well, {1..12,26}, {1..12,14}
  share one image; the two tight families have different images). This is a *confirmation* of
  opus-S399's triage axis 3 (mean-type statistics cannot see the floor) — logged as such.
- The faithful residue is the **deficit matrix** on the pair hyperbolas ds = b²−a²: the crossing
  moduli of the dominance functor are exactly the (D,s) frame's pair-sum and pair-difference
  lattices. The tournament functor did not detect the floor; it **rederived the rung lattice**.
- **Slack parity = s parity = straddle parity class** (14D even ⟹ slack ≡ s mod 2): odd-slack
  rungs (s = 14D−1: 27, 41, 55, 69, 83 …) force MIXED-parity active pairs — strict dominance
  deficit 1/(2ds) ≠ 0; the tight floor s=14D admits BOTH-ODD balanced straddles, and both tight
  families indeed straddle on (1,13), both odd, deficit exactly 0. The 2-adic layer ("the binding
  obstruction is 2-adic, not 7-adic", phase 6) is now a *coordinate on the rung lattice itself*.

### 3b. The blocker mirror palindrome (HYP-7965, resolved)

LEM-020's involution τ↔1−τ acts on THM-1240's centered-spoke blocker data exactly:

- **Set-valued palindrome law:** p ↦ q−p bijects the nearest-integer phase choices at gap k of
  carrier c onto those at gap c−1−k with EQUAL blocker sets. (A first draft claimed same-gap tie
  choices agree — refuted by the referee at c=3, d=12: {8} vs {10}; the mirror pairs choices
  ACROSS gaps, not within. Lesson kept in-file.)
- **Central column law:** for odd c, odd d, the central gap's spoke phase is exactly τ=1/2 and
  the blocker set is exactly the EVEN pack speeds — LEM-020's depth-1 2-adic layer sitting
  inside Wall B's machinery. For even d the central tie-pair is mirror-self and harmless.
- **Parity localization:** with any mirror-equivariant gauge, any gap predicate satisfies
  #{k : Q(k)} ≡ Q(central) mod 2 (c odd), ≡ 0 (c even). An odd total of odd-cycle gaps is
  FORCED onto the central gap. Verified across 8 families × all carriers, zero failures
  (`lrc14_blocker_mirror_palindrome_macmini_S123.py`).
- Empirical bonus: among all tested families only the **deep well** has fully-saturated small
  carriers (c=1,2,3,4 all gaps full-blocker; its two odd-cycle gaps at c=4 pair exactly as the
  law demands) — the covering-min extremal is extremal for blocker saturation too.

This answers opus-S399 §5's shaped question ("does the involution pairing typecheck on blocker
cycles?") — YES, and the fixed locus is computable: τ=1/2, even speeds, depth 1. The named road:
covering kills LEM-020 depths 1–3 but NOT depth 4; the palindrome now says **mirror-parity
obligations of the whole comb program localize onto the 2-adic tower**, so a Wall-B parity
argument needs exactly the depth-≥4 layers that covering leaves alive.

## 4. Typed next moves (extending the S399 practice: every residual named A / B / A-half)

1. **[Wall A] Effectivize δ via the rung-realizability map.** The metagraph-precise program:
   complete the "flip spectrum of the AP pole" — for each rung (D, s=14D−slack) in (1/14,3/41],
   either realize it or kill it (cage congruences HYP-7920 + G-K finiteness THM-1289 + gate-tower
   witnesses THM-1286 + slack-parity/straddle-parity from §3a as a new necessary condition:
   odd-slack realizers need a mixed-parity straddle, i.e. cannot be parity-homogeneous — one new
   clean kill condition per rung). Assembly line exists (death-star + kps).
2. **[Wall B] Quotient the comb program by the mirror.** All blocker enumeration halves; parity
   invariants localize to the central column = even speeds; the germ-lift obligation should be
   tested for mirror-oddness — if any lift-debt count is odd in total, it must be paid at τ=1/2,
   where the blockers are even speeds with their own 2-adic stratification (LEM-020 depth ≥ 4
   live). One codex-shaped session: rerun THM-1248's address compression on the palindrome
   quotient.
3. **[A-half] The n=12 sporadic banks** inherit the same slack-parity coordinate (s odd there
   too on the relevant rungs); check whether the H5/H6 face census can use the mixed-straddle
   necessary condition.
4. **[Engineering, dual mandate]** The deficit law is a shippable pairwise-fairness primitive:
   "which of two periodic processes is farther from sync more often" has a closed form, ties
   exactly on the dyadic classes — a natural addition to the tournament_tda / mod_rank library
   arc, and independently a clean OEIS/short-note candidate.

## 5. Scope honesty

Nothing here claims progress on the analytic content of Wall A (the effective isolation) or
Wall B (the germ lift). The session's two results are instruments: one refutation-with-structure
(the dominance functor is blind at threshold 1/2 but rederives the rung lattice and yields the
slack-parity coordinate) and one exact symmetry (the palindrome + central column + parity
localization). Both are proved at referee level (exact arithmetic, frozen outputs); the deficit
law's Fourier proof is written at sketch level in the script docstring and deserves a short
formal write-up before any THM claim. The convergence narrative (§1–2) is sourced to today's
canon and the three prior syntheses, not re-derived.

**Files:** `04-computation/lrc14_measure_dominance_tournament_macmini_S123.py` (+.out),
`04-computation/lrc14_blocker_mirror_palindrome_macmini_S123.py` (+.out), HYP-7960, HYP-7965.
**Cross-links:** THM-1289, HYP-7920/7921/7922/7930, THM-1240/1248/1250, LEM-020, THM-640,
THM-633/1230/1235/1269/1286, boxeph-S110 (MSG-1223), opus-S399 §5/§7.3/§8.6, MISTAKE-190.
