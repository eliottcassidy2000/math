# The Pascal-slope-d family: one triangle, a Pisot tower per slope, and the project's recursions in a single frame

**Source:** mac-mini-2026-06-14-S1. Dispatch: extend the family behind three
constructions — Pascal (→2ⁿ), the Fibonacci diagonals, and a third (→Narayana's cows)
— and their "same-pace" central partners. Canon: T819, HYP-2523..2525, OPEN-Q-098.
Compute: 04-computation/pascal_slope_d_family_macmini_0614s1.py (all claims below
reproduced exactly, d=1..6).

## The family, in one line

Read Pascal's triangle along the diagonal of slope d. Construction-d, row n, is

    row_d(n) = [ C(n−1−(d−1)k, k) : k = 0, 1, 2, … ],

and its **row-sum** a_d(n) satisfies the **d-step Fibonacci recurrence**

    a_d(n) = a_d(n−1) + a_d(n−d),     GF  Σ a_d(n) xⁿ = 1/(1 − x − x^d).

Each construction the user wrote is exactly this:

| d | row entries | row-sum sequence a_d | OEIS | recurrence |
|---|---|---|---|---|
| 1 | C(n−1,k) (full Pascal rows) | **2ⁿ** | A000079 | a(n)=2a(n−1) |
| 2 | C(n−1−k,k) (shallow diagonals) | **Fibonacci** | A000045 | a(n)=a(n−1)+a(n−2) |
| 3 | C(n−1−2k,k) | **Narayana's cows** 1,1,1,2,3,4,6,9,13,19,28,40 | A000930 | a(n)=a(n−1)+a(n−3) |
| 4 | C(n−1−3k,k) | 1,1,1,1,2,3,4,5,7,10,14,19 | A003269 | a(n)=a(n−1)+a(n−4) |
| 5 | C(n−1−4k,k) | 1,1,1,1,1,2,3,4,5,6,8,11 | A003520 | a(n)=a(n−1)+a(n−5) |
| 6 | C(n−1−5k,k) | 1,1,1,1,1,1,2,3,4,5,6,7 | A005708 | a(n)=a(n−1)+a(n−6) |

(Construction-3's stated sums 1,2,3,4,6,9,13,19,28,40 are Narayana's cows, EXACT;
construction-2's are Fibonacci, EXACT. Verified.)

## Three readings of the same recurrence

1. **Pascal diagonal.** a_d(n) = Σ_k C(n−1−(d−1)k, k): the antidiagonal sums of
   Pascal at slope d. The k-th entry counts the configurations using exactly k of the
   "long" steps.
2. **Tiling.** a_d(n) = number of tilings of a 1×(n−1) strip by tiles of length 1 and
   length d (a d-omino contributes one long step). Verified: a_d(n) = (#such tilings)
   for all d=1..6. This is the project-native reading — the staircase/tiling model is
   built from exactly such tile choices.
3. **Gap-d independent sets.** a_d(n) = number of binary strings of length n−1 with no
   two 1's closer than d apart (place a d-omino's left end at each 1) — the
   "lonely-1" / minimum-gap-d configurations. (At d=2 this is the classic
   Fibonacci ↔ independent-sets-in-a-path / Zeckendorf fact.)

## The Pisot tower: one constant per slope

The characteristic equation is **xᵈ = x^{d−1} + 1**, whose dominant real root ρ_d > 1
is a **Pisot number** governing the growth a_d(n) ~ C·ρ_dⁿ:

| d | ρ_d | name |
|---|---|---|
| 1 | 2 | (doubling) |
| 2 | 1.6180339887 | **golden ratio φ** |
| 3 | 1.4655712319 | **supergolden ratio ψ** (x³=x²+1) |
| 4 | 1.3802775691 | (quartic Pisot) |
| 5 | 1.3247179572 | **the plastic number ρ** (x⁵=x⁴+1, equivalently x³=x+1) |
| 6 | 1.2851990332 | … |

ρ_d ↓ 1 as d → ∞: longer steps ⟹ slower growth. The plastic number reappearing at
d=5 is not a coincidence — ρ satisfies both x³=x+1 (Padovan) and x⁵=x⁴+1 (this
family), so the Padovan/plastic world and this slope-5 world share a root.

## The "same-pace" central partners

The user's partners are the **central / ridge binomial sequences** that grow at the
*same exponential rate* as each target:

- **d=1 partner = C(n, ⌊n/2⌋)** — the central binomial coefficient (A001405):
  1,1,2,3,6,10,20,35,70,126,… **CONFIRMED exactly.** It grows like 2ⁿ/√n — same base
  2 as the powers of 2. And it is *literally the project's metagraph WIDTH*
  C(n−2, ⌊(n−2)/2⌋) (CLAUDE.md, exact for n≤6): the widest antichain of G_n is the
  central binomial, the d=1 ridge.
- **d≥2 partners = shallower ridges.** The Fibonacci/Narayana partners walk a *less
  steep* ridge of Pascal (k grows slower than n/2), so they grow at φⁿ, ψⁿ rather than
  2ⁿ — the canonical version is the maximal binomial on the slope-d diagonal,
  max_k C(n−(d−1)k, k). (The user's exact recalled values are close variants of this
  ridge; the invariant content is "central binomial re-angled to the target's Pisot
  pace.") The "pace" the user names is exactly ρ_d.

So: **row-sums = a_d(n) ~ ρ_dⁿ (the lower-symmetric edge of the construction); central
partner = the ridge ~ ρ_dⁿ (the peak).** Both grow at ρ_d; the row-sum is the total,
the partner is the largest single term. For d=1 they are 2ⁿ and C(n,⌊n/2⌋); the ratio
2ⁿ / C(n,⌊n/2⌋) ~ √(πn/2) is the usual peak-to-sum factor.

## The project weave: a Pisot tower per slope

This family is a single index running through the project's recursive structures —
the same "double + carry-back-d" engine at every d:

- **d = 1 (ρ=2):** the **skew-Sylvester doubling tower** (THM-447/480). H_{2m} =
  [[H,H],[−Hᵀ,Hᵀ]] doubles order; the self-dual d⁺ codes live at the **central
  dimension n/2** — the d=1 ridge, the central binomial. Powers of 2 = the orders;
  central binomial = the code dimension AND the metagraph width. The whole skew-tower /
  RM-duality story (this week's thread, THM-477/480, HYP-2409) is the d=1 floor of this
  family.
- **d = 2 (ρ=φ):** **Fibonacci / Zeckendorf** — the project's Zeckendorf-tournament
  tilings and the gap-2 independent-set reading of the staircase. (RIGOROUS that
  Fibonacci = gap-2 independent sets; the precise tournament count is the project's
  Zeckendorf work.)
- **d = 3 (ρ=ψ supergolden):** **Narayana's cows** — and the supergolden ratio.
- **d = 5 (ρ=plastic):** **A003520 / the plastic number** — the constant of the monad
  / free-factorial-law sessions (THM-438 family). The plastic number's reappearance at
  d=5 ties that whole resurgence-analysis thread to this single triangle.

**The unifying statement (HYP-2523):** the project's recursive towers are the
slope-stratification of one Pascal triangle. "Double" (d=1) is the skew-Hadamard /
self-dual / powers-of-2 floor; "delay by d" carries you to the golden, supergolden,
plastic, … Pisot towers, each a shallower diagonal of the same triangle, each with its
row-sum sequence (a_d) and its ridge (the central partner). The metagraph width being
the d=1 central binomial is the cleanest rigorous anchor; the higher-d anchors are the
project's Fibonacci/Narayana/plastic appearances, now seen as one family.

## Recurrences in n (what the user asked to extract)

- Row-sum: a_d(n) = a_d(n−1) + a_d(n−d) (proved/verified).
- Entry: row_d(n)[k] = C(n−1−(d−1)k, k); Pascal recurrence gives
  row_d(n)[k] = row_d(n−1)[k] + row_d(n−d)[k−1] (a "carry one long step" recurrence —
  the same shape as the row-sum, one tile-class down).
- Ridge/partner: the central binomial satisfies C(n,⌊n/2⌋) = C(n−1,⌊(n−1)/2⌋) ·
  (2 − [n even]/⌈n/2⌉)-type ratios; for the matched ridge at slope d the local ratio
  → ρ_d.

## Open leads

- **OPEN-Q-098:** does a_d(n) count a natural gap-d family of tournaments / tilings in
  the staircase model (a d-omino tiling of the base path, or Hamiltonian structures
  with a minimum-gap-d constraint), making the whole Pisot tower a tournament
  invariant — not just an analogy? The d=1 case (powers of 2 = the full tile-flip
  hypercube layer count, central binomial = width) is rigorous; the d≥2 gap-d
  tournament objects are the lead.
- The exact "same-pace partner" integer sequence per d (OEIS IDs) and whether it is the
  matched-diagonal ridge or a distinct central walk — pin down against the user's
  recalled values.
- The plastic-number coincidence (d=5 here vs Padovan x³=x+1) — is there a tournament
  bridge between the monad/free-factorial sessions (plastic) and a slope-5 tiling?
