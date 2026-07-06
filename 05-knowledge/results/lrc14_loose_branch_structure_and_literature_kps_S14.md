# The loose branch: its precise open structure, the consecutive-rung lemma, and the literature verdict

**kind-pasteur-2026-07-05-S14 (HYP-4167).** After the fleet converged (klein HYP-4151
rigidity, mac-mini HYP-4152 AP-base + descent, kps HYP-4157 peeling/critical, opus HYP-4156
"one rigidity three names"), this note pins down (a) an elementary rung lemma that
generalizes the descent's engine, (b) the *exact* two-piece open structure of the loose
branch, and (c) what the LRC(13) literature does and does not give — checked directly.

## (1) The consecutive-rung lemma (elementary, proved, no citation)

> **LEMMA.** If `{1,2,…,k} ⊆ B` (with `k ≤ 11`) and `M(B) < 2/25`, then `(k+1)` divides
> some element of `B`.

*Proof.* Evaluate at `t = 1/(k+1)`. Each `i ∈ {1,…,k}` has `‖i/(k+1)‖ = min(i,k+1−i)/(k+1)
≥ 1/(k+1)`. Each other element `X` has `‖X/(k+1)‖ ≥ 1/(k+1)` **unless** `(k+1)∣X` (in which
case it is `0`) — because `2(k+1)/25 < 1` for `k ≤ 11`, so the only residue within `2/25` of
an integer is `0`. If no element were divisible by `(k+1)`, every runner would be
`≥ 1/(k+1) ≥ 1/12 > 2/25` far at `t = 1/(k+1)`, giving `M(B) ≥ 1/12 > 2/25` — contradiction. ∎

This is the `c=1` core of mac-mini's AP-base mechanism (HYP-4152), isolated and generalized
to every `k` in one line, and it needs *no* LRC citation. The dilated version
(`c·{1,…,k} ⊆ B`) is mac-mini's CRT free-rider (their `c` witnesses `j ≡ 1 mod (k+1)` spread
the extra runner's residues so `< c` fall in the danger band). Verified: for `k = 7,…,11`,
`{1,…,k} ∪ {extras}` is **never** in the gap (0 hits over thousands of configs each).

**Why the descent does not, by itself, close the crux.** Applying rungs `k = 11, 10, …`
shows a gap-violator contains no *consecutive-from-1* run `{1,…,k}` unless it also carries a
multiple of `k+1` that extends it. But a gap-violator need not contain `{1,…,k}` at all
(it can be a dilated / "AP-free" tuple), so the descent constrains but does not empty the
gap. It is a genuine engine, not a termination.

## (2) The loose branch splits into exactly two open pieces

The dichotomy `M(B) = 1/13` (tight, `B = c·{1..12}`) **or** `M(B) ≥ 2/25` (loose) is
*equivalent* to the conjunction of:

- **(U) Extremizer uniqueness:** `M(B) = 1/13 ⟹ B` is a dilated AP `c·{1,…,12}`.
  (This is LRC(13) *equality-case* rigidity — the unique minimizer.)
- **(G) The second-value gap:** `M(B) > 1/13 ⟹ M(B) ≥ 2/25`.
  (Equivalently `σ₂ = 2/25`; the empty spectral gap `(1/13, 2/25)`.)

Both are needed; neither is trivial. (U) is the equality case of the *proven* LRC(13); (G)
is the *second* value, a strictly harder object than the conjecture's bound.

## (3) Literature verdict (checked directly, July 2026)

**Sungkawichai–Trakulthongchai, "Eleven, twelve, and thirteen lonely runners"
(arXiv:2604.23906)** — the cited LRC(≤13):
- Their proof is **computational** (sieving + a polynomial method), verifying the *bound*
  `M ≥ 1/(k+1)` for `k ≤ 12`. It does **NOT** characterize tight/extremal configurations,
  does **NOT** prove AP-uniqueness (U), and says **nothing** about the second value (G).
  → **The loose branch is genuinely open; it cannot be discharged by citation.** (Good to
  know: the owner's "cite LRC(≤13)" gives the bound only; (U)+(G) are the fleet's to prove.)
- Two usable pieces do appear: (i) their **continuity lemma** — *"if `v` is not tight, then
  there is `T` and `ε>0` such that every `t ∈ (T,T+ε)` is a witness time"* — is exactly
  opus-S87's "extremality persisting over an **interval**" frame, now with a literature
  anchor; (ii) their **polynomial method** proves any `(u_1,…,u_k) ≡ (1,…,k) mod p` with
  `gcd=1` proper when `k+1` and `p > k²+k` are odd primes — a direct tool for the
  **AP-residue** side of the rigidity (klein's "residues form the AP mod p").

**Jain–Kravitz, "Relative Lonely Runner spectra" (arXiv:2411.12684)** — a lead worth mining:
they show the *relative* Lonely-Runner spectrum (over subtori) has **"very rigid arithmetic
structure"** and admits **"explicit finite characterization"** for 2-dimensional subtori.
This is the closest existing machinery to the second-value / gap rigidity (G). If the
12-runner max-min spectrum near the floor can be cast as (or bounded by) a relative spectrum
over a low-dimensional subtorus, their finite-characterization theorem may hand the gap
directly. **Recommended fleet action:** read Jain–Kravitz and test whether the gap
`(1/13, 2/25)` is a relative-spectrum void in their sense.

## Status

- **Proved:** the consecutive-rung lemma (elementary).
- **Clarified:** loose branch = (U) extremizer uniqueness + (G) second-value gap, both open,
  both *not* citable from S-T (confirmed).
- **Leads:** S-T continuity (interval frame) + polynomial method (AP-residue); Jain–Kravitz
  relative spectra (rigid finite structure — the strongest external lead for (G)).

Connects to klein HYP-4151, mac-mini HYP-4152, kps HYP-4157 (peeling/critical), opus HYP-4156.
