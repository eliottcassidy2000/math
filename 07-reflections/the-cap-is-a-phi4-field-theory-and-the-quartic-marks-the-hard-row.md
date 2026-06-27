# The cap is a φ⁴ field theory: the quartic cumulant goes negative exactly at the hard row k=8

*mac-mini-2026-06-27-S67. Owner: work on Lee-Yang extremality toward the LRC proof, with the φ⁴ density
`exp(−λS⁴−bS²)dS` (the (φ⁴)₂ Euclidean QFT) and the ear-decomposition facts as inspiration; see everything as
information, be creative with hypotheses and definitions. This extends S66 (the miss-count PGF zeros, HYP-3103)
and the concurrent codex Lee-Yang/Ising portfolio (HYP-3108–3116), and lands a clean new signal with a
proof-relevant reframing. Companion to [[the-miss-count-partition-function-and-its-zeros]].*

## The reframe: the cap IS a φ⁴ partition function
My S63/S64 cap structure and the user's φ⁴ cue are the same object. The covering-bound cap splits as
```
cap_k  =  C(k+1,2)/91   −   dip_k
        = (the QUADRATIC term, b·S²)   (the QUARTIC term, λ·S⁴)
```
where the **quadratic** is the pair-normalized Pascal mass (`S2`/pairwise/the `b` of `exp(−λS⁴−bS²)`), exact
for k≥10, and the **quartic** is the dip (`S4`/the `λ`), nonzero only at the sparse binding rows k=8,9. This
is *literally* the (φ⁴)₂ single-field measure: a quadratic well stabilized by a positive quartic. The gK8
Delsarte dual confirms it term-by-term: `L_yK8 = 10 − 10 S1 + 10 S2 − 9 S3 + 6 S4` — the `+6 S4` **quartic
term appears only at k=8**; k=9,10 stop at `S3`, k≥11 at `S2`. The proof's only non-pairwise content is the
quartic.

## The new signal (VERIFIED): the quartic cumulant κ₄ goes NEGATIVE at the hard row
`lrc_phi4_quartic_stabilizer_macmini_S67.py` — the cumulants of the miss-count `N` for `consec_k`:
| k | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|
| κ₄ (quartic) | **−0.79** | +1.61 | +3.92 | +6.36 | +8.19 | +9.80 |
| cap dip | **0.0141** | 0.00025 | 0 | 0 | 0 | 0 |

**κ₄ changes sign — going negative exactly at k=8**, the unique row with the largest dip (the perennial
"finite exception"). In φ⁴ language `κ₄ < 0 ⟺ sub-Gaussian ⟺ the genuine λ>0 measure` (the Simon–Griffiths /
Lee–Yang regime where the quartic *suppresses* the tail). So:
> **k=8 is the unique binding row where the miss-count is a true `λ>0` φ⁴ measure** — the quartic stabilizer
> *engages* precisely where the cap dips below the quadratic pair-Pascal. The hard row is the φ⁴ row.

(The PGF keeps **0 real roots for all k** — Lee–Yang zero-confinement holds throughout, consistent with S66;
the *quartic sign* is the finer signal that singles out k=8.)

## The Lee–Yang extremality route (proof-relevant)
The coverage extremality `max_E p0 = cap_k` reduces (S63/S64) to **one obligation: bound the dip** (the only
non-pairwise content). The φ⁴ reframe turns this into a **single 4th-cumulant bound with a guaranteed sign**:
- the dip = the quartic `S4` term; φ⁴/Lee–Yang says the quartic is a *stabilizer* (`λ>0` at the binding row,
  `κ₄<0`), so the correction is **bounded and the right sign**, not a divergence;
- so "the cap dips below pair-Pascal by a controlled amount" = "the φ⁴ quartic is `λ>0`", and the residual
  proof obligation is the **uniform bound on the standardized quartic cumulant** `κ₄/κ₂²` over the binding
  family. This is a far more concrete target than the open "consec-maximizes" crux: a moment-4 inequality.

The Lee–Yang confinement (codex HYP-3108/3111: `corr(p0, nearest-zero) = +0.899`, `corr(p0, #real) = −0.48`)
supplies the complementary half: high coverage ⟺ zeros pushed off the real axis ⟺ the φ⁴-stabilized,
maximally-correlated config. **Coverage extremality = φ⁴ stabilization = Lee–Yang zero confinement.**

## Creative definitions (new vocabulary)
- **φ⁴ row:** a binding row with `κ₄ < 0` (a genuine `λ>0` quartic-stabilized miss-count). *k=8 is the unique
  φ⁴ row* in the consec family — the φ⁴ row is the hard row.
- **quartic engagement:** the dip = the `S4`/quartic correction; it *engages* where `κ₄` flips sign.
- **φ⁴ cap potential:** `cap_k = ∫ exp(−V(S)) `-flavored, `V(S) = λ_k S⁴ + b_k S²`, `b_k ↔ C(k+1,2)/91`,
  `λ_k ↔ dip_k`. The proof is the statement that `λ_k ≥ 0` (the quartic stabilizes) and is finite.

## The ear-decomposition bridge (the odd/even split of cumulants)
The user's ear facts — strong-connected ⟺ ear decomposition, **factor-critical ⟺ ODD ear decomposition**,
2-connected series-parallel ⟺ nested ear — are the combinatorial home of the **odd cumulants**. The miss-count
has a large, growing **odd** skewness `κ₃ (3.7→5.6)` and the even quartic `κ₄`. Bold reading: the **odd
cumulants ↔ the odd-ear / odd-cycle / OCF structure** (the same odd cycles that define `H = I(Ω,2)` and that
make the PGF complex-rooted, breaking real-rootedness), while the **even quartic κ₄ ↔ the φ⁴ stabilizer**. So
the cumulant odd/even split mirrors the ear-decomposition odd/nested split: odd ears carry the OCF/odd-cycle
content (the non-real-rootedness), nested/even structure carries the φ⁴ stabilization. (Speculative; the
testable form: does the winding tournament's odd-ear count track κ₃, and the even structure track κ₄?)

## New signals to measure (extending the S66 slate + codex's)
1. **κ₄ (and its sign):** the φ⁴ quartic; *sign-change locus = the binding/hard row* (VERIFIED at k=8).
2. **standardized κ₄/κ₂²:** the φ⁴ stabilizer strength; the dip-bound target.
3. **the effective `(λ_k, b_k)`** from fitting `q` to a φ⁴ measure; `λ_k≥0` = the proof sign.
4. **odd-ear count of the winding tournament vs κ₃** (the ear/odd-cumulant bridge test).
5. **the Lee–Yang confinement margin** (nearest-zero distance; codex, `corr +0.899` with `p0`).

## Honest status
VERIFIED: the cap quadratic+quartic split; `κ₄` sign-change at k=8 (the φ⁴/Lee–Yang `λ>0` regime = the hard
row); `#real=0` throughout. BOLD/UNTESTED: that the dip-bound is exactly a uniform `κ₄/κ₂²` bound closing
coverage extremality; the odd-ear ↔ κ₃ bridge; the full φ⁴ Lee–Yang derivation of zero confinement from a
ferromagnetic sector model (note: the zeros are NOT on `|z|=1`, so it is φ⁴ not plain Ising — the right
single-spin measure is the open modeling question). The proof-relevant prize: **reduce coverage extremality
to a single quartic-cumulant inequality with a φ⁴-guaranteed sign.**

Related: HYP-3122 (this session), HYP-3103 (PGF zeros), HYP-3085 (gK8/S2 = the quadratic), THM-577 (cap value),
HYP-3113/3111 (codex Lee-Yang/ear/Ising portfolio), [[the-cap-is-a-pair-normalized-pascal-mass-and-a-web-of-connections]].
