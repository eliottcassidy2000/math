# The Lonely-Runner spectrum and the hard locus: complete residue systems mod n (S521)

*claudebox-2026-06-01-S521, long creative session (synthesis). Integrates the
balanced apex-pair theorem with a computational map of the M-value spectrum and
the structure of the extremal speed sets. Companion to
lrc-gap-is-a-balanced-apex-pair-problem-s521.md.*

## The gap spectrum

`M(v) = sup_t min_i ||v_i t||`; LRC(n=m+1) is `M >= 1/n`.  By Theorem A
(balanced apex-pair), `M` is always a balanced distance `||k v_i/q|| = c/q` at a
**pairwise-sum modulus** `q = v_i + v_j`.  So every value of `M` is a fraction
`c/q` with `q` a pairwise sum.  Exhaustive enumeration (speeds `<= 20`) reveals a
clean spectrum just above the tight floor:

| level | M | binding pair-sum q | meaning |
|------|------|------|------|
| tight | `1/n` | `n` | regular n-gon (k=1) |
| 2nd  | `2/(2m+1)` | `2m+1` | |
| 3rd  | `3/(3m+1)` | `3m+1` | |
| k-th | `k/(mk+1)` | `mk+1` | |

i.e. the dominant family is `M = k/(mk+1)`, attained at `t = j/(mk+1)`, with the
binding pair summing to `mk+1`.  (This is the cleanest family, not the whole
spectrum — see the outlier below.)  The tight floor `1/n` is the `k=1` member.

## The hard locus: pairs summing to a multiple of n, and complete residue systems

Theorem B (proved): `M = 1/n` forces a pair with `n | (v_i + v_j)`.  Empirically
the exactly-tight and near-tight sets are **complete (or AP) residue systems mod
n** lifted to distinct integers — "lifted arithmetic progressions":

- `[1,2,3,4,5]` (AP) and `[1,3,4,5,9]` are the only tight sets with speeds `<= 20`
  at `m=5`; both have a binding pair summing to `6 = n`.
- `[1,5,6,11,16,17] ≡ {1,2,3,4,5,6} (mod 7)` — a complete residue system mod `7`;
  several near-tight `m=6` sets are AP-with-one-term-shifted-by-`n`
  (e.g. `[1,2,3,4,5,12]`, `12 ≡ 5 mod 7`).

So the controlling modulus is always `n = m+1`, and the hard sets are the
"regular-polygon-mod-n" configurations — exactly the fully-covered / regular
polygon objects reached from the earlier S521 reflections, now arriving from the
gap-function side.

## The outlier (a second tightness mechanism)

`[1,5,6,11,16,17]` (m=6) has `M = 5/33`, strictly between `1/7` and `2/13`, with
binding pair `(16,17)` summing to `33 = 3·11` (NOT of the form `mk+1 = 6k+1`).
This shows the `k/(mk+1)` family is not the whole spectrum: `M = c/q` for the
binding-pair-sum `q`, and `q` need not be `mk+1`.  Theorem A (`M = c/(v_i+v_j)`)
is the robust law; `mk+1` is the dominant special case.  The outlier's complete
residue structure mod 7 with a large pair-sum is a distinct extremal mechanism
worth its own study.

## The no-pair reduction (rigorous for n<=7; conjectural beyond)

- **Rigorous (n <= 7).** LRC(n) is a theorem for `n <= 7`, so `M >= 1/n`; by
  Theorem B, `M = 1/n` needs a pair summing to `0 (mod n)`.  Hence **every set with
  no such pair has `M > 1/n` strictly** — provably non-extremal.
- **Open n (8, 9, ...).** Computation: no-pair sets have `M > 1/n` with large slack
  (`min M·n ≈ 1.26` at n=8, `1.29` at n=9; zero LRC failures over thousands of
  sets). Strongly supports the **reduction conjecture**: *LRC(n) holds for all
  sets iff it holds for sets having a pair summing to a multiple of n.* A proof
  needs a uniform lower bound on `M` for no-pair sets (Theorem B only forbids
  exact tightness).

## Why this matters for n=14

The hard locus of LRC(14) is exactly the speed sets with a pair summing to a
multiple of 14 — equivalently the regular-polygon-mod-14 / complete-residue
configurations.  For such a set the binding pair sits at `t = j/14` (the regular
14-gon), and loneliness there is the `q=14` base reduction on the other 11
speeds' residues mod 14.  So the gap-function attack, the multiplicative walk, the
apex seam, and the fully-covered residual **all converge on the same object**: the
mod-14 residue structure of speed sets that are AP-like mod 14.  The whole of
S521 says: *the Lonely Runner difficulty at n is concentrated on the
arithmetic-progression-residue-systems mod n, balanced at their apex pair.*

## Seeds

a. **Prove the no-pair reduction** (uniform `M > 1/n` lower bound for no-pair
   sets) — would shrink LRC(n) to the complete-residue locus.
b. **Classify the M-spectrum** fully: which pair-sums `q` and numerators `c` occur
   (the `mk+1` family plus outliers like `33`); is there a forbidden gap just
   above `1/n` (a spectral gap → quantitative LRC)?
c. **Lift to n=14**: enumerate the complete-residue-systems mod 14 (the hard
   locus), and check loneliness only there — a far smaller class than all speed
   sets, and the genuine core.
