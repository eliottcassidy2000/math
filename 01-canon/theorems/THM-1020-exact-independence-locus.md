---
id: THM-1020
title: THE EXACT-INDEPENDENCE LOCUS FOR PAIRS — ρ(a,b) equals the independence value (2λ)² EXACTLY iff M ∣ 2a or M ∣ 2b, where M = g/λ (M = 13g in the 1/13-convention, 14g in the 1/14-convention) and g = gcd(a,b); the mechanism is the REFLECTION SYMMETRY fold(r) = fold(M−r) of the folded identity THM-965, which makes the two folds cancel precisely when S ≡ ±D (mod M); verified 3000/3000 exact — and it DISSOLVES the repo's recorded "Sidon-far exact hit at (77,143,169)" from a numerical coincidence into a structural fact
status: PROVED (from THM-965's folded identity + the fold's r ↔ M−r symmetry); verified 3000/3000 on random pairs; the (77,143,169) triple checked pair-by-pair against the law
source: opus-2026-07-17-S356 (owner: use the wall bound on the dense core; go back through the 169 work as a depth-first search for inspiration)
depends_on: [THM-965 (the two-variable folded identity — the deviation IS the two folds), THM-1012 (the independence constant as a floor), S181 (the dilation-invariant scalar saw(S) = Σ_pairs [ρ − 4/169])]
scripts: 04-computation/dense_core_169_opus_S356.py -> 05-knowledge/results/dense_core_169_opus_S356.out
---

# THM-1020 — the exact-independence locus

**Theorem.** For positive integers a ≤ b with g = gcd(a,b) and M = g/λ,

> **ρ(a,b) = (2λ)² exactly ⟺ M ∣ 2a or M ∣ 2b.**

At λ = 1/13 the independence value is 4/169 and M = 13g; at λ = 1/14 it is
1/49 and M = 14g.

*Proof.* THM-965 writes the pair overlap as the independence value plus
fold(S mod M) − fold(D mod M), where S = a+b, D = b−a and fold(r) = r(M−r).
Since fold is symmetric under the reflection r ↔ M−r, the two folds cancel
exactly when S ≡ D or S ≡ −D (mod M) — i.e. when M ∣ (S−D) = 2a or
M ∣ (S+D) = 2b. ∎ (Verified 3000/3000 on random pairs.)

**What it dissolves.** The repo has long recorded (77, 143, 169) as a
"Sidon-far triple with exact hit 1.000, deviation +0.000000" — a striking
numerical coincidence, all three pairs landing exactly on independence
despite carrying three DIFFERENT gcds (11, 1, 13). It is not a property of
that triple. Every pair satisfies M ∣ 2b:

| pair | g | M = 13g | 2b mod M |
|---|---|---|---|
| (77, 143) | 11 | 143 | 0 |
| (77, 169) | 1 | 13 | 0 |
| (143, 169) | 13 | 169 | 0 |

77 = 7·11, 143 = 11·13, 169 = 13² all live inside the 1001 = 7·11·13
system, which is exactly why each pair's b is a multiple of its own 13g.
The coincidence is the fold reflection, seen through the 7-11-13 lattice.

**Corollary (the dilation-invariant scalar).** S181's
saw(S) = Σ_pairs [ρ(a,b) − 4/169] vanishes identically on any set all of
whose pairs lie on this locus — saw({77,143,169}) = 0 exactly. Combined
with THM-1012, saw is now bounded on both sides: below by
−2λ·Σ_pairs (a/b), and it vanishes precisely on the locus above.
