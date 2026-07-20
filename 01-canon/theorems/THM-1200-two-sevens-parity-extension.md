---
id: THM-1200
title: THE TWO SEVENS, AND THE PARITY EXTENSION — this corpus contains two unrelated sevens: the TOURNAMENT seven (Paley on 7 vertices, QR₇, H7, the heptagon; 620+ mentions), where 7 is a prime ≡ 3 mod 4, and the LRC seven (2λ = 1/7, S₁ = 13/7, ĥ zeros at 7|n, the k=7 arity ceiling), where 7 is n/2 = 14/2. They are NOT the same constant. The LRC seven is a PARITY artifact: ĥ(m) = sin(2πm/n)/(πm) vanishes iff n | 2m, so zeros recur every n/2 steps for EVEN n but every n steps for odd n, and the arity ceiling k < n/2 degenerates exactly only when n is even. LRC(14)'s coincidences are exact because 14 is even. Separately, THM-1195's tent/independence coincidence is confirmed to hold for EVERY n (both equal 1/(2n)), so it explains LRC(n)'s tightness in general, not 14's in particular
status: the parity facts are exact algebra, verified at n = 12…17; the tent/independence identity is verified as an exact rational equality at n = 10, 12, 13, 14, 20. The archaeology is a survey of existing canon, and it establishes that my THM-1170 DUPLICATED prior work — HYP-2059 (opus-S557) and THM-401 (opus-S590) — which is recorded here with credit rather than presented as new
source: opus-2026-07-17-S389 (owner: see how past work in the repo has treated 7 and extend)
depends_on: [HYP-2059 (the pinch lemma, S557 — prior art for THM-1170), THM-401 (pair-sum modulus C = 2n−1, S590), THM-1170 (the duplicate), THM-1195 (withdrawn; see THM-1201), MISTAKE-173, MISTAKE-174]
scripts: 04-computation/seven_synthesis_opus_S389.py -> 05-knowledge/results/seven_synthesis_opus_S389.out
---

# THM-1200 — two sevens, and why LRC(14)'s are exact

## The archaeology

Searching the corpus for structural sevens turns up two families that have
nothing to do with each other:

| | the TOURNAMENT seven | the LRC seven |
|---|---|---|
| appears as | Paley on 7 vertices, QR₇, H7, heptagon | 2λ = 1/7, S₁ = 13/7, ĥ zeros, k=7 ceiling |
| mentions | 620+ (Paley), 230 (H7), 12 (QR₇) | 38 (13/7), 199 (mod 7) |
| what 7 *is* | a prime ≡ 3 (mod 4) | **n/2 = 14/2** |

They coincide numerically and not structurally. Nothing in the LRC seven
requires primality; nothing in the Paley seven requires evenness.

## Prior art I duplicated

**HYP-2059 (opus-2026-06-02-S557), the pinch lemma:** the loneliness radius
is attained at a pair-sum time t = m/(v_a+v_b). **THM-401 (opus-S590)** then
proved the modulus identity C = 2n−1 on top of it, via a Farey-neighbour
argument.

My **THM-1170 (S383)** rederived the same critical-point structure 800
sessions later, without finding the earlier work. It is recorded here as a
duplicate. What THM-1170 adds is the *measured* material — the 169 = 13²
candidate bound, the ~15.5% witness density, and the tight families landing
at 14 = 1+13 — not the structural claim, which is HYP-2059's.

**A near-miss.** I briefly believed I had a counterexample to THM-401: 4 of
25 families whose optimum sat at a difference-only denominator. It is an
artifact of denominator reduction — for V = {9,15,16,…,51} the optimum 1/6
appears at q = 6 (a difference), but 1/6 = 4/24 and 24 **is** a pair sum, so
pair-sums attain the identical value with deficit 0.00000000. **THM-401
stands.** See MISTAKE-173.

## The extension: the LRC seven is a parity effect

ĥ(m) = sin(2πmλ)/(πm) with λ = 1/n vanishes iff n | 2m, i.e.

> m ≡ 0 mod n/gcd(2,n)

so the zeros recur every **n/2** steps when n is even, and every **n** steps
when n is odd — twice as dense for even n. Likewise the density coefficient
1 − 2kλ vanishes at k = n/2, an integer **only** for even n:

| n | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|
| ĥ zeros every | 6 | 13 | **7** | 15 | 8 | 17 |
| n/2 | 6 | 6.5 | **7** | 7.5 | 8 | 8.5 |
| arity ceiling exact? | yes | no | **yes** | no | yes | no |

**So LRC(14)'s sevens are exact because 14 is even.** At odd n the same
quantities exist but never degenerate exactly — the arity ceiling is never
attained, and the Fourier vanishing is half as dense. That is a structural
difference between even and odd LRC, and it predicts that the machinery of
THM-1155/1165/1180 transfers to even n and weakens at odd n.

## THM-1195 generalised

**Correction (codex-S78): withdrawn.**  The arithmetic equality described in
this section is true, but one side came from the false THM-1195 cell upper
bound.  It therefore does not compare two valid estimates and does not explain
LRC tightness.  See THM-1201 and MISTAKE-174.  The parity analysis in the
preceding sections is independent of this correction.

The tent threshold λ/2 = 1/(2n) and the independence expectation
(1/2)/((n−1)+1) = 1/(2n) are **equal for every n** — verified as exact
rationals at n = 10, 12, 13, 14, 20. So the coincidence THM-1195 found is
not special to 14: it says LRC(n)'s floor sits exactly where the geometric
and probabilistic estimates cross, **for all n**. That is the general reason
the conjecture is knife-edge rather than a fact about 14.


---

**AMENDMENT (opus-2026-07-20-S401, THM-1380 §5).** The claim that the two
sevens are *structurally* unrelated is **too strong at the point `n = 14`**: the
extremal maximizer set is exactly `(ℤ/14)*` (φ(14)=6 points, 3 antipodal orbits), and
`14 = |D₇|`. What survives is the **general-`n`** statement — the LRC seven is `n/2`,
fixed by parity, and needs no primality, while the Paley seven needs `p ≡ 3 (mod 4)`
prime; at general even `n` they diverge. Read THM-1200 as about the family, not the
point. Credit: kind-pasteur-2026-06-28-S31av. Not withdrawn — scoped.
