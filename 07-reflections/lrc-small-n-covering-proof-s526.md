---
source: oracle-2026-06-01-S526
status: RESULT — LRC(n=3) proved by the new methodology; n>=4 reduced to an explicit character sum
tags: [LRC, covering, fourier, roots-of-unity, legendre, n3-proof, permutohedron, A000016]
---

# Proving small LRC with only the new methodology: n=3 closes; n>=4 is one character sum away

**Prompt (user):** try proving cases smaller than n=14 using only the new
methodology.

**Result.** The covering reformulation (S525), evaluated with the roots-of-unity
characters that underlie this whole thread (S522/S523), gives a **complete,
elementary proof of LRC(n=3)** and reduces every larger n to bounding one
explicit character/resonance sum. n=2 is trivial; **n=3 is proved**; n=4..7 are
verified and reduced but not closed — the residual is, provably, the same wall as
n=14.

## The methodology in one identity

Covering (S525): `B_s = { t in [0,1) : ||s t|| < 1/n }` is open of measure `2/n`;
the safe set `SAFE = [0,1) \ U_i B_{s_i}`, and **LRC(n) for a speed system holds
iff SAFE != empty**. The safe indicator `g(x)=1[||x||>=1/n]` has Fourier
coefficients `g_0 = 1-2/n`, `g_k = -sin(2pi k/n)/(pi k)` (k!=0). By orthogonality,

```
|SAFE| = |∩_i S_{s_i}| = Σ_{ (k_1..k_m): Σ k_i s_i = 0 }  Π_i g_{k_i}
       = (1 - 2/n)^m  +  (resonance corrections),     m = n-1.
```

The main term `(1-2/n)^{n-1}` is exactly opus-S524's "independence" value
(`(6/7)^13 = 0.1348` at n=14). LRC(n) is the statement that the corrections never
drive `|SAFE|` below 0 (or, at tightness, leave a boundary point). The characters
`g_k` are the n-th roots of unity sieved by `sin(2pi k/n)` — the same roots-of-
unity substrate as the round-tournament necklace (S523).

## n=3 — a complete proof

Two runners, speeds `a<b`, `gcd(a,b)=1`, threshold `1/3`. Here `m=2`, so the only
resonances are the **2-term** ones `k a + j b = 0`; with `gcd(a,b)=1` they are
`k=bM, j=-aM`. The n=3 character is `sin(2pi k/3) = (sqrt3/2) chi(k)` where
`chi = ( . /3)` is the Legendre symbol mod 3 (`chi(1)=+1, chi(2)=-1, chi(0)=0`).
Summing (`Σ_{M!=0, 3∤M} 1/M^2 = 8 pi^2/27`) gives the **closed form**

```
|B_a ∩ B_b| = 4/9 + (2/9) * chi(a)chi(b)/(a b),     equivalently
|SAFE|       = 1/9 + (2/9) * chi(a)chi(b)/(a b).
```

(Verified against numeric integration for all coprime pairs up to 12.) Now
`chi(a)chi(b) in {-1,0,+1}` and `ab >= 2`, so

```
|SAFE| = 1/9 + (2/9) chi(a)chi(b)/(ab) >= 1/9 - (2/9)/(ab) >= 1/9 - 1/9 = 0,
```

with equality iff `ab = 2 and chi(a)chi(b) = -1`, i.e. iff `{a,b} = {1,2}`.

- If `{a,b} != {1,2}`: `|SAFE| > 0`, so the closed safe set has positive measure,
  hence is nonempty — a lonely time exists. 
- If `{a,b} = {1,2}`: `|SAFE| = 0`, but `t = 1/3` has `||1/3|| = ||2/3|| = 1/3 >= 1/3`,
  so SAFE `= {1/3, 2/3}` is nonempty.

Either way a lonely time exists. **LRC(n=3) holds for every speed pair, with
`{1,2}` (the AP / regular-polygon extreme, S522) the unique tight case.** ∎

The proof is nothing but: covering reformulation + the mod-3 quadratic-residue
character (the A_2 / Paley-3 character of the methodology). The independence term
`1/9` is *exactly* compensated by the resonance correction, never over-compensated.

## Why n=3 closes and n>=4 does not (yet)

n=3 is the case `m=2`: the safe set is a **pairwise** intersection, so its measure
is a single 2-term resonance sum, which collapses to one Legendre character and is
bounded in closed form. For `n>=4` (`m>=3`) `|SAFE|` is an `m`-fold intersection
whose resonance sum contains **3-term and higher** resonances `k_a a + k_b b +
k_c c + ... = 0`; these do not factor into a single character and are exactly the
"multi-way correlations" opus-S524 could not bound. The main term stays positive
(`-> e^{-2} ~ 0.135`), so a proof needs only a *bound* on the higher-resonance
correction — but that bound is, at n=14, LRC@14 itself.

## n=4..7: verified and located (`lrc_small_n_covering_proof_s526.py`)

Exact scan of primitive speed sets (speeds bounded): **0 LRC failures** for
n=4,5,6,7. The AP `{1,...,n-1}` is always a tight set. Honest caveat: the AP is
**not the unique tight set** for n=5,6 (`{1,3,4,7}` at n=5; `{1,3,4,5,9}` at n=6
are also tight). So the S525 sub-target "wall-only => AP" must be refined: there
are several margin-`1/n` sets, though the AP remains the canonical one.

## Scorecard (only this methodology)

| n | runners | status by this methodology |
|---|---------|----------------------------|
| 2 | 1 | trivial: `t=1/(2s)` gives `||1/2||=1/2`. PROVED. |
| 3 | 2 | **PROVED** (covering + Legendre-mod-3 character; closed form above). |
| 4 | 3 | reduced to bounding the 3-term resonance correction; 0 failures in scan. OPEN. |
| 5,6,7 | 4,5,6 | reduced + verified by scan; multi-way resonance correction OPEN. |
| 14 | 13 | same residual, larger; localized to the AP wall (S525). OPEN. |

So the methodology genuinely *proves* the bottom of the ladder (n<=3) and turns
every larger case into one explicit, finite-to-state object — the higher-resonance
character sum `Σ_{Σ k_i s_i=0, |k|>=3 support} Π g_{k_i}` — whose sign control is
the whole game. The next honest target: bound the **3-term** resonance sum (n=4),
the first case the pairwise n=3 trick cannot reach.

## Anchor
`04-computation/lrc_small_n_covering_proof_s526.py` (+ `.out`): n=3 closed form
verified; `(1-2/n)^{n-1}` = opus independence term; n=4..7 scan (0 failures, AP
tight, tight-set list). Builds on S525 (covering), HYP-1998/2000 (round/block),
S522 (roots of unity). The closed form `|B_a∩B_b| = 4/9 + (2/9)chi(a)chi(b)/(ab)`
is THM-worthy (proves LRC n=3).
