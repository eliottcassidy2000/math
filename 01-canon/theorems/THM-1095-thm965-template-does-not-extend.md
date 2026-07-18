---
id: THM-1095
title: THE THM-965 TEMPLATE DOES NOT EXTEND TO k=3 — δ·m₇ is NOT a function of the residues (a,b,c) mod 14: within a single residue class it varies by 150–195% AND CHANGES SIGN (+0.029115 at (1,2,3) against −0.026993 at (1,2,17)), which truncation cannot produce, while the k=2 control varies only 0.8–8.6% consistent with exact determination plus tail error. So the 3-body term is not expressible as (residue data) / (multiplicative minimum), the shape that made THM-965 a FORMULA rather than a computation; higher-dimensional Dedekind sums carry reciprocity laws giving an ALGORITHM, but per-family evaluation is already cheap (THM-1065), so an algorithm buys nothing for uniformity. Combined with THM-1070, THM-1085 and THM-1090, the uniform Bonferroni-ledger route is blocked and should be abandoned
status: measured; the sign reversal within a residue class is robust (an order of magnitude, at small speeds where truncation is least) and the k=2 control behaves as the known closed form predicts. This refutes the THM-965 TEMPLATE (residues × multiplicative minimum); it does NOT prove that no closed form of any shape exists
source: opus-2026-07-17-S372 (owner: work the higher-dimensional Dedekind sum evaluation)
depends_on: [THM-1090 (which named evaluation as the only remaining route), THM-965 (the k=2 formula whose template this tests), THM-1080 (m₇, the correct normalisation), THM-1065 (per-family evaluation is already cheap — why an algorithm does not help)]
scripts: 04-computation/dedekind_evaluation_opus_S372.py, dedekind_normalized_opus_S372.py -> 05-knowledge/results/
---

# THM-1095 — the evaluation route, tested and blocked

THM-1090 concluded that δ(S) must be evaluated rather than estimated, and
named higher-dimensional Dedekind sums. This file tests whether that
evaluation can take the form the ledger actually needs.

## What the ledger needs, precisely

THM-965 is a **formula**, not a computation:

> mu(D_a ∩ D_b) = [4ab + fold_M((a+b) mod M) − fold_M((b−a) mod M)] / (196ab),  M = 14g

equivalently **δ·ab is determined entirely by residues mod 14g**. That is
what allowed a uniform bound over infinitely many families — the sawtooth
floor. THM-1065 already established that per-family *evaluation* is cheap
(8191 subset measures in seconds), so an evaluation **algorithm** adds
nothing. Only a formula in the parameters helps.

## The test, and a correction to my first attempt

My first run was wrong twice and I record both. The k=2 control was
**vacuous** — every residue class it displayed had a ≡ 0 mod 14, hence
7 | a, hence δ = 0 by the S368 criterion, so it showed only trivial
zeros. And the normalisation was wrong: I used δ·abc, but at k=2 the
natural factor ab **is** m₇, so the correct analogue is δ·m₇ (THM-1080:
δ ~ 1/(π^k m₇)). Corrected:

**Control (k=2, 7-free residues), δ·ab within a residue class:**

| residues | relative spread |
|---|---|
| (1,1) | 8.62% |
| (1,2) | 3.71% |
| (1,3) | 2.62% |
| (1,4) | 0.80% |

Small and consistent with exact determination plus truncation tail — the
closed form predicts exactly 0, and the residual is the N = 300 cutoff.

**The k=3 test, δ·m₇ within a residue class:**

| residues (1,2,3) | δ | m₇ | δ·m₇ |
|---|---|---|---|
| (1,2,3) | +0.029115 | 1 | **+0.029115** |
| (1,2,17) | −0.003374 | 8 | **−0.026993** |
| (1,16,17) | +0.014355 | 1 | +0.014355 |
| (15,16,17) | +0.007614 | 2 | +0.015228 |

Relative spread **192.7%**, and across the four classes tested 153.7%,
177.2%, 193.7%. Critically the quantity **changes sign** within a single
residue class, at small speeds where truncation error is smallest. A tail
cutoff cannot flip a sign by an order of magnitude.

## Conclusion

**δ at k=3 is not a function of (residues mod 14, m₇).** The THM-965
template does not extend. Note carefully what this does and does not
show: it refutes that specific shape, not the existence of a closed form
involving richer invariants (GL-equivalence class of the lattice,
continued-fraction data). Zagier's higher-dimensional Dedekind sums do
satisfy reciprocity laws — but reciprocity yields a **descent algorithm**,
and by THM-1065 algorithms do not address the obstruction, which is
uniformity over infinitely many families, not per-family cost.

## Recommendation: abandon the ledger route

Four consecutive independent attempts now fail, and they fail
consistently:

| session | approach | outcome |
|---|---|---|
| THM-1070 | containment / fragmentation | valid, looseness grows ~5×/level |
| THM-1085 | absolute values on the lattice | valid, constant 2.4→19→51 with k |
| THM-1090 | coset decomposition, keeping signs | identical to the absolute bound |
| THM-1095 | closed-form evaluation | template provably does not extend |

The uniform Bonferroni ledger — B₇ built slot by slot as B₅'s S₂ slot was
— should be **abandoned**. What survives and remains useful: the
seven-moduli sieve (THM-1035/1040, ~89% kill rate), the exact mod-7
independence criterion (THM-1075/1080), the AP decoupling and its
limiting profile (THM-1080), and the k=2 sawtooth floor itself. Those are
tools; the ledger was a plan, and the plan does not survive contact with
k=3.
