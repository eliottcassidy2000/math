---
id: THM-1095
title: THE THM-965 RESIDUE-TIMES-MULTIPLICATIVE-MINIMUM TEMPLATE FAILS EXACTLY AT k=3 — the same residue class contains delta(1,2,3)*m7=61/2058>0 and delta(1,2,17)*m7=-148/5831<0, from exact rational centered moments and exact product-ordered m7 enumeration. Thus delta*m7 is not a function of the speed residues modulo 14. Richer closed forms and all possible uniform signed bounds are not ruled out.
status: PROVED EXACT for the stated template refutation — two independent exact intersection engines give delta(1,2,3)=61/2058 and delta(1,2,17)=-37/11662; exhaustive product-ordered resonance enumeration proves m7=1 and 8. The broader recommendation to abandon the Bonferroni ledger is strategic evidence, not a no-go theorem.
source: opus-2026-07-17-S372 experiment; codex-2026-07-18-S67 exact rational audit
depends_on: [THM-1092 (exact centered-moment/lattice identity), THM-1093 (residue-coset motivation), THM-965 (the k=2 formula whose template this tests), THM-1080 (m7 normalization), THM-1065 (per-family evaluation cost)]
scripts: 04-computation/dedekind_evaluation_opus_S372.py, 04-computation/dedekind_normalized_opus_S372.py, 04-computation/kfold_resonance_exact_referee_codex_S67.py -> 05-knowledge/results/
---

# THM-1095 — the residue-times-`m7` template fails exactly at `k=3`

THM-1093 identified generalized higher-dimensional Dedekind sums as a natural
evaluation target. This file tests one precise form: whether the normalized
triple term `delta*m7` is determined by the speed residues modulo `14`, as the
pair term is after the appropriate gcd normalization.

## What the ledger needs, precisely

THM-965 is a **formula**, not a computation:

> mu(D_a ∩ D_b) = [4ab + fold_M((a+b) mod M) − fold_M((b−a) mod M)] / (196ab),  M = 14g

equivalently **δ·ab is determined entirely by residues mod 14g**. That is
what allowed a uniform bound over infinitely many families — the sawtooth
floor. THM-1065 already established that per-family *evaluation* is cheap
(8191 subset measures in seconds), so an evaluation **algorithm** adds
nothing. Only a formula in the parameters helps.

## The test, its two corrections, and the exact referee

The first exploratory run was wrong twice, and both errors remain useful
warnings.  Its `k=2` control used residue classes with `7|a`, so it saw only
the forced zero locus.  It also normalized triples by `abc` rather than the
multiplicative minimum `m7`.  The corrected truncated experiment used
`delta*m7` and found a sign reversal, but a finite lattice cutoff alone could
not prove that reversal.

The exact referee removes the cutoff.  For any finite speed set `A`, it
constructs the rational danger arcs, computes every intersection in two ways
(iterated interval intersection and an independent endpoint-cell sweep), and
uses

```text
delta(A) = integral product_(a in A)(1_(D_a)-1/7)
```

to obtain the centered moment by finite subset expansion.  As controls it
recovers exactly

```text
delta(1,r)*r = 6/49,5/49,4/49,3/49 for r=15,2,3,4,
delta(a,b)*a*b = 5/49 for (a,b)=(1,2),(1,16),(15,16),
delta(a,b)*a*b = 6/49 for (a,b)=(3,5),(3,19),(17,19).
```

Thus the apparent variation in the truncated `k=2` control was cutoff error:
the tested lifts in each residue class agree exactly.

The decisive triple values are

```text
delta(1,2,3)  =  61/2058,
delta(1,2,17) = -37/11662.                                  (1)
```

The product-ordered resonance enumeration checks every full-support 7-free
integer vector with coordinate product at most `36`.  It finds

```text
m7(1,2,3)=1,   witnessed by (-1,-1,1),
m7(1,2,17)=8,  witnessed by (-1,-8,1).                       (2)
```

Because every vector of product below the displayed witness occurs in that
finite enumeration, (2) proves minimality.  It can also be seen elementarily
in the second row: a product below eight forces `|n_3|=1`, while
`|n_1 n_2|<8` gives `|n_1+2n_2|<17`, contradicting
`n_1+2n_2+17n_3=0`.

Both speed triples reduce to `(1,2,3)` modulo `14`, but (1)-(2) give

```text
delta(1,2,3)*m7  =  61/2058 > 0,
delta(1,2,17)*m7 = -148/5831 < 0.                            (3)
```

The sign reversal is therefore exact; no truncation-tail assertion is needed.

This is not an isolated numerical accident.  Exact interval arithmetic gives,
for the four values `c=3,5,9,11`,

```text
delta(1,2,c)    = (70-3c)/(686c) > 0,
delta(1,2,c+14) = (14-3(c+14))/(686(c+14)) < 0,
m7(1,2,c)       = (c-1)/2,
m7(1,2,c+14)    = (c+13)/2.
```

For the first two identities, expand the centered product as

```text
delta = mu_123 - (mu_12+mu_13+mu_23)/7 + 2/343.
```

The exact arc intersections are `mu_123=1/(7d)` and `mu_12=1/14`;
the other two pair masses are each `1/(7d)` for `d=c` and each
`3/(7d)` for `d=c+14`.  Substitution gives the displayed formulas.
The witnesses `(-1,-(d-1)/2,1)` give the stated multiplicative minima,
and product-ordered exhaustion proves minimality.  Hence four distinct
residue classes exhibit the same exact sign reversal.

There is a stronger failure of sufficiency.  The rows `(1,2,3)` and
`(1,16,17)` have the same residues **and** `m7=1`, but

```text
delta(1,2,3)   = 61/2058,
delta(1,16,17) = 643/46648.
```

So even `delta` itself is not determined by the pair
`(speed residues modulo 14, m7)`.

## Conclusion

**At `k=3`, `delta*m7` is not a function of the speed residues modulo `14`.**
The THM-965 template in that precise form does not extend.  This says nothing
against a formula involving richer invariants such as the integral
equivalence class of the resonance lattice or continued-fraction data.
Generalized Dedekind reciprocity may still supply a descent algorithm; the
point is only that the two proposed inputs, residues and `m7`, are
insufficient.

## Strategic implication (not a no-go theorem)

Four consecutive independent attempts now fail, and they fail
consistently:

| session | approach | outcome |
|---|---|---|
| THM-1070 | containment / fragmentation | valid, looseness grows ~5×/level |
| THM-1085 | absolute values on the lattice | valid, constant 2.4→19→51 with k |
| THM-1093 | coset decomposition, keeping signs | exact reindexing; one sample estimate improves on the absolute bound, but no uniform inequality follows |
| THM-1095 | residue/`m7` closed form | this precise template provably does not extend |

These results justify deprioritizing the particular slot-by-slot ledger built
from absolute values or from the proposed residue/`m7` lookup.  They do **not**
prove that every uniform signed Bonferroni inequality is impossible.  A richer
invariant, a reciprocity descent, or cancellation across several ledger slots
could still work.  What survives and remains useful is the seven-moduli sieve
(THM-1035/1040), the exact mod-7 independence criterion (THM-1075/1080),
the AP decoupling and limiting profile (THM-1080), the exact coset split
(THM-1093), and the `k=2` sawtooth floor itself.

## Reproducibility

The frozen exact audit is
`04-computation/kfold_resonance_exact_referee_codex_S67.py` with SHA-256
`bea1539313883f80af0098ab7ca14c35d29c9f6e29104084d92b33845530fd64`.
Its byte-for-byte output is
`05-knowledge/results/kfold_resonance_exact_referee_codex_S67.out` with
SHA-256
`93dc829a2c3eca980d0d9543e72441e69ccff0552253e9f9141fd9fc66710dbe`.
