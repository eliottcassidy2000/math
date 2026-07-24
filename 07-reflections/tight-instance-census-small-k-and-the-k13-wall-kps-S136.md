---
source: kind-pasteur-2026-07-23-S136 (Opus 4.8)
status: CENSUS + method validation. Worked the S135 continuations (modulus-K rigidity; completeness of the
  k=13 tight locus) and, when those saturated, the emergent move: EXHAUSTIVE tight-instance classification
  for small k. Produces a census table, a validated multi-lift search method, an independent rediscovery of
  the Goddyn-Wong construction, and an honest bound on what my k=13 search can and cannot claim.
  Directly targets opus-S4's "the sole wall is the tight locus" (OPEN-Q-108).
tags: [lrc, lonely-runner, tight-instances, census, classification, goddyn-wong, modulus-rigidity, exhaustive]
related: [kps-S135, kps-S134, opus-S4 (OPEN-Q-108), THM-518, Cusick, Goddyn-Wong]
---

# Tight-instance census for small k, and what the k=13 wall actually looks like

opus-S4 established the endgame shape: the Fejér/variational route certifies `gap > 1/14` for everything with
`gap > 1/14` at practical degree, so **the sole wall is the tight locus** (configs with `gap` exactly `1/(k+1)`).
This session characterises that wall.

## 1. Modulus-K rigidity lemma — survives every test
> **Conjecture (modulus-K rigidity).** Every tight instance is witnessed at `τ = a/K`, `gcd(a,K)=1`;
> equivalently **no speed of a tight instance is divisible by `K`**.

Tested for `k=13, K=14`: all single-replacements of `{1,…,13}` by multiples of 14 (14,28,…,98) and all 3510
two-replacements containing 14 → **zero tight instances**. Also verified directly: both known k=13 tight
instances admit witnesses only at `m=1` (i.e. `q=14`), with `a ∈ {1,3,5,9,11,13}` = the units mod 14, never at
`q=14m, m≥2`. Status: strongly evidenced, unproved. (Proof sketch target: a speed `≡0 mod K` kills `τ=a/K`, and
the `q=Km` witnesses appear to always admit a strictly better `τ`.)

## 2. The k=13 tight locus: two instances, across ~12 independent searches
```
T1 = {1,…,13}                    T2 = {1,…,11,13,24}      (12 accelerated by 2)
```
Searches all returning ONLY these: single-replacement from T1/T2 (values ≤300); 2-replacement (<40); BFS depth 2
(≤90); `{1..12}∪{W}`, W ≤ 400; `{1..11}∪{A,B}`, ≤120; `{1..10}∪{A,B,C}` sampled; **all subset-accelerations by
factors 2,3,4 and mixed powers of 2**; **acceleration axes scanned to W = 2200** (past `2^11=2048`, only W=24
survives); residue-pattern + multi-lift (≤4 lifts); global simulated annealing (speeds ≤40, never even reaches
the locus). **T2 is exactly the Goddyn–Wong construction** ("accelerate a speed slightly less than n by a
suitable integer factor"; their k=7 example `{1,2,3,4,5,7,12}` is `{1..7}` with `6→12`) — independently
rediscovered here before I found the citation.

## 3. NEW — exhaustive census for small k (the emergent result)
| k | K=k+1 | factorisation | # tight (up to dilation, within range) | non-canonical members |
|---|---|---|---|---|
| 3 | 4 | 2² | 1 | — |
| 4 | 5 | prime | **2** | `{1,3,4,7}` |
| 5 | 6 | 2·3 | **2** | `{1,3,4,5,9}` |
| 6 | 7 | prime | 1 | — |
| 7 | 8 | 2³ | **3** | `{1,2,3,4,5,7,12}`, `{1,4,5,6,7,11,13}` |
| 8 | 9 | 3² | 1 | — |
| 9 | 10 | 2·5 | 1 | — |
| 13 | 14 | 2·7 | **2** (this work) | `{1,…,11,13,24}` |

**The count is NOT monotone in k and is NOT explained by primality alone** (K=5 prime → 2, but K=7 prime → 1).
It is governed by the arithmetic of `K`, in a law not yet identified. *Caveat: counts are within the searched
speed range (TOP = 14–26); larger-speed instances may exist, so these are lower bounds.*

## 4. The k=7 exotic instance — and the bias it exposed in my k=13 search
`{1,4,5,6,7,11,13}` (k=7) is **not** an acceleration of the canonical instance: its residues mod 8 are
`{1,3,4,5,5,6,7}` (miss 2, duplicate 5) and it requires **two** elements lifted by `K` (3→11 and 5→13).
All my earlier k=13 searches lifted at most **one** element — a genuine blind spot. I built a multi-lift search,
**validated it as a control** (it recovers both non-canonical k=7 instances), and applied it to k=13 with up to
4 lifts over all miss-one/duplicate-one residue patterns (68 640 configs): **no new k=13 tight instance**.

## 5. Honest status
- k=13 tight locus = `{T1, T2}` is **strongly evidenced but NOT proved complete.** The searches are extensive and
  now include the previously-missing exotic direction, but remain biased toward the canonical neighbourhood, and
  the literature's `2^{n-2}` "barrier" family was never located for k=13 (it must use multiple large speeds).
- If `{T1,T2}` *is* complete, the LRC(14) endgame residual is **two explicit configurations** — a remarkably
  small wall, and both have `gap = 1/14` verifiable exactly. That is the encouraging reading of opus-S4's
  OPEN-Q-108.

## 6. Next
1. **Prove modulus-K rigidity** (§1) — self-contained, would cut the search space enormously.
2. **Find the arithmetic law** behind the §3 census (why K=8 → 3 but K=9 → 1?). This predicts the k=13 count and
   would either confirm `{T1,T2}` or tell us exactly what is missing. *Most interesting open thread.*
3. Locate the `2^{n-2}` barrier family explicitly (multiple large speeds) and test it at k=13.

Files: `/tmp/{lemma14,bigfam,anneal,accel,bigW,small_n,residue,multilift,multilift3,k7full,k8,table}.py`.
