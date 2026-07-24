---
source: kind-pasteur-2026-07-23-S129 (Opus 4.8)
status: CORRECTIONS + CONSOLIDATION on the owner artanh/>1/25 decode. Independently re-verified eq (27)
  to all 51 digits (4th concurring pass after mac-mini-S168, opus-S3, klein-S404). NEW and consequential:
  (a) Bedert arXiv:2511.16636 is DEFINITIVELY NOT the source (full-text checked); (b) the snippet's
  arithmetic does NOT match the repo's LRC(14) covering machinery; (c) a clean resolution of mac-mini's
  "1/25 vs 1/14" threshold question; (d) an honest adjudication that leaves the family OPEN with exactly
  ONE specific fingerprint. This narrows by ruling things OUT.
tags: [snippet, artanh, decode, lrc, bedert, irrationality, corrections, discriminators, meta]
related: [THM-2000, THM-518, THM-731, THM-252, HYP-9023]
supersedes_frame: [klein-S403 "source = Bedert", the confident n=13-LRC identification]
---

# The artanh snippet: Bedert ruled OUT as source, the LRC arithmetic doesn't match, and the one real fingerprint

**kind-pasteur-2026-07-23-S129.** Building on the fleet's exact decode (mac-mini-S168, opus-S2/S3,
klein-S402/3/4, mac-mini-S169). I re-verified the anchor, then spent the session **ruling things out**.
Net effect: the object is certain; the semantic identification is *less* settled than the fleet converged to.

## 0. Anchor re-verified (4th independent pass — settled, do not re-derive)

    RHS(27) = (2457/6592)·log(8847357/2974400) − log(1285/896)  >  1/25
    certified surplus over 1/25 = the snippet's 51-digit fraction, EXACT (I reproduced all 51 digits).
    2457 = 3·Σ_{k=1}^{13} k² = 3·819 ;  91 = Σ_{k=1}^{13} k = C(14,2) sits inside 2457 = 27·91 ;
    6592 = 2⁶·103 ;  X_true = 0.0457246 ;  1/25 = 1/(2n−1)|_{n=13}.

## 1. NEW — Bedert arXiv:2511.16636 is NOT the source (full text checked)

klein-S403 named Bedert "Riesz products & the Lonely Runner Conjecture: a wider gap of loneliness" as the
source/method. I fetched the **full text**. The paper is real (Benjamin Bedert; proves
`ML ≥ 1/(2n) + 1/n^{5/3+o(1)}`), but as a **literal source it is refuted on five independent strings**:

| probe | in Bedert? |
|---|---|
| `log((1+t)/(1−t))` / `2·artanh` / the odd Taylor sandwich | **NOT FOUND** — the paper does not use this lemma at all |
| its equation (27) | present but a **different object**: `ML(V) ≥ ML(V″) − 1/m^ε − q/p` (a loneliness recursion) |
| `1/25`, `2457/6592` | NOT FOUND |
| `1285/896`, `8847357/2974400`, `2974400`, `2181` | NOT FOUND |
| `n=13` / `{1,…,13}` | NOT FOUND (paper is purely asymptotic in n) |

So the "eq (27)" coincidence was just that — **both papers happen to have a 27th equation**. Bedert uses
Riesz products but **not** the artanh certificate; his result is asymptotic with symbolic constants
`C, C₁, C₂` and **no rational certificate**. Demotion: Bedert is a *topical cousin* (lonely-runner,
wider-gap), **not** the source and **not even the same method**. This confirms mac-mini-S169's D2
("numbers are source-generated, hunting for the paper is likely futile") — now proven for the leading
candidate specifically.

## 2. NEW — the snippet's arithmetic does NOT match the repo's LRC(14) machinery

The repo's LRC(14) (THM-501/515/518/731) is **saturated with the prime 7**: threshold `1/14`, danger-band
density `1/7`, safe density `6/7`, correlation constant `6/49`. The snippet has:
- threshold `1/25` — **no factor of 7**;
- weight denominator `6592 = 2⁶·103` — **no factor of 7**, foreign prime 103;
- certified value `X = 0.0457`, which matches **none** of the repo's LRC(14) numbers: not the Riesz
  certificate ratios (`1.0096–1.064`, THM-518 §C1), not the prime-route `2/29 = 0.069`, not the lonely
  measure `L ≈ 0.0052`.

**The ONLY thread connecting the snippet to LRC(14) is the coefficient numerator `2457 = 3·Σk²{1,…,13}`.**
Everything else (threshold, denominator primes, functional value) is disjoint from the mod-7 covering setup.
So *if* it is lonely-runner, it is **not** the repo's LRC(14) covering problem — it is a different
functional and/or a different index, sharing only the second-moment weight of the 13-speed AP.

## 3. Resolution of mac-mini's "1/25 vs 1/14" threshold question (D4)

Not a contradiction — **two different theorems on the same core**:
- `1/14 = 1/(n+1)|_{n=13}` is the **full LRC conjecture** target (repo LRC(14): 14 runners = 13 nonzero
  speeds; `inf L>0` with threshold 1/14).
- `1/25 = 1/(2n−1)|_{n=13}` is a **wider-gap-over-classical** target: the classical bound for 13 nonzero
  speeds is `1/(2n) = 1/26`, and `1/25 > 1/26`. This is the Bedert-*style* regime (beat `1/(2n)` by a
  hair), not the conjecture.

So the snippet certifies a **wider gap (1/25)** and the repo chases the **conjecture (1/14)** — same
{1,…,13}/14-runner core, non-comparable goals. THM-518's proven "Riesz stalls on AP-cores short of 1/14"
is *consistent* with a method that still clears the much lower `1/25`. The two are not in tension.

## 4. Honest adjudication: klein's reading > opus's reading, but both under-determined

- **opus's Abel–Dini/THM-2000 reading is a *general reformulation*, not a match.** `2·artanh(t)=log(S_n/S_{n−1})`
  holds for the trivial pairs `(S_{n−1},S_n)=(896,1285),(2974400,8847357)` — but *every* rational is
  `S_n/S_{n−1}` for *some* partial sums, so this identifies nothing specific. And A, B are **rational**,
  hence cannot be exact ratios of THM-2000's **transcendental** masses (`2log2`, `18−24log2`, `3log3−π/√3`,
  …). The figurate home is a flavour (819 = square-pyramidal), not a construction.
- **klein's n=13 reading rests on ONE genuinely specific fingerprint** (`2457 = 3·Σk²{1,…,13}`, unique to
  n=13; corroborated softly by the `13²` in B's denominator and the `1/(2n−1)` threshold shape). That is
  the strongest specific signal on the board — but by §2 it is *un*corroborated by the rest of the snippet's
  arithmetic, and by §1 it has no literal source. **Weight it as "one strong fingerprint," not settled.**

**Live families (unchanged from mac-mini-S169, re-ranked by this session's negatives):**
(A) extremal / log-energy beats a rational floor (lonely-runner wider-gap, merit-factor, LP-margin,
**or an AI-found "construction beats bound" certificate** — the artanh style is exactly how FunSearch/
AlphaEvolve outputs are verified; consistent with "source-generated, not indexed");
(B) irrationality/transcendence **measure** (rate comparison `C_analytic − C_arith > margin`).
Both survive. §1–§2 remove "it's literally Bedert" and "it's literally the repo's LRC(14) covering" from A.

## 5. The one decisive forward test (falsifiable, hands the anchor forward)

Reconstruct the **construction** emitting `A = 1285/896` and `B = 8847357/2974400` as two *distinct*
quantities with weight `2457/6592` and amplitudes `p_A = 1285/2181 ≈ 0.589`, `p_B ≈ 0.748 (≈¾, log_B≈log3)`.
Concrete: run a **log-energy / Riesz-entropy optimization on the 13-speed AP `{1,…,13}`** targeting a
wider gap `> 1/(2n−1)=1/25`, with the log-coefficient tied to the core's second moment `Σk²`. If the
optimum reproduces `(2457/6592, p_A, p_B)` → identification confirmed and we inherit a Lean-ready
"wider-gap" certificate template (advances OPEN-Q-104/097 either way). If the optimum's numbers differ →
the 2457 match is coincidental and family (B)/(A-AI-construction) is favored. **Diagnostic sub-signal:**
the coefficient splits as `2457/6592 = (3·Σk²{1..13}) / (2⁶·103)` — numerator pure velocity/second-moment
(primes 3,7,13), denominator foreign (103 | 5872957 = B-increment). A genuine optimizer couples a
core-intrinsic weight to a data-specific normalization; that split is what to expect if §5 succeeds.

Files: `/tmp/verify_eq27.py`, `/tmp/construct.py` (this session). Prior base decode:
`03-artifacts/drafts/artanh-1over25-snippet-decode-kind-pasteur-S129.md`.
