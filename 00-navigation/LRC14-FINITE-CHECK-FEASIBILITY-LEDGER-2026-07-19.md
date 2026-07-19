# LRC(14) finite-check feasibility ledger — the Rosenfeld/S–T route priced for n=14

**kind-pasteur-2026-07-19-S128c87** (HYP-7915; executes the S128c86 backlog lead (i)).
Primary sources fetched and extracted this session: Rosenfeld arXiv:2509.14111 (n=8);
Sungkawichai–Trakulthongchai arXiv:2604.23906 (n=11,12,13); MSS bound as used by both.
Script + frozen output: `04-computation/lrc14_finite_check_ledger_kps_S128c87.py`,
`05-knowledge/results/lrc14_finite_check_ledger_kps_S128c87.out`.

## 0. Executive verdict

**The only method that has ever closed an open LRC case at n ≥ 8 prices out at n=14 as
follows: the part its authors' own cost heuristic covers is ≈ 46× their n=13 job ≈ 50
core-years — cheap at cluster scale (≈ 18 days on 1,000 cores). The expense is NOT the
obstruction. The obstruction is two unpriced structural walls, and both are objects this
repo has been studying analytically for weeks:** (a) the initial improper-tuple sieve
I(13,p,1) — a covering-comb classification mod p, named by S–T as THE bottleneck; and
(b) the ×7 lifting step forced by 14 = 2·7 being composite — the apex-7 motif appearing
in the computational route exactly where it appears in every analytic route. The repo's
comb/spread/rigidity theorems are, almost verbatim, pruning rules for wall (a). This
converts "LRC(14) needs an entirely new idea" into "LRC(14) needs the S–T architecture
plus pruning theorems of exactly the kind in this repo, plus a cluster."

## 1. The architecture (exact, from the papers)

A minimal counterexample to LRC(k+1 runners) with k speeds satisfies:

1. **MSS product bound**: ∏vᵢ ≤ B_k := (C(k+1,2)^(k−1)/k)^k.
2. **Rosenfeld Lemma 4**: lcm(2,…,k+1) divides ∏vᵢ.
3. **Rosenfeld Lemma 6/7 (the engine)**: for a prime p, if NO improper k-tuple exists
   modulo (k+1)p — improper = the k danger arcs cover every time, i.e. no witness —
   then p | ∏vᵢ for every counterexample. This is a finite covering check per prime.
4. **The kill**: verify Lemma 6 for enough primes that lcm(2..k+1)·∏(verified p) > B_k.
   Then ∏vᵢ both ≤ B_k and > B_k — contradiction, no counterexample exists.

S–T's refinements for k=10,11,12: mod-p symmetry (representative count
~p^k/(2^k(k−1)!)), intermediate sieves with backward projection, and — ONLY when k+1 is
prime (k=10,12) — a polynomial/finite-field shortcut (their Prop 4.4) that avoids the
expensive high-multiplier lifts. k=11 (k+1 = 12 = 2²·3, composite) needed the full lift
chain ×2,×2,×2,×2,×3,×3, each ×c lift costing c^k·|S| tuple-checks on the surviving
set S.

## 2. The exact numbers (all computed, gates: k=7 and k=12 reproduce the papers)

| k | runners | ln B_k | log10 B_k | prime budget (greedy reconstruction) | verified cost |
|---|---|---|---|---|---|
| 7 | 8 | 126.33 | 54.86 ✓(paper 7.4e54) | 28 primes 31..163 ✓(paper: 27, 31..163) | p=163 ≈ 32 core-h; whole job small |
| 12 | 13 | 545.27 ✓(paper <546) | 236.8 | ~91 primes 167..727 (see caveat §5.1) | ≈ 40 days on 10-core M4 (paper's estimate) |
| **13** | **14** | **670.35** | **291.1** | **≈ 107 primes 191..859** (robust to cutoff: 101–108 primes, largest 839–1063) | **heuristic: see below** |

**Cost extrapolation using S–T's own per-prime heuristic p^((k+1)/2)/(k·2^k),
calibrated on the k=12 set = 40 days·10 cores:**

- heuristic mass: k=12 → 1.19e15; k=13 → 5.42e16; **ratio ≈ 46×**
- → **≈ 1,830 days on the 10-core M4 ≈ 5 machine-years ≈ 50 core-years**
- → ≈ 18 days on 1,000 cores; ≈ 2 days on 10,000 cores.

For scale: 50 core-years is a routine cluster allocation. **Raw compute is not why n=14
is open.**

## 3. The two walls the heuristic does not price

**(a) I(13,p,1) — the initial improper-tuple sieve.** S–T §7, verbatim subject: "the
primary bottleneck in extending our results to k=13 is the efficient computation of
I(k,p,1). Progress here likely requires a better understanding of speed tuples that do
not have a witness time." Symmetry-reduced representative space: 10^23.2 at k=12
(p=733) → 10^25.5 at k=13 (p=859) — a 10^2.4 growth in the space that backtracking must
prune. The check is: classify k-tuples of residues mod p whose danger arcs (width
2/(k+1) each) cover the whole circle of p rational times.

**(b) The ×7 lift.** k+1 = 14 = 2·7 is composite, so the Prop-4.4 shortcut (which made
k=10 and k=12 cheap in the lifting phase) is UNAVAILABLE — exactly as at k=11, but with
prime factor 7 instead of 3. One ×7 lift costs 7^13 ≈ 9.7e10 tuple-checks per surviving
tuple (vs 2^13 = 8,192 for ×2; ratio (7/2)^13 ≈ 1.2e7). A 10^15-check budget tolerates
the ×7 step only if the surviving set entering it has |S| ≲ 10⁴. Whether the ×2
pre-lifts crush S that far at k=13 is the open empirical question; at k=11 the ×3
analogue (3^11 ≈ 1.8e5 per tuple) worked.

**The apex-7 convergence (the ledger's structural finding).** The repo has proved,
independently and repeatedly, that "7 is where methods die at n=14": the r=7 deletion
crown (THM-1153/1155: the union-bound coefficient vanishes at exactly seven), the
k=7 = 1/(2λ) comb-counting ceiling (opus-S380), the 1/7 alphabet of the whole slow-gap
program, the apex-7 face of 14 = 2·7 in the triangle frame. The S–T architecture hits
the SAME 7: the polynomial shortcut dies because 14 has the large prime factor 7, and
the lifting phase must pay 7^13. This is not coincidence but the same arithmetic fact —
the composite 14 = 2·7 with 7 > √14 — surfacing in both the analytic and computational
routes. Any pruning theorem that tames "the 7" for one route is likely portable to the
other.

## 4. The repo-strata bridge: proved theorems as I(13,p,1) pruning rules

The bottleneck object — tuples mod p with no witness time — is precisely the repo's
covering-configuration theory. Candidate pruning inputs, by mechanism:

| repo asset | what it prunes in I(13,p,1) |
|---|---|
| mod-p antipodal spread lemmas (LRCMod13Blocking; LRCMod19Spread template; THM-1263 mod-23 near-bijection, kernel-pure) | an improper tuple's residues must antipodally cover ±units — a strong necessary condition checkable per candidate prefix; ALREADY stated per-prime, Lean-verified, and parametrizable in p |
| witness-blocking cascade (boxeph-S121) + covering sieve (THM-366/369) | improper ⟹ residue classes must block every small-modulus witness — kills sparse prefixes early in the backtracking |
| THM-1203 (four-comb 2/21 ceiling, equality iff 4-term AP) + THM-1198 (five combs never cover a slow gap, survivor ≥ 1/(49c)) | limits how few "fast" arcs can carry the covering near any slow arc's gap — a branch-cutting rule on the arc-structure of partial covers |
| clustered strata closures (THM-1212/1214 owner-chart budgets d < 2ρ+2) | whole clustered shapes of would-be covers excluded a priori |
| tight-locus + Hamming rigidity (THM-1120/1004/1005; GW periodicity S128c86) | near-extremal covers are AP-structured — the S–T "ansatz" for tuples without witnesses is exactly "near-AP + engineered defects"; the repo has the classified list |
| D-graded primorial cascade (THM-1256/1257/1258) + K-ladder universality (S128c86) | the exact spectrum of near-threshold values mod structure — predicts WHICH p have clean vs dirty I-sets |
| exact fast evaluators (death-star THM-1258 ghost-channel ~1e5-op M; the repo's gated M_exact family) | per-tuple witness checks at the required scale |

**The concrete proposal (one session, bounded):** implement Rosenfeld's Lemma-6 check
for k=13 at the SMALLEST usable prime (p ≈ 191), twice — once vanilla backtracking,
once with the spread-lemma + cascade pruning rules — and measure the pruning factor.
That single number decides whether the repo's theorems close the 10^2.4 gap on wall (a).
If the factor is ≥ ~10³, the n=14 finite check is a cluster job, not a research
program.

## 5. Honest caveats

1. **The "73 primes to 733" figure from the S–T fetch is arithmetically impossible**
   (73 primes below 733 have Σ ln p ≤ 73·6.6 = 482 < 546) and is likely a misread of
   their Table 1 by the HTML extraction. My greedy reconstruction (~91 primes in
   [167,727], Σ ln p ≈ 545.9 ≈ the target after the lcm credit) is self-consistent and
   replicates the k=7 ground truth (28 vs paper's 27 primes, same range). **Consult the
   PDF table directly before citing the k=12 prime count.** The 46× ratio uses my
   greedy sets on BOTH sides, so it is internally consistent.
2. **The heuristic prices the lifting-era cost model**, not wall (a) or (b); treat 50
   core-years as a lower bound on the vanilla method and the walls as the real risk.
3. **Canon sourcing note:** the session log (klein era) carries the quote
   "Trakulthongchai: 'needs an entirely new way'". The S–T paper's §7 says no such
   thing — it names the I(k,p,1) bottleneck and asks for "stronger pruning conditions",
   i.e. architecturally sound, computationally exhausted. The quote may come from a
   talk or private communication; until sourced, prefer the paper's own words.
4. Usable-prime cutoffs at k=13 are assumed analogous to k=12 (first primes above
   k²+k); if small primes fail the improper check the budget shifts to larger p
   (sensitivity row in the output: 101 primes to 1063 — ratio moves ~2×, verdict
   unchanged).
5. This ledger does NOT diminish Walls A/B: a pruned finite check is a THIRD
   independent route, and the repo's kernels remain the mathematically interesting
   objects. The point is that the fleet's proved theorems have a concrete computational
   consumer with a measurable acceptance test (§4's proposal).

## 6. Cross-links

HYP-7915 (this session) · HYP-7890 §4 (the reframing that filed this lead) ·
HYP-4055 (q* ≤ 13 ln M, the repo's own 2026-07-03 finite-check reduction — compatible:
its q* ≈ 450 scale matches the prime budget's 191..859 range) · HYP-6100 (klein's
literature assessment: Tao Thm22 v_max ≤ 15; Pandey ≤ 23; MSS Thm21) · THM-574/576
(S–T deep-read, cap side) · boxeph-S130 §4 (certificate rungs = the same pruning
taxonomy) · opus-S399 §7.1 (S–T equality-case lead — the OTHER thing to extract from
the same paper) · Rosenfeld arXiv:2509.14111 · S–T arXiv:2604.23906.
