# LRC(14) finite-check feasibility ledger — the Rosenfeld/S–T route priced for n=14

> **SCOPED COST SNAPSHOT:** retain its estimates and architecture, but take
> mathematical status from [`CURRENT-FRONTIER.md#lrc14`](CURRENT-FRONTIER.md#lrc14)
> and literature scope from [`../05-knowledge/reference/CORE-PAPERS.md`](../05-knowledge/reference/CORE-PAPERS.md).

**kind-pasteur-2026-07-19-S128c87** (HYP-7921; executes the S128c86 backlog lead (i)).
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


## 4b. ACCEPTANCE TEST RESULT (S128c88, HYP-7930 — run same day, owner-directed)

The §4 proposal was executed at p = 43, 61, 71 (complete) and p = 191 (budgeted),
with a 9-gate harness (dual consistency; (1..13) improper mod every tested p; GW
improper; k=7/p=31 replication; A=B=C cover-count equality at completed primes).

**Verdict 1 — the pruning bet is REFUTED.** Node factors: naive-lex/MRV = 40× / 194× /
257× (growing), but MRV/(MRV + comb-capacity + slow-gap-run cuts) = **1.00× / 1.03× /
1.06×** — nowhere near the 10³ threshold. Mechanism: analytic coverage upper bounds are
subsumed by the search's exact uncovered-set bookkeeping; they can only win by seeing
future structure, which top-r capacity bounds don't. Level-1 search is NOT where repo
theorems pay.

**Verdict 2 — level-1 improperness is generic, and the ansatz describes the wrong
level.** Minimal covers: 1,280 / 25,711 / 260,568 at p = 43/61/71 (×20 per prime
step); the full known-family ansatz classifies 10 / 0 / 11 of them. The near-AP
ansatz is a property of tuples surviving the WHOLE lift tower, not of level 1.
(Caveat: the danger fraction (2·⌊p/14⌋+1)/p is inflated at small p — 0.163 at 43 vs
0.141 at 191 — so growth should soften at scale.)

**The redirect (what §4's table should have said):** level-1 at p=191 is a
compiled-code triviality (measured Python MRV rate 310k nodes/s; ~10¹¹–10¹² nodes
extrapolated ⟹ hours–days single-core in C, per prime, embarrassingly parallel —
consistent with §2's 50-core-year aggregate). The expensive object is the LIFT TOWER
fed by a huge level-1 sea. The repo's leverage is therefore **lift-level kill
theorems**: witness-quantization cages (opus-S400's HYP-7920 is exactly this at k=12
— every computational certificate strictly above 1/13, forcing AP-congruence mod ~90
primes at once), spread lemmas at the LIFTED moduli 2ᵃ7ᵇp (the LRCMod19/23Spread
template, re-aimed), and the primorial gate laws (which p×lift combinations admit
clean kills). The ×7 lift (§3b) remains the structural risk these theorems must
address. Scripts + frozen outs: lrc14_I13p1_acceptance_test_kps_S128c88,
lrc14_I13p1_minimal_covers_kps_S128c88.


## 4c. THE k=13 CONDITIONAL TIGHTNESS CAGE (S128c89, HYP-7940 — the first lift-level
kill theorem instance, built same day as the redirect)

The §4b redirect's item (a), executed. Full statement in HYP-7940; script
`lrc14_k13_cage_kps_S128c89.py` (+frozen out). Summary: CONDITIONAL on the un-run
k=13 sieve (D = {1,2,4,8} × 107 primes 191..859) terminating with terminal classes
⊆ {AP, GW}-dilates, every primitive 13-family with M < 1/14 + 1/48104 (the weakest
14-coprime rung, 491/6872) and max speed ≤ 281,577 IS {1..13} or {1..11,13,24} —
tight-locus rigidity + micro-gap emptiness at the k=12 cage's full strength
(258,276), despite n=14's two-family tight locus. The repair that makes this
possible: the c-free separator J = S₂S₆/S₄² (J(AP) ≠ J(GW), separator integer
divisible by only one caging prime, 443) converts the T=2 pigeonhole (which
collapses H0 to 489) into a two-stage elimination (stage-1 forcing height 1.03e8 ≫
stage-2 tower 281,577).

Two structural additions to §3's wall analysis: (a) **the ×7-lift quantization
degeneracy** — at 14 | lp the certificate rung is 1/14 exactly (zero margin), so the
lifts forced by composite 14 are precisely where the cage's strictness dies; the
apex-7 wall's third appearance. Cage rungs must come from the 14-coprime grid only.
(b) **the terminal table is unconditionally computed**: T = 2 over the entire known
ladder — every non-tight family is captured by the sieve's very first primes
(3/41 at q=193; 2/27, 3/40, 4/53, 5/66, 1/13-tier, deep-well all at q=191), so the
sieve's grind is concentrated entirely on the two tight families and whatever
exotic survivors exist (Wall A's shadow: each new terminal family would cost one
more J-separator stage, not a collapse).

## 5. Honest caveats

1. **The "73 primes to 733" figure from the S–T fetch is arithmetically impossible**
   (73 primes below 733 have Σ ln p ≤ 73·6.6 = 482 < 546) and is likely a misread of
   their Table 1 by the HTML extraction. My greedy reconstruction (~91 primes in
   [167,727], Σ ln p ≈ 545.9 ≈ the target after the lcm credit) is self-consistent and
   replicates the k=7 ground truth (28 vs paper's 27 primes, same range). **RESOLVED
   same day: opus-S400's independent S–T proof-mining session extracted "~90 primes in
   [167,733]" — confirming the ~91 reconstruction over the fetched 73.** The 46× ratio
   uses my greedy sets on BOTH sides, so it is internally consistent.
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

HYP-7921 (this session) · HYP-7890 §4 (the reframing that filed this lead) ·
HYP-4055 (q* ≤ 13 ln M, the repo's own 2026-07-03 finite-check reduction — compatible:
its q* ≈ 450 scale matches the prime budget's 191..859 range) · HYP-6100 (klein's
literature assessment: Tao Thm22 v_max ≤ 15; Pandey ≤ 23; MSS Thm21) · THM-574/576
(S–T deep-read, cap side) · boxeph-S130 §4 (certificate rungs = the same pruning
taxonomy) · opus-S399 §7.1 (S–T equality-case lead — the OTHER thing to extract from
the same paper) · Rosenfeld arXiv:2509.14111 · S–T arXiv:2604.23906.
