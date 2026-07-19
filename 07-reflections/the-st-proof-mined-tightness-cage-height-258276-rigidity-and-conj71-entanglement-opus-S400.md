# The S-T proof, mined: the tightness cage, height-258276 rigidity, and the Conjecture-7.1 entanglement

**Instance:** opus-2026-07-19-S400 (owner: "mine the Sungkawichai-Trakulthongchai equality case,
spend a full session on it" — executing boxeph-S114's proposal, headlined in the S399 synthesis).
**HYP-7920.** Scripts: `04-computation/lrc13_st_cage_height_bound_opus_S400.py`,
`04-computation/lrc13_eleven_even_branch_and_conj71_probe_opus_S400.py` (+ frozen `.out`s).

---

## 0. The verdict, in three layers

boxeph-S114 proposed: *"if the proof pins M(C)=1/13 ⟺ AP, HYP-4382 follows and the residual
collapses — the most plausible route, living inside a proof not in this repo."*

1. **Statement level: NO** — and canon already knew. kind-pasteur-S1 (2026-07-05, HYP-4096)
   read the paper and recorded: it "does *not* characterize extremizers (it only cites
   Goddyn–Wong for the word 'tight')." The paper's own §7 confirms: the sieve identifies
   (1,…,k) as tight but "no direct analysis of tightness structure emerges from the sieving."
   My S399 synthesis headlined this lead without that grep — MISTAKE-189 logged (§7 below).
2. **Proof level: YES, substantially** — the proof's *internal quantization* pins far more
   than its statements. Mining the witness denominators yields a **tightness cage** (every
   near-tight 12-family is dilated-AP-congruent modulo every sieve prime not dividing its
   product), which converts — via a power-sum/CRT forcing argument — into
   **tight-locus rigidity and micro-gap emptiness up to height 258,276** (prior state:
   exhaustion to height ~200). §§2–5.
3. **Bridge level: two live two-way bridges** — their named k=13 bottleneck is the repo's
   core expertise (§6.1), and their concluding Conjecture 7.1 is quantitatively entangled
   with the repo's stability-gap program, with probe data produced this session (§6.2).

Nothing here closes HYP-7310 (the n=12 AP-uniqueness inverse) at unbounded height. The height
obstruction kps-S1 diagnosed remains the wall. What changed: the bounded-height zone is now
~1300× taller, proof-grade instead of exhaustion-grade, and the unbounded-height problem
acquires a new necessary-condition battery.

## 1. Their proof, decoded (extraction basis)

Paper: **"Eleven, twelve, and thirteen lonely runners"**, Sungkawichai–Trakulthongchai,
arXiv:2604.23906 — proves LRC(k) for k ∈ {10,11,12} (speeds; = 11,12,13 runners). Taken as
TRUE per the owner's LRC(≤13) policy. Read via arXiv HTML extraction in four targeted passes;
**a direct PDF verification pass is the named follow-up** (§8 lists exactly what to re-check).

Architecture (k=12, threshold 1/13):
- **Properness** (Def 2.1): a residue class v mod l·p is proper iff [∃ witness t ∈ (1/(lp))ℤ
  with ‖tvᵢ‖ ≥ 1/13 for the whole class] OR [gcd branch: gcd(l, all-but-one coordinates) > 1].
- **Sieve** (Prop 3.1, §5): per prime p, lift by c=2 through depths l ∈ {1,2,4,8}, eliminating
  proper classes; **terminal state: only equivalence-classes of (1,…,12) remain** (equivalence
  = permutation + per-coordinate sign + global unit scaling mod p). The c=13 lift is never
  run (Remark 3.2: (1,…,k) is lonely only at t=s/(k+1), so I(k,p,l)=∅ would need 13|l —
  "extremely expensive"); instead:
- **Polynomial method** (Prop 4.1→4.4): for k+1 and p > k(k+1) both odd primes, every class
  over the AP-mod-p stratum is (k,p,k+1)-proper, via witness t = s/13 + r/p whose coordinates
  land in [Cᵢ/13, (Cᵢ+1)/13), Cᵢ ∈ {1,…,11}. Three cases: all-zero mod 13 (gcd branch);
  no-zero mod 13 (pure t=s/13 witness — **non-strict**, value exactly 1/13 iff the mod-13
  residues cover all six antipodal pairs); mixed (strict, the r/p part forces interiority).
- **Assembly** (§6): J(12,p)=∅ for ~90 primes p ∈ [167,733] (Table 1, ln∏ > 547) ⟹ any
  counterexample's product divisible by ∏P₁₂ > B₁₂ (Lemma 2.6, from Malikiosis–Santos–
  Schymura), contradiction.

## 2. The one-line quantization, and the cage

**Quantization (repo-elementary):** a witness t = a/(lp) gives every ‖wᵢt‖ ∈ (1/(lp))ℤ, so a
proper-via-witness class guarantees M ≥ ⌈lp/13⌉/(lp), which is **strictly** > 1/13 whenever
13 ∤ lp. For their k=12 pipeline, 13 ∤ lp always (l ∈ {1,2,4,8}, p ∈ [167,733]); the weakest
rung over the pipeline is at lp = 2·733 = 1466 (=8·733 reduced):

> **cage threshold = ⌈1466/13⌉/1466 = 113/1466 = 1/13 + 3/19058 ≈ 1/13 + 1.574·10⁻⁴.**

**Derived cage (conditional on extraction facts F1–F3 of §8):** let W be a primitive
12-family with M(W) < 113/1466. Then for every p ∈ P₁₂ with p ∤ w₁⋯w₁₂, either
(a) {±wᵢ mod p} = c·{±1,…,±12} as multisets for some unit c (the dilated-AP class), or
(b) eleven of the twelve wᵢ are even (the gcd branch, l ∈ {2,4,8} ⟹ shared factor 2).
*Proof:* W's class mod p is not eliminated-by-witness at any pipeline depth (a witness would
force M ≥ 113/1466); so it is either gcd-branch-proper (⟹ (b), a global integer condition) or
survives to the terminal state ⟹ AP-equivalent (⟹ (a)). Witness values are preserved under
their §5.1 equivalences (unit scaling reparametrizes t; signs and permutations are free). ∎

## 3. The branch lemma (unconditional, proved this session)

**Lemma (11-even branch):** any 12-family with ≥ 11 even speeds and one odd speed has
M ≥ 1/12. *Proof:* write the evens as 2u₁,…,2u₁₁ and the odd speed b. LRC(12) (11 speeds,
settled) gives τ with min‖uᵢτ‖ ≥ 1/12. The two lifts t = τ/2 and t = τ/2 + 1/2 have identical
even-runner distances ‖2uᵢt‖ = ‖uᵢτ‖, while the odd runner's two positions differ by exactly
b/2 ≡ 1/2 (mod 1), so one lift has ‖bt‖ ≥ 1/4. Hence M ≥ min(1/12, 1/4) = 1/12. ∎
(The THM-760 sheet-dodge mechanism at c=2, one rung down; exact verification 300/300 random
families in the part-2 script. All-12-even is non-primitive.)

Since 1/12 > 113/1466, **branch (b) never occurs below the cage threshold**: the cage's
conclusion for primitive W with M(W) < 113/1466 is (a) at every non-dividing sieve prime,
unconditionally on branches.

## 4. From cage to rigidity: the power-sum forcing and H₀ = 258,276

Eliminating the unknown dilation c_p: for m = 1,…,12 define the integers
R_m(W) := P₂(W)^m·S₂ₘ − S₂^m·P₂ₘ(W), with P₂ₘ(W) = Σwᵢ^{2m}, S₂ₘ = Σ₁¹²i^{2m} (even power
sums kill the signs, symmetry kills the permutation, the homogeneous pairing kills c_p).
Cage ⟹ p | R_m(W) for every caging prime. If some R_m ≠ 0 then |R_m| ≥ ∏(caging primes)
≥ e^547/733^ℓ with ℓ ≤ 12⌊ln H/ln 167⌋ (H = max W; primes ≥ 167 dividing the product).
Exact-arithmetic binary search (script, frozen): the forcing inequality holds for all m up to

> **H₀ = 258,276.**

R_m = 0 for all m = 1..12 pins the multiset {wᵢ²} = λ·{i²} (Newton, char 0), so wᵢ = t·σ(i)
with t rational, and primitivity gives t = 1. Newton-step sanity: AP and dilates all-zero ✓;
0/2000 random non-dilates all-zero ✓.

**Derived theorem (conditional on F1–F3 only):** every primitive 12-family W with
M(W) < 113/1466 and max(W) ≤ 258,276 is {1,…,12}. Equivalently:
- **tight-locus rigidity to height 258,276**: M(W) = 1/13 ⟹ W = {1,…,12} for max ≤ H₀;
- **micro-gap emptiness**: no primitive W with max ≤ H₀ has M(W) ∈ (1/13, 113/1466).

Prior state (kps-S1/mac-mini-S47): exhaustive verification to height ~200, plus the
elementary necessary conditions (N1) w₁₂ ≥ 12w₁, (N2) w₁₂ ≤ 12w₁₁. This is a ~1300× height
extension, by argument rather than enumeration, and the micro-gap statement is the first
gap-emptiness result of any width at 12 speeds.

## 5. What it feeds (and what it does not)

- **Wall A wiring.** THM-1017's open half concerns covering 13-families V with M(V) < 1/13
  and their 12-cores W. LRC(13) gives M(W) ≥ 1/13; the inverse theorem's hard branch is
  exactly W *near-tight* (the fattening room vanishes as M(W) → 1/13). The cage says: below
  height 258,276 that branch is **only the AP itself** — i.e., the inverse theorem's
  conclusion holds on the near-tight branch at bounded height; and at unbounded height the
  near-tight branch carries the full congruence battery (a). Any future Wall-A argument can
  cite: "W near-tight (M < 113/1466) ⟹ W ≡ c_p·(±perm)AP mod every non-dividing p ∈ P₁₂."
- **The mod-13 face** is recovered *en passant* from their Prop 4.4 case analysis (all
  residues nonzero + all six antipodal pairs covered, else M ≥ 2/13) — matching the repo's
  boxeph mod-13 machinery; independent convergence, no new content.
- **NOT delivered:** HYP-7310 at unbounded height; the (1/13, 2/25) gap beyond the micro-gap
  (certificate rungs at large p enter the gap — e.g. 8/103 — so the cage cannot reach 2/25);
  anything at 13 speeds/n=14 directly (k+1 = 14 composite kills Prop 4.1 there, as
  kps-S31ag/klein-S151 established). The height obstruction kps-S1 named is intact — the
  cage RAISES the fence, it does not remove it.

## 6. The two bridges

### 6.1 Their bottleneck is the repo's expertise
§7: "The primary bottleneck in extending our results to k=13 is the efficient computation of
I(k,p,1). Progress here likely requires a better understanding of speed tuples that do not
have a witness time in an ansatz." The repo owns exactly this theory: which families escape
denominator-q witnesses is the certificate-rung ladder (boxeph-S130), the q₀-stratification
(THM-1105/1210), the near-AP Hamming banks, the (D,s)/slack frame, and the detection-floor
discipline. A concrete joint target: characterize I(13,p,1) structurally (the level-1
improper classes at 13 speeds) using the repo's rung/stratum theory — both a contribution to
their program and a test of the repo's frame on their data. Note their k=11 case (k+1 = 12
composite) sieved to EMPTY because the AP's witness denominator 12 = 2²·3 is reachable by
cheap c ∈ {2,3} lifts; at k=13 (k+1 = 14 = 2·7) the analogous lift is c = 7 — expensive but
not 13-expensive; the bottleneck is the raw size of I(13,p,1), which structural pruning
attacks directly.

### 6.2 Conjecture 7.1 is entangled with the stability gap
Their concluding Conjecture 7.1: ∃D such that every d ≥ D gives every non-tight coprime v a
witness in (1/d)ℤ. Probe (part-2 script, exact arithmetic, k=13 speeds, threshold 1/14): the
slack-1 family {1,…,11,13,36} (M = 3/41, the repo's THM-1230 gap witness) has **961 BAD
denominators in [20,1500], 235 of them ≥ 1000, persisting to the scan edge**; the ladder
families' BAD sets shrink as M rises (m=5: largest BAD 729; m=8: 478); the loose control
{1..12,15} dies at 155. Mechanism: witness-interval width scales like (M − 1/14)/w_max, so
good denominators below ~1/width are only the lucky (Ostrowski-aligned) ones. Hence:
**D(k) < ∞ requires the witness-width infimum over non-tight families to be positive — which
at k=13 is exactly the repo's "(1/14, 3/41) emptiness + 3/41 isolated" question (THM-1268),
and at k=12 the (1/13, ·) gap question.** If the gap is populated with M → floor⁺,
Conjecture 7.1 fails at that k; if the spectrum above the floor is discrete, D(k) is
computable from the second value. The repo's Wall-A program is, verbatim, the study of their
concluding conjecture. (Evidence, not a refutation; logged in HYP-7920.)

## 7. Honesty ledger

- **MISTAKE-189 (logged):** S399 headlined this lead as "unclaimed" without grepping canon
  for the statement; HYP-4096 had already answered the statement-level half on 07-05. The
  session recovered because the proof-level half was genuinely unmined — but the rule
  (MISTAKE-183: grep the STATEMENT and its constants first) fired again, on the agent who
  had just written a synthesis praising it.
- **Extraction fidelity:** all of §§2–5 is conditional on three facts read via HTML
  extraction (F1: pipeline depths l ∈ {1,2,4,8} for k=12, no c=13 lift, P₁₂ ⊂ [167,733]
  with ln∏ > 547; F2: Def 2.1's two properness branches as stated; F3: terminal S₅ ⊆
  AP-equivalents per prime, equivalence = perm/sign/unit). The quantization line, branch
  lemma, power-sum forcing, and H₀ computation are unconditional given F1–F3.
- The derived theorem is filed as **HYP-7920 (DERIVED, verification-pending)**, not a THM.
  Repo rule: never claim proved what you haven't seen proved — I have seen extractions, not
  the paper. One session upgrades it: pull the PDF, verify F1–F3 and Table 1's exact prime
  list (which also sharpens H₀ — the e^547/733^ℓ bound is conservative), then promote to a
  THM-with-citation and hand the quantization + branch lemmas to the Lean crew (both are
  kernel-ready shapes; the branch lemma needs only LRC(12)-citation + two-sheet arithmetic).

## 8. PDF verification checklist (for the follow-up session)

(i) Table 1: exact P₁₂ list and the ln∏ > 547 line. (ii) §5.2's k=12 lift diagram: confirm
only c=2 lifts (depths {1,2,4,8}) and per-prime uniformity. (iii) Def 2.1 verbatim (gcd
branch form; any third properness clause?). (iv) Prop 4.4's proof: the three-case split and
the witness interval [Cᵢ/13, (Cᵢ+1)/13) interiority for r ≢ 0. (v) §5.1 equivalence relation
verbatim. (vi) Whether elimination ever uses any mechanism other than Def-2.1 properness.
Each of F1–F3 maps onto these; any deviation changes constants (threshold, H₀) but the
derivation template survives unless a 13 | lp witness stage exists for k=12 (it must not, by
Remark 3.2's own economics).

## 9. Cross-links

boxeph-S114 (the proposal, MSG-1823) · HYP-4096 / kps-S1 (statement-level verdict + N1/N2 +
the height diagnosis) · HYP-7920 (this) · MISTAKE-189 · THM-760 (sheet dodge) · THM-1017 /
HYP-7310 / HYP-4382 (Wall A) · THM-1105/1210/1215-corrected, boxeph-S130 rung ladder (the
I(13,p,1) bridge) · THM-1230/1268 (the 3/41 witness and gap forcing; Conj-7.1 entanglement) ·
kps-S31ag, klein-S151 (composite-14 wall) · S399 synthesis §7 lead 1 (executed here) ·
scripts + frozen outs (S400).
