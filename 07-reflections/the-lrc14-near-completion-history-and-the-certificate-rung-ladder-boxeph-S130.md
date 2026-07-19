# The LRC(14) near-completion history, and the certificate-rung ladder (think 2/19)

**boxeph-2026-07-19-S130** · owner-directed review: *"see the history of the work on the
14-runner Lonely Runner Conjecture, see how we have gotten close — where agents thought the
proof was complete — and the similar mistakes and repeated work; look for the insights and
clarification it provides; think 2/19."*

Method: three very-thorough Explore sweeps (session log 72,690 lines; MISTAKES.md ~100
LRC-era entries; the finish-map/frontier document sequence + current PROOF-MAP), plus one
exact computation (`04-computation/lrc14_certificate_rung_profile_boxeph_S130.py`,
HYP-7880). This builds on the only prior review of this kind — the owner-commissioned
**death-star-S58 synthesis** (2026-07-18, SESSION-LOG:1241), which named five failure
genera; the new content here is the *arithmetic mechanism underneath the genera* (§4), the
full episode chronicle (§2), and the forward levers it yields (§7).

---

## 1. The one-sentence verdict

**Nobody ever pushed a false QED to canon; what recurred — at least fourteen times — is a
true statement of the form "LRC(14) is closed modulo X" in which X, examined later, turned
out to contain the whole conjecture.** The residual was never eliminated; it was renamed.
And the renaming has an arithmetic explanation (§4): every tool that "closes 95%" is a
bounded-denominator certificate; the surviving 5% is the class of families that
block-or-saturate *every* such certificate, which is exactly the additively-structured
(Freiman/near-AP) class where the extremals live — so each new tool rediscovers the same
core under a new name.

## 2. The chronicle: fourteen near-completion episodes

(Full sourcing: SESSION-LOG.md line numbers and commits as cited; verified spot-checks on
lines 1241, 9762, 2255.)

| # | date | claim | undone by | mechanism |
|---|------|-------|-----------|-----------|
| 1 | 07-11 | FINISH-MAP: "proving EITHER route finishes LRC(14)" (klein-S258) | superseded 07-13; 3 internal walk-backs (klein-S263/S265/S266) | Route B window unbounded, not finite; HYP-6150 false vs deep well |
| 2 | 07-13 | "any power-saving Q_s=O(r^{2-ε}) suffices" (klein-S281; FINISH-MAP-07-13 §3) | klein-S314 (THM-883) + boxeph-S96/97 (THM-886(V)) | at the peel w=d the family is maximally resonant: Q_s=Θ(r²) SHARP; Cauchy–Schwarz infinitely lossy |
| 3 | 07-13 | "LRC(14) = a single equidistribution inequality; prove EITHER" (klein-S284) | overtaken by #2 + the covering crux | the [A] face rested on #2's premise; the [B] face became the open problem |
| 4 | 07-14 | "covering case closed MODULO one bounded finite-check band" (mac-mini-S104; echoes mac-mini-S58 "closed modulo THM-527-A", mid-June) | July 18–19: the covering case is the entire open problem | the "band" was not bounded; samplers couldn't see the dilated-AP/lcm stratum (death-star-S58: "demoted to sample statistic") |
| 5 | 07-16 | "the LRC(14) mathematics in this decomposition COMPLETE; remaining = bookkeeping + Lean" (mac-mini-S113, THM-878/879) | same entry (k-uniform O(1) REFUTED, ~log growth) + #2's refutation chain | O(r) held only at the k=13 interval core; the uniform premise was false |
| 6 | 07-16 | "LRC(14) endgame fully enumerated / analytic ledger empty" (boxeph-S47–S49 certificate v2 + manifest) | July 18–19 reopening | the "ten consumable items" grew back into an open inverse theorem |
| 7 | 07-18 | THM-995(X) empirical covering floor M≥1/9 | boxeph-S85 (MISTAKE-160) | sampled floor sat ABOVE the proved covering minimum 14/183 — the sampler missed the measure-zero dilated-AP stratum |
| 8 | 07-18 | "compact covering ⟺ LRC(14)" equivalence chain (boxeph S86/S88/S94/S113/S114) | codex-S74 (MISTAKE-170, THM-1099/1149) + codex-S67 (MISTAKE-166) | a SUFFICIENT residual (1/13 floor, stronger than the 1/14 target) promoted to an equivalence; INV omitted Covering(2..14); literal INVcov refuted by the dilation 2·{1..13} |
| 9 | 07-18 | "density row closes iff A=0" (boxeph-S98) | boxeph-S99 ("OVERCLAIM — A=0 is the O(R) tail") → S100 ("|S|=o(R) is FALSE") | the sufficient target was false; only the explicit O(R) bound closes SEPARATED far elements |
| 10 | 07-18 | THM-1029 "finitely many candidate cores" (death-star-S56 cont15) | same agent, next continuation (cont16) | "the clean lemma is FALSE… reduces to the alignment wall = LRC(14)" |
| 11 | 07-18 | "sharp horn certifies the r=6 covering stratum" (death-star-S58, THM-1123/1132) | self-qualified in-entry + codex tail audit (MISTAKE-168) + MISTAKE-175/182 | midpoint cells / a min branch / an AP slice promoted to uniform horn facts; centre-hitting refuted (117 non-proportional hitters) |
| 12 | 07-19 | death-star S58f→g→h→i: "closes [missed-modulus / binding / fully-clustered / 95% of spread cores]" | each entry's own tail | each closure re-spawns "sole residual = [renamed Freiman core]" — HYP climbs 7744→7746→7748→7750, the object never closes |
| 13 | 07-19 | THM-1230/1235/1240: the (1/14, 3/41) stability interval | opus-S397/S398 amendments: "NOT settled… I did not settle it, and I want to lead with that" | the (D,s) candidate list is infinite without a bound on D; 8.5M random families are "scale, not coverage" (missed both known near-floor families) |
| 14 | 07-19 | THM-1203 §12 residue list carried as "the actual remaining passage" (codex-S77 → THM-1245 §0) | kind-pasteur-S128c85 (THM-1251) | the list was STALE one session after it was written — its own author's THM-1214 closed the stratum; a residual list is a claim about canon and must be re-checked against canon |

Older instances of the same shape: mac-mini-S13 "L7 CLOSED / ALL PASS" (reduced, not
closed); the mid-June "criterion C holds universally" false target (HYP-2580a, refuted by
S\* = {1,2,3,5,7,…,13,38,42}, M=2/23 — see §4 for why S\* is the paradigm case);
MISTAKE-116: the finite covering system {2..Q0} is complete for NO finite Q0.

**Repeated work** (distinct from repeated *mistakes*): kind-pasteur + death-star both
re-proved codex's THM-1203 (MISTAKE-183); opus re-derived the Pinch Lemma "800 sessions
later" (THM-1170 vs HYP-2059/THM-401, and nearly filed a court case against the proved
THM-401 over a denominator-reduction artifact — MISTAKE-173); boxeph-S120 re-derived the
same Pinch Lemma independently two days later; boxeph-S123 re-derived THM-633's ladder;
kind-pasteur's THM-1153 = THM-633 and its |B|=1/343 = codex's THM-1183; kind-pasteur re-derived
THM-523's q-witness (MISTAKE-158); klein-S228 duplicated kps's THM-712 formalization;
mac-mini re-derived a construction *the same agent had refuted three sessions earlier*
(MISTAKE-099). Root cause, stated twice in canon and confirmed by its author working one
session apart (MISTAKE-183/158): **searches ran on the METHOD, never on the STATEMENT —
"2/21 alone would have found it in one grep."** Plus concurrency: same-owner-prompt
collisions produced THM-882/883/884, THM-735 (triple), THM-1238→1245, HYP-7305/7315, and
a still-live HYP-7750 double-assignment (opus-S389's archaeology vs death-star-S58i's
residual) — the reserve-the-id-first rule was logged FOUR times (MISTAKE-052→053→057→058)
and never fleet-adopted (~18 duplicate IDs live in MISTAKES.md today).

## 3. The taxonomy (what actually goes wrong)

From the full MISTAKES.md sweep (~100 LRC-era entries, 2026-06-15 → 07-19):

- **~45%: instrument-blind extremizers** — a sampled/bounded/grid/limit computation
  promoted to a universal law, false because the true extremizer is an *arithmetic* family
  (dilated AP, lcm speed, deep well, mod-p pole, sweep-boundary corner) that the instrument
  structurally cannot generate. Canon tracks this as the 073→100→119, 095/096/098→101/102
  →104→110→112/113/125, and 137→139→140→141→145 chains — it recurred a dozen+ times AFTER
  being first logged. The three sharpest rules it produced: *"a sampled floor above a
  proved minimum is a red flag"* (160); *"a null result from a search is worthless without
  a POSITIVE CONTROL"* (162); *"always include near-dilate/coherent seeds — a standing
  seed-battery requirement"* (140); *"scan orbit representatives, never a coordinate slice
  — third time dilation has bitten"* (156b).
- **~14 closure overclaims** — the "modulo X" failures of §2, plus route proposals never
  instantiated (105, 157b: *"when a search keeps finding worse, try to prove it can always
  find worse before proposing the bound"*).
- **~11 incomplete exact computations** — missing breakpoints (pair sums s_i+s_j, half-turn
  cusps q=2v), coarse grids masking sub-cell sets, Fourier truncation: *"incomplete
  breakpoints… fabricate tights — the most dangerous direction for a census"* (086).
- **~10 wrong-population tests** — non-covering families contaminating covering claims
  (089, 167), reduced vs represented denominators (173, oracle pinch-M, 114b), q≤14
  excluded as "classical" exactly where the witness lives (172).
- **~8 vacuous-but-kernel-pure Lean** — quantifier slips (∀m for ∃m: MISTAKE-186; ∀ over a
  TYPE instead of reachable values: 136a "TWICE now"; unsatisfiable KCL hypothesis: 146)
  and certificates placed where their cap is impossible (154, 155a). **This is the live
  July failure surface** — the three newest mistakes (184, 185, 186) are all in this
  cluster or its analytic sibling.
- **~7 stale-canon/repeat-work**, **~9 ID/process**, **~7 wrong-tool persistence** (soft
  Weyl/large-sieve for a lower bound — every such tool caps at O(r²) or is refuted, HYP-6445;
  Koksma–Hlawka vacuous on near-APs, 127b/128; descent/sieve/measure provably too weak on
  compact cores, 160), **~6 limit-as-finite-bound** (decorrelated M→∞ facts used as finite
  majorants: 080/082/137a/145/184/185).

Hypotheses that were **relied on before refutation**: HYP-6994 & HYP-7017 (uniform Q_s —
the 07-13 finish-map premise), HYP-6445 (offdiag≤0), the finite-covering-completeness
belief (MISTAKE-116), HYP-2580a (criterion C universal), HYP-7355's evidence base (refuted
twice; the conjecture itself now stands UNTESTED at n=14).

## 4. The mechanism: the certificate-rung ladder (this is where 2/19 comes in)

The five death-star-S58 genera and the §3 taxonomy are shadows of ONE arithmetic
mechanism. It is easiest to state through the general-p spread lemma (the S115 mod-13 /
S126 mod-19 template, one line at every prime):

> **M(V) < 2/p and p divides no speed ⟹ the residues {v mod p} hit ALL (p−1)/2 antipodal
> unit-pairs of Z/p.** (At t=b/p all distances are multiples of 1/p; M<2/p forces some
> v·b ≡ 0,±1; 0 is excluded, and b ↦ ±b^{-1} sweeps the pairs.)

**RUNGS.** Every elementary tool that ever "closed 95%" of anything here is a
*bounded-denominator certificate*: one rational time t=a/q whose success is a residue
condition mod q — the sieve (1/q, q≤13), pair-blocking (2/13), the spread lemmas (2/17,
2/19, 2/23), the pigeonhole floor (2/29), the dilated sieve (1/13 at nd), rational-time
floors (d_k/k), the mod-25 floor. By the Pinch Lemma (M = D/s at a pair-sum time),
**attainable loneliness values and certificate values are the same discrete grid** — the
spectrum ladder 1/13, 1/12, 2/23, 2/21, 2/19, … IS the certificate ladder. Progress
happens in rung-jumps; every tool strands at a rung boundary; coverage percentages are
statements about the *generic* stratum between rungs.

**BLOCKING.** A family defeats a mod-p certificate only by divisibility engineering:
one multiple of p blocks the whole prime (AP13 blocks 13 with 13; the deep well with
182=14·13; **S\* blocks 19 with 38=2·19**); a *dilation factor* blocks a prime globally
for all speeds at once — the cheapest block, which is exactly why "dilation bites" keeps
refuting scans (2·{1..13} vs INVcov; MISTAKE-156b "third time"); an *lcm element* blocks
or covers many moduli simultaneously — why every bounded-denominator raw-count conjecture
died against one lcm speed (MISTAKE-157b, 096, 116). Otherwise the family must SATURATE
the certificate: full antipodal spread. The HYP-7880 computation makes this exact for the
live targets:

- every 13-set in the unsettled window (1/14, 3/41) is covering-{2..13} AND
  block-or-full-spread at p=17, 19, 23 **simultaneously** (slack at 23: TWO);
- every 12-set in the uniqueness gap (1/13, 2/25) is a **mod-23 NEAR-BIJECTION** — 12
  residues onto all 11 antipodal pairs, exactly one collision ({1..12} doubles {±11};
  {1..11,24} doubles {±1}). Nobody has used p=23; canon contains 2/23 only as an attained
  M value;
- **S\*** — the family that killed "criterion C is universal" and taught MISTAKE-076's
  lesson — is fully explained: it defeats every single-removal arc criterion, and its
  rescue ("a global witness, not any single-removal arc") is *literally the mod-23
  certificate*: unhit pair {±6}, maximizer numerator a = 4 = 6^{-1} mod 23, M = 2/23
  exactly. The criterion-universality overclaim failed because the certificate lattice is
  plural, and the sampler saw only one sheet of it.

**WHY THE RESIDUAL IS INVARIANT.** Being extremal (M near the floor) MEANS blocking or
saturating every certificate that would otherwise lift M above the target. So the
residual class of any certificate-family is always the same object: CRT-thin (measure
zero to every sampler — the ~15× "sampling artifact" genus), divisibility-rich and
additively rigid (13 speeds must service ~15 constraints — 12 covering conditions + 3
spread conditions — so only highly-composite / near-AP configurations fit: the Freiman
core), and *renamed by each tool that newly fails on it*: compact core (S86), non-dilated
core (S113), alignment wall (death-star cont16), covering-core gap (S58h),
rational-time-evasive core (S58i), escape tail (HYP-4667), bound-D obstruction (THM-1240).
Seven names, one object. The MAX-not-MEAN genus is the same mechanism one level up:
averaging (Q_s, second moments) integrates across rungs and cannot see the single
saturated rung that carries the extremal family — which is exactly how the 07-13
finish-map premise died (Q_s = Θ(r²) at the resonant peel, CS infinitely lossy).

**THE NEAR-MISS ARITHMETIC.** The pigeonhole prime for 13 moving runners is 29 (14
antipodal pairs > 13 speeds): every family without a multiple of 29 has M ≥ 2/29 — 96.6%
of 1/14. This has been canon since THM-518 ("the prime route lands at 2/29, short by
3.4%"), and MISTAKE-093 fixed the last-mile window as [2/29, 1/14), width 1/406. **The
entire conjecture lives in a sliver of width 0.0025 above the free pigeonhole floor** —
which is why "one lemma away" is a recurring *sensation*: every certificate family
genuinely does get within a few percent, and the remaining sliver genuinely does contain
the whole problem, every time. 2/19 is the mirror rung on the other side: the largest
prime certificate that clears the ENTIRE gap regime with room to spare (2/19 > 2/25 >
3/38 > 1/13), which is why the mod-19 spread condition (S126–S129, now the ledger rung
`gap_regime_mod19_spread`) constrains every family below it — and why 23 (slack 1–2) is
the sharpest unused pin between them.

## 5. The residual's trajectory (07-11 → today): renamed and sharpened, never replaced

How the "remaining" item read in each planning document:

- **07-11** (FINISH-MAP, klein-S258): "the ENTIRE remaining content is one statement about
  the covering class… two independent routes — proving EITHER finishes LRC(14)"; "both
  FINITE/BOUNDED in principle." Self-corrected mid-document: Route B's window is
  UNBOUNDED (klein-S263).
- **07-13** (FINISH-MAP, klein-S284): "one equidistribution/cancellation inequality…
  Route A is the softer face; ANY power-saving suffices." Premise refuted 07-16/18.
- **07-14** (opus-S285 frontier): "the binding low-M residue patterns — THE ONE GENUINELY
  ANALYTIC PIECE LEFT"; same day the guardrails note "THM-724 and THM-726 titles are
  stronger than the honest scope in their bodies."
- **07-16** (frontier assessment, boxeph-S23): "what remains analytic is exactly two named
  uniformities: HYP-6994 and the compact-core exposure bound"; top recommended move: "the
  decidable compact-core sweep… could CLOSE [B] without any new analysis." **Both halves
  of my own assessment were mis-aimed**: HYP-6994 was refuted within days (klein-S314 +
  boxeph-S25), and the compact core later closed by one elementary line (M ≥ 1/(2ρ),
  death-star-S58h) — because the crux was never there: "the wall is SPREAD, not
  comparable" (S56 synthesis). A §2-style episode in miniature, mine.
- **07-18** (death-star-S56 unified synthesis): crux = THM-1017 inverse theorem ("M<1/13
  covering ⟹ 12 non-max speeds form a dilated AP" = HYP-7310/HYP-4382 = Tao n=12);
  soft-route refutation recorded in place.
- **07-19** (today, operative decomposition across S58h/S58i + ledger): kernel = covering
  M<1/13 ⟹ AP core. Partition: near-AP (Hamming ≤ 6) DONE via 2/25 rigidity;
  fully-clustered (ρ < 6.5) DONE via M ≥ 1/(2ρ); spread far cores 94.9% DONE via the
  rational-time floor M ≥ max_k d_k/k — **whose class-infimum lands on 3/29, the p=29
  rung again** (D=3, s=29: even death-star's newest certificate bottoms on the ladder of
  §4); the 5% "rational-time-evasive" residual is *named to be* HYP-4382, the
  Freiman/12-set-uniqueness core. The mod-19 spread lemma (S126–S129) is a proved
  necessary condition ON that residual — it sieves but does not remove it (3/38 is
  defeated only cross-modulus; S124/S126).

**DONE items that later reopened**: the covering-min 14/183 "PROVED" (07-13) was
qualified within 24h and its global form is STILL the subject of the one active
substantive court case (CASE-convergent-not-covering-min, GRANT-recommended, opus
conceding the earlier "14/183 CONFIRMED" was a restricted-family scan); the 07-11
"density floor CLOSED" is still OPEN in the PROOF-MAP as the (A′) per-k tail-minimality
lemma; the tight locus {AP} was corrected to {AP, GW} by court case (MISTAKE-100); the
formal premise INVcov was refuted by the dilation 2·{1..13} (THM-1158).

**Threshold migration**: 14/183 (07-11/13) → 1/13 compact floor (07-18) → the direct
(1/14, 3/41) interval (07-19, opus). Three constants, one object: the n=12/near-floor
AP-uniqueness core identified on 07-11 is exactly today's "5%". Every intermediate
"remaining" was this object wearing the current tool's name.

**A live staleness warning**: `LRC14-PROOF-MAP.md` has NOT yet absorbed the four newest
frontier items (the mod-19 rung, the clustered floor, the rational-time floor, the
D-bound reduction) — they live only in Lean/theorem files/INDEX. That is precisely the
channel through which MISTAKE-183-genus staleness propagates; the next PROOF-MAP editor
should fold them in before anyone works "the remaining passage" from the map.

**The one monotone meter.** Against fourteen prose walk-backs, the Lean corpus never
retreated: 569 LRC\*.lean modules today (LRCLadderD1 2026-07-06 → LRCDilatedSieve /
LRCAPCoreBridge 07-18 → LRCMod19Spread 07-19), including the kernel-checked reduction
LRC14 ≤ LRC(≤13) + INV (S105–S109, hypothesis-repaired by codex MISTAKE-166), THM-1214's
clustered r=5 closure, the c=7/c=8 wall lines, and the mod-13/19 rungs. Kernel-purity is
not sufficient (MISTAKE-186: vacuity passes the kernel silently) — but it is the only
asset class in this project that has never been walked back. Prose "95%" claims reset at
least fourteen times; the theorem-object ledger only grew.

## 6. The rules the history writes (consolidated, each traceable to ≥2 incidents)

1. **Grep canon for the STATEMENT and its CONSTANTS before working any inherited
   residue.** A residue in a letter is a claim about canon; check it against canon
   (MISTAKE-183, applied successfully one session later by its own author; 158, 134, 099).
2. **A sampler can never validate a closure claim** — only a blocked-stratum enumeration
   can (list who blocks each certificate, then classify the blockers). Mandatory: positive
   controls (162), near-dilate seed batteries (140), orbit representatives not slices
   (156b), covering predicate imposed up front (089). A sampled floor above a proved
   minimum is a contradiction, not evidence (160).
3. **Type the residual in every "closed modulo X" claim**: X ∈ {bounded finite check (state
   the bound and who verified it), named citation, OPEN}. The fourteen episodes of §2 all
   had X mis-typed: "finite" checks that were unbounded (116, klein-S263), "bookkeeping"
   that was the inverse theorem, "equivalent" residuals that were strictly stronger (170:
   *"write every arrow with its threshold, extraction hypothesis, and reverse
   realization"*).
4. **Kernel-pure ≠ non-vacuous**: exhibit a witnessing family for every hypothesis before
   wiring (186, 136a, 146, 154, 155a). Prefer dist(·,ℤ) < c to hand-rolled quantifiers.
5. **Never use a limit/decorrelated bound as a finite majorant** (080, 082, 137a, 184,
   185) and never aim an equidistribution UPPER bound at a concentration LOWER-bound
   target (S95's Weyl lesson; 127b/128's Erdős–Turán lesson).
6. **Reserve IDs before writing** — logged four times (052→058), still unadopted, ~18
   duplicate IDs live including HYP-7750 twice this week. Grep theorems dir + INDEX tail +
   origin log; renumber losers by first-push.
7. **Measure progress by the Lean ledger and by the residual OBJECT'S definition getting
   sharper** — never by coverage percentage. The correct current statement of the residual
   is not "5%": it is a *definition* (the class blocking/saturating all bounded
   certificates) plus a conjecture (that class = dilated APs + the known killers).

## 7. Where the 2/19 lens points next (logged as HYP-7880 + backlog lead)

- **Use p=23.** The near-bijection pin (slack 1 on 12-sets, 2 on 13-sets) is strictly
  sharper than the wired mod-19 rung and costs one Lean file by the LRCMod19Spread
  template. Impose it (with covering{2..13} + mod-17/19 spread) as a FINITE CRT
  feasibility system mod lcm(8,9,7,11,13,17,19,23) on: (i) the 3/38 depth-minimal target
  (HYP-7782), (ii) opus's D-rung candidates — inside (1/14, 3/41) the rung equation
  forces s uniquely per D: **D=4 ⟹ s=55 (the 4/55 mediant), D=5 ⟹ s=69, D=6 ⟹ s=83** —
  so "bound D" = "bound the rung", and each rung is a finite residue problem; (iii)
  death-star's spread-far-core killer tuples (S58i): the killers must complete the
  antipodal covers the 8-speed core leaves open — a per-prime codimension count that
  could carve the 5% directly.
- **The blocking-budget heuristic** (conjecture, honest status OPEN): 13 speeds must
  service ~15 simultaneous certificate constraints; quantify the budget (how many
  constraints one speed can service as a function of its divisor structure) and prove a
  counting lower bound forcing near-AP structure. This would be the Diophantine→additive
  bridge (S104's missing Half 1) in finite, checkable form.
- **Do not re-aim**: order-2/scalar invariants, soft-Weyl at lower bounds, descent on
  compact cores, single-modulus kernels at 3/38 (S126: the obstruction is irreducibly
  cross-modulus), finite Q0 covering systems (116). Each is provably dead in canon.

## 8. Cross-links

Builds on: death-star-S58 synthesis (SESSION-LOG:1241) · THM-518 (2/29 diagnosis) ·
MISTAKE-093 (the 1/406 window) · the-resonance-fill-profile (S81) · S121 spectrum ladder ·
S123 determinant strata · S126 mod-19 spread · MISTAKES 073–186 (esp. 160, 162, 170, 183,
186) · HYP-7880 (this session's data) · `LRC14Ledger.lean` rung `gap_regime_mod19_spread`.
Files: `04-computation/lrc14_certificate_rung_profile_boxeph_S130.py`,
`05-knowledge/results/lrc14_certificate_rung_profile_boxeph_S130.out`.
