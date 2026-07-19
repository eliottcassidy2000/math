# The strategy cross-walk, the two-vocabulary problem, and the Goddyn–Wong periodicity

**kind-pasteur-2026-07-19-S128c86** · owner prompt: *"synthesize the history of the repo's
work on the 14-runner lonely runner conjecture extensively and look for common patterns or
mistakes and slight similarities or differences between strategies we have pursued, and
opportunities to leverage big picture understanding for critical math insights and
reframing of assumptions."*

**Concurrency note (lane split).** This prompt was fanned out. boxeph-S130 pushed first
(episode chronicle: fourteen "closed modulo X" events + the certificate-rung-ladder
mechanism + the p=23 lever — read it first, it is excellent); opus-S399 reserved a
history-synthesis reflection (six parallel readers, pending at this writing); death-star-S59
pushed the single-far absorption atlas (new math: the single-far first-gap stratum is EMPTY
at every N in 6..13 except N=7 and N=13). This session's lane, chosen to complement rather
than duplicate: (i) the **strategy cross-walk** — how the decompositions relate, where they
redundantly cover, and what each one's kernel is; (ii) the **two-vocabulary problem** — the
same near-floor spectrum theory built twice, 13 days apart, with a measurable cost caught
and corrected here; (iii) the **Goddyn–Wong periodicity** — an exact computation showing
n=14's "unique non-rigidity" is the second member of an arithmetic progression whose first
member is already proved, with the reframing that follows. Method: four parallel readers
(session log 72,690 lines; MISTAKES.md ~190 entries; the 16 frontier docs; the THM-1000+
atlas), the deep-read of all 19 recent inbox messages, and one gated exact computation
(`04-computation/lrc14_ladder_realization_crossN_kps_S128c86.py`, HYP-7890).

---

## 1. The shape of the history (compressed; sources for the long form)

Fifty days, ten phases, ~19 distinguishable strategy families (the full timeline with
session cites is in this session's reader output, folded below into what it implies; the
episode-level chronicle is boxeph-S130 §2; the frontier-document evolution — eight renamings
of "the one remaining item" — is tracked doc-by-doc in the 07-08→07-18 sequence).

The arc has **three regimes**, and the transitions are the story:

1. **Existential/witness regime** (06-16 → 07-09): find a good period / density floor.
   Output that held: the denominator sieve, THM-527's skeleton, the μ_{1/7} floor with all
   six legs k=8..13 closed (THM-651 shifted-tent at k=8; THM-653/656/660 window–energy–PZ
   at 9..13), the good-period dichotomy (LEM-012/013), the coarse reduction to LRC(≤13)
   (kps-S52/S53). Its residual — "realization/Part A" — was never discharged; the effort
   pivoted instead.
2. **Analytic-inequality regime** (07-09 → 07-16): one signed cancellation
   (hB5 / OffLine ≤ f(E3) / Q_s power-saving). Killed by its own sharpness results:
   Q_s = Θ(r²) at the resonant peel (THM-883/886), the Delsarte-LP impossibility
   (THM-1185), the "irreducibly signed" verdicts (klein-S223/S279). The correct lesson,
   proved rather than felt: **soft/averaged instruments cannot see the extremal stratum.**
3. **Inverse/rigidity regime** (07-16 → now): classify the near-extremal families
   themselves. Everything currently live is a face of this: the n=12 AP-uniqueness complex
   (HYP-7310 + crown collapse + AP-extraction, with 3/38 as boxeph's depth-minimal
   target), the comb/killer geometry (clustered r=5/r=6 CLOSED by THM-1212/1214; the
   six-comb phase stalk and r=7 crown live with codex), and the determinant/interval
   thread ((1/14, 3/41) emptiness = bound D = bounded primitive near-floor speeds,
   THM-1268/1269, ex-opus-1240/1245 after the renumber).

**The one-sentence history:** every strategy closed its generic stratum quickly, then hit a
rigidity kernel; the kernels were renamed at least eight times (boxeph-S130 counts the
"closed modulo X" episodes at fourteen); and the fleet only recently began proving
*blindness theorems* — measure-blindness (THM-1185), translation-invariance blindness
(THM-1225), detection floors (HYP-7870 IV) — that explain in advance which instruments
cannot work.

## 2. The cross-walk: eight stratifications, three kernels, and a redundancy worth using

The fleet has cut the same space of 13-speed families along at least eight axes:

| axis | split | closed side(s) | kernel left |
|---|---|---|---|
| covering | misses a q ≤ 14 / covers all | non-covering: sieve t=1/q (THM-366/369) | covering families |
| scale | multi-scale (clusterable) / single-scale | coarse reduction → LRC(≤13) (kps-S52/S53, GREEN) | single-scale compressed |
| diameter | prim-diam ≤ 75 (k=13; ledger per k) / spread | intersection ledger (kps-S59/S60) | spread families |
| additive energy | high (near-AP) / low (dissociated) | both ends: tent variance (THM-656) + LEM-013 | the middle |
| Hamming | radius ≤ 2 (n=13: THM-1004/5), ≤ 6 faces (codex) / far | near-AP rigidity | far-from-AP with structure |
| determinant | D=1 (classical sieve ≡ q₀ ≤ 14) / D ≥ 2 | branch 1 incl. every extremal (THM-1210/1215) | branch 2: q₀ > 14 needs pair D/s ≥ 1/14 |
| clustering (killers) | clustered r ≤ 6 / spread killers | THM-1212/1214 (owner charts) | spread far cores |
| favorable shape | near-AP ∪ fully-clustered (ρ<6.5) ∪ spread far core | first two (THM-1004-6; M ≥ 1/(2ρ), S58h) | spread far covering cores (5% residual, S58i) |

Three observations the table forces:

**(a) The kernels coincide in pairs but are NOT all one object.** The covering, Hamming,
clustering, and favorable-shape kernels all name the same thing — a spread, covering,
non-AP-core family near the floor — i.e. the **Freiman/E3 inverse** (HYP-4382/7310 family).
boxeph-S130 §4 is right that every bounded-denominator certificate strands on it. But the
**six-comb phase stalk** (codex: oriented germ lift on THM-1248/1250's located tree) and the
**branch-2/bound-D** statement (opus) are *different* obligations: the first lives inside a
hypothetical cover's geometry with no 12-core in sight; the second is a compactness claim
about primitive near-floor speeds. Calling everything "the same wall" is comfortable and
wrong; the honest count today is **three kernels** (inverse/rigidity, six-comb phase,
bound-D), of which the first is the deepest and the other two might be independently
closable.

> **Post-pull adjudication note (integrating opus-S399, pushed mid-session).** opus's
> concurrent synthesis counts **2.5 walls** — A = the inverse complex, B = the phase
> transport, half = n=12 sporadic — and files bound-D INSIDE Wall A via D = M·s
> ("bounding D ⟺ bounded primitive near-floor speeds ⟺ the inverse statement in (D,s)
> coordinates"). The delta with this session's count of three is exactly one judgment
> call, worth recording precisely: bound-D FOLLOWS from the inverse theorem, but it is a
> strictly weaker COMPACTNESS statement (bounded speeds, not AP structure), and weaker
> statements sometimes admit proofs their parents don't — e.g. death-star's clustered
> floor M ≥ 1/(2ρ) closed a stratum the full inverse route had circled for days. So:
> treat bound-D as **Wall A's separable outer shell** — count it inside A for
> residual-typing (opus's rule), but permit dedicated attacks on the shell that do not
> route through AP-extraction. Both syntheses agree Wall B stays distinct until a
> reduction is written, and agree on protecting it from vocabulary re-aiming.

**(b) At least three decompositions are COMPLETE — so LRC(14) has at least three
independent sufficient kernels.** Branch-1 + branch-2 = everything (opus THM-1210); the
k-leg witness architecture was complete modulo realization (boxeph-S1 banner, 07-09); the
comb program (all carrier cardinalities + clustered strata + slow-gap covers) is complete
modulo the six-comb stalk and spread r=5. Nobody has compared the three kernels' *relative
hardness* in one place, and sessions keep choosing a kernel by momentum rather than by
comparison. The bound-D kernel is the newest and least explored; the six-comb stalk has the
most machinery attached; the inverse theorem has the most failed attacks. **Rule proposal:
when a session picks a wall, it should say which of the three kernels it serves and why
that one.**

**(c) Redundant closures are an unexploited asset.** Several strata are closed twice by
different proofs (e.g. near-AP by Hamming rigidity AND by the pair-sum competitor's Schur
deficit; clustered by owner charts AND by the 1/(2ρ) floor; k=9..13 legs each by 2–3
instruments). Where two proofs of the same stratum exist, their *difference* is information
about what the open kernel needs — the instrument features that appear in BOTH working
proofs (exact residues, located pairs, budgets) vs those in neither (means, doubling,
unlocated measure). This is how THM-651 was found (the tent kink sits BELOW threshold —
a located, signed feature) and it is the cheapest heuristic we have for kernel work.

## 3. The two-vocabulary problem, with a measured cost (and a correction landed)

The sharpest process finding of this synthesis: **the repo built the same near-floor
spectrum theory twice in 13 days, in disjoint vocabularies, and the two copies never
cross-cited until death-star-S59's verify-first check this morning.**

- **2026-07-06 vocabulary (binder congruences):** opus-S117/S118/S119 + mac-mini-S26..S30 —
  HYP-4506 (first-gap emptiness is non-monotonic in N; members at N=6,7,13), HYP-4572
  (the trichotomy for F(N) = {1..N-2, N, 3(N-1)}), **HYP-4516 (the definitive mod-30
  binder gate: F(N) attains the mediant 3/(3N+2) ⟺ N ≡ 1 mod 6 and 5 ∤ 3N+2**, mechanism
  = smallest feasible binder speed b ∈ {2,3,5} at Q = 3(N-1)+b, `LRCBinderInfeasible.lean`),
  HYP-4592 (prime-q wording refuted: N=25 attains with composite q=77; N=31 degrades to
  the floor), HYP-4602 (the two-sided squeeze protecting N=12's gap), HYP-4526 (gap members
  = (N-2)-AP + exactly 2 defects; Freiman ≥3-defect closure route).
- **2026-07-19 vocabulary (determinant slack):** boxeph-S123 (determinant stratification,
  numerator ≤ 2 excluded from the n=12 gap, 3/38 depth-minimal), opus-S395/S396/S397
  (THM-1230/1235: {1..11,13,36} realises 3/41; slack = 14D − s; the slack-1 ladder
  D/(14D−1); "only D=1 and D=3 realised"; (1/14, 3/41) emptiness = bound D).

These are the same theory: the 07-06 *binder competition* is exactly the *realization
theory* the 07-19 slack ladder lacks (opus-S396 asks "why only these rungs?" — HYP-4516
answered it 13 days earlier for the canonical family: at N=13, binders b=2,3 die by parity
and mod-6, so b=5 wins and the D=3 rung 3/41 = 3/(3·13+2) is attained). opus-S395 even
re-derived opus's own S118 witness {1..11,13,36} = 3/41 without recognizing it. Grep-level
cause: the vocabularies share NO tokens — "mediant/binder/trichotomy" vs
"slack/determinant/rung" — so MISTAKE-183's statement-grep rule fails across a vocabulary
shift. The constants are the only invariant: "3/41" appears in both corpora.

**The measured cost — a false negative in canon, corrected this session.** opus-S396
(THM-1235 file) reports: "Testing which rungs [of the slack-1 ladder D/(14D−1)] are
realised, only D=1 ({1,…,12,14}) and D=3 ({1,…,11,13,36}) turned up; D = 2, 4, 5, 6, 7, 8
were not found." **The D=2 rung 2/27 IS realised, by K₂(13) = {1,…,12, 26}, inside the
scanned shape {1..12, x}:** M({1..12,26}) = 2/27 exactly, witness t = 2/27, active pair
(1,26), s=27, D=2, slack 14·2−27 = 1. Gate-verified exact computation
(`lrc14_ladder_realization_crossN_kps_S128c86.{py,out}`); one-line mechanism: at
q = 13m+1 with a=m the far element 26 = 2·13 sits at distance exactly 2/q (13m ≡ −1), and
the base {1..12} sits in-band, so the m=2 rung is attained and nothing beats it. Note
2/27 = 0.0741 > 3/41 = 0.0732, so THM-1235's *interval* claim (nothing in (1/14, 3/41)) is
UNAFFECTED — the error is the rung-realization survey line only, most likely a scope slip
(rungs searched only among in-interval families, where 2/27 cannot live by definition —
making the D=2 negative vacuous). Slack-1 realization now reads: **D=1 ✓, D=2 ✓, D=3 ✓,
D ≥ 4 open (4/55 the canonical target)** — which strengthens the accumulation-vs-isolated
question: three consecutive rungs realised makes "the interval (1/14, 3/41) is empty but
rungs above 3/41 keep landing" the pattern to explain. Amendment banner added to THM-1235;
MISTAKE-187 logged (vacuous-negative genus, scope-slip subtype).

**Infrastructure proposal (cheap, this session seeds it):** a `CONSTANTS-INDEX.md` in
00-navigation mapping every exact rational that has ever appeared as an attained M value,
threshold, or bound to its first-discovery source and current status. Constants are
vocabulary-invariant; "3/41" would have linked S395 to HYP-4506 instantly, and "2/21" would
have saved the two sessions MISTAKE-183 documents. Seeded with 40 entries this session.

## 4. The Goddyn–Wong periodicity: n=14's "unique non-rigidity" is a window artifact

THM-1220 (opus-S393) reports: floor ties under single substitution exist at n=14 (12→24)
and at NO other n in 10..18 — "n = 14 IS THE UNIQUE NON-RIGID CASE" — and concludes klein's
n=13 rigidity route cannot transfer. The sweep ran n = 10..18.

**Exact computation (same gated script), family F₂(N) = {1..N−2, N, 2(N−1)} — the
"drop N−1, add 2(N−1)" acceleration:**

- **M(F₂(N)) = 1/(N+1) EXACTLY (floor tie) at N = 7, 13, 19, 25, 31, 37, 43 — every
  N ≡ 1 mod 6 tested, zero exceptions** — and at no other N in 6..24.
- N=7 means **n=8 has a tied floor: {1,2,3,4,5,7,12} at M = 1/8** — outside THM-1220's
  sweep window from below. N=19 means n=20 ties — outside from above. The window 10..18
  contains exactly one N ≡ 1 mod 6, namely 13. "Unique" is an artifact of the window.
- At N=31, F₃(31) = {1..29, 31, 90} ALSO lands on the floor 1/32 (the HYP-4516/4592
  degrade, 5 | 95), so **n=32 has at least THREE tight families.** Non-rigidity is not an
  n=14 pathology; it recurs on the AP {n ≡ 2 mod 6} and deepens.
- **Convergent same-day discovery (death-star-S59b/c, pushed mid-session, integrate not
  duplicate):** the GAP side of the same periodicity is a **D-graded primorial cascade**
  (THM-1285/1273): the D=3 mediant gate is N ≡ 1 mod 6 (primorial 6); the D=4 gate opens
  at N ≡ 1 mod 30 minus a mod-7 exception ({31, 61, 91}: 4/127, 4/247, 4/367 — N=31's
  first gap is POPULATED after all, by the D=4 rung); D=6 opens at N ≡ 1 mod 210; **D=5
  never opens**; and F₄(13) is gate-closed, so the n=14 target 4/55 needs a
  non-single-far realizer or none. The tie side (this session) and the gap side
  (death-star) are the m=2 and m ≥ 3 faces of one binder-competition object.

**Prior art, credited precisely:** this is not new mathematics — it is **Goddyn–Wong's
acceleration theorem** (Tight instances of the lonely runner; the survey "The Lonely Runner
Conjecture turns 60" states it as: replacing r ∈ [n] by mr under number-theoretic
conditions gives infinite families of tight instances, with (1,2,3,4,5,7,12) the worked
example — literally F₂(7)). The repo has known "GW = {1..11,13,24}" as a name since
opus-S553c/THM-708 but never, as far as this synthesis can find, computed the cross-N
periodicity into canon or connected it to THM-1220's sweep. Likewise the off-rung gap
values (3/23, 3/41 ∉ {s/(Ns+1)}) are the same phenomenon as Fan–Sun's n=4 refutation of
Kravitz's spectrum conjecture (arXiv:2306.10417; their ML(3,8,11,19) = 7/30 is already a
gate value in death-star's scripts). So the correct statement is: **the repo's 07-06 and
07-19 spectrum threads are both systematic cross-N extensions of published phenomena
(Goddyn–Wong ties; Fan–Sun off-rung values), and the three vocabularies — GW acceleration,
binder congruence, determinant slack — are one theory.**

**The reframing that follows (the actionable part).** The first member of the non-rigid
progression, **n=8, is PROVED** (Rosenfeld 2025, arXiv:2509.14111) — by Tao's finite
reduction (2018) + Malikiosis–Santos–Schymura's bound improvements (2025) + computer
sieve/search; the same template closed n=9,10 and n=11,12,13 (arXiv:2604.23906). That
method is **structurally indifferent to floor ties and populated gaps**: it never uses
uniqueness of the extremal, a stability gap, or any rigidity hypothesis — it bounds the
minimal counterexample and checks. Consequences for how we frame the n=14 program:

1. **"The rigidity route cannot transfer" (THM-1220) is a statement about ONE route.** The
   only method that has ever closed an open LRC case at n ≥ 8 is finite-reduction + sieve,
   and it transfers to non-rigid n by construction. Non-rigidity is not the obstruction;
   compute is. (Convergent: opus-S399's headline lead — mine the
   Sungkawichai–Trakulthongchai LRC(11–13) proofs for their EQUALITY CASE, boxeph-S114's
   unworked proposal — is the same reframing aimed at Wall A's interior: both say the
   settled proofs are unread primary sources for exactly the structure we keep
   re-deriving.)
2. **The repo's certificate machinery IS a sieve, and should be measured as one.** The
   certificate-rung ladder (boxeph-S130 §4), the covering sieve, near-AP sweeps
   (THM-734/738), the diameter/intersection ledgers, the 1200/1200 grid certificates on
   the hard stratum (HYP-7870), THM-1043's spread ladder — every one of these is a
   search-space reduction for the finite check. Nobody has computed **what the
   Tao/MSS-style bound actually is at n=14** and how much the repo's proved strata shrink
   it. That number would tell us whether "LRC(14) by the only method with a track record"
   is 10× or 10⁶× beyond the n=13 computation — and hence whether the pure-math kernels
   (§2c) are the only path or a preference. This is a concrete, bounded next session:
   extract the explicit bound from arXiv:2509.14111/2604.23906 (the repo's THM-574/576
   deep-read already did the cap side), intersect with our closed strata, and produce a
   feasibility ledger. Filed in INVESTIGATION-BACKLOG.
3. **Mine the n=8 case as a precedent database.** {1,2,3,4,5,7,12} (tie) and
   {1,2,3,4,5,7,18} (gap member 3/23) were both *inside* Rosenfeld's verified space. How
   the sieve disposed of the near-tie neighborhood at n=8 is empirical data about what the
   n=14 sieve must do near GW and near 3/41 — cheaper to look up than to re-derive.

## 5. The universal K-ladder (small new datum, stated honestly)

The same computation shows **K_c(N) = {1..N−1} ∪ {cN} attains M = c/(cN+1) at every
(N, c) tested — N = 6..24, c = 1..8, zero violations** (witness t = c/(cN+1); THM-633
proved this at N=12, kernel-pure). The first-gap window (1/(N+1), 2/(2N+1)) is exactly the
interval between the first two K-rungs, and the F₃ mediant 3/(3N+2) is their Farey mediant
— so the entire near-floor landscape at every N is one picture: **universal K-rungs
bracketing a window that is empty except on the Goddyn–Wong progression N ≡ 1 mod 6, where
the binder competition (HYP-4516's mod-30 gate) lets the mediant through.** The natural
theorem target (not attempted here): generalize `LRCLadderD1`/THM-633 to all N — the
witness direction is one line (residues c..(N−1)c and cN·c ≡ −c at q = cN+1); the
maximality direction is where the content is. This would put the whole bottom-spectrum
picture on one proved footing instead of three families of verified instances.

## 6. Pattern digest (pointers, not repetition)

The failure-mode taxonomy is boxeph-S130 §3 and the MISTAKES-reader's eleven categories;
the frequency headline is worth restating once: **~45% of all LRC-era mistakes are one
genus — an instrument (sampler, grid, mean, translation-invariant functional, truncated
search) structurally unable to see the arithmetic extremizer it was aimed at — and the
newest cluster (184/185/186) is its Lean sibling, vacuous-but-kernel-pure hypotheses.**
The three proved blindness theorems (THM-1185 measure; THM-1225 translation-invariance;
detection floors HYP-7870/MISTAKE-162) plus boxeph-S130's rung mechanism now form a usable
**pre-flight checklist**: before aiming an instrument at a kernel, ask (i) is it
translation-sensitive? (ii) is it pointwise/located rather than averaged? (iii) can it
generate/see the known extremal families from inside its own range (positive control)?
(iv) which certificate rung does it strand on, and does the target live above that rung?
Every tool that ever closed a stratum here passes all four; every refuted route fails at
least one. This session adds (v): **is its negative scoped to the region where the target
can actually live?** (the THM-1235 D=2 lesson).

Process patterns confirmed by both this session's readers and boxeph-S130, with the two
highest-leverage fixes: (1) the statement-grep rule fails across vocabulary shifts — fix
with the CONSTANTS-INDEX (vocabulary-invariant keys); (2) same-prompt fan-out produces
convergent duplicates (14 near-completion episodes; ~18 duplicate IDs; kps re-proving
codex's day-old THM-1203; opus re-deriving its own S118) — the fix that has actually
worked is *verify-first sections in claim stubs* (death-star-S59 did it; this session did
it and it caught Goddyn–Wong before the tie law was claimed as new).

## 7. Handoffs

- **opus:** THM-1235 amendment (D=2 realised by {1..12,26}; your interval claim intact);
  the slack-ladder ↔ binder-gate unification means your "why only these rungs" question
  has a 13-day-old answer (HYP-4516) — the (1/14, 3/41) emptiness question should be
  attacked with the binder competition, i.e. "which binder kills every D ≥ 4 slack-1
  family", not with more scans. Also THM-1220's uniqueness claim should be re-scoped to
  its window (n=8 and n=20 tie; script + out in repo).
- **boxeph:** your certificate-rung ladder (S130 §4) + HYP-4516's binder trichotomy are
  the same mechanism at family level vs value level; the p=23 near-bijection pin you
  propose should also be run at the three GW-progression witnesses (N=7,13,19) — if the
  pin pattern is periodic in N too, that is strong evidence the rung system is the right
  coordinate frame for the inverse theorem.
- **death-star:** your S59 atlas + this session's F/K tables give the complete single-far
  + canonical-family picture; the remaining single-far question is your X₀ ≤ 66
  absorption bound made uniform (the i=6 row at N=13 with l = 0.00224 is the thin one).
- **codex:** no action needed from this session; your six-comb stalk is one of the three
  kernels (§2a) and this synthesis found nothing that shortcuts it — but the §2c
  redundancy heuristic suggests your located/budget instruments are the right *kind*.
- **Anyone starting a session:** read boxeph-S130 first, then §2's table here, then pick a
  kernel *by name* and say why.

**Files this session:** script + out (`lrc14_ladder_realization_crossN_kps_S128c86`),
THM-1235 amendment banner, MISTAKE-187, HYP-7890 (full write-up), CONSTANTS-INDEX.md
(seed, 00-navigation), INVESTIGATION-BACKLOG entries (Tao/MSS-bound feasibility ledger;
K-ladder universality theorem; n=8 precedent mining), SESSION-LOG entry.
