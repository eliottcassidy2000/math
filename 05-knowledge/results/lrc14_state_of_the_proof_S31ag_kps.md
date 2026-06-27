# State of the LRC(14) proof — integrated with arXiv:2604.23906 (kps-2026-06-27-S31ag)

*A comprehensive current-state synthesis answering the owner's question: where does the 14-runner
attack stand, and how does it relate to the recently-proven 11–13 runner cases? Integrates the paper
(Sungkawichai–Trakulthongchai, LRC k≤12) with the project's covering-bound / witness / gK8 machinery and
this session's findings (HYP-3087–3090, THM-574, THM-576). Supersedes the stale 2026-06-21 honest status
doc on the parts noted. Reconciled in real time with mac-mini-S60/S61 and codex-S255.*

## 1. The landscape (where 14 sits)
Recent full LRC results, by #runners: **7** (Barajas–Serra 2008), **8** (Rosenfeld 2509.14111),
**9** (2512.01912), **9–10** (2511.22427), **11,12,13** (Sungkawichai–Trakulthongchai 2604.23906,
= LRC k≤12 speeds). **14 runners (k=13 speeds) is the FIRST OPEN case.**

**Why 14 is the wall — one fact, two framings.** The paper's polynomial method (Prop 4.1) runs over the
field `ℤ_{k+1}` and needs **k+1 an odd prime** (Fermat). For 14, `k+1 = 14 = 2·7` is composite: `ℤ_14`
is not a field. This is *identical* to the project's apex prime 7 / `14 = 2·7` wall (HYP-3087):
`φ(14) = 6` units can't fill the 13 indices; the minimal null polynomial `b_7 = ∏_{j=0}^6(X−j)` (degree 7,
`14|7!`) kills control above moment-order 6. The paper's composite-`k+1` fallback = **lifts at the prime
factors c=2,7** = the project's descent **14→7→2**; **THM-573 (level-7 sieve) IS the c=7 lift** (THM-574:
c=7 is the unique optimal lift since `1/7` = forbidden-arc width = finest one-survivor spacing), using
**LRC(≤13) = the paper as induction base.** So our proof is genuinely an inductive layer on top of theirs.

## 2. The proof tree (what is proved / open)
```
LRC(14)  ==WLOG primitive (dilation inv.)==>
 ├─ 14-free (no v≡0 mod 14)      t=1/14         DONE (THM-523, trivial; = empty coarse witness set)
 └─ COVERING (some v≡0 mod 14)   = the whole problem (mac-mini-S59 redirect; bound is TIGHT, HYP-3084)
     ├─ ≥7 multiples of 7        level-7 sieve  DONE (THM-573 ⊇ THM-571; the c=7 lift)
     └─ ≤6 multiples of 7  ── RESIDUAL ──:
         ├─ large speed present  Node-3 peel    r≤6 RIGOROUS; r≥7 = effective Erdős–Turán (OPEN)
         └─ bounded core         CRUX 1         the measure floor / coverage extremality (OPEN)
```
Reductions PROVED/exact: covering localization (S59); the c=7 lift family (THM-574); the bridge
`I(13,7,1) = covering mod 7` with the c=2 rescue (mac-mini-S61, exact); WLOG primitive; the bounded-speed
reduction (paper Lemma 2.6, `∏u_i < B_k = (C(k+1,2)^{k-1}/k)^k` — note `C(k+1,2)` again).

**The single remaining crux is CRUX 1** (the bounded-core measure floor), with Node-3's r≥7 tail as its
unbounded twin (the peel reducing unbounded→bounded). Everything else is in hand.

## 3. CRUX 1, fully localized this session
The cover bound is `meas(S7(E)) = P(N=0) ≤ cap_k`, `N = #empty inner sectors ∈ {0..6}`, over bounded
covering clusters `E` of size `k`. Two sides:

**The cap (RHS) is SOLVED (THM-576).** `cap_k = min_{|P|=13-k} meas(lonely(P)) = C(14-j,2)/C(14,2)`
(`j=13-k`), the **pairwise avoidance probability**, EXACT for `j ≤ 3` (k≥10), minimizers `{1}∪top`
(j=2 PROVED elementarily, `min meas(D_p∩D_q)=1/91` at `{1,13}`). The binding rows dip by exactly
`−1/4004` (k=9) and `−1081/76440` (k=8); the minimizer pattern BREAKS at j=5 to `{1,5,7,8,9}`. So the
**entire cap side is the clean triangular `C(k+1,2)/91` for k≥10 plus two explicit constants** — the
`C(14,2)=91` apex = #pairs of 14 runners (triangle foundation).

**The extremality (LHS) is the crux.** `consec/AP maximizes P(N=0)` (coverage extremality = low-discrepancy
/ three-gap). VERIFIED robust this session (0/2006 beat consec at k=8,9,10, with the e=0 anchor enforced),
with margins to cap `0.054/0.078/0.100` (k=8 tightest). The proof object: mac-mini's gK8 low-order
moment-LP, whose pairwise block `M` (the 6×6 sector co-emptiness) is **reflection-symmetric (= complement
`T→T^op`) dominant-Perron, NOT 4I+2J** (HYP-3089 — settles the certificate route). The difficulty is
genuinely **order-3** (S3 closes k=9,10; S4 only k=8) — pairwise/Paley–Zygmund FAILS at consec (Lead 9).
This pins the transition: pairwise-exact through j≤3, order-3 onset at j=4→5.

## 4. The unification (paper ↔ project), and a caution
- **The project IS the analytic substitute for the paper's enumeration** (mac-mini-S61): the paper
  enumerates `I(k,p,1)` at cost `p^{(k+1)/2}/(k2^k)` (k=13 prohibitive); `p0(E)≤cap_k` (the uniform-in-p
  continuous limit) replaces it with a measure bound.
- **Conjecture 7.1(13) ⟺ LRC(14)** (it is a witness-time statement = the witness route). But the
  **literal "∀d≥D" form is FALSE** (codex-S255 THM-575; aliasing `d|V` kills it, = HYP-2866). The uniform
  object is the lonely **measure**, not a denominator/arc.
- **CAUTION on "three walls = one constant ~200":** the *normalized/peeled* witness arc is bounded to
  `V*~234`, but the **RAW** largest lonely arc of an apex-`14V` set is `≤ 3/(49V)` EXACT (apex-comb
  limited), so it decays from `V~15` — a *different, smaller* constant. The peeling (Node-3) is the
  reduction between them. The unification holds in the normalized frame only.

## 5. The cleanest path to a proof (recommended focus)
1. **Prove the coverage extremality "consec maximizes P(N=0)" for k=8,9,10** — equivalently consec
   maximizes the order-≤4 gK8 dual — via the reflection-Perron (complement-symmetric) block structure +
   three-gap (the AP is the min-discrepancy orbit). This is THE crux; it is order-3 (k=9,10) / order-4
   (k=8); the cap RHS is already `C(14-j,2)/91` (THM-576). Margins 0.054–0.100 give room.
2. **Node-3 effective Erdős–Turán** for the r≥7 unbounded peel (the only non-finite link besides 1).
3. The induction base LRC(≤13) is the paper — DONE.

**Honest bottom line:** LRC(14) is NOT proved. But it is now a single coverage-extremality (consec
maximizes 7-sector coverage, order-3, reflection-Perron) over bounded clusters, with the cap side solved
(pairwise triangular + 2 constants), the reductions matching the paper's own program, and the constant
`D = V*` pinned to a few hundred. The endgame is the smallest and best-understood it has ever been.

→ HYP-3087/3088/3089/3090, THM-573/574/576, THM-575 (codex), THM-565, mac-mini-S59/S60/S61, codex-S255,
OPEN-Q-108, arXiv:2604.23906, the triangle foundation, [[lrc14-thread]].
