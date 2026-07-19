# CONSTANTS-INDEX — exact rationals ↔ first source (LRC thread)

**Purpose (kind-pasteur-S128c86, HYP-7890).** Vocabulary drifts; constants don't. The
statement-grep rule (MISTAKE-183) fails across vocabulary shifts (binder/mediant vs
slack/rung cost a canon error, MISTAKE-188; "2/21" cost two duplicated sessions;
"Hamming" vs "single-far" cost 30 minutes, MISTAKE-187). **Before working ANY exact
value — attained M, threshold, bound, measure — grep THIS file for it.** If your value
is here, read the source before deriving anything. If it isn't, add it with your session
ID when you freeze it. One line per constant: value — what it is — first source /
proof status.

## Attained loneliness values (13 speeds, n=14 unless marked)

- `1/14` — the floor; attained by {1..13} (AP) and {1..11,13,24} (Goddyn–Wong; =GW; THM-1115/1120; F₂(13) acceleration). Tight locus = exactly these two up to dilation (THM-1120).
- `3/41` — first value above the floor known; {1..11,13,36}; the N=13 mediant 3/(3N+2). First: opus-S118/HYP-4506 (2026-07-06); re-derived opus-S395/THM-1230 (07-19). Gap (1/14, 3/41) empty in all searches (THM-1235/1240, ~12,400 families; NOT a theorem).
- `2/27` — slack-1 D=2 rung; ATTAINED by {1..12,26} = K₂(13) (kind-pasteur-S128c86, MISTAKE-188 correction to THM-1235). Also the Farey neighbour of 1/14.
- `3/40` — k=3 Kravitz rung s/(13s+1); attained by {1..12,39} = K₃(13) (S128c86 table).
- `4/53` — attained by {1..11,13,48} = F₄(13) (opus-S395 ladder m/(12m+5)) and K₄(13)={1..12,52} (S128c86).
- `1/13` — {1..12,14}, {1..12,15} etc. (slack-1 D=1); also the 12-speed floor (below).
- `14/183` — covering minimum; the deep well {1..12,182}; 183 = Φ₆(14). THM-724/726 (07-11..13); unique covering strict-interior family (death-star S58f).
- `8/105` — covering-2..13 (not 14) witness populating the (1/14,1/13) covering stratum (MISTAKE-161b, death-star-S57).

## Attained values (12 speeds, n=13 core / Tao n=12 uniqueness)

- `1/13` — floor; unique primitive attainer {1..12} (verified to {1..20}, boxeph-S121/S122; Tao n=12 uniqueness = CRUX (C) still open in general).
- `2/25` — second value; {1..11,24}; klein Hamming-1 rigidity equality case (THM-1004); k=2 Kravitz rung; `LRCMod25Floor`/`LRCLadderD1` green.
- `3/37` — {1..11,36} = K₃(12) (THM-633 ladder m/(12m+1), kernel-pure).
- `3/38` — the DEPTH-MINIMAL open gap value for (1/13, 2/25); mediant; NOT attained anywhere tested (kps-S12 1.5M bases in [1,26]; boxeph-S123/S124/S125 obstruction mapped: cross-modulus needle-covering; open).
- `2/23` — tightest covering competitor {1..13}\{6} (boxeph-S122, covering-rigidity margin 3/299); also S* = {1,2,3,5,7,…,13,38,42} rescue value via mod-23 (boxeph-S130 §4).
- `1/156` — raw spectrum gap 1/12 − 1/13 (boxeph-S121).

## Cross-N landmarks

- `3/23` — N=7 (n=8) first-gap member {1,2,3,4,5,7,18} = F₃(7); mediant 3/(3·7+2). HYP-4506; single-far-unique (death-star-S59 atlas).
- `5/33` — N=6 first-gap member {1,5,6,11,16,17}; bordered dilated AP, order 3 (HYP-4506/4602).
- `1/8` tie — N=7 floor tie {1,2,3,4,5,7,12} = F₂(7) = Goddyn–Wong's worked example; n=8 is NON-RIGID and PROVED (Rosenfeld). S128c86.
- `1/20`, `1/26`, `1/32`, `1/38`, `1/44` ties — F₂(N) ties the floor at every N ≡ 1 mod 6 tested (N = 19, 25, 31, 37, 43); N=31 also F₃ → THREE tight families (HYP-4516 5|95 degrade). S128c86.
- `7/30` — Fan–Sun n=4 off-rung value ML(3,8,11,19); refutes Kravitz spectrum conjecture at n=4 (arXiv:2306.10417); used as evaluator gate.
- `c/(cN+1)` — the universal K-ladder {1..N−1}∪{cN}; attained at every (N,c) tested, N ≤ 24, c ≤ 8 (THM-633 proved at N=12; S128c86 data elsewhere).
- `3/(3N+2)` — the F₃(N) mediant; attained ⟺ N ≡ 1 mod 6 and 5 ∤ 3N+2 (HYP-4516 mod-30 binder gate; exception N=31 degrades to floor, HYP-4592).

## Thresholds, measures, bounds (LRC(14) machinery)

- `2/21` — the continuum four-comb BAD-measure ceiling, equality iff 4-term AP (THM-1203, codex-S77, Lean-certified). Re-proved redundantly by two agents (MISTAKE-183).
- `1/343` — |B| the bad-set/slack-simplex volume (codex THM-1183; grid artifact 0.003367 corrected in MISTAKE-176).
- `477/1078` — μ_{1/7}({1..13}) exact (opus-S130 three-gap; = rhoGlobFloorRat(13)).
- `3/4` — the k=8 leg floor μ_{1/7}(E₈) ≥ 3/4 (THM-651 shifted-tent, kind-pasteur-S73).
- `m_P = 0.0565…` — quantitative bar for hlarge; honest T_k bars in MISTAKE-123.
- `31807/194040` — min |S(P)| over 495 eight-cores, unique at (1,2,3,5,7,8,9,11) (THM-1251; replaces measured "0.164").
- `1/70` — shortest longest-component ℓ(P) floor, unique core (1,2,6,7,8,9,10,11) (THM-1123/1251).
- `6/(7k1)` — complete k1-gap width (eroded start complex; THM-1251: erosion threshold k1 ≥ 204, bottleneck core (1,2,3,5,7,8,9,11), K* = 203.865).
- `2/29` — pigeonhole prime floor (no multiple of 29 ⟹ M ≥ 2/29); THM-518; the "last mile" window [2/29, 1/14) has width 1/406 (MISTAKE-093).
- `2/19` — largest prime-certificate rung clearing the whole n=12 gap regime (mod-19 spread lemma, boxeph-S126..S129, `gap_regime_mod19_spread` in LRC14Ledger).
- `15/154` — seven-speed Hunter/safe-set floor (THM-1221, codex-S82; crosses 2/21 by 1/462).
- `1/(2ρ)` — clustered floor M ≥ 1/(2ρ), ρ = v_max/v_min (death-star-S58h; M < 1/13 forces ρ > 6.5).
- `1/(49c)` — five-comb survivor floor in every c-slow gap (THM-1198, codex-S76).
- `4/55` — the canonical D=4 slack-1 target; unique least-denominator fraction in (1/14, 3/41) (THM-1260, ex-opus-1240). OPEN; F₄(13) is binder-gate-closed so any realizer is non-single-far (death-star-S59b, backlog vii).
- `D = M·s` — the determinant identity at the active pair (THM-1261, ex-opus-1245); LRC(14) ⟺ s ≤ 14D (THM-1205).
- `4/127`, `4/247`, `4/367` — D=4 slack-1 rungs ATTAINED at N = 31, 61, 91 ({1..29,31,120} etc.; binder 7); the D-graded primorial cascade: D=3 gate mod 6, D=4 gate mod 30 (minus mod-7 exception), D=6 opens at N ≡ 1 mod 210 (6/1271 at N=211), D=5 NEVER (death-star THM-1256/1257, HYP-7900).

**Maintenance rule:** add on freeze; correct in place with your session ID; never delete —
strike through with the correction pointer.

## Tower rung values (canonical-D single-far gap members; death-star S59b-S59e)
- 4/127 — M(F_4(31)) = M({1..29,31,120}), D=4 slack-1 rung, witness 55/127 — THM-1256 (proof-backed exact via THM-1258)
- 4/247 — M(F_4(61)) = M({1..59,61,240}), witness 70/247 — THM-1257 (predicted-then-found)
- 4/367 — M(F_4(91)) = M({1..89,91,360}), witness 53/367 — THM-1257
- 6/1271 — M(F_6(211)) = M({1..209,211,1260}), witness 115/1271 — THM-1257 (out-of-sample #2)
- 7/16183 — M(F_7(2311)) = M({1..2309,2311,16170}), witness 12449/16183 — THM-1258 (out-of-sample #3; window 9.4e-8)
- 9/270287 — M(F_9(30031)) = M({1..30029,30031,270270}), witness 31799/270287 — THM-1270 (out-of-sample #4; window 5.5e-10)
- 1/212 — M(F_5(211)) degrade-to-floor (composite binder 9); a tight family at N=211 — THM-1257
- 10/5105119 — PREDICTED M(F_10(510511)) (binder 19; next rung, unrun) — THM-1270 lead
