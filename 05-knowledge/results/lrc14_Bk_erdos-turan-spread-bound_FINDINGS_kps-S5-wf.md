# LRC(14) S3: the erdos-turan-spread-bound angle -- FINDINGS (kps-2026-06-18-S5-wf)

Angle: prove the spread bound `B(k)`: `spread(E) > B(k) => mu(E) >= mu_min^bdd(k)`, the missing
piece of the uniform floor lemma (THM-527-G / THM-528-E / OPEN-Q-108).
`mu(E) = meas{x in [0,1): maxgap{frac(e_i x): e in E} > thr}`, `0 in E`, `k=|E|`.

## 0. A RIGOROUS exact mu engine (the prompt's tool was NOT rigorous)
`lrc14_Bk_erdos-turan-spread-bound_kps-S5-wf.py` computes mu EXACTLY: on each order-cell the
cyclic gaps are AFFINE in x, so `{maxgap>thr} = UNION_gaps {affine_gap > thr}` (exact, no missed
gap=thr crossings). VERIFIED vs canon (mu_3..mu_6, mu_13) and vs numeric over 400 random E
(max err 3e-5 ~ sampling noise). CORRECTION: canon THM-527-C lists "mu_7 = 13/35"; that is the
PERFORATED minimizer mu({0,2,3,4,5,6,8})=13/35, NOT consecutive. The CONSECUTIVE value is
mu_7 = 83/210 = 0.39524 (numeric-confirmed). (THM-527-F already uses 0.395 as the consecutive value.)

## 1. The EXACT Weyl deviation identity (PROVED; independently re-derived = concurrent kps-S5)
Since all e_i are integers the orbit `x -> (e_i x mod 1)` is a closed curve, so by Parseval
> **mu(E) = F(k) + sum_{m != 0, m.e = 0} ghat(m)**,
`F(k)=ghat(0)=` iid ceiling `= sum_{j>=1}(-1)^{j+1} C(k,j)(1-2j/7)_+^{k-1}` (F(13)=3132376013/13841287201).
The deviation is the Fourier mass of `g` on the RELATION LATTICE `Lambda(E)={m: m.e=0}`.
VERIFIED: structured/dense E -> large negative dev (consec/perforated dev ~ -0.40); Sidon-like
large-spread E -> mu ~ F(k) (k=7 near-Sidon mu ~ 0.83-0.86 ~ F(7)=0.80). So small mu <=> many
SHORT relation vectors <=> dense bounded-spread ruler. (Concurrent kps-S5 pushed this to a smooth
minorant G with explicit kernel psi-hat; same lattice structure.)

## 2. The 2/7 measure has NO clean spread bound -- minimizer drifts to spread ~2-3k, mu keeps dropping
EXACT per-spread minima (mu_exact, exhaustive where noted):
| k | mu_min^bdd(2/7) | minimizer spread | note |
|---|---|---|---|
| 4 | 19/21 = .9048 | 3 (=k-1) | consecutive |
| 5 | 9/14 = .6429 | 4 | consecutive |
| 6 | 4/7 = .5714 | 5 | consecutive |
| 7 | 13/35 = .3714 | 8 (~1.1k) | perforated {0,2,3,4,5,6,8}, exhaustive to s=21 |
| 8 | 71/220 = .3227 | 11 (~1.4k) | perforated |
| 9 | **811/4095 = .1980** | **18 (~2.0k)** | symmetric perforated, EXHAUSTIVE s=18,19,20,21 |
CORRECTION TO CONCURRENT kps-S5: their "bounded-spread min" search (capped spread ~14) recorded
k=9: 164/735=.2231 (s12), k=10: 468/2695=.1737 (s14), k=11: 409/2548=.1605 (s14). ALL are
UNDERCUT at spread ~2k: k=9 -> 811/4095=.1980 (s18, exhaustive); k=10 -> ~.1326 (s20);
k=11 -> ~.1214 (s20). The true minimizer spread is ~2k, NOT ~1.3k. Their "stretch-adversary STABLE"
verdict is an artifact of insufficient spread; their own `true_floor_k13` (correctly) shows the drop.
Annealing trajectories (validated on k=9 exact: annealer hits the exhaustive minimizer at s18):
 - k=9 (DONE, full U-shape): .283(s8) .223(s12) **.198(s18=MIN)** .259(s27) .270(s36) .312(s54)
   .314(s72). FLOORS at .198 at spread ~2k, then RISES back toward F(9)=.569. So the 2/7 floor
   IS POSITIVE for k=9, at BOUNDED spread ~2k.
 - k=13: .179(s12) .132(s16) .084(s26) .068(s39) .060(s52) -- still falling at 4k; by analogy with
   k=9 the U-min is expected at spread ~2.5-3.5k (=~33-45) but I did not reach the turnaround
   exactly. The 2/7 floor for k=13 is thus presumably positive but at LARGER spread (~3k=39) and
   LOWER value (~.06) than canon's "spread <= ~30 / mu_min ~ .11" estimate. The minimizer-spread
   grows ~2.5k, exceeding the canon "O(k) <= ~30" claim for k=13 (where ~2.5k ~ 33).

## 3. THE 2/7 MEASURE IS THE WRONG OBJECT (integrating concurrent HYP-2592 / THM-530)
The via-max criterion (cluster gap > 2/7) is SUFFICIENT but NOT necessary; `rho*_{2/7}=0` on
explicit admissible families (concurrent, exact), yet those sets are LONELY. The object the
reduction ACTUALLY needs is the **global-witness 1/7 measure** `mu_{1/7}(E)=meas{x:maxgap>1/7}`.

## 4. THE 1/7 MEASURE HAS A CLEAN SPREAD BOUND -- CONSECUTIVE MINIMIZES (the deliverable)
- **k <= 7: mu_{1/7}(E) = 1 for ALL E (PROVED, pigeonhole).** k<=6: 6 cyclic gaps sum to 1 so
  maxgap >= 1/6 > 1/7 always (identically 1). k=7: =1 a.e. (the only way all 7 gaps <=1/7 is all
  EXACTLY 1/7, a measure-zero x-set). Verified exhaustively (distinct value = 1.0 over all shapes).
- **k = 8, 9: CONSECUTIVE is the GLOBAL minimizer of mu_{1/7} (EXHAUSTIVE over all spreads to ~2k).**
  Every larger spread strictly RAISES mu_{1/7}. So `B_{1/7}(k) = k-1` (trivial spread bound).
- **k = 13: consecutive robustly minimizes (adversarial check).** The 2/7-minimizer perforated
  shapes give mu_{1/7} ~ 0.95-0.99 (FAR above consecutive 477/1078=.4425); no single-drop or
  random spread-2k..4k shape beats consecutive. THE TWO MEASURES HAVE OPPOSITE EXTREMIZERS:
  perforation HELPS create 2/7-gaps (lowers mu_{2/7}) but HURTS 1/7-netting (raises mu_{1/7}).
- Consecutive mu_{1/7} values (EXACT, cross-validated vs concurrent mu17_crudefloor, ALL MATCH):
  k=8 691/735, k=9 247/294, k=10 38/49, k=11 1381/2205, k=12 13823/24255, k=13 477/1078=.4425.

## 5. NET (honest)
- PROVED: the rigorous mu engine; the Weyl deviation identity; mu_{1/7}=1 for k<=7 (pigeonhole).
- VERIFIED (exact/exhaustive where stated): consecutive minimizes mu_{1/7} (k=8,9 exhaustive;
  k=13 adversarial); the 2/7 bounded-spread minima + the ~2k-minimizer-spread correction to
  concurrent kps-S5.
- OPEN (the residual, now ISOLATED and EASIER): a PROOF that consecutive minimizes mu_{1/7} for
  all 8<=k<=13 (then mu_{1/7}(E) >= mu_{1/7}({0..k-1}) >= thr_k is a single explicit check, and
  with the concurrent union bound this closes the k>=8 branch of the global-witness floor). This
  is the WEAK 1/7 spread bound flagged OPEN by HYP-2592 -- strictly easier than the 2/7 B(k),
  which (per section 2) likely does NOT have a clean spread bound at all. LRC(14) NOT proved.
