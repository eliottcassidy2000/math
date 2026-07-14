---
id: HYP-6830
title: Scale uniformization of Claim B — proved sheet regime, refuted raw-fragmentation complementarity, and a peel-relative splice
status: PARTIALLY PROVED / REFUTED AS STATED — THM-761 proves the large-scale regime; raw fragmentation is not controlled by divisor-packet scale; a peel-relative splice and the r>=7 tiling residue remain open
source: opus-2026-07-14-S299
depends_on:
  - THM-755   # capped envelope, v* = r_P / (pi |G'_P|)
  - THM-760   # r=1 coprime sheet dodge
  - THM-761   # multi-exception sheet covering bound (this session)
  - HYP-6780  # v*(cP) = c v*(P): the scale covariance that killed raw-height bands
related: [THM-756, THM-757, THM-758, HYP-6785, HYP-6815, HYP-6820, HYP-6835]
---

# HYP-6830 — scale uniformization of the ≥4-far covering endgame

## The proposed splice, with its failed bridge isolated

For every covering 13-speed family V and every scale c ≥ 2, write V = cP ⊔ W
(canonical: P = {v/c : c|v}, r = |W|). Define c*(V) = the largest c with |P| ≥ 7
(equivalently r ≤ 6), or 1 if none. Then:

1. **Large scale (PROVED, THM-761):** if c*(V) ≥ 43 (or the exact per-(r,c)
   criterion fires, or the gcd budget Σ g_a(⌊c/(7g_a)⌋+1) ≤ c−1 holds), V is closed:
   M(V) ≥ 1/14, by a free witness sheet. No enumeration, no witness search.
2. **Small scale (OPEN — coordinate incomplete):** if c*(V) ≤ 42, then every dilated
   sub-structure inflates the capped-envelope band edge by at most the SAME bounded
   factor (`v*(cP) = c·v*(P) ≤ 42·v*(P)`, HYP-6780 used positively). This local
   covariance remains exact, but `c*(V)≤42` does not by itself bound the core's
   fragmentation or produce a globally bounded normalized domain. TO PROVE, with a
   richer state: every ≥4-far covering family with
   c*(V) ≤ 42 either safe-peels (THM-753), fires the capped envelope, or lies in an
   explicitly bounded normalized band family. This is the corrected, quantifier-honest
   replacement for THM-758 Claim B's "finite check" and for the REFUTED q≤25 finish
   (codex-S3, family 26·{1..12} ∪ {339}, first witness q=27).
3. **Raw fragmentation/divisibility complementarity (REFUTED):** there is no
   function `B` with `r_P≤B(c*(P))` for all twelve-cores. A scale-free family with
   `c*(P)=1` and `r_P→∞` is given below. Dilation is one fragmentation mechanism,
   but a single high-frequency runner cuts a fixed safe interval into arbitrarily
   many pieces without creating a seven-runner divisor packet.

## Exact infinite-family falsifier

For every prime `N>11`, put

```text
P_N = {1,2,...,11,N},
V_N = P_N union {1092 N},       1092 = lcm(12,13,14).
```

This family has all the scope properties needed to refute the proposed bridge.

1. **No high-support scale.** Among `{1,...,11}`, at most five entries share any
   nontrivial divisor (the five even entries attain the maximum). The prime `N`
   shares no divisor with them, and adjoining `1092N` can add at most one member to
   any old divisor packet. Hence the largest divisor packet has size five in `P_N`
   and six in `V_N`; in particular `c*(P_N)=c*(V_N)=1`.
2. **Unbounded exact fragmentation.** The base good set for `{1,...,11}` contains

   ```text
   J = [1/14, 13/154],             |J| = 1/77.
   ```

   Inside `J`, the `N`-runner removes the disjoint open danger teeth centered at
   `k/N` with radius `1/(14N)`. The number of full teeth in `J` is `N/77+O(1)`;
   every full tooth separates two components of `G'_{P_N}`. Thus `r_{P_N}` is
   unbounded although `c*(P_N)=1`. Exact interval arithmetic gives component counts
   `18,22,38,72` for `N=101,211,503,1009`.
3. **It occurs inside the covering endgame.** Every `V_N` is primitive because it
   contains `1`. Speeds `2,...,11` carry moduli `2,...,11`, while `1092N` carries
   `12,13,14`; hence `V_N` is covering. This is not an LRC counterexample: the
   top peel fires THM-755's capped-envelope test in every audited instance. It is a
   counterexample specifically to using `c*` as a sufficient fragmentation
   coordinate.

The surviving quantity is therefore peel-relative, not component-count absolute.
For a proposed peel `v`, the cap sees

```text
kappa(P;v) = r_P / (v |G'_P|) = pi v*(P) / v,
```

and not `r_P` alone. A replacement splice should retain at least this normalized
load, the divisor-support profile `c -> #{p in P:c|p}`, and enough endpoint-owner
data to distinguish dilation copies from high-frequency tooth insertion. It remains
open whether those data admit a finite or recursively compact quotient; no
replacement theorem is claimed here.

**Ratio study (opus-S300, independent confirmation + the measured constants).** The
stress battery `lrc14_regime2_complementarity_stress_opus_S300.py` (+ .out)
independently refuted `r_P <= B(c*)` empirically (scale-free cores: median r_P
doubles with height, 60 -> 1700 over H = 50 -> 1600) and measured the peel-relative
invariant at the first peel, `rho(P) = v*(P)/max(P) = kappa(P; maxP)/pi`:

- `rho` is **scale-invariant on dilates** (9.334... at every `c` for `c*{1..12}`) and
  `O(1)` on every tested family: spread scale-free 2.0–3.0, partial dilates 1.9,
  deep-well shape 1.22, GW 12-core 8.06;
- an adversarial hill-climb over `c* <= 42` cores **converged back to {1..12}**:
  measured max `rho = 9.335` at the interval shape itself — the same extremal that
  carries the covering-min and the H-band corners;
- the codex falsifier `P_N = {1..11, N}` has `rho < 1` for large `N` (band EMPTY —
  its top peel indeed fires the envelope), so it kills the raw-`r_P` bridge while
  CONFIRMING the ratio coordinate;
- since `r_P <= Sum(P)` always, `rho <= 12/(pi |G'_P|)`: the remaining proof
  obligation is a **|G'| floor off the classified tight families** (the mac-mini B5
  stability lane) — that alone converts `rho = O(1)` from measured to proved and
  makes the regime-2 band domain bounded in NORMALIZED (peel-relative) coordinates.

## The residual inside regime 1 (named by THM-761; r=7 stratum CLOSED by THM-767)

- **r = 7: CLOSED above the shape bound (THM-767, opus-S300 — the event pierce).**
  For `7|c`-compatible strata the bad-sheet counts are CONSTANT (= c/7) off a finite
  event set; at every event moment in the closed core-safe set the total count drops
  to `c-1`, freeing a sheet with closed clearance `>= 1/14`. Route (a) below was the
  right one, sharpened: not the whole safe set — the SWITCHING TIMES. The realized
  c=7 wall instance is pierced at all 203 tested core-safe event moments. Explicit
  sufficient condition: `w_max > 7*Sum(P)`. Remaining inside r >= 7: the r >= 8
  alignment residue (single-event pierce fails structurally; realized), `7 ∤ c`
  decks (s-threshold), and the bounded-shape residue `w_max <= 7*Sum(P)`.
- The **KCL absorption law** (THM-767(4)) replaces the Newman-route sketch in (b):
  maintained exact tilings need `sum over mirror partners (14*gcd | w_a+w_b) of
  gcd(w_a,w_b) >= w_a` per exception — violated by 1399/1400 random 7-sets; the
  KCL-feasible packets are the rigid arithmetic corners of the deck.
- **Deep gcd entanglement** (Σ g_a over budget): recursive descent c → c/g;
  bookkeeping to write; termination is clear (c strictly decreases).

## The recursion (the structural content, from the S299 reflection)

The sheet residual at scale c IS an inhomogeneous discrete lonely-runner instance on
Z_c: runners = exceptional residues w_a mod c, offsets = w_a·t0/c, radius 1/14,
lonely time = free sheet. The tight case (r = 7, arcs tile) is the 7-clock partition
(THM-754) one level down: tight instances are tilings at every level. The underlying
object is a pointed circle with burned arcs, self-similar under scale descent; the
free sheet is the next level's basepoint (the observer-lens principle survives
descent).

## Probes filed (testable, not claimed)

- **FI cubic probe** (arXiv:2605.29035 structural analogy): on the 8,260-family band
  bank, test whether an exact third-moment (cubic) functional of the good-set
  indicator decides the 19 direct-L closures that the quadratic disc certificate
  (THM-731/732) misses — if yes, the band protocol becomes one uniform inequality
  (strengthen-then-deduce, the Frank–Ivanisvili shape: bulk by low moments, sharp
  constant one moment deeper).
- **Rédei deck-parity probe:** does the sheet deck carry a parity statement (odd
  number of free sheets under a suitable weighting) refining bare existence? Aim at
  equality/boundary structure only (guardrail C18).

## Verification and tooling

- THM-761 battery: 04-computation/lrc14_multi_exception_sheet_bound_opus_S299.py
  (+ .out): 50,964 exact counting instances, exact failure sets, end-to-end 13-speed
  closures including the q25-refutation family, r=7 wall realization.
- Terminal certificate in the library: `sheet_certificate(speeds, c)` in
  lrc14_certificates.py (self-test 15/15) — speed-only, O(r), exact.
- Exact falsifier and covering extension:
  `04-computation/lrc14_affine_slope_suspension_codex_S2.py` (+ stored `.out`).
- The natural next sessions are: (i) formulate and stress-test a peel-relative
  replacement for the false `r_P≤B(c*)` bridge; (ii) prove that its residual has a
  bounded or recursive normalized atlas; (iii) classify the `r≥7` deck residue.
