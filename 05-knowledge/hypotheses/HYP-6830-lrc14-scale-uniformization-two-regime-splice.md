---
id: HYP-6830
title: Scale uniformization of Claim B — proved sheet regime, refuted raw-fragmentation complementarity, and a peel-relative splice
status: PARTIALLY PROVED / REFUTED AS STATED — THM-761 proves its exact sheet-budget regime; raw fragmentation is not controlled by divisor-packet scale; a peel-relative splice and the remaining r>=7 residue stay open
source: opus-2026-07-14-S299
depends_on:
  - THM-755   # capped envelope, v* = r_P / (pi |G'_P|)
  - THM-760   # r=1 coprime sheet dodge
  - THM-761   # multi-exception sheet covering bound (this session)
  - THM-771   # exact seven-owner defect and corrected reduced-winding event pierce
  - THM-779   # exact scale-free f=4 transverse-fragmentation falsifier
  - THM-780   # invariant-level good-set-state transverse-tooth cap
  - THM-777   # rho bridge and exact bounded-height safe-mass floor
  - HYP-6780  # v*(cP) = c v*(P): the scale covariance that killed raw-height bands
related: [THM-756, THM-757, THM-758, HYP-6785, HYP-6815, HYP-6820, HYP-6835, MISTAKE-146]
---

# HYP-6830 — scale uniformization of the ≥4-far covering endgame

## The proposed splice, with its failed bridge isolated

For every covering 13-speed family V and every scale c ≥ 2, write V = cP ⊔ W
(canonical: P = {v/c : c|v}, r = |W|). Define c*(V) = the largest c with |P| ≥ 7
(equivalently r ≤ 6), or 1 if none. Then:

1. **Sheet-budget regime (PROVED, THM-761):** for a chosen decomposition
   `V=cP union W`, the exact per-`(r,c)` criterion or the mixed-gcd budget
   `sum_a g_a(floor(c/(7g_a))+1)<=c-1` gives a free witness sheet and
   `M(V)>=1/14`.  In particular `c>=43` is uniform when all exceptions are
   coprime to `c`.  The numerical condition `c*(V)>=43` alone is not a theorem
   hypothesis.

   For example `V={100,200,...,900,50,150,1,3}` has `c*(V)=100`, but at that
   canonical scale the exception gcds are `(50,50,1,1)` and the THM-761 budget
   is `130>99`.  This does not say that `V` violates LRC or lacks another
   routing scale; it shows why the exact gcd sidecar cannot be omitted.
2. **Complementary regime (OPEN — coordinate incomplete):** if `c*(V)<=42`, then
   every high-support dilated packet whose scale is counted by `c*` inflates its
   capped-envelope band edge by at most that bounded factor
   (`v*(cP)=c*v*(P)<=42*v*(P)`, HYP-6780 used positively). Low-support packets can
   still carry arbitrarily large common factors. The local covariance remains exact,
   but `c*(V)<=42` does not by itself bound fragmentation, nor has it produced a
   globally bounded normalized domain. TO PROVE, with a
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

## Exact infinite-family falsifier (THM-779)

For every prime `N>110`, put

```text
P_N = {1,2,...,9,15,110,N},
V_N = P_N union {1092 N},       1092 = lcm(12,13,14).
```

This family has all the scope properties needed to refute the proposed bridge.

1. **No high-support scale.** Among `{1,...,9,15,110}`, at most five entries
   share any nontrivial divisor (the four small even entries together with `110`
   attain the maximum). The prime `N` shares no divisor with them, and adjoining
   `1092N` can add at most one member to any old divisor packet. Hence the largest
   divisor packet has size five in `P_N`
   and six in `V_N`; in particular `c*(P_N)=c*(V_N)=1`.
2. **Unbounded exact fragmentation.** The base good set for
   `{1,...,9,15,110}` contains

   ```text
   J = [1/14, 111/1540],           |J| = 1/1540.
   ```

   Inside `J`, the `N`-runner removes the disjoint open danger teeth centered at
   `k/N` with radius `1/(14N)`. The admissible center interval has length
   `N/1540-1/7`, so it contains at least `N/1540-8/7` integers. Every resulting
   full tooth separates two components of `G'_{P_N}`. Thus `r_{P_N}` is
   unbounded although `c*(P_N)=1`. Exact interval arithmetic gives component counts
   `66,104,174,310` for `N=211,503,1009,2003`.
3. **It occurs inside the covering endgame.** Every `V_N` is primitive because it
   contains `1`. Speeds `2,...,9` carry their own moduli, `110` carries `10` and
   `11`, and `1092N` carries `12,13,14`; hence `V_N` is covering. Moreover it is
   literally in the first open far-count chart: its nine small speeds are
   `{1,...,9}` and its four far coordinates are `15,110,N,1092N`. This is not an
   LRC counterexample: THM-779 proves the top peel fires THM-755's capped-envelope
   test for every prime `N>110`. It is a counterexample specifically to using `c*` as a
   sufficient fragmentation coordinate, now internal to the four-dimensional
   object rather than imported from another stratum.

The intended `f>=4` lane does not rescue the raw claim.  Put
`B={1,...,8,21,22,23}`, `Q_N=B union {N}` for odd `N>23`, and adjoin a larger
multiple `W` of `360360`.  Then `Q_N union {W}` is primitive, covering, has
five speeds above `14`, and has `c*=1`.  Since
`|t-1/12|<1/1932` is uniformly safe for `B`, the same comb argument gives
`r_{Q_N}>=N/966-13/7`.

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

## A terminal transverse face (THM-780)

The falsifier mechanism itself admits a general theorem. For a fixed base `B`,
put `mu=|G'_B|` and let `r_B` be the number of components of `G'_B`. After
adjoining frequency `N`,

```text
|G'_{B union {N}}| >= 6mu/7-2r_B/(7N),
r(B union {N}) <= N+r_B.
```

Thus the old good-set state `(mu,r_B)`, rather than every endpoint, is sufficient
for this transition. A single marked safe interval of length `L` gives the more
portable corollary

```text
|G'_{B union {N}}| >= 6L/7-2/(7N),
r(B union {N}) <= N+sum(B).
```

Consequently, for every fixed base with `mu>0`,
`liminf_(N->infinity)|G'_{B union {N}}|>=6mu/7>0`. A single transverse tooth
frequency cannot cause safe-mass collapse; any global decay sequence must vary
the base or use correlated iterated/multiscale insertions.

A proportional peel `aN` closes beyond an explicit rational threshold
whenever `(333/106)*a*6L/7>1`. This proves that transverse wall proliferation
alone is a terminal face: it can make raw components unbounded while keeping
the theorem-facing load `r/(aN|G'|)` bounded. For THM-779, `L=1/1540`,
`a=1092`, and the exact crossing is `11734415/9278<1265`.

This does not prove the global splice. It sharpens its negation space. A truly
unresolved sequence must couple increasing wall frequency to safe-mass collapse
relative to its component state, to a subcritical peel rate, or to an owner
alignment that evades the proposed descent. Those three channels must remain
distinguishable.

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
- THM-779 proves that the entire chart-native prime family has an empty top band;
  the other high-frequency falsifier is compatible with the same candidate on its
  audited range. Neither fact proves the ratio-coordinate conjecture globally;
- THM-777 proves the bridge `rho <= 12/(pi |G'_P|)`, the unconditional tail
  `|G'_P|>=1/(91 max(P))`, and by exact census the sharp candidate floor
  `7/858` for every primitive 12-core with `max(P)<=18`. Its explicitly listed
  scale/tooth samples and five seeded searches find no lower value. The
  global floor is explicitly **CONJECTURAL**. If proved, it would give
  `rho<469` uniformly and bound the normalized regime-2 band; the census and
  adversarial battery alone do not do so.

The three scale statements now form a precise logical staircase. THM-779 shows
that component topology is not compact modulo `c*`. THM-780 shows that a fixed
positive safe-mass state survives transverse refinement strongly enough for a
proportional peel. THM-777 shows that a global positive safe-mass floor would
compactify the remaining ratio coordinate. What is missing is exactly either
that shape-space floor or a classification/descent of every sequence on which
safe mass tends to zero.

## Assumption challenge and tournament scope

The falsifier was found by rejecting the default runner/divisor-packet
vertices. Alternate vertices considered were runners, maximal divisor packets,
base safe components, individual tooth intervals, endpoint events, residues,
Fourier modes, and peel obligations. The operative incidence is bipartite:
an `N`-tooth cuts a named base component. A binary orientation of runner pairs
does not preserve how many cuts land in the same component, their widths, or
the cap load, so there is no clean proof tournament on runner vertices here.

The companion HYP-6815 audit instead uses candidate representations as
tournament vertices. Its pairwise observable is the fraction of
proof-critical row pairs separated, its two gauges prioritize predicate
retention or compression, and the declared tie Hamiltonian path is printed in
the stored output. That tournament is diagnostic: it confirms that scale-only
and residue-only representations lose metric topology. The exact carrier for
this claim is the component/tooth incidence with owner and peel sidecars, not
the orientation itself.

## The residual inside regime 1 (named by THM-761; unramified r=7 high reduced winding closed by THM-771)

- **r = 7, all effective grids unramified: CLOSED above a reduced-shape bound
  (THM-771).** Put `g_a=gcd(w_a,c)`, `C_a=c/g_a`, and `u_a=w_a/g_a`. If `7|C_a`
  for every owner, its count is `c/7` off its endpoint coset and `c/7-g_a` on it.
  The endpoint mesh is `1/u_a`; the six-speed core supplies a closed safe interval
  of length at least `1/(7 max(P))`. Thus `max u_a>=7 max(P)` forces a free sheet.
  The realized c=7 wall is pierced at `t0=1/8`. Remaining inside r>=7: r>=8,
  ramified effective grids, and bounded reduced winding.
- The promoted **KCL absorption law is withdrawn** (MISTAKE-146). Strict endpoint
  equality is safe for both exiting and entering owners, so it tears rather than
  transfers a bad sheet. The theorem-facing replacement is the exact incidence
  defect `F=Q+Omega-sigma`; the primitive c=21 row has `Q=0`, `Omega=sigma=12`,
  and no free sheet on its static chamber. Mirror arithmetic remains diagnostic only.
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
