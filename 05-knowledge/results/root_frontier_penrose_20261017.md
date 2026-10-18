# A Penrose frontier for the ROOT transport: basin currents, and the termination climb near a fixed scale

2026-10-17, opus-2026-10-17-S42. Owner prompt: keep working on uniform ROOT energy; read two new preprints
(OpenAI, *The Quasi-Riemann Hypothesis*, Oct 5 2026; OpenAI, *Conformal flow and the Riemannian Penrose inequality
with minimizing frontiers*, Sep 27 2026); pattern-match their ideas liberally onto our work. This note records the
transfer that produced a concrete new result. The broader dictionary is in
[preprint_connections_solfeggio_digits_20261017.md](preprint_connections_solfeggio_digits_20261017.md).

**Typing.**
- PROVED (one line): the frontier inequality `R(N)^2 <= F_Y E_Y(N)`.
- FINITE-EXACT: the stopped energies for every `N < 2^23` (`Y = 30, 100, 1000, 10^4`), every `N < 2^21`
  (`Y = 3, 10, 10^5`), every `N < 2^26` (`Y = 30, 100, 1000, 10^4, 10^5`) and every `N < 2^28`
  (`Y = 100, 1000, 10^4`), and every `N < 2^30` (`Y = 100, 1000, 10^4`, int64-guarded). The `R(N)` invariant is
  checked at every `N`. The S42 audit reproduced every octave statistic to `2^28` with independent code.
- FINITE-EXACT: the scale profiles and the cross-term decomposition (exact integer arithmetic at the listed `N`).
- ANALOGY: the Penrose dictionary.
- OPEN (CONJECTURE): HYP-9290, the frontier energy cone.

**Prior art.** clock_duality_oct15's terminal mesh is the closest method
([collatz_terminal_mesh_20261015a.md](collatz_terminal_mesh_20261015a.md): T1, terminal extension never worsens the l1
bound; T2, Gram bound over unknown terminal labels). It already stops sources at finite terminals and bounds the current
without knowing their labels. What is new here:
- a fixed low frontier (first seed `<= Y`);
- the l2/Cauchy-Schwarz form of the bound;
- the basin form, i.e. Lemma F as Cauchy-Schwarz over first-entry basins (section 2);
- exact evidence that the stopped energy stays `O(N)` while the large climbs happen at small seeds;
- a partial scale collapse;
- the conditional power-saving route.

## 1. What the Penrose proof does, and the transfer

Bray's conformal flow, in the preprint's all-dimension version:
- **Area is preserved.** The flow `g_t = u_t^(4/d) g_0` keeps the frontier area `A` fixed.
- **Mass decreases.** It drives the ADM mass down by `m'(t) = -2 μ(t)`, where `μ = m - c >= 0` is mass minus capacity
  (the mass-capacity inequality, proved by **doubling the exterior with an odd reflection** `ψ = ±φ` and applying
  positive mass).
- **The end limit.** It normalises coordinates by `x = e^(-2t/d) y`. It never asks the singular frontier to converge.
  It replaces the frontier by its outer-minimizing hull, keeps the area, and finishes with an *enclosing-sphere*
  comparison at the largest limiting radius.

The ROOT transport of THM-4667 has the same skeleton, with the same factor 2:

| conformal flow | signed ROOT transport |
|---|---|
| preserved area `A` | preserved invariant `Σ_u c_d(u) f(u) = R(N)` |
| mass `m(t)` | energy `E_d(N) = Σ_u c_d(u)^2` |
| `m' = -2μ` | `E_(d+1) - E_d = 2 Cross_d` |
| mass-capacity `μ >= 0` | **Lemma F**: `E_Y >= R^2/F_Y`, where the capacity is the minimal energy carrying the invariant on `F_Y` seeds |
| monotone mass, `m' = -2μ <= 0` | dissipation `Cross_d <= 0` (ideal regime, `d <= 12`) |
| doubling by odd reflection `±φ` across a seam | sibling pairs `(s, 2s+1)` merged with opposite signs (statement 8) |
| normalised coordinate `x = e^(-2t/d) y` | band normalisation `M_d = 4(3/2)^d` of clock_duality's R8 |
| Schwarzschild equality `m = (A/ω)^(d/k)/2` | `E_∞ = R(N)^2` only when the unrooted energy vanishes |
| singular minimizing frontier (codimension 8) | ROOT, where coherent merging spoils monotonicity |
| outer-minimizing hull + enclosing spheres | **a frontier at seed height `Y` + Cauchy-Schwarz** (this note) |

**Pointwise Penrose inequality.** The pointwise analogue of the Penrose inequality, `R(N)^2 <= E_0(N)`, is **false**:
at `N = 7910`, `R^2 = 9.08 N` but `E_0 = 4.86 N`. Monotone mass fails because the termination merges have positive cross
terms; in the conformal flow, this is like negative scalar curvature behind the horizon.

**What the Penrose argument suggests instead.** Penrose needs positivity only *outside* an outer-minimizing frontier, and
handles the frontier by a comparison. Here that means stopping the transport at a frontier and comparing there.

## 2. The frontier inequality

For `Y >= 1`, stop each source `x` at the first seed `u <= Y` on its odd-transport path. Sources `<= Y` stop at once.
Write:
- `c_Y(u) = Σ_(x stops at u) σ(x -> u) c_0(x)` for the frontier coefficients;
- `E_Y(N) = Σ_u c_Y(u)^2` for the frontier energy;
- `F_Y = #{seeds <= Y} ≈ 3Y/8`.

The odd transport preserves the `f`-weighted sum at every step, and stopping only freezes coordinates. Hence:

**Lemma F (PROVED).** `R(N) = Σ_(u <= Y) f(u) c_Y(u)`, and therefore

    R(N)^2 <= F_Y · E_Y(N).

- `E_Y` is the energy of all actual hitting contributions. Nonhitting sources remain in a moving complement. Define f=0 off the ROOT basin; only initially ROOT-supported sources are needed for Lemma F. Their finitely many arrivals suffice for the invariant, without asserting that the moving complement disappears.
- Its cross terms include moving/moving merges above and at the frontier, and moving/frozen arrivals. Only transport after freezing is omitted.
- The final coherent merging near ROOT is replaced by the factor `F_Y`. This is the analogue of the enclosing-sphere
  comparison, and the Cauchy-Schwarz form of clock_duality's label-free terminal bound.
- Equality holds iff `c_Y(u) = (R/F_Y) f(u)` on every frontier seed, i.e. `E_Y = R^2/F_Y`: the maximally spread state.
  This is the analogue of the capacity, the minimal energy carrying the invariant on `F_Y` seeds.

**Basin form (PROVED; checked by `root_transport_basin_currents_check_20261018.py`, 1522 seeds, 0 violations).**
- Along any path the transport sign is `σ(x -> u) = f(x) f(u)`.
- So `c_Y(u) = f(u) R_u(N)`, where `R_u(N) = Σ_(m odd <= 4N, first seed <= Y of m is u) f(m) k_N(m)` is the ROOT current
  restricted to the **first-entry basin** of `u`.
- Hence `R(N) = Σ_u R_u(N)` and `E_Y(N) = Σ_u R_u(N)^2`.

Lemma F is Cauchy-Schwarz over the basins. HYP-9290 asks that the basin currents have square-root size in mean square.
For fixed `Y` this is equivalent, up to the factor `F_Y`, to every basin current being `O(sqrt N)`.

The termination climb is the basin currents aligning in sign. Every single basin current stays small
(`root_transport_basin_maxima_20261018.py`): `max_u R_u(N)^2/N <= 0.78` (`Y = 100`) and `<= 0.22` (`Y = 1000`) in every
octave `2^10-2^27`, with no growth. Meanwhile their sum reaches `R^2/N = 9.08`.

## 3. The data

The method is `root_transport_frontier_energy_20261017.py`. The stopped vector changes in only three source entries
from `N` to `N + 1`, so every `N` costs `O(1)`. `R(N) = H(4N) - H(N)` is verified at every `N`.

**Maxima of `E_Y(N)/N` by octave, every `N < 2^23`:**

| octave | `2^12` | `2^13` | `2^14` | `2^15` | `2^16` | `2^17` | `2^18` | `2^19` | `2^20` | `2^21` | `2^22` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `Y = 100` | 2.35 | 1.68 | 1.80 | 1.65 | 1.85 | 1.22 | 1.60 | 1.48 | 1.24 | 1.46 | 2.03 |
| `Y = 1000` | 2.08 | 1.92 | 1.75 | 1.41 | 1.42 | 1.44 | 1.38 | 1.45 | 1.50 | 1.52 | 1.31 |
| `R(N)^2/N` (unstopped) | 9.08 | 7.53 | 2.55 | 3.36 | 4.23 | 2.77 | 3.92 | 1.66 | 1.56 | 0.98 | 6.06 |

- In octave `2^12`, where `R^2/N` reaches 9.08, the `Y = 1000` frontier energy never exceeds `2.08 N`. In octave `2^22`
  (`R^2/N` up to 6.06) it stays below `1.31 N`.
- Lemma F with `Y = 1000` gives `|R(N)| <= 28 sqrt(N)` for every `N` in `[2^12, 2^23)`. This is lossy by about 9
  against the truth, but it is a square-root bound obtained without the termination regime.

**Out of sample: every `N < 2^30`** (`root_transport_frontier_energy_big_20261017.py`, memory-lean, and its
int64-guarded copy `..._huge_20261018.py` for `2^28-2^30`). Sixteen tests were pre-registered in HYP-9290 before the
scans, and all passed. Maxima below are rounded up.

| octave | `2^23` | `2^24` | `2^25` | `2^26` | `2^27` | `2^28` | `2^29` |
|---|---|---|---|---|---|---|---|
| `max E_100/N` | 2.09 | 1.59 | 2.63 | 1.39 | 1.80 | 1.37 | 1.20 |
| `max E_1000/N` | 1.34 | 1.24 | 1.45 | 1.30 | 1.25 | 1.34 | 1.23 |
| `max E_10000/N` | 1.41 | 1.42 | 1.39 | 1.27 | 1.27 | 1.31 | 1.24 |
| `max R(N)^2/N` | 3.08 | 2.64 | 4.37 | 7.49 | 7.38 | 3.88 | 1.65 |

So on the stated range `2^12 <= N < 2^30` (not the untested finite-prefix calibration):
- `E_1000(N) <= 2.08 N`;
- `E_100(N) <= 2.63 N`;
- `|R(N)| <= sqrt(375 · 2.08 N) < 28 sqrt(N)` by Lemma F (with `Y = 1000`, for `N >= 2^12`).

At six arg-max points in octaves `2^23-2^25` (P6), merges above `Y = 1000` remove `3.57-3.64 N`, and arrivals add at
most `0.246 N`.

**Certified closed frontiers** (clock_duality's R6; `root_transport_closed_frontier_energy_20261017.py`, P11-P12
pre-registered, passed).
- Setting: stop on first entry into the forward-closed sets `V_100` (74 seeds) and `V_1000` (627 seeds). There the
  stopped and unstopped transports agree outside `V`.
- Result: over `2^12 <= N < 2^26`, `E_(V_1000)/N <= 2.47` and `E_(V_100)/N <= 2.56`, octave profiles close to the
  height frontiers'.
- Observation (finite, `N < 2^26`): the frontier energy shows no dependence on the return kernel that the height
  frontier hides. The octave maxima of `E_(V_1000)` and `E_1000` differ by at most 18.8%.

**Scale profiles** (`root_transport_frontier_profile_20261017.py`). These show `E_Y/N` as the frontier height `Y` comes
down from `4N`.

| `Y` | 4N | `10^6` | `10^5` | `10^4` | 1000 | 300 | 100 | 30 | 10 | 3 | 1 (`= R^2/N`) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `N = 7910` | 4.86 | | | 3.87 | 1.92 | 1.88 | 1.78 | 2.40 | 4.19 | 9.29 | 9.08 |
| `N = 66463` | 4.86 | | 4.12 | 1.81 | 1.25 | 1.20 | 1.29 | 1.25 | 1.46 | 3.14 | 4.23 |
| `N = 2^20` | 4.86 | 3.28 | 1.71 | 1.47 | 1.37 | 1.33 | 1.22 | 1.14 | 0.96 | 1.06 | 1.53 |
| `N = 6904411` | 4.86 | 1.81 | 1.48 | 1.36 | 1.17 | 1.12 | 1.37 | 4.09 | 6.40 | 7.23 | 6.06 |

- The energy falls as the frontier comes down, to about 1.1-1.9 near `Y ≈ 100-1000`. At `N = 2^20` it falls further, to
  0.96 at `Y = 10`.
- In these four profiles every large climb happens below `Y ≈ 30-100`; it is the merging of a few dozen frontier
  currents into ROOT. The S42 audit (a8) added two out-of-sample profiles:
  - at the `R^2/N = 7.48` record (`N = 133,970,036`), the climb is again below `Y = 100`;
  - at the `E_100` record (`N = 48,251,574`), the profile is 1.29, 1.38, 1.42, 1.60, **2.62**, 1.31, 1.90, 3.09, 1.33
    for `Y = 10^4, 3000, 1000, 300, 100, 30, 10, 3, 1`. That is a transient bump of `+1.02 N` between `Y = 300` and
    100, which cancels below 30 while `R^2/N` stays 1.33.
  - So the climb is not always confined below 30.
- **Partial scale collapse.** For `N/Y` up to about 10, the profile depends essentially on `N/Y` alone. For example, at
  `N/Y ≈ 7`, `E_Y/N = 1.81` for both `N = 66463` and `N = 6904411`. The octave means agree to about 5%:
  - at `N/Y ≈ 6`, 1.91 (`Y = 1000`) against about 1.85 (`Y = 10^4`, interpolated);
  - at `N/Y ≈ 12`, 1.74 against 1.65.

  Beyond `N/Y ≈ 20` the curves separate into plateaus whose height depends on `Y`:
  - `Y = 10^4` descends slowly: octave means 1.61, 1.57, 1.51, 1.45, 1.44, 1.37 for `N/Y` = 20-630, then 1.34, 1.32,
    1.29, 1.20, 1.21, 1.26, 1.19 to `2^29`;
  - `Y = 1000` has octave means 1.185-1.345 for `2^15-2^22`, then 1.16, 1.04, 1.28, 1.21, 1.10, 1.18, 1.02 for
    `2^23-2^29`;
  - `Y = 100` has octave means 0.83-1.71 for `N/Y >= 20` (1.71 in octave `2^12`, 1.69 in `2^25`). It is noisy, with
    38 frontier seeds.

  So the scale-free part (the ideal regime seen in scale rather than depth) is visible only for `N/Y <= 10-20`. Beyond
  it the octave means of the stopped energy stay in about `1.0-1.35 N` at effective depths about 35, far beyond the
  banded ladder (`d <= 12`, `E_12 = 1.479`). There is no monotone trend to `2^29`.

**Where the climb comes from** (`root_transport_frontier_decomposition_20261017.py`). The stopped process is run
synchronously in depth, and every energy increment is split into three kinds:
- `MM_hi`, two moving currents merging above `Y`;
- `MM_lo`, two moving currents arriving at the same frozen seed in the same step;
- `MF`, a moving current arriving at a frozen seed that already carries current.

At each finite depth, `E_frozen + E_moving = E_0 + 2(MM_hi + MM_lo + MF)` exactly. In the completed finite runs below the moving vector vanishes, giving the displayed E_Y identity. Values below are `/N`, with the factor2 included.

| `N` | `Y` | `E_0` | merges above `Y` | same-step at frontier | arrivals | `E_Y` | `R^2` |
|---|---|---|---|---|---|---|---|
| 7910 | 30 | 4.856 | -3.237 | -0.115 | +0.896 | 2.401 | 9.080 |
| 7910 | 1000 | 4.856 | -2.739 | -0.366 | +0.164 | 1.916 | 9.080 |
| 66463 | 1000 | 4.857 | -3.267 | -0.103 | -0.240 | 1.248 | 4.226 |
| 1048576 | 100 | 4.857 | -3.593 | -0.003 | -0.037 | 1.224 | 1.533 |
| 6904411 | 30 | 4.857 | -3.653 | -0.031 | **+2.915** | 4.088 | 6.055 |
| 6904411 | 100 | 4.857 | -3.643 | -0.037 | +0.195 | 1.372 | 6.055 |
| 6904411 | 1000 | 4.857 | -3.547 | -0.039 | -0.101 | 1.170 | 6.055 |

- **Merges above the frontier are dissipative.** At every computed `(N, Y)` with `N >= 6Y` they remove 54-76% of the
  initial energy `E_0 ≈ (34/7) N`. Almost all of it happens in the first third of the stopped depth.
- **Arrivals.** Arrivals against frozen seeds are the only positive term. They are small at `Y >= 100` (at most `+0.24 N`
  in all 12 cases), and they *are* the climb at `Y = 30`. At `N = 6904411`, the arrivals at seeds below 30 add
  `+2.92 N`, two thirds of it in the middle third of the depth.

The termination regime of THM-4667 is therefore an arrival phenomenon: currents that have already come down meet the
current frozen at a few small seeds, with the same sign.

## 4. What this changes in the plan

1. **The termination regime is (in the data) a small-seed problem.** In the computed profiles the large coherent
   climbs happen among seeds below about 100, apart from transient bumps such as the one at `N = 48,251,574`. HYP-9288's third item (shell partial sums) is therefore paid by Lemma F at a constant cost
   `sqrt(F_(Y_0))`, provided the frontier energy above `Y_0` is controlled.
2. **What the frontier cone needs.**
   - By the decomposition, the frontier cone follows from two statements. The implication (i) + (ii) ⇒ (A) is one
     line, as HYP-9290 (B) states. Whether (i) and (ii) are also *necessary* is unknown.
     - (i) moving/moving merges, including simultaneous frontier arrivals, satisfy `MM_hi + MM_lo <= 0` uniformly in finite depth;
     - (ii) arrivals at frozen seeds, for a cutoff `Y>=Y_0`, stay `O(N)` uniformly in depth. Individual landing seeds can be smaller than Y_0.
   - Statement (i) is the stopped-process form of HYP-9288 (C), ideal dissipation. It is not the same statement.
     - The stopped process cuts every path at first entrance, so its coefficients above `Y` differ from the unstopped
       ones by the currents that dipped below `Y` and came back.
     - S41's old matching test omitted ideal-only seeds. The [whole-support repair](root_sparse_boundary_module_20261017i.md) proves O(d^2) dyadic error support and a growing-depth window, but does not identify stopped and unstopped moving vectors. R5 there retains the first-return kernel explicitly.
   - Statement (ii) is new. It asks that a current coming down to the frontier be incoherent with what is already
     frozen there. In the data this holds at `Y >= 100` and fails at `Y = 30`.
   - Neither statement involves the coherent final merging.
3. **Conditional power saving.** Take `Y = N^θ`. Lemma F gives `|R(N)| <= sqrt(F_Y E_Y) ≪ N^((1+θ)/2)`, provided
   `E_(N^θ)(N) ≪ N`.
   - HEURISTIC depth scale: about `(1-θ) log_(4/3) N`. This is not a deterministic stopping-time bound; one step can overshoot far below N^θ. Uniform stopped-flow and arrival estimates remain necessary.
   - With clock_duality's horizons (P3/O3, `l1`, depth `m/10`) one would need their `l2` analogue to depth
     `(1-θ)·2.4 log_2 N`. Comparing those heuristic horizons does not by itself remove the stopped-flow, arrival or coverage obligations.
4. **Relation to the octave cone.** HYP-9288 (B) bounds the unstopped traces. HYP-9290 bounds the stopped energy.
   - No implication is known in either direction.
   - HYP-9290 at a fixed `Y_0` gives the pointwise bound `|R(N)| <= sqrt(F_(Y_0) C N)`.
   - HEURISTIC: if the frontier coefficients fluctuate like independent Gaussian sums, the running maximum of `E_Y/N`
     still drifts by `O(sqrt(log log N / F_Y))`. That is why HYP-9290 also has a log log form.

## 5. HYP-9290 (filed)

**Frontier energy cone.** [HYP-9290](../hypotheses/HYP-9290-frontier-energy-cone-of-the-stopped-root-transport.md): there
are `Y_0` and `C` with `E_(Y_0)(N) <= C N` for all `N`. The initial Y_0=100,C=2.4 proposal was surpassed by the incoming E_100/N>2.62 record at N=48,251,574. The large-source test now uses C=3. Any all-N C must also be at least10, from N=1 and Y>=3. A weaker form allows
`C log log N`. Consequence (Lemma F): `|R(N)| <= sqrt(F_(Y_0) C N)`.

Sixteen tests were pre-registered before the corresponding scans, and all sixteen passed (HYP-9290). The audit
verified the order of commits for P1-P14.
- P1-P6: `2^23 <= N < 2^26`.
- P7-P10: `2^26 <= N < 2^28`.
- P11-P12: the closed frontiers.
- P13-P14: plateau statistics, read descriptively after the audit.
- P15-P16: `2^28 <= N < 2^30`.

The independent audit is `04-computation/experiments/s42_audit_20261017/REPORT.md`.

## Reproduction

```
python 04-computation/experiments/root_transport_frontier_energy_20261017.py 22 30 100 1000 10000
python 04-computation/experiments/root_transport_frontier_energy_20261017.py 20 3 10 100 1000 10000 100000
python 04-computation/experiments/root_transport_frontier_energy_big_20261017.py 25 30 100 1000 10000 100000
python 04-computation/experiments/root_transport_frontier_energy_big_20261017.py 27 100 1000 10000
python 04-computation/experiments/root_transport_frontier_decomposition_20261017.py 9666088,32532811,47930131,8429844,32184394,37675352 1000
python 04-computation/experiments/root_transport_frontier_profile_20261017.py 7910 66463 1048576 6904411
python 04-computation/experiments/root_transport_frontier_decomposition_20261017.py 7910,66463,1048576,6904411 30,100,1000,10000
python 04-computation/experiments/root_transport_frontier_energy_huge_20261018.py 29 100 1000 10000
python 04-computation/experiments/root_transport_closed_frontier_energy_20261017.py 25
python 04-computation/experiments/root_transport_frontier_statistics_20261017.py 27
python 04-computation/experiments/root_transport_basin_currents_check_20261018.py
python 04-computation/experiments/root_transport_basin_maxima_20261018.py 27 100
python 04-computation/experiments/root_transport_basin_maxima_20261018.py 27 1000
```
