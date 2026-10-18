# A Penrose frontier for the ROOT transport: the termination climb lives below a fixed scale

2026-10-17, opus-2026-10-17-S42. Owner prompt: keep working on uniform ROOT energy; read two new preprints
(OpenAI, *The Quasi-Riemann Hypothesis*, Oct 5 2026; OpenAI, *Conformal flow and the Riemannian Penrose inequality
with minimizing frontiers*, Sep 27 2026); pattern-match their ideas liberally onto our work. This note records the
transfer that produced a concrete new result. The broader dictionary is in
[preprint_connections_solfeggio_digits_20261017.md](preprint_connections_solfeggio_digits_20261017.md).

**Typing.**
- PROVED (one line): the frontier inequality `R(N)^2 <= F_Y E_Y(N)`.
- FINITE-EXACT: the stopped energies for every `N < 2^23` (`Y = 30, 100, 1000, 10^4`), every `N < 2^21`
  (`Y = 3, 10, 10^5`), every `N < 2^26` (`Y = 30, 100, 1000, 10^4, 10^5`) and every `N < 2^28`
  (`Y = 100, 1000, 10^4`). The `R(N)` invariant is checked at every `N`.
- FINITE-EXACT: the scale profiles and the cross-term decomposition (exact integer arithmetic at the listed `N`).
- ANALOGY: the Penrose dictionary.
- OPEN (CONJECTURE): HYP-9290, the frontier energy cone.

**Prior art.** clock_duality_oct15's terminal mesh is the closest method
([collatz_terminal_mesh_20261015a.md](collatz_terminal_mesh_20261015a.md): T1, terminal extension never worsens the l1
bound; T2, Gram bound over unknown terminal labels). It already stops sources at finite terminals and bounds the current
without knowing their labels. What is new here:
- a fixed low frontier (first seed `<= Y`);
- the l2/Cauchy-Schwarz form of the bound;
- exact evidence that the stopped energy stays `O(N)` while the climb is confined below a fixed scale;
- the scale collapse;
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
| mass-capacity `μ >= 0` | dissipation `Cross_d <= 0` (ideal regime, `d <= 12`) |
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

## 3. The data: the climb lives below a fixed scale

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

**Out of sample: every `N < 2^28`** (`root_transport_frontier_energy_big_20261017.py`, memory-lean; ten tests
pre-registered in HYP-9290 before the scans, all passed).

| octave | `2^23` | `2^24` | `2^25` | `2^26` | `2^27` |
|---|---|---|---|---|---|
| `max E_100/N` | 2.08 | 1.59 | 2.62 | 1.38 | 1.79 |
| `max E_1000/N` | 1.34 | 1.24 | 1.45 | 1.30 | 1.25 |
| `max E_10000/N` | 1.40 | 1.42 | 1.38 | 1.26 | 1.27 |
| `max R(N)^2/N` | 3.08 | 2.63 | 4.36 | 7.48 | 7.37 |

So on the stated range `2^12 <= N < 2^28` (not the untested finite-prefix calibration):
- `E_1000(N) <= 2.08 N`;
- `E_100(N) <= 2.63 N`;
- `|R(N)| <= sqrt(375 · 2.08 N) < 28 sqrt(N)` by Lemma F (with `Y = 1000`, for `N >= 2^12`).

At six arg-max points in octaves `2^23-2^25` (P6), merges above `Y = 1000` remove `3.58-3.63 N`, and arrivals add at
most `0.245 N`.

**Scale profiles** (`root_transport_frontier_profile_20261017.py`). These show `E_Y/N` as the frontier height `Y` comes
down from `4N`.

| `Y` | 4N | `10^6` | `10^5` | `10^4` | 1000 | 300 | 100 | 30 | 10 | 3 | 1 (`= R^2/N`) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `N = 7910` | 4.86 | | | 3.87 | 1.92 | 1.88 | 1.78 | 2.40 | 4.19 | 9.29 | 9.08 |
| `N = 66463` | 4.86 | | 4.12 | 1.81 | 1.25 | 1.20 | 1.29 | 1.25 | 1.46 | 3.14 | 4.23 |
| `N = 2^20` | 4.86 | 3.28 | 1.71 | 1.47 | 1.37 | 1.33 | 1.22 | 1.14 | 0.96 | 1.06 | 1.53 |
| `N = 6904411` | 4.86 | 1.81 | 1.48 | 1.36 | 1.17 | 1.12 | 1.37 | 4.09 | 6.40 | 7.23 | 6.06 |

- The energy falls as the frontier comes down, to a minimum of about 1.1-1.9 near `Y ≈ 100-1000`.
- Every large climb happens **below `Y ≈ 30-100`**, at every `N` from `2^13` to `7·10^6`. It is the merging of a few
  dozen frontier currents into ROOT.
- **Partial scale collapse.** For `N/Y` up to about 10, the profile depends essentially on `N/Y` alone. For example, at
  `N/Y ≈ 7`, `E_Y/N = 1.81` for both `N = 66463` and `N = 6904411`. The octave means agree to about 5%:
  - at `N/Y ≈ 6`, 1.91 (`Y = 1000`) against about 1.85 (`Y = 10^4`, interpolated);
  - at `N/Y ≈ 12`, 1.74 against 1.65.

  Beyond `N/Y ≈ 20` the curves separate into plateaus whose height depends on `Y`:
  - `Y = 10^4` is still descending slowly (octave means 1.61, 1.57, 1.51, 1.45, 1.44, 1.37 for `N/Y` = 20-630);
  - `Y = 1000` is flat at 1.20-1.34 from `N/Y ≈ 50` to 6300;
  - `Y = 100` is flat at about 0.8-1.35, noisy with 38 frontier seeds.

  So the scale-free part (the ideal regime seen in scale rather than depth) is visible only for `N/Y <= 10-20`. The
  plateau is a property of the frontier itself. Either way the stopped energy stays near `1.2-1.5 N` at effective
  depths about 35, far beyond the banded ladder (`d <= 12`, `E_12 = 1.479`).

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

- **Merges above the frontier are dissipative.** At every computed `(N, Y)` with `N >= 6Y` they remove 54-75% of the
  initial energy `E_0 ≈ (34/7) N`. Almost all of it happens in the first third of the stopped depth.
- **Arrivals.** Arrivals against frozen seeds are the only positive term. They are small at `Y >= 100` (at most `+0.24 N`
  in all 12 cases), and they *are* the climb at `Y = 30`. At `N = 6904411`, the arrivals at seeds below 30 add
  `+2.92 N`, two thirds of it in the middle third of the depth.

The termination regime of THM-4667 is therefore an arrival phenomenon: currents that have already come down meet the
current frozen at a few small seeds, with the same sign.

## 4. What this changes in the plan

1. **The termination regime is a fixed-scale problem.** In the data, all coherent climbing happens among seeds below
   about 100. HYP-9288's third item (shell partial sums) is therefore paid by Lemma F at a constant cost
   `sqrt(F_(Y_0))`, provided the frontier energy above `Y_0` is controlled.
2. **What the frontier cone needs (HEURISTIC reduction, not an implication).**
   - By the decomposition, the frontier cone follows from two statements:
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
   - Neither implies the other.
   - HYP-9290 at a fixed `Y_0` gives the pointwise bound `|R(N)| <= sqrt(F_(Y_0) C N)`.
   - HEURISTIC: if the frontier coefficients fluctuate like independent Gaussian sums, the running maximum of `E_Y/N`
     still drifts by `O(sqrt(log log N / F_Y))`. That is why HYP-9290 also has a log log form.

## 5. HYP-9290 (filed)

**Frontier energy cone.** [HYP-9290](../hypotheses/HYP-9290-frontier-energy-cone-of-the-stopped-root-transport.md): there
are `Y_0` and `C` with `E_(Y_0)(N) <= C N` for all `N`. The initial Y_0=100,C=2.4 proposal was surpassed by the incoming E_100/N>2.62 record at N=48,251,574. The large-source test now uses C=3. Any all-N C must also be at least10, from N=1 and Y>=3. A weaker form allows
`C log log N`. Consequence (Lemma F): `|R(N)| <= sqrt(F_(Y_0) C N)`.

Ten tests were pre-registered before the corresponding scans: P1-P6 for `2^23 <= N < 2^26`, P7-P10 for
`2^26 <= N < 2^28`. All ten passed (HYP-9290).

## Reproduction

```
python 04-computation/experiments/root_transport_frontier_energy_20261017.py 22 30 100 1000 10000
python 04-computation/experiments/root_transport_frontier_energy_20261017.py 20 3 10 100 1000 10000 100000
python 04-computation/experiments/root_transport_frontier_energy_big_20261017.py 25 30 100 1000 10000 100000
python 04-computation/experiments/root_transport_frontier_energy_big_20261017.py 27 100 1000 10000
python 04-computation/experiments/root_transport_frontier_decomposition_20261017.py 9666088,32532811,47930131,8429844,32184394,37675352 1000
python 04-computation/experiments/root_transport_frontier_profile_20261017.py 7910 66463 1048576 6904411
python 04-computation/experiments/root_transport_frontier_decomposition_20261017.py 7910,66463,1048576,6904411 30,100,1000,10000
```
