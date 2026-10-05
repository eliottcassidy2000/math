# Information versus dimension for the Syracuse law: pressure bounds and the integer-source gap (2026-10-05)

**Current scope correction, 2026-10-05:** Proposition1 below now distinguishes
an exact valuation cylinder from a coarse one, and coefficient contraction
from actual descent. The spectrum upper bound is retained; equality and
the diffuse-regime lower bound remain conditional on unproved asymptotic
inputs. The script uses float64 and a finite valuation cutoff, so its data
are **VERIFIED finite numerical evidence**, not exact rational computations. Computable
points are not random for a computable atomless measure, but that is not a
no-go theorem for all measure arguments. The
[atomic-prefix follow-up](collatz_effective_prefix_mass_20261005.md) constructs
a computable measure with a positive atom at every target integer and gives
an exact mass-one coverage criterion. It does not prove the needed decay.
This targeted correction is not a full audit of the inherited spectrum,
local-dimension, or dimension-of-E_inf claims.

**Session:** opus, `pascal-lyapunov-20261005` (fourth task), 2026-10-05.
Owner's directive: read [arXiv:2408.05121](https://arxiv.org/abs/2408.05121) (Jan Reimann, *Information vs
Dimension -- an Algorithmic Perspective*) -- universal prefix-free machines, the Levin--Schnorr/Gács theorem for
computable measures, the effective/pointwise mass distribution principle, multifractal measures -- and say how it
meshes with the repo's Vitali-set thread and with the structure of Collatz the repo is trying to capture.
**Pattern:** the owner's arXiv analogy-bridge pattern: structural analogy (shape of the hard direction), reconciled
with the frontier, not subject overlap.
**Inherits (cited):** [`collatz_directions_20260926.md`](collatz_directions_20260926.md) Proposition 1 (`E_inf`
closed, box dimension `h*`); the [five-mirrors note](collatz_five_mirrors_20260929.md) (`3^(h*-1)` with `h* =
h(log_3 2)` the binary entropy, PROVED standard); [THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md)/[THM-4499](../../01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md)
(thin divergence `N(X) <= K X^(h*) (log X)^a`); the [Fourier experiments](collatz_fourier_experiments_20261004.md)
E2 (`E[rho_n^q]` bounded iff `q < 2`, the `t^-2` tail, the maximal atom at `-1`); the
[Cauchy criticality note](procgen_cauchy_20260925_cauchy_schwarz_criticality.md) (the `L^inf` dimension `log_3
2`; "the divergence of the `L^2` norm is a bulk multifractal effect of the `s = 2` tilt, odd frequency `3/4`");
[THM-4263](../../01-canon/theorems/THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport.md)
(density transport iff uniformly integrable fibre weights); [THM-4517](../../01-canon/theorems/THM-4517-no-root-uniform-positive-density-of-collatz-basins.md)
(no root-uniform density; the three `3x-1` basins); the LRC "Vitali wall" ([THM-406](../../01-canon/theorems/THM-406-covering-depth-master-object-factorial-moments-and-spectral-identity.md)
M2 and the [anti-Poisson atlas](../../07-reflections/anti-poisson-coimage-atlas-s604.md)); the
[Pascal-tower note](collatz_pascal_tower_lyapunov_20261005.md) (pressure, Jensen gap, the `0.1%` of frequencies).
**Paper, read (CITED):** Theorem 4.11 (Gács; `x` is `mu`-random iff `K(x|n) >= -log mu[[x|n]] - c`), Theorem 4.16
(effective mass distribution principle), Theorem 4.1 (Eggleston: `dim_H D(p) = H(p)/log N`), Corollary 4.13
(`dim_H(x) = liminf K(x|n)/n`), Theorem 4.18 (Shannon--McMillan--Breiman for random points), Proposition 5.8
(`mu`-random `x` has `dim_H(x) = delta_mu(x)`, the local dimension), Theorem 5.12 (algorithmic decomposition of the
multifractal spectrum; the universal semimeasure `M` as the perfect multifractal `f_M(alpha) = alpha`).

**Status: CITED (the paper; the repo's `h*`, `E_inf`, thin-divergence and `L^p` results); PROVED (elementary:
Proposition 1, the moment bound of Theorem 3, the identity of Proposition 4, the dimension-`0` remark of section 5);
VERIFIED (finite numerical) (`collatz_syracuse_multifractal_20261005.py`: the Syracuse law to level `16`, `43` million residues,
partition functions at eleven `q`; local dimensions at six integer points); CONJECTURED (equality in the
spectrum bound for `q > 2`); the bridge itself is a typed analogy. Collatz OPEN; nothing here is a Collatz step.**

---

## 0. What is new, in one screen

1. **The dictionary (section 1).** Reimann's identity "entropy = complexity = dimension" is already realized in the
   repo at one point: the coefficient-survivor set `E_inf` has the cited box-counting exponent
   `h* = H(1/log_2 3) = 0.94996`. This matches Eggleston's frequency-set entropy
   at the critical shortcut odd-bit frequency `log_3 2`; identification of
   Hausdorff dimension for the all-prefix set needs a separate argument; the thin-divergence bound `N(X) <= K X^(h*+o(1))` is its counting form (the
   mass distribution principle read on the integers below `X`). The paper's pointwise principle "a point that
   supports measure decaying like `2^(-sn)` has dimension at least `s`" has an exact Collatz counterpart:
   **coefficient descent at odd step `k` is the inequality `2^(-S_k) < 3^(-k)`** between the Haar mass of the
   point's 2-adic parity cylinder and the Haar mass of its 3-adic image cylinder (Proposition 1): descent is a
   comparison of two adic masses along the orbit, the binary clock against the ternary clock.
2. **The Syracuse law has a pressure upper bound and a proposed spectrum (section 3).** With `Z_n(q) = sum_c mu_n(c)^q`
   over `Z/3^n`, superadditivity gives **`tau(q) <= min(q - 1, log_3(2^q - 1))` for `q >= 1`** (PROVED); the two
   bounds cross exactly at `q = 2` because `sum_a 4^-a = 1/3 = 3^-1` (the Parseval = incoherent coincidence of the
   Fourier note in multifractal form), so the proposed `D_q = 1` for `1 < q <= 2` is conditional on the asymptotic moment bounds, and `D_q <= log_3(2^q - 1)/(q - 1)` beyond, with
   `D_inf <= log_3 2` (equality is proposed). The Legendre point at `q = 2` is the valuation-one frequency `p = 3/4`,
   through the exact identity `2 - H(3/4) = (3/4) log_2 3` (Proposition 4): the `s = 2` tilt of the Cauchy note
   and the energy-weighted step law `3 4^-a` of the ridge analysis are this point. Numerical partition functions to level
   `16` sit below the bound with `1/n` corrections (`tau_16(2.5) = 1.367` against `1.400`); equality is CONJECTURED.
3. **Hot spots are the negative cycles (section 4).** The local dimension of the Syracuse law at the 3-adic shadow of
   an integer is bounded by its parity cylinder: `log_3 2 times (halvings per odd step)`, i.e. `0.631, 0.946, 0.992`
   for `-1, -5, -17` and `1.26` for the trunk; measured at level `14`: `0.61, 0.87, 0.92` for the three negative
   cycles and `1.05, 1.05, 1.11` for `1, 5, 7`. The "earthquake hot spots" of the pasted paragraph are the
   `3x-1` cycles seen from the `3x+1` law; the displayed positive points have finite-level estimates near one; this does not establish their limiting density.
4. **The Vitali wall in algorithmic form (section 5).** Every integer is a computable point of `Z_2`: its parity
   vector has `K(x|m) = O(log m)`, effective dimension `0`, and it is `mu`-random for no computable `mu` without
   atoms. The paper's hard direction -- from almost everywhere to every random point, via complexity -- stops
   exactly at the points Collatz is about. This is the repo's LRC "Vitali wall" (moment and measure methods see
   only initial moments; the residual lives beyond every moment order) in its sharpest form: the whole
   multifractal spectrum, every `L^q` bound (P7), THM-4263's transport and the typical rates of the Pascal tower
   are statements about Haar-random points, and `Z^+ cap E_inf = empty` is a statement about computable ones.
   This is a limitation of an automatic typical-point upgrade, not of all measure arguments. The atomic-source reroute in the current correction retains every integer.

---

## 1. The dictionary

| Reimann (CITED) | repo object | map / preserved | lost / sidecar | status |
|---|---|---|---|---|
| Eggleston: `dim_H {x : digit frequencies -> p} = H(p)` | `E_inf` has box dimension `h* = H(log_3 2)` (directions note, Prop. 1); `3^(h*-1)` is the no-descent rate (five-mirrors) | parity vector of `x in Z_2`; the constraint "ones-frequency `>= log_3 2` at every prefix" is the no-descent condition | `E_inf` is defined by all prefixes, Eggleston's set by the limit; Hausdorff = box for `E_inf` is expected from the mass distribution principle with the uniform measure on its cylinders, not checked here | CITED + PROVED (repo) |
| mass distribution principle, counting form | THM-4476/4499: integers below `X` in the thin-divergence set number `<= X^(h*)` up to logs | cover by 2-adic cylinders of depth `log_2 X` | the polynomial prefactor (`(log X)^a`) | CITED |
| Gács/Levin--Schnorr: `mu`-random iff `K(x|n) >= -log mu[[x|n]] - c`; pointwise principle `-log mu[[x|n]] >= sn` gives `dim(x) >= s` | coefficient contraction at step `k` iff `2^(-S_k) < 3^(-k)`: the parity cylinder is lighter than the image cylinder (Proposition 1) | the two Haar measures along one orbit; `-log_2` of the 2-adic cylinder is the binary clock `S_k`, `-log_3` of the 3-adic cylinder is the ternary clock `k` | the paper compares one measure with the uniform one; here two measures on two spaces, linked by the map | PROVED (elementary) |
| Shannon--McMillan--Breiman for random points, Prop. 5.8: `dim(x) = delta_mu(x)` | Haar-a.e. 2-adic `x` has ones-frequency `1/2 < log_3 2` (descends) and the Syracuse law has local dimension `1` at its image (the density is finite a.e. under P7) | typical points | nothing about computable points (section 5) | CITED + SUPPORTED |
| multifractal spectrum `f_mu(alpha)`, generalized dimensions | `tau(q) <= min(q-1, log_3(2^q-1))`, `D_inf = log_3 2`, hot spots at `-1, -5, -17` (sections 3-4) | partition functions of `mu_n` | equality for `q > 2` | PROVED bound + VERIFIED (finite numerical) + CONJECTURED equality |
| universal semimeasure `M` as the perfect multifractal `f_M(alpha) = alpha` | the Syracuse spectrum sits below the line: `f(alpha) = alpha H(1/(alpha log_2 3)) <= alpha` on its concentrated branch | the deficit `1 - H(p)` of the tilted valuation law | no analogue of universality | ANALOGY |
| Theorem 4.16: high effective dimension gives a dispersing measure through the point | no analogue for integers (dimension `0`) | -- | -- | the wall |

## 2. Proposition 1: coefficient contraction is a mass comparison (PROVED, corrected)

For odd `x` with accelerated orbit `x_0 = x, x_i = U(x_(i-1))`, valuations `a_i` and `S_k = a_1 + ... + a_k`, the
exact first `k` odd valuations determine one class modulo `2^(S_k+1)`.
Its Haar mass **conditional on an odd source** is `2^(-S_k)`; unconditional
binary Haar mass is `2^(-S_k-1)`. Modulo `2^S_k` one has only the coarse
condition, with the final valuation at least the prescribed one. On the exact cylinder,
`x_k = (3^k x + B)/2^(S_k)` with `B >= 0`. The image `x_k mod 3^k` is a 3-adic cylinder of Haar mass `3^(-k)`.

**Proposition 1 (repaired).** The affine **coefficient** contracts iff
`3^k < 2^(S_k)`, equivalently `2^(-S_k) < 3^(-k)`. These are respectively
the conditional binary source mass and the full ternary target-cylinder mass.
When the coefficient contracts, actual descent additionally requires
`x > B/(2^(S_k)-3^k)`. The coefficient-survivor set E_inf imposes failure
of coefficient contraction at every prefix; it must not be identified with
actual non-descent without retaining exceptional representatives.
*Proof.* Subtract x in the affine formula. The carry B cannot be dropped.
The minimal hostile is the valuation word(2) at x=1: its coefficient is3/4,
but U(1)=1. The exact word cylinder is1 mod8; the coarse cylinder1 mod4
also contains5, whose actual valuation is4. THM-4512 bounds the exceptional
representatives; it does not erase them. □

**Reading against Theorem 4.11/4.16.** Reimann's principle bounds the dimension of a point by how fast one
measure decays along it; Proposition 1 compares Collatz coefficients through of two measures (the Haar measures of
the two adic completions) along the orbit of one point, with the ratio `3^k 2^(-S_k)` as the running "dimension
ratio". For a Haar-random `x` the ratio tends to `0` (`S_k ~ 2k`); on `E_inf` it stays `>= 1`; the set where it is
critical (`S_k ~ k log_2 3`) has dimension `h* = H(log_3 2)` by Eggleston, which is the repo's constant. Nothing in
this is new to the repo except the typing; it is the exact place where the paper's idea lands.

## 3. The multifractal spectrum of the Syracuse law (PROVED bound, VERIFIED (finite numerical), CONJECTURED equality)

Let `mu_n` be the law of the `n`-th accelerated iterate modulo `3^n` for a Haar-random odd 2-adic start (the repo's
Syracuse law; density `rho_n = (2/3) 3^n mu_n` on 3-units), `Z_n(q) = sum_c mu_n(c)^q`, `tau(q) = liminf
-log Z_n(q)/(n log 3)`, `D_q = tau(q)/(q-1)`.

**Theorem 3 (bound).** For every `q >= 1`, `tau(q) <= min(q - 1, log_3(2^q - 1))`. Hence `D_q <= 1` for `q > 1`,
`D_q <= log_3(2^q - 1)/(q - 1)` for `q >= 2`, and `D_inf <= log_3 2`.
*Proof.* `mu_n(c) = sum_(w: Y_n(w) = c) 2^(-S_n(w))` over the parity words `w = (a_1, ..., a_n)`. For `q >= 1`,
`(sum_i m_i)^q >= sum_i m_i^q`, so `Z_n(q) >= sum_w 2^(-q S_n(w)) = (sum_a 2^(-qa))^n = (2^q - 1)^(-n)`, giving
`tau(q) <= log_3(2^q - 1)`. Jensen on the `(2/3) 3^n` residues of positive mass gives `Z_n(q) >= ((2/3)
3^n)^(1-q)`, i.e. `tau(q) <= q - 1`. □ Besides the normalization equality at `q = 1`, the two bounds cross for `q > 1` exactly at `q = 2`, where `2^2 - 1 = 3 = 3^(2-1)`:
`sum_a 4^-a = 1/3` is the incoherent (word-energy) rate and `3^-1` the Parseval rate, the same coincidence that makes
H1's margin `1.3%` in the Fourier note. For `q < 2` the second bound is the smaller one and the repo's P7
(`E[rho_n^q]` bounded for `q < 2`, SUPPORTED to level `14`) would imply it is attained for `1 < q < 2`: `D_q = 1`; the finite moment data alone do not prove this. For `q > 2` the
word-energy bound is the smaller one.

**Conjecture (equality).** `tau(q) = log_3(2^q - 1)` for `q >= 2`: in the concentrated regime each residue is carried
by few parity words at the `q`-tilted scale. Then `D_2 = 1` with `E[rho_n^2]` linear (the observed `+0.31` per
level), `D_2.5 = 0.9335`, `D_3 = 0.8856`, `D_4 = 0.8217`, `D_inf = log_3 2 = 0.6309` (the atom at `-1`:
`mu_n(-1) 2^n -> 1.462`, so the local dimension there is exactly `log_3 2`), and the Legendre transform gives
`alpha(q) = 2^q ln 2 / ((2^q - 1) ln 3)` on `[log_3 2, 4/(3 log_2 3)]` with `f = q alpha - tau`.

**Numerical partition functions to level `16`** (`collatz_syracuse_multifractal_20261005.py`; `43` million residues at
`n = 16`, valuations to `60`, total mass `1` to twelve digits; per-level slope `-log(Z_n/Z_(n-1))/log 3`):

| `q` | bound `tau(q)` | slope `n = 14` | `15` | `16` | bound `D_q` |
|---|---|---|---|---|---|
| 0.5 | `-0.5000` | `-0.4978` | `-0.4981` | `-0.4983` | `1` |
| 1.5 | `0.5000` | `0.4855` | `0.4868` | `0.4880` | `1` |
| 1.9 | `0.9000` | `0.8546` | `0.8575` | `0.8601` | `1` |
| 2.0 | `1.0000` | `0.9430` | `0.9463` | `0.9492` | `1` |
| 2.5 | `1.4003` | `1.3586` | `1.3630` | `1.3669` | `0.9335` |
| 3.0 | `1.7712` | `1.7356` | `1.7401` | `1.7440` | `0.8856` |
| 4.0 | `2.4650` | `2.4257` | `2.4297` | `2.4332` | `0.8217` |
| 6.0 | `3.7712` | `3.7260` | `3.7291` | `3.7317` | `0.7542` |
| 8.0 | `5.0439` | `5.0038` | `5.0060` | `5.0078` | `0.7206` |

Every slope is below its bound and rising, with increments of `1/n` type (`q = 2`: the slope is `1 - (log
E[rho_n^2]/E[rho_(n-1)^2])/log 3 ~ 1 - 1/(n log 3)`, as the linear growth of `E[rho^2]` predicts). The `q = 2.5`
gap `0.033` closes at `0.004` per level; a `1/n` extrapolation overshoots the bound, a geometric one lands on it. So
the data are consistent with equality and do not decide it; the repo's ridges (several parity words on one
residue) would make `tau` strictly smaller, i.e. the law more multifractal than the word-energy bound.

**Proposition 4 (the critical frequency).** `2 - H(3/4) = (3/4) log_2 3` exactly, and `p = 3/4` is the tangency
point of `(2 - H(p))/(p log_2 3)` with `1`. *Proof.* `2 - H(3/4) = 2 + (3/4)(log_2 3 - 2) - 1/2 = (3/4) log_2 3`;
the derivative condition `-p H'(p) = 2 - H(p)` with `H'(p) = log_2((1-p)/p)` holds at `p = 3/4` since
`(3/4) log_2 3 = 2 - H(3/4)`. □ In words: at the `q = 2` tilt the words that carry the second moment have
valuation-one frequency `3/4` (mean valuation `4/3`, the step law `3 4^-a`); their 3-adic local dimension is
`4/(3 log_2 3) = 0.841` and the dimension of that set of words is `H(3/4)/((3/4) log_2 3) = 0.682`. The Cauchy
criticality note's "`s = 2` tilt with odd frequency `3/4`" and the ridge analysis's "energy-weighted mean step
`4/3`" are this Legendre point; the identity makes it exact.

## 4. Hot spots: where the Syracuse law concentrates (VERIFIED (finite numerical))

The local dimension of `mu` at the 3-adic shadow of an integer `m` is bounded above by its parity cylinder:
`alpha(m) <= log_3 2 times (mean halvings per odd step along the orbit of m)`. Measured `-log_3 mu_n(m mod
3^n)/n` at level `14`, with the cylinder bound:

| point | orbit structure | bound | `alpha_14` |
|---|---|---|---|
| `-1` | fixed, valuation `1` | `0.631` | `0.606` (the maximal atom at every level; `1/n` correction `-log_3 1.46/n`) |
| `-5` | 2-cycle, word `(1, 2)` | `0.946` | `0.873` |
| `-17` | 7-cycle, `11` halvings | `0.992` | `0.917` |
| `1` | trunk, valuation `2` | `1.262` | `1.049` |
| `5`, `7` | enter the trunk | -- | `1.045`, `1.108` |

The three negative cycles are the hot spots (local dimension below `1`, mass denser than Haar), in the order of
their non-descent; the positive integers tested have finite-level dimension estimates about `1`; no finite-density limit is proved by those values. The
measured values at `-5` and `-17` are below the single-cylinder bound: other parity words land on their residues
too, which is the coherent landing the ridge analysis measures. In the pasted paragraph's image, the fault system
is the `3x-1` sheet with its three attractors (THM-4517's basins `0.327 / 0.325 / 0.348`), and the earthquake hot
spots of the `3x+1` law are its cycles.

## 5. The Vitali wall, algorithmic form (PROVED remark + typed reading)

**Remark 5 (repaired proof).** For every integer m its parity vector x(m)
is computable, so `K(x(m) restricted to n) <= K(m)+O(log n)` and its
effective dimension is0. If mu is computable and mu{x(m)}=0, computably
search for prefixes with cylinder mass below2^-j; these cylinders form a
Martin-Lof test covering x(m). Thus the point is not mu-random. Mere
divergence of minus-log mass, compared with an O(log n) complexity upper
bound, would not by itself prove this last conclusion. For a computable
atomic mu with positive mass at x(m), the point is mu-random instead. □

An upgrade from almost everywhere to every point random for a fixed
computable atomless measure does not directly cover the specified integer
itineraries. The cited density, moment, and Fourier results retain their
individual hypotheses and scope; they do not imply a blanket prohibition
on uniform estimates or a change of measure. In particular the target-law
finite moments do not by themselves prove asymptotic moment boundedness.
The coefficient-stopping target `Z^+ cap E_inf = empty` is necessary for
Collatz, while actual endpoint exceptions and the terminal cycle remain
separate obligations. The negative point -1 is one explicitly observed hot
spot, not the only computable point a measure can see. The two-sheet guard
also has integer consequences: it certifies common futures for the sources
satisfying it, without giving universal root coverage.

**What a crossing would need (repaired).** A typical-point theorem for the
fixed atomless measure does not automatically settle a specified integer.
One possible reroute is a computable atomic source measure, under which
mass-one root coverage implies every integer is covered; the necessary
decay estimate is still open. Another is a uniform arithmetic bound that
does not require randomness. Eventual periodicity alone permits additional
cycles and is weaker than entry into the specified root cycle. No actual
Vitali set or nonmeasurable selector has been constructed here. The older
"Vitali atom" and "Vitali covering" connections are finite-observer analogies.

## 6. Reproduction

```bash
python 04-computation/experiments/collatz_syracuse_multifractal_20261005.py 16 60    # 72 s, 1.2 GB at level 16
python 04-computation/experiments/collatz_syracuse_multifractal_20261005.py 13 60    # 2 s
```
Output `.out` beside the script; the local dimensions of section 4 are the six atoms printed by the same recursion
at level `14`.

## 7. Verdicts

| claim | status |
|---|---|
| the cited box-counting exponent of E_inf matches the Eggleston frequency entropy; the all-prefix Hausdorff lower bound is not established here | CITED counting input + scoped analogy |
| coefficient contraction iff `2^(-S_k) < 3^(-k)`; actual descent also needs the carry threshold (Proposition 1) | PROVED, corrected |
| `tau(q) <= min(q-1, log_3(2^q-1))`, crossing at `q = 2` by `sum 4^-a = 1/3` (Theorem 3) | PROVED |
| `2 - H(3/4) = (3/4) log_2 3`; the `q = 2` Legendre point is frequency `3/4` (Proposition 4) | PROVED |
| numerical partition functions to level `16` below the bound with `1/n` corrections | VERIFIED (finite numerical) |
| equality `tau(q) = log_3(2^q - 1)` for `q > 2` | CONJECTURED |
| hot spots `-1, -5, -17` at local dimensions `0.61, 0.87, 0.92`; positive integers at about `1` | VERIFIED (finite numerical) (level `14`) |
| integers have effective dimension `0`; the a.e.-to-random upgrade cannot reach them (Remark 5) | PROVED (trivial) + typed reading |
| Collatz | OPEN |
