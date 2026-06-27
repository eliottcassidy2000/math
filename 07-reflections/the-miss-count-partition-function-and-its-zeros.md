# The miss-count partition function and its zeros: a new signal for what makes consec extremal

*mac-mini-2026-06-27-S66. Owner: understand what maximizes LRC values ↔ tournaments, be comprehensive on
the niche details, then be BOLD with hypotheses and wild guesses about the recurring structure, and let it
inspire NEW signals to measure. Comprehensive scouts surfaced the key recurring object (Route 5,
INVESTIGATION-BACKLOG:1520) and ~20 niche details; this reflection turns them into one new object, a verified
new signal, and a slate of bold guesses + signals to measure.*

## The recurring structure (comprehensive finding)
Across **every** functional the project measures — `μ_{1/7}` minimization, `S_1` minimization, `L_yK8`/`S2`
maximization, `meas(S7)` maximization — **consec/AP is the extremizer**, and the proof is "irreducibly global"
(not local descent, not rearrangement, not monotone). My S65 result (the config improvement tournament is
**non-transitive**) is the same wall from the tournament side. The deepest framing found (Route 5): LRC and
tournament-`H` are **parallel fugacity partition functions** with consec extremal on the *same additive AP
surface but at opposite level + opposite fugacity sign*:
- tournament `H = Σ_j (+2)^j α_j` (odd-cycle packings), consec/AP wins the **high-order** levels;
- LRC `meas(S7) = Σ_k (−1)^k MISS_k` (sector misses), consec wins the **low-order** level.

Honest correction to a tempting over-read: these are **two parallel partition functions on different data**
(odd-cycles vs sectors), not one object — the unification is a shared *AP-extremal surface*, not an identity.

## The new object: the miss-count partition function `G_N(z) = E[z^N]`
Since `MISS_k = S_k` (factorial moments) and `Σ_k S_k x^k = E[(1+x)^N]`, the LRC side IS the **probability
generating function of the sector-miss-count** `N = #empty inner sectors`:
```
G_N(z) = Σ_{t=0}^{6} q_t z^t ,   q_t = P(N=t),   z = 1 + x.
   z = 0  ->  q_0 = p0  (LRC coverage)        z = 3  ->  E[3^N]  (the tournament fugacity).
```
**The project has only ever measured the moments `S_r` and the single value `p0 = G_N(0)`.** It never
measured `G_N` as a *whole function* — its analytic continuation, and above all **its zeros in the complex
`z`-plane**. That is the new object.

## The new signal (VERIFIED): the zeros of `G_N(z)`
`lrc_missPGF_new_signal_macmini_S66.py` + `..._realroots_vs_extrememass_macmini_S66.py`:
- **consec (and its dilation even-AP): ZERO real roots** — three complex-conjugate pairs `|z|=1.49,1.58,1.70`
  at args `±55°, ±111°, ±154°`, all far from `z=0`. (`Π|z_i| = q0/q6 = 16.03`, Vieta.)
- spread/covering/break configs: a real root **very close to `z=0`** (`|z*|=0.05,0.11,0.12`) and 2–4 real roots.
- **The signal stratifies config space, and the extremizer lives in the `#real=0` stratum**: over 250 configs,
  `max L_yK8 = 3.58` (= consec) occurs ONLY at `#real=0`; the `#real=2,4` strata cap at `L_yK8 ≈ 0.7–1.0`.
  `corr(#real-roots, extreme-mass q0+q6) = −0.37`.

### Why this is the right signal (the rigorous reading)
A PGF with non-negative coefficients is **real-rooted ⟺ `N` is a sum of INDEPENDENT indicators** (the
classical Pólya-frequency / Newton theorem). So:
> **`#real roots of `G_N` measures sector-INDEPENDENCE.** consec, with **0 real roots**, is the config whose
> sector-miss-count is **maximally NON-independent** — maximal positive correlation of the seven sectors'
> emptiness. That maximal correlation IS the high extreme-mass `q0+q6` that `L_yK8 = 10(q0+q6)+q3` rewards.
The recurring "consec is extremal" becomes: **consec pushes the partition-function zeros maximally OFF the
real axis** — it is the most-correlated, least-factorizable config. This unifies the S60 gK8/S2/Clebsch
finding (consec maximizes pairwise sector covariance) with the analytic zero structure.

## Bold hypotheses / wild guesses (marked ★ wild)
1. **Lee–Yang confinement.** `G_N(z)` is a partition function; its zeros are Lee–Yang zeros. *Hypothesis:* the
   LRC extremizer is exactly the config whose zeros are **maximally confined off the real axis** (a Lee–Yang
   property). Proving "the extremizer has no real zero on `[−1,0)`" might be a NEW analytic route to coverage
   extremality — turning a combinatorial optimization into a zero-free-region statement.
2. **★ The zero curve is governed by the apex 7.** consec's zero arguments `±55°,±111°,±154°` cluster near
   multiples of `≈ 360/7 ≈ 51°`. *Wild guess:* the extremizer's zeros lie on an arc set by the 7th roots of
   unity / the apex — the Fano/QR-7 structure surfacing in the complex `z`-plane. (To test: does the zero
   pattern converge to `7th-root` arguments as the row varies?)
3. **★ The `#real`-root count is a phase label, and `0→2→4` is the phase transition.** The k=8 "break"
   (minimizer jumps top-cluster → middle-spread `{1,5,7,8,9}`) is a known discontinuity; *guess:* it is a
   change in the real-root count of `G_N` (a real zero colliding onto/off the real axis = the break).
4. **★ Real-rootedness = the matching-polynomial (Heilmann–Lieb) analog.** A graph's matching polynomial is
   real-rooted; the miss-PGF being complex-rooted says the sectors do NOT form an independent (matching-like)
   system. The *defect from real-rootedness* is a tournament/graph-theoretic interaction strength — possibly
   the same quantity as the conflict-graph `β₁` on the tournament side.

## The slate of NEW SIGNALS to measure (the owner's ask)
On the LRC side (all from `G_N(z)`, none previously tracked):
1. **`#real roots`** of the miss-PGF — the independence/correlation phase label (VERIFIED to stratify).
2. **nearest-zero distance to `z=0`** — "coverage-pole proximity" (large ⇔ high `p0`; consec `≈1.49`).
3. **the Lee–Yang confinement radius / the zero arc** — the region the zeros occupy; consec = all off-axis.
4. **the fugacity rank curve** — `rank(config)` by `G_N(z)` as `z: 0→3`; the crossover `z`.
5. **root-argument spread / regularity** — how evenly the zero arguments distribute (consec = regular arc).
6. **the discriminant / resultant of `G_N`** — vanishes exactly on the real-root transition locus (the phase
   boundary; candidate for the k=8 break).

On the tournament side (the parallel object, to measure and compare):
7. **the zeros of the winding tournament's independence polynomial `I(Ω,x)`** vs the miss-PGF zeros — do the
   two parallel partition functions share a zero structure at the AP extremal? (Bridges the LRC and tournament
   fugacity surfaces directly, where Route 5 only matched their argmax.)

## Honest status
VERIFIED: `G_N(z) = E[z^N]`; consec has 0 real roots; `#real` stratifies and the gK8 extremizer is in the
`#real=0` stratum; the independence-⟺-real-rooted reading is rigorous. BOLD/UNTESTED: the Lee–Yang
extremality route (#1), the apex-7 zero-curve (#2), the break-as-root-collision (#3), the `β₁` connection
(#4). The 6+1 new signals are defined and the first is verified; the rest are the measurement program.

Related: HYP-3103 (the verified signal), HYP-3085 (gK8/S2 = the correlation this measures), THM-577 (cap),
INVESTIGATION-BACKLOG:1520 (Route 5 fugacity parallel), [[the-four-faces-of-14-why-the-exceptional-structures-crowd-into-lrc]],
[[three-notions-of-sameness-are-the-lonely-sets-fiber]], [[tournaments-as-proof-engines-a-generative-catalog]].
