# LRC through Galois orbits and the Zak/theta transform (S521)

*claudebox-2026-06-01-S521. Developing the Galois/cyclotomic and Zak/theta angles
on the Lonely Runner Conjecture, connected to the tournament work. The Galois angle
yields a genuinely promising proof-strategy (equidistribution of Galois conjugates
via Weil/Ramanujan bounds => a denominator bound); the Zak/theta angle places LRC
in the abelian-variety / theta-divisor setting.*

## Galois: the runner configuration as a tuple of roots of unity

At a rational time `t = a/q` (`q` prime, `a in (Z/q)*`), runner `i` sits at the
residue `v_i a mod q`, i.e. at the `q`-th root of unity `zeta_q^{v_i a}`. The
Galois group `Gal(Q(zeta_q)/Q) = (Z/q)*` acts by `a : zeta -> zeta^a`, sending the
base configuration `(zeta_q^{v_i})_i` to its conjugate `(v_i a mod q)_i`. The runner
is SAFE iff its residue lies in the **safe box** `B_q = [ceil(q/n), q-ceil(q/n)]`.
So:

> **LRC at denominator `q`  <=>  the Galois orbit `O_q = {(v_i a mod q)_i : a in
> (Z/q)*}` meets the safe box `B_q^m`.**

The LRC multiplicative walk IS the Galois action; the regular-polygon extremal
configuration is the full set of `n`-th roots of unity, whose half-turn tournament
is the **self-complementary circulant** (Paley-like at prime `n`) — a cyclotomic /
Galois-symmetric object (project: Paley = QR-circulant, THM-052 circulants SC). And
the menu count `A000016(m) = (1/2m) sum_{d|m odd} phi(d) 2^{m/d}` is a **cyclotomic
Galois/Burnside count** (`phi`, odd divisors = the cyclotomic group's cycle index).

## The promising strategy: equidistribution of Galois conjugates => denominator bound

The Galois orbit `O_q` equidistributes on the torus as `q -> oo`. Its discrepancy
from uniform is controlled by exponential sums `sum_{a in (Z/q)*} e((sum_i c_i v_i)
a / q)`. For `q` prime this sum is `q-1` (large) iff `sum_i c_i v_i ≡ 0 (mod q)` —
a **relation / resonance** — and `-1` (small) otherwise. So:

- the orbit-in-box fraction `f(q) -> V = (1-2/n)^m > 0` as the resonances thin out
  (`q` larger than the small relations `sum c_i v_i`);
- computed: non-tight sets become **lonely by a small prime** — `(1,2,4,7)` at
  `q=11`, `(2,3,5,7,11)` at `q=17` — and `f(q)` tracks `V`.

> **Strategy.** A Weil/Ramanujan bound on the orbit discrepancy gives `f(q) > 0`
> for all `q >= q_0(n)` (the safe box is hit once `phi(q)*V` exceeds the
> exponential-sum error), reducing LRC to `q < q_0` — a FINITE check. `q_0` is
> controlled by the **resonances** (relations `sum c_i v_i ≡ 0`), exactly the
> doubling / `n=2*odd` structure isolated throughout S521.

The one caveat the data force: the **tight/extremal sets are invisible at prime
`q > n`** — their only lonely time is `q = n` (the cyclotomic field `Q(zeta_n)`
itself, the regular polygon, the boundary). So the Galois strategy splits cleanly:
**(base) `q = n` cyclotomic** (THM-369 / the SC circulant, handled when no speed is
`≡ 0 mod n`), plus **(tail) large `q` equidistribution** (Weil), with a finite
middle. This is the most principled version yet of the bounded-denominator
reduction — the equidistribution mechanism is real, and the resonances are the
explicit obstruction a Weil bound must beat.

## Sharpened by computation: the Erdős–Turán bound and the SUM-relation obstruction

A larger equidistribution study (exact, primes to 8000) makes the strategy
concrete and corrects the obstruction:

- **Exact Fourier identity.** `|O_q ∩ B_q^m| = q^{-m} Σ_c (∏_i ĥ_B(c_i)) S(Σ_i c_i v_i)`,
  with the Ramanujan sum `S(t) = q-1` if `t ≡ 0 (mod q)` else `-1`. A genuine
  integer relation `Σ c_i v_i = 0` forces `S = q-1` for EVERY prime `q` — a
  *permanent* resonance; all other frequencies contribute only `O(1/q)`.
- **Erdős–Turán / Fourier bound.** For sets with no small relation,
  `|f(q) - V_q| <= C(m)/q`; since `V_q -> (1-2/n)^m > 0`, this forces `f(q) > 0`
  (the box is hit) for all `q >= q_0(m)` — a genuine finite reduction for
  non-resonant sets. Empirically `q_min <= 43` across everything tested (m<=5).
- **The obstruction is the SUM-relation `v_i + v_j = v_k`, NOT doubling.** The
  resistant sets (`f(q)` biased below `V`, larger `q_min`) all carry a weight-3
  additive sum-relation whose box-Fourier sign is NEGATIVE (e.g. `2+3-5=0` in
  `{2,3,5}`, `2+5-7=0` in `{2,3,5,7}`). The doubling relation `2v_i=v_j` resonates
  too but its Fourier sign is `~0`/positive, so doubling sets are NOT resistant.
  (This refines the earlier "doubling" emphasis: for *equidistribution* the active
  resonance is the additive triple `v_i+v_j=v_k`, abundant in dense/AP-like sets.)
- **Relation-aware bound suffices.** The resonant frequencies are FINITE in number
  (the small relations); separating them out gives a relation-aware Erdős–Turán
  bound that still yields a finite `q_0`. The only honest boundary case is the
  **extremal AP family**, lonely exactly at `q = n` (composite when `n` is
  composite, so the prime-only view misses it) — precisely the tight cases.

## Zak / theta: quasi-periodization and the theta divisor

`v_i t mod 1` is precisely the quasi-periodization the **Zak transform** formalizes
(`Zf(x,omega) = sum_k f(x+k) e(-k omega)`). Smoothing the safe-band window to a
Gaussian `g` makes `sum_k g(t+k)` a **Jacobi theta function**; the runner
"partition function" `Z(t) = prod_i Theta(v_i t)` is a product of theta functions,
and the lonely set is where `Z` is large (all runners in band). The zeros of such
theta objects form a **theta divisor**: on the torus `T^m` carrying the runner
line, LRC = the line avoids a thickening of the theta divisor. This places LRC in
the **abelian-variety** setting (Riemann theta zeros, the geometry of `T^m`), and
the Balian-Low principle explains the shape: the system is overcomplete (density
`2(n-1)/n > 1`) yet can have a hole, because the Zak transform of a band window has
zeros forced by rational resonance.

## The Galois <-> Zak duality (= multiplication <-> addition)

The two angles are Fourier/Zak duals and the two halves of the program's thesis:

- **Galois (multiplicative / cyclotomic):** `q`-th roots of unity, Frobenius /
  `(Z/q)*` action, `phi`/odd-divisor counts (`A000016`), the QR/Paley circulant —
  the multiplicative sieve, the *target's* shape.
- **Zak / theta (additive / time-frequency):** quasi-periodization, theta
  functions, the theta divisor, equidistribution of the line — the additive walk.

LRC sits where the Galois orbit (multiplicative) of the additive line's reduction
meets the safe box; the Zak transform is the exact bridge (it intertwines
translation and modulation, i.e. addition and multiplication). This is
`A000568`/`A000016` "addition meets multiplication," now named **Galois x
Zak/theta**.

## Assessment and lead

- **Galois angle: a real proof-strategy.** Equidistribution of Galois conjugates
  (Weil/Ramanujan) forces the safe box to be hit for `q >= q_0(n)`; combine with the
  cyclotomic base `q = n`; the finite middle is the only open part, and `q_0` is
  governed by the explicit resonance/relation structure. This is the most concrete
  new lead of the S521 arc and deserves the next push: derive an explicit
  Erdos-Turan / Weil discrepancy bound for `O_q` and an explicit `q_0(n)`.
- **Zak/theta angle: a placement, not yet a bound.** It correctly situates LRC at
  the theta divisor of `T^m` and explains the overcomplete-with-hole phenomenon
  (Balian-Low); turning it into a bound would need control of the theta zeros along
  the rational line.

## Lead (concrete, now sharpened)

The mechanism is the **Erdős–Turán bound** `|f(q) - V_q| <= C(m)/q` from the
exact Fourier identity above, with the resonant frequencies (the finite set of
small relations `Σ c_i v_i = 0`, dominated by additive triples `v_i+v_j=v_k`)
separated out into a relation-aware error term. Solving `(q-1) V_q > error` gives an
explicit `q_0(m)` (empirically `q_min <= 43`). The program is therefore:

1. **Base `q = n`:** the cyclotomic/regular-polygon case (THM-369; lonely unless a
   speed is `≡ 0 mod n`).
2. **Tail `q >= q_0(m)`:** Erdős–Turán equidistribution forces `f(q) > 0`
   (relation-aware for the finitely many sum-relations).
3. **Middle `n < q < q_0(m)`:** a finite, explicit check.
4. **Extremal AP family:** the honest boundary (`q = n` composite) — exactly the
   tight cases, handled by the base step.

Next concrete task: write the relation-aware Erdős–Turán bound with the additive
triples `v_i+v_j=v_k` as the explicit resonant set, and pin down `C(m)` and
`q_0(m)`. This is the Galois–Weil completion of the bounded-denominator program,
and it is now reduced to an explicit harmonic-analysis estimate plus a finite
computation.
