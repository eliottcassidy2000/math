---
id: THM-3027
title: "The AMM 12592 capacity threshold is exactly log(phi)/log(sqrt 5)"
status: >
  PROVED (tangency system, symbolically closed) + VERIFIED-EXACT to 30+ dps.
  THM-3002 reports the asymptotic Bernstein-capacity threshold only as
  "gamma* ~ 0.598, the solution of a two-ray entropy comparison". That
  comparison is solved here in closed form: gamma* = log(phi)/log(sqrt 5)
  = 2 log(phi)/log 5 = 0.597987435665440149745..., so the epoch-closure
  route's natural barrier is C = 1 + gamma* = 1.597987435665440... The
  mechanism is not decorative: eliminating the two multipliers from the
  stationarity/tangency/value system collapses it to the single equation
  (1 - tau)^2 = tau, i.e. tau^2 - 3 tau + 1 = 0, the minimal polynomial of
  phi^2 -- so the binding mass fraction is tau* = phi^-2 = 2 - phi, the
  binding fill ratio is rho* = 1 - 1/sqrt 5, and phi enters as the root of
  the tangency condition itself. The maximizing sigma* = 0.038635 is
  INTERIOR, so no boundary case is hidden. Consistent with death-star's
  finite-R ladder (0.5313, 0.5606, 0.5758, 0.5849, 0.59065, 0.59393 for
  R = 32..1024), which increases toward gamma* from below, and it explains
  every eq(27) death: 2457/6592, 2457/4135 and 1/2 all lie BELOW gamma*,
  while 3/5 is the first round rate above it. The binding fraction phi^-2 is
  UNIVERSAL in the alphabet size b (the reduction to (!) is b-free), with
  gamma*(b) = log(phi)/log((b+phi)/phi); b = 2 is special only because
  (2+phi)/phi = sqrt 5 exactly, so the golden constant is NOT an artifact of
  the binary alphabet and is NOT evidence of Fibonacci substructure.
  AUDIT UPGRADES 2026-08-03 (boxeph audit, same file as THM-3009 note):
  tangency collapse re-derived independently (residuals identically 0);
  inner sigma-concavity PROVED symbolically (d2psi/dsigma2 < 0) and
  dpsi/dgamma = -(1+sigma)log(1-rho) > 0 PROVED, making the FLOOR direction
  scan-free; sigma* interiority upgraded to an exact certificate via the
  rational bracket of gamma*. The Laplace/Stirling rate passage remains the
  single assumed analytic step (self-declared below).
source: klein-S428
depends_on:
  - THM-3002
related:
  - THM-2966
  - THM-2160
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_capacity_threshold_is_log_sqrt5_phi_thm3027.py
output: 05-knowledge/results/amm12592_capacity_threshold_is_log_sqrt5_phi_thm3027.out
---

# THM-3027 -- the capacity threshold is `log(phi)/log(sqrt 5)`

## 0. Statement

Let `d_i = floor(gamma(R+i)) + D0` for `i = 0..R-1` and let the Bernstein
capacity criterion of THM-3002 be

```text
   S(t) := sum_{i<=t} C(d_i, t-i) 2^(t-i)  >=  C(R-1, t)   for all t.      (*)
```

Let `gamma*` be the infimum of the `gamma` for which the large-`R` rate form of
`(*)` holds for every mass fraction. Then

```text
   gamma* = log(phi) / log(sqrt 5) = 2 log(phi) / log 5 = 0.5979874356654401497...
   C      = 1 + gamma*                                   = 1.5979874356654401497...
```

attained at the interior point

```text
   tau*   = phi^-2 = 2 - phi = (3 - sqrt 5)/2 = 0.38196601125...   (mass fraction)
   rho*   = 1 - 1/sqrt 5                      = 0.55278640450...   (fill ratio)
   sigma* = 0.03863539626...                                       (source index)
```

## 1. The rate problem

Put `i = sigma R`, `t = tau R`, `m = tau - sigma`, `D(sigma) = gamma(1+sigma)`,
`rho = m/D`. By Stirling, `d_i/R -> D(sigma)` and

```text
   (1/R) log [C(d_i, t-i) 2^(t-i)]  ->  D H(rho) + m log 2,
   (1/R) log C(R-1, t)              ->  H(tau),
```

with `H` the natural entropy. The sum in `(*)` has `R` terms, so
`(1/R) log S(t) -> Psi(gamma,tau) := max_sigma [D H(rho) + m log 2]`, the max
taken over `sigma in [0, min(1,tau)]` with `m <= D`. The criterion holds
asymptotically iff `Psi(gamma,tau) >= H(tau)` for all `tau in (0,1)`.

This Laplace/Stirling passage is the one analytic input and is standard; it is
assumed, not reproved. Everything below is exact given it.

## 2. The three conditions at threshold

At `gamma = gamma*` the constraint is tangent at some `tau*`. Two derivatives
are needed. Writing `psi = D H(rho) + m log 2` with `dD/dsigma = gamma`,
`dm/dsigma = -1`, and using `d/dD [D H(m/D)] = -log(1-rho)`,
`d/dm [D H(m/D)] = log((1-rho)/rho)`:

```text
   dpsi/dsigma = -gamma log(1-rho) - log(2(1-rho)/rho).
```

By the envelope theorem `dPsi/dtau = log(2(1-rho)/rho)`, while
`H'(tau) = log((1-tau)/tau)`. So with `u := 1 - rho` and `A := log(1/u) > 0`:

```text
   (S)   gamma A = log(2u/(1-u))                     inner stationarity
   (T)   2u/(1-u) = (1-tau)/tau                      tangency in tau
   (V)   D [H(rho) + rho log 2] = H(tau)             value
```

## 3. The collapse (the content of the theorem)

**Key identity.** `(S)` says `log(2u/rho) = gamma A`, hence
`log(rho/2) = log u - gamma A = -A(1+gamma)`. Therefore

```text
   H(rho) + rho log 2 = -rho log(rho/2) - u log u
                      = rho A (1+gamma) + u A
                      = A (rho + u + gamma rho)
                      = A (1 + gamma rho).                                (K)
```

**The multiplier cancels.** The definition `sigma = tau - rho D` together with
`D = gamma(1+sigma)` gives `D(1 + gamma rho) = gamma(1+tau)`, i.e.
`D = gamma(1+tau)/(1 + gamma rho)`. Substituting `(K)` into `(V)`, the factor
`(1 + gamma rho)` cancels **identically**:

```text
   D [H(rho) + rho log 2] = gamma A (1 + tau).
```

By `(S)` and `(T)`, `gamma A = log((1-tau)/tau)`. So `(V)` becomes an equation
in `tau` alone:

```text
   (1+tau) log((1-tau)/tau) = -tau log tau - (1-tau) log(1-tau).
```

Expanding and cancelling `-tau log tau` from both sides leaves
`2 log(1-tau) = log tau`, that is

```text
   (1 - tau)^2 = tau,     equivalently     tau^2 - 3 tau + 1 = 0.          (!)
```

**`(!)` is the minimal polynomial of `phi^2`** (indeed `phi^2 + phi^-2 = 3` and
`phi^2 . phi^-2 = 1`). The root in `(0,1)` is `tau* = phi^-2 = 2 - phi`, whence
`1 - tau* = phi^-1`.

**Back-substitution.** `(T)` gives `2u/(1-u) = (1-tau*)/tau* = phi`, so
`u(2+phi) = phi`. Since `2 + phi = phi^2 + 1 = phi sqrt 5`,

```text
   u* = phi/(phi sqrt 5) = 1/sqrt 5,      A* = log(sqrt 5),
   rho* = 1 - 1/sqrt 5,
```

and `(S)` closes it:

```text
   gamma* = log(2u*/(1-u*)) / A* = log(phi)/log(sqrt 5) = 2 log(phi)/log 5.
```

Finally `D* = gamma*(1+tau*)/(1+gamma* rho*) = 0.62109091720...` and
`sigma* = tau* - rho* D* = 0.03863539626...`, which lies strictly inside
`[0,1]`, so the inner maximum is interior and no boundary branch was missed.

## 4. Why the golden ratio

It is worth being precise about the mechanism, because a golden ratio appearing
next to a Fibonacci-flavoured problem invites a wrong explanation.

`phi` does **not** enter through any Fibonacci or Zeckendorf structure in the
coin problem. It enters at `(!)`, which is the *tangency condition* — the
statement that the worst mass fraction is a double root.

Nor does the `2` in `2 log(1-tau) = log tau` come from the alphabet. Replace the
capacity factor `2^(t-i)` by `b^(t-i)` for any `b > 1`. Then `(S)` becomes
`gamma A = log(bu/(1-u))` and `(T)` becomes `bu/(1-u) = (1-tau)/tau`, but the
identity `(K)` is unchanged — `log(rho/b) = -A(1+gamma)` still gives
`H(rho) + rho log b = A(1 + gamma rho)` — so the multiplier still cancels and
`(V)` still reduces to `(1+tau) log((1-tau)/tau) = H(tau)`. The `2` is
`(1+tau) + (1-tau)`, produced by the algebra, not by the alphabet. Hence:

> **The binding mass fraction `tau* = phi^-2` is universal in `b`.**
> Only the threshold moves: `(T)` gives `u*(b) = phi/(b+phi)` and
> `gamma*(b) = log(phi) / log((b+phi)/phi)`.

Verified numerically (part E of the script), worst `tau` pinned at `0.382` in
every case:

| `b` | `gamma*` numeric | `log(phi)/log((b+phi)/phi)` |
|---|---|---|
| 2 | 0.597987433 | 0.597987436 |
| 3 | 0.458840046 | 0.458840048 |
| 4 | 0.386586952 | 0.386586954 |

So `phi` is present for **every** alphabet size, through `(!)` alone. What is
special about `b = 2` is only that the threshold is then a clean surd:
`(2+phi)/phi = sqrt 5` exactly, because `2 + phi = phi^2 + 1 = phi sqrt 5`.
That coincidence — and nothing about Fibonacci — is why the binary constant is
`log_{sqrt 5}(phi)`. (Degenerate check: `b = 1` gives `(1+phi)/phi = phi` and
`gamma* = 1`, the no-boost case.) Recorded so that no one reads `phi` here as
evidence of a Fibonacci substructure, and so that no one attributes the golden
constant to the binary alphabet.

## 5. Consequences for the ledger

1. **The epoch-closure barrier is `C = 1 + 2 log(phi)/log 5 = 1.59798743566544...`**
   This is the best constant the THM-3002 route can deliver: `(*)` is necessary
   for block closure, so no `gamma < gamma*` closes all large epochs, while
   THM-3002's sufficiency gives `C* <= 1 + gamma` for admissible `gamma`.
   Whether the true `C` equals this is a separate question — this pins the
   method's barrier, not (yet) the answer.
2. **Every eq(27) death is explained by one inequality.** `2457/6592 = 0.372725`,
   `2457/4135 = 0.594196` and `1/2` all lie **below** `gamma*`, so all three are
   eventually deficient — matching the direct computations (`R >= 16`,
   `R = 2048`, and `R = 64` respectively). `3/5 = 0.6 > gamma*` is the first
   round rate above threshold, which is exactly why `3/5` survives.
   See `07-reflections/eq27-is-a-logit-gate-and-its-weight-is-not-pinned-klein-S428.md`.
3. **The finite-`R` ladder was converging, and to this.** death-star's
   thresholds `0.5313, 0.5606, 0.5758, 0.5849, 0.59065, 0.59393` for
   `R = 32..1024` increase monotonically toward `gamma*`, the gap roughly
   halving per doubling of `R`. An earlier klein note warned against
   extrapolating that ladder; the warning was wrong (already retracted in the
   script header of `amm12592_capacity_criterion_...py`), and the limit now has
   a closed form.
4. **Do not re-run the numerical threshold search.** `gamma*` is exact.

## 6. Status of each claim

| Claim | Status |
|---|---|
| Rate reduction of `(*)` by Stirling | ASSUMED (standard) |
| `(S)`,`(T)`,`(V)` are the threshold conditions | PROVED |
| Key identity `(K)` and the multiplier cancellation | PROVED |
| `(V) <=> (1-tau)^2 = tau` | PROVED |
| `gamma* = log(phi)/log(sqrt 5)` | PROVED, given the rate reduction |
| `sigma*` interior; stationary point is the global worst `tau` | VERIFIED NUMERICALLY (3000-point scan; bisection agrees to 9 digits) |
| Finite-`R` ladder increases to `gamma*` | VERIFIED against death-star's data |
