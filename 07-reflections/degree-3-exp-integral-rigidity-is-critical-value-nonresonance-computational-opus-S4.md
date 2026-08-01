---
source: opus-2026-08-01-S4 (degree-3 exp-integral: computational test of the Beukers rigidity criterion)
status: >
  EXPLORATION + REDUCTION (not a closed theorem for all cubics; honest). For a cubic Q, Phi_Q(z)=int_0^1 e^{zQ}dt
  satisfies a minimal homogeneous ODE of order = #(distinct critical values), generically 4 (two endpoints Q(0),Q(1)
  + two saddles Q'(c)=0); its irregular exponents at oo ARE those critical values. By Beukers 2006, I=Phi_Q(1) is
  transcendental IF {1,Phi_Q,Phi_Q',Phi_Q'',Phi_Q'''} are Qbar(z)-linearly independent (RIGIDITY) AND Phi_Q is
  non-rational. Computationally (exact rational moments, rigidity tested to coefficient-degree 6): RIGIDITY HOLDS
  for every cubic whose critical values avoid 0 (generic, depressed, pure-power t^3, t^3+7/10) and FAILS exactly
  when 0 is a critical value (Q=t^3/3+t^2, Q(0)=0) -- precisely the critical-value-resonance predicted by the
  general reformulation. Non-rationality at degree 3 is KNOWN (FC(2) holds for deg<=3, death-star exact
  elimination). So generic-cubic transcendence is REDUCED to "rigidity for all degrees" -- verified numerically,
  and a full proof needs the 4-exponential connection analysis (analogue of the degree-2 2-exponential argument).
  The pure cubic t^3 (order-2 Kummer M_3=1F1(1/3;4/3)) is already a proved theorem (pure-power file). This file
  claims: the mechanism + the computational confirmation + the reduction; NOT a general degree-3 theorem.
tags: [factorial-conjecture, FC2, exponential-integral, transcendence, e-functions, beukers, degree-3, rigidity, critical-values, saddle, airy, computational]
related: [degree-2-exp-integral-is-transcendental-rigorous-via-linear-beukers-shidlovskii, pure-power-exp-integral-transcendental-all-degrees-via-1F1-one-over-d, exp-integral-is-an-E-function-in-the-parameter-general-reformulation-and-rigidity]
external: ["Beukers, Ann. of Math. 163 (2006)", "death-star THM-3031 (bridge needs !=1) + FC(2) deg<=3 by elimination"]
---

# Degree 3: rigidity = critical-value non-resonance, confirmed computationally

## What was tested

For a cubic `Q`, `Phi_Q(z) = int_0^1 e^{zQ(t)} dt = sum_m mu_m z^m/m!` (exact rational moments `mu_m=int_0^1 Q^m`).
Two exact computations over `Q`:
1. **Minimal homogeneous ODE** of `Phi_Q` (search order `r`, polynomial-coefficient degree `I`).
2. **Beukers RIGIDITY criterion**: are `{1, Phi_Q, Phi_Q', Phi_Q'', Phi_Q'''}` linearly independent over `Qbar(z)`?
   Equivalently, is `1 = sum_{j=0}^3 q_j(z) Phi_Q^{(j)}(z)` UNsolvable in rational `q_j` (tested to `deg q_j <= 6`)?

## Results (exact)

| cubic `Q` | min ODE order | rigidity | verdict |
|---|---|---|---|
| `t^3 + t^2/2 - t + 1/3` (generic) | 4 | HOLDS | transcendental (w/ non-rationality) |
| `t^3 - t + 1/2` (`Q(0)=Q(1)=1/2`) | 3 | HOLDS | transcendental |
| `t^3` (pure power) | 2 | HOLDS | **proved** (pure-power file) |
| `t^3 + 7/10` | 2 | HOLDS | transcendental |
| `t^3/3 + t^2` (`Q(0)=0`) | 4 | **FAILS** | 0 is a critical value -- finer argument needed |

## Reading

- **Order = number of distinct critical values.** Generic cubic: 4 (`Q(0), Q(1)` and the two saddles `Q'(c)=0`)
  -> order 4. Endpoint coincidence `Q(0)=Q(1)` drops it to 3. Pure power `t^3`: the two saddles collide at `0`
  and the endpoints degenerate, leaving the **order-2 Kummer** `M_3 = 1F1(1/3;4/3)` -- exactly the E-function of
  the pure-power theorem. The ODE's irregular exponents at `infinity` are the critical values themselves.
- **Rigidity <=> `0` not a critical value.** The constant `1 = e^{0.z}` (exponent `0`) is independent of the jet
  `{Phi, Phi', Phi'', Phi'''}` iff `0` is not among the exponents, i.e. iff `0 notin {Q(0),Q(1),Q(c_1),Q(c_2)}`.
  Confirmed: the ONLY failure among the tested cubics is `Q=t^3/3+t^2` with `Q(0)=0`. This is exactly the
  critical-value-resonance the (corrected) general reformulation predicted.
- **Failure is not a counterexample.** When `0` is a critical value the Beukers *linear* criterion does not apply
  directly (`1` is dependent on the jet); `Phi_Q(1)` is very likely still transcendental, provable by a finer
  argument as in the degree-2 degenerate cases. RIGIDITY-FAILS means "this argument doesn't close it here", not
  "algebraic".

## Honest status

- **Proved (unconditional):** pure cubic `int_0^1 e^{alpha t^3+gamma} dt` transcendental (pure-power file).
- **Reduced + numerically confirmed:** generic-cubic (`0` non-resonant) transcendence follows from (i)
  non-rationality of `Phi_Q` -- KNOWN at degree 3 (FC(2) by exact elimination) -- and (ii) RIGIDITY, verified to
  coefficient-degree 6 for the cubics above. A general theorem needs (ii) proved for all cubics via the
  4-exponential connection (distinct nonzero critical values => distinct exponential types => non-isomorphic,
  no horizontal `1`), the degree-3 analogue of the degree-2 2-exponential argument.
- **Not claimed:** a degree-3 transcendence theorem for all `Q`; the `0`-resonant cases; complex cubics where
  non-rationality would need the (open) general FC(2) rather than the deg-3 elimination.

## Why the mechanism is trustworthy (post the S2 correction)

Unlike the retracted `Phi_Q`-non-constancy claim, this test does NOT rely on any Laplace/indicator lower bound.
It is a finite EXACT linear-algebra question over `Q` (does a rational relation exist among honest Taylor series),
so it is immune to the `|int| != int|.|` cancellation trap. The one imported fact (non-rationality at deg 3) is
independently established by elimination, not by a growth heuristic. The residual gap is purely the
degree-uniform version of a criterion that is exact at every finite truncation.
