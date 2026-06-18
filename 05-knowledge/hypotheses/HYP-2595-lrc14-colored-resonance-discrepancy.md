# HYP-2595: LRC(14) Colored Resonance Discrepancy Bound

**Status:** OPEN proof program with a rigorous Fourier/color-resonance identity
and exact stress evidence.

HYP-2594 reduced finite placement to bounding

`Delta(P,E,V) = V*Sigma(P,E) - actual_count(P,E,V)`.

The raw interval bound uses the total micro-component count `K`, but this is
far too pessimistic.  The sharper target is:

`Delta(P,E,V) <= 8*(k + c_GP) + 1`,

where `k=|E|` and `c_GP` is the number of components of the small-speed safe
set `G_P` (with `c_GP=1` when `P` is empty).

## Fourier/Color-Resonance Identity

Let `h(x)=1_{||x||>=1/14}` with boundary values taken as `1/2` for the
Fourier identity.  For a fixed color `b`, define

`F_b(x)=prod_{p in P} h(px) * prod_{e in E} h(ex-b/14)`.

The half-boundary colored count is

`N(P,E,V)=sum_b sum_{t=0}^{V-1} F_b(t/V+b/(14V))`.

Then

`N - V*Sigma = V * sum_{ell != 0} sum_b Fhat_b(ell V) exp(2*pi*i*ell*b/14)`.

Expanding each `h` factor, a Fourier monomial can survive the color sum only
if

`sum_p p*n_p + sum_e e*n_e = ell*V`

and

`sum_e n_e == ell (mod 14)`.

This is the colored resonance condition.  It is the first rigorous reason raw
endpoint count `K` should overestimate discrepancy: most endpoint events do
not satisfy the compatible `V`-resonance plus color congruence.

The closed CRT count differs from the half-boundary count only on explicit
binding-boundary residues, which can be counted separately.

## Evidence

Script: `04-computation/lrc14_colored_resonance_deficit_codex.py`.
Output: `05-knowledge/results/lrc14_colored_resonance_deficit_codex.out`.

Named/adversarial rows:

- The HYP-2594 stored max-deficit row has `K=4234`, but actual deficit is
  only `71.867`.
- A `P=empty`, `k=13`, `Emax=660` row has `K=11536`, open deficit `128.347`,
  closed actual deficit `94.347`, and boundary bonus `34`.
- A boundary-bonus row has open deficit `75.457` but closed actual deficit
  `-44.543`, because `120` binding-boundary grid points are recovered.
- A deterministic constant-100 failure
  `P=(1,2,11)`, `E=(0,84,293,301,355,416,485,665,843,886)`, `V=1203`
  has actual deficit `135.435`, so a universal constant `100` is false,
  but `8*(k+c_GP)+1 = 177` clears it.

Reproducible stress bank:

- `149` total rows (`9` deterministic plus `140` random covering rows).
- `0` zero actual witnesses.
- Bound `100`: `1` violation.
- Bound `7*(k+c_GP)+1`: `0` violations, max pressure `0.989069`.
- Bound `8*(k+c_GP)+1`: `0` violations, max pressure `0.866530`.
- Bound `8*k*c_GP+1`: `0` violations but much larger median size.
- Raw `K`: `0` violations but median size `862`, much looser.

The old covering hard-core tower `{1,...,11,13} union {84m}` is tame in this
colored lens: for `m=1`, actual deficit is only `492/385`; for `m=5`, the
actual deficit is already negative.  Thus the colored discrepancy spike is not
coming from that principal tower.

The candidate-bound tournament ranks raw `K` most robust but too large; among
small useful candidates, `8*(k+c_GP)+1` beats `8*k*c_GP+1` and the constant
bound.

## Why This Would Nearly Close The Placement Layer

Across all `P subset {1,...,13}`, the maximum `c_GP` is `32`, attained at
`P=(10,11,12,13)`.  Since `k<=13`, the candidate bound gives

`Delta <= 8*(13+32)+1 = 361`.

Combining this with the HYP-2593 structured-bank `Sigma` floor

`Sigma >= 14249/28028`

would certify every large-`V` row once

`V > 361 / (14249/28028)`, i.e. `V >= 711`.

Thus the proof would reduce the colored CRT placement layer to a finite exact
check for `V<711`, assuming the uniform `Sigma` floor is also promoted from
evidence to theorem.

## Proof Route

1. Formalize the half-boundary Fourier identity above for step functions via
   Fejer approximation or direct finite Fourier-Stieltjes endpoint accounting.
2. Bound the surviving color-compatible resonances by `O(k+c_GP)`.
3. Bound the closed-boundary correction by the same order, or show it only
   improves the lower bound for actual witnesses.
4. Combine with HYP-2593/HYP-2594 and exact finite checking below the resulting
   cutoff.

This is not a proof of LRC(14), but it turns the remaining colored arithmetic
discrepancy question into a precise resonance lemma.

See also `HYP-2594`, `HYP-2593`, THM-527/528, and `OPEN-Q-108`.
