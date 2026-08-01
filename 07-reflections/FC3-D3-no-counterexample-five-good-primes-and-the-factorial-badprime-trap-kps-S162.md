---
source: kind-pasteur-2026-07-24-S162 (Opus 4.8)
status: RESULT (strong finite-field evidence at D=3) + a caught trap. Pushed the FC(3) integer-leak
  ideal-emptiness test to D=3 via a direct F_q point search (evaluate all 6 leaks over the whole grid F_q^5).
  CAUGHT a bad-prime trap: q<=18 give FACTORIAL-VANISHING artifacts (fake "counterexamples": q=7->374, q=13->10)
  because T=A!P!Q!*g kills terms mod small q. Good primes q>18 are faithful. Method VALIDATED on D=2
  (proved empty over Q -> 0 zeros over every good q). D=3: 0 common F_q-zeros of all 6 leaks for
  q=19,23,29,31,37 -> strong evidence NO FC(3) counterexample at D=3 (cyclic-weight family). Not an exact proof
  (extension-field points not excluded); the exact ideal-emptiness needs an F4/F5 CAS.
tags: [factorial-conjecture, finite-fields, point-search, bad-primes, groebner, transversality]
related: [kps-S159, kps-S161]
---

# FC(3) at D=3: no counterexample over five good primes (and the factorial bad-prime trap)

## 1. The method: F_q point search for the leak ideal
The FC(3) cyclic-weight leaks `leak_j = L(f^{3j})/(3j)!` are integer polynomials in `c1..c5`. A counterexample at
D=3 is a common zero of `leak1,...,leak6` (`P=5`, need `P+1` leaks). Emptiness of that variety = no counterexample.
Direct test over a finite field: evaluate all 6 leaks at **every** point of `F_q^5` (chunked over `c0` to bound
memory) and count common zeros. `0` common zeros over a good prime = no `F_q`-rational counterexample.

## 2. The bad-prime trap (caught)
First pass over `q=7,13` reported `374` and `10` "counterexamples" -- **artifacts.** The moment table
`T(A,P,Q)=A!P!Q!*g` has factorial factors, so for `q<=18` many `T ≡ 0 (mod q)` (e.g. `A!≡0 mod 7` once `A>=7`),
which **truncates the leaks** and inflates the mod-`q` variety. `q<=18` are **bad primes**; the reduction is not
faithful. **Only `q>18` (no factorial vanishing, `A!,P!,Q! != 0` for degrees `<=18`) is valid.** (Same
detection-floor discipline as kps-S151: a signal from a degenerate reduction is not a signal.)

## 3. Validation on D=2 (proved empty over Q)
kps-S161 PROVED `V(leak1,leak2,leak3)=varnothing` at D=2 (Groebner `=(1)`). The `F_q` search must therefore return
`0` for every good prime -- and it does:
> **D=2, `q=19,23,29,31,37`: `0` common zeros each.** Method validated against the exact result.

## 4. D=3 result
> **D=3, `q=19,23,29,31,37`: `0` common `F_q`-zeros of all six leaks `leak1..leak6`, every prime.**
No `F_q`-rational FC(3) counterexample for five independent good primes. Combined with kps-S159 (generic
transversality / no family, D<=5) and kps-S161 (exact empty at D=2), this is **strong evidence that the
cyclic-weight FC(3) leak variety is empty at D=3** -- i.e. no counterexample, isolated or family.

## 5. Honest limits
- **Not a proof.** Zero `F_q`-points for good `q` is strong evidence but not exact emptiness over `Qbar`: a
  counterexample defined over a number field with a large Galois orbit could evade small primes (its reduction is
  `F_{q^k}`-rational, not `F_q`). Five good primes at `0` makes this unlikely (by Chebotarev an actual point would
  surface as an `F_q`-point for a positive density of primes), but it is not a certificate. Strengthening:
  more primes, or `F_{q^k}` point-counts, or the zeta/point-count of `{leak1..leak5}` (0-dim) then `leak6`.
- **Exact D=3** (Groebner `=(1)` over `Q`, or a Nullstellensatz certificate) remains at the tooling frontier:
  sympy Buchberger does not finish on degree-18 polys in 5 vars, python-flint exposes fast `nmod_mat` but no
  Groebner, and the Macaulay matrix at the needed degree is memory-prohibitive in dense form. A real F4/F5 engine
  (msolve / Magma / Singular `std`) would settle it; not available here.

## 6. Status
- D=2: **exact**, empty (no counterexample, isolated or family).
- D=3: **strong finite-field evidence** (0 over 5 good primes, method validated), empty; not exact.
- D<=5: generic transversality (no family).
The FC(3)-in-the-cyclic-weight-family no-counterexample claim now has exact confirmation at D=2 and validated
finite-field confirmation at D=3. Frontier: an F4/F5 run for the exact D=3 certificate.

Files: `/tmp/{fq3,disc_d2}.py`. Builds on kps-S159/S161.
