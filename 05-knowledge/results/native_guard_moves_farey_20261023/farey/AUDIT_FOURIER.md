# Independent audit of the exact ternary mod-eight observer

23 October 2026. **VERIFIED:** the stated algebraic proof and its native
scope. **FINITE-EXACT:** the final normal and optimized receipt replays
both pass 17,225 checks. Both maintenance nits below were repaired and the
receipt now matches its frozen guard input.

Audited files: [`../fourier/THEOREM.md`](../fourier/THEOREM.md) and
[`../fourier/check.py`](../fourier/check.py). This is a proof audit, separate
from rerunning the producer.

## Proof findings

1. **Coefficient ring and inversion.** In the quotient by
   `w^2+w+1`, the product `(w-1)(w^2-1)=3` holds exactly. Since 3 is a unit
   modulo eight, both factors are units in this commutative ring. Thus
   even a nonzero zero-divisor Fourier coefficient cannot be annihilated
   by a nontrivial root difference. Orthogonality gives the displayed
   inverse because `3^r` is a unit. Independently enumerating multiplication
   by each root difference on all 64 ring elements gives a permutation.

2. **Support and the full translation stabilizer.** Translating the input
   by `v` multiplies Fourier coefficient `k` by `w^(k dot v)`. Fourier
   uniqueness and the preceding unit assertion imply that equality of the
   whole translated record is equivalent to `k dot v=0` for every nonzero
   supported coefficient. Hence `K=(span S)^perp`, including tuple-valued
   records by union of supports. No field assumption on the coefficient
   ring is being smuggled into this step. Linear algebra is over F3 only
   for the frequency and translation coordinates.

3. **Minimal observation.** If `f=g(Ca)`, every vector in `ker C` is a
   universal period, so `ker C` is contained in `K`. This proves the lower
   bound on rank; choosing a basis of `span S` attains it. Surjectivity of
   the resulting linear map gives exactly `3^(r-rank B)` lifts per
   observation. Equal output values can unite several such cosets, so the
   report correctly declines to make them all universal moves.

4. **Sharp move budget.** The low-weight translations are symmetric and
   their span is exactly their reachable subgroup: coefficients two can
   be executed by repeating a permitted addition. Connected components
   are the cosets of that span. A spanning set of weight at most `b` has
   a basis of the same weight bound, proving the alternative minimax-basis
   formula. This is a cost in parameter coordinates only. Actual vertex
   count, duration and electrical energy need separate sidecars.

5. **Chart changes, pins, and residual states.** For the convention
   `f_new(a)=f_old(Ma+c)`, the support transforms by `M^T`; phases from
   the offset are units. Rank is invariant but coordinate footprint need
   not be. Pin substitution and complete residual-table comparison are
   exact finite operations. The report correctly limits residual-state
   minimality to the fixed variable order and makes no efficient general
   compiler claim.

6. **Modulus-nine hostile.** The nonzero constant-three function on F3
   has identically zero unnormalized transform in the corresponding
   modulus-nine ring. This directly defeats this inversion mechanism at
   nine; it does not assert that all possible encodings modulo nine fail.

## Native scope and reproduction

The original order-eight lane reads the frozen owner map and harmonic
potential, reconstructs normalized masks from each actual path and source
matrix, and checks both endpoint classes and the phase. Converse reverses
the observer path as well as every arc; the merged endpoint classes and
phase remain unchanged. This is correctly distinguished from swapping the
two ends of an oriented complement line, which negates the phase.

The two-phase coordinate fixture has record values `(1,1)`, `(1,5)`,
`(0,0)` for relative phases zero, one, two. Its diagonal kernel and coupled
cost two are genuine conclusions for that fixture. They do not identify
all Hamiltonian paths with one automorphism orbit or imply that all equal
phase values have a universal kernel move. The separate validity component
prevents invalid rows from masquerading as valid zero phases.

Two maintenance nits were reported and then verified as repaired: the
electrical dependency link now names `guards/README.md`, and the checker
explicitly bounds each target ID below `8**4` and verifies its exact
reconstruction from four base-eight digits. The actual merged IDs are
below 3,528 and therefore fit; neither nit invalidated this fixture.

Commands:

```text
python -B 05-knowledge/results/native_guard_moves_farey_20261023/fourier/check.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/fourier/check.py
```

Both commands pass 17,225 checks, including the frozen-receipt comparison.
For provenance, an independent recursive comparison had isolated the sole
initial mismatch to `native/guard_receipt_sha256` while that input was
being updated. The final receipt includes the refreshed checksum and the
new digit-bound controls. No theorem correction is requested.
