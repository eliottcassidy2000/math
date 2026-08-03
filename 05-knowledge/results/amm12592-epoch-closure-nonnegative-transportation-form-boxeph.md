# AMM 12592: the epoch-closure problem in nonnegative transportation form

**Status: VERIFIED-EXACT (identity + R=8 witness replay). boxeph, 2026-08-03.**
Companion observation to THM-3002 (*) / THM-3026; no new closure is claimed.

## Statement

In THM-3002's representation problem (*) for the epoch `[R, 2R)`, substitute
`delta_{i,k} = binom(d_i,k) - 2 f_{i,k}` per block. The Lucas box
(`|delta| <= binom`, `delta == binom mod 2`) maps bijectively onto

```text
f_{i,k} integer,   0 <= f_{i,k} <= binom(d_i, k),
```

and, because `sum_k binom(d,k) B_{d,k} = (x + (1-x))^d = 1`, the identity (*)
`q^{R-1} = sum_i p^i Delta_i` becomes

```text
sum_{i<R} p^i F_i(p)  =  T_R(p) := [ (1-p^R)/(1-p) - (1-p)^{R-1} ] / 2,   (**)
F_i = sum_k f_{i,k} B_{d_i,k},
```

with `F_i` exactly the [0,1]-valued decided-tree polynomials of THM-2966 —
the extractor's probability weights, restored at epoch scale. Verified: the
substitution round-trips the reconstructed R=8 floor witness exactly and
(**) holds as a polynomial identity in `Z[p]`.

## Why this form is useful

1. **Parity clock made scalar.** `T_R in Z[p]` iff every `binom(R-1,j)` is
   odd iff `R` is a power of two (Lucas) — THM-3002 4b's depth-free mod-2
   clock is literally the integrality of the target vector.
2. **The LP relaxation is the Hall/ARCH analysis.** (**) asks for an integer
   point in a capacitated transportation-style polytope (variables
   `f_{i,k} >= 0` with upper bounds, one linear equation per coefficient of
   `T_R`). THM-3009's (ARCH) cuts and THM-3024's Hall cuts are exactly
   feasibility conditions for its real relaxation. Hence the entire
   construction-vs-floor gap ("does the golden floor profile close at every
   R") is an **integrality-gap** question for this polytope family.
3. **Where total unimodularity fails.** The Bernstein-to-power conversion
   couples each `f_{i,k}` into many coefficient rows with signs
   `(-1)^j binom(k,j)`; equivalently, in the deficit-flow picture the
   elementary move multiplies by `(p+q)` and copies mass to BOTH lattice
   children with the same coefficient. This coupled copy is what breaks the
   network-matrix structure, so integer feasibility does NOT follow from LP
   feasibility with margin — the redistribution obstruction of THM-3026
   (row overshoots 49/276/1541 at R=8/16/32) is the visible face of that
   integrality gap.
4. **Sharp target for the all-R theorem.** A proof that the (**) polytope
   with profile `d_i = floor(gamma(R+i))`, `gamma > gamma* = log_5(phi^2)`,
   contains an integer point for all large `R` gives `C* <= 1 + gamma` for
   every such gamma, hence `C* <= 1 + gamma*` by infimum — which with
   THM-3024's floor pins `C* = log_5(5 phi^2)`. The floor profile itself
   never needs to close.

## Files

Verification inline in this session's close-out; substitution replay against
`04-computation/amm12592_floor_witnesses_R8_R16_R32.json` (R=8 entry).
