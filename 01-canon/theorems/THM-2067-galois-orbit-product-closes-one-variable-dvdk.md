---
id: THM-2067
title: "Galois orbit-products close the one-variable constant-term nonvanishing input"
status: >
  RESERVED FOR ACTIVE PROOF AUDIT. The proposed proof replaces the monodromy
  step in THM-1605 by the transitive Galois action of the irreducible
  polynomial X^M-tR(X) over C(t). The orbit-product contradiction is proved
  below; the remaining audit is to verify every interface with THM-1550's
  small-root-product criterion and to rule out hidden separability or base-
  field assumptions before promoting this file to PROVED.
source: codex-2026-07-21-dvdk-galois-orbit
depends_on:
  - THM-1550
related:
  - THM-1605
  - THM-1630
  - THM-2022
---

# THM-2067 -- Galois orbit-products and one-variable nonvanishing

This namespace is reserved because it may remove the only cited mathematical
premise from the THM-2022 proof of NC2/GMC(2). The target is deliberately
narrow: reprove the one-variable statement actually consumed there, not the
stronger asymptotic theorem of Duistermaat--van der Kallen.

Let `R in C[X]` have `R(0) r_d != 0` and degree
`d=M+N`, where `M,N>=1`, and put

```text
Phi(X)=X^M-tR(X) in C(t)[X].
```

The proposed proof has three components.

1. `Phi` is irreducible: viewed in `C[X,t]` it has degree one in `t`; any
   factor independent of `t` would divide both `X^M` and `R(X)`, whose gcd is
   one because `R(0)!=0`. Gauss then passes irreducibility to `C(t)[X]`.
2. The following orbit-product lemma is purely algebraic. If a separable
   irreducible degree-`d` polynomial over a field has transitive Galois group
   on its roots, its total root product is a nonzero constant, and a proper
   nonempty root subset has product `c t`, then multiplying that identity over
   the subset orbit gives

   ```text
   (c t)^|O| = (product of all roots)^eta,
   ```

   because transitivity makes every root occur in the same number `eta` of
   orbit subsets. The right side is constant and the left side is not.
3. THM-1550 identifies vanishing of every positive-power constant term with
   the identity `Pi(t)=c t` for the product of the `M` roots tending to zero.
   Since `0<M<d`, the orbit-product lemma forbids that identity.

Known: the orbit count, Vieta constant
`product(Phi-roots)=(-1)^d r_0/r_d`, and `c!=0` are elementary. Still being
audited: the exact field in which THM-1550's germ identity is interpreted and
the passage from its chosen small-root subset to the splitting field over
`C(t)`. No DvdK theorem, saddle asymptotic, or numerical-semigroup
noncancellation is being assumed.
