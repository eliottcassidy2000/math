---
id: THM-854
title: Path-relative kappa equivariance and the two-path holonomy obstruction
status: RESERVED STUB - all-size identity proved; exact n=5..8 witness census in progress
source: codex-2026-07-15-S15
depends_on: [THM-781, THM-796, THM-849, THM-852]
related: [THM-549, THM-793, HYP-6905]
---

# THM-854 - path-relative kappa holonomy

Reserved for the correction to THM-852's proposed involutive-realizer
consequence.  The all-tile flip is relative to the marked Hamiltonian path
`P`; its actual equivariance law is

```text
sigma kappa_P = kappa_(sigma P) sigma,
```

not commutation with a fixed `kappa_P`.  Thus a direct realizer
`sigma T=kappa_P T` satisfies

```text
sigma^2 T=kappa_(sigma P) kappa_P T
             =Flip_(E(P) symmetric_difference E(sigma P)) T.
```

The exact finite census and the corrected scope of the conditional odd-coset
lemma are being recorded here.  What is already known is that the square
condition used in THM-852 fails for every black direct realizer through
`n=8`; what remains before promotion is a standalone replay and a precise
statement of which marked-path data the node quotient destroys.
