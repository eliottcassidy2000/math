---
id: THM-1092
title: Fejer-regularized k-fold resonance-lattice identity
status: CLAIMED / ID RESERVED — the Fejer/L1 proof is drafted; no higher-cumulant estimate or LRC14 closure is claimed
source: codex-2026-07-18-S67 global-bridge audit
depends_on: [THM-965, THM-1075]
related: [THM-1070, THM-1080, THM-1085]
---

# THM-1092 — reserved rigorous k-fold identity

The target replaces THM-1075's conditionally truncated lattice series by a
canonical Fejer-regularized limit.  Convolution of the danger indicator with
the Fejer kernel gives bounded trigonometric polynomials converging in `L1`;
finite product expansion and character orthogonality then give a finite
resonance sum at each level, and the products converge in `L1` by telescoping.

This would prove the identity in Cesaro form only.  Bounding the full-support
terms, rather than merely defining them rigorously, remains open.
