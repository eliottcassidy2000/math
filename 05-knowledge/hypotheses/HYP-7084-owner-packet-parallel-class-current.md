# HYP-7084 — the owner packet as a parallel-class current

**Status:** CLAIMED / IN PROGRESS (codex-2026-07-16-S19).

This note reserves the exact class-circle reformulation of `HYP-7083`.  On the
one-miss set of `F=E union {t}`, write

```text
q=floor(14{tx})=2c+eta,   eta in {0,1},   u=2c mod 7.
```

Then `sec(tx)=c` while `sec(2tx)=u+eta`.  The two possible target classes are
therefore adjacent vertices `u,u+1` of the `THM-913` parallel-class circle:
the even half-sector is the degenerate pair `(c,c)` of sum class `u`, and the
odd half-sector is the boundary pair `(c,c+1)` of sum class `u+1`.

Let `mu_s^eta(u)` be the exact mass on which `F` has unique miss `s` and the
owner half-sector has class coordinate `(u,eta)`.  Put

```text
M_s(u)=mu_s^0(u)+mu_s^1(u),
J_s(u)=mu_s^0(u)-mu_s^1(u).
```

The claimed exact decomposition is

```text
Delta_F(2t)
 = 1/2 sum_s (M_s(s)+M_s(s-1)) - p1(F)/7
   + 1/2 sum_s (J_s(s)-J_s(s-1)).
```

Thus the parity remainder is the same-label trace of the discrete divergence
of an oriented current on `C_7`.  This quotient preserves the exact finite
owner correction and its miss-label incidence.  It destroys slow-wall order,
the additive relation lattice of `E`, and most crossing information in the
book drawing; the `THM-913` circulant crossing energy is therefore only a
candidate norm until a comparison inequality is proved.

Reserved verifier/output:
`04-computation/lrc14_owner_packet_parallel_class_current_codex_S19.py` and
`05-knowledge/results/lrc14_owner_packet_parallel_class_current_codex_S19.out`.
The immediate tasks are to verify the decomposition fraction-exactly on the
`HYP-7083` bank and test whether any circulant energy controls the diagonal
divergence with a constant strong enough for `0.097`.
