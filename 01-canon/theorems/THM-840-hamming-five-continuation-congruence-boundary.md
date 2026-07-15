---
id: THM-840
title: Exact endpoint geometry is insertion-Markov while the Hamming-five handoff quotient is not
status: RESERVED — abstract congruence lemma proved; the explicit H0/H1 insertion counterexample is being independently replayed
source: codex-2026-07-15-S10 continuation
depends_on: [THM-822, THM-828, THM-832]
related: [THM-837, HYP-6820]
---

# THM-840 — the Hamming-five continuation-congruence boundary

This namespace is reserved for the distinction between static reconstruction
and deterministic continuation.

For an observation `O:X->Y` and a family of operations `T_a:X->X`, an update
map `F_a` satisfying

```text
O(T_a x)=F_a(O(x))
```

exists if and only if

```text
ker(O) subset ker(O after T_a).                           (1)
```

For pure LRC comb insertion, the ordered exact open-endpoint word is such a
state because

```text
E(S union {u})=E(S) intersect {t:||ut||>1/13}.            (2)
```

The proposed sharp counterexample says that THM-822's coarse `H0=H1`
observation is not a congruence even for one named insertion.  The rows with
replacement speeds `(14,15,16,17,18)` and `(27,28,29,30,31)` share `H0=H1`.
After inserting labelled replacement `(q,u)=(6,19)`, the first row gains
handoffs `(6->2,k=3),(6->3,k=2)`, while the second gains
`(6->3,k=2),(5->6,k=3)`.

The abstract lemma is immediate, but this file remains reserved until the
handoff convention and the displayed edge sets have an independent exact
replay.  Static three-face injectivity, trivial Cech cohomology, and a
deterministic update quotient are logically different properties.
