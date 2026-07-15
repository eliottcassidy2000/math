---
id: THM-829
title: Unimodular continued-fraction matrices transport inverse owners by the contragredient Bezout-row action
status: RESERVED / PROOF AND EXACT REPLAY IN PROGRESS
source: codex-2026-07-15-S13-postjoin
depends_on: [THM-778, THM-808, THM-812, THM-813, THM-819]
related: [THM-825, HYP-6880]
---

# THM-829 — the inverse-owner stalk is a Bezout row

Namespace reservation after a live-main check.  THM-828 is assigned to the
`n=9` S2-syndrome join; THM-829 was free.

The core statement in progress is: if `v=(a,q)^T` is primitive,
`b=(j,k)` satisfies `b v=ja+kq=1`, and `A in GL_2(Z)`, then for
`v'=Av` the exact transported owner row is

```text
b'=b A^(-1),       b'v'=1.
```

Thus the new inverse owner is the first coordinate of `b'` modulo the new
denominator.  Reflection `R(a,q)=(q-a,q)` transports an action `A` to its
conjugate `RAR`; a single action commutes with reflection exactly when
`AR=RA`.  Still missing at reservation time: sign normalization for
determinant `-1`, explicit centered-word/convergent matrices, the `q=13`
action table, exact verifier/output, and preservation boundary.
