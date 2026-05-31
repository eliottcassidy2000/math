---
id: HYP-1865
status: OPEN
source: codex-2026-05-31-S397
related:
  - HYP-1833
  - HYP-1858
  - HYP-1859
  - HYP-1860
  - HYP-1861
  - HYP-1862
  - HYP-1863
  - HYP-1864
---

# HYP-1865: n=16 endpoint debt has an Egyptian split ledger

## Statement

At `n=16`, any primitive Lonely Runner open-cover counterexample with a
`16`-gate should induce an Egyptian split ledger across dyadic endpoint layers.
The local endpoint-protection inequalities should be compressible to a finite
bookkeeping system whose atoms are:

```text
16/N = 1/b + 1/c             Egyptian split or unsplit residue
|p*(16m+eps)-16*a*u| < u     LRC endpoint protection residue
q0(x)=x-1/x                  anti-reciprocal counterterm
```

The conjectural proof target is:

```text
Unsplit Egyptian residues must be paid by even breaker speeds,
and those breakers force either positive LRC gap or a private endpoint leaf.
```

## Evidence

`lrc_egyptian_q_bridge_s397.py` builds a coordinate ledger for residues
`1..15`.  It compares:

- two-term Egyptian splitability of `16/c`;
- product-sum target arities;
- exact endpoint count protected by residue `c` for owner `u=16`;
- exact one-coordinate half-turn missed-cell count in the normalized
  `n=16` scalar quotient;
- membership in the exact nine-residue lower cover of a `16`-gate.

The standout row is `c=15`:

```text
16/15 = 1 + 1/15,
c = -1 mod 16,
half-turn missed cells = 128.
```

This is the best one-coordinate half-turn defect from S393.  It is the
anti-residue analogue of the S394 counterterm `q0=x-1/x`.

The best support-2 half-turn pair is:

```text
(10,15), missed=160.
```

Here `15` is splitable and anti-reciprocal, while `10` is an unsplit even
breaker residue.  S397 reads this as a finite `q`-pattern: one coordinate
carries unsplit debt, the other supplies the reciprocal counterterm.

S397 also audits the exact `16`-gate lower cover

```text
(1,3,5,7,8,9,11,13,15)
```

plus five of the six even breakers

```text
(2,4,6,10,12,14).
```

Every simple five-breaker choice remains positive-gap.  The closest among
those six choices omits `14` and still has `gap/th=0.060606` with `12`
unprotected endpoints.  High `16`-multiple probes tied to the `(10,15)` pair
also remain positive-gap with many exposed endpoints.

## Interpretation

Egyptian fractions answer:

```text
Can reciprocal mass k/N split cleanly into two unit fractions?
```

The master criterion says this is equivalent to:

```text
exists d | N^2 such that d == -N mod k.
```

LRC endpoint protection asks the approximate/congruence version:

```text
Can endpoint (16m+eps)/(16u) be strictly swallowed by protector p?
```

equivalently:

```text
|p*(16m+eps)-16*a*u| < u.
```

Both are residue-splitting problems.  HYP-1865 proposes that a complete `n=16`
proof can turn endpoint protection cycles into an Egyptian-style split ledger.
The missing mass in the ledger is a `q`-counterterm: it cannot disappear; it
must reappear as positive gap or as a peelable private endpoint.

## Predictions

1. Any labelled endpoint cycle from THM-365/THM-367 should have a nonnegative
   Egyptian defect, with equality only for imprimitive fan systems.
2. Primitive gcd-breaker speeds should contribute positive defect in the
   Egyptian ledger.
3. The pair `(10,15)` should remain the canonical two-coordinate obstruction:
   `10` as unsplit even debt and `15` as anti-reciprocal counterterm.
4. A future branch-and-bound proof should prioritize unsplit residues
   `{2,3,4,5,6,7,9,10,11,13}` differently from split residues `{8,12,14,15}`,
   not merely by dyadic valuation.

## Sources

- `04-computation/lrc_egyptian_q_bridge_s397.py`
- `05-knowledge/results/lrc_egyptian_q_bridge_s397.out`
- `07-reflections/lrc-egyptian-q-bridge-s397.md`
