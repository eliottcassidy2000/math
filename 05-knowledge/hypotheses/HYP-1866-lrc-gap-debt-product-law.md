---
id: HYP-1866
status: OPEN
source: codex-2026-05-31-S398
related:
  - THM-357
  - THM-365
  - THM-366
  - THM-367
  - HYP-1842
  - HYP-1844
  - HYP-1857
  - HYP-1858
  - HYP-1859
  - HYP-1861
  - HYP-1867
  - HYP-1868
---

# HYP-1866: LRC gap-debt product forbids simultaneous zero

## Statement

For a Lonely Runner speed set `V` at denominator `n=|V|+1`, define the
Archimedean gap size

```text
Gap_A(V) = max_gap(V) / (1/n).
```

Define `Debt_2(V)` as the exposed endpoint debt, first by the raw count of
unprotected forbidden endpoints and then, in the refined version, by their
dyadic or denominator-quotient layers.

The hypothesis is that every genuine repair branch, after the small-denominator
boundary skeleton has been closed, satisfies a positive product bound

```text
Gap_A(V) * Debt_2(V) >= c(n) > 0.
```

This is not a universal identity over all speed sets.  Boundary-only tight
systems have

```text
Gap_A = 0,  Debt_2 > 0,
```

and are already certified by endpoint debt.  The product law belongs to the
positive-gap debt-export branch: when a repair shrinks the visible real gap, it
pushes debt into a denominator layer.  A counterexample to LRC must instead
force both coordinates to zero:

```text
Gap_A = 0,  Debt_2 = 0.
```

So the target proof is an adelic obstruction: the real absolute value and the
2-adic/denominator absolute value cannot both vanish on the same endpoint core.

## Evidence

The S398 audit computes exact product rows for the known debt-export ladders.

```text
n=14, d=7:   Gap_A=5/924,  Debt=84,   product=5/11
n=14, d=14:  Gap_A=5/1848, Debt=168,  product=5/11

n=16, d=2:   Gap_A=1/33,   Debt=34,   product=34/33
n=16, d=4:   Gap_A=1/66,   Debt=68,   product=34/33
n=16, d=8:   Gap_A=1/132,  Debt=140,  product=35/33
n=16, d=16:  Gap_A=1/264,  Debt=280,  product=35/33

n=18, d=9:   Gap_A=1/176,  Debt=176,  product=1
n=18, d=18:  Gap_A=1/352,  Debt=352,  product=1
```

The `n=16` row has a meaningful phase step:

```text
(35/33) / (34/33) = 35/34.
```

So the conserved quantity is not simply "raw count times gap" across every
layer boundary.  The better reading is that raw endpoint count is the first
2-adic size proxy.  A final invariant probably needs a layer weight that
charges the transition from the `34/33` plateau to the `35/33` plateau.

The boundary branch confirms the logical split:

```text
n=14 initial: Gap_A=0, Debt=6
n=16 initial: Gap_A=0, Debt=8
n=18 initial: Gap_A=0, Debt=6
```

Those are not failures of the product route; they are immediate witnesses.

## Proof Route

The proof should have three regimes.

```text
positive gap       => lonely interval exists;
boundary debt      => lonely endpoint exists;
gap-debt product   => no repair branch can converge to (0,0).
```

The missing formal object is a weighted endpoint-debt norm.  Candidate weights:

1. raw unprotected endpoint count;
2. first exposed quotient layer `den(endpoint)/gcd(den(endpoint),n)`;
3. dyadic depth for pure `n=16` systems;
4. labelled-cycle slack from THM-365;
5. local capacity deficit from THM-367.

The desired theorem would say that every sieve-complete, all-endpoint-repaired
attempt either has a positive visible gap, has exposed endpoint debt, or carries
positive weighted product.

At `n=16`, HYP-1859's global dyadic debt-flow inequality is the most concrete
test case.  The product law says that closing the `16`-gate branch cannot simply
delete the old odd-unit witnesses.  It must export their mass to descendants,
and the exported mass has positive norm.

## Cayley-Dickson Reading

The Cayley-Dickson analogy becomes literal here.  Doubling does not preserve
the same visible object; it preserves a norm-like obstruction after changing
coordinates.  The `16`-gate kills the old unit endpoints the way a zero-divisor
move kills one component, but the obstruction reappears as endpoint debt in a
lower quotient layer.

The `n=16` phase tax `35/34` is therefore useful.  It may be the first place
where the raw endpoint count is not the true norm and the weighted dyadic debt
must be corrected.

HYP-1868 gives the proposed correction geometry: exposed endpoints form a
finite frontier in the Bruhat-Tits tree of `PGL_2(Q_p)`.  On the pure `n=16`
dyadic branch, `Gap_A*2^h` stays constant and the `35/34` tax is exactly a
normalized frontier-mass jump from `17` to `35/2`.

## Sources

- `04-computation/lrc_gap_debt_product_s398.py`
- `05-knowledge/results/lrc_gap_debt_product_s398.out`
- `07-reflections/lrc-gap-debt-product-s398.md`
- HYP-1859.
- THM-357.
- THM-365.
- THM-366.
- THM-367.
- HYP-1867.
- HYP-1868.
