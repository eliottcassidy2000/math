# Independent audit: gcd-fibre transfer

**Verdict:** proof audit passes. No unrestricted LRC(14) result is claimed.
The audit code imports only Python's standard library, not the producer's
code. Tests are exact rational arithmetic and remain enabled under `-O`.

## Mathematical checks

- The `d` phases form the `h=d/gcd(d,v)` translated grid with multiplicity
  `gcd(d,v)`. Closest distance to a half-integer is
  `||h/2-u*s||/h`, yielding F1. For `h=1` this is precisely `||u*s||`.
- The nearest-integer choice of `m`, followed by `u*j=m mod h`, reaches the
  maximizing phase modulo one. Negative `m` and rounding ties are harmless.
  The chosen `0<=j<h<=d` is a native lift, including endpoint base times.
- For positive `a`, odd `h` imposes `||u*s||>=a`; even `h` instead imposes
  `||u*s||<=h(1/2-r)`. The parity reversal cannot be suppressed.
- The closed interval criterion follows directly from
  `n+lo<=uR` and `uL<=n+hi`. Its ceiling/floor test and
  `s=max(L,(n+lo)/u)` construction preserve singleton intersections.
  It costs O(k) rational operations, not O(k) bit operations.
- Every proper orbit has order at least the smallest prime factor of `d`.
  Thus the advertised radius bound gives all nonmultiple attachments;
  multiples reduce exactly to the undilated problem. Nonempty base-safe
  sets are necessary for the stated rejection-set identity.
- Safe-fibre sharpness holds also for `h=1`: `C={1},s=1/2,d=1,v=2`
  gives attached clearance zero. For odd `h` the choice `v=2`, and even `h`
  the choice `v=1`, attains the threshold at a core-safe point.
  This does not establish optimality of the global rejection-law threshold.
- Independent tail argument for each bank core: a closed safe interval of
  length `L>0` cannot fit inside a connected open bad gap of a speed `v`
  when `v*L>=2r`. Therefore `ceil(2r/L)` is a valid all-speed cutoff;
  equality is valid because endpoints are safe. The audit recomputes all
  28 maximal-width cutoffs and every finite attachment below each cutoff.
- The resulting all-dilation corollary is valid for these inherited cores.
  Dividing a dilated core by its gcd gives its existing normalized shape;
  there are no new primitive core families. A full attached row may still
  be primitive. No arbitrary-core claim follows.

## Reproduction

```
python3 session120_research/audit_lrc/check_independent.py > session120_research/audit_lrc/normal.json
python3 -O session120_research/audit_lrc/check_independent.py > session120_research/audit_lrc/optimized.json
cmp session120_research/audit_lrc/normal.json session120_research/audit_lrc/optimized.json
```

The receipts contain exact declared counts. Finite checks support the audit;
they do not replace the all-input proofs above or extend the bank's scope.
