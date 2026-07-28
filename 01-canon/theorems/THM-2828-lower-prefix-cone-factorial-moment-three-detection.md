---
id: THM-2828
title: "Lower-prefix cone factorial moment-three detection"
status: RESERVED / UNPROVED EMPTY STUB
source: root/audit-2809-2026-07-28
depends_on: []
related:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - HYP-8765-gmc2-radial-channel-return-tower
---

# THM-2828 -- lower-prefix cone factorial moment-three detection

**RESERVED / UNPROVED EMPTY STUB.**

Proposed scope: extend THM-2824 from one lower interval
`f_b-f_a` to every nonzero prefix-cone direction

```text
U=sum_(i=0)^(b-1) lambda_i(f_(i+1)-f_i),   lambda_i>=0,
V=f_c-f_b,                                  b<c.
```

The intended mechanism is:

1. THM-2824 makes the mixed atomic determinant linear and nonnegative in
   every `lambda_i`;
2. an exact positive formula for
   `L((f_(i+1)-f_i)(f_(j+1)-f_j)(f_(k+1)-f_k))`
   gives strict cubic orientation for every nonzero `U`; and
3. the same binary quadratic-divisibility identity then forces
   moment-three detection on the full plane `span_C{U,V}`.

In coefficient coordinates this is the cone of mean-zero lower
polynomials whose proper prefix sums are nonpositive.  It allows
arbitrarily many occupied lower factorial slots, but it does not claim
arbitrary radial coefficients, a second arbitrary cone direction, or
general GMC2.

Until the cubic tensor formula, cone theorem, exact companion, and
independent hostile audit are installed, this reservation is not a proved
result or proved dependency.
