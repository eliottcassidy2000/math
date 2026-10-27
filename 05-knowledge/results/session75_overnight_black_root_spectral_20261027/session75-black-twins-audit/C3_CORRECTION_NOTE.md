# Source-bound correction for HYP-9421 (1) and (2)

**Status: REFUTATION of the specified clauses, not a new classification.**
This separate note does not change the frozen Phase2 theorem or audit files.

## Exact incoming source

Commit `87d61f67673a6b10cbafff1bf9e1377f3c41e56c`, file
`05-knowledge/hypotheses/HYP-9421-black-descent-exceptions-tc-spectrum-and-the-native-forcing-tail.md`,
is copied byte-for-byte as `inputs/HYP-9421-at-87d61.md`. Its exact repository
path, commit and SHA256 appear in this packet's `SOURCE_BINDINGS.json`.

Clause (1) claims that the only strong tournaments with no black minimum-cut
exit word are C3, R5, and the doubled family E(M;Z). It expressly states that
E has only even orders and predicts no exception at order11.
Clause (2) combines a universal upper bound lambda+1 with the assertion that
this extra move is needed exactly for members of E(M;Z). These must be tracked
as two logically separate assertions.

## Already-audited counterexample

The order11 odd-middle family member has named rows

```
(2032,2025,2019,2021,1006,448,640,768,576,1056,496).
```

It is strong and black-active, with lambda=1 and delta=2. Its complete minimum
outgoing-shore family is masks992 and2016: respectively {5,6,7,8,9} and
{5,6,7,8,9,10}. Each complete decreasing cone consists of the same nine BLUE
Hamilton words. Exhaustion of all1,215 Hamilton words finds1,116 black words,
and every black one-step neighbor remains strong. The literal two-black route

```
P=(2,1,0,5,6,7,8,9,10,4,3)
Q=(9,8,7,6,5,0,1,3,4,2,10)
```

reaches a nonstrong target, with a strong intermediate. Thus the exact black
nonstrong distance is2=lambda+1. This is not a pure-blue source or a failure
restricted to one arbitrarily chosen minimum shore. It cannot be in the
incoming E(M;Z) family because its order is odd.

The source-bound excerpt `inputs/odd11-audited-head.json` contains the exact
previously accepted source, complete-head counts and named R/K endpoints.
It is extracted from the independent Phase2 receipt, whose SHA256 and frozen
audit manifest are recorded in that excerpt. No new search range is used here.
The accepted structural odd-family proof also gives such examples at every
odd order >=11; the order11 member alone is sufficient for both corrections.

## Required status split

- **HYP-9421 (1): REFUTED as stated.** The E-only exception classification and
  explicit prediction of no order11 exceptions fail. A replacement complete
  exception classification remains open; merely adding this family is not a
  proof of completeness.
- **HYP-9421 (2), universal upper bound d_black <= lambda+1: OPEN / unchanged.**
  The witness attains the proposed bound and does not refute it.
- **HYP-9421 (2), extra move needed exactly for E(M;Z): REFUTED.** This odd
  example also needs the extra move and is not in that even-order family.
- Clauses (3), (4), (5), the incoming full order10 census, and all other claims
  are outside this correction's certification scope.

The separate same-family core-distance2 result does not settle global
black-core distance=delta. That global question remains open.
