---
id: THM-357
name: lonely-runner-endpoint-protection-trichotomy
status: PROVED
date: 2026-05-30
session: codex-2026-05-30-S360
depends_on:
  - THM-355
  - HYP-1794
  - HYP-1802
results:
  - 05-knowledge/results/lonely_runner_endpoint_protection_s360.out
---

# THM-357: Lonely Runner Endpoint-Protection Trichotomy

## Statement

Let

```text
V = {v_1, ..., v_k}
```

be a finite set of distinct positive integer speeds, let `n=k+1`, and set
`delta=1/n`.  Define the open forbidden set on the time circle:

```text
F(V) = union_{v in V} { t in R/Z : ||v t|| < delta }.
```

Every forbidden endpoint lies in the finite quotient

```text
Q(V) = n * lcm(V),
```

and has the form

```text
e = (n*m + eps)/(n*v)  mod 1,
eps in {-1,+1},  m=0,...,v-1.
```

Call such an endpoint `e` **protected** when it is strictly contained in at
least one forbidden interval, equivalently when

```text
e in F(V).
```

Then exactly one of the following alternatives holds.

1. **Positive gap.**  `measure(F(V)) < 1`.  The complement contains a
   positive-length open interval, and every point of that interval is a lonely
   witness.

2. **Boundary-only tight case.**  `measure(F(V)) = 1`, but at least one
   forbidden endpoint is unprotected.  Every unprotected endpoint is a lonely
   boundary witness.

3. **Full open cover.**  `F(V)=R/Z`.  Equivalently,

   ```text
   measure(F(V)) = 1
   and every forbidden endpoint is protected.
   ```

Consequently, a reduced Lonely Runner counterexample is exactly a
full-measure endpoint-protection certificate: the forbidden interval union has
measure `1` and every endpoint is protected.

For an endpoint

```text
e = (n*m + eps)/(n*v_i),
```

a speed `v_j` protects `e` exactly when there is an integer `a` such that

```text
| v_j*(n*m + eps) - a*n*v_i | < v_i.
```

Thus the certificate is a finite arithmetic incidence hypergraph on the
quotient `Z/Q(V)Z`.

## Proof

The set `F(V)` is a finite union of open intervals in `R/Z`.  Its boundary is
contained in the finite endpoint set displayed above, because the boundary of

```text
||v t|| < 1/n
```

is obtained from

```text
v t = m +/- 1/n  mod 1.
```

If `measure(F(V)) < 1`, the finite interval decomposition of the circle has at
least one complementary cell of positive length.  Every point in such a cell
satisfies

```text
||v_i t|| >= 1/n
```

for all speeds, so it is a lonely witness.  This is the positive-gap case.

Assume now that `measure(F(V)) = 1`.  Since the decomposition has only
finitely many endpoints, the complement of `F(V)` cannot contain any
positive-length interval; otherwise the measure would be less than `1`.
Therefore any nonempty complement is a finite set of boundary points.  A
boundary point is outside `F(V)` exactly when no forbidden interval strictly
contains it, i.e. exactly when it is unprotected.  Such a point satisfies the
closed inequalities and is a boundary lonely witness.

It follows that `F(V)` covers the whole circle if and only if the full-measure
case has no unprotected endpoint.  That is precisely the endpoint-protection
certificate.

Finally, for the integer criterion, write

```text
e = (n*m + eps)/(n*v_i).
```

The condition that `v_j` strictly protects `e` is

```text
|| v_j e || < 1/n.
```

Equivalently, for some integer `a`,

```text
| v_j*(n*m + eps)/(n*v_i) - a | < 1/n.
```

Multiplying by `n*v_i` gives exactly

```text
| v_j*(n*m + eps) - a*n*v_i | < v_i.
```

This proves the finite arithmetic form.

## Computed Audit

`04-computation/lonely_runner_endpoint_protection_s360.py` checks the
trichotomy exactly over rational arithmetic and verifies that the direct
strict-containment test agrees with the integer inequality above.

In the bounded primitive scans inherited from S357:

```text
k=3, max_speed=24
k=4, max_speed=24
k=5, max_speed=20
k=6, max_speed=16
k=7, max_speed=14
```

every full-measure case was boundary-only and no open-cover certificate was
found.  The known tight examples have unprotected boundary witnesses, with
first witness `1/(k+1)`, while the full endpoint graph lives naturally at the
larger quotient `Q(V)=n*lcm(V)`.

## Use

THM-357 separates the topology from the conjecture.  The positive-gap and
boundary-only alternatives are already witnesses.  The unresolved Lonely
Runner content is now the absence of the third alternative:

```text
no primitive speed set V admits a full-measure all-protected endpoint graph.
```

This is the LRC analogue of THM-355's quotient-gap transport viewpoint.  A
missing open cell is a positive quotient gap; an unprotected endpoint is a
boundary quotient gap; and a counterexample would have to erase both by making
every endpoint an incident protected boundary.

## Related

- HYP-1794: Lonely Runner as a quotient-gap certificate.
- HYP-1802: endpoint-protection obstruction.
- THM-355: quotient gaps are zero transport rows and columns.
- `04-computation/lonely_runner_residue_probe_s356.py`.
- `04-computation/lonely_runner_tight_scan_s357.py`.
- `04-computation/lonely_runner_endpoint_protection_s360.py`.
