# Independent divisor/F1/defect audit

2026-10-26. **PASS** for the elementary claims in the root's
`session66-root/DIVISOR_F1_BRIDGE.md`, with explicit semantic boundaries.

The free prime monoid and its monoid ring are a legitimate carrier, not a
field structure on inverse aliquot labels. The displayed polynomial deletes
only the top monomial: prime evaluation gives `s(n)`, augmentation at one
gives `tau(n)-1`. For n=1 both give zero. Swapping prime generators 2 and3
keeps the unweighted chain but changes s(4)=3 to s(9)=4. The fiber relation
fails multiplicative congruence by s(2)=s(3), s(4)!=s(6).

Both strict divisibility inequalities follow by the stated distinct-divisor
injections, also at n=1. They imply antichain fibers, not finite fibers:
the prime fiber is infinite. No density conclusion follows.

For N>1, the interior excludes 1 and N and is complement-stable. The proper
set excludes only N and is not stable. At N210 these sizes are14 and15;
all seven listed complement pairs and s(210)=366 check. Of all pairs the
unique smallest numerical gap is the pair14,15.

For any involution with f fixed vertices and b swaps, the unordered edges
split into C(f,2) fixed-endpoint edges, b forced reversed edges, and
M=b(f+b-1) two-element edge orbits. Their polynomial factors respectively
are2,2z,2(1+z²). This proves exactly

    sum_T z^D(T)=2^[C(f,2)+b+M] z^b (1+z²)^M.

Its minimum is b. The increasing numerical orientation under divisor
complementation reverses every edge, so has the maximum C(v,2), not the
minimum. At210 the polynomial is `2^49 z^7(1+z²)^42`, minimum7 and maximum91.
This is about a specified involution, not arithmetic Frobenius or a native
orientation forced by arithmetic. Complement-invariant trace ties cannot
orient the swapped endpoints.

Primary-source scope checked: [Deitmar, Schemes over F1](https://arxiv.org/html/math/0404185v7),
section1/Theorem1.1, defines F1 as the trivial monoid, monoid-ring base
extension, and its adjunction to forgetting ring addition. This supports
the stated compatible construction and no additional F1 fiber structure.

## Finite exact verification

`audit_divisor.py` independently enumerates the **entire defect distribution**
for every fixed-point/swap type on at most6 vertices (including f>1), and
matches every coefficient of the polynomial. It also checks divisor
polynomial evaluation for N1..600, divisibility inequalities for n1..100
and k2..30, and the root's explicit minimizing210 orientation edge by edge.
The frozen root checker is replayed with normal and optimized Python.

```
PYTHONDONTWRITEBYTECODE=1 python3 session66-audit-aliquot/audit_divisor.py > session66-audit-aliquot/divisor-normal.json
PYTHONDONTWRITEBYTECODE=1 python3 -O session66-audit-aliquot/audit_divisor.py > session66-audit-aliquot/divisor-optimized.json
cmp session66-audit-aliquot/divisor-normal.json session66-audit-aliquot/divisor-optimized.json
```

Both independent runs pass with byte-identical output. Suggested wording
repairs (use “evaluation” rather than prime-weighted “augmentation”; show
s(4)!=s(6) explicitly) were applied by the root. No mathematical failure
or unresolved blocker remains. No producer source was edited by this auditor.
