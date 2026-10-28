# Twenty as the unique five-divisor / unit-pair overlap

2026-10-28. **PROVED:** elementary classifications below. **FINITE-EXACT:**
two independent divisor-sum algorithms agree at every source through1000001.
**CITED / INHERITED:** near-perfect terminology and the1001/343 uniqueness
mechanism. No novelty claim for classical divisor formulas or square repunits.

The user's refined idea is an overlap of representations. For20,

    s(20)=1+2+4+5+10=22=1+1+20,
    20*1=20,      20+1=21=s(20)-1.

The literal proper-divisor list has five members; the compressed expression
has three summands and uses the unit pair(1,20). Retain which representation
is meant. The following is a precise version of the proposed cubic rule:

    C(n): there exist positive integers a,b with ab=n and s(n)=1+a+b.

**Theorem P1.** For every positive integer n,

    C(n) iff tau(n)=4 or s(n)=n+2.

The first branch consists exactly of prime cubes and distinct-semiprime
products; the second uses(a,b)=(1,n). These branches are disjoint.

**Theorem P2.** **20 is the unique integer with five proper divisors that
also satisfies C(n)**, equivalently the unique such integer with s(n)=n+2.
This proves the user's proposed extra rarity in a specified all-size sense.

The formulas also distinguish nearby examples:

| n | proper-divisor count | s(n) | factor-pair compression |
|---:|---:|---:|---|
| 8 | 3 | 7 | (2,4), a prime-cube case |
| 6 | 3 | 6 | (2,3) |
| 16 | 4 | 15 | none |
| 12 | 5 | 16 | none |
| 18 | 5 | 21 | none |
| **20** | **5** | **22** | **(1,20)** |
| 28 | 5 | 28 | none |
| 24 | 7 | 36 | none |
| 104 | 7 | 106 | (1,104) |
| 464 | 9 | 466 | (1,464) |
| 650 | 11 | 652 | (1,650) |

Thus16 in the user's additional list is a genuine boundary: it is a prime
power, with s(16)=1+2+4+8, but does not satisfy C. The whole aliquot row15
contains both16 and33; only33 has the literal three-divisor form. Likewise
row22 contains exactly20 and38, with different compression types.

[THEORY.md](THEORY.md) proves P1–P4, explains the affine row through12,20,28,
connects the overlap to near-perfect numbers, and recovers the two different
square-sum mechanisms through1001. It also explains how the earlier
Liouville/Möbius work applies to the exponent-pattern change
`7*11*13 -> 7^3 -> 3*19`.

The established1001 statement is correct with its exact restriction:
among squarefree products of three primes, only1001 feeds a prime cube whose
own aliquot edge has square total. The native path is1001->343->57; its
second edge has total400. Separately1001+1001000=1001^2 is the universal
construction m+m(m-1)=m^2 at m1001, not an aliquot edge from1001.
The short orbit of1000001 passes through21->11->1->0, as independently
recomputed. Full unbounded aliquot dynamics are not inferred from these facts.

Reproduce from the repository root:

```text
python -B 05-knowledge/results/aliquot_divisor_overlap_20261028/verify.py --replay
```

The producer uses additive divisor sieving; the audit uses least-prime-factor
multiplicative recurrences and imports none of its functions. Hashes compare
all1000002 entries of both sigma and tau arrays, including index0. Normal and
optimized receipts agree. Explicit divisor lists and factor-pair searches
through5000 provide a third bounded check. This is a same-agent proof review
with independent implementations, not an independent-agent audit.
