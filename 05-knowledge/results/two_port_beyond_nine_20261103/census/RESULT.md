# Every labelled tournament through ten is reconstructed by two ports

**FINITE-EXACT.** The complete named sidecar (H,S,E,C) determines every
labelled tournament of order10. Together with the preceding complete census
through9, this excludes every two-port collision through10. A counterexample
to universal two-port reconstruction must have at least11 vertices.

This is not an all-order theorem. The smaller sidecar is not asserted to
determine the full two-port response for arbitrary order; the additional
joint two/three-path data is retained in [the response engine](../response/RESULT.md).

## Inherited universe, new observer

The complete class bank is the existing
[order-ten augmentation packet](../../consecutive_flags_order10_20261028/GENERATION.md).
Its alternate full-child canonicalization established uniqueness and labelled
mass2^45, including every automorphism order. This session does not claim a
new enumeration of A000568 and does not re-enumerate that bank.

The bank represents9,733,056 classes as increasing attachment addresses to
191,536 nine-vertex parents. Both compressed and decompressed source hashes
are checked before decoding. Edge bits are lexicographic i<j, with1 meaning
i->j; the nine attachment bits give the NEW vertex's outgoing neighbours.
The full run requires exactly9,733,056 increasing addresses and end of stream.

The observer records H, each named start/end count S,E, and ordered two-path
covers C_uv ending the first path at u and starting the second at v. H=sum S
is redundant at a fixed order. Session27's actual two-port query battery
recovers this sidecar, so its injectivity suffices for the claimed result.

## An exact endpoint sieve makes the larger census practical

For every class, first store the **full sorted list of ten pairs(S_u,E_u)**.
No probabilistic hash or approximate equality is used. Sort these80-byte
keys and partition the bank into equal-key groups.

If a group has one class and all its vertex pairs are distinct, its sidecar
cannot be isomorphic to another class's sidecar, and every sidecar automorphism
fixes all ten vertices. This is a proved sufficient shortcut.

Every other class receives full C computation and exact sidecar canonicalization:
refine using the invariant vertex colours and ordered C-edge multisets, then
enumerate **all** permutations inside the remaining colour cells. Keep every
minimizing labelling, through the minimum and maximum underlying tournament
codes. Distinct endpoint-key groups cannot have isomorphic full sidecars, so
comparison of canonical full keys is needed only within each equal-key group.

Two checks are essential: different input classes have different full keys;
and every minimizing labelling produces the same underlying tournament code.
The latter says Aut(sidecar) is contained in Aut(tournament); the reverse
containment holds by construction. Class separation plus this automorphism
condition gives **labelled** injectivity, exactly as in the prior order-nine
proof. Class separation alone would be insufficient.

## Complete receipt

| Quantity | Exact value |
|---|---:|
| Tournament classes | 9,733,056 |
| Endpoint-multiset groups | 9,726,272 |
| Groups containing multiple classes | 493 |
| Largest endpoint group | 49 |
| Classes with repeated vertex endpoint pairs | 227,764 |
| Full ordered-cover canonicalizations | 227,772 |
| Canonical leaves visited in those cases | 3,272,365 |
| Largest individual leaf search | 362,880 |
| Full-sidecar class collisions | 0 |
| Extra full-sidecar automorphisms | 0 |
| Sum10!/Aut(sidecar) | 35,184,372,088,832 = 2^45 |

The full automorphism histogram agrees with the inherited bank, including
both classes of order81. There are no exact S=E endpoint-balanced classes
at order10; this finite observation does not prove their absence at all even
orders. A stronger parity shortcut is refuted at order8 in the audit lane.

Arithmetic uses unsigned32-bit path/cover entries. Every endpoint count is
at most9!, and the total number of ordered two-path covers is at most
9*10!=32,659,200<2^32. Every nonnegative partial sum and product contributes
to such a count. Adjacency codes need45 bits; group/mass/census counters use
unsigned64-bit integers. No floating-point comparison is used.

## Independent controls and reproduction

The [full tensor engine](../response/RESULT.md) independently enumerates all
10! vertex permutations and their admissible cuts on six specified order-ten
controls. Its endpoints and ordered C agree with this census. It also checks
refined symmetry against all10! labellings, retaining all minima; transitive
and order81 examples supply distinct symmetry controls. The complete order10
conclusion still relies on the proved sieve/canonicalization and inherited
complete bank, not on extrapolating those six controls.

From the repository root:

```sh
python3 05-knowledge/results/two_port_beyond_nine_20261103/census/run.py
```

The frozen complete output is `full10.json`. This is one full C++-O3 census,
not two different C++ optimization audits. The wrapper's bounded100000-class
normal/optimized replays are separate controls and are labelled incomplete
in their own receipts. An existing C++17 compiler named g++ is required;
nothing is installed. Temporary binaries and decoded banks are removed by
the wrapper. The explicit `--limit` option is for bounded controls only and
never emits `complete_bank:true`.
