# THM-4922: first creators of odd Hamiltonian-path counts - converse fibers, the two kinds, and genealogy

kind-pasteur-2026-11-07-S16. Status by section:
- **PROVED**: sections 1 and 4.1.
- **FINITE-EXACT**: orders <= 9 (exploratory) and order 10 (pre-registered Q1-Q5,
  [PREREGISTRATION.md](PREREGISTRATION.md), pushed in 9e07033ca4 before the order-10 run).
- **CONJECTURE**: HYP-9512.

Owner prompt (2026-11-07): merge the hypotraceable and hypohamiltonian work, and follow the structure that first creates
each odd value except 7 and 21 as the order grows. The conjecture is that each value is created by a single isomorphism
class or by exactly two classes that are converses. The odd values would then split into two kinds, and larger
structures might be generated from smaller ones.

**Data.** `classgen.c` is an independent class generator: vertex extension plus individualization-refinement canonical
forms, with no nauty. It reproduces A000568 (1, 1, 2, 4, 12, 56, 456, 6880, 191536) and klein's level-resolved census
(`../path_spectrum_two_angles_metagraph_20261106/data/census_n.txt`) exactly at orders 3-9, as (h, c3, strong) class
counts. For every class it records:
- H, c3 and strong connectivity;
- self-converse (SC) status and |Aut|;
- the class and H of every vertex-deleted subtournament.

`order10.c` extends the census to all 9,733,056 classes of order 10, for the registered tests.

Credited prior work:
- klein's least-order census (n(h), strong windows, N(n) ⊆ Ss(n));
- THM-791 (SC tournaments of even order have H = 1 mod 4);
- THM-344 (H = 63 occurs at order 8 in exactly two classes);
- THM-4870 (the Thomassen merge: non-critical vertices, no hypohamiltonian or hypotraceable tournament);
- THM-1370 (7 and 21; the insertion lemma).

## 1. Elementary facts (PROVED)

Let n(h) be the least order of h, and call the classes of order n(h) with H = h the **first creators** of h.

1. **Strong, and closed under converse.** Every first creator is strong. Otherwise H is a product of the counts of its
   strong components, and removing a one-vertex component would realise h at a smaller order. The first creators are
   closed under converse, since H(T^op) = H(T). So they are s SC classes plus p converse pairs.
2. **Endpoint deletion.** Deleting the last vertex of a Hamiltonian path is an injection
   HP(T) -> union over v of HP(T - v), so H(T) <= sum over v of H(T - v).
3. **Genealogy above a threshold.** If h > n max H(n - 2) with n = n(h), then every first creator T of h has a vertex v
   such that T - v is a first creator (of H(T - v), at order n - 1).

   *Proof.* By 2, some v has H(T - v) >= h/n > max H(n - 2). Hence H(T - v) is not attained at order n - 2, and its
   least order is n - 1. QED.

## 2. The conjectured dichotomy is false

| | |
| --- | --- |
| values with n(h) <= 9 (h >= 3) | 1,520 |
| first creators: exactly one SC class, (s, p) = (1, 0) | 44 values |
| first creators: exactly one converse pair, (s, p) = (0, 1) | 25 values |
| everything else (several classes) | 1,451 values; fiber sizes up to 2,458 |

- **Smallest counterexamples.**
  - h = 23: two converse pairs at order 6.
  - h = 29: two SC classes and one pair.
  - h = 33 and h = 37: one SC class and two pairs.
  - Refining by the least level c3, or by the number of Hamiltonian cycles, does not repair it: 5,691 values still have
    more than two classes at their least order and least level, out of 7,663 values with n(h) <= 10, from klein's
    census.
- **Where the dichotomy holds.** It holds at the small values and near the top of each order's range, where few classes
  exist.
  - One SC class: 1, 3, 5, 11, 13, 41, 61, 85, 87, 125, 127, 147, 155, 157, 171, 175, 189 (= max H(7)), 629, 637,
    2891, ..., 3357 (= max H(9)).
  - One converse pair: 19, 25, 27, 31, 35, 39, 43, 47, 63, 601, ..., 633, 2885, ..., 3109.

## 3. The two kinds and the splitting pattern

Let a value be **blue** if some first creator is SC, and **black** if all of them are converse pairs. Through order 9
there are 1,149 blue and 371 black values. The split is governed by the parity of the least order and by h mod 4:

| least order | h mod 4 | blue | black |
| --- | --- | --- | --- |
| even | 1 | 89 | 43 |
| even | 3 | 0 | 124 |
| odd | 1 | 521 | 103 |
| odd | 3 | 538 | 101 |

- The empty cell is THM-791 (KNOWN): at even order every SC class has H = 1 (mod 4), so a value 3 (mod 4) first
  created at an even order is black.
- **Observation (FINITE-EXACT, n <= 9, not found stated in the repository).** At odd orders the SC classes are split
  exactly evenly between H = 1 and H = 3 (mod 4): 1/1, 4/4, 44/44 and 1376/1376 at n = 3, 5, 7, 9. At even orders all
  of them are 1 (mod 4), as THM-791 says.
- The SC-only values cluster at the top of each order: median h/max H(n(h)) is 0.918, against about 0.55 for black and
  mixed values.

## 4. Genealogy: larger first creators grow from smaller ones (the hypohamiltonian merge)

For a first creator T of order n, let k(T) be the number of vertices v for which T - v is a first creator at order
n - 1. In THM-4870's language, these are the vertices whose deletion stays in the class.

| order | first creators | k = 0 (orphans) | k = n (fully critical) |
| --- | --- | --- | --- |
| 7 | 328 | 0 | 78 (24%) |
| 8 | 5,457 | 2 | 2,118 (39%) |
| 9 | 163,676 | 0 | 82,553 (50%) |
| 10 | 8,647,232 | 0 | 5,433,631 (62.8%) |

- **Orphans.** Through order 9 the only orphans are the two first creators of h = 63, a converse pair at order 8
  (THM-344's two classes). 63 is the isolated hole of order 7 (klein), and all its vertex deletions have values already
  present at order 6. Every other value has a first creator built from a first creator by adding one vertex. This is
  the owner's "generate larger structures from smaller ones", now an exact statement.
- **No hypo analogue.** A hypohamiltonian graph is non-Hamiltonian although every vertex deletion is Hamiltonian, and
  no tournament is hypohamiltonian or hypotraceable (THM-4870). Likewise, through order 9 no class that is not a first
  creator has every deletion a first creator.
- **Proved part.** Section 1.3 proves the genealogy for every value above n max H(n - 2): 1701 at n = 9, 6610 at
  n = 10. The data show it holds far below that threshold, except at 63.

## 4.1 Order 10 (pre-registered)

The full order-10 census (`order10.c`, [order10.out](order10.out), 750 s) has:
- 9,733,056 classes (A000568);
- 8,647,232 first creators;
- 6,143 new values, matching klein's |N(10)|.

| Prediction | Outcome |
| --- | --- |
| Q1: no orphan first creator at order 10 (63 stays the only orphan value) | **CONFIRMED** (k = 0 never occurs; the minimum is k = 2, for 41 classes) |
| Q2: no hypo analogue | **CONFIRMED** (0 classes) |
| Q3: full criticality above 50.4% | **CONFIRMED** (62.8%) |
| Q4: the dichotomy for at most 5% of new values | **CONFIRMED** (55 of 6,143, all in [14675, 15745]; max H(10) = 15745 has a unique SC class) |
| Q5: mixed fibers are the majority | **FAILED**: SC-only 5, black 3,605, mixed 2,533 |

**Why Q5 failed: THM-791.** At even order every SC class has H = 1 (mod 4), so the 3,067 new values that are 3
(mod 4) are all black. Among the 3,076 values that are 1 (mod 4), 2,538 (82.5%) are blue, close to the odd-order rate
(84%). The prediction was registered before THM-791 was connected to the split, so this failure is honest. The lesson
is in HYP-9512.
