# THM-4922: first creators of odd Hamiltonian-path counts - converse fibers, the value colours, and genealogy

kind-pasteur-2026-11-07-S16. Status by section:
- **PROVED**: section 1 (as corrected), the proved part of section 4, and the odd-order SC balance (section 3; proof
  by the audit).
- **FINITE-EXACT**: sections 2, 3, 4 and 4.1.
- **CONJECTURE**: HYP-9512.

The independent audit is PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md); MISTAKE-735). Pre-registration:
[PREREGISTRATION.md](PREREGISTRATION.md) (9e07033ca4).

Owner prompt (2026-11-07): merge the hypotraceable and hypohamiltonian work, and follow the structure that first
creates each odd value except 7 and 21 as the order grows. The conjecture is that each value is created by a single
isomorphism class or by exactly two classes that are converses. The odd values would then split into two kinds, and
larger structures might be generated from smaller ones.

**Data.**
- `classgen.c` is an independent class generator: vertex extension plus individualization-refinement canonical forms,
  with no nauty. It reproduces A000568 and klein's level-resolved census at orders 3-9.
- `order10.c` covers all 9,733,056 classes of order 10.
- The audit rebuilt every class through order 10 with its own canonical form, by a different method. It matches
  everywhere, including all 6,143 order-10 value lines.

**Credited prior work.**
- **Session224** (`session224_first_occurrence_20261107`, 2054513960, same owner prompt, pushed first): the refutation
  of the dichotomy at h = 9, the (s, p) birth ledger through order 8, and the strong-component structure of minimal
  realizations.
- **klein's least-order census:** n(h), and the first-creator fiber sizes for n(h) <= 10 in `least_order_10.txt`.
- **THM-4889 and klein section 4:** the counting identity H(T) <= sum over v of H(T - v).
- **THM-791 and Session29 converse congruences:** the even-order black cell.
- **63:** THM-344, THM-4748 and the Session24 multiplicities.
- **THM-1370:** 7 and 21, and the insertion lemma.
- **No hypohamiltonian or hypotraceable tournament:** Rédei and Camion, as merged in THM-4870's source packet
  `thomassen_merge_hypohamiltonian_groups_20261104.md`, section 2.

## 1. Elementary facts (PROVED)

Let n(h) be the least order of h. The classes of order n(h) with H = h are the **first creators** of h.

1. **No one-vertex strong component, and closed under converse.**
   - A first creator of h > 1 has no one-vertex strong component. H is the product of the counts of the strong
     components, and deleting a one-vertex component would realise h at a smaller order.
   - For the same reason, each strong component C is itself a first creator of H(C). Replacing it by a smaller
     tournament with the same count keeps H (Session224 section 3).
   - So a first creator is either strong, or an ordinal sum of at least two strong first creators, each of order at
     least 3. The second case never occurs through order 11. Order 10 was recomputed by the audit; at order 11, klein's
     census_11 has no non-strong class with a new value. Excluding it in general is HYP-9481 (b), which is OPEN.
   - The first creators are closed under converse, since H(T^op) = H(T). So they are s SC classes plus p converse pairs.
2. **Endpoint deletion (THM-4889).** Deleting the last vertex of a Hamiltonian path is an injection, so
   H(T) <= sum over v of H(T - v).
3. **Genealogy above a threshold.** If h > n max H(n - 2) with n = n(h) >= 3, then every first creator T of h has a
   vertex v such that T - v is a first creator.

   *Proof.* By 2, some v has H(T - v) >= h/n > max H(n - 2). So H(T - v) is not attained at order n - 2, nor, by
   padding with a sink (S(m) ⊆ S(m + 1)), at any smaller order. T - v attains it at order n - 1, so T - v is a first
   creator. QED.

## 2. The conjectured dichotomy is false

| | |
| --- | --- |
| values with n(h) <= 9 (all h >= 1) | 1,520 |
| first creators: exactly one SC class, (s, p) = (1, 0) | 44 values |
| first creators: exactly one converse pair, (s, p) = (0, 1) | 25 values |
| everything else (two SC classes, or three or more classes) | 1,451 values (8 with two SC classes); fiber sizes up to 257 (h = 1503); 2,458 at order 10 (h = 7799) |

- **Smallest counterexamples.**
  - h = 9 and h = 15: two SC classes at order 5. Session224 proved h = 9 first.
  - h = 17: two SC classes at order 6.
  - Eight values with n(h) <= 9 fail this way: 9, 15, 17, 45, 145, 621, 3075, 3105.
  - h = 23: two converse pairs at order 6.
  - h = 29: two SC classes and one pair.
  - h = 33 and h = 37: one SC class and two pairs (33 is Session224's recovered H = 33 fiber).
- **Refinements do not repair it.** Out of 7,663 values with n(h) <= 10 (klein's census, recomputed by the audit):
  - at the least order and least level c3, 5,691 values still have more than two classes;
  - at the least order and least number of Hamiltonian cycles, 1,901 do.
- **Where the dichotomy holds.** Mostly at small values and near each order's maximum. At order 7 it also holds at 61,
  85, 87, 125 and 127, which lie at 0.32-0.67 of max H(7).
  - One SC class: 1, 3, 5, 11, 13, 41, 61, 85, 87, 125, 127, 147, 155, 157, 171, 175, 189 (= max H(7)), 629, 637,
    2891, ..., 3357 (= max H(9)).
  - One converse pair: 19, 25, 27, 31, 35, 39, 43, 47, 63, 601, ..., 633, 2885, ..., 3109.

## 3. The value colours and the splitting pattern

A value is **blue** if some first creator is SC, and **black** if all of them are converse pairs. These value colours
are distinct from the metagraph's line colours. Through order 9 there are 1,149 blue values (1,148 with h >= 3, as in
the table) and 371 black values.

| least order | h mod 4 | blue | black |
| --- | --- | --- | --- |
| even | 1 | 89 | 43 |
| even | 3 | 0 | 124 |
| odd | 1 | 521 | 103 |
| odd | 3 | 538 | 101 |

- **The empty cell is THM-791 (KNOWN).** At even order every SC class has H = 1 (mod 4).
- **Odd-order balance (PROVED; proof by the audit, audit section 4).** At every odd order n >= 3, the SC classes split
  exactly evenly between H = 1 and H = 3 (mod 4): 1/1, 4/4, 44/44 and 1376/1376 at n = 3, 5, 7, 9.
  - An SC tournament of odd order has an anti-involution sigma0 fixing one vertex. All of them are conjugate under
    Aut, so SC classes correspond to centralizer orbits of the tournaments anti-invariant under a fixed sigma0.
  - By OCF, H = 1 + 2F (mod 4), where F counts sigma0-symmetric odd cycles through the fixed vertex.
  - Reversing the internal arcs of the pairs whose path-count parity is odd is an equivariant involution.
  - It changes F by an odd number. That number is the count of all directed paths of a pair tournament, odd by Rédei.
- **Blue shares.** The blue share among values not forced black is 0.72, 0.66, 0.84 and 0.825 at orders 7-10.

## 4. Genealogy: larger first creators grow from smaller ones (the hypohamiltonian merge)

For a first creator T of order n, let k(T) be the number of vertices v for which T - v is a first creator at order
n - 1. In THM-4870's language, these are the non-critical vertices: their deletion stays in the class.

| order | first creators | k = 0 (orphans) | k = n (fully non-critical) |
| --- | --- | --- | --- |
| 7 | 328 | 0 | 78 (24%) |
| 8 | 5,457 | 2 | 2,118 (39%) |
| 9 | 163,676 | 0 | 82,553 (50%) |
| 10 | 8,647,232 | 0 | 5,433,631 (62.8%) |

- **Orphans.** At orders 4-10 the only orphans are the two first creators of h = 63 at order 8. These are the Fibonacci
  fan F_6 and its converse (THM-4748), THM-344's two classes. All their vertex deletions have values already present at
  order 6.
  - Among the values h >= 5 with n(h) <= 10, every other value has a first creator built from a first creator by adding
    one vertex.
  - The 3-cycle, at order 3, is also an orphan, since its deletions have H = 1.
  - This is the owner's "generate larger structures from smaller ones": finite-exact through order 10, proved above
    n max H(n - 2) (section 1.3), conjectured in general (HYP-9512 (A)).
- **No hypo analogue.** No class of order 3-10 that is not a first creator has every deletion a first creator. At
  order 2 the unique class trivially is one. Compare: no tournament is hypohamiltonian or hypotraceable.
- **Proved part.** Section 1.3 covers values above the threshold: 1,701 at n = 9 (671 of the 1,200 new values lie above
  it) and 6,610 at n = 10 (4,360 of 6,143).

## 4.1 Order 10 (pre-registered)

The full order-10 census (`order10.c`, [order10.out](order10.out), 750 s) has 9,733,056 classes, 8,647,232 first
creators and 6,143 new values. The audit's own order-10 run matches every value line.

| Prediction | Outcome |
| --- | --- |
| Q1: no orphan first creator at order 10 | **CONFIRMED** (blind; k = 0 never occurs, the minimum is k = 2, for 41 classes) |
| Q2: no hypo analogue | **CONFIRMED** (blind; 0 classes) |
| Q3: full non-criticality above 50.4% | **CONFIRMED** (blind; 62.8%) |
| Q4: the dichotomy for at most 5% of new values | **CONFIRMED, not blind** |
| Q5: mixed fibers are the majority | **FAILED**: SC-only 5, black 3,605, mixed 2,533 |

Q4 was not blind: the registration quoted klein's order-10 maximum fiber (2,458), and klein's fiber sizes bound the
count by 57 (0.93%). The final count is 55 of 6,143, all in [14675, 15745]; max H(10) = 15745 has a unique SC class.

**Why Q5 failed: THM-791.** The 3,067 new values that are 3 (mod 4) are all black. Among the 3,076 that are 1
(mod 4), 82.5% are blue.

**Erratum for the frozen registration.** Its data item 1 has two errors. The smallest counterexample is 9, not 23. And
2,458 is the order-10 maximum fiber; the order-<= 9 maximum is 257.
