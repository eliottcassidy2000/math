# Pre-registration: first creators of odd Hamiltonian-path counts at order 10

kind-pasteur-2026-11-07-S16. Reserves **THM-4922** and **HYP-9512**.

## The owner's conjecture

The owner prompt (2026-11-07) asks to merge previous hypotraceable and hypohamiltonian work, and to follow, for every
odd number except 7 and 21, the structure that first created it as the order grows. The conjecture is that each value
is first created by a single isomorphism class or by exactly two classes that are converses of each other. The odd
values would then split into two kinds, and the splitting pattern should be studied. Larger structures might be
generated from smaller ones.

## Definitions

- The least order n(h) is the first order at which some tournament has H = h.
- A **first creator** of h is a class of order n(h) with H = h. All first creators are strong, the set of them is closed
  under converse, and so it is s SC classes plus p converse pairs.
- k(T) is the number of vertices v such that T - v is itself a first creator, of its own value, at order n(h) - 1.

## Data before this registration (orders <= 9, exploratory)

The class census comes from an independent generator, `classgen.c`, using individualization-refinement canonical
forms. It reproduces A000568 and klein's level-resolved census exactly at orders 3-9, as (h, c3, strong) class counts.

1. **The dichotomy fails.** Of the 1,520 values with n(h) <= 9, only 69 have (s, p) = (1, 0) or (0, 1). The
   smallest counterexample is h = 23, with two converse pairs at order 6. Fiber sizes reach 2,458.
2. **Kinds.** 57 values are SC-only, 371 NS-only, and 1,092 mixed.
3. **Genealogy.** Of the 169,501 first creators at orders 4-9, all but 2 have k(T) >= 1. The two exceptions are the
   converse pair of first creators of h = 63 at order 8, which are orphans. The value 63 is the only one with no
   first-creator lineage.
4. **Full criticality.** k(T) = n holds for 24%, 39% and 50% of the first creators at orders 7, 8 and 9.
5. **No hypo analogue.** No class that is not a first creator has every vertex deletion a first creator.

## Predictions for order 10 (`order10.c`, run after this file is pushed)

- **Q1.** No orphan first creator at order 10: every first creator has k(T) >= 1. So 63 stays the only orphan value
  through order 10.
- **Q2.** No hypo analogue at order 10.
- **Q3.** The fraction of first creators with k(T) = 10 exceeds 50.4%, the order-9 value.
- **Q4.** The dichotomy holds for at most 5% of the values new at order 10.
- **Q5.** Mixed fibers, with at least one SC class and at least one converse pair, are the majority of the values new
  at order 10.

Any failure is reported as a failure.
