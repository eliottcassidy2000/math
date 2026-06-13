# THM-008b: General mu Triviality Bound (generalizes THM-008)

**Type:** Theorem
**Certainty:** 5 — PROVED (elementary + exhaustive verification)
**Status:** PROVED
**Last reviewed:** opus-2026-03-05-S2
**Generalizes:** THM-008 (extends from 3-cycles at n<=5 to all L-cycles with L >= n-2)
**Disputes:** none
**Tags:** #mu-triviality #n5-mystery #open-q-001 #claim-a

---

## Statement

**Theorem (mu triviality bound).** Let T be a tournament on n vertices, v a vertex, and C a directed odd L-cycle through v. Then:

```
mu(C) = 1   whenever   L >= n - 2.
```

In particular:
- At n <= 5, mu(C) = 1 for ALL directed odd cycles C through v (both 3-cycles and 5-cycles).
- At n <= 7, mu(C) = 1 for all directed 5-cycles C through v.
- The first n at which mu(3-cycle) can exceed 1 is n = 6.
- The first n at which mu(5-cycle) can exceed 1 is n = 8.

---

## Proof

By definition, mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2).

The conflict graph Omega(T-v)|_{avoid C\{v}} is built on directed odd cycles of T-v that are vertex-disjoint from C\{v}. Such cycles must use vertices from the "available" set:

```
Available = V(T-v) \ C\{v}
|Available| = (n-1) - (L-1) = n - L
```

An odd cycle requires at least 3 vertices. If |Available| < 3, i.e., n - L < 3, i.e., L > n - 3, i.e., L >= n - 2 (since L is odd and n-2 could be even, but the condition L > n-3 suffices), then no odd cycles exist among the available vertices. The conflict graph is empty, and I(empty, 2) = 1.

**For n = 5:**
- 3-cycles (L=3): |Available| = 5 - 3 = 2 < 3. So mu = 1.
- 5-cycles (L=5): |Available| = 5 - 5 = 0 < 3. So mu = 1.

**For n = 6:**
- 3-cycles (L=3): |Available| = 6 - 3 = 3 >= 3. mu CAN exceed 1.
- 5-cycles (L=5): |Available| = 6 - 5 = 1 < 3. So mu = 1.

QED.

---

## Verification Record

| n | Cycle length | mu values observed | Matches theorem |
|---|-------------|-------------------|-----------------|
| 5 | 3 | {1: 7680} | Yes |
| 5 | 5 | {1: 3840} | Yes |
| 6 | 5 | {1: 15200} (sample) | Yes |

---

## Corollary: Resolution of OPEN-Q-001

At n = 5, since mu(C) = 1 for all cycles C through v, Claim A simplifies to:

```
H(T) - H(T-v) = 2 * (N_3 + N_5)
```

where N_L = number of directed L-cycles through v. This is a pure cycle-counting identity with no mu weights, which explains why the per-path identity (a 3-cycle counting mechanism) can "see" the full RHS: the 5-cycle contributions have the same trivial weight (mu=1) as the 3-cycle contributions.

At n = 6, mu(3-cycle) can exceed 1, introducing non-trivial weights that the per-path identity's 3-cycle counting cannot capture. This is why the per-path identity breaks at n = 6.

---

## Notes

- The "n=5 mystery" (OPEN-Q-001) is not a deep structural coincidence. It's a consequence of a simple vertex-counting argument: at n=5 there aren't enough vertices left for the mu computation to be non-trivial.
- The real frontier for understanding Claim A is n = 6, where mu weights become non-trivial for 3-cycles.
- THM-008 (upstream) proved the n<=5 case for 3-cycles. This generalization (opus-2026-03-05-S2) extends to all cycle lengths.
