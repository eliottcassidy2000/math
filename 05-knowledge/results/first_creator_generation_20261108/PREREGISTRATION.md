# Pre-registration: first-creator genealogy at order 11, and generation of all tournaments from first creators

kind-pasteur-2026-11-08-S17. This continues THM-4922 and HYP-9512 (kind-pasteur S16). Reserved IDs: **THM-4942** and
**HYP-9532**. This file is frozen once pushed. Results go to README.md and PROOF.md, with an erratum section if needed.

Owner prompt (2026-11-08): run the order-11 test. For every odd number besides 7 and 21, keep track of the structure
that first created it as the order grows. See how all tournaments can be generated from this minimal set if you get
creative enough with the operations. Relate the objects used in the proofs that 7 and 21 never appear to fundamental
concepts and constructions.

## Definitions

- S(n) is the set of Hamiltonian-path counts H of tournaments of order n.
- A **first creator** is a class of order n with H not in S(n - 1).
- An **extension pattern** of an order-n class T is one of the 2^n labelled ways to add a vertex v. Write
  x_a = [a -> v].
- A pattern is **fertile** if H(T + v_x) is not in S(n), that is, the extension is itself a first creator.

## Tool: `ext.c`

`ext.c` evaluates H on every extension through the quadratic attachment form

    H(T + v_x) = H(T) + sum_a x_a L_a - sum_{a<b} x_a x_b M_ab,
    W[a][b] = sum_{S: a in S, b not in S} f[S][a] g[V\S][b],
    L_a = f[V][a] - g[V][a] + sum_{b != a} W[a][b],   M_ab = W[a][b] + W[b][a].

Here f[S][a] (respectively g[S][a]) counts the Hamiltonian paths of T[S] ending (respectively starting) at a. The
patterns are visited in Gray order with O(n) updates per step. For a first creator T, the engine also carries the forms
of the n deletions T - u, which give every H(T' - u) for T' = T + v_x.

(Degree two in the attachment variables is folklore: each Hamiltonian path uses at most two arcs at v; cf. Session189's
response polynomials and the insertion calculus of THM-012b and THM-1900.)

## Data seen before this registration

**Validation of `ext.c`.**
- Extending every class of order 7, 8 and 9 reproduces S(8), S(9) and S(10) exactly. S(10) equals klein's 7,663
  values.
- Brute-force subset-DP checks of every pattern found 0 failures: all 456 classes of order 7, 300 classes of order 8
  and 160 of order 9, including every deletion count.

**Genealogy at the value level (HYP-9512 (A)).** The only new value not reached from a first creator one order down
is 63, at order 8. That is the known orphan, so the gate passed. None is unreached at orders 9 and 10.

**Hypo pairs (HYP-9512 (B)).** There are 0 pairs at orders 8, 9 and 10. A hypo pair is a first creator T and a pattern
x such that T + v_x is not a first creator but every deletion of it is.

**Fertility at orders 3-9.**
- No first creator is sterile: every one has a fertile pattern.
- Every first creator has at least half of its patterns fertile, with one exception. The two first creators of 63
  (order 8) have only 77 of 256.

**Share of new values reached only from first-creator parents.** These are the values new at order n + 1 for which
every (class, pattern) realization from order n has a first-creator parent: 85/243 = 0.350 at order 8, 510/1200 =
0.425 at order 9, and 2886/6143 = 0.470 at order 10.

**Generation census at orders 3-9 (`gen.c`).** Every class is classified, in this priority order:

| Category | Meaning |
| --- | --- |
| FC | a first creator |
| EXT | not FC, but some deletion is a first creator |
| SUM | otherwise, not strong (an ordinal sum) |
| SUB | otherwise, a nontrivial module (a substitution) |
| REST | strong, prime, not FC, with no first-creator deletion |

The counts are:

| order | non-first-creators | EXT | SUM | SUB | REST |
| --- | --- | --- | --- | --- | --- |
| 8 | 1,423 | 1,162 | 216 | 33 | 12 |
| 9 | 27,860 | 23,855 | 2,489 | 1,010 | 506 |

At orders 3-7, REST is empty.

The REST classes:
- Every REST class at order 8, and 505 of the 506 at order 9, has a **pendant** vertex: in-degree 1 or out-degree 1.
- The single exception is an order-9 class with H = 385. It is self-converse with trivial automorphism group and
  arc-strong connectivity 2: two strong 4-tournaments in ordinal sum, plus a middle vertex that splits both.

**Also consulted.**
- klein's `least_order_11.txt`: S(11), n(h) and the strong windows. It defines the new values at order 11 but has no
  deletion or parent data.
- mac-mini's THM-4930 checkpoint (self-converse colours). It is a different question.

## Not yet seen

- the order-10 generation census;
- order-10 fertility;
- every order-11 extension datum.

## Decisive runs

1. `ext.c` over all 9,733,056 classes of order 10 and all 1,024 patterns each (9.97 x 10^9 extensions).
2. `gen.c` over all 9,733,056 classes of order 10.

## Predictions

| ID | Prediction | Expected |
| --- | --- | --- |
| R1 | HYP-9512 (A) at order 11: every value with n(h) = 11 is H(T + v) for some first creator T of order 10. So no new orphan value appears at order 11. | CONFIRMED |
| R2 | HYP-9512 (B) at order 11: 0 hypo pairs. No order-11 tournament whose H lies in S(10) has every vertex deletion a first creator of order 10. | CONFIRMED |
| R3 | No first creator of order 10 is sterile. | CONFIRMED |
| R4 | Every first creator of order 10 has at least 512 fertile patterns out of 1,024. The only exceptions so far are the orphan creators of 63, and order 10 has no orphan. | CONFIRMED |
| R5 | The share of new order-11 values reached only from first-creator parents lies in [0.48, 0.56]. | CONFIRMED |
| R6 | REST at order 10 is between 2% and 6% of the non-first-creators (0.84% at order 8, 1.82% at order 9). | CONFIRMED |
| R7 | At least 99% of REST classes at order 10 have a pendant vertex. | CONFIRMED |
| R8 | At most 100 order-10 classes lie outside FC, EXT, SUM, SUB and the pendant extensions (1 at order 9). | CONFIRMED |
| R9 | EXT is at least 85% of the non-first-creators at order 10 (81.7% at order 8, 85.6% at order 9). | CONFIRMED |

A failure of R1 or R2 refutes HYP-9512 (A) or (B) at order 11. The others are exploratory shares.

Files: `ext.c`, `gen.c`, `dump10.c` and `merge.py`, in this folder.
