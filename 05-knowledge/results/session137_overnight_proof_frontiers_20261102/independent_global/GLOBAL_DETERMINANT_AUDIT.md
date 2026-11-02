# Independent audit: global Hamiltonian determinant inequality

**Verdict: ACCEPT the all-order proof, weighted extension, equality classifications,
and the unweighted cycle-defect strengthening.** This is an independent mathematical
review with separate exact implementations, not a formal proof-assistant certificate.
The computations below are controls for the proof, not finite evidence substituted
for universal quantifiers.

Baseline inspected read-only: `bb7cb2d1f61a14632389814469317131d897b640`.
Reviewed final manuscript: Session137 `discriminant/GLOBAL_PROOF.md`, 10,860 bytes,
SHA-256 `f7e79e9aa935533f2a0444ce2bf6b220f5ff1cb29fea16adc1fa19a7cc71ba6f`.
The final text, including its added Sections 5–6, was read and accepted.
No Git mutation, fetch, or publication was
performed by this reviewer. All reviewer files are in this independent packet.

## 1. Accepted statements and indispensable scope

For a **nonempty** tournament with adjacency matrix A, put

```
H = number of directed Hamiltonian vertex words,
u = det(I+A),             v = (-1)^n det(A),
D = det(I+A-A^T)/2^(n-1), s = 1^T(I+A-A^T)^(-1)1.
```

The accepted conclusion is

```
H >= u+|v| = max(1,s)D >= D.
```

For ordinary zero/one tournaments, equality in the first inequality holds exactly
for transitive tournaments and C3. Equality `H=D` holds exactly for transitive
tournaments. This closes the stronger residual inequality in the inherited
THM-1950 route, not the separate Hamiltonian value-spectrum problem.

The same inequality holds for real **weighted tournament matrices** with zero
diagonal, `0<=A_ij<=1`, and `A_ij+A_ji=1`. Here H is the sum of products of arc
weights along Hamiltonian words. Its complete equality classification is:

* order one: equality;
* order two: equality `H=u+|v|` for every edge probability;
* order at least three: equality only for deterministic transitive tournaments
  or the deterministic three-cycle.

In the weighted setting too, `H=D` holds exactly for deterministic transitive
tournaments. No integer-cycle defect statement below is exported to weighted
matrices.

For the empty object use `H=u=v=1`. The strong inequality is **false** there
(`1<2`), while the weak induction assertion `H>=u` is true. This exception is
load-bearing and explicitly retained in every recurrence audit.

## 2. Independent proof review

### Algebraic positivity is not circular

Let `M=I+A-A^T`, `J=11^T`. The identities

```
M+J=2(I+A),       M-J=-2A^T
```

and the determinant lemma imply `u+v=D>0`, `u-v=sD>0` for every nonempty
induced submatrix. Indeed M is identity plus real skew, and for `z=M^-1 1`,
`s=1^Tz=||z||^2>0`. These are algebraic facts established before the H bound.
For the empty remainder the same two quantities are 2 and 0, nonnegative.

### Determinant signs and endpoint orientation checked independently

Add x to parent T. An old directed path P starts at an x-outneighbor and
ends at an x-inneighbor. Its closure through x is one distinguished cycle,
with no rotation factor. If P has ell old vertices, that cycle contributes
`(-1)^ell` to `det(I+A)`. In `(-1)^n det(A)`, every spanning-cover cycle
has sign -1, independent of its length. Thus, writing E/O for even/odd ell
and `R=T-V(P)`, the exact identities are

```
u_child = u_T + sum_E u_R - sum_O u_R,
v_child = -sum_(E union O) v_R,
H_child = H_T + 2 sum_E H_R.
```

The H identity uses the inherited odd-cycle collection formula, not a
determinant approximation. In the weighted case every term has the product
weight of the old path and its two incident x-arcs. Length-one P must then be
allowed: it represents a weighted two-cycle and occurs in the determinant
recurrences but not in OCF. Those additional O terms are beneficial, not an
omitted error term.

For either `epsilon=+1,-1`, subtraction gives

```
H_child-u_child+epsilon*v_child
 = H_T-u_T
   + sum_E [2(H_R-u_R)+u_R-epsilon*v_R]
   + sum_O [u_R-epsilon*v_R].
```

Induction on **only** `H>=u`, starting with the empty object, makes every
summand nonnegative. Both signs simultaneously prove `H>=u+|v|` for the
nonempty child and propagate the weaker induction assertion. No strongness,
feedback bound, finite census, or unsupported sign assumption is needed.

Weighted OCF follows by taking expectations in independent random tournaments:
each Hamiltonian path and vertex-disjoint odd-cycle collection uses distinct
unordered edges. A naive expectation argument for the determinants would be
invalid, because a fractional matrix can have nonzero two-cycle products
`A_ij A_ji`; the proof correctly retains these separately.

### Strictness and equality

In a nontransitive binary tournament of order at least four, choose x on a
triangle. The other two vertices give an E term whose remainder is nonempty;
`u_R-epsilon*v_R>0` for both signs. C3 and transitive cases are computed directly.

For a weighted matrix of order at least three, any fractional edge x-i gives
a positive length-one path weight `A_xi A_ix`. Its remainder has order at
least one, so the corresponding O term is strictly positive for both signs.
At order two, direct calculation gives

```
H=1, u=1-p(1-p), v=-p(1-p), D=1-2p(1-p).
```

These observations prove the equality classifications above, including the
fractional order-two hostile.

## 3. Unweighted cycle-defect strengthening

Let `C_all` count all literal directed simple cycles (all lengths at least
three, modulo rotation), and let `C_proper` exclude Hamiltonian cycles. Then

```
H-u >= C_all,
H-u-|v| >= C_proper.
```

Equality in either statement holds exactly when there are no two
vertex-disjoint directed triangles. If there are two, each corresponding
excess is at least two.

Independent proof check: the recurrence for `H-u` contributes at least one
per new cycle, even for empty remainder. In each signed recurrence, each
proper cycle through x has nonempty remainder and contributes at least one
because `u_R±v_R` is a positive integer. Cycles avoiding x are counted by
the parent weak defect. Two disjoint triangles yield an E term with cyclic
remainder and additional `2(H_R-u_R)>=2`. Conversely, when no two triangles
are disjoint, no two directed cycles are disjoint. The single-cycle expansions
give `H-u=C_all` and `v=-C_Hamiltonian`, proving both equalities.

## 4. Exact independent replay receipts

No producer module is imported by these checks.

| Check | Exact universe / consequence | Paired runs |
|---|---|---|
| `global_gap.cpp` | 34,380 tournaments: every labelled order1..6; 128 fixed-seed samples at each order7..10; exact strong `s<1` witness. 42,532 root insertions, 4,748,109 literal old paths, 2,342,526 induced full-matrix determinant checks, 654,234 empty-complement terms. All three insertion identities, both gaps, positivity, ordinary equality. | O2 / UBSan byte-identical |
| `audit_weighted_gap.py` | All edge probabilities in `{0,1/2,1}` through order4; sixteen fixed-seed rational matrices at each order5..7. 808 matrices, 3,292 insertions, 70,745 positive-weight paths, 3,742 two-cycle terms, 19,870 empty-complement terms, 14,670 induced determinant checks. Weighted identities and full equality classification. | normal / `-O` byte-identical |
| `audit_defects.cpp` | Every labelled order1..6 plus 128 fixed-seed examples at each order7..10: 34,379 tournaments, 15,931 with two disjoint triangles, 1,147,738 literal cycles and 1,024,879 proper cycles. Both defect bounds, exact equality locus, and excess at least two. | O2 / UBSan byte-identical |

The C++ determinants are **full-matrix modular Gaussian elimination**, not the
producer's Bareiss, Schur, or small residual formulas. At order at most ten,
Hadamard bounds the absolute determinant by `10^5`, below half the prime
`1,000,003`; symmetric-residue reconstruction is therefore an exact integer
computation, not a probabilistic modular test. Weighted determinants use
independent `Fraction` Gaussian elimination. Path counts use subset DP;
insertion paths and cycle counts use separate literal recursions.

Output SHA-256:

```
7ea06c6757380f62310fb5768bd87c0b906b0a20aa3d82072b079b4e7711b6cc  global_gap.normal.json
7ea06c6757380f62310fb5768bd87c0b906b0a20aa3d82072b079b4e7711b6cc  global_gap.ubsan.json
74c80b7bedbbdfb2363dd83600fa70d694af5e592ae1900ede43f5297e30554c  weighted_gap.normal.json
74c80b7bedbbdfb2363dd83600fa70d694af5e592ae1900ede43f5297e30554c  weighted_gap.optimized.json
```

The defect hashes and all current source hashes are in `GLOBAL_SHA256SUMS`.
`GLOBAL_FILES.txt` lists the complete portable first-checkpoint subset.

## 5. Reproduction

From this directory, with Python 3 and a C++17 compiler:

```
g++ -std=c++17 -Wall -Wextra -Werror -O2 global_gap.cpp -o /tmp/global_gap
/tmp/global_gap > /tmp/global_gap.json
cmp /tmp/global_gap.json global_gap.normal.json
g++ -std=c++17 -Wall -Wextra -Werror -O1 -fsanitize=undefined global_gap.cpp -o /tmp/global_gap_ubsan
/tmp/global_gap_ubsan > /tmp/global_gap.ubsan.json
cmp /tmp/global_gap.ubsan.json global_gap.normal.json
python3 -B audit_weighted_gap.py > /tmp/weighted_gap.json
python3 -B -O audit_weighted_gap.py > /tmp/weighted_gap.optimized.json
cmp /tmp/weighted_gap.json weighted_gap.normal.json
cmp /tmp/weighted_gap.optimized.json weighted_gap.normal.json
g++ -std=c++17 -Wall -Wextra -Werror -O2 audit_defects.cpp -o /tmp/audit_defects
/tmp/audit_defects > /tmp/defects.json
cmp /tmp/defects.json defects.normal.json
```

All Python validity checks use explicit exceptions, not removable assertions.

## 6. Correction lineage and boundaries

The false THM-1950 supermultiplicativity sentence was independently refuted:
strong rows `[30,8,122,80,2,91,19]` have `(u,v,D,s)=(5,1,6,2/3)`; the ordinal
square has `(25,1,26,12/13)`, so `26<6*6`. The exact SCC composition law and
its invariant reduction survive. Also `C3=>C3=>C3` has `s=9/7>1`, refuting
the overbroad claim that every reducible response lies below one.

The new proof repairs and surpasses that reduction; it does not use either
false sentence. Neither H-spectrum completeness, LRC(14), efficient exact
Hamilton counting, nor external literature priority follows from this audit.
