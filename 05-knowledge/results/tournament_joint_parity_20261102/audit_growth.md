# Independent audit: the sharp attachment-distribution loss at order five

**PASS. FINITE-EXACT exhaustive sharpness through order four; explicit
order-five counterexample independently reconstructed by named joins.**
The general endpoint formula and its mean consequence below are PROVED.
This audit imports neither `joint_audit.py` nor any repository algorithm.

## 1. Objects and preserved source law

The producer's `joint_audit.py` enumerates all labelled tournaments at orders
0..6, quotients by an exact degree-compatible canonical form, and compares
the 2^n one-vertex attachment distributions of equal-H ordinary classes.
Its canonical masks 2 and 32 at n=5 are independently reconstructed here as

```
U = (0) -> (1) -> (2,3,4 cycle 2->3->4->2),
V = (0) -> (1,2,3 cycle 1->2->3->1) -> (4).
```

The explicit vertex transports into the producer's gauge are, respectively,

```
U: (0,1,2,3,4) -> (4,3,1,0,2),
V: (0,1,2,3,4) -> (4,1,3,2,0).
```

Every arc is checked after this transport. There is no implicit redefinition
of source vertices. Both graphs have H=3 and automorphism order 3; their
triangle rotations are the only automorphisms.

For each named source, b=(b_0,...,b_4) runs through all 32 equally weighted
binary vectors, with b_i=1 meaning i->x and b_i=0 meaning x->i for a new
vertex x. This is the uniform law on **labelled attachment vectors**, not
a fresh uniform law on child isomorphism classes, attachment orbits, or
Hamilton paths. Duplicate target classes retain their original mass.

## 2. Independent exact results

The checker literally enumerates all 6! vertex orders of each of the 64
children. It also counts every ordered spanning k-path cover by choosing a
vertex order and k-1 cut locations; every uncut adjacency must be a source
arc. Thus the counts do not use the producer's subset DP or binomial formula.

Both ordered cover profiles (k=1..5) are

```
[3,42,162,240,120].
```

The attachment histograms are:

| child H | number of U attachments | number of V attachments |
|---:|---:|---:|
|3|4|4|
|5|6|6|
|9|4|6|
|13|3|6|
|15|4|2|
|17|3|0|
|23|3|0|
|25|3|0|
|27|1|2|
|29|0|6|
|33|1|0|

Both first-moment sums are 432. The second-moment sums are 8064 for U and
8640 for V. Under the stated probability law the means are both 27/2,
whereas the variances are 279/4 and 351/4. Even the exact support and its
maximum differ: U reaches 33 and does not reach 29, while V reaches 29 and
does not reach 33. The equal ordinary H and equal automorphism orders do not
repair this loss.

**Sharpness.** The audit exhausts every labelled tournament at orders
0,1,2,3,4, and every attachment vector, using literal permutation path
counts throughout. At fixed order, any two sources with equal H have equal
attachment histograms. Therefore order five is the first failure, including
for the stronger quotient retaining the full cover profile, since that
profile includes H. The only repeated-H pair of nonisomorphic order-four
classes is the converse pair 1->C3 and C3->1; converse preserves the entire
attachment histogram by complementing the attachment vector. This explains
the small boundary instead of treating it as an unexplained numerical match.

## 3. Why the cover profiles agree and why the next distribution does not

The equality is already predicted by
[THM-3166, falling-factorial order-join transform](../../../01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md).
If p_T(k) counts **unordered** spanning path covers and
Q_T(t)=sum_k p_T(k)(t)_k, that theorem proves

```
Q_(A->B)(t)=Q_A(t)Q_B(t).
```

Multiplication is commutative, so this complete profile forgets SCC order.
U and V have the same SCC multiset {1,1,C3}. Ordered cover counts are k!
p_T(k), hence they also agree. The present counterexample identifies a
specific observable that does not descend through that inherited quotient:
the law, and even the support, of H after arbitrary one-vertex attachment.
No novelty claim is made for the order-join transform or its loss of order.

The endpoint sidecar explains the failure exactly. Let s_i and e_i count
Hamilton paths starting and ending at vertex i. Let C_ij count ordered
spanning two-path covers (P,Q) with last(P)=i and first(Q)=j. Then i!=j and
for every attachment vector b,

```
H(T+x)
 = sum_i [s_i(1-b_i)+e_i b_i]
   + sum_(i!=j) C_ij b_i(1-b_j).                         (1)
```

**Proof.** In a child Hamilton path, x is either first, last, or internal.
Deleting x in the first two cases gives a source Hamilton path with one
specified endpoint. Deleting internal x gives a unique ordered two-path
cover with the two displayed boundary endpoints. Conversely each such
source object produces a child path precisely under the indicated one or
two attachment-bit conditions. These three cases are disjoint. QED.

For independent fair attachment bits,

```
E[H(T+x)] = H(T)+P_2(T)/4,                               (2)
```

where P_2=sum_ij C_ij is the ordered two-cover count. Indeed every endpoint
condition has probability 1/2 and each two-endpoint condition has probability
1/4. U and V therefore have the same mean 3+42/4=27/2.

Squaring (1) instead introduces pairs of conditions. Their joint probability
depends on whether they mention the same named vertex, impose the same bit,
or conflict. The totals H and P_k forget this joint incidence. A pair of
requirements b_i=1 and b_i=0 is impossible, while requirements on distinct
vertices can coexist; replacing all incidences by scalar cover counts
changes the second moment.

Here the start vectors agree, but the end vectors already differ:

```
s_U=s_V=[3,0,0,0,0],
e_U=[0,0,1,1,1],        e_V=[0,0,0,0,3].
```

The checker retains the full C matrices, verifies (1) at every attachment,
and prints them. These are sufficient sidecars for this one-vertex response.
No all-order minimum carrier claim follows, and the Session146 order-ten
kernel correction is unaffected.

## 4. Audit of producer scope and reproduction

The producer's canonicalizer searches every permutation within equal-score
cells, and orders those cells by score. This is exact: an isomorphism must
preserve score, so every possible isomorphism is represented. Its reported
stabilizer multiplicity is the number of orders giving the canonical image,
equal to the automorphism order. The class-orbit masses independently agree
with n!. Its order-0..6 F_n and self-converse polynomials agree with the
separate augmentation census in `spectrum_probe.py`. Neither census uses the
other's representatives.

The sharp failure supports a precise negative conclusion: the scalar H,
automorphism order, and even full spanning-cover profile do not determine
the arbitrary-attachment H distribution. It does **not** refute the exact
SCC convolution, a recurrence supplied with named endpoint data, or the
inherited scalar A000568 recurrence with its explicit coefficients.

Run the independent audit with:

```sh
python3 -I -B audit_growth.py > /tmp/growth-audit.out
python3 -I -B -O audit_growth.py > /tmp/growth-audit-optimized.out
cmp /tmp/growth-audit.out /tmp/growth-audit-optimized.out
```

Both runs produced identical output. All checks use explicit exceptions.
No source filters are imposed in the smaller-order sharpness universe;
all joins, named vertices, ports, source weights and attachment bits are
specified explicitly.
