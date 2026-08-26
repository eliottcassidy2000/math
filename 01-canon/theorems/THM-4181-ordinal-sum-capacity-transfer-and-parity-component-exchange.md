---
id: THM-4181
title: "Ordinal-sum capacity transfer and parity component exchange"
status: >
  PROVED exact ordinal-sum capacity factorization, rooted two-path-cover
  parity identities, even-order component-exchange bijection, and arbitrary
  block-tensor G_+ formula + FINITE-EXACT strict positive ordinal remainder
  on 242,060 ordered block-class presentations with A of orders one through
  seven and no-sink B of orders three through seven +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. The all-order ordinal remainder,
  the order-at-least-four no-sink/no-source sign law, the strong residual,
  and the order-eleven asymmetric bank remain OPEN.
source: codex-frontier-synthesis-creative-20260826az
depends_on:
  - THM-4177-root-split-johnson-current-and-source-sink-boundary-census
related:
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4167-tournament-exposure-capacity-deletion-support-moment-and-parity-holonomy
  - THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer
  - THM-4172-multideletion-support-tomography-and-same-parity-johnson-holonomy
script: 04-computation/tournament_ordinal_sum_capacity_transfer_thm4181.py
output: 05-knowledge/results/tournament_ordinal_sum_capacity_transfer_thm4181.out
independent_audit_script: 04-computation/tournament_ordinal_sum_capacity_transfer_thm4181_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_ordinal_sum_capacity_transfer_thm4181_independent_audit.out
script_sha256: 735f1cadb699cb44b190a84ab7284fe8aedc3ecc68109a416302926a29c36099
output_sha256: b4bd9d46d2964566f6326d52caf3648237713876779e50c9a469c5fb589f5922
independent_audit_script_sha256: 2be285497e4e930279c97df904874753f38a266bbdd0fbaa5941e13fa46c2089
independent_audit_output_sha256: a6f97f326d159b288ea398b7674605497ffe4d16deadbd9ea9d3abdc387ec5f3
gentourng_sha256: 89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone endpoint/subset-convolution implementation checks the
  proved sidecar and layer identities on all 532 tournament-class
  representatives through order seven, directly reconstructs 24 small
  ordinal capacities, and scans all 242,060 declared ordered block-class
  presentations with the right factor restricted to no-sink orders three
  through seven. Normal, optimized, and fixed-hash streams byte-match.
independent_audit: >
  ACCEPT. A separate C++ referee uses literal permutation and labelled-
  tournament enumeration, not gentourng or subset DP. It checks 33,867
  labelled sidecars through order six, 65,600 component exchanges, 374
  literal ordinal transfers, 55,350 labelled remainder controls, and
  90,166 labelled dominance presentations. Clang O0/O3 and ASan/UBSan
  streams byte-match.
---

# THM-4181 -- ordinal-sum capacity transfer and parity component exchange

**PROVED exact transfer and component exchange + FINITE-EXACT ordinal
remainder census + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4177 reduces a minimum counterexample to the no-sink sign law to a
source-and-sink-free tournament, but it does not reduce it to a strong
tournament. This theorem computes exactly what happens at a nontrivial
strong-component cut. The cross transfer has rank at most two, and genuine
rank-two examples occur: contracting
a block while retaining only its Hamilton count or its capacity tensor loses
two rooted path-cover parity states. Those states are the required sidecar.

The resulting ordinal remainder is strictly positive in a large exact
finite universe, but its all-order sign is not proved. Thus the strong-core
reduction in Section 7 is **CONDITIONAL**, and neither THM-4177's open sign
law nor THM-4169's order-eleven asymmetric bank is closed here.

## 1. Conventions and object types

For nonempty tournaments `A,B`, let `T=A triangleright B` be their ordinal
sum: every vertex of `A` dominates every vertex of `B`. Use `H(empty)=1`,
and let `c_ij(X)` be THM-4177's integer exposure capacity of the **actual
tournament** `X`. Put `W(X)=sum_(i<j)c_ij(X)`.

Because both factors are nonempty, `A triangleright B` has a sink if and
only if `B` has a sink; dually, it has a source if and only if `A` has a
source. Thus Section 6's no-sink filter on `B` is exactly the intrinsic
no-sink sector of the ordinal-sum child, not an auxiliary presentation
condition.

For a vertex `a` of a tournament `X` and `epsilon in {0,1}`, define

```text
U_a^epsilon(X)
 =sum_(P directed simple path in X starting at a,
        |V(P)| == epsilon mod 2) H(X-V(P)),               (1)

V_a^epsilon(X)
 =sum_(P directed simple path in X ending at a,
        |V(P)| == epsilon mod 2) H(X-V(P)).               (2)
```

Directed simple paths in `(1)--(2)` include the one-vertex, zero-arc path.
Those singleton terms realize the direct cross-arc contribution in `(4)`.

The letters `U,V` record the distinguished path's start or end. They are
not capacity coordinates and are not retained by contracting the block to
one quotient vertex.

No deletion card is used below. For example, `c(T)|_A` means the literal
restriction of the actual parent tensor `c(T)` to edges inside `A`, while
`c(A)` is the actual capacity tensor of the factor. Equation `(3)` relates
these distinct objects. The finite census counts ordered block-class
**presentations** `(A,B)`. It is not a labelled-child count, a rooted-orbit
count, or an unrooted child-isomorphism count.

## 2. Exact ordinal capacity transfer

> **Theorem 1 (rank-two ordinal transfer).** For `T=A triangleright B`,
>
> ```text
> c_aa'(T)=H(B)c_aa'(A),                         a,a' in A,
> c_bb'(T)=H(A)c_bb'(B),                         b,b' in B,       (3)
>
> c_ab(T)=2[U_a^0(A)V_b^0(B)+U_a^1(A)V_b^1(B)], a in A,b in B.  (4)
> ```
>
> In particular the `A x B` capacity matrix has ordinary rank at most two.

### Proof

Use first the exposed-word interpretation of capacity. If the marked gap
lies inside `A`, an exposed word cannot contain an unmarked transition from
`B` back to `A`; the marked endpoints are both in `A`, so such a transition
could not be hidden in the marked gap. Therefore all `A` vertices precede
all `B` vertices. Restriction gives an exposed `A` word and an arbitrary
Hamilton path of `B`, and concatenation is the inverse. This proves the
first identity in `(3)`; the second is identical.

For `(4)`, use THM-4177's odd-path capacity formula. A directed simple path
between `a in A` and `b in B` can only run from `a` to `b`. It crosses the
ordinal cut exactly once and decomposes uniquely as `P_A` followed by
`P_B`, where `P_A` starts at `a` and `P_B` ends at `b`. Its arc length is
`|V(P_A)|+|V(P_B)|-1`, which is odd exactly when the two vertex counts have
the same parity. Moreover

```text
H(T-V(P_A)-V(P_B))
 =H(A-V(P_A))H(B-V(P_B)),                         (5)
```

because Hamilton paths of an ordinal sum factor, including an empty factor.
Substitution into the odd-path formula proves `(4)`. QED.

The bound is attained. For `A=B` equal to the transitive two-vertex
tournament, the cross matrix is

```text
[[2,4],
 [2,2]],
```

whose determinant is `-4`.

The source, target, and loss of the quotient map are now explicit:

```text
source:       the actual factor tournaments A,B,
target:       the actual tensor c(A triangleright B),
preserved:    internal capacities, after H-scaling,
lost by bare contraction: rooted parity states U and V,
needed sidecar: the two coordinates in (1)--(2).
```

For a general strong quotient substitution, a Hamilton path may enter the
same module in several sections. Then the two-state sidecar is no longer
complete: one needs section-counted path-cover data. Thus `(4)` is an exact
ordinal-cut transfer, not a claimed transfer theorem for arbitrary modules.

## 3. Path-cover mass and component exchange

Write `U^epsilon(X)=sum_a U_a^epsilon(X)` and similarly for `V`.

> **Theorem 2 (two-path-cover identities).** Every nonempty tournament
> satisfies
>
> ```text
> U^0+U^1=V^0+V^1=W+H.                              (6)
> ```
>
> If `|X|` is odd, then
>
> ```text
> U^1-U^0=V^1-V^0=H.                                (7)
> ```

### Proof

An object counted by `U^0+U^1` is an ordered pair `(P,Q)` in which `P` is a
nonempty directed path and `Q` is a Hamilton path of the complementary
vertex set. If `Q` is nonempty, concatenate `P|Q` and mark the intervening
gap. This is an exposed word, and splitting an exposed word at its marked
gap is the inverse. Those objects contribute `W`; the objects with
`Q=empty` are exactly the `H` spanning paths. This proves the first equality
in `(6)`. Recording the end of `P` instead of its start proves the second.

When `|X|` is odd, component exchange `(P,Q) <-> (Q,P)` reverses the parity
of the distinguished component on every nonspanning object. Those terms
cancel in `U^1-U^0`; an unmatched spanning `P` has odd vertex count and
contributes `H`. The endpoint version is identical. QED.

There is a different, capacity-layer form when the tournament order is even.
Recall THM-4177's odd length layer `c^[ell]`, where `ell` is the number of
arcs of the distinguished path.

> **Theorem 3 (even-order component exchange).** Let `n=|X|` be even and
> let `ell` be odd with `1<=ell<=n-3`. Expand each capacity atom as
>
> ```text
> (R,S),
> R = a directed ell-arc simple path,
> S = a chosen Hamilton path of X-V(R),                    (8)
> ```
>
> retaining the factor-two weight from the odd-path formula. Then
> `(R,S) |-> (S,R)` is a weight-preserving involutive bijection from
> expanded layer-`ell` atoms to expanded layer-`n-ell-2` atoms. The two
> endpoint chords are vertex-disjoint.

### Proof

The complement path has `n-ell-1` vertices and hence `n-ell-2` arcs. That
number is odd because `n` is even and `ell` is odd. The upper bound on
`ell` ensures that both components have at least two vertices, so both have
endpoint chords. After the swap, `R` is a Hamilton path of the complement
of `S`; swapping again recovers the original object. Both expanded atoms
have the same factor-two weight, and their vertex-disjoint components give
disjoint endpoint chords. QED.

Consequently

```text
sum_(i<j)c^[ell]_ij
 =sum_(i<j)c^[n-ell-2]_ij.                              (8a)
```

This is a total layer-mass identity. The exchange changes endpoint chords,
so `(8a)` is not a coordinatewise tensor equality; the primary audit's 116
checks are exactly mass checks.

The boundary `ell=n-1` is deliberately excluded: the complement is empty
and supplies no chord to exchange. This is exactly the boundary on which
THM-4177 found the fifteen negative order-eight self-layer carriers. The
bijection does not by itself prove that all cross-layer polarizations are
nonnegative.

## 4. Exact block-gate formula

This subsection is valid for arbitrary symmetric tensors; only the
specialization through `(3)--(4)` is capacity-specific. On the oriented
complete graph `A triangleright B`, let `x,y,z` be supported respectively on
edges inside `A`, edges inside `B`, and cross edges `a->b`. Put

```text
p_a=sum_b z_ab,                 q_b=sum_a z_ab,
Z=sum_a p_a=sum_b q_b,          Z_2=sum_(a,b)z_ab^2.      (9)
```

For `x`, let `W_x` be its total edge mass, `d_a^x` its unsigned degree, and
`o_a^x=sum_(a->a')x_aa'`. For `y`, define `W_y,d_b^y` and
`r_b^y=sum_(u->b)y_ub`.

> **Theorem 4 (generic ordinal block formula).**
>
> ```text
> G_+(x+y+z)
>  =G_+(x)+G_+(y)+W_xW_y
>   +sum_a p_a(W_x-d_a^x+4o_a^x)
>   +sum_b q_b(W_y-d_b^y-4r_b^y)
>   +1/2[Z^2+Z_2+3sum_a p_a^2-5sum_b q_b^2].             (10)
> ```

### Proof

Apply THM-4177's polarization kernel. An `x` edge and a `y` edge are always
disjoint, giving `W_xW_y`. For a cross edge `a->b`, the `x` edges avoiding
`a` contribute `W_x-d_a^x`; the edges leaving `a` have a common tail and
coefficient `+4`; the edges entering `a` are mixed and contribute zero.
This gives the first linear sum. Dually, the `y` edges avoiding `b`
contribute `W_y-d_b^y`, while the edges entering `b` have a common head and
coefficient `-4`, giving the second sum.

For the cross tensor alone,

```text
D(z)=1/2[Z^2-sum_a p_a^2-sum_b q_b^2+Z_2],
C(z)=sum_a p_a^2-sum_b q_b^2.                            (11)
```

Thus `G_+(z)=D(z)+2C(z)` is the last line of `(10)`. QED.

For the actual ordinal capacity, take

```text
x=H(B)c(A),
y=H(A)c(B),
z_ab=2[U_a^0V_b^0+U_a^1V_b^1].                          (12)
```

Since `G_+` is homogeneous quadratic,

```text
G_+(x)=H(B)^2G_+(c(A)),
G_+(y)=H(A)^2G_+(c(B)).                                  (13)
```

Define the exact ordinal remainder

```text
R_+(A,B)
 =G_+(c(A triangleright B))
   -H(B)^2G_+(c(A))-H(A)^2G_+(c(B)).                     (14)
```

Equations `(10)--(13)` are a proved finite-dimensional formula for `(14)`.
They do **not** prove its sign.

## 5. Why the pointwise cancellation fails

The negative kernel terms suggest two coordinate dominances:

```text
(L) c_uv(T)>=c_vb(T)  for u->v in A and b in B,
(R) c_au(T)>=c_uv(T)  for u->v in B and a in A.           (15)
```

Both implications are false. Use the lexicographic pair-bit convention
`(0,1),(0,2),...,(0,n-1),(1,2),...,(n-2,n-1)`.

In the complete declared representative universe used in Section 6,

```text
1<=|A|<=7,  3<=|B|<=7,  B has no sink,                  (16)
```

ordered by child order, block orders, labels, and coordinates, the first
`(L)` failure already occurs at

```text
A=B=C3, labels 101,101:
c_01(T)=6 < c_1,3(T)=10.                                 (17)
```

The lexicographically first `(R)` failure occurs one child order later:

```text
A=one vertex, label empty,
B=111111100111111 (order six, no sink):
c_0,1(T)=22 < c_1,6(T)=28.                               (18)
```

Within that first failing block order `(|A|,|B|)=(1,6)`, the sharp slack is
`-8`; with label-first tie-breaking a witness is

```text
B=111111101111110: 18 < 26.                              (19)
```

If `|A|>=2` is imposed, the first failure has child order eight with
`A=1` (the transitive two-vertex tournament) and the same order-six factor
as `(18)`. If `|A|>=3` is imposed, the first failing block order is
`(3,5)`: `A=C3`, `B=1111100111`; its lexicographically first negative row
is `20<24`, while the same presentation contains the sharper row `20<30`.
Thus the global singleton witness and the strong-left-block stress control
are distinct statements.

These first-failure statements are relative to the exhaustively declared
universe and ordering `(16)`. Five-by-five stress controls are

```text
(L): A=B=1100111111,               88-730=-642,
(R): A=1110101111, B=1111100111,   90-150= -60.          (19a)
```

Thus a proof of `R_+>0` must aggregate coordinates or retain additional
path-pair multiplicity. The failed inequalities cannot be promoted into an
edgewise injection.

## 6. Finite exact ordinal remainder census

**FINITE-EXACT.** The primary universe contains one `gentourng`
representative of every tournament class `A` of orders one through seven
and every no-sink tournament class `B` of orders three through seven:

```text
A class counts:       1,1,2,4,12,56,456,       total 532,
B no-sink counts:         1,2,8,44,400,         total 455.
```

Hence there are exactly `532*455=242,060` ordered block-class
presentations. On all of them,

```text
R_+(A,B)>0.                                             (20)
```

There are no zero rows. The exact minimum is
`R_+(one vertex,C3)=120`. Again, this is not a child-orbit count. The same
unrooted child may admit different strong-component cuts, and the census
intentionally retains the ordered factors.

The primary engine also checks all 532 sidecar identities, 31,543 odd-path
layer coordinates, 116 even component-exchange masses, and 24 direct small
transfers. Its presentation-stream digest is

```text
1726ec4ca437f4417e6d0696f70303eee088d765aeeb6627e8e72e9b08b0c55a.
```

The independent literal-permutation referee does not use `gentourng` or the
primary subset DP. It checks all 33,867 labelled tournaments through order
six, 65,600 even component exchanges, 374 literal direct transfers, and
55,350 labelled ordinal-remainder controls through `|A|<=4,|B|<=5`; its
minimum is again `120`. It separately derives the stated first dominance
failures and first-block sharp slacks from 90,166 labelled ordinal-sum
presentations of total order at most eight.

These exact finite facts are evidence for, not a proof of, the all-order
sign of `(14)`.

## 7. Conditional strong-core reduction and open boundary

Consider the explicit hypothesis

```text
(OS+)  R_+(A,B)>0 for every nonempty A and every no-sink B
       with |B|>=3.                                      (21)
```

> **Conditional corollary.** If `(OS+)` holds, a minimum-order
> counterexample to THM-4177's no-sink `G_+>0` conjecture is strong.

Indeed THM-4177 already makes a minimum counterexample source-free and
sink-free. If it is not strong, its strongly connected components are
linearly ordered. Let `A` be the first component and let `B` be the union
of the remaining components. Source-freeness forces `|A|>=3`, and
sink-freeness forces the terminal component, hence `B`, to have no sink and
order at least three. Both factors are smaller no-sink tournaments.
Minimality gives nonnegative factor gates, with the directed three-cycle
handled by its exact zero gate. Equations `(13)--(14)` and `(OS+)` then give
`G_+(T)>0`, a contradiction. Taking converses gives the analogous
conditional reduction for `G_-`.

But `(OS+)` is **OPEN**. Therefore all of the following remain open:

```text
|T|>=4 and no sink   => G_+(c(T))>0,
|T|>=4 and no source => G_-(c(T))>0,
the all-strong residual of those implications,
the full order-eleven rational asymmetric bank,
every exact Johnson-coset and actual-response consequence.               (22)
```

The even component exchange does not cover spanning atoms, the two
coordinate dominances are false, and finite positivity through factor order
seven supplies no multiplicity-preserving all-order injection. This is the
exact remaining obstruction.

## 8. Replay

Primary exact audit:

```bash
python3 -B \
  04-computation/tournament_ordinal_sum_capacity_transfer_thm4181.py
python3 -O -B \
  04-computation/tournament_ordinal_sum_capacity_transfer_thm4181.py
PYTHONHASHSEED=4181 python3 -B \
  04-computation/tournament_ordinal_sum_capacity_transfer_thm4181.py
```

Independent literal audit:

```bash
clang++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_ordinal_sum_capacity_transfer_thm4181_independent_audit.cpp \
  -o /tmp/thm4181_independent_O0
/tmp/thm4181_independent_O0

clang++ -std=c++20 -O3 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_ordinal_sum_capacity_transfer_thm4181_independent_audit.cpp \
  -o /tmp/thm4181_independent_O3
/tmp/thm4181_independent_O3

clang++ -std=c++20 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_ordinal_sum_capacity_transfer_thm4181_independent_audit.cpp \
  -o /tmp/thm4181_independent_san
/tmp/thm4181_independent_san
```

The three Python streams byte-match the frozen primary output. Clang O0/O3
and ASan/UBSan byte-match the frozen independent output. **QED for the exact
transfer, path-cover identities, component exchange, generic block formula,
and the stated finite universes only.**
