# One fixed center sees a quadratic seam that the shell forgets

**Status: PROVED exact one-center law; finite-exact separation.** This
does not resolve Chan–Pak. It exhibits a strictly stronger
strong-self-converse observation than shell iteration on a specific
pair of equal-profile hosts.

## Two-path kernel

Let U be a nonempty even-order self-converse tournament with a specified
anti-involution `sigma`. Add a vertex z fixed by sigma, and from each
sigma-pair choose exactly one vertex that **z beats**. Equivalently,
write `y_v=1` if `v → z`, with `y_sigma(v)=1−y_v`. The extension is
self-converse; if U is strong, the extension is strong because z has
a nonempty in-neighborhood and a nonempty out-neighborhood in U.

Let `E_i` count Hamiltonian paths of U ending at i. Let `A_ij` count
ordered pairs of nonempty disjoint directed paths covering U, the
first ending at i and the second starting at j. Set
`B_ij=A_(i,sigma(j))`. Then B is symmetric, and

```text
H(U+z) = 2 sum_i E_i y_i + sum_(i,j) B_ij y_i y_j.
```

To prove this, cut a Hamiltonian path at z. If z is last, the
contribution is `sum E_i y_i`; the anti-involution identifies this
with the z-first contribution. Otherwise both sides are nonempty,
and a pair with endpoints i then j contributes exactly when
`y_i(1−y_j)=1`. Replace j by sigma(j) in the sum. To see symmetry of
B, reflect and reverse both old paths and exchange their order. The
first now ends at j, the second begins at sigma(i), and the operation
is an involution.

Fix a reference choice for every pair and write `t_i ∈ {0,1}` for the
decision to flip pair i from that reference. Substituting
`y_v ∈ {t_i,1−t_i}` above gives a multilinear polynomial of degree at
most two. Thus the center counts are recovered from
`1+m+binom(m,2)` queries (empty flip set, singletons, pairs). The
off-diagonal interaction matrix with entries

```text
M_ij = H({i,j}) − H({i}) − H({j}) + H(empty),  M_ii=0
```

has invariant rank: a change of reference conjugates it by a diagonal
sign matrix and a permutation of pairs permutes its rows and columns.
In fact if `y(t)=y(0)+sum_i t_i w_i`, where each `w_i` is a signed
difference of matched unit vectors, then `M_ij=2 w_i^T B w_j` for
`i≠j`. This is the off-diagonal part of the two-path kernel on the
anti-invariant degree of freedom, not a general profile reconstruction.

## The orientation-neutral first moment

Choose the center's signs independently and fairly on the mate pairs.
Then $\mathbb E[y_i y_j]$ is 1/4 across pairs, 1/2 if i=j, and zero
if j=sigma(i). The latter entries of B already vanish, since
$A_{ii}=0$. Also $\sum_{i,j} B_{ij}=Z_U(2)-2H(U)$, as this counts
all ordered two-path covers with both paths nonempty. Consequently

```text
E[H(U+z)] = H(U)/2 + Z_U(2)/4 + tr(B)/4.
```

Thus the average center observation carries exactly one extra scalar,
the reflected endpoint trace, once the slot profile is known. Each
$B_{ii}=A_{i,\sigma(i)}$ counts spanning paths using the mate edge
$i\to\sigma(i)$ in the orientation of that edge consistent with
the direction shown: glue the two paths across that edge. Its value
depends only on the other arcs of the host, since the split cover
does not use the mate edge. This is a direct endpoint interpretation
of the missing statistic.

The same argument is graded. Let $R_U(r)$ count assignments to two
ordered nonempty distinguished path slots, with first endpoint $i$
and second start $\sigma(i)$ summed over $i$, and $r-1$ additional
labelled path slots (which may be empty). Then for positive $r$,

```text
E[Z_(U+z)(r)] = r [ Z_U(r)/2 + (Z_U(r+1)+Z_U(r-1)+R_U(r))/4 ].
```

Choose the color of z (r possibilities), then split that slot at z.
If its two old segments are empty, their total old weight is
`Z_U(r−1)`. If exactly one segment is nonempty, either side is
permitted with probability 1/2, producing `Z_U(r)−Z_U(r−1)`.
If both are nonempty, use allowance 1/4 except at matching ends,
where it is 1/2. The unrestricted two-nonempty weight is
`Z_U(r+1)−2Z_U(r)+Z_U(r−1)`, with additional quarter weight for
`R_U(r)`. This sums to the displayed formula. Thus full
center-and-shell observations obtain an entire reflected-endpoint
polynomial beyond the principal specialization, with `R_U(1)=tr(B)`.

## The two six-vertex hosts

Use outgoing masks, where bit j is set exactly when the vertex beats j:

```text
U=(58,4,1,6,14,30)     sigma_U=(2,5,0,4,3,1)
V=(38,8,2,5,15,30)     sigma_V=(3,5,4,0,2,1)
```

They are strong and both have principal specialization
`r^6+8r^4+8r^2`. The reference for z is to beat the lower-index
member of each mate pair; list those pairs by their lower index. An
exact subset-DP evaluation (or direct enumeration) gives the
degree-two polynomials

```text
U: 49 − 16 t1 + 56 t2 + 18 t3 + 4 t1 t2
V: 41 − 14 t1 + 34 t2 + 46 t3 + 8 t1 t2 − 8 t1 t3 + 8 t2 t3.
```

The interactions already distinguish the hosts: U has one nonzero
unordered pair of interactions (matrix rank 2), while V has three
(rank 3). Stronger, the two sets of eight Hamiltonian counts are
disjoint. Their ascending lists are

```text
U: 33,49,51,67,93,105,111,123
V: 27,41,65,69,75,87,115,129.
```

Thus a common symmetry-preserving rule—adjoin one fixed center and
let it beat one vertex per mate pair—distinguishes the matching
specialization fiber without relying on any alignment between the
two hosts' pairs.

Their mean observations are 79 and 76 respectively, so even the
orientation-neutral average separates them; the corresponding
reflected endpoint traces are 58 and 46.

The literal reflected-endpoint polynomials are R_U(r)=6r⁴+38r²+14
and R_V(r)=6r⁴+30r²+10. Direct enumeration verifies each at five
positive integers, enough because two nonempty distinguished slots
leave at most four old vertices for ordinary slots. The companion
checker also independently verifies the graded identity at r=1,2,3
on both hosts using directed path enumeration and subset partitioning;
normal and optimized Python agree.

Exact check: `python3 check_center_probe.py`.
