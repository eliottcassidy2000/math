# Decode the total before choosing the tournament

This note extends the [first triangle terrace factorization](https://github.com/eliottcassidy2000/math/blob/8051b02275faab24bffba7fbed5cc8d3afc29bac/05-knowledge/results/chanpak_target_first_terrace_density_20261114/TERRACE_DENSITY.md)
and the [two-wing corridor](https://github.com/eliottcassidy2000/math/blob/26519368143bbd88e824d8aba2c86ef2c0adffb0/05-knowledge/results/session337_rooted_line_surfaces_20261114/ZERO_CORRIDOR.md).
A terrace word records the nonanchor gaps along a backward chain
in descending score order. Its positive widths are (a_1,ldots,a_k),
and the numbers of intervening zero gaps (q_1,ldots,q_{k-1}).
Normalized words begin and end positive. They represent actual
isomorphism classes; reversing the word gives the converse.
Write (D_y=\operatorname{diag}(1,y)), (y_i=2^{a_i}), and
(C_q=\bigl(\begin{smallmatrix}e_q&f_q\\f_q&g_q\end{smallmatrix}\bigr))
for the shared zero-corridor kernel. The exact ordinary count is
[
 (1,1)D_{y_1}C_{q_1}D_{y_2}\cdots C_{q_{k-1}}D_{y_k}(1,1)^t.
]
The kernel counts zero choices conditioned on the two neighboring
positive selections, so the expression retains the whole count.

## Exactly when a corridor forgets a degree of freedom

The dream research observation is
[
 \det C_q=-T_{-q-4},
]
where the tribonacci sequence on all integer indices has
(T_0=0,T_1=T_2=1). To prove it, the entries of (C_q)
are bilinear pairings of two fixed input vectors and two fixed output
functionals across a common three-dimensional tribonacci transfer.
The determinant is therefore a coefficient of its second exterior
power. Concretely, use the Session 337 map
[
 M(E,O,V)=(E+2V,E+V,O),
]
with inputs ((1,0,0)), ((0,1/2,1/2)), and output
functionals (E+2V), (E+O+V). Then
(C_q=L M^q U), using those two functionals as rows
and inputs as columns (the offdiagonal entries agree).
The characteristic polynomial of (M) is (x^3-x^2-x-1)
and its determinant is 1. The characteristic polynomial of
the exterior square is (x^3+x^2+x-1): its eigenvalues are the
reciprocals of those of (M). The determinants thus follow
(d_{q+3}=d_q-d_{q+1}-d_{q+2}) from (0,-2,3),
as do the displayed negative-index terms.

[Mignotte and Tzanakis](https://doi.org/10.4064/aa-57-4-357-364)
prove that the only integral zeros of this sequence are
(0,-1,-4,-17), in this indexing; see also
[Bilu et al.](https://arxiv.org/pdf/2401.06537).
Thus the only reset lengths are **0 and 13**, not merely the only
ones seen in a search. The rank-one factorizations are
[
 C_0=(1,1)^t(1,1),\qquad C_{13}=(47,56)^t(47,56).
]
Both give a genuine factorization of the ordinary path count at
that seam, with positive odd factors greater than one.
In particular no prime-count terrace witness can contain either.
For 13, a two-wing count is
((47+56\cdot2^a)(47+56\cdot2^b)); an isolated positive
gap flanked by two such seams contributes (2209+3136\cdot2^a).

Every finite modulus sees false resets: the tribonacci triple map
is invertible modulo any modulus (m), so there is a period (p>0)
for the simultaneous coefficient triples of (C_q) modulo (m).
Take (q=13+jp) with (j>0). Then (C_q\equiv C_{13}\pmod m),
but (C_q) is invertible over \(\mathbb Q\). The same works for
any finite list of moduli by taking their least common multiple.
This distinguishes a genuine determinant zero from any finite
congruence certificate that only reads the coefficient residues.

For an odd prime (p), if (p\nmid T_{-q-4}), the seam cannot
turn a nonzero boundary vector into zero modulo (p).
If (p\mid T_{-q-4}), use the actual input vector to test
its membership in the kernel of (C_q\bmod p).
A common (p) divisor of the two state entries can never disappear
under subsequent integral transfers and forces (p\mid H).
In particular, a prefix collapsing modulo (p<N) can be discarded
for a prime target (N). A singular matrix alone does not imply
that a particular prefix collapses.

## Decode a frame with its outer widths unspecified

Fix the zero seams, the interior positive widths, and the two
boundary rows (u,v\in\{(1,1),(47,56)\}). The frame readout is
[
 N=A+B x+C y+Dxy,\qquad x=2^a,\quad y=2^b,\quad a,b\ge1.
]
Here (A,B,C,D>0), all at most (D). To see this, each kernel
has nondecreasing rows and columns in the two boundary indices
((e_q\le f_q\le g_q)). The first and last kernels preserve
these inequalities under multiplication by nonnegative matrices
in the relevant directions, and the boundary flags satisfy
(1\le1) and (47<56). No common interior orientation is assumed.
Consequently
[
 D2^{a+b}\le N\le D(1+x)(1+y)\le(9/4)D2^{a+b}.
]
Test only (s=\lfloor\log_2(N/D)\rfloor) and (s-1).
For each, test the integer quadratic
[
 Bx^2-(N-A-D2^s)x+C2^s=0
]
for roots (x=2^a\) with (1\le a<s), then verify the actual
readout. The factor view of the same fiber is
[
 (Dx+C)(Dy+B)=DN+BC-AD.
]
The correction (AD-BC) is the signed product of the seam
determinants and the interior widths' powers of two, times the
four positive boundary entries. A 13-zero reset therefore erases
exactly the correction; a nonreset preserves a nonzero one.

For two wings with both flags ((1,1)), iterating (q) while
(g_q\le N/4) gives a full inverse in (O(\log N)) arithmetic
tests, one log factor fewer than the previous scan. For arbitrary
length, finding the interior frame itself is the expensive part.
The result promises two exponent sums per fixed frame, not
(O(\log N)) for the whole terrace.

## Select a canonical isomorphism-class witness on the whole terrace

Build from the exact target, not a list of tournaments. Cut at all
reset seams; the adjacent components share the boundary flag
corresponding to that reset. A reset-free component may be one
positive gap (decode with one divisibility/power-of-two test)
or a nonreset chain (use the two-exponent inverse on each interior
frame). Enumerate frames with outer widths one in the pruning
bound. All kernel entries are nondecreasing with (q) and eventually
grow; increasing an interior width increases the bound; appending
one more positive site and corridor strictly increases the minimal
readout. Hence all three loops terminate without imposing an
arbitrary size cutoff. Equivalently, the compound lower bound
in the previous note supplies (O(\log N)) total width and
zero budget, even before cutting.

For a prescribed (N), use states ((N,l,r)) where (l,r)
are the two boundary flags. Either choose one reset-free component
of value (N), or choose a proper divisor (d) as the leftmost
component and append the state ((N/d,t,r)) at a reset of
type (t) between them. Both factors are odd and at least three:
the first-state coordinate after each kernel is odd, and its
second coordinate is made even by the next positive width, so
both allowed outer flag vectors preserve oddness. They have
positive entries and cannot reduce a readout. The value strictly
decreases at each recursion. Compare witnesses by number of
vertices, then lexicographically by the normalized gap word.
Additivity of the gap cost under concatenation makes it safe to
retain just the best witness in every state. The final flags are
both ((1,1)). This produces a deterministic first-terrace class
(including a deterministic choice between converse classes if they
tie), or a proof that no terrace witness has the desired count.

Example: 573 is first decoded there at order 12, with one primitive
representative ((1,0,0,3,0,1)). 385 has no terrace witness even
though the repo has two-hub witnesses. 25281 decodes at order 18
using the 13-zero reset. See `python decoder.py 385 573 25281`
and `python replay_decoder.py`.

**Boundary:** the decrease above is valid once one knows a legal
seam/frame for the target. No arbitrary odd input is guaranteed
such a choice. A global least-order or lexicographic rule selects
uniquely *conditional on a witness*; it cannot by itself prove
existence. A sufficient proof route would use several verified
endpoint states and prove that every target outside a finite seed
set has at least one physically realizable decreasing move. The
13-zero reset and the determinant test show exactly which steps
cannot possibly be used for primes in this family.

## Exact independence versus fast loss of information

Give every positive-gap bit its Boltzmann weight (2^{a_i}) if
selected and reweight every neighboring pair by (C_{q_i}).
Normalize the resulting weights by the actual count (H). When
only two positive-gap bits (b_i,b_j) are observed, their 2-by-2
joint frequency table is, up to normalization and strictly positive
row/column scales, the matrix product of the intervening kernels
and interior diagonal weights. Its determinant has the sign of
(prod_{i\le h<j}\det C_{q_h}). For binary variables this is
also the sign of the covariance. Therefore their exact statistical
independence occurs if and only if a 0- or 13-zero reset separates
them. The rank-one decomposition is the corresponding independence
of the whole left and right choices, not just of one pair.

There is nevertheless uniform *approximate* memory loss at every
seam. If (t>0) is the ratio of the two coordinates entering a
corridor, the outgoing ratio after the next positive width is
(y(f+gt)/(e+ft)). Its derivative with respect to (log t),
measured in outgoing log-ratio, is
[
 \frac{(eg-f^2)t}{(f+gt)(e+ft)}.
]
The maximum absolute value over (t>0) is
[
 \frac{|\sqrt{eg}-f|}{\sqrt{eg}+f}.
]
This follows by minimizing the denominator divided by (t),
whose variable terms are (ef/t+fg t).
The boundary inequalities (e\le f\le g\le2e) hold from the
first three kernel triples by the positive tribonacci recurrence.
Consequently (eg/f^2\in[1/2,2]), and the log-ratio contraction
is at most (3-2\sqrt2<0.172), uniformly even in the positive
width (y), which merely adds a log-ratio translation. The
bound is attained as a derivative norm at (q=1).

Unlike exact factorization, this contraction does not by itself
select integer targets. It explains why interval computations
of word pressure can be sharp: conditional projective states
quickly forget remote boundaries, even though exact projective
state differences can be arithmetically decisive.
