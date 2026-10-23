# The native Hamiltonian incidence has exactly two components

Status: **PROVED**, with constructive **FINITE-EXACT** controls through native
order six. This strengthens the finite component counts in
[`RESTRICTIONS.md`](RESTRICTIONS.md). It is about the unmarked incidence
coupling; source-specific pinned metagraph connectivity remains a different
obligation.

## 1. Objects and the legitimate move

Put `m=n-1`, `M=binom(m,2)`. A graph word `G` on compressed positions
`0,...,m-1` realizes an `n`-vertex tournament with Hamiltonian path
`(0,...,n-1)` by turning the native nonconsecutive arc `(i,j+1)` backward
exactly when `{i,j}` is an edge of `G`. Consecutive native arcs stay forward.

The incidence graph has ordinary compressed graph classes on one side and
native tournament classes on the other. We can therefore:

1. rearrange the compressed graph's vertex positions, retaining its source
   isomorphism class;
2. reorder the actual native vertices along another Hamiltonian path,
   retaining its native tournament class;
3. recompress along this new path, reaching another source graph class.

Only step 2 traverses a source-to-source incidence adjacency. Its entire
Hamiltonian word, including the junction to the unreordered suffix, must
be checked. Preserving the *number* of crossing backward arcs alone would
not suffice.

**Theorem.** For every `n>=3`, the native incidence graph has exactly two
components: the empty compressed graph with the transitive native target,
and all nonempty compressed graph classes with all nontransitive native
targets. Every nonempty `G` with `e` edges is at source-adjacency distance
at most `2(e-1)` from the one-edge graph class.

The singleton component was proved in `RESTRICTIONS.md`. We prove the second
component by descent in `e`, using four cases with an actual boundary guard.

## 2. A guarded strong-four reduction

Suppose a compressed vertex `u` has

\[
                    2\le d_G(u)\le m-2.
\]

Choose two neighbors `a,b` and a nonneighbor `w`. Put the compressed
vertices in positions `(a,u,b,w,...)`. Thus the first compressed triple
contains edges `01,12`, and possibly `02`, while `13` is absent.

The native first four vertices have backward free arcs `02,13`, and
possibly `03`. If the compressed triple is a path, use the native order

\[
                       (2,0,3,1,4,5,\ldots).
\]

If it is a triangle, use

\[
                       (2,3,0,1,4,5,\ldots).
\]

Both first-four orders are Hamiltonian and have precisely one backward
internal arc. Both end at old native vertex `1`. The next junction
`1 -> 4` is forward **because the compressed edge `13={u,w}` is absent**.
This is the needed source-specific boundary guard. The old suffix is then
unchanged.

The first four vertices remain before every outside vertex, so the total
number of backward arcs crossing the block is unchanged. Consequently the
compressed edge count drops by one for the path triple and two for the
triangle triple. Neither move creates the zero word, since it retains a
nontransitive native tournament.

## 3. The exceptional graph shapes also descend

If no vertex has degree between `2` and `m-2`, every degree belongs to
`{0,1,m-1}`. The graph is a matching, a spanning star, or a complete graph:
without a universal vertex it is a matching; with exactly one universal
vertex it is a spanning star; with two universal vertices every remaining
vertex has degree at least two and must also be universal. The same list
covers `m=2,3` directly.

**Matching with at least two edges.** Put two edges in compressed positions
`01` and `23`. The native backward arcs include `02` and `24`. Rotate the
first three native vertices to obtain

\[
                       (1,2,0,3,4,\ldots).
\]

The internal triple is a directed cycle; the exit `0 -> 3` is forward
because compressed `02` is absent in a matching. The edge count is
unchanged, but the new compressed graph contains adjacent edges `01` and
`13`. Its edge count is at most `m/2`, so it cannot be a spanning star or a
complete graph. In fact its vertex `1` has degree exactly two. It therefore
admits the guarded reduction in section 2 on the next move. This is a
one-step preparation, not a descent assertion for the rotation itself.

**Spanning star with at least two edges.** Put its center in compressed
position zero. The native arcs among `1,...,n-1` are all forward, and every
native vertex at least two points to zero. The Hamiltonian order

\[
                         (1,2,\ldots,n-1,0)
\]

has exactly one backward arc, namely the old `0 -> 1`. Thus the star drops
directly to the one-edge class, with no external boundary to check.

**Complete graph.** The native tournament points forward on consecutive
pairs and backward on every other pair. Use descending pairs, ascending
within each pair:

\[
 (n-2,n-1,n-4,n-3,\ldots),
\]

ending in the singleton zero if `n` is odd. For example the orders are
`(4,5,2,3,0,1)` at `n=6` and `(3,4,1,2,0)` at `n=5`. Within pairs the
arcs are consecutive and forward. Between pairs the label decreases by
three, and into the final singleton by two, so these arcs are also forward.
All nonconsecutive original pairs now point forward unless within a pair,
where they do not occur. The only backward arcs are the original
consecutive pairs separated across blocks. Their number is

\[
                  \lceil n/2\rceil-1
                  =\lfloor(n-1)/2\rfloor < M
\]

for `n>=4`. This is a strict descent. The case `n=3` is already the
one-edge terminal.

Every strict drop costs at most two source-adjacency moves, including a
possible matching preparation. Induction on `e` reaches the one-edge class
in at most `2(e-1)` moves. Every nontransitive native target is incident to
at least one such source, so the claimed component statement follows.

## 4. Quantitative survivors and limits

The nonempty source component has radius at most `2(M-1)` about the
one-edge class and diameter at most `4(M-1)`. It is an all-size structural
result, not an inference from the finite census.

It also supplies an explicit, deliberately weak, native spectral bound.
Give the `Z=2^M-1` nonzero words equal weight and use the exact native
coupling `P`. Its stationary source masses are `r_i/Z`. If two source
states share a target with column mass `c_t`, their edge conductance is

\[
 \pi_iP_{ij}=\frac1Z\sum_t\frac{W_{it}W_{jt}}{c_t}
                    \ge\frac1{Z^2},
\]

since positive incidences are positive integers and `c_t<=Z`. Choose the
descent path from every source to the one-edge base, of length at most
`R=2(M-1)`. For any real source function `f`, Cauchy--Schwarz along these
paths gives

\[
 \operatorname{Var}_\pi(f)
 \le\sum_i\pi_i(f_i-f_{\rm base})^2
 \le R\sum_{\{i,j\}\text{ used}}(f_i-f_j)^2
 \le RZ^2\mathcal E_P(f).
\]

Thus, for `n>=4`, the native nontransitive component has gap at least

\[
                  \boxed{\frac1{2(M-1)(2^M-1)^2}}.           \tag{1}
\]

No uniform-in-`n` positive bound is claimed. The exact order-four gap is
`7/15`, so the auxiliary `8/9` bound is still false for this component.

There is an additional complement corollary. Pairing compressed source
graphs by complementation identifies the empty class with the complete
class. These lie in the two different components just proved. Consequently
the source-complement-paired native incidence is connected at every order
`n>=3`; converse-merging targets preserves connectedness. At `n=4` its
exact gap is `8/15`, again not `8/9`.

The endpoint of a source-to-source move is an actual Hamiltonian
presentation of the same native tournament *for that move*. Repositioning
the compressed source before the next move may change the native
tournament. Fixed native source, indexed ports, shared attachment data and
electrical mod-8 guards are therefore not claimed to survive the whole
descent. This does not finish the guarded native wiring problem.

## 5. Exact receipt

`check_native_descent.py` visits every nonzero labelled compressed graph for
native `n=3..6`: `1,7,63,1023` starting words. It implements the stated
descent, verifies every actual Hamiltonian junction, independently checks
both source-graph and target-tournament isomorphism under the full
permutation group, verifies lossless recompression and the radius bound,
and includes every exceptional case. No numerical spectral computation or
optimization-sensitive assertion is used.

```text
python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/coupling/check_native_descent.py
python -O -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/coupling/check_native_descent.py
```

The first write run passed 30,881 checks; receipt replays add one equality
check. The maximum observed descent lengths are `0,1,3,7`, respectively.
