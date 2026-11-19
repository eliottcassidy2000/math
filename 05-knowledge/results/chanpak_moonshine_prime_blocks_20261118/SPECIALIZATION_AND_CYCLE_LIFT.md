# Specialization and physical cycle lifting at a marked seam

**Proved.** The full Hamilton count at every transitive block width
can agree while the underlying odd-cycle packing polynomials differ.
For a tight marked cut we can also identify the two extreme
exponential modes of the packing polynomial exactly. These give
a practical distinction between a scalar certificate and a stronger
cycle-level certificate for the remaining $38\cdot2^j+1$ family.

## What transfers from BSD and Hodge

The repository's [BSD/Hodge note](../../../07-reflections/bsd-hodge-the-polygonal-ladder-and-what-merges-kps7.md)
distinguishes two tasks: diagnose an object at a distinguished value of
a generating function, then ask what genuine geometric objects carry
that invariant. The [Selmer-growth session](../bsd_selmer_growth_20261013a.md)
uses a rank plus a growth profile, insisting that finite receipt data
not be confused with a full arithmetic realization. The [cycle-lifting
session](../../../07-reflections/hodge-conjecture-cycle-lifting-lrc14-codex-s232.md)
similarly warns against treating a compatible formal class as an
actual cycle certificate.

The analogy is precise about its limit. The [BSD rank statement](https://www.claymath.org/wp-content/uploads/2022/05/birchswin.pdf)
uses an order of vanishing at $s=1$, and the [Hodge conjecture](https://www.claymath.org/wp-content/uploads/2022/06/hodge.pdf)
asks for rational combinations of algebraic cycle classes. Neither is
being reduced to tournaments. The theorem below instead computes,
for literal finite tournaments, what a distinguished specialization
forgets, and what still has to be lifted.

## The exact two-visit cycle profile

For any tournament $T$, let $Z_T(z)$ enumerate collections of
vertex-disjoint directed odd cycles, weighted by $z$ per cycle;
unused vertices are free. The specialization $Z_T(2)=H(T)$ is the
[Grinberg–Stanley identity](https://arxiv.org/abs/2307.05569).

Take a quotient with marked vertex $x$, outneighbors $O$ and
inneighbors $I$, whose return-edge matching number is at most one
([marked cut criterion](MARKED_VISIT_MATCHING.md)). Replace $x$ by
a transitive block of width $m\ge1$. Define three polynomials on the
unexpanded quotient:

- $U(z)=Z_{T-x}(z)$.
- $V(z)=\sum_P z Z_{T-x-V(P)}(z)$ over directed outside paths $P$
  from $O$ to $I$ containing an odd number of vertices.
- $W(z)$ is the same sum for even order paths.

A packing uses the block in at most one run in total: each block
run in a cycle requires a distinct outside return path from $O$
to $I$, and disjoint returns would contain disjoint return edges.
If the unique outside path has odd order, it needs an even nonzero
subset of the block to make an odd cycle. If it has even order, it
needs an odd subset. Every chosen subset has a single increasing order.
The rest of the packing lies on the outside complement of that path.
Consequently
$$
Z_{T(m)}(z)=U(z)-V(z)+2^{m-1}(V(z)+W(z)).
$$
Thus the remaining ray at seed width $m=1$ amounts exactly to
$$
U(2)-V(2)=1,\qquad V(2)+W(2)=38\cdot2^j.
$$
These are equations on **one physical marked quotient**. In
particular $U(z)-V(z)=1$ is sufficient for the intercept equation
but is too strong as a universal search constraint. The condition
actually required is
$$
U(z)-V(z)-1\in z(z-2)\,\mathbb Z[z].
$$
The factor $z$ is the empty packing; the factor $z-2$ is the
Hamilton specialization. A nonzero multiple is a genuine cancellation,
not permission to replace the tournament by a signed sum of graphs.

### Two literal equal-count rays which separate at the cycle level

Both marks below have in- or outdegree one, so the theorem applies.
The masks are outgoing, with vertex 0 at the low bit.

| Masks | Mark | $U-V$ | $V+W$ |
|---|---:|---|---|
| 4,5,48,7,43,75,159,63 | 0 | $1$ | $36z+4z^2$ |
| 0,21,73,67,45,15,179,63 | 7 | $1+2z-z^2$ | $34z+5z^2$ |

Both have $1+44\cdot2^m$ Hamiltonian paths at every width. At the
polynomial level their difference is
$(1-2^{m-1})z(z-2)$ in table order. They agree even as entire
polynomials at width one, but differ at every width from two onward.
Exact cycle enumeration computes the table, and independent path
subset DP checks the specialization through four widths.

## The two extreme poles of a tight cut

Now suppose the return-edge matching number is $\nu\ge2$, and
there are exactly $2\nu$ outside vertices. Let $M$ count the
maximum matchings. For positive widths the cycle-packing polynomial
is a finite combination of $\lambda^m$, with integer bases
$1-\nu\le\lambda\le1+\nu$. In the generating function in block
width, nonzero modes correspond to poles at $1/\lambda$.
The coefficients of the two extreme modes are exactly
$$
[ (1+\nu)^m ]\,Z_{T(m)}(z)=M(z/2)^{\overline\nu},\qquad
[ (1-\nu)^m ]\,Z_{T(m)}(z)=M(-z/2)^{\overline\nu},
$$
where $a^{\overline\nu}=a(a+1)\cdots(a+\nu-1)$ and the left sides
mean coefficients in the unique expansion in distinct nonzero
exponential bases, not polynomial coefficients in $m$.

For the proof, fix a compressed packing with $r$ distinct block
runs. Label the runs by the immediately following outside vertex.
Assign each block vertex to one of the $r$ runs or leave it unused.
Use inclusion-exclusion to make every run nonempty and a parity
character per marked cycle to enforce odd order. The base of each
term is $1+\sum_{i\in S}\epsilon_{c(i)}$, where $S$ is the set of
enabled runs, $c(i)$ is the cycle containing run $i$, and
$\epsilon_c=\pm1$. We have $r\le\nu$ by the matching bound. The
extreme bases require $r=\nu$, all runs enabled, and all
characters positive or all negative. With only $2\nu$ outside
vertices the return paths then have two vertices each and exhaust
the exterior, hence they are the edges of a maximum matching.
A permutation of those labeled edges says which return follows
which around a cycle. For a permutation with $c$ cycles the parity
filter contributes $2^{-c}$ and the cycle weight is $z^c$.
Each outside cycle length is even before inserting clones, so
the all-negative choice adds a factor $(-1)^c$. Sum over the
permutations using $\sum_{\pi\in S_\nu} t^{c(\pi)}=t^{\overline\nu}$
(by inserting the last label). This proves the formulas.

In particular at $z=2$ the negative extreme vanishes, while the
positive one is $M\nu!$. The negative coefficient need not vanish
at other values: at $\nu=2$ it is $M z(z-2)/4$. The replay shows
this on rows (6,28,24,17,1) marked at 0, where $M=2$. This is a
literal special-value cancellation. More generally every negative
mode must vanish at $z=2$ because the Hamilton run formula only
has positive exponential bases, but the tight-cut result exhibits
the parity/grouping mechanism and its precise factor.

**Next use.** Search the defect-14 star hosts with the scalar
equations on $U,V,W$ and track the excess multiple of $z(z-2)$
explicitly, rather than discarding either all cancellations or all
cycle data. The two equal rays prove why either shortcut loses
information. This does not yet place the $38$ ray or close Chan–Pak.

Run python3 cycle_lift_replay.py.
