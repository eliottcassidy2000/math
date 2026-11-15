# Neményi's inverse direction: geometry, equilibrium, compatibility

This is the continuum-mechanics Neményi, not the multiple-comparison
test. In [his survey](https://www.sciencedirect.com/science/article/pii/S0065215608703004),
one chooses part of the geometry or response before solving for the
remaining boundary data. [Sáenz's five-constant theorem](https://doi.org/10.1007/s00707-008-0118-3)
gives a 2–5 real-parameter result under its plane elasticity
hypotheses. That theorem is not a graph theorem and contributes
no bound of five to the tournament problem.

A concrete Cartesian control explains the useful distinction.
With no body force and straight Cartesian principal trajectories,
the shear stress is zero. Equilibrium then gives
stress_xx=f(y), stress_yy=g(x), with two apparently arbitrary
functions. For homogeneous isotropic plane elasticity an Airy
stress potential is biharmonic. Here zero shear permits a
separated potential, so compatibility is f''(y)+g''(x)=0.
Thus f and g are quadratic with opposite quadratic coefficients:
there are five parameters before imposing the boundary tractions.
One may restrict to a patch where the two principal values
are unequal to avoid a degenerate trajectory direction.
This example distinguishes the *geometry* of trajectories
from the *physical admissibility* of the stresses.

## An exact equilibrium object for the tournament

For a fixed tournament W write H for its Hamilton path count,
s_i and t_i for the full path start/end counts at i, and U_ij
for the number of full paths using directed arc i→j. Each
vertex contributes outgoing usage except when last, and incoming
usage except when first:
```
sum_j U_ij = H−t_i;     sum_j U_ji = H−s_i;
div U = s−t.
```
The difference of responses with the same H and endpoint data
has zero row and column sums on the outgoing/incoming bipartite
lift. It is an invisible circulation at this observation level.
But equilibrium alone admits nonphysical solutions. For instance,
on a tournament with an arc a→b and a directed triangle on
the other three vertices, set the nominal H=1, start a and
end b, and place unit usages only on that arc and that triangle.
Every displayed vertex equation holds, yet no spanning path
has those arc usages. The cross-cut usage between the two parts
is zero. More generally, if A is a nontrivial vertex subset
and s(A)+t(A)=u, actual path traffic has cross-cut sum
at least max(u,2H−u). Each path crosses at least twice
unless its endpoints lie on opposite sides, and at most
min(u,2H−u) paths can have opposite-side endpoints.

There is an exact lifted compatibility test for **a multiset**
of paths. In the subset-state DAG take states (S,i), i in S,
and edges (S,i)→(S∪{j},j) when i→j and j is not in S.
Add edges from source to singleton states and from full states
to sink. Nonnegative integer flows decompose into source-to-sink
paths and project to arc/end marginals. This still represents
an arbitrary multiset, not the physical unit-weight sum over
all Hamilton paths. In the physical response each lifted edge
has weight equal to the product of the number of legal prefixes
at its tail and the number of legal suffixes at its head.
Both multiplicities are fixed by the actual tournament.

## The requested load is a joint split count

Add v beating only u of W. In a full path of the new
tournament, v first contributes s_u possibilities; v last
contributes H−t_u. With v inside, split at v: the suffix
must start at u, and the prefix is a nonempty path avoiding u
and covering its complementary vertices. It automatically
ends in an in-neighbor of v. Write C_u for the resulting
number of ordered nonempty two-path covers with second
path starting at u. Then the exact equation is
```
H(new) = H+s_u−t_u+C_u = H+2E_u.
```
The latter equality is the
[odd-cycle emptying identity](../chanpak_ternary_emptying_20261115/EMPTYING_EXCLUSIONS.md).
Thus the target demands
`C_u=N−H−s_u+t_u`, a **joint** observable that edge
conservation and ordinary undirected Hamilton counts do
not impose. Root admissibility is also needed: u must lie
in the first SCC of W if the augmented tournament is
to be strong.

For any chosen order there is a finite exact inverse with
the desired geometry imposed first. Choose binary variables
for the free frame slots, fix the outlet arcs, and use
the exact subset/last dynamic program of the augmented
tournament. Each term multiplies an arc bit a by a path count x
bounded by (|S|−1)!. In an integer formulation it can be
replaced by z with 0<=z<=x, z<=B a and z>=x−B(1−a).
Enforce the terminal sum N. This admits exactly the physical
framed witnesses, unlike the balance relaxation, but it is
still exponential in order. The
[complement-star family](COMPLEMENT_STAR_INVERSE.md)
is a constructive semi-inverse specialization: the geometry
fixes a graph class, the response has a closed binary law,
and the remaining compatibility can be decoded at a
prescribed integer without a general tournament search.

## A graph-side boundary equation with the opposite sign

Fix an undirected graph R and add x with binary contacts b_i
to its vertices. Count ordered Hamilton paths. Its exact count
after the insertion is
```
sum_i A_i b_i + sum_(i<j) B_ij b_i b_j,
```
where A_i counts full paths of R with endpoint i (one end
position per occurrence; for a singleton count both).
B_ij is the number of ordered nonempty two-path covers of R
whose first path ends and second starts at i,j, or j,i.
Removing the new vertex from an interior position gives exactly
such a split. If i and j are adjacent in R, concatenating at
that edge also gives the interpretation of B_ij as the number
of full paths of R using edge ij; if they are not adjacent
this full-path surrogate misses the entire interior contribution.

For the directed tournament's apex readout, the exact
[recontact energy](../chanpak_tail_return_20261115/RECONTACT_ENERGY.md)
instead has interactions −W_ij b_i b_j for nonnegative W_ij:
marking both endpoints makes a prospective odd-cycle
crossing unavailable. This sign and the ordered split data
are where inversion of the undirected path statistic stops.
