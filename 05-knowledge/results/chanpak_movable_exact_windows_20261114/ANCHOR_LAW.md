# A common anchor from different tilted readers

Use the [bilateral profile](MOVABLE_TRIPLES.md)
(h,f,P,Q,S,T,K) of a marked core. Here f counts ordered
two-path covers, S counts full core paths ending in the
outcontact set U, and K is the core count with the apex.
The two transitive wings with addresses x,y have exact count

    H = fxy+Qx+Py+h(g(x)+g(y)−2)+K+T p(x)+S p(y).

All the witnesses here have x,y>0 so are strong. As in the
[finite-support defect](../chanpak_f1_physical_carry_20261114/ORDINAL_DEFECT.md),
K−P+S=h−E(U), where E(U) counts rejected complementary
pairs whose two adjacent endpoints are contacts. Each
unordered pair has an even contribution E_ij; write
E(U)=Σ_{i<j in U} E_ij.

## Compress the reader, keep the exact height

Take x=a u, y=b u−1, with a,b positive even and u≥1.
At a fixed h,f and tilt the only exposed linear coordinate
needed is W=bP+aQ. The entire ordinary count is

    H = f a b u² + (W−fa)u + h(g(au)+g(bu−1)−2) + h−E(U).

Hence a fiber at one exact W with E(U) moving in steps of
two produces a window at every u, not just at one residue.
For primitive tilt a₀:b₀, if two different (P,Q) agree in
W then (ΔP,ΔQ)=(a₀ v,−b₀ v) for nonzero integer v. Thus
both components of the primitive tilt are at most f−h for
any genuine compressed switch in this bank. More generally
at fixed parities and equal h,f, a specified exact difference
between two cores imposes one affine equation in x,y unless
both P,Q already match. Independent full-profile movement
and compensation along a given slope are different mechanisms.

## Synchronize different source scales and reader widths

If a,b are powers of two at least two, then

    g(au)=g(u)+(a−1)u+p(u),      g(bu−1)=g(u−1).

The complete common-height signature is now

    D=fab,   C=bP+a(Q−f+h)−h,   δ=h−E(U),
    H(u)=D u²+C u+h(g(u)+g(u−1)+p(u)−2)+δ.

This allows f, a and b to change together at fixed (h,D,C).
Here is a literal five-window with h=9, D=5376, C=627.
Rows are outgoing masks of the puncture; U is an apex mask.

| a | b | Rows | U | f | P | Q | δ |
|---:|---:|---|---:|---:|---:|---:|---:|
| 4 | 8 | (28,25,10,16,4,31) | 45 | 168 | 129 | 60 | −59 |
| 4 | 8 | (28,25,10,16,4,31) | 46 | 168 | 127 | 64 | −57 |
| 2 | 8 | (124,121,58,0,8,88,28) | 74 | 336 | 129 | 129 | −55 |
| 4 | 8 | (28,25,10,16,4,31) | 54 | 168 | 125 | 68 | −53 |
| 4 | 8 | (28,1,2,6,14,31) | 43 | 168 | 124 | 70 | −51 |

The middle row uses a different source cover count and different
left width; all five are read at the same absolute center.
For u=1 they give 5953,5955,5957,5959,5961. The earlier
[pure complementary chart](../chanpak_tight_mod_six_20261114/TIGHT_COMPLEMENT.md)
rejects all five, by its complete width inverse.
At each u the candidate counts cover all five odd phases mod 10.
If the prescribed target is in that window, choice of its phase
makes the discrepancy at most 8<10, hence zero.
The other [tilted certificates](DIAGONAL_WINDOWS.md) realize
the same heights using different puncture charts.

For any fixed dyadic signature, recover u by monotone search.
If F(u)=g(u)+g(u−1)+p(u)−2, then F(u+1)−F(u)≥−u by
1≤g(z)≤z+1. Also C≥−af+ah−h, D+C>0 and 2D>h.
Thus H strictly increases at each positive u. The positive
product term implies (D/2)u²≤H, giving a finite height bound
and bit-length binary search. Even freely varying dyadic a,b
over a fixed core bank gives only O(sqrt(X)) individual counts
up to X: the number of possible u is
O(sqrt(X)/(sqrt(fab))), and the dyadic sum converges.
This counts the tilted chart only, not the independent-wing charts.

## An exact small phase move and the target-selection problem

On the tilt x=X, y=R X−1 with X=M u and powers of two M,R≥2,
make a phase move to x=X+r, y=R(X−r)−1, for an integer
−M/2≤r<M/2. Its exact change in H is independent of u:

    −f r(Rr+1) + r(Q−RP) + h d(r) + T p(r),
    d(r)=2g(t)−t+p(t)−2,  t=max(r,−r−1).

To see the g term for r≥0, the two low M-block words of
X+r and X−r−1 are complements, with high prefixes u,u−1;
g(R(X−r)−1)=g(X−r−1). Their address cross terms cancel in
the change except for the reverse address of r. The
[complement identity](../chanpak_absolute_reflection_20261114/COMPLEMENTARY_READER.md)
then gives the displayed d. The negative phase interchanges
the two g arguments and gives the reflection r→−r−1.
The remaining terms are the literal change of product, contact
and parity terms. Tested phase moves do not imply a wider
contiguous window.

A clean remaining exact problem is visible already at x=2,y=1.
Put w_i equal to the start-minus-end response of vertex i on
complementary paths of the core. Then W=f−h+Σ_{i∈U}w_i
for tilt a=b, and the ordinary target equation becomes

    (N−4f−h)/2 = Σ_{i∈U} w_i − Σ_{i<j in U} E_ij/2.

Both off-diagonal and linear coefficients come from the same
physical core. The right side is a submodular quadratic
set function: adding a contact decreases the marginal effect
of every other contact by E_ij/2≥1. This gives precise target
tests for a chosen core and a way to look for a core with
a long interval of attainable levels. There is also a structural limit on a bank of unit interactions.
Make a graph joining contacts i,j only when E_ij=2. It cannot
contain K₄. Indeed the [two-contact minimum-frame classification](../chanpak_f1_physical_carry_20261114/CONTACT_PAIR_DEFECT.md)
shows that E_ij=2 survives deletion of other vertices: the
remaining order and its allowed two mixed types are inherited.
If four contacts had six such pairs, on the restricted four
vertices their combined E would be 12. But when all four are
contacts, E counts every ordered nonempty two-path cover,
at least 2⁴−2=14 by assigning each proper nonempty subset and
its complement a Hamiltonian path. This is a contradiction.

For a nontransitive core, every pair with E_ij=2 must also hit
every directed triangle, since its complementary remainder is
transitive. On r≥3 contacts there are therefore at most 3r−6
unit pairs (fix any one triangle). In a transitive core the
same bound holds: directly, an interval of separation d gives
E_ij=2^(d−1)+[d=1], so at most 2r−3 pairs can be unit.
All other coefficients E_ij/2 are at least two. This rules out
large independent unit-penalty pools, but does not rule out
balanced cancellation of larger coefficients.

Submodularity alone
does not prove interval coverage or guarantee the needed
choice of f,h at every target.

Run **python3 replay_relocation.py** for independent literal
window checks, nonoverlap with the old chart on the named
examples, and the exact phase and shared-height identities.
