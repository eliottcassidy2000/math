# Many distinct groups at one order inside one complete path-observer fiber

**Status: PROVED deduction; FINITE-EXACT permutation checks.** This is a
consequence of the canonical substitution group law already used in
[LIFTING.md](LIFTING.md) and [FIBER_ACTION.md](FIBER_ACTION.md), not a claim
of literature novelty. All isomorphisms are ordinary, arc-preserving ones.

## Construction and hidden-equality check

Fix h≥2, L=3^h and n=4L. The leaves are the two four-vertex tournaments
A0=singleton→C3 and A1=A0* on identically addressed labels. Every upper
quotient is the same oriented C3. As proved in the lifting note, independently
reversing leaf modules preserves every labelled subset Hamiltonian-path
count and every subset path-cover count. Thus **all witnesses below belong
to one full labelled observer fiber**, not merely a fiber of an aggregate.

Take the minimal pair from that note:

    D0=A0, E0=A1,
    D(j+1)=C3[Dj,Dj,Ej],    E(j+1)=C3[Ej,Ej,Dj].

Their upper stabilizers are all trivial and their groups are C3^(3^j).
They are distinct: their A1-leaf counts differ by one (Dj has one fewer),
and these canonical leaf types are isomorphism invariants.

For each r∈{1,…,h}, put a=3^(r−1) and define

    S(r)=C3[D(r−1),D(r−1),D(r−1)],
    S(j+1)=C3[S(j),Dj,Ej]       (r≤j<h).

The requested witness is T(h,r)=S(h). It has exactly **one upper node with
nontrivial quotient stabilizer**, at height r. This does not mean that there
is only one order-three automorphism. The initial active node has three
identical children. At every subsequent node S(j) has group order
3^(3^j+1), whereas Dj and Ej have order 3^(3^j); hence S(j) cannot equal
either child type. Dj and Ej are themselves distinct as just proved. Thus no
accidental upper rotation appears. This induction also avoids relying on
an unproved visual distinction between child words.

## Group theorem and proof

Writing elementary abelian groups additively, the canonical group law gives

    Aut(T(h,r)) ≅ ( (F3^a)^3 ⋊ <σ> ) × F3^(L−3a),
    σ(x0,x1,x2)=(x2,x0,x1),       |σ|=3.                 (1)

The three identical D blocks admit coherent identical-label transports, so
the extension splits with the displayed genuine block permutation. All
later nodes have trivial quotient stabilizer and merely add direct factors.
In particular every group in (1) has order 3^(L+1).

Let B=F3^(3a), and R=F3^(L−3a). For (b,t)∈B⋊C3,

    (b,t)(c,u)=(b+σ^t c,t+u).

A central element must have t=0, since σ and σ² act nontrivially on B.
It must then have b∈ker(σ−I), the diagonal copy of F3^a. Thus

    Z(Aut(T(h,r))) ≅ C3^(L−2a).

Since B, R and C3 are individually abelian, the derived subgroup is
im(σ−I)⊆B. Explicitly it is the space of triples whose sum is zero,
which has dimension 2a (also in characteristic three). Hence

    Aut(T(h,r))' ≅ C3^(2a).

These center orders are strictly decreasing as r increases, proving **h
pairwise nonisomorphic groups of the same order inside one full subset-H
fiber and one automorphism-order cell**. This is a lower bound of h, not a
classification of all abstract group types in that cell.

All groups have exponent nine. Indeed (b,1)^3=(Nb,0), where
N=I+σ+σ² sends (x0,x1,x2) to the diagonal triple with entry x0+x1+x2.
It is nonzero for some b, and its image has exponent three; elements with
t=0 have exponent at most three. More quantitatively,

    #{g:g³=1}=3^L+2·3^(L−a),
    #{g:ord(g)=9}=2·3^L−2·3^(L−a).                    (2)

This gives another abstract-group distinction without relying on centers.

They also all have **nilpotency class exactly three**. Put M=σ−I over F3.
Then M³=0, and M²=I+σ+σ²=N has rank a>0. Commuting the derived
subgroup with the full group applies M again (up to sign or invertible
factor), so the lower central series has

    γ2=im M,   |γ2|=3^(2a),
    γ3=im M²,  |γ3|=3^a,
    γ4=1.

Thus neither order, exponent nor nilpotency class distinguishes these groups;
center size or either nontrivial lower-central-series size does.

## Explicit permutations, not only abstract symbols

Use vertices (k,q), 0≤k<L, q∈{0,1,2,3}, encoded as 4k+q.
For each leaf k let c_k cycle (k,1)→(k,2)→(k,3)→(k,1), fixing all other
vertices. This is valid for both A0 and its converse. Since each embedding
places S first, the active subtree occupies the first 3a leaves. Define

    ρ(ia+k,q)=((i+1 mod 3)a+k,q)
      for i∈{0,1,2}, 0≤k<a;
    ρ fixes leaves k≥3a.

These L+1 actual permutations generate precisely (1). They preserve arcs:
leaf rotations do so internally and modules have uniform outside arcs; ρ
cycles three literally identical blocks and preserves their C3 quotient.
They satisfy c_k³=ρ³=1, [c_k,c_l]=1 and
ρ c_(ia+k) ρ^−1=c_(((i+1) mod 3)a+k). Their normal forms c^bρ^t are unique,
as the quotient block action detects t and disjoint leaf actions detect b.
The canonical decomposition law rules out further automorphisms.

At h=2 (n=36), the two groups both have order 3^10=59,049 and exponent 9:

| Active height r | Center order | Derived order | Elements of order 9 |
|---:|---:|---:|---:|
| 1 | 3^7=2,187 | 3²=9 | 26,244 |
| 2 | 3³=27 | 3^6=729 | 37,908 |

## Reproducible independent checks

From this directory run:

    python3 check_group_shapes.py > check_group_shapes.out
    python3 -O check_group_shapes.py > check_group_shapes.optimized.out
    cmp check_group_shapes.out check_group_shapes.optimized.out

The standard-library checker uses explicit guards, not Python assertions.
It independently canonicalizes cyclic child types and counts active nodes
for every proposed witness through h=7; compares literal cyclic-triple sets
for all witnesses through h=3; checks every generator against the actual
arc matrix; and enumerates generated permutations at h=1 and both h=2
cases. It computes centers by commuting with every generator, and derived
groups from all generator commutators, verifying their normality. It checks
all element orders and formula (2), and verifies that the centers and derived
groups have exponent three. It also generates γ3 from commutators of the
derived group with the generators, checking its order 3^a and centrality.
No 36! permutation search is used. Enumeration
establishes the generated group independently; exhaustiveness as the full
automorphism group uses the cited canonical decomposition law, not a hidden
claim of unrestricted automorphism search. The all-height result is proved
above, not extrapolated from the finite checks.

## Independent audit

A separate implementation closed the two literal degree-36 permutation groups, checked every generator against tournament arcs, computed centers and derived groups, verified their normality and all element orders, and checked the third lower-central term and its centrality. Both normal and optimized runs agree.
