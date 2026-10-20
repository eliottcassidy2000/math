# Two uniform repairs in the closed (n-3,2) endpoint layer

**PROVED:** the two marked families below, for n>=7.
**FINITE-EXACT:** all32,144 stated presentations through n10.
The remaining ports and arbitrary black-core connectivity are **OPEN**.

## Object and exact exceptional ports

Inherit the K convention and blue anti-involution test from
[cyclic endpoint detours](../kernels_games_detours_20261020f/detours.md).
Let r=n-2 and let V have a blue Hamiltonian word (v_1,...,v_r).
Adjoin a,z and fix j in {2,...,r}. Set
```
z->a,
a->v_i except v_j->a,
v_i->z except z->v_(r+1-j),
P=(a,v_1,...,v_r,z).                                     (1)
```
P is closed and blue, with scores d(a)=r-1=n-3, d(z)=2.
Conversely every closed blue word with these endpoint scores has exactly
this form. The additional incoming arc at a identifies v_j; anti-reversal
forces the reflected outgoing arc at z. No middle port is omitted or merged.

K_P(T) complements all nonpath edges and keeps path edges. Blueness is
unchanged across that edge. A zero endpoint-score excess is not an iff test.

## A low-degree cyclic-cut lemma

For a closed blue P on n>=7 vertices whose last score is2, its one-step
cyclic cut Q=(v_1,...,v_r,z,a) is black.

If Q were blue, the two positional anti-involutions would compose to the
cyclic shift by two, an automorphism. At odd n this acts transitively, forcing
the tournament regular with score2, hence n5. At even n it has two cycles of
length s=n/2. An automorphism of a tournament cannot have an even cycle:
a suitable power reverses an opposite pair. Thus s is odd. Its induced
tournament is regular, of score (s-1)/2. The score2 vertex also has an
outgoing edge across the two parity orbits, namely the closing path edge.
Therefore (s-1)/2+1<=2, forcing s<=3 and n<=6. Contradiction.

The boundary is real: the exact checker retains closed blue six-vertex
examples with a blue cyclic cut. This repairs the old shift-by-one shortcut
by retaining the actual shift-by-two permutation.

## Family A: open middle and one excluded port

Suppose v_1->v_r and j!=r-1. Put b=v_1 and W=(v_2,...,v_r).
The inherited guarded endpoint ledger is
```
Q1=(b,W,z,a),   Q2=(z,W,b,a),   Q3=(a,z,W,b).              (2)
```
Its direction guards are z->a, b->v_r and v_2->z. The last holds precisely
because j!=r-1. The old XOR proof applies without changing a single middle
arc and gives the literal target K_P(T). Q1 is black by the cyclic-cut
lemma; Q2 always has endpoint-score excess -1.

For Q3 inspect it at B=K_P(T). Anti-reversal would require the arcs
z->b and a->v_r to agree. They do not: in B,
```
[z->b]=[j!=r],       [a->v_r]=[j=r].
```
Both are complemented nonpath comparisons; this is exactly the exceptional
port information. Thus Q3 is black, even when its endpoint-score test ties.
All three moves are authenticated for every r>=5 in this family.

## Family B: closed middle without its near-end chord

Suppose v_r->v_1, v_(r-1)->v_1, and j!=r. Write b=v_1,
c=v_r, W=(v_2,...,v_(r-1)), x=v_2, y=v_(r-1). Use
```
Q1=(c,W,b,z,a),   Q2=(z,c,W,b,a),   Q3=(a,z,b,W,c).        (3)
```
The inherited no-chord XOR ledger is still directed: the only extra failed
join could be b->z, and j!=r is exactly its exclusion. Blue symmetry supplies
c->x, and the displayed premise supplies y->b.

Put d=d_V(c). Closure and c->x imply d>=2. The first two endpoint excesses
are d-1 and -1, hence black. The final excess is d-2, so only d2 needs more
than scores. In that case the outneighbors of c within V are exactly b,x.
At B=K_P(T), the inneighbors of c are exactly {b,x,y}, whereas the
outneighbors of a are {z,b,v_j}. Reversal of Q3 sends b to v_(r-2).
For r>=5 this vertex is none of b,x,y. Hence Q3 cannot be blue.
Again the endpoint is the complete labelled K_P(T).

## Controls and remaining branches

The generator takes every fixed-path blue middle of orders3..8 and every
j=2..r. This is32,144 presentations, with no isomorphism or strongness
filter. Every case satisfies the specified endpoint pattern. The two proved
families cover20,580 of these presentations at n7..10. All successive paths,
their actual positional blackness and the entire target matrix are checked.

The same script tries the old third (closed-chord) ledger and retains
failures, rather than interpreting a failed template as a disconnected edge.
Its legal instances all succeed at n8..10 in this finite universe, but
no universal conclusion is drawn. The open-middle port j=r-1, closed
no-chord port j=r, and the general closed-chord branch remain next targets.
No arbitrary replacement of an endpoint, converse identification, or clone
contraction is licensed by these results.

