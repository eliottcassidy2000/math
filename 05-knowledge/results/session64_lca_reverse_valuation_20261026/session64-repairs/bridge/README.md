# A degree-and-port certificate for cyclic-shift black detours

2026-10-26. **PROVED, all finite orders in the stated scope.** The three
blackness tests in a cyclic-shift replacement can be reduced to one
endpoint-degree mismatch and one permutation/port test. This is an
explicit sufficient certificate, not a guaranteed-existence theorem for
the required Hamiltonian word and not full black-core connectivity.

## 1. Literal source, target and candidate word

Let T be a tournament with a directed closed blue Hamiltonian word
`P=(0,1,...,z)`, where `z=n-1` and n>=3. Blue means that the actual
positional reversal is an anti-automorphism. Let

```
k=d_T^+(z),        I=N_T^-(0) minus {z},
B=K_P(T),          J=I union {1}.
```

Since P is blue, `d_T^+(0)=n-1-k`; since it is directed and closed,
`|I|=k-1` and 1 is not in I. The target has

```
N_B^+(0)={1,z} union I,       |J|=k.
```

Choose a directed Hamiltonian word W of `B[V minus {0}]`, starting at 1
and ending at z. Set

```
Q1=C=(1,2,...,z,0),       Q2=(W,0),       Q3=(0,W),
g=r_Q3 composed with r_P.
```

All words and permutations act on the original labels. No clone or
module is contracted.

## 2. Sufficient criterion

The three words are a directed black detour from T to its exact labelled
B whenever

1. `d_T^+(1) != k`; and
2. either g has even order, or `g(J) != J`.

Both guards are explicit. For a permutation, even order is equivalent
to at least one even cycle in its cycle decomposition.

The conclusion concerns all finite n for which these hypotheses hold.
It does not assume that every source has a W satisfying the guards.
The ordinary W=(1,...,z) is always available, but gives Q3=P and fails
the second guard, so Hamiltonian existence alone is insufficient.

## 3. Legality and the literal target

Let `S=K_C(T)`. Since

```
E(C) xor E(P)={01,0z},
```

S agrees with B on all pairs away from 01 and 0z. Thus W is directed in S,
and `z->0` is retained by C. This makes Q2 directed. In S, `1->0` holds;
Q2 omits 01, so its K move makes `0->1` true. Every W edge is retained,
so Q3 is directed at the next source U.

Moreover `E(Q2) xor E(Q3)={01,0z}`, giving

```
E(Q1) xor E(Q2) xor E(Q3)=E(P).
```

Three full complements and this edge XOR imply
`K_Q3(K_Q2(K_Q1(T)))=K_P(T)` pair by pair, on literal vertices.

## 4. The two endpoint tests collapse to the same mismatch

An anti-automorphism pairs the first and last vertices of its word;
their outdegrees must sum to n-1. For Q1 in T this sum is

```
d_T^+(1)+d_T^+(0)=d_T^+(1)+n-1-k.
```

It is unequal to n-1 by guard 1, so Q1 is black.

At S the degree of 0 is k-1: its original incoming ports I become its
outneighbors, while both special pairs01 and 0z point into0. Since 1 is
first in C, only its outgoing path arc is retained, giving

```
d_S^+(1)=n-d_T^+(1),      d_S^+(0)=k-1.
```

The endpoint sum for Q2 is therefore `n+k-1-d_T^+(1)`, again unequal
to n-1 by exactly the same guard. Thus Q2 is black, regardless of the
interior order of W.

## 5. The final word is controlled by a literal port set

Blue/black status of a word Q is unchanged across its own K edge:
its unordered path-edge set is invariant under its positional reversal,
so K complements both sides of every anti equation by the same bit.
Consequently Q3 is blue in U iff it is blue in B. P is also blue in B.

If Q3 were blue there, the composition `g=r_Q3 r_P` of these two anti
maps would be a tournament automorphism of B. Because Q3 and P have
first vertex0 and last vertexz, g fixes 0 and z. Therefore g must preserve
`N_B^+(0) minus {z}=J`. The second alternative in guard 2 excludes this.

No tournament automorphism has even order. Indeed, an automorphism of
even order has a nontrivial involution power; that permutation exchanges
some pair of vertices, which reverses their tournament arc, impossible
for an automorphism. The first alternative in guard 2 likewise excludes
Q3 being blue. Hence all three moves are black.

## 6. Exact hostile boundaries: the guards are not iff statements

The independently replayed 20 n12/n14 certificates all satisfy both
sufficient guards. Of the six cyclic n11 certificates, five satisfy the
degree guard. All 26 cyclic n11/n12/n14 certificates fail port-set
invariance; ten also have even permutation order.

The remaining cyclic n11 certificate is a positive equality hostile:

```
T=(126,156,440,1008,1248,1474,1926,1793,1555,1079,15),
k=4,        d_T^+(1)=4.
```

Its first cyclic word is nevertheless black. An exact census of all
1-to-10 Hamiltonian words in `B minus {0}` gives 132 candidates: 127 make
both remaining moves black and five make Q3 blue. Thus degree equality
is outside this sufficient certificate, not a connectivity obstruction.
The supplied successful word
`W=(1,5,6,7,2,3,4,8,9,10)` creates the exact intermediate clone block
(5,4), retained in the separate order11 certificate.

Conversely, for the genuine order-seven hostile

```
T=(6,4,8,19,103,79,15),
B=(122,124,121,112,32,64,16),
```

the degree mismatch holds, but every Hamiltonian word of B is blue.
Hence no choice of W can pass guard 2. Failure here is a real black-target
obstruction, established by an exact independent HP census, not inferred
from failure of this sufficient menu.

The n11 source whose stored certificate uses a noncyclic first word also
has 122 endpoint-fixed W candidates, 119 of which give cyclic three-black
repairs. The noncyclic discovery was therefore not a cyclic-bridge
obstruction. The bounded block-permutation search had simply found a
different short certificate first.

## 7. What this changes and what remains

This is an actual operation-level connection from endpoint degrees and
labelled ports to the blackness of a three-step K detour. It exposes the
remaining search coordinate precisely: find an endpoint-fixed Hamiltonian
word W whose reflection product breaks J or has even order. It does not
replace W with an isomorphism class, a degree sequence or a contracted
clone core.

The separate exact small-order repairs establish the sharp closed
endpoint-four threshold n>=8. This all-size sufficient lemma may be used
at other endpoint degrees, but guarantees no all-source repair there.
Open marked paths and unrestricted black-core connectivity remain open.
The independent bridge checker/review supplies finite controls; the
all-size conclusion rests on the explicit proof above.
