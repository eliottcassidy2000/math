# A minimal extra bit for a genuine port-legality test

Status: **PROVED** for the specified one-backward-bridge joins. The statement
concerns future H and the existence of a Hamilton path starting at the port;
it is not a complete graph or metagraph carrier.

## P1. An exact multiplicative coordinate

Let (T,a) be any nonempty pointed tournament. Let B be a disjoint nonempty
tournament with an entrance b and a designated new port q. Form T star B by
directing all cross edges T->B, except the single edge b->a. Write S(T,a)
for the number of Hamilton paths starting at a. Define kappa(B;b,q) to count
ordered spanning covers of B by:

- one path starting at q and ending at b;
- one further path, possibly empty.

Then

    S(T star B,q) = kappa(B;b,q) S(T,a).               (P1)

Proof: a Hamilton path starting in B must enter T, and the only B->T edge
is b->a. It can do so only once. Thus its restriction to T is a Hamilton
path starting at a; the initial and final B segments are exactly the two
paths counted by kappa. Conversely join any such cover through any such
T path. If the second B segment is nonempty, its first vertex is not b,
so the required T->B exit edge exists. This is a bijection, including an
empty final B segment and q=b. No assumption about T's source language enters.

For each of the65 native blocks, the retained reversed spine runs from q
to b. Hence kappa>=1. Exact independent cover enumeration computes all65
factors, whose distinct values are

    {1,2,3,4,5,6,7,8,9,10,11,12,13,18}.

All1270 literal pointed joins in the declared tests agree with(P1), including
every tournament through order4 at every port with four selected blocks,
and all65 blocks on a start-zero and a start-positive source.

## P2. The coarsest joint observer

Retain the inherited h,x,y and one bit theta=1 if S(T,a)>0, else0. For every
bank append, h,x,y update by the existing matrix and theta is unchanged.
Therefore this is a closed carrier for every future observation

    (Hamilton-path count, existence of a Hamilton path starting at the port).

It is the coarsest such observation equivalence on nonempty pointed sources.
Indeed the empty, a, aa contexts read

    h, h+2x, h+2x+2y,

so distinct triples are separated by those contexts. Distinct theta values
are separated already by the current legality observation. Conversely equality
of all four retained values gives equality after every future bank word by
the exact update laws. The carrier consists of THREE unbounded integer
coordinates and ONE bit; it is not a four-state finite automaton.

The extra bit is necessary. The literal five-vertex source and four-vertex q1
source in CONTEXT_ENTRY have equal(h,x,y)=(3,3,3), but start counts0 and1.
No further bank append can turn the former into a port-preserving bank source:
zero is absorbing, while every bank source has its retained Hamilton spine
starting at the port. Observer-image entry therefore does not repair literal
source legality in this example.

## P3. If the full start count is observed, four linear coordinates are sharp

For a nonempty source the vector(h,x,y,S) evolves by block-diagonal matrices
diag(P_l,kappa_l). The empty formal root needs its own initialization: compute
S on the first actual block, then use(P1). It must not be treated as a real
zero-start tournament. `port_guard.json` records both first-block counts and
subsequent factors.

The rational linear dimension is exactly4 for future(H,S). Three actual
bank sources U,UB,UQ give independent h,x,y directions, as in ARITHMETIC A5.
The two literal equal-triple sources with S=0 and1 supply a fourth independent
direction. Conversely H after empty,a,aa together with S directly reads all
four coordinates. This gives a rank4 lower bound and the displayed4-coordinate
upper bound.

For an arbitrary block with kappa=0 the bit update becomes
theta'=theta AND(kappa>0); zero remains absorbing. Neither the bit nor S
classifies all other legality predicates, isomorphism types, source orders,
or possible re-markings. The literal word and clock still accompany every
construction in this packet.
