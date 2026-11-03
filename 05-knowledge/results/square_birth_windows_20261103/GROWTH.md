# Edge birth, cycle freedom and the prime comparison

Status: **PROVED elementary deductions / FINITE-EXACT controls / CITED Ogg
classification.** Birth and Pell laws and the small path banks are inherited,
not new claims of priority. See [RECOVERY.md](RECOVERY.md).

## Count, root and parent are different coordinates

The older square-sum packet already proves that adding integer n creates

    D(n)={s^2-n : n<s^2<2n},
    b(n)=|D(n)|=floor(sqrt(2n-1))-floor(sqrt(n)).

The strict inequalities exclude a zero parent and a loop. Every old edge
persists. Between changes of the root interval, every named parent moves
down by one each step. A root enters when 2n-1 first reaches its square,
and exits when n reaches its square. The inherited Pell analysis classifies
simultaneous entry/exit events, including 9 and 25. Equal b(n) does not
detect these events.

| n | Roots of new edge sums | Actual old neighbors | b(n) |
|---:|---|---|---:|
| 8=9-1 | 3 | 1 | 1 |
| 9 | 4 | 7 | 1 |
| 10=9+1 | 4 | 6 | 1 |
| 24=26-2 | 5,6 | 1,12 | 2 |
| 25 | 6,7 | 11,24 | 2 |
| 26 | 6,7 | 10,23 | 2 |
| 27 | 6,7 | 9,22 | 2 |
| 28=26+2 | 6,7 | 8,21 | 2 |

The arithmetic reflection 24 <-> 28 about 26 does not preserve this port
stencil: 25 is a root-replacement event, whereas 26 is not. On the fixed
root chart 25 through 28, translating n by one translates each parent by
minus one. This is an actual local arithmetic map; it is not a graph
isomorphism between Q_n and Q_(n+1), which have different vertex sets.

## Why the first nontrivial path appears at 15

The recovered Session147 proof gives Q12 as three paths, with vertex sets
{1,3,6,8,10}, {4,5,11,12}, {2,7,9}. Vertex 13 joins the first two components;
14 joins the third, producing a tree with leaves 8,9,10. Adding 15 joins
1 to 10 through 15 and creates the cycle 1--3--6--10--15--1. The forced
endpoints 8 and 9 determine the path

    8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.

Its reverse is the only other directed choice. It leaves exactly edge
1--3 unused. The path extends at 16 and 17; adding 18 creates a third leaf
and destroys Hamiltonicity. Thus edge additions need not preserve a spanning
path when they also introduce a new required vertex. The singleton Q1 is
the trivial earlier exception to 'first at 15'.

The user's displayed sets are the components of Q18[{1,...,12,16,17,18}],
not the first 15 integers. They are paths on 6 and 4 vertices and a tree
on 5 vertices whose vertex 7 has degree 3. The omitted vertices 13,14,15
are essential connectors and cycle data. This exact identification is
inherited from Session147 and checked again here.

## Unused edges and the natural oriented abstraction

**Proposition (any finite graph).** If G has a spanning path P and m edges,
P leaves m-n+1 edges unused. Since G is connected, this is its cycle-space
dimension. Each unused edge, together with its unique segment in P, gives
a fundamental cycle; these cycles form a binary basis.

For a specified direction of P, the orientations extending P form a cube
of dimension m-n+1. Its coordinates are the unused edges. Different path
cubes may overlap, so summing their sizes counts incidences, not distinct
orientations. A basis for one path does not identify bases for all paths.
This is standard spanning-tree cycle algebra applied to the current guard.

**Proposition (any finite graph).** Directed Hamilton paths are in bijection
with acyclic orientations that admit a Hamilton path: orient every edge
forward in the path order. The order is unique, because reachability along
the path forces the relative order of every pair. Conversely an acyclic
orientation with a directed Hamilton path has every edge forward in that
order. Equivalently, an acyclic orientation has a unique topological order
iff it has such a path: without an arc between consecutive vertices in
a unique order, those two vertices could be interchanged.

**Complement boundary.** For the forward orientation of P, K_P reverses
every extra edge. If its endpoints occur at positions i<j, its reverse
closes a directed cycle along the P segment from i to j. Therefore K_P
preserves acyclicity iff G is exactly the path (including the singleton).
The five-cycle in Q15 is the first nontrivial prefix hostile.

Orienting instead by numerical birth order is a DAG too, but it has no
directed Hamilton path for n>=2: such a path would have to be 1,2,...,n,
and edge 1--2 is absent. Thus the useful order is a path order, not the
age order. Constructing it already requires finding a Hamilton path.

At Q25/Q26/Q27 the acyclic active counts are exactly 20/24/70, against
4,864/11,232/62,952 arbitrary active guard orientations. These equalities
are checked by topological elimination on every active orientation.
The smaller acyclic family cannot serve as a complement-closed source.

## When most edges really are unused

At 15, a path leaves 1 of 15 edges unused; at 25 it leaves 8 of 32 unused.
Only edge 1--15 of Q25 lies in no path at all. These are different quantifiers.
At 27 every edge belongs to some Hamilton path, although each path omits
10 of the 36 edges. At 28 every path omits 11 edges, and 7--9 is the only
edge omitted by all paths. The exact lists are in `growth.json`.

**Proposition.** The first n>=2 for which E(Q_n)>2(n-1) is 56. The inequality
holds for every n>=56. There is a Hamilton path at 56, supplied in `path56.json`.

Proof. Sum the exact b(j) for j<=56: E55=108=2(55-1), E56=111>110;
the complete smaller ledger has no positive excess. For 19<=j<=24,
b(j)=2. For j>=25 put k=floor(sqrt(j))>=5. Then

    2j-1 >= 2k^2-1 >= (k+2)^2,

since k^2-4k-5=(k-5)(k+1)>=0. Hence b(j)>=2 and the excess
E(Q_n)-2(n-1) is nondecreasing from 19 onward. This proves the whole
infinite tail, not just the finite scan. At 56 the supplied path uses
55 edges and leaves 56 unused.

Summing b(j)=(sqrt(2)-1)sqrt(j)+O(1) also gives

    E(Q_n) = (2/3)(sqrt(2)-1)n^(3/2) + O(n).

Incoming [Session148 growth](../session148_square_sum_20261103/growth/GROWTH.md)
sharpens the remainder to **O(sqrt(n))** using its exact cumulative formula.
With K=floor(sqrt(n)) and M=floor(sqrt(2n-1)), that formula is

    E(Q_n)=(M-K)n-K+K(K+1)(2K+1)/6
           -M(M+1)(2M+1)/12+ceil(M/2)/2.

It follows by summing the square-entry thresholds; the two order-n terms
cancel on expansion. This stronger incoming result is credited, and the
formula is checked against direct edge enumeration here through n=300.

Whenever a Hamilton path exists, the used fraction is asymptotic to
3(sqrt(2)+1)/(2sqrt(n)), tending to zero. Edge abundance does not by itself
prove Hamiltonicity; incidence and separator obstructions remain.

## The first 25 primes, and the 15 supersingular primes

Three objects must be distinguished:

1. Q_k on integer values 1 through k, possibly indexed by the label p_k.
2. Q_(p_k), an integer prefix of prime order.
3. The square-sum graph on actual values p_1,...,p_k.

The [ledger](LEDGER.md) records all three through k=25. At k=25 their
edge counts are 32,258,24. In the third graph the component sizes are
22,1,1,1, with isolated primes 37,43,67. Its 24=n-1 edges do not form a
spanning tree: the graph has cycle rank 24-25+4=3 and is disconnected.
The prime-rank bijection is not adjacency-preserving: prime edge 2--7
maps to ranks 1--4, while integer edge 1--3 maps to prime values 2--5.
Neither resulting sum is a square.

For any graph on prime values with square-sum edges, a bipartition is

    A={2} union {odd primes congruent to 1 mod 4},
    B={primes congruent to 3 mod 4}.

Two odd primes must sum to an even square, hence to 0 mod 4; an edge at
2 has the other prime 3 mod 4. A spanning path therefore needs part sizes
differing by at most one. This is inherited prime-square-sum structure.

The Ogg set is 2,3,5,7,11,13,17,19,23,29,31,41,47,59,71: precisely the
primes for which all supersingular j-invariants are defined over F_p.
This published classification concerns X_0(p)^+; it is not the shorter
unquotiented X_0(p) genus-zero list. See Duncan--Ono,
[The Jack Daniels Problem, introduction and (1.1)](https://arxiv.org/html/1411.5354v2).
Its square-sum bipartition is

    A={2,5,13,17,29,41},
    B={3,7,11,19,23,31,47,59,71}.

Sizes 6 and 9 forbid a Hamilton path. Thus the equality '15 objects' cannot
identify the first nontrivial Q15 path with a path through all Ogg primes.
A possible deeper map would have to change the predicate or retain extra
data; no such map is claimed here.

PRIMEGAME gives another precise distinction. Permuting prime registers
induces a bijection on positive integers by permuting valuations. It
preserves divisibility and conjugates the ordered FRACTRAN transition when
all fractions and inputs are renamed consistently. It need not preserve
addition: swapping primes 2 and 3 sends square-sum edge 1--3 to nonedge
1--2. Reading numerical prime-output exponents also requires retaining the
decoding map. The register-machine graph and the square-sum graph share
possible labels, but their intrinsic operations are different.
