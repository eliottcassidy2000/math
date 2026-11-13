# Cyclic frames: the choice of origin and its parity

**Proved.** This refines THM-4729's closed/open split and THM-4741's
rotation count. We count tournaments together with a marked directed
Hamiltonian cycle, modulo isomorphisms preserving that directed cycle.
We do not identify different cycles of the same tournament unless an
automorphism identifies them. In particular, this is not an enumeration
of the unmarked merged nodes.

There is an intrinsic origin-shift action on **line instances** of the
merged metagraph, before passing to its node support. Each complement
line has a unique closed endpoint (THM-4729). Interpret that tiling as
a directed Hamiltonian cycle rooted at its first vertex. Rotate the
root, take the resulting closed tiling, and recover the line by
complementing its free tiles. This conjugates cyclic shifts of closed
frames to an action on line instances. All lines in an orbit share the
merged node at their closed end; the open end need not stay in one
node. The action is not claimed to preserve individual line colors or
the unmarked node graph. The closed endpoint, not an arbitrary tiling
of the pair, is essential to defining the action.

For a concrete order-four check, start with 0→1→2→3→0 and
chords 2→0, 3→1. The complementary open endpoint at root 0 is
transitive (all arcs point forward in its frame), whereas at root 1
the open endpoint is nontransitive (rows of outgoing bit masks
14,4,8,2, which contain the cycle 1→2→3→1). Thus the open merged
node does not stay constant under an actual change of origin.

## 1. Origins form a torsor; blue origins detect a reflection

Write the cycle vertices as v_i, indexed by Z/n with successor i+1. Let
R be the subgroup of shifts i↦i+r that preserve all arcs. Let A be the
set of offsets a for which the reflection i↦a−i reverses all arcs.
Two reflections differ by a shift, and composing such a reflection
with an arc-preserving shift gives another one. Thus A is empty or a
coset of R. The order of R is odd: an even-order cyclic rotation group
would contain a nonidentity involution, which interchanges some pair
of vertices and hence cannot preserve the arc between them.

Root i means the Hamilton path starting at v_i and following the cycle.
Its reflected path reverses the index j to 2i−1−j. Hence the rooted
closed tiling is blue exactly when 2i−1 belongs to A. The equation
2i−1=a has exactly one solution for odd n. For even n, any a in A is
odd, since an even offset reflection fixes two vertices and cannot
reverse the arc between them. An odd a gives exactly two solutions.
Therefore an unrooted cycle with a reflection has exactly b|R| blue
origins, where b=1 for odd n and b=2 for even n. Two origins give the
same rooted tiling exactly when their difference lies in R. There
are thus exactly b **distinct blue tilings** in a reflection-symmetric
cycle orbit, or zero if A is empty.

Changing origin preserves the existence of blue origins, but does not
in general preserve the color of the chosen tiling. At n=4, the unique
cyclic orbit has four distinct rooted tilings but only two are blue.
"Self-converse"
here means with the directed cycle marked: an isomorphism to the
converse carries the cycle to its reverse as a directed cycle in the
converse. This is precisely existence of a reflection in A; an
anti-automorphism of the underlying tournament need not preserve the
chosen cycle.

## 2. Exact census by rotation stabilizer

Fix the directed cycle arcs i→i+1. For n≥3 let D=n(n−3)/2 count
unordered non-cycle chords and let F=floor((n−3)/2). There are 2^D
labelled chord orientations. For an odd divisor d of n let R_d be the
unique rotation subgroup of order d. It acts freely on chords: a
nonidentity shift fixing an unordered chord would swap its endpoints
and have order two. Therefore exactly 2^(D/d) orientations are fixed
by R_d. There are no rotation stabilizers of even order.

Now impose the rooted blue condition, i.e. the anti-reflection
j↦−1−j. It acts on the D/d rotation-orbits of chords. On a fixed
orbit, indexing chords by the d rotation elements makes its action
j↦c−j, so there is exactly one geometrically fixed chord. Conversely
each geometrically fixed chord belongs to a fixed orbit. A fixed chord
has its two endpoints swapped by the reflection: at odd n the
reflection has exactly one fixed vertex, at even n it has none. There
are (n−1)/2 swapped pairs at odd n, of which one is a cycle edge, and
n/2 swapped pairs at even n, of which two are cycle edges. Thus there
are F fixed orbits of chords. Each fixed orbit allows both orientations
of its representative chord, because endpoint exchange and direction
reversal cancel; choosing that chord as representative propagates
consistently by the rotations. Each pair of exchanged orbits permits
one free bit, because the other is forced. Consequently the count fixed
by R_d and blue at root 0 is 2^((D/d+F)/2).

For odd r|n, define
```
U_r = sum_{r|d|n, d odd} mu(d/r) 2^(D/d),
B_r = sum_{r|d|n, d odd} mu(d/r) 2^((D/d+F)/2).
```
Divisor Möbius inversion counts frames with exact rotation stabilizer
order r, respectively all and those blue at root 0. There are rU_r/n
unrooted directed-cycle-marked classes with stabilizer order r. There
are B_r/b self-converse classes among them: by Section 1 each such
rotation orbit of chord masks contributes b distinct masks blue at
the fixed root 0.

## 3. The first split and its large-size rate

Direct enumeration of chord assignments, independent of inversion:

| n | Unrooted cyclic classes | Reflection-self-converse |
|--:|--:|--:|
|3|1|1|
|4|1|1|
|5|8|8|
|6|88|16|
|7|2344|256|
|8|131072|1024|

For every n≥6 the number of all cyclic classes is at least 2^D/n,
whereas the number of reflection-self-converse ones is at most
2^((D+F)/2). The strict inequality n<2^((D−F)/2) holds at n=6
and propagates: D(n+1)−D(n)=n−1, and F grows by at most one,
so the right-hand side grows by a factor at least 2^((n−2)/2),
larger than (n+1)/n. Thus the split at n=6 persists at all larger
sizes.

For any nonidentity rotation its fixed orientations number at most
2^(D/3) (zero for even order), so Burnside gives the total number
as 2^D/n times (1+O(n 2^(−2D/3))). The self-converse count is
exactly 2^((D+F)/2)/b by summing B_r. The fraction therefore is
asymptotic to (n/b)2^(−(D−F)/2). This gives a micro/macro transition
with a proved global threshold and decay scale. It does not give
an all-order cutoff for scalar Hamilton-count omissions.

## 4. What is shared with Goldbach, and what is not

The swap p↔N−p acts on the set of prime pairs for a given N. If
R_2(N) counts ordered representations of N as two primes, the
unordered count is (R_2(N)+1_{N/2 is an integer prime})/2.
The fixed point equation is again a doubling equation, now 2p=N
on the integer line. For odd N the swap has no fixed point and
the only possible unordered pair is {2,N−2}; for even N
the arithmetic allowance for a diagonal pair does not force
any off-diagonal witness. The common mechanism is the midpoint
equation under an involution, not any correspondence between
prime pairs and Hamilton paths. The cyclic frame admits changes
of origin as actual relabelings, but an integer translation of
both Goldbach summands changes N and does not preserve primality.

## 5. The induced parity law at merged nodes

For a self-converse tournament T, consider its automorphism orbits
of directed Hamiltonian cycles. Any anti-automorphism induces the same
involution on this set, after reversing the cycle order; any two such
maps differ by an automorphism. Write I(T) for the number of fixed
cycle orbits. Section 1 shows that the closed blue incidence at its
merged node is exactly b I(T). Non-self-converse nodes have no blue
incidences. An incidence counts a loop end separately; closed and open
ends are distinct framed tilings.

The total number of blue incidences at a self-converse node is odd by
THM-643. Thus for even n the open blue incidence is odd. For odd n,
I(T) has the parity of the total number of cycle orbits. Each such orbit
has odd size because the tournament automorphism group has odd order,
so the open blue incidence is congruent to 1+HC(T) modulo 2, where HC
counts directed cycles modulo rotation, not modulo relabelling.

There is a literal secondary graph at even n: pair each blue closed
frame with the frame cut at its antipodal origin. This stays blue,
has no fixed masks (the half-turn would be an order-two automorphism),
and preserves the closed merged node. Turn the pair of blue lines into
one edge between their two OPEN endpoints. The resulting multigraph
has odd degree at every self-converse node, zero at all others; loops
count twice. This is a canonical reorganization of existing blue
incidence, not an enumeration of all merged nodes.
