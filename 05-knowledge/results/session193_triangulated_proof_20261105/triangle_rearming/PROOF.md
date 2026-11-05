# A sharp interface budget for breaking tournament profile content

Session 193, 2026-11-05. **PROVED elementary orbit theorem and constructions;
FINITE-EXACT certificates; independently reviewed.** The unrestricted
Hamilton-path spectrum conjecture with exceptions `{7,21}` remains **OPEN**.
This packet proves when genuinely nonuniform attachments can escape the
intact-module divisor, and how much interface they must supply. It does not
prove a literal ternary return operation.

## 1. Inheritance and the new question

The closest mechanism is Session 193 [Session 193 rooted-clock profiles](../profile_substitution/PROOF.md),
Proposition 2: keeping a tournament S as an intact substitution factor
preserves its complete path-cover content g(S). Its exact rooted-clock
operators and transported two-level fan therefore cannot send a clock with
g(S)>1 to the literal digits `3H(S)+2` or `3H(S)+4`. The Session 193 induction
packet supplies an actual H41-to-H127 old-arc reset, but no indefinitely
reusable return.

The hostile example here is the strong nine-vertex clock
`C3[C3,C3,C3]`: H=3159 and profile content 81. The corrected near miss is
that destroying module homogeneity is not by itself enough to remove its
arithmetic obstruction. Even an exterior vertex with arbitrary nonuniform
incidences to all three triangle modules leaves every ordered path-cover
count divisible by three.

The least-used sidecar is the set of old modules touched by *used path
adjacencies* to the exterior. There are at most twice as many such modules
as new vertices. This converts the full-profile obstruction into a sharp
quantitative interface budget.

Concept board: intact-module content; actual path adjacencies; independent
triangle rotations; split-interface births; old-arc resets; affine rank four.
The anchor is a genuinely physical escape. The niche is a sharp lower bound
on its interface. The wildcard is that the sharp one-vertex splitter has
H41, reconnecting the budget to the finite reset frontier. No identification
with either H41 presentation in the induction packet is assumed.

## 2. Complete profiles and the interface-budget theorem

For a finite nonempty tournament T, let pc_T(j) count spanning covers by j
unordered nonempty directed paths, and set

    F_T(j)=j! pc_T(j),
    P_T(d)=sum_j binom(d,j) F_T(j),
    g(T)=gcd_j F_T(j)=gcd_(d>=1) P_T(d).

Thus F counts ordered nonempty path tuples, and P(d) counts paths in d
labelled bins, with empty bins allowed. In particular P(1)=H(T). Equality of
the gcds follows in both directions from the displayed binomial expansion
and its integer finite-difference inverse. Values d=1,...,|T| suffice.
Ordered or labelled covers matter: the group action below need not be free
on *unordered* path collections.

A subset M is a module if every vertex outside M has the same orientation
to all vertices of M. A triangle module induces C3 and has its independent
cyclic rotation of order three.

**Theorem 1 (arbitrary split-interface budget).** Let T contain m pairwise
disjoint triangle modules. Let U be any tournament obtained by adjoining k
new vertices, preserving every old arc. Every old/new and new/new arc may
be chosen arbitrarily. Then

    3^max(m-2k,0) divides F_U(j) for every j,
    and therefore 3^max(m-2k,0) divides g(U).             (1)

T need not be strong, and its triangle modules need not cover all vertices.
Neither the external incidence pattern nor the order of the births is
restricted. The statement concerns every path-cover coordinate, not merely H.

**Proof.** Consider an ordered spanning cover of U by j paths. Call an old
triangle module *exposed* if some used consecutive path adjacency joins a
vertex in it to a new vertex. Each new vertex has at most one predecessor
and one successor in its path. Consequently at most 2k old modules are
exposed.

Partition the covers according to their exact exposed-module set E. Within
one such class, independently rotate every unexposed old triangle. These
rotations preserve every old/old adjacency because they are independent
module automorphisms of T. They preserve used old/new adjacencies because no
unexposed module participates in any such adjacency. New/new adjacencies
are unchanged. The class E is preserved as well.

The resulting group `(C3)^(m-|E|)` acts freely on the ordered covers. Indeed,
a nonidentity element moves some old vertex, and that named vertex occurs
in one fixed position of an ordered path word. Thus it cannot fix the
ordered tuple. Every orbit in this class has size `3^(m-|E|)`, which is
divisible by `3^max(m-2k,0)`. Summing over E proves the first statement.
The second follows by taking a gcd. QED.

The same argument works for any prime p when each of m disjoint modules has
a nontrivial automorphism of order p: replace every 3 in (1) by p. A group
element need not move every vertex of its module. It is enough that a
nonidentity element moves at least one named vertex, since all vertices
occur in a spanning ordered cover. The sharp constructions below are for
p=3 only.

**Corollary 2 (old-arc reset budget).** Suppose U additionally changes old
arcs, and A is the set of old triangle modules containing an endpoint of
any changed old arc. Then

    3^max(m-|A|-2k,0) divides g(U).                    (2)

In particular, if f old arcs differ from T, then

    3^max(m-2f-2k,0) divides g(U).                     (3)

**Proof.** The m-|A| unaffected triangles still have their independent
old-graph rotations. Apply the same argument to them. Each changed old arc
has endpoints in at most two of the original modules, so |A|<=2f. QED.

Here f counts arcs between the original vertices; all incidences involving
new vertices are already unrestricted. This is a statement about the final
labelled graph and allows any intervening order of births and reversals.
Deleting old vertices is outside the statement.

## 3. The exponent is sharp, at every size

Use outgoing adjacency bitmasks on the named vertices 0,...,n-1. Start with

    C = C3 = [2,4,1],
    C -> C = [58,60,57,16,32,8].

Adjoin vertex6 beating exactly old vertices0 and3. Its outgoing mask is 9.
The resulting seven-vertex tournament is

    B = [58,124,121,16,96,72,9].

It is strong and its exact profile values are

    H(B)=41,             P_B(2)=712,            g(B)=1. (4)

Thus one vertex breaks two separate triangle-module divisors. It is not a
uniform substitution wrapper around the old six-vertex tournament.

For a single leftover triangle, adjoin vertex3 with sole old out-neighbor0:

    D = [2,12,9,1],
    H(D)=5,              P_D(2)=32,             g(D)=1. (5)

Both (4) and (5) are verified by independent literal permutation-and-cuts
profiles and vertex-subset path counts. They are finite input identities,
not extrapolations.

For an ordered join, binwise concatenation gives

    P_(X->Y)(d)=P_X(d)P_Y(d).                         (6)

**Theorem 3 (sharp partial repair).** For every m>=1 and every
`0<=k<=floor(m/2)`, the ordered join of m triangle modules has an induced
extension by exactly k vertices whose profile content is

    g=3^(m-2k).                                      (7)

Its H and second bin evaluation are

    H=41^k 3^(m-2k),
    P(2)=712^k 12^(m-2k).                            (8)

**Construction and proof.** Pair the first 2k consecutive triangles and
perform the physical B birth inside each pair. Give each new vertex the
ordered-join incidences to every other block. This preserves every original
arc. The final graph is

    B -> ... -> B -> C -> ... -> C,
       k copies         m-2k copies.

Equation (6) gives (8). The gcd of those two displayed values is exactly
`3^(m-2k)`, since 712=8*89 and the only common prime is the indicated power
of three. Hence g divides this number. Theorem 1 gives the reverse
divisibility, proving (7). QED.

**Theorem 4 (exact minimum primitive completion).** The minimum number of
added vertices needed to turn the ordered join of m triangle modules into
an induced extension with g=1 is exactly `ceil(m/2)`.

**Proof.** Theorem 1 requires 2k>=m. Write m=2a+r, r in {0,1}. The graph
`B^a -> D^r`, with all products interpreted as ordered joins, uses a+r
births and has

    H=41^a 5^r,          P(2)=712^a 32^r.

These two integers are coprime, so g=1. QED.

The statements minimize births for *profile content one*, not for a
specified target H or minimum tournament order. Most joins here are not
strong. That is permitted and explicitly separate from the strong example
in the next section.

## 4. An exact strong two-birth escape from content 81

Let S=C3[C3,C3,C3], with the three consecutive labelled triangle modules.
Its rows and two physical births are as follows:

| Stage | Rows | New outgoing mask | H | P(2) | g |
|---|---|---:|---:|---:|---:|
| S | `[58,60,57,464,480,456,135,263,71]` | none | 3159 | 48924 | 81 |
| S+x | `[58,572,569,464,992,968,135,775,583,73]` | 73 | 13869 | 246600 | 9 |
| S+x+y | `[58,572,1593,1488,2016,1992,1159,1799,1607,1097,3]` | 3 | 51703 | 1106588 | 1 |

The first new vertex x beats exactly old vertices0,3,6. The second y beats
exactly old vertices0,1 and loses to all other existing vertices, including
x. Each birth preserves all existing labelled arcs. All three stages are
strong, and full exact ordered profiles are in `receipt.json`.

One birth can never achieve content one from S, irrespective of its
attachment mask: Theorem 1 with m=3,k=1 preserves a factor of three. The
last row has coprime H and P(2), so two births suffice. Thus **two is the
exact minimum** number of arbitrary births needed to remove all profile
content from this particular strong clock.

This also demonstrates the distinction between splitting all old modules
and consuming their arithmetic obstruction. The first birth is nonuniform
to all three triangle modules, yet the content is still 9.

This route is not a ternary construction: the desired children of 3159 are
9479 and 9481, and already the first listed birth has larger H. No claim is
made here that some other two-birth route reaches those values. Coprime
content is a necessary gate for the old intact-clock strategy, not sufficient
arithmetic coverage or a target-preserving return.

## 5. The geometry connection and the global stopping point

The sharp B and D blocks are in the common-affine-chart alternating-rank-four
class. For an attachment word c_1,...,c_l on two apexes, use K=3l+1 and points

    a=(1,0,K,0),  b=(1,0,0,K),
    x_i=(1,3i, -1 if c_i has bit2 else 1,
                   1 if c_i has bit1 else -1).

Here bit1 means the least significant bit and bit2 the next bit. The form is
`omega(u,v)=u1 v2-u2 v1+u3 v4-u4 v3`. Every pairing is nonzero and realizes
the prescribed two-apex tournament. For B use word10310 with vertex order
`(0,3,1,2,4,5,6)`; for D use word22 in its original order; for C use word2.
The verifier checks every sign as an exact integer. Full coordinates are
retained in the receipt.

The previously proved Session 192 affine symplectic substitution compiler,
used in [Session 193 rooted-clock profiles](../profile_substitution/PROOF.md), then realizes arbitrary ordered
joins of these blocks in the same rank-at-most-four class. Thus all of
Theorems 3 and 4 have rank-four realizations. The construction preserves
labelled incidences, not one fixed numerical coordinate configuration
through time. No rank-four certificate is asserted here for the specific
11-vertex strong splitter in Section 4.

Conversely the standard iterated cyclic substitutions

    S_1=C3,       S_(t+1)=C3[S_t,S_t,S_t]

are in that same affine rank-four class and contain `m=3^(t-1)` disjoint
triangle modules on their leaves. Every induced extension with H not
divisible by three must therefore add at least

    ceil(3^(t-1)/2)

vertices. If it also changes f old arcs and adds k vertices, it necessarily
satisfies `f+k>=ceil(3^(t-1)/2)`. In particular the literal ternary targets
`3H(S_t)+2` and `3H(S_t)+4` cannot have a bounded-size birth/old-arc repair
uniformly over these clocks. This is a lower bound for the specified
label-preserving operations, not an obstruction to reconstructing an
unrelated target witness.

This is the exact link to the earlier results: the complete clock profile
supplies the arithmetic invariant; the path interface measures how much a
non-module operation can remove; old-arc resets are charged by the modules
they actually touch; the symplectic class supports both the obstruction and
sharp escape family. Geometry alone does not supply the missing target.

The unrestricted `{7,21}` proof, universal prime-lane closure, and a reusable
literal two-digit return remain **OPEN**. The separately saved exploratory
root-moving operator is not used in any theorem in this packet.

## 6. Reproduction and scope

From the workspace root run

    python3 verify.py > normal.stdout
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout

The standalone verifier uses only the Python standard library and explicit
exceptions. Its finite checks include:

- three complete seed profiles, independently by permutation cuts and
  least-name path-cover decomposition;
- all 11 labelled quotient tournaments through order three, with C3 factors,
  checking every ordered-profile coordinate;
- all 512 one-vertex attachment masks of the nine-vertex strong wreath,
  checking H and P(2), with 16 full-profile and literal-path controls;
- every sharp partial-budget example through four old triangle modules and
  all four minimum primitive-completion examples;
- full profiles on the three stages of the strong two-birth escape, with all
  induced-subset path counts independently checked by enumeration of literal
  simple directed words;
- all 30 exact pairings in the three affine rank-four seed configurations.

For large joins, the verifier checks the two explicit evaluations whose gcd
certifies the upper bound, while the all-order orbit theorem supplies the
lower bound. It does not silently replace a complete-profile claim by H
alone. The universal statements follow from the group action and physical
join construction, not from the finite control range.
