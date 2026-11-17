# The minimum Chan–Pak proof: atomic cores, inherited factors, and typed descent

**Proved.** Ordinal composition reduces Chan–Pak to two families of
strong classes, one for each odd prime p other than 7 and one for
each 7p with prime p other than 2,3,7. The three remaining
arithmetic atoms have literal witnesses. A separate uniform-module
theorem gives the **exact gcd** of every possible future strong
count containing a fixed block, not merely a sufficient divisibility
rule. For sevenfold visible symmetry, the first (7p) with a
prime quotient and a nontrivial automorphism is 7·593.
The quotient by that very symmetry has 5 paths, not 593. This
identifies precisely where an arithmetic tree is not a tournament
descent tree.

## 1. One strong isomorphism class per atom is enough

Let M be the positive odd integers excluding 7 and 21. This is
the maximal admissible arithmetic monoid, not an assumption that
the physical spectrum equals it. The SCCs of a tournament are
totally ordered, and a Hamiltonian path visits them as a
concatenation of full paths. Thus ordinal sums multiply H.
The exact atoms of M are:

* all odd primes other than 7;
* 7p for prime p ≥ 5, p ≠ 7;
* 49, 63, 343.

No displayed atom has a proper split in M. Conversely factor
N = 7^a 3^b times a product of other primes. For a ≥ 2, represent a
using 2 and 3 (49 and 343); for a=1 use a 7p if available, or 63
when the only other prime is 3 and b≥2. The cases b=0,1 without
any other factor are precisely the two holes. The literal strong
classes for 49,63,343 are in [the replay](atoms_replay.py).
The exclusions of [7 and 21](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
are inherited. No universal existence proof for either infinite
prime family is asserted here. The prior [corrected finite certificate](../../../session203/coverage-recertification/README.md)
supplies both families wherever the target count is at most 4,294,967,295,
but it gives no asymptotic existence theorem.

This proves two important minimality statements. First, within
ordinal multiplication no smaller numeric basis can suffice.
Second, any witness for an atom is one strong core plus ordinal
singleton padding; proving existence of just one unrooted strong
class per atom suffices. If at order k there are s_h(k) strong
classes for atom h, the full order-n fibre has exactly
sum_k (n−k+1)s_h(k) classes. The self-converse fibre keeps only
self-converse cores and equal left/right padding. It is not
necessary to classify all tournaments, prove uniqueness of first
creators, or prove that all atom cores are self-converse.

**Highly composite counts.** [ATOMS.md](ATOMS.md) gives the
exact arithmetic factorization-length equations. For example
63² = 49·3^4 has two versus five atomic factors; the largest
ratio of lengths for any fixed integer is exactly 5/2. The
possible gaps between adjacent lengths over all targets are
exactly 1,2,3. When a=0 or a≥2 and the non-7 prime factors are
already realized, the maximum number of nonsingleton SCCs
among all physical witnesses is b+sum_{p≠3,7}v_p(N)+floor(a/2).
This follows by refining each SCC count into arithmetic atoms;
equality uses known 49 and 343 and ordinary prime factors, not
a 7p theorem. Thus highly composite targets yield many
compositions but require fewer existence proofs, not more.

## 2. The exact context-independent divisor of a block

For a nonempty tournament B let c_r(B) be the number of ordered
covers by r nonempty directed paths, and let P_B(q) allow empty
slots. Let D(B) be the gcd of H(T) over all **strong**
tournaments containing B as a proper homogeneous module
(all outside contacts are uniform). Then

    D(B) = gcd_{1≤r≤|B|} ((r−1)! c_r(B)).

**Divisibility.** In a quotient word through a mark X for B,
write the r−1 internal return runs between visits to X.
They are nonempty, labelled by disjoint sets of exterior
vertices, and can be freely permuted without changing any
needed boundary arc. Thus the coefficient a_r of c_r(B) in
any quotient read is divisible by (r−1)!. This is stronger
than just taking gcd(c_r(B)): the exterior itself supplies
a factorial.

**Sharpness with strong test contexts.** Two triangular readers
C3[X,1,1] and C3[X,TT_2,1] read 3c_1+c_2 and 5c_1+3c_2.
Any common divisor is odd, so their combination 4c_1 recovers
divisibility of c_1 and then c_2. For r≥3 take r−1 pairs
(U_i,V_i) with X→U and V→X, matched U_i→V_i, all other
V_j→U_i, and arbitrary arcs within U and within V.
The quotient is strong (each pair lies on an X-triangle).
Any return consumes at least two vertices; at r visits
all helpers are consumed in matched two-vertex returns.
Hence its highest coefficient is exactly (r−1)! and
there are no higher coefficients. Subtract lower terms
inductively, proving the gcd equals the displayed integer.

This is an exact composition classification, or **divisibility
profile**. When B sits uniformly in T and T sits uniformly
in a larger tournament, B is still uniform. Therefore
D(B) divides D(T): a forced factor cannot be repaired by
nested uniform wrapping. The only way out is to change
the block or allow nonuniform contacts. D(B) describes a
guaranteed divisor, not which multiples can be placed.

### The finite symmetry types at 7

For any odd prime ell, the factorial factors and the free
permutation of the r slots make all terms r≥ell automatically
multiples of ell. The remaining c_r vanish mod ell exactly
when P_B(j)=0 for j=1,…,ell−1 (triangular binomial inversion).
The [odd-cycle packing polynomial](https://arxiv.org/abs/2307.05569)
has parity P_B(−q)=(-1)^|B|P_B(q), so it suffices to check
j=1,…,(ell−1)/2. In particular 3 needs just H, 5 needs
(H,P(2)), and 7 needs **(H,P(2),P(3))**. For ell^e the exact
cutoff is the last r with v_ell(r!)+v_ell((r−1)!)<e, since r!
divides c_r. This is a bounded test independent of B's size.

Put Z_7(B)={j in {1,2,3}:7 divides P_B(j)}. An ordinal sum
multiplies the whole slot polynomial, so its Z_7 is the
**union** of its constituents' sets. All eight subsets
occur physically; the replay gives literal profiles,
and the three singleton types generate the eight by union.
This is an exact small compositional catalogue for
hereditary sevenfold divisibility, not a seedbank for H.
The same prime-power valuations add coordinatewise when
ordinal sums are composed.

Three separating examples:

| B description | (H,P(2),P(3)) | forced gcd D(B) |
|---|---:|---:|
| rigid strong 49, rows (118,108,120,49,98,64,8) | (49,812,6363) | 7 |
| rigid strong 49, rows (2,20,17,7,105,15,47) | (49,872,6843) | 1 |
| strong 427, rows (8,25,51,116,225,195,135,15) | (427,6496,49851) | 1 |

Thus even identical H and identical trivial automorphism group
do not determine the factor a block contributes in the future.
The 427 core forces 7 in **every** at-most-two-visit exterior,
but the three-visit matching quotient above gives 125259 ≡1
mod 7. This demonstrates why the third coordinate is needed.
The asymmetric eight-vertex example (4,49,90,99,169,197,147,15)
has plain cover gcd 9 and exact future gcd 27. The universal
factorial is the explanation, not an extra graph symmetry.

## 3. Visible symmetry has its own sharp law

Automorphisms act freely on full ordered Hamiltonian paths,
so their order divides H. On a strong tournament of order n,
[Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
gives H ≥ 5^((n−1)/3). A strong p-path core at prime p≥5
therefore has n<p, and is rigid. For 7p with p≥11,
1+3 log_5(7p)<p; consequently the group is trivial or C7.
This is a restriction on the atom witness type, not a
proof of existence in either type. The same order bound and
D(B) | n!(n−1)! imply D=1 on strong prime-count cores p≥5,
and D is 1 or 7 on the 7p cores with p≥11. Visible symmetry
and inherited divisibility are different tests within those limits.

**Composition of symmetry.** Fix a prime order permutation
g, order ell, acting as automorphism. On each nontrivial
orbit choose cyclic labels in F_ell so g adds one. Fixed
vertices see each orbit uniformly. Between two moved orbits,
arcs are specified by a difference set D_ij in F_ell.
Join i,j in a graph when D_ij is neither empty nor full.
An independent shift t_i on each orbit preserves the
tournament exactly when t_i=t_j along each such edge:
a nonempty proper subset of a prime-order cyclic group
has no nonzero translation period. Thus pure translations
are exactly C_ell^c, where c is the number of connected
components of the interaction graph, and ell^c divides H.
For H=7p, p≠7, if there is a nontrivial sevenfold rotation
then its orbit interaction graph must be connected. Independent necklace symmetries are impossible
on that atom; they have to be synchronized.

For a normalizing map fixing each orbit setwise, the common
linear multiplier a and offsets t_i must satisfy
a D_ij+t_j−t_i=D_ij on each interface and preserve each
internal difference set. On a nonuniform edge the offset
difference is unique if it exists. Global offsets exist iff
these prescribed differences sum to zero around every cycle.
Consequently the affine multiplier can fail inside an orbit,
at a single interface, or through an incompatible cycle of
phases. This statement is for maps fixing the orbits setwise;
arbitrary orbit permutations are not silently included.
At ell=7 the multiplier 2 cycles (1,2,4) and (3,6,5).
If it extends with no orbit permutation, it gives an
order-three automorphism (since 1+2+4=0 mod 7), forbidden
on a 7p count coprime to 3.

**Why the first nonrigid prime quotient is 593.** Each moved
seven-cycle induces one of the two cyclic seven-tournaments:
interval (175 paths) or Paley (189 paths). Their entire
ordered-cover lists agree from c_2 onwards, starting at
c_2=1638. Their exact D values are 7 and 63. In particular
a homogeneous Paley7 never belongs to a prime 7p reading:
it contributes 63 in every uniform exterior. For two or
more moved orbits, odd-cycle-pack positivity gives
H(T)≥175²: every pack formed separately inside the orbits
is also a pack of T. More generally H(T) is at least the
product of the H values of any partition into blocks.
Equality is possible exactly when an SCC is never split:
any crossing triangle adds a positive term, and connected
components of the triangle hypergraph are precisely SCCs
(overlapping triangles give strong sets, and a vertex outside
one such component sees the entire component uniformly,
otherwise it would form a triangle crossing components).

For p<4375 there is therefore at most one moved seven-orbit.
It is a homogeneous module B. If there is no exterior,
the quotient numbers 25 and 27 are composite. If the strong
quotient has three paths it is the cyclic triangle, and the
interval read gives 175+1988=7·309, also composite.
For any strong quotient on at least four vertices, the
marked vertex X has at least two two-visit words. Start
from a Hamiltonian cycle cut as X,V,X. Along V the contacts
start with an outneighbor and end with an inneighbor.
An in-to-out transition splits V into the inner run and
a suffix after X, giving a second word. If there is no
such transition, at least two consecutive initial outcontacts
or two final incontacts occur. Move the first to a suffix,
or the last to a prefix, respectively, to get a second word.
The quotient has at least five paths, so
H ≥ 5·175+2·1638 = 4151 = 7·593. This is attained by
inflating vertex 0 of the strong four-tournament with
rows (2,4,9,3) by interval7; 593 is prime. Its group has
order 7 and the vertex quotient has H=5, not 593.

The tiny p=5 exception to n<p is also rigid: a five-cycle
automorphism orbit is a cyclic five-tournament of H=15.
Two such orbits force H≥225. A single orbit with any
exterior is a module; a strong quotient triangle through
it gives a subtournament with H≥15+2·15>35, and no
exterior gives only 15. For p=7, H=49 cannot even contain a
seven-orbit. Ordinary strong H=3 is the cyclic triangle.

The prior [Mersenne, Collatz and F_21 work](../seven_twentyone_mersenne_openai_math_20261006.md)
already distinguishes 7=111 in binary from 21=10101, with
21 a direct odd Collatz predecessor of 1. The exact shared
arithmetic is the clock 1→2→4→1 mod 7 and the factor
7→21→63 by multiplication by 3 (the first admissible
term is 63). Those are not proofs of a common descent.
The new link is precise: the order-three part of the old
Frobenius symmetry can force a compositional factor, and
that factor is incompatible with a prime 7p quotient.

## 4. What a decreasing constructor would have to do

A quick constructor exists for every allowed *composition* of
available atoms: factor N, choose the arithmetic factors in
section 1, then concatenate witnesses. This needs only
O(log N) output vertices once the atomic witnesses are
available: Busch bounds every strong order by
1+3 log_5(H). It does not provide a quick generator for
missing atomic counts. Exact search at an atom is finite
(bounded strong order), but naive labelled search costs
2^O((log N)^2); it is not a fast universal recipe.

There is a genuine decreasing tree **on witnessed graphs**:
any vertex of a nonsingleton strong core lies in a triangle,
so deleting it strictly lowers H (it removes a positive
odd-cycle term). At least one deletion has H≥H(T)/n,
since each full path can be charged to the path obtained
by deleting its last vertex, and there are n choices for
that last vertex. Thus for an atom, a physical descent
candidate always lies in [H/(1+3 log_5 H),H). But a
decreasing deletion of an *existing* witness is not a
way to produce that witness from its count.

The exact [compensated ×7 reader](../fermat_flip_stern_20261117/FACTOR_SEVEN.md)
gives a typed alternative edge: a cyclic target
ab+2uv on typed smaller profiles (u,R=a+u),(v,S=b+v)
is greater than both R and S. Indeed S≥2v gives
ab+2uv≥v(R+u)>R, and symmetrically it is >S. Thus any
successful choice decreases a numeric rank and can be
lifted constructively. A concrete inverse even selects
divisors of a shifted prescribed prime: see [ATOMS.md](ATOMS.md)
for the exact factorization and the verified mixed-core
landing p=292235682641783. Its two factors are actual
cover differences, not arbitrary divisors.

For a literal one-outgoing-edge version, take a state to be a
list of pending marked subproblems. Always work on the
largest rank first, breaking ties by a fixed order. An ordinal
factorization peels off an atom; a cyclic lift replaces one
entry by its two smaller marked profiles. The ordered multiset
of ranks strictly decreases (equivalently, the ordinal sum of
one omega^rank per entry, arranged in decreasing rank, falls).
Thus one can have a deterministic drain even when construction
has two children. A scalar edge to only one child would discard
the other prerequisite.

The minimal missing lemma is a **total target-dependent
choice of one physical typed step for each uncovered atom**.
Such steps would give a terminating tree by strict decrease.
The known scalar spectrum, visible symmetry, and the
finite compositional divisibility types remove impossible
edges, but none guarantees a child interface with the
correct absolute count. The graph can drain to 1 only
after that totality is proved; no Collatz conjecture is
needed or inferred.

Run `python replay.py` and `python atoms_replay.py`.
The optional counts in the first replay use literal
subset path counting and independent cover convolution;
the matching test quotients are explicit, not searches
over exteriors.
