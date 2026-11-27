# The third ternary chain: two exact holes, two SCC banks, and explicit local shadows

2026-11-27. Read the multiply-by-three question alongside the
[older address and atom inventory](../residues_monoid_tower_tree_20261116/MONOID.md),
[strong ternary births](../session302_ternary_ray_birth_20261113/README.md),
and [even-exponent self-converse center](../chanpak_self_converse_core_20261118/encoder/RAY.md).
The Moore 7-to-21-to-63 constructions in the preceding sessions were
cycle lengths, not tournament Hamiltonian path counts.

## Two views of the same hierarchy

Every odd H has a unique expression H=a*3^k with a odd and
3-free. Ordered by their roots, the chains start at

    1, 5, 7, 11, 13, 17, 19, 23, ...

The (2j+1)st root is 6j+1, the (2j+2)nd is 6j+5.
Every level is a dilated copy of the root set. Under
multiplication the two coordinates act independently:
(a,k)(b,l)=(ab,k+l). Thus three, the smallest nonunit,
is precisely the shift in k.

There is a complementary additive tree. The three children of
an odd H are 3H-2, 3H, 3H+2. In the address t=(H-1)/2 these
are 3t, 3t+1, 3t+2, i.e. appending the ternary digit 0,1,2.
Start at t=0, allowing leading zeros for equal-depth reads.
The multiplicative chains are the constant-middle lines; k
is exactly the number of trailing ternary 1s. Stripping them
recovers the 3-free address. At depth j below H, the full block
is the 3^j odd integers from 3^j(H-1)+1 to 3^j(H+1)-1,
centered at 3^jH. Algebraically the two-sided offset is
2 times a balanced ternary word of length j.

On the third chain the addresses are 10,101,1011,... in base
three. The first cell at 7 is (19,21,23), with attained
flanks and missing center. At the next scale the middle
cell is (61,63,65), with all three attained. This is a literal
micro/macro relation on the integer labels; it does not claim
there is a single tournament edit mapping between the two cells.

## Exact and asymptotic versions of Chan–Pak on the chains

Let R be the tournament Hamiltonian-path spectrum. Its closure
under ordinal join is a product of path sets, and a directed
triangle has exactly three paths. Thus
h(a)=min{k:a*3^k in R} (infinity if empty) gives exactly a
final segment [h(a),infinity) on every chain, and
h(ab)<=h(a)+h(b). The known forbidden values 7 and 21 and the
eight-vertex witness for 63 give h(7)=2. A factor carrying 7 in
an ordinal factorization of 63 cannot have count 7 or 21. In
particular the first hit is one strong component after removing
singletons, despite its two numerical 3 factors.

[Chan–Pak](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf)
is equivalent to h(a)=0 for every other root, not just to
finiteness of the h(a). The known atom reduction separates
ordinary prime roots and the 7p roots (with 49,63,343 already
in hand). If all prime roots other than 7 are eventually shown
attained at height zero, all the chains are occupied from
height two at latest, and roots with v_7 different from one are occupied from height zero:
for v_7>=2 use 49 and 343 to manufacture every exponent >=2.
Only the first two strata of roots 7b, with b coprime to 21,
would remain in question. More directly, h(7b)<=2+h(b),
and each attained 7p with p|b lowers the first-hit upper bound
for that whole divisibility branch when the other factors of b
are attained roots. These are uses of the inherited atom
reduction, not a new proof of the prime lanes.

For a quantitative version put
D_j(X)=#{a<=X: a odd, 3-free, h(a)>j}. The exact number of
missing odd values up to X is sum_{j>=0} D_j(X/3^j).
For example, zero relative density of missing roots is
equivalent to zero relative density of all missing odd
values. One direction follows from the j=0 term, and the
other by bounding the tail j>=J by the count of multiples
of 3^J. More generally, if delta_j is the relative density
of D_j in the root set, then the relative density of all
missing odds is sum_{j>=0} 2 delta_j/3^(j+1).
So a density-one result on the root set transfers to the
entire spectrum, but neither direction alone excludes
additional isolated gaps.

## Exact SCC factorization of the unit chain and the third chain

Keep two variables x for tournament order and z for the
chain exponent; count isomorphism classes before merging
converse pairs. Let P(x,z) enumerate non-singleton strong
components with counts 3^j, j>=1. Let Q(x,z) enumerate strong
components with counts 7*3^j. Let F_a(x,z) be the all-tournament
series on chain a, including the empty word of count one.
Uniqueness of the order of strong components gives

    F_1 = 1/(1-x-P),       F_7 = Q F_1^2,
    Q = (1-x-P)^2 F_7.

This is a two-variable packaging of the inherited SCC
counting law. It is an exact deconvolution on the whole
chain, not just a leading asymptotic. Each word on the
third chain has exactly one component whose count is
divisible by seven, with two independent pure-3 contexts.
To retain only self-converse words, let Q_sc mark only
self-converse strong centers. Such a word has that unique
special component in the center and one freely chosen
context on the left, reflected and conversed on the right:

    F_7^sc(x,z) = Q_sc(x,z) F_1(x²,z²).

The series for converse-merged classes is (F_7+F_7^sc)/2.
The nontrivial point for the micro/macro reading is that
hiding the contexts can change both the nth-order count
and its converse type, but deconvolution recovers the core
inventory coefficient by coefficient.

For concrete next strata write Q_j=[z^j]Q. The
[all-order census](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
gives P_1=x³, P_2=2x^5 and Q_2=2x^8, with the two Q_2
components a converse pair. Thus the exact ordinary order
series are

    H=63:   2x^8/(1-x)^2
    H=189:  Q_3/(1-x)^2 + 4x^11/(1-x)^3
    H=567:  Q_4/(1-x)^2 + (2x³ Q_3+8x^13)/(1-x)^3
                              + 6x^14/(1-x)^4.

For example at H=189 exactly 4 binom(n-9,2) classes on n
vertices have more than one nontrivial strong component,
interpreted as zero for n<11. None is self-converse, since
one triangle cannot balance around the central core. At
both 189 and 567 all self-converse words are single strong
centers with mirrored singleton padding: the possible
centers needed by a nonempty pure-3 context at those
exponents are absent. The Paley order-seven center and the
known order-nine 567 center witness the two cases.

## A physical dyadic family shadows both holes at every ternary precision

The [one-apex binary accumulator](../tournament_computational_interfaces_fractran_20261127/README.md)
uses x as the binary address and c=(H-1)/2 as the count
register. Appending 0 adds ceil(x/2) to c and appending 1
does not change c. Take the sparse contact word

    110 1^b 0,     b>=1.

Its 1 first and 0 last make the tournament strong: every
spine vertex lies on a directed route from the apex's first
victim to its last attacker. The two zeros contribute 2
and 7*2^(b-1), so its exact count is

    H_b = 5+7*2^b,      H_{b+3}=8H_b-35.

For every k>=1 choose b with 2^b=2/7 mod 3^k. This is a
single residue class modulo 2*3^(k-1), since 2 generates the
unit group modulo 3^k; indeed 4 has order 3^(k-1) by cubing
and the valuation of 4-1, and 2 also supplies the sign
modulo 3. Then H_b=7 and H_{b+3}=21 modulo 3^k.
For k>=2 they even have exactly the same respective
3-adic valuations as the two forbidden values. The
examples are not on the third chain, since every H_b is
5 modulo 7. At k=3 one gets H_3=61 as a physical shadow
of 7 modulo 27 and H_6=453 as its partner shadow for 21.
Increasing b by full periods gives infinitely many strong
shadows at each precision. Appending triangles ordinally
lifts their matching suffixes up the middle lines.

Thus even an explicit two-zero family can look arbitrarily
close to both holes simultaneously in the 3-adic read,
with a fixed edit of three ones between them. The actual
strong tournament outputs never hit either integer. This
is a concrete test distinguishing a finite ternary
observation from exact selection at a prescribed height;
it gives no obstruction to any unknown allowed root.

Run node replay.js for exact modular lifting, direct
Hamiltonian path counts of the sparse family and the
63-to-189 join, and finite checks of the address decoding.
