# Cycle-cut obstruction for strong black deficit two

Status: PROVED scoped obstruction, independently audited; not a complete strong near-purity classification. 2026-10-20.

Let T be a self-converse strong tournament, h=H(T)/|Aut(T)| its Hamiltonian-path orbit count, and b its blue orbit count. Blue means that positional reversal of the path is an anti-automorphism of T. Assume h-b=2.

## 1. Exact cycle-cut orbit lemma

Fix any directed Hamiltonian cycle C=(v_0,...,v_{n-1}), indexed modulo n, and write c(v_j)=v_{j+1}. Let

    D = Aut(T) intersect <c>,   d=|D|,   q=n/d.

Every tournament automorphism has odd order: an even-order element has a nonidentity involutory power, and an involution reverses the arc of any transposed pair. Thus d is odd. As a subgroup of the cyclic group <c>, D=<c^q>.

Cutting C before v_i gives P_i=(v_i,...,v_{i-1}). Any automorphism taking P_i coordinatewise to P_j is necessarily c^(j-i). Consequently these n cut paths meet exactly q Aut(T)-orbits, indexed by i modulo q. No possible automorphism outside <c> can identify two additional cuts: the full ordered words uniquely specify that map.

The positional reversal of P_i is

    r_i(v_j)=v_(2i-1-j).

If at least one r_a is an anti-automorphism, then r_i is an anti-automorphism exactly when r_i r_a=c^(2(i-a)) is an automorphism. Thus the blue cut-orbit count is exactly gcd(2,q). If none is anti, it is zero. In either case, the number of black cut orbits is at most the global number h-b=2.

The blue criterion here is literal positional reversal, not an unjustified requirement that a chosen global converse identification fix the path. Equivalently, a path orbit fixed by the induced converse/reversal action has anti-automorphic positional reversal: the required orbit-identifying automorphism composes with the chosen anti-isomorphism to give precisely that reversal.

## 2. Odd order

For odd n, q is odd. If no cut is blue, the bound q<=2 forces q=1. But q=1 makes c an automorphism, identifying T with a circulant in the cycle coordinates; every map j -> a-j is then an anti-automorphism. Therefore the no-blue case is impossible.

Hence q-1<=2, so each cycle has either:

- q=1, with its unique cut-path orbit blue; or
- q=3, with one blue and two black cut-path orbits.

In the second case 3 divides n. Moreover, all q=3 cycles belong to a single Aut(T)-orbit: each uses both of the only two black path orbits, and a closing path determines its directed Hamiltonian cycle uniquely. If such a cycle exists, every bypass path is blue. Here a bypass path has its first vertex beating its last; a closing path has its last vertex beating its first. Both conditions are preserved by automorphisms and by converse/reversal.

In particular, if 3 does not divide n, the successor of EVERY Hamiltonian cycle is an automorphism, and all cycle cuts are blue.

## 3. The prefix obstruction needs only cycle-cut purity

Source: ../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md, especially sections 2, 4, and 5.

That proof gives, in every odd circulant of order n>=7, a directed Hamiltonian path whose positional reversal is not an anti-automorphism. Crucially, each displayed witness is also a CLOSING path:

- F1 starts at 2 and ends at n-1; the return increment is 3, in S by the F1 hypotheses.
- F2 starts at k-1 and ends at n-1; the return increment is k, in S by the F2 hypotheses.
- Both F3 variants start at 1 and end at n-1; the return increment is 2, in S by the F3 hypotheses.
- The unit-step coordinate changes in section 4 preserve cycle closure and the violated-arc certificate.
- The four order-seven witnesses in section 5 have return increments 4,4,3,2, respectively, each in the displayed connection set.

The prefix lemma itself only needs the standard cycle cut and the displayed rerouted cycle cut to be blue. The three-pattern and unit-reindexing arguments therefore run unchanged under the weaker assumption that every Hamiltonian-cycle cut is blue. They do not need blue reversals of arbitrary bypass paths.

Consequently no odd circulant of order at least seven has every Hamiltonian-cycle cut blue.

Combining this with section 2 proves:

    A strong self-converse tournament with h-b=2 and odd n>=7
    must satisfy 3 | n.

This is an all-order exclusion in the complementary odd congruence classes, not a full classification. The small orders must be handled separately; no claim about all strong order-five tournaments is inferred from the circulant witness theorem.

## 4. Even order: exact surviving alternatives

For even n, q is even, and the blue cut-orbit count is either zero or two. Since q minus that count is at most two, q is either 2 or 4.

Since d is odd, n and q have the same power of two dividing them. Therefore:

- If 8 divides n, deficit two is IMPOSSIBLE.
- If n is 2 modulo 4, every cycle has q=2. Its two cut-path orbits may be BOTH BLUE or BOTH BLACK. Rotation by two does not itself force a reversal to be anti. Claiming that all cuts are blue in this case would be a gap.
- If n is 4 modulo 8, every cycle has q=4, with exactly two blue and two black cut-path orbits. Every Hamiltonian cycle therefore contains both global black orbits, so all Hamiltonian cycles lie in one Aut(T)-orbit. Its stabilizer has size d=n/4, giving HC(T)=4|Aut(T)|/n when cycles are counted up to cyclic rotation. Every bypass path is blue.

For n=2 modulo 4, if any all-black-cut cycle exists then all such cycles lie in one Aut(T)-orbit and every bypass path is blue. If no such cycle exists, all black paths are bypasses. These alternatives are not excluded here.

## 5. Universal lower bounds beyond deficit two

The cycle-cut argument also gives a general lower bound for any strong
self-converse tournament of even order n, without assuming h-b=2:

    h-b >= 2^v2(n)-2.

Indeed d is odd, so q>=2^v2(n), and at most two cut-path orbits are blue.
For odd n>=7 let p(n) be the least prime divisor of n. There is at least
one Hamiltonian cycle whose successor is not an automorphism: take any
cycle; if its successor is an automorphism, this gives a circulant
presentation, and the closing-prefix obstruction supplies a cycle with
a nonblue cut. That latter cycle cannot have automorphic successor,
because the inverse reflection of a circulant is always anti.
For a nonautomorphic successor q>1, and q divides n. Hence

    h-b >= q-1 >= p(n)-1.

The first inequality allows the possibility of zero blue cut orbits.
In particular a prime-order strong self-converse tournament at n>=7
has black deficit at least n-1. These inequalities are lower bounds,
not claims that equality is always attained.

## 6. Honest remaining gap

The preceding argument leaves possible strong deficit-two atoms at odd orders divisible by three and at even orders congruent to two or four modulo eight. It does not prove that the unique strong S4 is the only strong deficit-two atom. Finite census exclusions in these residual classes are separate evidence and cannot close this all-order gap. The pure-blue atom theorem alone is insufficient: one must rule out a single exceptional reflection-paired black path orbit pair.
