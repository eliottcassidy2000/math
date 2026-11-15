# Pair-lines: a smooth-core obstruction and an exact inverse for vertical fibres

**Proved:** all-order small-core obstructions for fixed readers; order and
congruence filters for actual pair-lines; and the fixed-H catalogue and generating
function. **Finite-exact replay:** sharp line-gap witnesses and labelled cycle-packing
checks. No exclusion from the unrestricted Chan–Pak spectrum is claimed.

## 1. Interfaces and the analogies that changed the question

Write h=H, p=P₂ (two ordered nonempty directed paths), and q=p+2h (two
ordered slots, empties allowed). For an ordinal join, every ordered slot
intersects each component in one possibly empty path, and conversely those
paths concatenate. Hence both h and q multiply, including across strong
components. A singleton contributes (1,2).

A seeded sample of older paths (lexicographic Git list, GNU shuf with a
stream of the string 47) included [scale credits](../scale_credit_and_tournament_compression_20261021a.md),
[one-sided walks](../meta_fibonacci_one_sided_walks_20261013.md), and an
audit of planar line fibres. Continuing from the latter into
[carrier factor rigidity](../planar_carrier_factor_rigidity_20261019.md)
suggested retaining prime support rather than catalogue size or the
scalar projection. The scale-credit analogy distinguishes an infinite
horizon in one variable from simultaneous growth in the hidden variable.
No Jacobian or algebraic-dynamical theorem is imported.

## 2. Finite prime rank obstructs fixed two-visit readers

Consider a finite bank of fixed two-visit readers

    N = Σ_{J⊆[m]} a_J ∏_{i in J} p_i ∏_{i not in J} h_i,

with a_J nonnegative integers and a_empty≥1. Each module is a nonempty
tournament, and each of its strong components has order ≤s(X), with s≥1.
We count **distinct output integers**, not presentations or classes.

For an SCC of size j≤s, h≤j! and q≤(j+1)j!, since a permutation and cut
specifies an ordered slot pair. Thus the prime factors of each h_i and
q_i are at most Y=(s+1)s!. Every nonsingleton strong component has h≥3
(it has a triangle). Since ∏h_i≤N≤X, the total number of nonsingleton
components is ≤log₃ X.

Call p_i active at a tuple if it occurs in a nonzero summand with positive
coefficient. Active p_i≤X. Also p_i≥2^{|T_i|}−2 by choosing any nonempty
proper vertex subset and a Hamiltonian path on each of the two sides.
For inactive p_i we can discard all singleton SCCs of its module (retaining
one singleton if it was transitive): h_i stays the same and p_i cannot
increase, by multiplicativity of q. Thus no zero monomial comes alive.
This deletion can be done for all inactive modules at once, and all
surviving nonzero monomials are unchanged. We can count presentations
using at most O(log X) SCCs per module. Each prime exponent of each
coordinate is then O(log X log Y), giving the bound

    exp(O(π(Y) log(log X · log Y))).

In particular, if s(X) ≤ (1−ε) log log X / log log log X for fixed ε>0,
then Y≤(log X)^{1−ε+o(1)}, and the support up to X is X^{o(1)}. A fixed
finite library gives polylogarithmic support. Strong block permutation
or perfect modular coverage cannot change this prime-exponent obstruction.

## 3. A larger threshold with a live top term

For a positive single-module line, or an m-module reader in which a_J>0
for every J, we get a stronger bound. It also holds when only those a_J
corresponding to the permitted nonsingleton sets J are required positive.
For nonsingleton modules p_i≥h_i by cutting spanning paths once, so
q_i≤3p_i; singleton modules give q_i=2,h_i=1. The live monomial at the
exact set J of nonsingletons therefore gives ∏q_i≤3^m X, and the all-H
monomial gives ∏h_i≤X. This applies to two-visit matrices with all four
entries positive, including the inherited hyperbolic kernel. It also
applies to the single-variable signed projection αq−βh with α>0 and
β<2α, since that quantity is bounded below by a positive multiple of q.

Let B_m(Z,Y) count m-tuples of Y-smooth positive integers with product
at most Z. For 0<σ<1 the elementary Rankin inequality is

    B_m(Z,Y) ≤ Z^σ ∏_{p≤Y}(1−p^(−σ))^(−m).

Indeed, dominate the indicator of product≤Z by the σ-power weight and
sum that weight over all tuples. The logarithm of the Euler product is
O_{σ,m}(Y^{1−σ}), by summing over all integers. If
s(X)≤(2−ε) log log X/log log log X, choose
1−1/(2−ε)<σ<1/2 (it suffices to take 0<ε<1).
Then Y≤(log X)^{2−ε+o(1)} and each tuple has at most X^{σ+o(1)}
possibilities. Counting h and q tuples independently gives
X^{2σ+o(1)}=o(X) distinct outputs. Losing physical compatibility
only enlarges the bound.

This rules out positive density, and therefore universal coverage, by
slowly growing libraries of small components in a fixed fully positive
bank. It is not a bound on total order, and it says nothing against
unrestricted or target-changing readers. For the census, rank misses by
the least **maximum SCC order** of a witness, not merely its total order:
the obstruction cannot be fixed by adding enough small components.

## 4. Hard shell constraints before the line test

Use the odd-cycle packing specialization of
[Grinberg–Stanley](https://arxiv.org/abs/2307.05569). For a vertex-disjoint
packing C of directed odd cycles, let c=|C|, u count uncovered vertices,
and d=Σ_{γ in C}(|γ|−1). The ordered slot-filling polynomial is

    f(x)=Σ_C 2^c x^(c+u)=Σ_C 2^c x^(n−d).

The empty term is x^n; every other term has even d≥2. Hence f(1)=h
is odd, f(2)=q, and:
* q≡(−1)^n h (mod 3);
* for odd n≥3, q≡4·(# directed Hamiltonian cycles) (mod 8);
* for even n≥4, q≡8·(# directed (n−1)-cycles) (mod 16).

Read powers 2^{u+2c} for the last two: at valuation 2 the only contributing
nonempty type is (u,c)=(0,1), valid only for odd order; at valuation 3 for
even order use (u,c)=(1,1). The empty term vanishes in the stated moduli.
At fixed n and h distinct q values are thus separated by multiples
of 12 (odd n), 24 (even n). The least nonzero gaps do occur:
at n=5 codes 40 and 76 have (h,q)=(15,108),(15,120);
ordinal singleton padding doubles q and gives the even bound at n=6.
Codes use 1 for i→j in lexicographic pair i<j.

The same formula gives a positive-coordinate envelope, sharpenable with
the directed triangle count τ:

    1 + τ/2 ≤ q/2^n ≤ 1 + τ/2 + (h−1−2τ)/16 ≤ (h+3)/4.

Indeed an individual triangle is the only kind of packing with d=2;
all others have d≥4. The final upper bound is attained precisely when
there are no longer odd cycles and no disjoint directed triangles,
e.g. in a transitive tournament or one C₃ padded by singletons. In
particular every physical (h,q) obeys

    ceil(log₂(4q/(h+3))) ≤ n ≤ floor(log₂ q),   n≤q/h−1.

The last bound uses cuts at all n−1 internal edges plus two empty-side
cuts of each spanning path. These are necessary, not sufficient:
(h,q,n)=(3,36,5) obeys the global envelope, cut bound and mod-12
check, but does not occur. The odd-cycle expansion forces exactly one
triangle if h=3. A nontrivial strong component larger than a triangle
contains at least two: any outside vertex mixed towards a fixed triangle
produces another; if all outside vertices are uniform towards it, then
strong connectivity forces both uniform types and an arc from the
beaten type to the beating type, producing another triangle. Thus
q=12·2^(n−3) at h=3, giving q=48 rather than 36.

## 5. All-order fixed-H inverse with multiplicities

Set ρ=q/2^n. Both h and ρ multiply under ordinal join, and a singleton
contributes (1,1). At fixed count h₀, all usable nonsingleton SCC letters
have h dividing h₀ and are found through order 1+floor(3 log₅ h₀)
by [Busch](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
Every ordered word of such letters with product h₀ yields a normalized
profile; its total order w is its unpadded order. Every n≥w is obtained by
inserting singleton factors. For h₀=1 take the empty word and n≥1. So for
each (h₀,ρ) the admissible n are empty or a whole ray n≥w_min.
This classifies all orders on H=h₀ by a finite computation.

A word of k nonsingleton SCC isomorphism classes of total order w
contributes binomial(n−w+k,k) classes at each n≥w, as the singleton
components are distributed in k+1 gaps. The converse reverses words
while replacing each letter by its converse. Let [a,b] be a basis of
a formal monoid algebra, with componentwise multiplication of indices.
Write W=Σ_C z^{|C|}[h(C),ρ(C)] over nonsingleton strong classes; let W_sc
be its subseries on self-converse C and W₂ replace each term by
z^{2|C|}[h(C)^2,ρ(C)^2]. The ordinary and self-converse series are

    A=1/(1−z−W),     S=(1+z+W_sc)/(1−z²−W₂).

The latter arranges an arbitrary left half, a self-converse middle or no
middle, and a forced mirrored half. (A+S)/2 counts converse-merged
classes, including the empty tournament. This refines the
[older fixed-H census](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
by retaining the normalized second coordinate. It does not claim a
faster unrestricted enumeration.

The practical restriction is that a fixed H cannot grow a richer normalized
fibre later. At prime H there is exactly one nonsingleton SCC per
presentation. A search for absolute placement must therefore let a strong
component grow with the target, or change the reader, rather than
padding each fixed-H seed. Which of those can eventually yield universal
coverage remains open.

## 6. The content budget for prime and 7p targets

There is a cheap multiplicative line sieve beyond the size screen.
For a module define g(T)=gcd(h(T),q(T))=gcd(h(T),p(T)). Every term of a
two-visit kernel contains one coordinate of each module, so

    ∏_i g(T_i) divides N.

This holds for signed coefficient identities too. On a module assembled
from strong components C_j, for any prime r the valuation of g is exactly

    min(Σ_j v_r(h(C_j)), Σ_j v_r(q(C_j))).

In particular it can grow across components even when every individual
component is primitive; a prime dividing h of one and q of another
spoils primitivity. The proper state carries both prime budgets, not
the componentwise gcds alone. For a prime target larger than every
atom coordinate, all modules must be primitive. For N=7p and p larger
than all atom coordinates, at most one module can spend a single factor
7, and no other overlap is allowed. Irrespective of how large the
atoms are, no C₃ may occur in any module for a target coprime to 3,
because q≡±h mod 3: indeed a whole 3-divisible h fibre is unusable in
these readers. This is a restriction on the reader, not the target's
unrestricted realizability.

For a concrete cross-atom obstruction, the order-four profile (5,32)
is primitive and the order-five strong profile (11,100) is primitive.
Their ordinal product is (55,3200) with content 5. So a reader into a
fixed prime target larger than all these coordinates must refuse this
seemingly eligible assembly. The rule suggests using a conflict graph
on the palette: connect pairs whenever either h prime support meets
the other q support (mark loops for nonprimitive atoms). Restrict the
prime corridor to independent sets, then use exact exponents to handle
the 7p cases. This cannot prove coverage, but changes which source
combinations are worth attempting.

## 7. Repairing a concrete band construction: joint dyadic, ternary and quinary units

The simpler inherited bank (5,32) and (11,100) gives
N=6·32^e·100^j·2^t−5^e·11^j through the fixed physical reader 6q−h.
If e,j≥1, the previous content obstruction makes every such N divisible
by 5. In particular the e=1 branch used to switch dyadic cosets is not a
prime-compatible branch. This does not refute its geometric-band theorem:
the issue is forced divisibility in the part meant to seed hard targets.

Replace the toggle by order-eight strong profiles (223,4096) and (599,8192),
codes 10760449 and 18975246. Adjoin the order-eight toggle (217,4096), code 27280129. The four
literal primitive components (including order-five code 548 with (11,100))
have no cross-atom conflict.
Their ordinal products and padding, sent through the same physical reader
6q−h, realize

    N(d,e,f,j,t) = 6·2^(12d+12e+13f+t)·100^j − 217^d·223^e·599^f·11^j,

with d∈{0,1,2,3}, e,f∈{0,1}, j≥1, t≥0. Every value is coprime to 2·3·5·11 and to the prime factors of selected toggles:
for prime factors of the centre use the nonzero offset, and conversely
for factors of the offset use the nonzero centre. This is not a primality
claim.

**Joint band theorem.** Fix exponents r≥3 and s,v≥1.
For every unit residue a modulo 2^r·3^s·5^v and every open I⊂(1,2), the above
family has Ω(K) distinct counts in 2^K I congruent to a, for all sufficiently
large K (constants may depend on the fixed parameters).

**Arithmetic choice.** Take t≥r with t≡0 modulo M=2·3^(s−1). The dyadic
subgroup generated by 11 consists exactly of the residues 1 or 3 modulo 8
and has order 2^(r−2); this follows inductively from v₂(11²−1)=3.
Both toggle bases are 7 modulo 8, so e+f parity determines the coset.
Within that coset the parity of j is determined already modulo 8.
Since 223≡1, 599≡−1, 11≡−1 modulo 3 and 3 divides the centre, f can be
chosen to meet a modulo 3, and then e is determined by the dyadic coset.
This fixes a j residue modulo 2^(r−2). To lift modulo 3^s, let j vary
through the required parity j=j₀+2u (j₀=0 or 1). The t restriction
makes 2^t=1 modulo 3^s. In the resulting expression as a function of u,
the centre difference between u and v has valuation ≥3+v₃(u−v) since
v₃(100²−1)=2; the offset difference has valuation 1+v₃(u−v) since
v₃(11²−1)=1. Thus the full difference has valuation exactly
1+v₃(u−v). The map on u modulo 3^(s−1) is a bijection onto the residue
class already prescribed modulo 3. The new toggle satisfies 217≡1 modulo 8 and modulo 3, so the choices
of e,f did not depend on d. Take j≥ceil(v/2), making the centre
vanish modulo 5^v. Because 217≡2 modulo 5 is a primitive root there
and 11≡1, there is a unique d∈{0,1,2,3} solving the residue modulo 5.
Before lifting the ternary residue, fix this d; its presence does not
change either valuation estimate above. Modulo 5^v the remaining
condition for j is exactly one residue modulo 5^(v−1), since
v₅(11−1)=1 gives the full subgroup congruent to 1 modulo 5.
Now combine the dyadic, ternary and quinary j conditions by CRT.

**Geometric choice.** Put B_j=log₂(6·2^(12d+12e+13f)·100^j) and set
t=K−floor(B_j), so the scaled centre is 2^{frac(B_j)}.
On the permitted progression for j, B_j modulo M is an irrational
rotation because log₂100 is irrational. Restrict B_j modulo M to
the interval with integer part K modulo M and fractional part in
an interior closed subinterval of log₂ I. This ensures the former
constraint t≡0 mod M and leaves a geometric margin. The interval has
positive length. Any segment of the progression of sufficiently many
consecutive indices hits that interval: a finite initial rotation net,
translated to the current segment, still has the same mesh. There are
only finitely many K residues modulo M, so a single mesh works. Among
positive j≤K/(4 log₂100) there are Ω(K) hits. For each, t≥3K/4−O(1)
and the offset is O(2^{K/4}), so the actual count lands in 2^K I.
If two such j give the same count, the power 2^min(t,t') divides
11^{|j−j'|}−1, impossible for large K: its valuation is 1 for odd
difference and 2+v₂(|j−j'|) otherwise. This proves distinctness.

The repair matters because it lifts geometry and three local arithmetic
coordinates without ever reintroducing the forced factor 5. Nevertheless,
the entire family still has O((log X)^2) distinct counts up to X:
for N≤X, its centre is at most 6X/5 by h≤q/2, bounding j,t by O(log X).
The absolute window still needs a changing core or reader to fill it.
