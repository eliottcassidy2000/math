# Independent audit of odd-character strata and frame families

**Verdict: the stated classifications pass.** This audit independently derives
the stabilizer and digit arguments and checks actual finite group actions. It
imports no producer code and makes no literature-priority, full-automorphism,
Ramsey-record, or circuit-lower-bound claim.

## 1. Fixed-character classification

For a nonzero residue of additive order m in Z/N, its unit-action orbit is exactly
S_m=(N/m)U(m), and its stabilizer is the kernel of U(N)->U(m). This reduction is
surjective by prime-power lifting and CRT. Thus a nonzero equivariant value exists
exactly when the character kills that kernel. An odd character then ensures that
the resulting nonzero orbit function is odd; an even character is incompatible
with oddness except at the zero function.

Every real character on U(p^a), p odd, kills the odd-order principal-unit kernel,
and the square quotient of U(p) has order two. Therefore a character on U(N) is
uniquely a product of local Legendre characters, with squarefree conductor f.
It kills the stratum stabilizer precisely when f divides m. This proves the
producer's admissible strata, independent choices and conductor formula.
Each supported orbit has phi(m) residues, so its support-polynomial contribution
is 1+2z^phi(m). This counts literal differences, not merely stratum names.

Full support includes every S_p. Hence its nontrivial odd conductor divides every
prime dividing N, forcing one prime p, conductor p and p=3 mod4. The converse is
the inherited prime-power frame theorem. The old source
`odd_prime_function_cycle_20261018d.md`, section3, already contains the 2^a
prime-power symbol count. The composite necessity and full support bookkeeping
extend that source; neither the old count nor elementary character theory should
be described as newly discovered here.

The unit stabilizer of a nonzero fixed-character symbol is exactly ker chi.
Its affine translation extension is an exhibited subgroup, not necessarily the
full automorphism group. Different characters cannot yield the same nonzero
symbol because one nonzero orbit value recovers chi(u) from its transport.
The zero symbol must be counted once in a union over characters.

## 2. Prime-power frames and the subgroup divisor

For x,y written in base p, the first differing low digit v gives
(y-x)/p^v modulo p equal to y_v-x_v. Independent nonzero digit scalings therefore
preserve the first-difference position and can change each Legendre frame sign
independently. This proves every frame is isomorphic to the all-positive
lexicographic Paley tower. It does not prove equivalence under a single ring
multiplier, and the ring Z/p^a is not the field F_(p^a).

At every node of the p-ary low-digit tree, translations and square multipliers
of the next digit act independently. The root choice and subtree choices can be
recovered from the resulting vertex permutation, so the iterated subgroup is
faithful and has order [p(p-1)/2]^((p^a-1)/(p-1)). Its action on ordered Hamiltonian
paths is free: an element fixing the entire ordered path fixes every vertex.
Thus this subgroup order divides H. No claim about the entire automorphism group
or an exact formula for H is needed. The prior ternary tower is credited to
THM-4754 and the lexicographic mechanism to THM-4755.

Independent subset-DP controls give H=3,189,3159,95095 for the tested orders
3,7,9,11, respectively. Their stated subgroup divisors are 3,21,81,55.
The 81 ternary-square automorphisms were also enumerated as literal permutations
and every arc checked, rather than inferring their faithfulness solely from a
scalar order formula.

## 3. Separate characters on separate strata

A strictly broader class allows a separate odd character chi_m of U(m) on every
additive-order stratum S_m. Full support is possible if and only if every prime
p dividing N satisfies p=3 mod4: prime strata prove necessity, while every divisor
m then has odd characters. If m has r distinct prime factors, exactly 2^(r-1)
real characters are odd, since evaluation at -1 is a nontrivial homomorphism on
the r-dimensional binary character group. Each has two independent initial signs.

Consequently the full labelled-symbol count is

    2^(sum_{m|N, m>1} omega(m))
      = 2^(sum_{p|N} a_p product_{q!=p}(a_q+1)).

There is no double counting: on a supported transitive orbit its function recovers
both its initial sign and its character. This broadening must not be conflated
with the original requirement of one global character.

For this broader full-support class, the stabilizer among unit multipliers is
exactly U(N)^2. The prime strata force every local Legendre value to be +1, and
this condition suffices for every chosen stratum character. A multiplier with
all local Legendre values -1 reverses every arc, since each odd stratum character
uses an odd number of those factors. These statements again concern unit
multipliers, not all automorphisms.

For N=21, literal enumeration gives sixteen distinct symbols and tests every
stratum covariance equation, square stabilizer, and simultaneous reversal coset.
The receipt records the exact fourth-trace distribution of the skew matrices:
15876,19908,23940,44100, each attained by four symbols. Here the isomorphism
classification **within the sixteen-symbol family** can be completed by retaining
the actual multiplier orbits, rather than inferring isomorphism from equal traces.

There are two possible characters on the order21 stratum, chi3 and chi7. For
fixed choice chi_j, the three stratum signs are epsilon3,epsilon7,epsilon21.
Unit multipliers realize arbitrary pairs of Legendre values (chi3(u),chi7(u))
by CRT. They change epsilon_j and epsilon21 by the same sign, so

    delta=epsilon21*epsilon_j

is invariant. Conversely, any choices of epsilon3 and epsilon7 can be reached,
and delta determines the resulting epsilon21. Thus there are exactly two unit
orbits for each character choice, four orbits total, each of size four. Character
choice on the order21 stratum is unchanged by multiplication, since the unit
group is abelian. Every multiplier orbit consists of isomorphic tournaments.

The independent checker explicitly forms these four orbits as sets of literal
symbols. Their four trace values are distinct. Since tr(S^4) is invariant under
any vertex permutation, no isomorphism can join different unit orbits. Therefore
there are **exactly four tournament isomorphism classes among these sixteen
stratum-character symbols**. This is not a census of arbitrary21-vertex
tournaments, and equal fourth traces are not asserted to characterize isomorphism
outside this finite classified family.

## 4. Review of the transport/eighths synthesis

The affine cycle theorem is a direct composition identity, with an abelian value
group and genuine automorphisms on every edge. Its weighted arithmetic use first
passes to Q; integrality is a further obstruction, not implied by rational
transport. The common product does not generally determine the inhomogeneous term.

The Boolean firewall is correct for n>=2: a nonzero translation partitions the
cube into pairs. Constant derivative zero gives even weight, and constant derivative
one gives weight2^(n-1), also even. An input permutation cannot turn any unbalanced
truth table into its complement. The n=1 exception is correctly retained.
The literal involution x->-x on F_2^n is the identity; no analytic oddness can be
silently assigned to it. Strong derivative obstructions are not circuit hardness.

The four eighth-character rows are orthogonal and the phase formula is correct.
Division by four requires the stated coefficient-ring hypothesis. Odd unit
differences modulo2^a connect precisely opposite parity classes, so the underlying
support graph is complete bipartite and has no odd cycle. These are exact maps and
hostiles, not a universal transfer from one meaning of oddness to another.

## 5. Reproduction and finite scope

Run `python3 check_independent.py` and `python3 -O check_independent.py` in this
folder. The outputs are byte-identical (`normal.json`, `optimized.json`).

The checker builds real characters by literal square-coset multiplication tables
for every odd N from3 through105, then builds the full signed action graph on
nonzero residues, with oddness imposed as a separate constraint. It does not use
producer orbit, conductor, or frame code. Small odd partial symbols are also
enumerated literally. Separate tests cover frame bijections, explicit tree
automorphisms, Hamiltonian subset DP, the N21 extension, and all odd-weight
Boolean tables through four variables.

The finite checks corroborate the displayed uniform proofs. They are not a proof
assistant certificate or independent human peer review.
