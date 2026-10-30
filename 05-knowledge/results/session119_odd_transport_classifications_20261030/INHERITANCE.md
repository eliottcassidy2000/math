# Session 119 inheritance: oddness is an action, not a shared adjective

**Recovered source statements, with their stated status and limits.** This is a
bounded inheritance pass, not a replay of every cited computation or a claim of
new solutions to Collatz, LRC, Ramsey, or unrestricted circuit lower bounds.

## 1. Keep four meanings separate

1. An **odd function** satisfies f(iota x)=-f(x) for a specified involution iota
   and a specified value group. Fixed points force 2f(x)=0, not always f(x)=0.
   Over a field of characteristic two the sign distinction disappears.
2. An **odd-weight Boolean function** has an odd number of ones in its truth
   table. This is not analytic oddness: on F_2^n, x=-x for every x. Negation on
   the Boolean cube needs a different named action, for example complementing
   every coordinate, and its covariance is an additional property.
3. An **odd cycle** has odd length. Directed tournament cycles, undirected graph
   cycles, cycles of a permutation, and arithmetic periodic orbits are different
   objects. Preserve the actual cycle word and orientation when relevant.
4. An **odd integer** is a unit modulo every power of two, but need not be a
   unit modulo an arbitrary odd modulus. Arithmetic action theorems must retain
   gcd/valuation strata, cycle order and the chosen multiplier.

## 2. Closest proved bridges

### Tournaments and anti-invariant sign functions

`01-canon/theorems/LEM-004-tournaments-are-odd-functions.md` proves a bijection
between circulant tournament connection sets on Z_n, n odd, and signs on the
nonzero differences satisfying sigma(-d)=-sigma(d). There are 2^((n-1)/2)
such symbols. For even n the self-opposite difference n/2 makes a sign assignment
impossible. General labelled tournaments are odd signs on ordered off-diagonal
vertex pairs, not necessarily functions on a difference group.

Its multiplicative prime-field classification gives the Legendre sign exactly
when p=3 mod4. Its higher-walk caveat is essential: oddness supplies cancellation
and the universal cherry term, not quasirandom spectral decay or a universal
Paley extremum. The rotation family is a canonical hostile to those overclaims.

### Translation, reflection, signed incidence, and odd-sheet descent

`05-knowledge/results/odd_prime_function_cycle_20261018d.md` is the closest
already-merged odd-number/odd-function/odd-cycle source. Its §2 already classifies
f(x+a)=-f(x) on Z_m with values in Z_q. If g=gcd(m,a), L=m/g, the number of
solutions is q^g for L even and gcd(q,2)^g for L odd. Reflection x->-x instead
has gcd(m,2) fixed points and a different solution count. Do not rediscover these
homogeneous formulas as new results.

Section 3 classifies all quadratic-unit-covariant signs on nonzero residues
modulo p^r: sigma(p^v a)=epsilon_v chi_p(a), one independent sign per valuation
layer. They are odd precisely when p=3 mod4, giving 2^r symbols. Nonunits cannot
be discarded: naive reduction modulo p leaves unoriented differences.

Section 4 already gives the signless-incidence kernel count q^b gcd(q,2)^o,
where b/o count bipartite/nonbipartite components. Section 5 proves odd-sheet
matching-cover descent for **fixed labels and uniformly replicated pins**, not
unrestricted colouring or multiplicity preservation. Section 6 separately proves
odd-degree square-class descent by the field norm. These are precise different
maps, not a single unrestricted transfer principle.

### Odd-cycle collections and two-adic Hamiltonian counts

`01-canon/theorems/THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md`,
building on `THM-002-ocf.md`, gives H(T)=sum_k 2^k alpha_k, where alpha_k counts
collections of k vertex-disjoint **literal directed odd cycles**, modulo cyclic
rotation. Its congruence modulo2^j retains only k<j.

Two current corrections are mandatory:

- The Oct26 correction in THM-466 says evaluation at2 does not recover the
  coefficient vector. At n6, profiles (1,2,1) and (1,4,0) both yield H=9.
- The Oct26 repair to `THM-209-independence-polynomial.md` forbids collapsing
  cycles to vertex supports. A five-vertex example has four triangles and three
  different directed5-cycles: literal conflict graph K7 gives H=15; collapsing
  supports gives K5 and the wrong answer11.

`01-canon/theorems/THM-4246-sun-reflection-orbit-two-adic-tower-and-odd-cycle-firewall.md`
gives a different expansion sum_d 2^d beta_d for coordinate-reflection fibres.
The native orbit graphs are hypercubes and are bipartite, so tournament odd-cycle
collections do not transfer to them merely because both formulas use powers of2.

## 3. Arithmetic and phase sidecars

`05-knowledge/results/collatz_mod6_20260917_floor_sums_odd_functions.md` proves
an exact floor-sum law for odd residue functions with the zero-fibre correction
Z_f(n). It expressly denies a found object-level map to Rédei cycles or Collatz
cycle counts. The even-modulus fixed point is retained correctly. Hostiles include
f=x^2 modulo3 and the false unrestricted formula for Z_7 at n=256.

`01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md`
already classifies -1 in <2> modulo odd m: all prime factors must have even
ord_p(2), with one common positive v_2 of that order. Thus3 and5 each alternate,
but15 does not;7 never alternates. The integer floor/ceiling maps are related by
translation, but that translation changes parity and does not conjugate the
accelerated odd-integer Syracuse maps. Negation exchanges the plus/minus maps.

`01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md`
retains unit and nonunit multiplication orbits for prime powers and connects
Wieferich order lifting to counts of invariant tournament symbols. The odd-order
multiplier count itself is inherited; a composite or inhomogeneous extension
must say what additional classification it establishes.

`05-knowledge/results/cycle101_tournament_collatz_dictionary_20261029.md` and
THM-4849 supply a real circulant/rational-Collatz dictionary with markings,
parity words and denominator cancellations retained. Two registered predictions
failed; the all-shape generalization is not automatic. A rational periodic word
is not evidence that every positive integer enters a cycle.

Session115's packet `session115_odd_eighth_classifications_20261030/` provides
complete phase, affine-group and one-reset-controller families. Its word-address
transport preserves positivity/products only under a congruence, and3->5 first
fails at the odd eighths. The same fractions can be different observables and
measures. Session117 proves exact costs for its single-point Boolean defects:
at most2n-2 unrestricted B2 gates, despite strong Walsh geometry. Session118
requires retained shared-DAG ancestry and coordinate maps; neither parity labels
nor free input-GL changes yield a circuit-complexity improvement.

## 4. Canonical hostiles and false shortcuts

- K4 modulo3: signless incidence has full rank4, while its Gram matrix has
  rank3. Real positivity cannot be imported into modular arithmetic.
- A triangle's double cover is a6-cycle: even covers can erase an odd obstruction.
- A genuine3-cover of K7 is3-colourable while the base is not6-colourable:
  unrestricted labels do not satisfy the fixed-fibre descent hypothesis.
- `THM-1380-free-Z2-but-no-odd-map.md`: LRC gap functions are even; free symmetry
  on maximizers does not supply the odd map required for Borsuk-Ulam. All-odd
  speed sets have the fixed maximizer t=1/2 and are easy, not hard LRC cases.
- `THM-3187-central-sign-parity-quotient-and-odd-observable-necessity.md`: every
  even-degree observable loses a central sign; a nonzero odd observable or fixed
  charged reference is needed. Necessity is not construction.
- `odd_symmetry_beyond_sylow_and_affine_circulant_rules_20261019.md`: all odd
  p-group results cannot be promoted to every odd-order group. A C15 action of
  type(3,5,1,1) has32 tournament types versus31 sign-even graph types.

## 5. Useful next theorem targets

These are proposed extensions to audit, not claims that major conjectures follow.

1. **Inhomogeneous signed-cycle equations.** Extend the inherited homogeneous
   alternating law to complete solvability, integer cokernels/Smith forms and
   explicit reconstruction for f(Px)+/-f(x)=b(x). Cycle sign holonomy, rather
   than the adjective odd, should classify the obstruction and its2-torsion.
   Multiplication actions then supply arithmetic cycle counts by gcd strata.
2. **Complete composite-modulus covariance.** For odd N, classify signs satisfying
   sigma(ud)=chi(u)sigma(d) on a specified collection of gcd strata. Stabilizers
   should force chi to descend to U(m) for each retained m=N/gcd(d,N). If all
   nonzero strata are retained and N has two different prime factors, descending
   simultaneously to each prime factor forces chi trivial by CRT, incompatible
   with oddness. This suggests a sharp prime-power-only classification extending
   the existing p^r theorem. Check independently before promoting.
3. **Keep the marked obstruction, not only its evaluation.** A common compiler
   can carry a cycle product/sign, a zero-fibre or a coordinate-change cost and
   provide reconstruction plus hostile tests. It must not identify directed
   cycle multiplicity with reflection-orbit dimension or spectral pseudorandomness.

The unifying mechanism worth pursuing is a specified action with a binary
character and a retained stabilizer/fixed-fibre obstruction. This already has
several exact consumers; none alone proves a general Collatz, LRC, Ramsey or
unrestricted circuit theorem.
