# Arithmetic decorations and the structure behind tournament quasirandomness

2026-10-29. **PROVED elementary deductions from CITED classical inputs;
FINITE-EXACT controls; OPEN extensions.** No literature-priority claim.

## 1. Fixed square-class decorations preserve prime-prefix quasirandomness

Let a_1,...,a_m be fixed positive integers. Let Q(a) be the vector of prime
valuations modulo two, and let M have columns Q(a_i). Its rank is r over F2.
For primes q congruent to 3 modulo 4, excluding divisors of the a_i, define

    c_i(q) = (1 - (a_i/q))/2 in F2.

**Theorem.** The possible colour vectors form C = image(M^t), with each
colour having limiting relative density 2^(-r) among these primes. For
any fixed labelled k-vertex tournament F and any specified colours
c_1,...,c_k in C, the proportion of ordered distinct prime k-tuples up to
X with those colours and Legendre tournament F tends to

    2^(-binom(k,2)) * 2^(-rk).                              (1)

Infeasible colours never occur. The prime-edge rule is p->q iff (p/q)=1;
quadratic reciprocity makes it a tournament on primes 3 mod4.

**Proof of the colour law.** Choose the distinct primes ell dividing the
a_i. The basic signs (ell/q) are independent uniform signs in the fixed
residue-class limit. Indeed, for an odd ell, reciprocity changes (q/ell)
only by the fixed sign (-1)^((ell-1)/2), since q=3 mod4. Exactly half the
units mod ell have each Legendre sign. If ell=2 occurs, q=3 or7 mod8 gives
the two signs. CRT combines all choices independently; the prime number
theorem in these FIXED arithmetic progressions gives their densities.
The map from basic bits z to colours is M^t z. All its fibres have the
same cardinality, proving the claim. No independence is asserted for
linearly dependent characters.

**Proof of the edge law.** Inherit [P1](../tournament_open_problem_routes_20261029/PROOFS.md):
for the signed prime matrix K, ||K||_2 <<_epsilon X^(1/2+epsilon), while
N=number of admitted primes is asymptotic to X/(2 log X). Put W=(J+K)/2,
with diagonal1/2. In the labelled-pattern product include the indicator
1_{c(q_i)=c_i} for each vertex. Telescope its binom(k,2) edge factors,
replacing each by1/2. After fixing the other k-2 vertices, the changed
factor is multiplied by f(x)g(y) with f,g in[0,1], even with the colours.
Its normalized contribution is bounded by ||K||_2/(2N). Thus

    pattern density = 2^(-binom(k,2)) product_i mu_X(c_i)
                      + O_k(||K||_2/N + 1/N).              (2)

The final term removes collisions; the finite omitted prime set costs
O_k(1/N). Taking the colour limits proves(1). This is a standard synthesis
of fixed-modulus PNT, CRT and the quadratic large sieve, with its joint
statement made explicit here.

**Concurrent credit and stronger bound.** Incoming966ab75d90 independently
proves and audits the same coloured-rooted theorem in
[Session107](../session107_quasirandom_arithmetic_20261029/characters/CHARACTER_STRUCTURE.md).
Its general fibre estimate is stronger than our fixed-colour formulation:
for sampling from arbitrary nonempty S_i of sizes n_i, the iid edge error
is at most(||K||_2/2) sum_(i<j)(n_i n_j)^(-1/2), with explicit overlap
collision correction. It also treats negative bases using an affine code.
These shared conclusions are credited jointly, not presented as exclusive
discoveries of this packet. The sigma-bank duality below is its additional
concrete application.

**What this adds.** Each fixed feasible colour class induces a quasirandom
tournament. Edges between finitely many such classes also have all fixed
random pattern frequencies. Taking the a_i to be distinct fixed primes
3 mod4 gives a named-root extension theorem: every attachment pattern to
those roots occurs at its expected density, jointly with any tournament
on the free vertices. This uses CRT/PNT in addition to P1's norm bound;
P1 alone did not assert a named-root estimate.

It does not handle a growing root bank, growing conductor, arbitrary thin
path-defined classes, exact finite response matrices, or complement/re-marking
closure. In particular it is not a uniform theorem for growing Sylvester
anchors. In every fixed nonempty colour class, the uncoloured twelve-vertex
81/27 hostile classes still have frequency ratio3: their orbit sizes differ
by3. See the inherited certificates linked in P1.

## 2. The six labels and the Hamming bank give literal dual codes

Recover the [incidence chart](../incidence_attractors_20261028/THEORY.md),
whose arithmetic atom bank is B={4,2,5,3,7,13}. There are two different
observers; they must not be confused.

* Direct tests (b/q): five independent prime coordinates2,3,5,7,13,
  while(4/q)=1. They give32 possible vectors.
* Divisor-sum tests (sigma(b)/q): the arguments are7,3,6,4,8,14,
  with square classes7,3,6,1,2,14. Their rank is3, hence exactly8 vectors.

Index the second set of bits by b. Its complete relations are

    c_3=0,    c_2+c_5+c_7=0,    c_4+c_7+c_13=0.             (3)

Each of the eight solutions has limiting prime density1/8. Modulo168,
the24 unit classes congruent to3 mod4 realize them exactly three times
each. Adding atom11 gives sigma(11)=12, hence c_11=c_2: this bit is free
because its square class is already represented. The three independent
bits may be chosen as c_2,c_4,c_7. Both constants and dependencies survive
conditioning; they are not random independent extra bits.

For the static task of answering every named character query in this bank,
three bits are also necessary: all eight distinct answer vectors actually
occur. This is minimality for that observer, not for arbitrary tournament
operations or future arithmetic updates.

For the inherited seven-prime bank {2,5,7,19,29,59,89}, the square classes
of p+1 are respectively3,6,2,5,30,15,10: all seven nonzero vectors of F2^3.
The old square-total subset code is ker M, the Hamming[7,4,3] code.
The prime-character pattern code is

    image(M^t) = (ker M)^perp,

the simplex[7,3,4] code: zero and seven weight-four words. Its eight words
occur with limiting frequency1/8. Every pair of distinct coordinates is
independent, but each Fano-line triple has parity zero identically. Thus
pairwise randomness and exact three-way dependence coexist in an actual
prime family. Equation(1) says the resulting constrained decorations still
coexist with quasirandom tournament edges.

This is an exact pairing between two operations on square classes: selecting
columns whose product is a square, and evaluating a quadratic character on
each column. It extends the repo's binary-matroid mechanism, not ordinary
graph cycle space: the Fano matroid is non-graphic. The identification with
square sigma(product b) additionally requires pairwise coprime atoms. The
old illegal selection{2,4} remains illegal; passing to characters does not
repair multiplicativity of sigma on noncoprime products.

## 3. Sylvester numbers give a completely nonrandom Jacobi tournament

Let s_1=2 and s_(j+1)=s_j^2-s_j+1. Induction gives

    s_(j+1)=1+product_(i<=j)s_i.

Hence the terms are pairwise coprime, and s_j=1 mod s_i for i<j. Every
s_i for i>=2 is3 mod4. General Jacobi reciprocity therefore gives

    (s_j/s_i)=1,       (s_i/s_j)=-1,       2<=i<j.            (4)

Orient s_i->s_j when(s_i/s_j)=1. This is exactly the reverse linear order:
every later term beats every earlier term. It is transitive at EVERY size,
including composite terms; no squarefreeness or prime infinitude is assumed.
Its cyclic-triangle density is zero, rather than the quasirandom value1/4.
Its two half-intervals have signed cut magnitude floor(n/2)ceil(n/2),
giving ||K||_2 >= sqrt(floor(n/2)ceil(n/2)), a linear obstruction.

This is a mechanism, not just a counterexample: every newly encountered
divisor becomes a permanently frozen future residue1. The recurrence pays
for Euclid coprimality by enforcing correlations. The terms grow doubly
exponentially, so their count at height X is O(log log X). The prime-only
subsequence, if infinite, is transitive too; its sparsity is outside P1's
large-subset guarantee. A Jacobi sign+1 for a composite modulus must not
be described as an iff test for being a quadratic residue.

The familiar reciprocal identity also follows exactly:

    sum_(i<=j)1/s_i = 1 - 1/(s_(j+1)-1).

It is a conservation law, not a mixing law. This Sylvester NUMBER sequence
must also be distinguished from the skew-Sylvester HADAMARD doubling tower.
Incoming [Euclid-offset syndromes](../session107_quasirandom_arithmetic_20261029/conditioning/EUCLID_OFFSET_SYNDROMES.md)
independently proves the prime-term transitivity and supplies the more
general factor-sign relation sum_j(e_j mod2)v(q_j)=v(c) when N=c at each
root modulus. Our Jacobi formulation extends transitivity to all odd terms.
The incoming13/139 table also records the direction guard: the factors
have equal signs when tested modulo3,7,43, but reversing numerator and
denominator makes their columns opposite because13 and139 differ modulo4.

### Why 13 and139 occur

1807=Phi_6(43)=13*139. If a prime p!=3 divides Phi_6(a), then a has exact
order6 modulo p: (a+1)Phi_6(a)=a^3+1, and a=-1 would force p=3.
Consequently p=1 mod6. Both13 and139 pass, and43 has order6 modulo each.
Equivalently4Phi_6(a)=(2a-1)^2+3 gives the Eisenstein discriminant-3.
For these two moduli the Sylvester orbit hits0 and then stays1. Splitting
alone does not force a hit:19=1 mod6 has the nonzero cycle below.

## 4. The19 cycle has a proved all-level clock

Put f(x)=x^2-x+1 and retain the original starting point2. Modulo19,

    2 -> 3 -> 7 -> 5 -> 2.

For G=f^4,

    G(2)-2 =1805=5*19^2,
    G'(2)=3*5*13*85=16575=7 mod19.

Using the reduced orbit gives the same derivative residue:
3*5*13*9=1755=7 mod19. The two integer products are NOT equal; only the
residue is used. The appearance of1755=(3511-1)/2 here does not identify
this nonlinear derivative with the doubling action at3511.

Write a lift of2 as2+19t. Modulo19^2 its return digit is t->7t.
The element7 has order3 modulo19. Therefore the full76-point preimage
of the four-cycle decomposes into one4-cycle and six12-cycles. The seed2
lies on the4-cycle at this level. At the next level,2+19^2 t has return
digit t->5+7t; its unique fixed digit is15, giving the4-periodic lift5417
modulo19^3. The integer seed2 instead has period12 there.

**All-size period theorem.** The orbit of2 modulo19^k is purely periodic,
of exact period

    L_1=L_2=4;       L_k=12*19^(k-3) for k>=3.              (5)

**Proof.** Exact modular iteration gives

    f^12(2)-2 =17*19^3 mod19^4,
    (f^12)'(2)=172 mod19^2, hence1 mod19.

These are short reproducible integer certificates, checked without floating
arithmetic. The minimal period at k=3 is12 by the affine digit calculation.
For completeness, the lifting step is elementary. If an integer polynomial
H satisfies H(x)=x+p^n u with p not dividing u, n>=2 and H'(x)=a=1 mod p,
Taylor iteration modulo p^(n+2) gives

    H^p(x)-x = p^n u*(1+a+...+a^(p-1)) mod p^(n+2).

Quadratic error terms vanish since2n>=n+2, while the sum is p mod p^2
for odd p. Thus the displacement valuation increases from n to n+1.
On the next digit the return acts as a nonzero translation, so its exact
period multiplies by p. The derivative stays1 mod p. Apply this repeatedly
to H=f^12, x=2, p=19, n=3. This proves(5) and pure periodicity. QED.

The unique19-adic4-periodic point above2 exists by Hensel lifting because
G'(2)-1=6 mod19 is invertible. The natural integer2 is not that point.
This cleanly separates an infinite compatible periodic lift from one fixed
integer source whose observed periods keep growing.

### An actual shared Collatz modulus, with different carries

The identity1805=2^11-3^5 recovers the earlier
[prime-power carry hostile](../collatz_local_prime_power_gates_20261020/LOCAL_CARRY_AND_PRIME_POWER_GATES.md):
the primitive parity word11101000100 has B=475, D=1805 and reduced
denominator19. Every prime of D divides B, but19^2 does not. The inherited
[circulant dictionary](../cycle101_tournament_collatz_dictionary_20261029.md)
also puts the two Paley11 parity cycles at denominators361 and1805.
These are literal consumers of the same modulus. No intertwining between
f and Collatz is supplied by the integer identity. Their carries, sources
and clocks differ, and the example specifically shows why depth matters.

## 5. Paley cancellation, Wieferich delay and Fermat rigidity

For a prime q=3 mod4 let K_xy=chi(y-x) on F_q. The classical character
correlation gives

    KK^t=qI-J,    K1=0,    ||K||_2=sqrt(q).                 (6)

Here is an elementary proof of the off-diagonal correlation. Translation
and scaling reduce it to sum_t chi(t(t-1)). The number of solutions to
y^2=t(t-1) is q plus this sum. Rewriting as
(2t-1-2y)(2t-1+2y)=1 gives exactly q-1 solutions, one per nonzero first
factor. Thus the sum is-1; the diagonal is q-1. This is exact orthogonality
inside one field, unlike the large sieve over many varying prime moduli.

At3511, ord(2)=1755=(3511-1)/2, so doubling generates the squares and
preserves the Paley tournament. It is nevertheless Wieferich. At1093,
ord(2)=364 and-1=2^182 modulo1093^2: doubling cannot preserve a tournament
connection set. The odd-order subgroup used in THM-4644 is essential.

Recover [the split/grow law](../fermat_wieferich_order_20261014e.md), not a
new theorem here. For d=ord_p(2), s=v_p(2^d-1),

    ord_(p^k)(2)=d*p^max(0,k-s).

Both1093 and3511 have s=2. Their special property is one delayed lift,
not absence of Paley mixing. The Paley operator in(6) and multiplication
on prime-power units are different operators. At a Fermat prime
p=2^(2^m)+1, d=2^(m+1), 2^d-1=p(p-2), so s=1. No Fermat prime is
Wieferich; this classical exclusion is credited to Crandall–Dilcher–Pomerance.

All circulant tournaments at a Fermat prime have only translation
automorphisms (inherited THM-4644 / its full slug in AUDIT.md), but this
rigidity does NOT imply randomness. The rotational connection set
{1,...,(p-1)/2} has two intervals of size floor(p/4) with every cross arc
pointing the same way. This is a macroscopic cut obstruction regardless
of its multiplier stabilizer. Conversely, the highly symmetric Paley
family satisfies(6). Symmetry size is not the mixing criterion.

There is also an exact comparison of the number sequences. Fermat numbers
F_j=2^(2^j)+1 satisfy product_(i<j) F_i=F_j-2, so later terms are2 modulo
earlier ones, whereas later Sylvester terms are1. For i>=2, F_i=1 mod8,
hence(F_j/F_i)=(2/F_i)=1. Reciprocity is symmetric on these1-mod4 moduli;
the result is not a tournament. All odd Sylvester terms are3 mod4, while
F_j=1 mod4 for j>=1: the two number sequences intersect exactly at3.
This assertion needs no classification of their prime terms. Also1093=1
mod4 is outside our prime-prefix tournament;3511=3 mod4 is inside it.

## 6. The Monster supplies a third mechanism through representation size

Gowers's **CITED** quasirandom-group theorem supplies a real connection.
If every nontrivial complex representation of a finite group G has degree
at least d, then its subset convolution operators have bounded nonconstant
singular modes. In particular, for A,B,C subsets of G,

    |A||B||C| > |G|^3/d  =>  some a*b=c.                    (7)

Use [Gowers, Lemma3.2 and Theorem3.3](https://arxiv.org/html/0710.3877).
The mechanism is that a nonconstant singular eigenspace is a nontrivial
representation: its dimension forces its energy to be spread out.
For the Monster d=196883 (Griess, source in AUDIT.md). The Baby Monster
has d=4371, also read from the small-degree2.B character table and descent
of its4371 constituent. These are quantitative finite statements, not an
asymptotic limit obtained by repeatedly enlarging one fixed sporadic group.

These operators are bipartite Cayley incidence operators. They are not
Monster tournaments: any nonidentity involution in a tournament automorphism
group would swap the endpoints of some arc, impossible. Monster and Baby
Monster have even order. The inherited7A normalizer quotient does act on
Paley7 by orientation OR converse, with its common C2 sign retained; it
does not make the whole Monster an automorphism group of that tournament.

The concurrent [Session107 group audit](../session107_quasirandom_arithmetic_20261029/groups/GROUPS_AND_RECOVERED_PRIMES.md)
independently reaches these conclusions and sharpens a useful boundary:
AGL(1,7) and7:3 have nontrivial cyclic quotients, hence complex minimum
degree1. The actual Monster-local-to-Paley map therefore does not transport
the Monster's group-quasirandomness parameter. Paley cancellation is a
property of its chosen connection set, even over an abelian additive group.

The exact moduli-space bridge is also inherited, via Aricheta: supersingular
rationality with level1,2,6 gives the audited Monster/Baby/Fi22 prime-support
filtration. Rational endpoints can still have isogeny edges exchanged by
Frobenius. That edge phase is not erased by either the filtration or(7).

### Give the six numerical labels an honest functional graph

For R=|M|/|B| define e_R(b)=max{j:b^j divides R} for every integer b>=2.
This extends prime valuations to composite bases without identifying13=4.
The cited group orders give

    4->2->5->3->7->4,       13->2,       11->1.

At1 the map is not defined; keep it terminal. Replacing R by the actual
involution-class index |M|/(2|B|) changes2->5 into2->4, leaving the2-cycle
4<->2 with5->3->7->4 and13->2 feeding it. This recovers the old central
factor-two sensitivity using a genuine integer function. It does not make
the valuation dynamics a group action or the supplied six labels a six-cycle.
The same bank's sigma-character interpretation in section2 is a separate,
precise map; equality of labels does not identify the two maps.

## 7. Fermat's little theorem and the coefficient3 in Collatz

Interpreting the prompt's3n=1 as the modular inverse condition and the
signed3n±1 maps: for p!=3, FLT gives3^(-1)=3^(p-2) mod p. This solves
3n=±1 locally; it does not select the multiplier3 uniquely. For the actual
odd3n+1 map, every image is nonzero mod3 and an inverse edge

    n=(q*2^a-1)/3

exists exactly when q*2^a=1 mod3. Since2=-1 mod3 this fixes the parity
of a; higher source digits require the full3-adic order2*3^(k-1).
For3n-1, replace the numerator by q*2^a+1, changing the congruence sign.
These are exact local statements, not statements about the number of
global integer cycles or convergence.

The source of the3 in Phi_6 is instead its discriminant-3 and sixth-root
order. These roles share multiplicative-order tools without being the same
dynamical system. For an odd valuation word, keep its entire affine carry:

    word(1,3): (9n+5)/16;      word(3,1): (9n+11)/16.

The same total exponents do not determine the map or its legal source
cylinder. Quasirandom fixed-pattern averages cannot erase that carry.

## 8. Proposed next problems

1. **Growing decorated primes.** Find an explicit rate of growth of the
   square-class rank/conductor for which(1) remains uniform. The proof
   above freezes both; Sylvester anchors are the hostile sparse limit.
2. **Sylvester factor selection.** Classify which selections of prime factors
   break the transitive Jacobi order. The relation(1807/q)=(13/q)(139/q)
   shows exactly what product signs forget; no random-factor theorem follows.
3. **Joint period and response lifts.** For arithmetic maps with a retained
   path/carry response, classify splitting, growth and nilpotent response
   layers together. The19 polynomial cycle and1093/3511 doubling towers
   are distinct positive controls; an equal scalar period is insufficient.

The common programme is to retain exact character relations first, then
bound the remaining coherent modes at a stated scale. Sections1–4 give
proved instances; no universal arithmetic convergence principle is claimed.
