# Paired prime rows, the two Collatz signs, and multiplication

Date: 2026-10-03 (America/Denver).

**Status:** PROVED for the elementary identities and gap statements below;
FINITE-EXACT for the bounded census; OPEN / CONDITIONAL for simultaneous-prime
infinitude and its Hardy–Littlewood prediction. No new theorem ID, priority
claim, or Collatz convergence result. The family-index interpretation of p
remains provisional, as in the preceding note.

This continues [the modular tiling and atom note](tiling_modular_atoms_20261003.md)
with the owner's prime/composite partition, 3n-1, and two-recursion question.

## Inheritance, portfolio, and concept board

Closest mechanisms:

- [THM-4470, pairing ladder](../../01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md)
  already proves the two sum-preserving consecutive pairings for 3n+1/3n-1.
- [THM-4505, small arithmetic graphs](../../01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md)
  and its [full note, Proposition F4](procgen_smallgraph_20260926_small_graphs_arithmetic.md)
  already identify a rectangular grid's fence count a+b+2ab with Sundaram's sieve.
- [THM-4524, selfie loop gauge](../../01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md)
  distinguishes the chosen path, six free bits at order five, and loop switches.
- [HYP-3244, interlocking tiling recursions](../hypotheses/HYP-3244-tiling-half-tiling-interlocking-recursions.md)
  concerns lift versus quotient; it is not a theorem identifying these arithmetic operations.

Canonical hostile: the accelerated 3n-1 map has the cycle 5 -> 7 -> 5.
Corrected near miss: MISTAKE-554 in `01-canon/MISTAKES.md` distinguishes the
codes of two different Collatz conventions from a conjugacy of one map.
Least-used relevant sidecar: the family sign in addition to y=N, followed by
the exact power of two removed from a Collatz step.

| Live object / lane | Operation and retained predicate | Loss / decisive test |
|---|---|---|
| Anchor: paired prime rows | y-projection with A/B tag; prime/composite | Dropping tag identifies different next atoms; N=2 |
| Odd multiplication | a star b=a+b+2ab; exact factorization | Untagged atom loses residue; four product formulas |
| Modular center | inverse of 4, two consecutive pair sums | Central entry is not an iterated orbit; test both signs |
| Niche: prime gaps | sieve by residue classes | Density need not fill holes; factorial construction |
| Wildcard: two recursion laws | interchange and affine composition | Missing cross terms/carry; all-one interchange hostile |
| Prime-preserving paired step | S+(4N-1), S-(4N+1) | Local admissibility is not prime infinitude; four-form sieve |

Cards used: “Search the statement before the method,” “Type every analogy
and every implication,” and “Find the hidden second coordinate in a nearly
true theorem.” The board's new connections are recorded below; no method
card promotion or manufactured tournament is needed.

## 1. Recover the geometry and the two tagged copies

The preceding note proves the mod-10 wedge, its boundary
5,0,5,0,5,4,3,2,1, and 15=6 free positions+4 path positions+5 loop positions.
Six binary tiles cover all 12 order-five tournament isomorphism classes;
they are 64 presentations, not a bijection to classes. Signed coordinate
reflections and transpose reconstruct the multiplication table, with smaller
boundary orbits and the zero axes treated separately.

Write F_k=(k,ceil(k/2),2k+1). The two triples over atom N>=1 are

    B_N=(2N-1,N,4N-1),    A_N=(2N,N,4N+1).

The corresponding odd-modulus central blocks are

    mod 4N-1: [[N,3N-1],[3N-1,N]],
    mod 4N+1: [[3N+1,N],[N,3N+1]].

Define P_A={N:4N+1 prime}, P_B={N:4N-1 prime}; let C_A,C_B be their
respective complements in the positive integers. There are two partitions
of N, not a single partition into P_A and P_B. They may overlap, and they
may both miss an atom. First terms:

    P_A: 1,3,4,7,9,10,13,15,18,22,...
    P_B: 1,2,3,5,6,8,11,12,15,17,...

N=3 lies in both (11,13); N=14 lies in neither (55,57).
The tag is sufficient to reconstruct the odd number from N; without it,
even primality is ambiguous (N=2 gives 7 and 9).

## 2. The exact composite sieve on the two rows

Odd multiplication transports under z(k)=2k+1 to

    a star b = a+b+2ab,    z(a star b)=z(a)z(b).

For positive a,b its image is exactly the indices of odd composites.
The converse follows by factoring any odd composite into odd factors >=3.
This is the Sundaram mechanism already present in THM-4505's fence model.
Adjoining index zero supplies the identity (odd number 1); 1 is neither
prime nor composite and is outside the owner's original list.

For positive r,s, multiplication on the z-labels gives the following
products of tagged families (the product notation in this display means
multiply z-labels and return the unique matching triple):

    A_r A_s = A_(4rs+r+s),    B_r B_s = A_(4rs-r-s),
    A_r B_s = B_(4rs-r+s),    B_r A_s = B_(4rs+r-s).

Consequently the first two polynomial images exhaust C_A, and the mixed
images exhaust C_B. This is an iff, by the factorization argument; the
factor pairs may repeat an output. Family letters multiply by the sign
character mod 4: AA=BB=A and AB=BA=B. Primality is irreducibility in the
full monoid with identity zero, not a property of this two-element quotient.

The owner's p, interpreted as two indices, is

    P(a,b)=(s,t)=(a+b,(2a+1)b),    a star b=a+t.

The prior note proves P is injective, so retaining both outputs loses no
input information. Retaining only their family letters does lose it.
Its outputs need not have the same family: P(1,2)=(3,6) has types (B,A).

For q(I,J)=(O,E)=(2(I+J)-1,4J), 1<=I<J, the output atom addresses are

    y(F_O)=I+J,    y(F_E)=2J,    address gap=J-I.

Thus q supplies an ordered B/A pair at two addresses, not a primality
criterion. For instance q(2,3)=(9,12) has z-labels 19 and 25.

## 3. The central blocks are the two Collatz pairings

Distinguish the shortcut maps from the fully accelerated odd maps:

    T_sign(n) = n/2 if n even, (3n+sign)/2 if n odd;
    S_sign(u) = (3u+sign)/2^v2(3u+sign), for positive odd u.

The modular-center rows are exactly

    (T_+(2N-1),T_+(2N))=(3N-1,N),
    (T_-(2N),T_-(2N+1))=(N,3N+1).

Both outputs sum to their pair's original sum, respectively 4N-1 and
4N+1. This explicitly connects the modular inverse-of-four calculation
to THM-4470; no primality assumption is needed. In particular, 3N-1 in
the B central block is the ordinary + shortcut on 2N-1, whereas 3N+1 in
the A central block is the minus shortcut on 2N+1. The sign cannot be
inferred just from the displayed affine expression without its input.

For the odd labels themselves, division by all powers of two gives

| Input | S_+ | S_- |
|---|---|---|
| 4N-1 (B) | 6N-1 | oddpart(3N-1) |
| 4N+1 (A) | oddpart(3N+1) | 6N+1 |

Here oddpart(m)=m/2^v2(m). The proof is the four factorizations
12N-2=2(6N-1), 12N-4=4(3N-1), 12N+4=4(3N+1), 12N+2=2(6N+1).
The boundary S_-(1)=1 and S_-(5)=7, S_-(7)=5 prevents any conclusion
that the minus map has only the trivial cycle.

The atom projection alone does not support either dynamics: 7 and 9 share
N=2, but S_+(7)=11, S_+(9)=7 have addresses 3 and 2; S_-(7)=5,
S_-(9)=13 have addresses 1 and 3. Even “family+prime” is insufficient:
7 and 23 are both B primes but S_+ sends them to prime 11 and composite
35. Under S_-, B primes 7 and 47 go to prime 5 and composite 35.
The strongest survivor is the full tagged-coordinate identity above.

There is a second exact multiplication connection. For the **unhalved**
Collatz map U_sign(n)=n/2 on evens and 3n+sign on odds, set
phi_sign(n)=2n+sign. Then

    phi_sign(U_sign(n)) = (phi_sign(n)+sign)/2 if n even,
                         3 phi_sign(n)          if n odd.

For sign + the domain is odd integers >=3; for sign - it is odd integers
>=1. These are conjugacies to explicitly different piecewise maps, not a
conjugacy between U_+ and U_-. Also 1 star n=3n+1. Thus odd multiplication
by 3, parity-dependent halving, and an affine change of coordinates are
literal operations in the same dictionary. Their finite identities do not
provide a global descent bound.

## 4. What can be proved about the gaps

**PROVED: both prime rows have arbitrarily long simultaneous empty blocks.**
Given L>=1, take K=(4L+2)! and N_j=K/4+j for 1<=j<=L. Then

    4N_j +/-1 = K+(4j +/-1)

is divisible by 4j +/-1, a proper divisor >=3. So all L consecutive
addresses lie in C_A intersect C_B. Multiples of K give arbitrarily far
such blocks. Neither prime row can uniformly fill the other's holes.
N=14 is the smallest both-composite atom; the factorial proof is a
construction, not a claim of minimality.

Both prime rows are infinite by elementary Euclid arguments. For any
finite list of primes 3 mod 4 with product R, 4R-1 has a new prime divisor
3 mod 4. For any finite list of primes 1 mod 4 with product R,
(2R)^2+1 has a new odd prime divisor p; the order of 2R modulo p is 4,
so p=1 mod 4. Thus the empty blocks also force unbounded gaps between
successive entries in either prime-address list and in their union.

**PROVED: successive composite addresses in each row have gap at most 3,
and this is sharp.** C_A contains 2,5,8,... because 4N+1 is then a
multiple of 3 greater than 3. C_B contains 4,7,10,... for the same reason
applied to 4N-1; N=1 is the prime exception 3. C_A starts at 2, C_B at 4,
so the progressions give the bound without an initial exception in the
consecutive-gap statement. Equality occurs at C_A's 2->5 and C_B's 4->7.
Consequently no prime row has three consecutive addresses beyond the B
initial block 1,2,3. These are gap statements about selected address lists,
not uniform bounds on gaps between primes.

Same-row prime gaps in the original odd numbers are exactly four times
their address gaps. There are two placements of twin primes:

    vertical: (B_N,A_N) -> (4N-1,4N+1),
    diagonal: (A_N,B_(N+1)) -> (4N+1,4N+3).

Their union accounts for all odd twin-prime pairs. Infinitely many
vertical pairs is a restricted twin-prime assertion, not equivalent by
itself to the usual twin-prime conjecture. Infinitely many in at least
one of these two placements is equivalent to that conjecture.

## 5. A precise conjecture joining primes and both Collatz signs

**OPEN conjecture (an instance of Dickson's prime-tuples conjecture):**
there are infinitely many N for which

    4N-1, 4N+1, S_+(4N-1)=6N-1, S_-(4N+1)=6N+1

are all prime. For N>1 these are four distinct numbers. The four forms
have no fixed prime obstruction: N=0 modulo any prime makes them +/-1.
The number nu(p) of forbidden residues is 0 at p=2, 2 at p=3 and p=5,
and 4 for every p>=7. For p>=7 their four roots +/-1/4,+/-1/6 are
distinct; collisions could only divide 2 or 10.

**CONDITIONAL:** the generalized Hardy–Littlewood prediction is

    count(N<=X) ~ S * integral_2^X dt /
      [log(4t-1) log(4t+1) log(6t-1) log(6t+1)],
    S=product_p (1-nu(p)/p)/(1-1/p)^4 > 0.

This is a one-variable prime-tuples problem; finite-complexity theorems
for multiple independent variables do not settle it. The primary source
used for the conjectural framework is Green–Tao,
[Linear Equations in Primes, v2, Conjectures 1.2/1.4 and the discussion of
Dickson on printed page 6](https://arxiv.org/pdf/math/0606088).
No novelty or unconditional prime-infinitude claim is made here.

FINITE-EXACT: for 1<=N<=100000, there are 125 such atoms, starting
1,3,18,45,87,357,528,942,1470,2523. For N=3 the four primes are
11,13,17,19; at N=18 they are 71,73,107,109.

This supplies a target-preserving connection: source = prime atoms with
two prime images; target = prime points of four affine forms; map = the
two specified opposite-sign odd steps. No information is lost with N and
the branch labels retained. It does not concatenate into one Collatz orbit,
since the two arrows belong to different maps. The next frontier is
simultaneous primality, not a missing algebraic identity.

## 6. Two recursion laws and the Eckmann–Hilton test

On nonnegative indices, addition and star both have unit 0. But their
interchange defect is exactly

    (a+b) star (c+d) - [(a star c)+(b star d)] = 2(ad+bc).

At a=b=c=d=1 the two sides before subtraction are 12 and 8. Hence
Eckmann–Hilton cannot identify these operations. On nonnegative indices
the defect vanishes iff ad=bc=0; generic strictly positive inputs never
satisfy it. This is a measured obstruction, not just a mismatch of names.
The standard theorem's requirement is two monoid structures with the
interchange law; see Batanin,
[The Eckmann–Hilton argument and higher operads, introduction](https://arxiv.org/pdf/math/0207281).

An exact positive replacement is the affine representation

    f_a(t)=a+(2a+1)t,    f_a composed with f_b=f_(a star b).

Thus the two orders of odd multiplication commute, while mixed addition
and multiplication retain the cross terms. The previous
[one-out-edge systems note, section 5](collatz_one_out_edge_systems_20260930.md)
discusses the parallel slope-versus-carry obstruction for Collatz words.
For the simple rational letters h(t)=t/2 and g_sign(t)=3t+sign,
g_sign(h(t))-h(g_sign(t))=sign/2: forgetting the translation makes the
slopes commute and loses the actual affine map. Parity legality must also
be retained before interpreting any word as a trajectory.

## Reproduction and finite scope

Run from the repository root:

    python3 04-computation/experiments/paired_prime_collatz_20261003.py

[Source](../../04-computation/experiments/paired_prime_collatz_20261003.py),
[exact output](paired_prime_collatz_20261003.out).
No random sampling, external data, or inherited filters. The universe is
every atom 1..100000; Eratosthenes covers every integer through 600001,
with independent trial division for 2..40001. Exhausting all unordered
odd factor pairs with smaller factor <=sqrt(400001) independently recovers
both composite rows (522157 factor pairs, 200000 labels). A second full
quartet census uses trial division directly for every atom through 100000
and agrees with the sieve result, including all exclusions.

The paired-step formulas and affine conjugacies are checked through 10000;
the factorial proof is instantiated for lengths 1..20; the interchange
defect is checked on {0,...,5}^4; q on all 1275 pairs 1<=I<J<=51.
Positive controls: N=3 prime quartet, mod-9/11 centers from the prior note,
and both twin placements. Hostiles: N=14 joint prime hole, mixed p/q output
types, same-atom non-descent, prime-to-composite steps, interchange at all
ones, and the 3n-1 cycle 5->7->5.

In the finite universe P_A/P_B have sizes 16900/16959, with 1934 shared
addresses. The longest both-composite run is 92566..92592 (length 27).
These finite records do not provide an asymptotic or any orbit convergence
statement. The universal assertions above rest on their displayed proofs.
