# Matched ternary clocks: norm-restricted primes and a closed carrier

2026-10-29 overnight continuation. **PROVED elementary consequences of
the classical Eisenstein norm formula / FINITE-EXACT supporting audit.**
The source here is an integral alternating form with two two-dimensional
blocks. It is not the earlier pair of native four-vertex tournaments.

## E1. An actual clock forces the Eisenstein grid

Let J=[[0,1],[-1,0]] and A=[[0,-1],[1,-1]]. Then A^3=I and A^T J A=J.
Fix nonzero integers d,e and the alternating form

    K(d,e,B)=[ dJ   B ],             T=diag(A,A).
            [-B^T eJ]

The clock T preserves K exactly when A^T B A=B. Solving those four
linear equations over Z gives precisely

    B(a,b)=[[a,b],[-a-b,a]],       a,b in Z.

Thus det B=Q(a,b)=a^2+ab+b^2, the norm on Z[(1+sqrt(-3))/2].
This is a literal intertwiner calculation, not a grid analogy. The two
coordinates and the matched choice of clock generator are retained.
Write z=de and t=z-Q(a,b). Exclude t=0, which is the singular boundary.
The preceding alternating Smith calculation gives

    coker(K) with its pairing = H_s + H_(|t|/s),
    s=gcd(d,e,a,b),       H_m=(Z/m)^2 with standard pairing.       (6)

In particular s^2 divides t. Different d,e with the same product can
have different s, so neither a determinant nor the norm fixes the group.

## E2. A complete all-order count of the attainable paired groups

Define r(N)=0 for N<0, r(0)=1 and, for N>0,

    r(N)=6 sum_(q|N) chi_(-3)(q),
    chi_(-3)(q)=0,1,-1 according as q=0,1,2 mod3.              (7)

This classical norm count is equivalently zero if an inert prime q=2 mod3
has odd exponent in N, and otherwise

    r(N)=6 product_(q=1 mod3) (v_q(N)+1).

The norm is Euclidean and has six units. The prime3 ramifies; primes1 mod3
split into two conjugate factors, with v_q(N)+1 choices; primes2 mod3 remain
prime and must occur to even exponent. Unique factorization gives (7).
These are classical arithmetic inputs, not new discoveries; see the direct
[course derivation, exercises4-17](https://categorified.net/23Winter3032/HW6.pdf).

Now fix positive d,e. For s dividing gcd(d,e), put D=d/s, E=e/s,
g=gcd(D,E). The exact number of labelled clock-preserving cross blocks
whose paired group is H_s+H_(s n), n>=1, is

    sum_(epsilon=+1,-1)
      sum_(h|g, h^2 divides N_epsilon)
        mu(h) r(N_epsilon/h^2),
    N_epsilon=DE-epsilon n.                                  (8)

Terms with negative N are zero. For N=0 every h^2 divides N and r(0)=1.
If s does not divide gcd(d,e), the answer is zero.

Proof: write a=sx,b=sy. Equation(6) requires
Q(x,y)=DE-epsilon n and gcd(D,E,x,y)=1. Mobius inclusion-exclusion
over h dividing g removes precisely the common content. Dividing x,y by h
divides their norm by h^2. The two signs cannot overlap because n>0.
Hence a group is attainable iff (8) is positive; this is a complete
finite divisor algorithm at every n, with no bounded-weight cutoff.
Counts are for labelled integer couplings, not orbits under re-marking.

This recovers the Liouville/Mobius distinction concretely: inclusion-exclusion
must remove common prime divisors, so the squarefree-supported Mobius
function occurs. Parity of the total prime-factor count is not a substitute.

## E3. Exactly which prime planes can be generated from d=e=3

For H_p with p prime, s=1 and (8) reduces to

    count(H_p)=r(9-p)+r(9+p).                                (9)

The primes2,3,5,7 occur. For every prime p>9, the exact criterion is:

    H_p occurs iff every prime q=2 mod3 has even valuation in p+9.

Such p must be7 modulo12: a positive Eisenstein norm is0 or1 modulo3,
and an even norm is divisible by4. This congruence is necessary, not
sufficient; p=31 has p+9=40 and is excluded. The first larger admitted
primes include19,43,67,103,139,163. Formula(9), rather than a finite list,
classifies all primes. It makes no density or infinitude claim about this
shifted norm-prime set.

Examples preserving the named clock:

- a=2,b=0 gives Q=4 and H5;
- a=2,b=4 gives Q=28 and H19;
- H11 would require Q=20, which is impossible.

Dropping the clock admits B=diag(1,-2), with Pfaffian11 and content1.
This is the cheapest explicit test that the new prime exclusion is an
action constraint, rather than a property of all two-port attachments.

For the requested composite examples:

- H25 requires norm34, whose inert primes2 and17 have odd exponents;
- H57 requires norm66, whose inert primes2 and11 have odd exponents;
- H27 has norm36 representations, but all have a,b divisible by3.
  Mobius content subtraction removes all six from the H27 count.
  They instead produce **H3+H9**.

Thus all three cyclic paired planes H25,H27,H57 are excluded in this
matched irreducible-clock model, although arbitrary weighted joins realize
them. The earlier native four-vertex C3 census permits H9/H27 because
its source representation and admissible matrices are different.

## E4. Smallest closed carriers for these actual cross-block operations

Allow cross attachments B(a,b) -> B(a+u,b+v), u,v arbitrary integers;
global sign reversal K -> -K; and block port re-markings R1,R2 in GL2(Z)
commuting with the dual clock A^T. The latter act on the canonical current
chart by B -> R1^-1 B R2^-T. Such R_i are the six Eisenstein units, all
of determinant1, so d,e and z are unchanged. The new a,b are the first
row of the transformed B. Global sign reversal sends (a,b) to(-a,-b)
and d,e to -d,-e, leaving z and gcd(d,e) unchanged. Singular updates
are reported as zero and stop a query branch.

For determinant observations, the exact minimal carrier is

    (z,a,b).                                                (10)

Sufficiency follows from det K=[z-Q(a,b)]^2 and the operation laws.
For necessity query translations along each coordinate axis. Writing
t=z-Q(a,b), L1=2a+b and L2=a+2b, the two response polynomials are

    f_i(w)=[t-L_i w-w^2]^2.

Their cubic coefficients are2L_i. The first quadratic coefficient is
L1^2-2t. Thus the responses recover L1,L2,t and then

    a=(2L1-L2)/3, b=(2L2-L1)/3, z=t+Q(a,b).

Five evaluations w=0,...,4 on each axis, sharing w=0, give nine queries.
Zeros are valid observations, so singular test values do not obstruct
the interpolation. Equal future responses are equivalent to equal (10),
relative to the named clock chart. This is a response-equivalence theorem,
not a claim that three real values require a particular number of bits.

If future observations include the **unmarked paired-group isometry type**,
the exact minimal carrier is instead

    (z,g,a,b),       g=gcd(d,e)>0.                           (11)

It is sufficient by s=gcd(g,a,b) in(6), and its operation laws above.
It is necessary because determinant responses first recover (10); the
legal cancelling attachment (-a,-b) then gives B=0 and the paired group
H_g+H_(|z|/g), from which g is its first paired invariant factor.
This target is nonsingular since d,e are nonzero. Thus the extra content
coordinate is exactly what the determinant query language misses here.

This does not retain named classes of ports in the cokernel or classify
the induced clock action on that finite group. For example, sources
(d,e)=(1,36) and(4,9) at B=0 have the same (11), but the first named
cokernel port has order1 versus4. Thus(11) is not a fully marked Selmer
carrier. The source clock and legal operation labels are preserved; the
declared output group is unmarked. Those are different requirements.

The hostile d,e=(1,9) versus(3,3), both at a=b=0, has identical (10)
and all future determinant responses, but groups H9 versus H3+H3.
The actual group query separates them. This is distinct from the tournament
automorphism81/27 hostile, whose transporter observations remain outside
this language.

Even the norm itself is too small a closed carrier: (a,b)=(1,0) and(0,1)
both have Q=1 and the same initial paired group when d=e=3. The SAME
marked attachment (u,v)=(1,0) changes their norms to4 and3, giving groups
H5 and H6, respectively. The missing coordinate is relative phase in the
Eisenstein lattice, which the cross term of Q reads.

There is also an exact cheapest missing operation. After cancelling B,
allow the single internal direction diag(J,0), which preserves the clock.
Its determinant responses at weights+1,-1 are (z+e)^2,(z-e)^2. Their
difference is4ze, so the already known nonzero z recovers e, and then
d=z/e. This restores the source factors which cross operations could not
read, and explains the relation with full invariant-coupling recovery.

## E5. A full ternary-precision recurrence and compatible CRT charts

Here change the observation language explicitly: observe ONLY determinants
modulo m, with no separate exact-singularity flag. Allow all formal cross
translations, including matrices singular over Z. For a residue state
(z,a,b) in (Z/m)^3, its complete response function is

    R_(z,a,b)(u,v)=[z-Q(a+u,b+v)]^2 mod m.

Every residue state has a nonsingular integer lift if desired, but detecting
exact singularity of a future integer lift would be an additional observable.
Let K(m) be the number of distinct complete response functions. Their
equivalence is closed under translations, sign reversal and the same clock
unit re-markings: these operations reindex the query arguments.

For m=3^k, k>=2, two states have the same complete responses exactly when

    z'=z mod m;
    a'=a mod (m/3), b'=b mod (m/3);
    if z is not1 mod3, also a'-b'=a-b mod m.                 (12)

For z=1 mod3 the last condition is absent. Consequently fibers have size9
in that stratum and size3 elsewhere. The exact all-precision count is

    K(1)=1,  K(3)=7,
    K(3^k)=7*27^(k-1) for k>=1.                           (13)

In particular every new ternary digit after the first multiplies the
carrier count by27. For the fixed source product z=9, the corresponding
phase-carrier count is 3^(2k-1) for k>=1, multiplying by9 per new digit.
This is unbounded precision growth, not a finiteness argument for Sha.

**Proof of(12).** Put L1=2a+b,L2=a+2b and t=z-Q(a,b). Expand

    R=Q(u,v)^2+2Q(u,v)(L1 u+L2 v)
      +(L1 u+L2 v)^2-2t Q(u,v)-2t(L1 u+L2 v)+t^2.

The difference of two response polynomials has degree at most3. Its
third forward differences along u and along u,u,v are respectively
12 times the difference of L1 and4 times the difference of L1+L2.
If all values agree modulo3^k, these give a'=a and b'=b modulo3^(k-1).
The second forward difference along u then gives t'=t modulo3^(k-1),
and hence z'=z modulo3^(k-1).

Write the remaining increments as3^(k-1)(epsilon_z,epsilon_a,epsilon_b).
Their square terms vanish modulo3^k since k>=2. Dividing the response
difference by2*3^(k-1) leaves the function over F3

    (z-r^2)[epsilon_z-r(epsilon_b-epsilon_a)],
    r=a-b+u-v.

Testing r=0,+1,-1 forces epsilon_z=0. It forces epsilon_a=epsilon_b
unless z=1 mod3, where both phase increments are free. This proves
necessity and sufficiency of(12). There are3^(3k-3) response classes
with fiber9 and six times that many with fiber3, proving(13) for k>=2.

At k=1, Q(a+u,b+v)=(a-b+u-v)^2. Source z=0 gives three translated
response functions, z=1 gives three, and z=2 gives one constant function.
This yields seven. The fixed z=9 claim follows directly at k=1 and
from the size-three fibers at higher precision.

For any modulus n coprime to6, every residue triple is distinguished:

    K(n)=n^3.                                              (14)

Indeed the values on the5-by5 grid determine the polynomial coefficients
because the Vandermonde denominators divide powers of2 and3. The cubic
and quadratic recovery in E4 then divides only by2 and3, which are units.

For coprime m,n, CRT is compatible with both states AND query arguments,
so the full response carrier is the product and

    K(mn)=K(m)K(n).                                        (15)

At modulus2, the parity response has a unique minority value whose
location recovers(a,b); that value also recovers z. Hence K(2)=8.
Squares modulo4 depend only on parity, giving K(4)=8 as well, with
eight states in each fiber. It follows, for example, that

    K(76)=8*19^3=54872,
    K(228)=8*7*19^3=384104.

These are exact complete-response counts for this arithmetic attachment
language, not the earlier125 tournament chart choices. At moduli8 and16,
the new finite census gives40 and320 classes; an all-order dyadic formula
is deliberately left OPEN here.

The5-by5 grid is sufficient for comparing complete functions over any
modulus, even when it does not recover polynomial coefficients: fifth
forward differences in either argument vanish, so the initial values
propagate by an integer recurrence without division. This distinction
is why the modular census is exact rather than just a probe sample.

## Prior niche work, scope and next step

The old [dyadic conductor audit](../dyadic_conductor_blue_fibres_20261021/review/evidence/INDEPENDENT_AUDIT.md)
has a sharp norm57 witness whose reduced class form is(25,20,52).
Its conductor40 condition and quadratic character at5 separate a norm
kernel from the square subgroup. Here57 itself is indeed a norm, for
example Q(1,7)=57, but H57 requires the DIFFERENT norm9+57=66. This
explicit shift prevents confusing the two arithmetic predicates. The
common lesson is to retain source, conductor or clock, and the actual map.

The portfolio is marked paired joins / conductor-square hostiles /
arithmetic local characters. This result fills the cross-block restriction
left by [FINITE_CLOCK_PORTS](FINITE_CLOCK_PORTS.md): allowing every invariant
coupling recovers G, but allowing only cross couplings leaves a genuine
source-factor ambiguity and a smaller closed carrier. Declaring the
operation set is therefore mathematically consequential.

The group order is still a square because the form is alternating. These
formal groups are not automatically Sha groups, and a norm-prime criterion
does not bound the prime support of Sha for a fixed arithmetic source.
An actual compatible descent/isogeny map is still required. No Sha-finiteness,
native tournament complement, moonshine or Collatz claim follows.

Next: prove or refute the dyadic precision pattern and supply a cross-source
arithmetic map. The internal-block probe above already reads d,e; named
cokernel classes and their induced clock action must be audited separately
before claiming a fully marked arithmetic carrier.
For native tournaments, repeat the legality test before transporting the
orbit-sum or Eisenstein coupling formulas.

## Reproduction

`python -B eisenstein_clock_joins.py` and the same with `-O` compare(7)
with a full norm-shell enumeration, (8) with direct labelled couplings,
the nine-query decoder with direct determinants, all clock-unit re-markings,
and the phase/content/clock hostiles. The JSON records exact finite bounds.
All-order claims use the displayed proofs; numerical prime scans are examples.
`clock_join_precision.py` separately enumerates the declared finite residue
universes, verifies(12) through modulus81, compares full residue-query banks
with the finite grid at small moduli, and checks the named CRT count76.
