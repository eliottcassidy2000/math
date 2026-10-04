# The sixth clock: Collatz rays, a golden sextic, and recursive prime depth

2026-10-04. **PROVED** for the identities, indexed ray families, number-field
statements, and all-height clock towers proved below. **FINITE-EXACT** for
the explicitly bounded controls. **CITED** for the global primitive-divisor
theorem. **OPEN** for constructing a terminating certificate for every
positive integer or assigning every negative integer to a certified cycle.
The proposed richer carrier is a research construction, not that missing proof.
No historical novelty claim is made.

The strongest connection is an actual object joining the owner's expressions:
if `lambda=(5*phi^4)^(-1/6)` and `t=1/lambda`, then the sextic field Q(t)
reduces modulo2 to F64. Its unit `t+1` has residue order63 and generates an
all-height binary tower with32 child cycles per parent. This gives precise
content to one proposed recursion; it does not identify its dynamics with
Collatz. The Collatz three-ray pattern and the cyclotomic19/76 connection
are separate exact maps, with their missing coordinates retained below.

## Inheritance and working board

The closest mechanisms are [the modulo-six braid and triadic odometer](arithmetic_braids_20260917_collatz.md),
[the audited 63/primitive-divisor comparison](collatz_mod6_20260917_zsigmondy_triad.md),
and [the first-hit/return-rank distinction](arithmetic_seams_20260921_primes.md).
The inherited [divisor-box classification](arithmetic_braids_20260917_divisors.md)
and [k-free extension](collatz_mod6_20260917_divisor_balance_family.md) already
prove the p^2*q*r and p^k*q*r families. They are not discoveries of this session.

Recent input includes [difference families and the 4*19^k phase tower](difference_families_20261004.md),
[rational-anchor integer certificates](rational_anchor_returns_20261004.md),
and the concurrently integrated [denominator-tree child selectors](denominator_fibre_recursion_20261004.md).
The golden phase theorem identifies periodic golden points with primitive
pairs modulo a denominator; passing its integer-realization filter remains
a separate requirement. The canonical hostile is the rational phase with
no integer root, already present above the denominator2->4 lift. Another
hostile here is `H(3)=5`: a residue row is not entirely a power-of-two ray.

The corrected near miss is a map convention: an early commentary in this
session assigned the owner's ordered exponents5,3,7 to the minus sheet.
The inherited convention is `H(n)=(3n+1)/2`, and this gives exactly the
requested source order3,5,1 modulo6. The least-used coordinates are
prime-power depth, the sequence-relative novelty index, and the integral
basis of the sextic field at the bad index prime3.

| Lane | Live object | Invariant or cheapest decisive test |
|---|---|---|
| Anchor | Direct Collatz return rays | Exact integer, valuation, source row, first root |
| Anchor | Recursive signed certificates | Guarded inverse steps and retained terminal certificate |
| Niche | Mersenne and golden clocks | Full prime ideal/valuation and exact return order |
| Niche | A marked divisor box | Exponent profile; direct F,S,U count |
| Wildcard | The lambda sextic | Minimal polynomial, maximal order, residue field |
| Boundary | Recursive phase lifts | Parent map, anchor, child selector, integer-realization gate |

## 1. The exception at six retains depth that prime support forgets

Start `x_0=0`, `x_(n+1)=2*x_n+1`. Then `x_n=2^n-1`. Thus

    x_2=3, x_3=7, x_6=63=3^2*7.

Every prime dividing x_6 has appeared earlier. The second power of3 is new
depth, however:

    ord_3(2)=2,                 ord_9(2)=6.                  (1)

There is no prime p with `ord_p(2)=6`: it would divide63 and hence be3 or7,
whose orders are2 and3. This is an elementary proof for the sixth term.
The classical Bang–Zsigmondy theorem says that coprime positive a>b have
a primitive prime divisor of `a^n-b^n` for n>6. Together with the first
six terms it gives exactly n=1 and6 as exceptions for this sequence.
This global step is **CITED**, from the introduction of the retrieved
[Voutier–Yabuta research paper](https://www.impan.pl/shop/en/publication/transaction/download/product/82701),
which recalls Zsigmondy's theorem; it is not proved by our finite scan.

Put `phi=(1+sqrt(5))/2`, `O=Z[phi]`. A dual event occurs in O:

    phi^3-1=2*phi,              phi^6-1=4*phi^3.             (2)

Since phi is a unit and X^2-X-1 is irreducible modulo2, (2) says that the
only prime ideal appearing at time6 is the prime (2), which already
appeared at time3. Its depth increases. In particular

    ord_(O/(2))(phi)=3,        ord_(O/(4))(phi)=6.           (3)

| Base | Old prime support | Old clock | Increased depth | New clock |
|---|---|---:|---|---:|
| 2 in Z | 3 | 2 | 3 -> 9 | 6 |
| phi in O | (2) | 3 | (2) -> (4) | 6 |

The two factors in6=2*3 exchange their roles. The shared predicate is
**increased return period without new prime support**, proved by (1)–(3).
There is no asserted conjugacy between the two rings. The inherited
prime-power arguments give the all-height clocks

    ord_(3^a)(2)=2*3^(a-1),
    ord_(O/(2^a))(phi)=3*2^(a-1),             a>=1.         (4)

The [prime-clock note](golden_prime_clocks_20261003.md) proves the second,
including the initially delicate binary lift. Taking radicals of ideals
would erase exactly the depth responsible for these recursions.

## 2. The three Collatz rays, including the offset

Use ordinary signed Collatz C, and on odd n use the single-halving expression

    H(n)=(3n+1)/2.

It is different from the fully accelerated odd map
`U(n)=(3n+1)/2^v2(3n+1)`. The owner's ordered powers correspond to H:

| Source row | H(n) | n, j>=0 | First three sources |
|---|---|---|---|
| 3 mod6 | 2^(5+6j) | (2^(6+6j)-1)/3 | 21,1365,87381 |
| 5 mod6 | 2^(3+6j) | (2^(4+6j)-1)/3 | 5,341,21845 |
| 1 mod6, root removed | 2^(7+6j) | (2^(8+6j)-1)/3 | 85,5461,349525 |

Each source reaches1 after one odd operation and its complete halving
chain. These are all the nonroot positive odd sources with U(n)=1.
They are sparse subsets of their residue rows: H(3)=5 is the minimal
counterexample to interpreting the table as every odd integer.

Before removing the root, the third exponent is1+6j, with source1,85,...
Removing `1 -> 2^1` shifts that row's initial exponent from1 to7. This
explains the single offset. It does not require a different Collatz rule.

All three rays lie in the single indexed family

    a_m=(4^m-1)/3, m>=1,
    a_1,a_2,a_3,a_4,... = 1,5,21,85,...,
    R(n)=4n+1,                 R^3(n)=64n+21.              (5)

The three rows are the classes of m modulo3. Lifting the exponent at3 gives

    v_3(a_m)=v_3(m).                                       (6)

Consequently the repeated3 in63 survives division by3 as the first
source21 divisible by3. This connects the clock to a precise Collatz
valuation, rather than just a matching occurrence of six.

For the all-depth offset, instead index by j=m-1>=0. Partition j into
its `3^s` residue classes. Removing j=0 changes the least index in exactly
one class, from0 to3^s; the other classes keep their least indices. This
is an exact recursive marked-branch structure. The inherited stronger
statement is that R is a single cycle on odd residues modulo `2*3^s`,
with

    v_3(R^t(n)-n)=v_3(t),        n odd, t>=1.              (7)

Indeed `R^t(n)-n=(4^t-1)(3n+1)/3`, and 3 does not divide3n+1.
Thus R has a genuine triadic odometer interpretation, already proved in
the braid note. Every such R-ray stays over the same accelerated target:
`U(R(n))=U(n)`.

**Important hostile: novelty depends on the retained sequence.** Although
63 has no new prime in all M_n=2^n-1, its factor7 *is new* in the subsequence
M_2,M_4,M_6=3,15,63 used by these rays. Likewise21 introduces both3 and7
relative to the earlier root sources1,5. A map from the full Mersenne
sequence to the root ray does not preserve primitive-prime novelty.
Retaining the full clock/first-hit index repairs the information loss;
calling all three occurrences the same novelty failure does not.

## 3. The corrected 63 identity leads exactly to19 and76

The arithmetic correction is `63=3*19+6`; `2*19+6=44`. It has a useful
cyclotomic decomposition. With `Phi_18(X)=X^6-X^3+1`,

    Phi_18(2)=64-8+1=57=3*19,
    63=Phi_18(2)+(2^3-2)=3*19+6.                           (8)

The failed new-prime stage at Phi_6(2)=3 is followed along the ternary
index lift6->18 by the new prime19: `ord_19(2)=18`. In the golden ring,

    Phi_9(phi)=7+10*phi,        N(Phi_9(phi))=19,
    Phi_18(phi)=5+6*phi,        N(Phi_18(phi))=19,
    Phi_9(phi)*Phi_18(phi)=19*phi^6.                        (9)

The norm here is `N(a+b*phi)=a^2+a*b-b^2`. Modulo19, the two golden roots
are5 and15. The first factor in (9) vanishes at5 and the second at15;
they select the two distinct primes above19. Factoring phi^18-1 using
its cyclotomic divisors gives

    phi^18-1
      =4*phi^3*Phi_9(phi)*Phi_18(phi)
      =76*phi^9.                                         (10)

This recovers the earlier `G^18-I=76G^9` for multiplication by phi, with
`G=[[0,1],[1,1]]`. The rational anchor -19/11 already has the scaled
3n-11 cycle19->23->29->19, whose closing edge is `3*29-11=76=4*19`.
Equation (10) now supplies an explicit cyclotomic connection to that
denominator. It does not turn a general golden phase into an integer cycle.

There is also an all-height synchronization. Direct calculation gives

    5=2^16 mod19,              15=2^11 mod19,
    v_19(2^18-1)=1.

For all k>=1, the elementary odd-prime lift therefore proves

    ord_(19^k)(2)=18*19^(k-1).                             (11)

This is exactly the common primitive phase period at denominator4*19^k
in the earlier difference-family theorem. After choosing an anchor v on
one such golden cycle, the assignment

    2^j mod19^k  ->  G^j*v mod(4*19^k)                    (12)

is a well-defined bijection of these two cycles intertwining multiplication
by2 with multiplication byG. Both cycles have the exact length in (11).
It preserves the clock and cyclic phase. It is not an additive or ring
map, and it does not preserve an integer Collatz itinerary. The chosen
cycle and anchor are necessary data.

## 4. The annotated lambda is cubic over the golden field

Keep the current expression distinct from the different earlier Higgs
normalization, as [the previous note's lambda section](difference_families_20261004.md#the-two-lambda-formulas-must-remain-distinct)
already does. For the current one,

    alpha=phi^4-1=1+3*phi=sqrt(5)*phi^2,
    alpha^2=5*phi^4,
    lambda=alpha^(-1/3) approximately0.554855513344434.     (13)

The positive cube root is meant. Set t=1/lambda. Since
`alpha^2-5*alpha-5=0`,

    t^6-5*t^3-5=0,
    5*lambda^6+5*lambda^3-1=0.                            (14)

These are degree6 minimal polynomials over Q, up to monic normalization
for lambda. To prove minimality, alpha has norm -5 in K=Q(sqrt5). It
cannot be a cube in K, because a cube would have a rational-cube norm.
Thus X^3-alpha is irreducible over K. Conversely
`phi=(t^3-1)/3` lies in Q(t), so `[Q(t):Q]=3*2=6`.

The earlier correct degree12 identity was not minimal:

    25*lambda^12-35*lambda^6+1
      =(5*lambda^6+5*lambda^3-1)
       *(5*lambda^6-5*lambda^3-1).                       (15)

Here six is a field degree decomposing as a cubic extension over a
quadratic field. It is not being equated with a dynamical period.
No Standard Model prediction or Higgs-mass derivation follows from (13).

### Integral basis and its necessary prime3 coordinate

For L=Q(t), the ring of integers is

    O_L=Z[phi,t],
    basis = {1,t,t^2,phi,phi*t,phi*t^2},
    disc(L)=3^6*5^5=2278125.                              (16)

Here is a proof of maximality, not just a discriminant guess. The relative
polynomial X^3-alpha has discriminant `-27*alpha^2`. Its only possible
bad prime ideals in O_K lie over3 or5. At the prime over5, alpha has
valuation1, so it is Eisenstein. At3, which is inert in K, put Y=t-1:

    Y^3+3*Y^2+3*Y-3*phi=0,

again Eisenstein. In each local extension the generating uniformizer has
degree3 and the residue field is unchanged. The valuations of the three
terms in a basis expansion are distinct modulo3, so an integral element
has integral coefficients. This proves local maximality at both possible
bad primes; elsewhere the order discriminant is a unit. Taking the
relative discriminant's norm and multiplying by disc(K)^3 yields (16).

By contrast the power order Z[t] has discriminant `3^12*5^5`, and
`[O_L:Z[t]]=27`. This matters: reducing its polynomial modulo3 gives
`(X-1)^6`, but it does **not** make the true residue field F3. The integral
coordinate phi survives, and the residue field is F9. The expression
`phi=(t^3-1)/3` cannot be reduced by division by3 at that prime.

| Rational prime | Number of primes above it | Ramification e | Residue degree f |
|---:|---:|---:|---:|
|2|1|1|6|
|3|1|3|2|
|5|1|6|1|

The last two rows follow from the Eisenstein proof; the first follows next.
Exact trace-matrix determinants and an independent PARI `nfinit`/`nfdisc`
calculation agree with (16) and the index27. The computation uses the
[documented PARI number-field routines](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html).

## 5. An exact bridge from lambda to the 63 nonzero elements of F64

Modulo2, the polynomial for t becomes

    X^6+X^3+1=Phi_9(X).                                   (17)

This is irreducible over F2. One proof starts in F4=F2(phi): there
`alpha=phi^2`, a nonidentity element of order3. It is not a cube in F4,
so X^3-alpha is irreducible there. A root t satisfies
`phi=t^3+1`, recovering the quadratic subfield. Hence its degree over F2
is6. Since the index27 is odd, the power order and maximal order agree
at2, proving

    O_L/(2) = F64,               |F64^*|=63=3^2*7.         (18)

Thus the annotated sixth root and the sixth term of2x+1 are connected by
an explicit number field and reduction map. They are not merely assigned
the same count. In this field write `b=t+1`. Direct polynomial reductions give

    ord(t)=9,
    b^9=t^5+t^2+t !=1,
    b^21=t^3+1=phi !=1,
    b^63=1,
    b^56=t.                                              (19)

The prime divisors of63 are3 and7, so (19) proves that b has exact order63.
Multiplication by b cycles through every nonzero element. The degree-six
field has subfields F4 and F8, and `F64^*` is cyclic of order9*7. The
golden subfield contributes its order3 subgroup, while t supplies order9.
Also b is a unit in O_L itself: its norm is f(-1)=1 for
`f(X)=X^6-5*X^3-5`.

An equivalent six-bit implementation uses b as the basis generator. Its
minimal polynomial is

    X^6-6X^5+15X^4-25X^3+30X^2-21X+1,

which reduces to the primitive binary polynomial `X^6+X^4+X^3+X+1`.
The recurrence

    s_(j+6)=s_(j+4)+s_(j+3)+s_(j+1)+s_j mod2

therefore visits all63 nonzero six-bit states before returning. The zero
state is fixed. This is a standard linear feedback shift-register
construction, here obtained from the specified lambda rather than an
arbitrarily chosen polynomial. The finite-state count is exact; its
periodicity is not a Collatz convergence claim.

This gives two different maps that must not be conflated. The scalar
sequence2^n-1 counts the nonzero elements of an n-dimensional binary
vector space. The element b generates a cycle on the nonzero elements
when n=6 in the particular field constructed here. It is not the map
`z -> 2z+1` inside that field: characteristic2 would make that map constant1.

There is also a useful loss test: reducing a polynomial *only in phi*
modulo2 stays in the four-element subfield. To access all64 elements,
one must retain the cubic coordinate t; the existing golden digit
register cannot acquire this extra information merely by being renamed.

## 6. A proved recursive tower with32 children per cycle

This construction has an all-height extension, rather than stopping at a
single field of64 elements. For a>=1 put

    A_a=O_L/(2^a),
    V_a={v in A_a: v mod2 !=0},
    B_a(v)=b*v, b=t+1.                                   (20)

Each member of V_a is a unit. Its exact period under B_a is

    L_a=63*2^(a-1),
    |V_a|=63*2^(6(a-1)),
    number of cycles=2^(5(a-1)).                          (21)

**Proof of the period.** The base order is (19). Work in the power basis,
which is valid locally at2. Direct exact reductions give

    b^63 =1+2*t^4                         mod4,
    b^126=1+4*(t^2+t^4+t^5)               mod8.           (22)

The coefficients after the first nonzero power of2 are units, because
their reductions are nonzero in F64. Starting with the second identity,
squaring `1+2^a*u`, a>=2, increases its 2-adic valuation by exactly1:
the new coefficient `u+2^(a-1)*u^2` remains a unit. This proves the order
formula in (21). If `b^j*v=v`, canceling the unit v makes this exactly the
order condition for b. The phase and cycle counts follow.

**Parent/child theorem.** Reduction `A_(a+1)->A_a` maps every child cycle
onto a parent cycle, with degree2. There are exactly32 children per parent.
Every parent phase has64 lifts, while each child cycle visits that fibre
twice. More explicitly choose an anchor v, a fixed lift v_tilde, and write
a lift as `v_tilde+2^a*w`, with w in F64. A parent return acts by

    w -> w+c,
    c=((b^L_a-1)/2^a)*v_tilde mod2 !=0.                   (23)

The child labels are the32 cosets of the one-dimensional F2-line {0,c}
in F2^6. To compute a label, choose a coefficient r where c_r=1 and
replace w by `w+w_r*c`; its rth coefficient becomes0. The remaining five
bits label the child. The anchor and chosen lift are part of this label.

This is the same *proved local mechanism* as the recent golden binary
and ternary trees, with a different dimension. Generally a nonzero return
translation on F_p^d has p^(d-1) orbits of length p. When a lift's return
has that form, it creates that many children. A zero return or nontrivial
linear part requires a different analysis; the formula is not automatic.

| Tower | Lift dimension d | Prime p | Child cycles | Period multiplier |
|---|---:|---:|---:|---:|
| R(n)=4n+1, triadic odometer |1|3|1|3|
| Golden denominator powers of2 |2|2|2|2|
| Golden denominator powers of3 |2|3|3|3|
| Golden denominator4*19^k |2|19|19|19|
| This sextic b=t+1 tower |6|2|32|2|

The dimension and return translation explain the different branch counts.
Matching the word “recursive” without those data would lose the mechanism.

### A concrete use: exact finite-ring Fourier separation

During integration, incoming commit79e9966a60 supplied
[all-period chart separation](periodic_chart_separation_20261004.md): two
distinct marked primitive growth charts repeated at least twice have
disjoint first-descent cylinders. A primitive complex character annihilates
both repeated words; if a common source forced them to differ only in their
last entry, that remaining coefficient could not vanish. Its
[coverage synthesis](subdivision_carriers_and_anchor_coverage_20261004.md#5-first-descent-clocks-turn-a-finite-atlas-into-a-general-theorem)
also proves that increasing this atlas's period cutoff cannot cover7.
That structural miss is another reason to retain a separate route grammar.
The incoming commit also repairs a denominator parenthesis in the older
rational-anchor proof; the repaired `d*2^t*n-z` is retained.

Our new register implements the Fourier mechanism over finite rings for
all lengths T>1 dividing63, with a necessary precision bound. For k>=1 set

    r_k=(2^(k-1))^(-1) mod63,       e_k=2^(k-1)*r_k,

where the inverse is taken in1,...,62 **before** multiplication (there is
no final reduction of the product). Thus e_k is divisible by2^(k-1) and
is1 modulo63. Define `tau_k=b^e_k mod2^k`. These elements have exact order63
and reduce compatibly: e_(k+1)-e_k is divisible by63*2^(k-1). Put

    zeta_(T,k)=tau_k^(63/T),
    E_(T,k)(a)=sum_(j=0,...,T-1) a_j*zeta_(T,k)^j.         (24a)

If a is a word of period p properly dividing T, then E(a)=0. Its geometric
sum has numerator `zeta^T-1=0`, and denominator `zeta^p-1` is a unit:
its reduction in F64 is nonzero. If two such repeated integer words agree
except for a last-entry difference D, then

    0=E(a)-E(a')=D*zeta^(T-1)
    implies D=0 mod2^k.                                 (24b)

At fixed precision this only detects D not divisible by2^k. With the
retained bound `|D|<2^k`, it proves D=0; doing this at every k has the same
conclusion. The smallest hostile is D=2, invisible modulo2 and visible
modulo4. Thus prime-power depth is operational information, not decoration.
The field F64 alone cannot reproduce the unrestricted integer conclusion.

For example T=21 supports both the inherited rational anchor word112
repeated7 times and the -17 anchor word1112114 repeated3 times. Both
annihilate (24a). Under the hypothetical shared-source condition, the
last-entry argument separates them exactly as in the incoming theorem.
The test supplies an exact finite-ring representation of that proof for
these lengths, not new coverage beyond the existing all-period theorem.
Even periods, or odd periods not dividing63, are outside this root-of-unity
register; equal-looking counts do not remove that restriction.

## 7. What survives of the p^2*q*r analogy

For a positive integer N>=2, let F be its number of proper nontrivial
divisors, S the number of those divisors that are squarefree, and U the
number that are prime. For distinct primes p,q,r and N=p^2*q*r,

    F=10,                 S=7,                 U=3,
    F=S+U.                                               (24)

The inherited classification is `N=p`, `p^3`, or `p^2*q*r`. The lattice
for the last family is the exponent box C3*C2*C2: a Boolean cube of8
vertices plus a four-vertex face. Removing that extra face's top vertex
leaves3 proper nonsquarefree divisors. This proves its balance mechanism.

The offset analogy can be made explicit as a **chosen combinatorial map**:
label the three branches by p,q,r, place p on the delayed branch, and add
the delay vector(1,0,0) to baseline exponents(1,1,1). The result is(2,1,1).
This preserves a distinguished coordinate. The prime labels are extra
choices; the map does not preserve the numerical Collatz operation.

There is also a literal completion using the new field size:

    63=3^2*7               gives (F,S,U)=(4,3,2),
    2*63=126=3^2*2*7       gives (F,S,U)=(10,7,3).         (25)

Thus63 itself fails the balance by1; adjoining the missing prime2 repairs
it. This addition of a prime axis, the cube-root coordinate giving F64,
and the delayed Collatz branch are three specified operations, not one
operation inferred from a matching shape.

The recursive boundary is informative. An exponent profile with one2
and r-1 ones has defect

    F-S-U=2^(r-1)-r-1.                                   (26)

It vanishes for r=3, but at the next triadic partition r=9 it is246.
Therefore the marked-branch construction recurses while the fixed
divisor equality does not. The useful inherited replacement is to retain
three prime axes, deepen p to p^k, and replace squarefree by k-free:

    N=p^k*q*r, k>=2:       F=4k+2, S_k=4k-1, U=3.         (27)

This is the correct all-depth divisor family. A fourth connection, also
inherited, is the squarefree density of the odd residue rows: the row3
mod6 has density6/pi^2 relative to that row, while rows1 and5 each have
9/pi^2. Of the three odd lifts of row3 modulo18, one is divisible by9;
rows1 and5 avoid3 already. The [squarefree braid note](arithmetic_braids2_20260917_squarefree_symmetry.md)
proves these constants. This marked exclusion is a prime-square local
condition; it is not the deletion of the Collatz root.

## 8. How to use these structures in an integer certificate carrier

An integer alone forgets a representation's clocks, valuation depth,
branch labels, and arrival witness. A useful object should retain these
as distinct coordinates. The new sextic supplies a particularly explicit
recursive arithmetic register:

    z=a0+a1*t+a2*t^2+a3*phi+a4*phi*t+a5*phi*t^2,
    a_i in Z, together with compatible z mod2^k.           (28)

Its additions and products are exact; the integral basis prevents the
prime3 loss. A **proposed** enlarged arithmetic fibre over an integer n is
the n-fibre in `Z x O_L`, projected by `(n,z)->n`, together with a finite
ordered route certificate when one has been constructed. The product ring
has a genuine ring projection to Z; a projection from O_L alone to Z
preserving1 would be impossible, because (14) has no integer root.

For a nontrivial coupling to the orbit, the register can store parity
history. Define the halved signed map `T(n)=n/2` for even n and
`T(n)=(3n+1)/2` for odd n. Then the lift

    (n,z) -> (T(n), b*z+epsilon(n)),
    epsilon(n)=n mod2 in {0,1}, b=t+1>2,                  (29)

projects to T. Starting at z=0 and retaining the word length, z stores the
finite parity word injectively: the largest differing power of b exceeds
the sum of all lower powers, since b>2. The six integer coefficients in
(28) store that exact history in a fixed algebraic basis. Reduction modulo
2^k gives compatible finite registers. The zero-input clock is precisely
the tower of section6; a driven parity word changes it, so one cannot
apply the zero-input period formula to the full lifted Collatz map.

This is more informative than a bare integer and is a meaningful history
encoding. It is not yet a termination certificate. Neither its unbounded
coefficients nor its finite phase cycles supply a decreasing rank.

A completed certificate must instead retain guarded inverse macros and
a terminal witness. For an odd certified target u with3 not dividing u,
choose a0=1 when u=2 mod3 and a0=2 when u=1 mod3. Then

    n_j=(2^a0*4^j*u-1)/3, j>=0,          U(n_j)=u.        (30)

Equation (30) stores an entire inverse ray using u,a0,j; attaching the
certificate for u gives an actual finite certificate for each n_j. A
dyadic sleeve2^e*n_j can be attached as well. The same rule works on the
negative sheet with a certified target in the cycles through-1,-5,-17.
This is the inherited route grammar, now with explicit room for the
triadic index and the sextic register. The root1 itself requires the
usual first-arrival convention; its self-return is not a new first hit.

The next decisive task is a rule that uses source-decodable data in these
registers to choose finitely many guarded macros ending at a certified
smaller target. A finite memory that records only arbitrary prefixes
does not supply this rule. The earlier rational-anchor compiler does
give all-height additional cylinders and remains a sound component.
Coverage of every integer remains **OPEN**, as does completeness of the
known negative cycles. No assignment to an unknown “respective cycle” is
silently assumed.

## Connection audit and reproduction

| Source -> target | Map/predicate retained | Lost data; required addition |
|---|---|---|
| Full Mersenne sequence -> Collatz root ray | M_(2m)/3=a_m, exact valuation | Prime novelty changes under subsequence; keep original index/clock |
| Base2 clock -> golden19 clock | Anchored bijection (12), exact period | Ring operations and itinerary; keep phase anchor and integer gate |
| Current lambda -> sextic -> F64 | (13),(14), reduction modulo2 | Real size and ramification at3; keep maximal integral basis |
| Sextic2-adic level -> next level | (23), 32 child labels, doubled period | Rotation and lift origin; keep chosen anchor |
| Three marked branches -> p^2*q*r | Chosen labels and one delayed coordinate | Canonical prime choice and dynamics; keep this typed as a model |
| Integer -> richer history/certificate | Projection (29), checked macros (30) | A total certificate constructor; still OPEN |

The concept board changed after each positive/hostile result: novelty
required an index sidecar; the cyclotomic19 gave an exact synchronized
clock; the lambda simplification produced a field; its index27 forced an
integral basis; the finite-field generator yielded an all-height tower;
the divisor hostile required k-free predicates; the driven-register test
separated memory from termination. These are the stopping boundaries for
this pass, rather than untested claims of one universal recursion.

Run from the repository root:

```text
python3 04-computation/experiments/sixth_clock_branches_20261004.py --pari
python3 -O 04-computation/experiments/sixth_clock_branches_20261004.py --pari
```

The [script](../../04-computation/experiments/sixth_clock_branches_20261004.py)
and [output](sixth_clock_branches_20261004.out) use exact integers and
rationals. There are no probabilistic primality tests. Universes: Mersenne
indices1..60;30 members of each displayed ray;100 root valuations; complete
triadic residue orbits through3^7; dual clock depths1..12; synchronized19
depths1..8; all4320 primitive phases at q76; all63 nonzero F64 elements;
all4032 primitive sextic phases modulo4; sextic order depths1..16; child
selectors at depths1..5; k-free divisor boxes for k2..8. The q4 sextic
census has32 cycles of length126. Norm identities and both trace
discriminants are computed exactly. PARI independently verifies the
minimal polynomial, maximal order, index, and prime decompositions.
Compatible odd roots and finite-ring Fourier tests cover precisions1..8,
all T in{3,7,9,21,63}, every proper divisor period, and last-entry changes
1,2,3,2^k, including the explicit two period21 anchor words.
The proposed history encoding is checked on all4096 binary words of
length12 and on12 steps from each signed source-200,...,200. Inverse macros
are replayed for five signed targets and26 ray indices each.
No independent-agent or formal proof-assistant audit is claimed.
Ordinary and optimized Python runs produced byte-identical output; linked
local targets and whitespace checks passed. The initial repository-wide
document check reported a pre-existing hypotheses-index budget violation;
incoming commit79e9966a60 repaired it. The integrated document check passes,
including the startup byte budget.
