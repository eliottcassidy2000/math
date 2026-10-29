# Ramanujan's surds: exact traces, changed CM fields, and the missing conjugates

2026-10-28. **CITED** for the displayed singular-value evaluations in the
user's quotation; **PROVED** for the elementary trace deductions and the
g18 obstruction; **FINITE-EXACT** for reduced forms; **NUMERICAL-VERIFIED**
for decimal values and independent product comparisons.

The quotation is Ramanujan, *Modular equations and approximations to pi*,
Quarterly Journal of Mathematics45 (1914),350--372, sections4--5. The
[primary transcription](https://ramanujan.sirinudi.org/Volumes/published/ram06.html)
and [PDF](https://ramanujan.sirinudi.org/Volumes/published/ram06.pdf) recover
the corrupted subscripts and powers: G_n^24 and g_n^24 are powers of class
invariants, not invariants indexed by24n. The restored definitions are

    q=exp(-pi sqrt(n)),
    G_n=2^(-1/4) q^(-1/24) product_(j>=0)(1+q^(2j+1)),
    g_n=2^(-1/4) q^(-1/24) product_(j>=0)(1-q^(2j+1)).

The supplied special values give

    g22²=1+sqrt2,       G37⁴=6+sqrt37,
    g130²=(2+sqrt5)(3+sqrt13)/2.

## M1. The near-integers are Pell--Lucas traces

For a root u of X²-tX+c, define A_r=u^r+(c/u)^r. Then
A_0=2,A_1=t,A_(r+1)=tA_r-cA_(r-1), so integer t,c give integer A_r.
This is the same trace mechanism as Lucas numbers, with the appropriate
trace and norm parameters. At c=-1, even powers pair u^r with u^(-r).

For u=1+sqrt2, t=2,c=-1, the recurrence gives A_12=39202. For
u=6+sqrt37, t=12,c=-1, it gives A_6=3111698. Thus the exact integer targets
are

    N22=64 A_12+24=2508952,
    N37=64 A_6-24=199148648.

The defining products independently explain their errors. Raising to24
and expanding through degree3 gives

    64g_n^24=q^(-1)-24+276q-2048q²+O(q³),
    64G_n^24=q^(-1)+24+276q+2048q²+O(q³).

Here276=binomial(24,2); the extra contribution to2048 is the first q³
factor. Inverting the series supplies4096q because64²=4096. Consequently

    64(g_n^24+g_n^(-24))=q^(-1)-24+4372q+96256q²+O(q³),
    64(G_n^24+G_n^(-24))=q^(-1)+24+4372q-96256q²+O(q³).

The constants therefore have an explicit product/trace origin. No Monster,
dimension, or tournament inference follows from the isolated number276.

The checker evaluates the products separately from the surds at120-digit
working precision. It obtains

    exp(pi sqrt22)=2508951.998257424467165529...,
    exp(pi sqrt37)=199148647.999978046551856766....

The errors below their respective integer targets are about0.00174257553
and0.0000219534481. These decimal strings are numerical checks, not interval
certificates. The exact integer targets and recurrence identities do not
depend on floating-point arithmetic.

## M2. Simple nested square roots do not imply a near-integer exponential

Let a=2+sqrt5 and b=(3+sqrt13)/2. Both have norm-1. The four conjugates of
(ab)^12 have sum A_12(4,-1) A_12(3,-1), an integer. Besides (ab)^12 and
(ab)^(-12), this sum contains (a/b)^12 and (b/a)^12. One of these is larger
than1. They are not exponentially small merely because the expression for
g130 is short.

The resulting full-trace target and the actual exponential are

    64 A_12(4,-1) A_12(3,-1)+24 =3599845752586520,
    exp(pi sqrt130)             =3599845752585248.574872003930....

The difference is about1271.42512799607. Its main term is exactly the
omitted-conjugate expression64((a/b)^12+(b/a)^12), numerically
1271.4251279960687.... Thus short radical syntax and excellent logarithmic
approximations to pi do not supply near-integrality of this exponential.
To use a trace as an integer approximation, retain every other conjugate
and a bound on their total. This is the arithmetic version of keeping the
joint component data before scalarizing.

### M2a. The same 5,13 allocation square gives a literal Pythagorean bridge

There is a precise connection with the older
[marked Gaussian/Pythagorean fiber](../endpoint_pythagorean_fibers_20261019.md),
not merely the repeated number 65. Fix Gaussian factors 2+i above 5 and 3+2i
above 13. Their four independent conjugation choices, modulo simultaneous
conjugation, give the two primitive triples with hypotenuse 65:

    (2+i)(3+2i)=4+7i  -> (33,56,65),
    (2+i)(3-2i)=8-i   -> (16,63,65).

In the g130 calculation fix a=2+sqrt(5), b=(3+sqrt(13))/2. The two independent
real conjugations act by a->-1/a and b->-1/b. On the even power (ab)^12,
simultaneous conjugation is inversion. Its four choices therefore give two
reciprocal-pair traces:

    K_same=(ab)^12+(ab)^(-12),
    K_opposite=(a/b)^12+(b/a)^12.

Matching the independent 5- and 13-conjugation bits gives an explicit
equivariant bijection between these two marked two-element quotients. It
preserves which prime coordinates are conjugated and the common-conjugation
quotient. Under this map the second Pythagorean allocation corresponds to
the cross-conjugate term that prevents near-integrality at 130.

This is a shared Boolean allocation carrier. It is not a field isomorphism:
the Gaussian construction lies in Q(i), while the real surds generate
Q(sqrt(5),sqrt(13)). The map does not preserve lengths, numerical traces,
or the square-sum positivity cone. Prime factors, chosen embeddings and
exponents remain required side data. The full trace is K_same+K_opposite;
discarding the other allocation is exactly the lost coordinate.

## M3. Doubling the radicand changes the CM field

The point tau=i sqrt(n) has endomorphism order Z[sqrt(-n)] of discriminant
-4n for the squarefree n in the quoted list. Replacing n by2n changes it
to discriminant-8n and changes Q(sqrt(-n)) to Q(sqrt(-2n)). These are
different imaginary quadratic fields. This is not a2-isogeny within one
fixed CM field: isogenous complex elliptic curves have isomorphic rational
endomorphism algebras. By contrast n->4n stays in the same field and changes
the order/conductor. Obtaining two evaluations by specializing one modular
equation does not make their fields identical.

Exact reduction of primitive positive binary quadratic forms gives:

| n | h(-4n) | h(-8n) | class group at -8n | classes per genus |
|---:|---:|---:|---|---:|
|23|3|4|C4|2|
|29|6|2|C2|1|
|31|3|8|C8|4|
|47|5|8|C8|4|
|71|7|4|C4|2|
|79|5|8|C8|4|

The complete form lists are in receipt.json. In each doubled case exactly
two classes are self-inverse (including the identity), checked by the
standard reduced-form boundary rules. A finite abelian2-group with exactly
two such elements is cyclic; hence the group column follows, not just its
order. The genus quotient A/A² retains only one bit in each row, while
classes per genus can still be2 or4.

This directly continues the repo's
[idoneal information theorem](../idoneal_information_channels_20261021/IDONEAL_INFORMATION_CHANNELS.md):
being a2-group is weaker than having exponent2. The C8 rows require a tower
of quadratic choices, and genus characters alone cannot choose the class.
Compare the inherited discriminant-1092 example with group C2³: it has the
same eight classes but its genus data is already complete. Equal class
counts do not identify composition laws or missing information.

### A decisive guard: class number a power of two is not sufficient by itself

Raw Weber values include normalization and level data. For example the
same primary paper gives g18=(sqrt2+sqrt3)^(1/3), while the order of
discriminant-72 has class number2. This g18 is not obtainable from rational
numbers by a tower of square roots.

Proof: put v=g18². Then v³=5+2sqrt6=epsilon, and sqrt6=(v³-5)/2 lies in
Q(v). The ring of integers of Q(sqrt6) is Z[sqrt6]. Its smallest real unit
of norm1 greater than1 is epsilon: for x+y sqrt6 with y>0, the equation
x²-6y²=1 first has a solution at y=2,x=5; smaller positive units would
require y<2. If epsilon were a cube in this field, its positive real cube
root would be an integral norm1 unit strictly between1 and epsilon. This
is impossible. Hence X³-epsilon is irreducible over Q(sqrt6), and
[Q(v):Q]=6, which cannot be a divisor of a power of two. Therefore g18
cannot belong to any iterated quadratic extension.

The missing data is the actual modular function and its multiplier/level
action, not merely a class-number table. [Schertz, *Weber's class invariants
revisited*, Theorems1--2 and Proposition2](https://numdam.org/item/JTNB_2002__14_1_325_0.pdf)
explicitly retains normalization powers, congruence conditions and eta
multipliers. We do not infer the rational minimal degree of raw G_n from
Ramanujan's cubic/quintic/septic wording alone. The listed surd evaluations
are inherited evaluations; this table is supporting structure, not a new
derivation of all six of them.

## M4. A useful arithmetic composition carrier

A proposed exact carrier for this lane should retain the quadratic order
(including conductor), the chosen modular function and level, its actual
Galois action, a selected root/embedding, and the bounds on discarded
conjugates used by the numerical target. These are separate obligations:

* recurrence traces can be integral while forgetting branch alignment;
* genus characters can be complete for C2³ and incomplete for C8;
* a class group can be a2-group while a raw modular value still needs a cube;
* simple surds can have large unwanted conjugates.

No minimality theorem for this proposed modular carrier is claimed. A
productive next target is an explicit, proved selector criterion for the
odd part of the modular multiplier action, followed by certified trace-error
bounds. Neither would by itself prove infinitely many Lucas prime values.
