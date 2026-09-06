# A universal odd-minor packet and a p31 ideal invisible to two earlier residue bits

**Status: PROVED ANALYTICALLY + EXACT POLYNOMIAL CERTIFICATE + INDEPENDENTLY AUDITED.** The universal packet
identity and the stated p31 determinantal ideal are proved below by algebraic
arguments with exact finite polynomial checks. Full-matrix banks are
FINITE-EXACT controls. No full p31 partition law, general-prime classification,
or external priority claim is made.

## 1. Inheritance and the new missing coordinate

The [eleventh full p11/p13 laws](overnight11_20260906_smith_prime_banks.md)
retain the minimum and next-weight minor bands, integral coefficient content,
and a divided companion. They show that a Deuring bit alone need not determine
the full observer. The [ninth theorem](overnight9_20260906_jets_deuring.md)
still gives the exact largest factor at prime-order multiplicity. The
least-used operation is the derivative/degree shift already identified in
the eleventh note. Here it turns a smaller-multiplicity Hermite cofactor
into an intermediate minor of a larger complete observer.

The concept board is: complete Hasse banks; Hermite inverse cofactors;
degree/derivative shifts; integral content; shallow weight bands; and
simultaneous unit packets. The map preserves a normalized determinant
polynomial up to a specified rational scalar. It does not preserve the
whole Smith form, nor the scalar's p-adic valuation. Those are essential
boundaries of the connection.

The cheap probe tested actual full Hasse matrices at the next primes. It
found a p31 hostile: a=3 and a=4 have the same metric, the same ordinary
Deuring bit, and neither forms an arithmetic progression modulo31. Both
largest factors have exponent47 at depth one, yet an intermediate ideal
and the actual kernel modulo31^43 differ. We pursue that one ideal rather
than extrapolating an entire full-partition law from the finite bank.

## 2. The p31 statement

Consider Hasse derivative orders0,...,15 at each of the three nodes

```
0, 31^e, 31^e a,       a,a-1 in Z_31^*.
```

The source module consists of polynomials of degree below48. Translation
and affine unit changes on complete banks are unimodular. The first16
Smith exponents are zero. Let E_r denote the sum of the first r exponents
after those zeros, equivalently the valuation of the full determinantal
ideal of order16+r. Put

```
kappa=[a mod31 in {3,11,15,17,21,29}].
```

**Theorem.** For every e>=1 and every admissible p-adic lift a,

```
E_29 = v31(D_45) = 631e+1+kappa.                     (1)
```

At e=0 the entire observer is unimodular and D_45 has valuation zero.
The exceptional six residues form precisely the cross-ratio orbit of3:
they are stable under a,1/a,1-a,1/(1-a),a/(a-1),(a-1)/a. Thus kappa is
intrinsic under permuting the normalized nodes. None of them is Deuring
supersingular or an AP class. The latter classes are2,16,30 modulo31.

In particular, a=3 and a=4 have D_45 valuations633 and632 at e=1,
while both have largest exponent47 by the ninth theorem. Their complete
finite matrices have kernel orders31^762 and31^761 modulo31^43. The last
statement is FINITE-EXACT, certified by direct full-matrix algorithms;
(1) is an all-depth and all-lift statement about one intermediate ideal.

## 3. Universal normalized minimum odd minors

For integers m>=s+1 and odd residual rank r=2s+1, define

```
T_s(a)=sum_(j=0)^s binom(s+j,j) binom(2s-j,s-j)
                         a^(s-j)(a-1)^j,
c_s=gcd of the integer coefficients of T_s,
q_s=T_s/c_s,
h=m-s-1.
```

After clearing the m rows at node0, the unscaled residual matrix has entries
`binom(c,j) x^(c-j)`, with x=1 or a, j=0,...,m-1 and c=m,...,3m-1.
The minimum r-minors use columns m,...,m+2s and every derivative order
h+1,...,m-1 at both nodes, plus order h at one of them. Up to signs their
two normalized determinant polynomials are

```
K_(m,s) a^(s^2)(a-1)^(s^2) q_s(a),
K_(m,s) a^((s+1)^2)(a-1)^(s^2) a^s q_s(1/a),          (2)
```

where the exact rational expression for the scalar is

```
K_(m,s)=c_s * product_(c=s+1)^(3s+1) [(c+h)!/c!]
              / {h! product_(j=1)^s [((j+h)!/j!)^2}. (3)
```

The scalar in (3) is integral: the determinant in (2) is an integer
polynomial, and its monomial, (a-1) power and primitive q_s have primitive
product by Gauss's lemma. The shift ratio before multiplication by c_s can
be rational and is not asserted to be a p-adic unit.

**Proof of (2).** First let m=s+1. These residual minors are full-observer
cofactors deleting the highest coefficient column and one order-zero row
at node1 or node a. The confluent determinant is, up to sign,
`a^(m^2)(a-1)^(m^2)`. The top coefficient of the order-zero fundamental
Hermite polynomial at node a is

```
[z^s] (a+z)^(-m)(a-1+z)^(-m)
 = (-1)^s a^(-m-s)(a-1)^(-m-s) T_s(a).
```

This follows by multiplying the two finite reciprocal Taylor coefficients
that can contribute to z^s. Multiplication by the confluent determinant
gives the first formula with K=c_s. Interchanging the two nonzero nodes
and scaling by a gives the reciprocal packet in the second formula; its
monomial exponents also follow directly from the same cofactor expansion.

Now simultaneously raise m,c,j by one. The exact entry identity
`binom(c+1,j+1)=(c+1)/(j+1) binom(c,j)` leaves c-j unchanged. Factoring
these rational column and row scalars out of the determinant, then repeating
h times, gives (3). The indices become exactly the two minimum row sets
described above. This proves the universal identity, not only its observed
agreement at selected multiplicities.

## 4. The exact p31 packet and all-lift attainment

For s14 the coefficient content is c14=19380. In increasing degree order,
q14 has coefficients

```
2070, -44505, 449995, -2838430, 12489092,
-40589549, 100591491, -193344684, 290017026,
-338353197, 302737071, -201824714, 94976336,
-28310254, 4044322.
```

These integers are computed directly from the displayed binomial sum,
not fitted from numerical determinants. Put q1=q14 and q2=a^14 q1(1/a).
At m16, formula (3) becomes

```
K=19380*44!/(15!)^3=23038504627568008043520,
v31(K)=1.                                             (4)
```

Coefficientwise, q2-a^2q1 is divisible by31. Let

```
R(a)=(q2(a)-a^2q1(a))/31 in Z[a].                     (5)
```

Complete exact substitution in the29 admissible residues gives common
zero set precisely3,11,15,17,21,29. The values of R at those residues,
in the same order, are11,3,17,20,14,28. Therefore at every lift of an
exceptional residue both packets are divisible by31 but their displayed
combination has valuation exactly one. Off those residues one packet is
a unit. Thus, for every admissible p-adic a,

```
min(v31(q1(a)),v31(q2(a)))=kappa.                     (6)
```

An individual packet may cancel more deeply. For example a=879 has
v31(q1)=3 and v31(q2)=1. The divided companion explains why such a lift
does not deepen the common ideal. Since a and a-1 are units, the monomials
in (2) cause no further valuation loss.

## 5. The entire next band has a characteristic31 rank obstruction

The least column sum uses16,...,44. The largest sum of29 derivative
orders from two copies of0,...,15 is239. Their difference is631, giving
the minimum factor31^(631e). Exactly two row sets attain that weight.
Every other minor has integral excess weight at least one.

There are six minors of excess weight one. Two raise the last retained
column44 to45, keeping a minimum row set. The other four keep the first29
columns and lower the row-order sum by one. Equivalently, the omitted three
derivative orders are either(0,0,2) or(0,1,1), with two node choices each.
There are no other excess-one choices because the row and column losses
are nonnegative integers and their sum must be one.

Every one of these six normalized polynomials has all coefficients divisible
by31. If all retained derivative orders satisfy j>=1, column31 is zero
modulo31: binom(31,j) is divisible by31 for1<=j<=15. In the remaining case
there is exactly one order-zero row, and every other row has j>=2. Both
columns31 and32 are then supported only on that one row. Their determinants
vanish identically over F31[a]. For column32 this follows from Pascal's
identity and divisibility of binom(31,j),binom(31,j-1). This is a polynomial
rank argument, not a sample of residue evaluations.

The minimum minors have exact common valuation631e+1+kappa by (4)--(6).
Every next-band minor has valuation at least632e+1, which is at least
631e+2 for e>=1. Every omitted higher band has valuation at least633e,
also at least631e+2. Thus no omitted minor lowers the ideal, including at
the shallow boundary e=1,kappa=1. One of the two minimum minors attains
the claimed valuation. This proves both directions of (1).

The structural zero-column argument replaces a large symbolic next-band
bank. It uses the actual retained derivative orders and columns; a generic
claim that the entire next band always has the needed content is not made.

## 6. Exact controls and scope

The standalone source uses integer coefficients, exact fraction-free
determinants and modular full-matrix elimination. It checks the universal
packet against literal minors at s0,...,6 and three multiplicities each;
recomputes q14 and its content; verifies every coefficient of (5); checks
all899 admissible second-digit lifts; and retains six individual-cancellation
hostiles. It independently enumerates all complements of three rows and
columns to certify the complete minimum and next bands, then verifies the
polynomial zero-column/one-row argument exactly. Literal29-by29 minimum
determinants provide additional controls at a2,3,4.

The118 full original Hasse matrices include every residue at depths0,1,2,4,
with distinct higher lifts, plus the a3/a4 precision controls. Their D45
valuations agree with (1), and their largest factors agree with the inherited
ninth theorem. The direct full e1 partitions differ at residual positions29
and30:43,43 for a3 versus42,44 for a4. Capping every factor at precision43
gives762 versus761. These exact finite calculations demonstrate the actual
observer consequence; they do not establish a full all-parameter p31 law.

```
python -B 04-computation/overnight13_20260906_jets_p31_intermediate.py
python -B -O 04-computation/overnight13_20260906_jets_p31_intermediate.py
```

The next open question is which smaller-multiplicity cofactor packets, after
the derivative shift and its content change, control further intermediate
ideals. Metric, Deuring and AP data together are already insufficient.

The [independent proof audit](overnight13_20260906_jets_p31_intermediate_audit.md)
accepts the universal identity and complete all-lift ideal proof, with44,617
normal/optimized gates and an independent full-matrix algorithm. Root read
the proof and audit. Primary gates97,294 also pass in both modes.

Source SHA256: `5dc0b2cb3af6e6e282c17eae5df81a2bc68117c6a29563e813cffd82fdc63525`.
Output SHA256: `efee820e2d67773205f14dbf3a7eed3d2361cdf67ac79c11241f6a34a2ffbda2`.

**Filing:** root integrated these independently audited artifacts in the thirteenth
checkpoint. Reproduction commands are relative to the repository root.
Transcript hashes refer to filed LF bytes; Windows CRLF captures were normalized.
