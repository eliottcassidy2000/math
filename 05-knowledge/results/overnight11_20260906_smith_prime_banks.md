# Two full prime-order Smith laws and a Deuring-blind intermediate ideal

**Status: PROVED by analytic reduction and exact finite polynomial
identities; independently audited.** Literal matrix banks are
**FINITE-EXACT controls**, not an extrapolation argument. No theorem
identifier, repository edit, or external priority claim is requested.

## 1. Inheritance, target, and the first decisive test

The observer consists of Hasse derivatives of orders 0,...,m-1 at three
nodes on polynomials of degree below 3m. We retain its entire p-local
Smith partition. The closest proved mechanisms are the
[eighth full p7 law](overnight8_20260906_jets_residue.md)
and the [ninth exact Deuring loss](overnight9_20260906_jets_deuring.md).
The latter determines the largest factor at prime-order multiplicity
m=(p+1)/2, including all higher lifts. It does not determine intermediate
determinantal ideals. The canonical hostile is the earlier ternary nested
four-node two-jet metric failure; the corrected near miss here is the
inference that an exact top factor automatically reconstructs the partition.

The least-used sidecar is the *next minor-weight band*. The p7 proof needed
only a one-digit gap, whereas some current ideals need two digits. At
depth e=1, checking only the smallest-weight minors would miss a genuine
lower-bound obligation. The concept board is complete jet banks, residual
weights, integral coefficient content, simultaneous unit packets, divided
companions, and actual finite-precision kernels. The LRC anchor retains
endpoint owners and the no-three-line niche retains actual permutations;
this wildcard retains intermediate ideals. These are methodological
comparisons, not theorem transfers.

The first literal p11 bank supported a one-bit full law. The next prime
immediately supplied a hostile: at p13, m7, the nodes (0,13,26) and
(0,13,39) have the same metric, the same Deuring bit, and the same largest
factor, but different intermediate factors. The positive mechanism below
both explains the failure and recovers full laws at p11 and p13.

## 2. Statements, including every depth and unit lift

Let p be 11 or 13, m=(p+1)/2, and assume the three pairwise p-adic depths
all equal e. Translation and a p-adic unit change of variable are
unimodular on complete Hasse banks, so normalize the nodes as

```
0, p^e, p^e a,       a,a-1 in Z_p^*.
```

The parameter may be any p-adic integer satisfying these unit conditions;
it is not restricted to a residue representative. At e=0 every factor
is a unit. In the following formulas e>=1 and the displayed lists are
already nondecreasing.

**Theorem A: p11, m6.** Put `sigma=[a mod11 in {2,6,10}]`. The full list
is six zero exponents followed by

```
e, 2e, 4e, 5e, 7e, 8e+1, 10e+1, 11e,
13e-1, 14e, 16e-1+sigma, 17e-sigma.                 (1)
```

Here sigma is precisely the Deuring bit. Thus it suffices for the entire
partition at p11, not only its largest factor.

**Theorem B: p13, m7.** Put `tau=[a mod13 in {2,7,12}]`. The full list
is seven zero exponents followed by

```
e, 2e, 4e, 5e, 7e, 8e, 10e+1, 11e+1, 13e, 14e,
16e-1+tau, 17e-tau, 19e-1, 20e.                    (2)
```

All admissible p13 parameters are Deuring-ordinary. The additional bit
tau records that one of the three normalized residue nodes is the average
of the other two. Indeed these three possibilities are a=2, a=-1, and
a=1/2. Consequently tau is intrinsic under node permutations and affine
unit transformations. No higher unit digit is needed in (1) or (2).

The sum of (1) is 108e and that of (2) is 147e, agreeing with the confluent
determinant. The last pair in (1) ties at e=1,sigma=1; the affected pair in
(2) ties at e=1,tau=1. Other shallow ties in the displayed lists are
retained automatically, with no unchanged-prefix assumption imposed.

In particular, at e=1 the p13 lists are

```
a=2: (0^7,1,2,4,5,7,8,11,12,13,14,16,16,18,20),
a=3: (0^7,1,2,4,5,7,8,11,12,13,14,15,17,18,20).    (3)
```

Both largest factors are 20. Modulo 13^16 their observer kernels have
orders 13^141 and 13^140 respectively, since a Smith factor p^s contributes
p^min(s,16) to the kernel. This is an actual precision-observer separation,
not merely a change in a selected minor polynomial.

Relative to the proved smaller-prime routes, 13 is the least odd prime
at which the Deuring bit fails to determine the full partition in this
precise three-node equilateral family with m=(p+1)/2. For p3 the two-jet
[THM-4429](../../01-canon/theorems/THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md)
applies; p5 is covered by the
[seventh three-jet law](overnight7_20260906_oddjets.md);
p7 by the eighth law; p11 by Theorem A. This is not an integer-diameter
minimality statement or a claim for arbitrary multiplicity or node count.

## 3. Elimination, weights, and a reusable shift operation

The m Hasse rows at zero give an identity on columns 0,...,m-1 and vanish
on the remaining columns. Clearing their entries in the other rows is
integral and leaves the residual 2m by 2m matrix

```
R_((x,j),c)=p^((c-j)e) binom(c,j) x^(c-j),
x in {1,a}, j=0,...,m-1, c=m,...,3m-1.              (4)
```

An r-minor is `p^(We) P(a)`, with integral P and
`W=sum(columns)-sum(derivative orders)`. The least weight uses the first
r columns and the r highest derivative orders. It has one row choice at
even r and two at odd r. Write w_r for this least weight. For r=1,...,14,
as far as the matrix size allows, these weights are

```
1,3,7,12,19,27,37,48,61,75,91,108,127,147.          (5)
```

Every other minor is at least one weight higher. If an ideal needs a
correction of two, the next band is also relevant: it consists exactly
of the maximal-row choices with the last column increased by one, and
the first-column choices with total row order reduced by one. There are
five such minors at even r and six at odd r in the ranges used below.
Every omitted minor then has W>=w_r+2.

There is a useful operation explaining why the same small polynomial
packets reappear. Increase m, every retained derivative order j, and
every retained degree c by one. The powers c-j are unchanged, and

```
binom(c+1,j+1)=(c+1)/(j+1)*binom(c,j).
```

Thus the normalized determinant polynomial changes by the scalar
`product_columns(c+1)/product_rows(j+1)`. Its residue-root packet before
content reduction has the same shape, while its integral p-content can
change. This is an exact derivative/degree shift, not an assertion that
changing the observer preserves its Smith form. The separate source
checks the operation on all ten shared minimal ranks at m6 and m7.

## 4. Complete minimal-weight identities and the extra band

For r=2s the unique minimal polynomial is

```
K_(p,r) a^(s^2)(a-1)^(s^2).
```

For r=2s+1 its two polynomials are

```
K_(p,r) a^(s^2)(a-1)^(s^2) q_s(a),
K_(p,r) a^((s+1)^2)(a-1)^(s^2) (-1)^s a^s q_s(1/a). (6)
```

The first choice retains one more high derivative at node1; the second
retains one more at node a. The needed primitive packets are

| s | q_s(a) |
|---|---|
| 0 | 1 |
| 1 | 2a-1 |
| 2 | 7a^2-7a+2 |
| 3 | (2a-1)(3a^2-3a+1) |
| 4 | 143a^4-286a^3+234a^2-91a+14 |
| 5 | (2a-1)(26a^4-52a^3+44a^2-18a+3) |

Here are all constants required for the proof:

| r | K_(11,r) | v11 | K_(13,r) | v13 |
|---|---:|---:|---:|---:|
| 1 | 6 | 0 | 7 | 0 |
| 2 | 126 | 0 | 196 | 0 |
| 3 | 2940 | 0 | 8232 | 0 |
| 4 | 185220 | 0 | 1037232 | 0 |
| 5 | 740880 | 0 | 11409552 | 0 |
| 6 | 40748400 | 1 | 1882576080 | 0 |
| 7 | 32016600 | 2 | 6409723320 | 1 |
| 8 | 124864740 | 2 | 116656964424 | 2 |
| 9 | 252252 | 1 | 1767529764 | 2 |
| 10 | 756756 | 1 | 42420714336 | 2 |
| 11 | not needed | — | 80047968 | 1 |
| 12 | not needed | — | 17153136 | 1 |

These are exact integer identities. They can be verified by a short
Laplace expansion along the node1 rows: every summand is a product of two
integer binomial determinants and one monomial in a. The standalone
certificate stores every coefficient, row set, and degree set. It also
checks each polynomial against a separate literal Bareiss determinant.

At p11, none of the paired q_s packets for s<=4 vanish simultaneously
at an admissible residue. At p13 the same holds for s<=4. These are exact
finite-field statements, checked on all residues by the certificate. As
a transparent example, the p11 s4 packet is, up to units, the common
quadratic a^2-a+1, whose discriminant -3 is a nonsquare modulo 11. At p13
the first s4 polynomial reduces to the nonzero constant one. The sole
exception needed here is the p13 s5 packet, treated in the next section.

For all ranks requiring a two-digit correction, the entire next-weight
band has at least one p in its integer coefficient content:

| p | r | Number with content valuation 1 | Number with valuation 2 |
|---|---:|---:|---:|
| 11 | 7 | 4 | 2 |
| 11 | 8 | 4 | 1 |
| 13 | 8 | 0 | 5 |
| 13 | 9 | 0 | 6 |
| 13 | 10 | 4 | 1 |
| 13 | 11 | 6 | 0 |

This table is an exact coefficient-divisibility certificate, not a
sample of evaluations. Together with the minimum-weight identities it
contains just 26 p11 polynomials and 40 p13 polynomials. All other minors
are handled by (5) and the weight gap; their expansions are unnecessary.

## 5. The p13 intermediate divided companion

At rank11, remove the unit monomials in a and a-1 from (6), and put

```
q1=(2a-1)(26a^4-52a^3+44a^2-18a+3),
q2=(a-2)(3a^4-18a^3+44a^2-52a+26).
```

Modulo 13 their common admissible zero set is exactly {2,7,12}, and
`q2=-a^2 q1`. The integral divided companion satisfies

```
(q2+a^2 q1)/13
 =(a-1)(a^2-a+1)(4a^4-2a^3-a^2-2a+4).              (7)
```

At a=2,7,12, the right side has respective residues 2,1,12. Hence it is
a unit on every p-adic lift of each exceptional class. Both q1 and q2
are divisible by 13 there, but their displayed combination has valuation
exactly one. Consequently

```
min(v13(q1(a)),v13(q2(a)))=tau                       (8)
```

for every admissible p-adic a. A lift may cancel an individual numerator
more deeply; it cannot increase this common ideal. Restoring K_(13,11)
adds its single factor of 13, giving correction 1+tau at this rank.
This is a new intermediate-ideal use of a divided companion; the ninth
Deuring proof concerns a different polynomial pair and the largest factor.

## 6. Lower bounds, attainment, and reconstruction

Write E_r for the valuation of the rth residual determinantal ideal, so
E_r is the sum of the first r exponents after the initial m zeros. For p11,
the preceding identities give

```
(E1,...,E10)
 =(e,3e,7e,12e,19e,27e+1,37e+2,48e+2,61e+1,75e+1). (9)
```

For p13 they give

```
(E1,...,E12)
 =(e,3e,7e,12e,19e,27e,37e+1,48e+2,
   61e+2,75e+2,91e+1+tau,108e+1).                 (10)
```

Here is the lower-bound argument including the shallow boundary. A target
correction zero follows from W>=w_r. For correction one the minimal
polynomials are all divisible by p and every other weight costs at least
e>=1. For correction two the minimal polynomials have valuation at least
two, the next band costs at least e+1>=2, and all omitted weights cost
at least 2e>=2. At exceptional p13 rank11 use (8) for the minimal band
and the last row of the next-band table. Thus every minor satisfies the
claimed lower bound, including e=1.

Every bound is attained by a minimum-weight polynomial. The finite-field
unit packets prove this for all lifts at the ordinary ranks; (8) proves
it at exceptional rank11. This establishes the ideals, not only upper
bounds on their valuations.

The ninth theorem supplies the last two ideals. Its Deuring polynomial
has roots exactly {2,6,10} at p11 and no admissible root at p13, so

```
p11: E11=91e+sigma, E12=108e;
p13: E13=127e,      E14=147e.                       (11)
```

Successive differences in (9)-(11) prove (1)-(2). At e=0 the determinant
is a p-adic unit and handles that boundary separately. All uses of the
weight gap therefore retain their stated e>=1 hypothesis.

## 7. Computation, scope, and reproduction

The source has three separate exact paths: symbolic Laplace coefficients,
literal integer Bareiss minors, and full p-local elimination of the
original Hasse matrices. Minimum-valuation pivots have p-integral row
multipliers; clearing the corresponding row and column is unimodular over
Z_(p), so this last path computes all Smith exponents directly. It does
not consult the predicted ideal list while eliminating.

The 310 full-matrix controls comprise every admissible a modulo p^2 at
e=1 for p11 and p13, every residue with a negative third-digit lift at
e=0,2,5, and signed affine/node-permutation controls. Their determinant
valuations and sorted exponents are checked at each instance. The exact
polynomial packet and divided identity, rather than the depth bank, prove
the all-lift statements. The final certificate is a standalone JSON file
containing all 66 retained polynomials and their explicit index sets.

```
python -B 04-computation/overnight11_20260906_smith_prime_banks.py
python -B -O 04-computation/overnight11_20260906_smith_prime_banks.py
```

Both runs pass **65,923 always-active gates**, with byte-identical LF
output. Only the Python standard library is required.

```
source SHA256:
ed102132d643cbdea698f3705576dced1ac51f6a26cd9a4659390da5e5a453ff
output SHA256:
1ef4498d01078e2085c65c025a73a00c5a342ade1aed188b74e9a45bfc18d621
certificate SHA256:
05f79d3dba323be53a1308ea174f61cc2cb79dff053c61cd4c93ec80a15751ec
```

The failed proposed connection was metric plus Deuring bit to the *entire*
partition. Its strongest survivor is the ninth exact largest-loss theorem,
and Theorem A shows that it even suffices at p11. At p13 the lost coordinate
is the rank11 arithmetic-progression packet; Theorem B restores it with a
single intrinsic bit and an all-lift divided companion. No general-prime
claim that these two bits always suffice is made. The precise next question
is which universal minimum-weight packets acquire simultaneous roots at
larger primes, and how much of the next weight band their content requires.

The [independent audit](overnight11_20260906_smith_prime_banks_audit.md)
accepts both full laws with260,363 exact gates and no repair. The
[66-polynomial certificate](../../04-computation/overnight11_20260906_smith_prime_banks_certificates.json)
lives beside its producer and referee; the producer regenerates those pinned bytes.

**Filing:** root integrated these audited artifacts in the eleventh checkpoint;
reproduction commands above are relative to the repository root. Outside-worktree
locations preserve author provenance, not the present reproduction location.
