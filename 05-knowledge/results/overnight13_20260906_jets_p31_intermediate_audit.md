# Independent audit: universal Hermite packets and the p31 intermediate ideal

**Status: PASS, proof-level independent audit.** The universal minimum-odd-minor
identity and the all-depth, all-lift p31 formula in the
[producer report](overnight13_20260906_jets_p31_intermediate.md)
are accepted. They may be promoted to **PROVED ANALYTICALLY + EXACT POLYNOMIAL
CERTIFICATE + INDEPENDENTLY AUDITED**. The specified complete Smith matrices
and kernel orders remain FINITE-EXACT controls. No full p31 partition law or
general-prime classification follows from this audit.

The frozen producer source SHA256 is
`5dc0b2cb3af6e6e282c17eae5df81a2bc68117c6a29563e813cffd82fdc63525`;
its output SHA256 is
`efee820e2d67773205f14dbf3a7eed3d2361cdf67ac79c11241f6a34a2ffbda2`.
Both hashes were checked. No producer or repository files were changed.

## 1. Scope and retained coordinates

The object is the complete order-zero-through-order-15 Hasse bank on
`0,31^e,31^e a`, with `a,a-1` units in Z_31. There are 48 coefficient columns,
for polynomial degrees below 48. Clearing the node-zero rows by integral
unit operations leaves 16 identity rows and the residual matrix

```
B_(x,j;c) = 31^(e(c-j)) binom(c,j) x^(c-j),
x in {1,a}, j=0,...,15, c=16,...,47.
```

For e>=1 all residual entries are divisible by 31, so exactly the first
16 Smith exponents are zero. The sum E_r of the first r residual exponents
is the valuation of the full determinantal ideal D_(16+r). This identifies
the actual consumer: E29 is v31(D45), not the largest exponent. At e=0
the three nodes are distinct modulo 31, the confluent determinant is a
unit, and every exponent is zero. The positive-depth formula is not applied
at that boundary.

Translation and affine unit changes are harmless only because the derivative
banks are complete. The proof uses unit monomials in a and a-1, so the
equilateral residue restrictions are necessary. It covers arbitrary p-adic
lifts, including normalized rational units, not only the printed integer
residue representatives.

The inherited ninth Deuring theorem determines the largest exponent and the
eleventh p11/p13 laws motivate retaining intermediate minor content. Their
role is correctly separated from the present proof. The new ideal formula
does not depend on assuming that a Deuring bit or an AP bit determines the
full Smith form. The literal a3/a4 hostile refutes precisely that shortcut.

## 2. The universal cofactor identity, with signs and scalar type

Take integers s>=0 and m>=s+1. For residual rank r=2s+1 put h=m-s-1.
The least column sum uses m,...,m+2s. The largest derivative-order sum uses
both copies of h+1,...,m-1 and one copy of h. Thus there are exactly two
least-weight minors. Their common weight is

```
(2s+1)(m+s)-[(2s+1)h+s(s+1)] = 3s^2+3s+1.
```

First set m=s+1. In the full 3m-by-3m Hasse matrix, delete the last
coefficient column and the order-zero row at one nonzero node. Clearing
the node-zero identity block gives precisely the two residual minors.
The full confluent determinant is `a^(m^2)(a-1)^(m^2)` in the usual node
and derivative order. The highest coefficient of the Hermite basis
polynomial for evaluation at node a is

```
[z^s](a+z)^(-s-1)(a-1+z)^(-s-1)
  = (-1)^s T_s(a)/[a^(2s+1)(a-1)^(2s+1)],

T_s(a)=sum_(j=0)^s binom(s+j,j) binom(2s-j,s-j)
                              a^(s-j)(a-1)^j.
```

This follows directly from the two reciprocal binomial series. It is a
rational-function identity; the finite numerical controls do not supply
its quantifiers. At node one the corresponding numerator is
`U_s(a)=a^s T_s(1/a)`. Multiplication by the confluent determinant and the
cofactor signs gives the following exact signs when rows are ordered first
at node one, then at node a, with increasing derivative order:

```
M1 = K_(m,s) a^(s^2)(a-1)^(s^2) q_s(a),
M2 = (-1)^s K_(m,s) a^((s+1)^2)(a-1)^(s^2) a^s q_s(1/a),
q_s=T_s/content(T_s).
```

M1 keeps the lower derivative order at node one; M2 keeps it at node a.
The producer's up-to-sign version is therefore correct, including odd s.

For general m the Hasse identity

```
binom(c,j) = [(c)_h/(j)_h] binom(c-h,j-h)
```

lowers every selected degree and derivative order by h and leaves c-j
unchanged. The column and row scalars give exactly the producer's formula
for K_(m,s). The ratio before multiplying by content(T_s) is rational in
general: at s=2,m=4 it is 560/3. It is not automatically a p-adic unit.
The final K_(m,s) is integral by Gauss's lemma: the determinant is an
integer polynomial and the displayed monomial times primitive q_s is
primitive. This proves both the universal transport and its precise loss
of scalar-valuation information.

The companion independently verifies the signs through rational Gaussian
determinants, rather than the producer's absolute-value Bareiss checks.
It also recomputes the reciprocal Taylor coefficient by multiplying the two
truncated rational series, independently of the determinant path.

## 3. Exact p31 packet and all-lift attainment

The direct integer expansion gives content(T14)=19380. Its primitive
coefficient list agrees coefficientwise with the producer. At m=16 the
shift gives

```
K=19380*44!/(15!)^3=23038504627568008043520,
v31(K)=1.
```

Write q1=q14, q2=a^14 q1(1/a). The independent coefficient calculation
checks the polynomial identity

```
q2-a^2 q1 = 31 R,             R in Z[a].
```

It also exhausts all 29 admissible residues, giving common root set
`{3,11,15,17,21,29}`. The R values at these roots, in that order, are
`11,3,17,20,14,28`, all nonzero modulo 31. Thus for every lift of a root,
both q1 and q2 are divisible by 31 but their displayed combination has
valuation exactly one. Their minimum valuation is exactly one. Off the
root set at least one is a unit. This is an all-lift argument, not a
conclusion from checking finitely many second digits.

The six roots form the full anharmonic orbit of 3. Direct evaluation of
`sum_(j=0)^15 binom(15,j)^2 a^j` verifies that every one is Deuring-ordinary;
none is in the AP set {2,16,30}. The node-permutation invariance asserted
for kappa is correct. No supersingular classification or external priority
claim is needed here.

The individual-cancellation hostile is retained: at a=879 the valuations
of q1 and q2 are 3 and 1. A single packet's valuation is not the ideal
valuation. Because a and a-1 are units, the monomial factors in the two
minors add no valuation. Consequently their exact common valuation is

```
631e+1+kappa,
kappa=[a mod31 in {3,11,15,17,21,29}].
```

## 4. Lower bounds for every omitted minor

The least column sum is 870 and the largest sum of 29 row derivative
orders is 239. A minor's exponent factor is therefore
`31^(e(631+u+v))`, where u is its nonnegative column-sum excess and v is
its nonnegative row-sum deficit. Polynomial coefficients are integral.
This gives a complete, disjoint classification by u+v.

There are two minimal row sets and one minimal column set. At row deficit
one the omitted derivative orders must be (0,0,2) or (0,1,1); each has
two choices of the nonduplicated node. At column excess one the last
column 44 is replaced by 45. Thus exactly six minors have total excess one.

For all cases retaining only derivative orders j>=1, column 31 vanishes
identically modulo 31. In the remaining row pattern there is exactly one
order-zero row and no order-one rows. Both columns 31 and 32 then vanish
on every other row, since their binomial coefficients are divisible by 31
for 2<=j<=15. These two columns have support in the same single row, so
the determinant is zero over F31[a]. This proves coefficient divisibility
by 31, rather than merely vanishing at all scalar residues.

Every excess-one minor consequently has valuation at least 632e+1, and
every higher-excess minor has valuation at least 633e. For e>=1 both are
at least 631e+2. The two minimum minors attain 631e+1+kappa, so no omitted
minor can lower the ideal. This handles the delicate shallow case e=1,
kappa=1 as well as every larger depth. The proof establishes both lower
bound and attainment for the complete determinantal ideal.

The full index universe has `binom(32,29)^2=24,601,600` minors. The proof
does not enumerate those determinant polynomials: only the two minimal
and six next-band structures need coefficient information. The companion
enumerates the 4,960 row and 4,960 column index sets separately to verify
the band classification. It does not claim a generic next-band divisibility
principle at other primes or ranks.

## 5. Independent complete-matrix controls and actual observer consequence

The companion materializes the original integer 48-by-48 Hasse matrix,
then uses fixed-precision elimination with a globally least-valued pivot.
It normalizes a pivot only by a unit and clears its column with integral
DVR multipliers. This differs from the producer's repeated common-layer
division and unit-peeling implementation. Precision is 48e+3, exceeding
the displayed exponents, with explicit pivot visibility and determinant
checks. No producer module is imported.

There are 32 fresh full-matrix controls: residues 2,3,4,11,16,29,30 at
depths 0,1,2,3, plus four negative or higher-digit lifts of residue 3.
They recover the full D45 formula and determinant valuation 768e.
The original full matrices at a=3 and a=4, e=1, agree outside residual
positions 29 and 30, where their factors are respectively

```
a=3: 43,43,
a=4: 42,44.
```

Both largest exponents are 47. Capping all 48 valuations at precision 43
gives sums 762 and 761, hence kernel orders 31^762 and 31^761. This confirms
the producer's actual finite observer consequence by a separate full-matrix
path. The kernel formula is `log_31 |ker mod31^b|=sum_i min(b,lambda_i)`;
using the largest exponent alone would lose the distinction.

These controls do not prove a full partition law at every depth or lift.
That limitation is stated correctly in the producer. The all-depth result
accepted here is E29 alone, together with the general minimum-band packet
identity.

## 6. Reproduction and promotion basis

The independent source checks 44,617 always-active gates, including:

* 84 signed rational determinant controls at s=0,...,6 and
  m=s+1,s+2,s+4, plus direct Hermite reciprocal coefficients;
* six literal 29-by-29 determinants at a=2,3,4, with rational elimination;
* every coefficient of the divided-companion identity, every admissible
  residue, 116 rational/higher-digit inputs, and the a879 hostile;
* the complete least two row/column bands and the polynomial rank argument;
* 32 independent full original Hasse matrices and the finite kernel split.

Run:

```
python -B 04-computation/overnight13_20260906_jets_p31_intermediate_audit.py
python -B -O 04-computation/overnight13_20260906_jets_p31_intermediate_audit.py
```

Both modes have identical LF output. Independent source SHA256:
`476d544142ccd82963149407d7aea1f7d5b4df6f53d1bca2768e5aec23739914`.
Independent output SHA256:
`857da5af360c621e0f19be27749ddbf96bb6e81ad97373ecf2462ee310292a90`.

Promotion rests on the universal cofactor/shift derivation, exact finite
coefficient and residue certificates, the all-lift unit argument, and the
complete omitted-band lower bound with attained minimum. The full-matrix
bank is independent corroboration, not a substitute for those proofs.
No repair to the frozen producer is requested.

## 7. One precise next implication, separate from the audited theorem

The cheapest next target is the adjacent even ideals

```
E28=588e+2,             E30=675e+1.
```

If these are proved, the audited E29 immediately yields the actual residual
factor pair

```
lambda29=43e-1+kappa,   lambda30=44e-kappa.
```

The proposed proof is small: rank 28 has one minimal minor, using rows
j=2,...,15 at both nodes and columns 16,...,43. Two derivative shifts from
the full multiplicity-14 determinant give scalar valuation two. Every
weight-one modification retains column 31 and only j>=1 rows, hence has
an additional factor 31; higher weight already suffices. Rank 30 has one
minimal minor, using rows j=1,...,15 twice and columns 16,...,45, with
scalar valuation one. Its ordinary weight gap suffices at every e>=1.

This is a separately identified two-factor closure target. It is not part
of the producer's current promotion and has not been installed as a new
theorem here. Even if obtained, it would identify this pair's kernel
contribution, not assert that all other factors are independent of kappa.
No further census is needed to attempt it.

**Filing:** root integrated these independently audited artifacts in the thirteenth
checkpoint. Reproduction commands are relative to the repository root.
Transcript hashes refer to filed LF bytes; Windows CRLF captures were normalized.
