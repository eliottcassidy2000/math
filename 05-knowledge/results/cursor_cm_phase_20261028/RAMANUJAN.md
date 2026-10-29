# The quoted radicals: a Fricke trace, weighted divisors, and a field change

**CITED** classical class-invariant evaluations and modular identifications;
**PROVED** elementary algebraic consequences below; **FINITE-EXACT** series,
forms and traces; **VERIFIED** certified decimal intervals. No new modular
identity or Moonshine theorem is claimed.

## R1. Recover the notation before comparing the numbers

The excerpt is from [Ramanujan, *Modular equations and approximations to pi*
(1914), sections1,4,5](https://ramanujan.sirinudi.org/Volumes/published/ram06.pdf).
With q=exp(-pi sqrt(n)), his definitions are

    G_n=2^(-1/4) q^(-1/24) product_(j odd>=1)(1+q^j),
    g_n=2^(-1/4) q^(-1/24) product_(j odd>=1)(1-q^j).

The garbled exponents in the prompt are g_n^24 and G_n^24, not subscripts
24n. The two longer evaluations are

    g_130 = sqrt((2+sqrt(5))*(3+sqrt(13))/2),

    G_65 = [((1+sqrt(5))/2)*((3+sqrt(13))/2)]^(1/4)
           * [sqrt((9+sqrt(65))/8)+sqrt((1+sqrt(65))/8)]^(1/2).

The values g_22=sqrt(1+sqrt(2)), G_37=(6+sqrt(37))^(1/4), and
g_58=sqrt((5+sqrt(29))/2) use the same normalization. The checker compares
all five radical expressions to their defining products to75 decimal places.
That is a numerical cross-check of cited identities, not their proof.

## R2. Why 24,64,276,4096,4372 belong together

Set P(q)=product_(j odd)(1-q^j)^24. In formal integer series,

    P(q)=1-24q+276q²-2048q³+11202q⁴+...,
    P(-q)=product_(j odd)(1+q^j)^24.

Here276=binomial(24,2); the q² term comes entirely from the first factor.
The normalization2^(-1/4), raised to24, gives1/64. Define

    t(q)=q^-1 P(q),
    F(q)=t(q)+4096/t(q)+24.

Then

    F(q)=q^-1+4372q+96256q²+1240002q³+...,
    4372=276+4096.

The eta product identity gives t=(eta(tau)/eta(2tau))^24 for q=exp(2pi i tau).
F is the classical Fricke Hauptmodul j_2^*=T_2A; see
[Zagier, *Traces of singular moduli*, section8](https://citeseerx.ist.psu.edu/document?doi=ad32ffe3debb0169cea3d102151a74b1edb76fed&repid=rep1&type=pdf)
and [Kim, *Borcherds products associated with certain Thompson series*](https://arxiv.org/pdf/math/0203076).
The exact connection to the quoted invariants is

    64(g_n^24+g_n^-24)+24 = F(q),
    64(G_n^24+G_n^-24)-24 = -F(-q).

The eta transformation under W_2:tau->-1/(2tau) interchanges t and4096/t.
The quotient loses exactly a reciprocal choice generically:

    t²-(F-24)t+4096=0.

On the positive real locus the two roots coincide at t=64; choosing t>64
selects the cusp branch. This is an exact two-valued algebraic quotient,
not a statement that every arithmetic trace loses only one bit. The
[ring-clock calculation](TRACE_CARRIER.md) supplies explicit larger fibers.

The deeper repo connection is recorded in
[MONSTER.md, section2](../cyclic_attachment_lattices_20261023/MONSTER.md).
The196883-dimensional Monster character restricts to2.B in dimensions
1,4371,96255,96256; a2A involution has trace4371. The additional trivial
dimension gives4372. The last constituent has negative central sign, so
it cannot be silently treated as a representation of B. This is established
representation data, credited there to
[Breuer, Magaard and Wilson, section2](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf).
It explains a shared coefficient, without identifying a tournament carrier
with a Monster module.

## R3. An exact arithmetic inverse behind the product

For P(q)=product_(d>=1)(1-q^d)^(a_d), write

    -q P'(q)/P(q) = sum_(n>=1) b_n q^n.

Formal differentiation gives b_n=sum_(d|n) d a_d. Mobius inversion recovers

    a_n = (1/n) sum_(d|n) mu(n/d)b_d.

In this example a_d=24 for odd d and0 for even d. Consequently

    b_n=24 sigma(oddpart(n)).

This is an actual map from the class-invariant product to weighted divisor
data. Both the complete product and the complete logarithmic derivative
retain every exponent; a finite coefficient truncation retains only a
finite initial segment. Replacing divisor weights by parity would lose
exactly what the inverse needs. The earlier Liouville/Mobius distinction
fits here: 1*mu is the convolution identity, whereas 1*lambda is the
indicator of perfect squares. They invert different observations.

There is also the formal identity

    P(q)P(-q)=P(q²).

It gives g_(4n)=2^(1/4) g_n G_n. The scale is4n because q(n)²=q(4n).
This is NOT the parameter change n->2n in the quoted comparison. Sign
switching, nome squaring and Fricke inversion are three distinct operations.

## R4. Integer traces explain the near integers; simple radicals do not suffice

For a norm-minus-one unit u with u+bar(u)=a, put S_k=u^k+bar(u)^k.
Then S_0=2, S_1=a and S_k=aS_(k-1)+S_(k-2). At even k,
bar(u)^k=u^-k, so S_k is the needed reciprocal trace.

| n | a | k | exact integer N | certified positive N-exp(pi sqrt(n)) |
|---|---:|---:|---:|---|
|22|2|12|2508952|0.00174257553283447...|
|37|12|6|199148648|0.0000219534481432335...|
|58|5|12|24591257752|0.000000177786758530424...|

The first/last integers are64S_k+24; the middle is64S_k-24. These integers
come from exact recurrences. The decimal enclosures use Machin's formula
with rational alternating-series bounds and outward-rounded Decimal
operations. They do not rely on the displayed asymptotic dots.

**Hostile130.** Put alpha=g_130². The two norm-minus-one factors in R1 give

    alpha+alpha^-1=6+sqrt(65),
    alpha^4-12alpha³-27alpha²-12alpha+1=0.

The reciprocal pair has lost only one branch, but its trace is still
quadratic irrational. With D_0=2,D_1=x,D_k=xD_(k-1)-D_(k-2),

    alpha^12+alpha^-12=D_12(6+sqrt(65))
      =28123794942082+3488327441280 sqrt(65).

Thus the trace responsible for F(q) is not an integer. The independently
certified interval gives

    exp(pi sqrt(130)) =3599845752585248.57487200393006006027... .

This explicitly refutes 'quadratic surds imply a near integer'. The repaired
criterion requires control of the other algebraic conjugates as well as
integrality of the relevant full trace. It does not classify all near integers.

The final incoming [modular comparison, M2-M2a](../reciprocal_lucas_carriers_20261028/MODULAR.md)
recovers those omitted conjugates explicitly. With a=2+sqrt(5) and
b=(3+sqrt(13))/2, the other reciprocal trace is (a/b)^12+(b/a)^12.
The FULL four-conjugate trace gives the integer3599845752586520 after the
64 multiplier and24 offset; the omitted contribution is about1271.425.
It also gives a precise marked bridge to the two Gaussian factorizations
(2+i)(3+2i)=4+7i and(2+i)(3-2i)=8-i, hence triples(33,56,65) and(16,63,65).
The map matches the two independent5/13 conjugation bits modulo simultaneous
conjugation. It preserves that allocation action, not a field or metric.

## R5. The easier g_(2n) changes the CM field

For the eta coordinates above, G_n uses tau=(1+i sqrt(n))/2 and g_(2n)
uses tau=i sqrt(2n)/2. For the squarefree odd n in this table, their
primitive quadratic discriminants are D_G=-n when n=3 mod4 and -4n
otherwise, versus D_g=-8n. Exact reduced primitive form counts give:

| n | D_G | h(D_G) | D_g | h(D_g) |
|---|---:|---:|---:|---:|
|23|-23|3|-184|4|
|29|-116|6|-232|2|
|31|-31|3|-248|8|
|47|-47|5|-376|8|
|71|-71|7|-568|4|
|79|-79|5|-632|8|
|65|-260|8|-520|4|
|37|-148|2|-296|10|

For the listed n=7 mod8, replacing D_G by the conductor-two discriminant
-4n happens to leave h unchanged; that does not make the orders identical.
The2-power orders in the doubled column explain why quadratic towers are
plausible in the corresponding class fields. A degree comparison by itself
does not prove a particular raw Weber invariant's radical expression:
normalizations, multiplier extensions and selected powers matter.
In particular the old phrase 'G_n is cubic/quintic/septic' must not be read
as an unqualified minimal polynomial degree of the raw normalized G_n over Q.

An explicit incoming hostile strengthens this warning:
[M3 of the reciprocal/Lucas modular note](../reciprocal_lucas_carriers_20261028/MODULAR.md)
uses g_18=(sqrt(2)+sqrt(3))^(1/3). The order discriminant-72 has class
number2, yet g_18² has degree6 over Q, proved using the fundamental positive
unit5+2sqrt(6). A raw Weber value therefore need not be a quadratic tower
even when the class number is a power of two.

The last row is the cheap hostile to a universal improvement: doubling37
introduces an odd class-number factor5. No universal monotonicity is claimed.

Two precise repo bridges survive:

- [THM-4825, half-shifted grid/class-group theorem](../../../01-canon/theorems/THM-4825-the-half-shifted-second-grid-is-a-ring-or-the-order-two-class-over-two-the-two-grids-are-the-class-group-iff-d-is-5-13-or-37.md)
  identifies the two-grid class group exactly at d=5,13,37, for d>1 and
  d=1 mod4. The classes are the principal class and the ramified ideal over2.
  This explains an algebraic role for the user's trio; the real unit in R4
  is not literally a point in the same imaginary norm lattice.
- [THM-4666, Kronecker ranking excess](../../../01-canon/theorems/THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md)
  gives excess p*h(-p) for the natural Paley ordering. Exact character sums
  at p=23,31,47,71,79 give69,93,235,497,395. Composite discriminants instead
  have zero character values on nonunits; those are missing relations, not
  arcs to orient arbitrarily. Nothing here proves the natural ordering optimal.

[cm_series.py](cm_series.py) computes the table by reduced forms, not by a
preloaded class-number list. Its two independent eta-product constructions
agree through q^32. See [cm_series.json](cm_series.json) for all exact data.
