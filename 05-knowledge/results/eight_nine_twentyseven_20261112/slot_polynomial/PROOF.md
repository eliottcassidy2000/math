# The slot polynomial: a genuine 8/9/27 precision boundary

2026-11-12. **PROVED:** the identities, exact residue observer, and padding
family below. **INHERITED:** the general binomial period law from
**THM-4027**, [sun-two-four-six-eight-universal-modular-solubility](../../../../01-canon/theorems/THM-4027-sun-two-four-six-eight-universal-modular-solubility.md).
The connection concerns finite arithmetic observers, not Collatz termination
or universal Hamilton-count coverage.

## 1. A polynomial with a literal meaning at every nonnegative integer

Let T be a nonempty n-vertex tournament, and let C_k(T) count ordered
spanning covers by k nonempty directed paths. Define

    F_T(x)=sum_(k=1)^n C_k(T) binom(x,k).

For integer q>=0 this counts ordered q-slot spanning path covers in which
empty slots are allowed. Choose the k occupied slots and then their ordered
nonempty paths. Equivalently, assign each vertex one of q colors and choose
a Hamilton path inside each induced color class.

The latter description, followed by the inherited odd-cycle expansion in
each color, gives

    F_T(x)=sum_(disjoint odd directed cycle packings P)
             2^|P| x^(n-sum_(C in P)(|C|-1)).             (1)

Every cycle in the packing must be monochromatic, while unused vertices
have arbitrary colors. This proves (1) first for all q>=0 and hence as a
polynomial identity. In particular F_T is monic with integer coefficients,
has no constant term, and obeys

    F_T(-x)=(-1)^n F_T(x),    F_T(1)=H(T).

These are the same retained cover coordinates as the inherited tournament
slot/odd-cycle polynomial. We do not assert a new OCF theorem.

Its integer fixed divisor is exactly the ordered-cover gcd

    h(T)=gcd_k C_k(T)=gcd_(q in Z) F_T(q).                (2)

One direction follows from the binomial expansion; the other from
C_k=Delta^k F_T(0), an integer combination of F_T(0),...,F_T(k).
Thus G_T=F_T/h(T) is integer-valued, though its power coefficients need
not be integers. Dividing power coefficients modulo a prime would lose
this distinction.

## 2. Bounded degree gives a finite precision observer

**Degree/precision theorem.** For every integer-valued polynomial G of
degree at most n>=1, every prime p and s>=1, the function G(x) modulo p^s
has period dividing

    p^(s+floor(log_p n)).                               (3)

Expand G in its integer binomial basis. For k<=n the denominator-safe
period theorem of THM-4027 gives period p^(s+floor(log_p k)); all these
divide (3). More directly, with P as in (3), Vandermonde gives

    binom(x+P,k)-binom(x,k)
       =sum_(j=1)^k binom(P,j)binom(x,k-j).

Since k<P and v_p binom(P,j)=v_p(P)-v_p(j)>=s, every term vanishes.
Generalized binomial coefficients at negative integer x are integral,
so the same proof covers all x in Z.

Apply this to G_T. Degree at most8 means its normalized mod3 observer
needs **at most two ternary digits**, a period dividing9. At degree9,
the general bound becomes27. The change is sharp for integer-valued
polynomials: binom(x,9) modulo3 has exact period27 by the same inherited
period theorem. This witness is not claimed to be a tournament profile.

Thus 8=3^2-1 and9=3^2 genuinely mark an observer-precision boundary.
It is distinct from the minimum tournament order permitting a context-gcd
gap. The content lane exhibits a formal degree7 cover polynomial with the
9/27 gap, so the latter minimum cannot follow from degree/parity alone.

## 3. The asymmetric core realizes exact two-digit memory

For the inherited core with outgoing rows

    (56,117,129,22,100,76,141,59),

the ordered covers are

    (513,6102,32508,99432,187200,216000,141120,40320),

and direct binomial conversion gives

    F(x)=x^8+34x^6+232x^4+246x^2.                       (4)

In particular F(1)=513, F(2)=7128 and F(4)=268128. Their gcd is9.
The input4 is a3-unit, so the loss of the extra factor3 is not confined
to nonunit inputs.

**Exact residue theorem.** For every integer x,

    F(x) mod27 = 18 if x=4 or5 mod9, and0 otherwise.    (5)

If 3 divides x, (4) is0 mod27. Otherwise put x^2=1+3t. Expanding gives

    F(x)=9x^2 (3t^3+37t^2+101t+57).

The parenthesis modulo3 is t(t+2), which is2 exactly at t=2 mod3 and
zero at t=0,1. Among unit residues modulo9, x^2=7 exactly at x=4,5.
This proves (5), including negative x.

Therefore F/9 modulo3 has **least positive period9**. Period9 follows
from (5), while periods1 and3 fail by comparing x=1 and4; any smaller
period together with9 would force their gcd to be a period. It records
the second ternary digit of the unit square class. The uniform-context
ideal27 is a different observer: its factorially weighted coordinates
annihilate precisely the raw-cover precision that this slot reading sees.

For completeness, on nonzero multiples of3 the valuation is exact:

    v_3 F(x)=2v_3(x)+1.                                (6)

Indeed the factor x^6+34x^4+232x^2+246 has valuation1: its constant246
has valuation1 and every other term has valuation at least2.

## 4. An infinite asymmetric family retains exactly9 versus27

Let T_m=C→TT_m be the ordinal sum of C and a transitive suffix of length
m>=0. Color each suffix vertex arbitrarily. Within each color its suffix
vertices have a unique order and append after that color's C-path. Hence

    F_(T_m)(x)=x^m F_C(x),    H(T_m)=513.               (7)

Every C_k(C) is divisible by9, so F_C(q) is for every integer q; (7)
proves9 divides h(T_m). Conversely the three values at1,2,4 have gcd9
for every m. The first anchors possible primes to3 and19; multiplication
by powers of2 cannot introduce either prime, and the other two already
remove19 and bound the3-adic valuation by2. Thus h(T_m)=9.

Every uniform context containing T_m is also a uniform context containing
the original C, so27 divides g(T_m). The general sharp content bound in
the companion content lane gives g(T_m)|h(T_m)^2/rad(h(T_m))=27.
Consequently

    h(T_m)=9,    g(T_m)=27    for every m>=0.            (8)

The core is strong and asymmetric, as inherited from its literal bank
record. Its SCC is distinguished from each suffix singleton, so an
automorphism fixes every suffix vertex and restricts to the identity on C.
Every T_m is therefore asymmetric. A strong claim is made only at m=0.
This is an infinite family at every order at least8, not an inference
from sampled paddings.

More generally, take r>=1 ordered copies of C and distribute m>=0
singleton strong components among the r+1 gaps. The polynomial is
`x^m F_C(x)^r`, independent of that distribution. The values at1,2,4
again give raw fixed divisor9^r. The first minimum3-adic evaluation is
still at4, so the first-minimum invariant in
[content/PROOF.md](../content/PROOF.md), Theorem A, bounds the context
divisor above by3*9^r. The values at0,...,3 are divisible by3^(3r),
while the fourth has valuation2r; all higher Newton terms gain at least
one factor3 from their excursion factorial. Therefore, for every r,m,

    H=513^r,       h=9^r,       g=3*9^r.             (9)

At fixed r,m there are exactly binom(m+r,r) unlabelled tournaments in
this constructed family: the unique ordered SCC decomposition recovers
the padding distribution. They are asymmetric, so their labelled
contribution at n=8r+m is `n! binom(m+r,r)`. This is a contribution to
the H513^r fibre, not the total number of tournaments in that fibre.
Its bivariate family generating function is

    sum_(r>=1,m>=0) binom(m+r,r) u^r z^(8r+m)
      = u*z^8 / ((1-z)*(1-z-u*z^8)).                (10)

This exact rational description makes the infinite padding family finite
to manage while retaining the context obstruction at every order.

## 5. The cross-problem connection and its boundary

The old Sun2/4/6/8 work supplied the exact binomial precision law despite
its global representation failure. Here that law makes a whole normalized
tournament observer finite. It preserves residue values, not the existence
of a tournament for an arbitrary integer-valued polynomial.

The history lane recovers a separate Collatz branch whose next coarse
observer requires one more ternary digit: dividing a known affine numerator
by3 consumes one input digit of precision. Formula(5) also depends on a
quotient by3, namely t=(x^2-1)/3, but these are distinct maps. Their shared
method is to retain the predivision residue and branch label. Neither the
equality of some small numbers nor the common modulus identifies their
orbits or proves termination.

The square-sum endpoints8 and9 instead arise from a removed fixed point
and an interval boundary in the9/16/25 reflection matchings. That supplies
the complementary geometric lesson used by the selective-port theorem:
retain the actual boundary or compound block before permuting it. The
history lane gives the exact maps and the original source paths.
