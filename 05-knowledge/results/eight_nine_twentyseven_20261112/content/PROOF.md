# Parity controls excursion content, and ordinal powers stabilize exactly

Date: 2026-11-12. **PROVED:** the parity-sensitive content bound, finite
first-minimum invariant, and all-power stabilization theorem below.
**FINITE-EXACT:** the inherited small-bank minimum-order control and the
named literal examples. No claim of full Chan–Pak realization is made.

## 1. Inheritance and the new invariant

For a nonempty tournament T on n vertices let C_k(T) count ordered covers
of its full vertex set by k nonempty directed paths. Put

    h(T)=gcd_(1<=k<=n) C_k(T),
    g(T)=gcd_(1<=k<=n) (k−1)! C_k(T).                    (1)

The second number is the exact gcd of Hamilton-path counts in **all
uniform one-module tournament contexts**. This is inherited from
[Session137, width, Theorems1–2](../../session137_overnight_proof_frontiers_20261102/width/PROOF.md)
and [Session293, MODULE_IDEAL.md](../../session293_cycle_returns_20261113/control/MODULE_IDEAL.md).
Permuting the k−1 actual nonempty exterior return excursions gives the
extra factorial; sharp triangular contexts prove equality of the ideals.

We retain the integer polynomial

    F_T(x)=sum_k C_k(T) binom(x,k).

Since C_k=k!p_k for integer unordered-cover counts p_k, F_T has integer
monomial coefficients. The inherited permutation/cut reversal identity is

    F_T(−x)=(−1)^n F_T(x),
    F_(T→U)(x)=F_T(x) F_U(x),
    F_T(1)=H(T).                                        (2)

See [Session139, all_moduli, §0](../../session139_incidence_moduli_and_stability_20261102/all_moduli/PROOF.md)
and Session137 §4 for the direct parity proof. Finite differences give
`C_k=Delta^k F_T(0)`, so h is also the fixed divisor
`gcd_(x in Z) F_T(x)`. The values x=0,...,n already generate that ideal.

Closest mechanism: the exact factorial-weighted context ideal. Hostile:
the asymmetric order-eight core has raw content9 but context content27;
automorphisms and the unweighted profile do not explain the extra factor.
Corrected near miss: ordinary unit-residue evaluation is not the full
context quotient at prime powers. Least-used sidecar: **where the raw
minimum p-adic valuation first occurs**, rather than its value alone.
Board: ordered path covers / excursion permutations / polynomial parity /
factorial plateaux / first minimum / physical ordinal powers.

## 2. The first-minimum index and the sharp parity bound

For a prime p dividing h, define

    a_p = v_p(h),
    kappa_p = min {k>=1 : v_p(C_k)=a_p},
    b_p = v_p((kappa_p−1)!),
    E(T) = product_(p|h) p^(b_p).                         (3)

If h=1, the empty product is1. A minimum exists since the cover profile
has finite support and C_1=H>0.

**Theorem A.** Every nonempty tournament satisfies

    h(T) | g(T) | h(T) E(T),
    E(T) | h(T)/rad(h(T)).                               (4)

In particular

    h(T) | g(T) | h(T)^2/rad(h(T)),                      (5)

and h and g have exactly the same prime divisors. If h is squarefree,
then g=h. For each prime with a_p>0,

    a_p <= v_p(g) <= a_p+b_p <= 2a_p−1.                 (6)

The raw divisibility h|g is immediate. First, without using parity, if
a prime p does not divide h, choose C_k not divisible by p. Since k!
divides C_k, necessarily k<p; thus (k−1)! is a p-unit and p does not
divide g. This proves equality of prime support. For p|h, the kappa_p
coordinate proves `v_p(g)<=a_p+b_p`.

It remains to prove the strict bound `b_p<=a_p−1`. The argument in fact
works for every nonzero integer-coefficient polynomial F with F(0)=0
and F(−x)=epsilon F(x), epsilon in {1,−1}; tournaments supply this parity.

Write a=a_p, kappa=kappa_p. Because k! divides C_k,

    b=v_p((kappa−1)!) <= v_p(kappa!) <= a.

Assume b=a. Then p does not divide kappa. The factorial valuation a
occurs on one plateau

    j=mp,mp+1,...,mp+p−1,       m>=1,

with kappa>mp. Every C_j with j<kappa has valuation greater than a by
the definition of kappa; every C_j beyond this plateau also has greater
valuation because j! divides C_j. Therefore the Newton coefficients
of the integer-valued polynomial G=F/p^a, reduced modulo p, can be
nonzero only at

    j=mp+1,...,mp+p−1.                                (7)

At least one is nonzero. Denote these coefficients by c_1,...,c_(p−1),
padding with zeroes beyond the degree if necessary. Since all their
indices exceed p−1, G(r)=0 mod p for r=1,...,p−1. Polynomial parity
forces G(−r)=0 mod p at the same arguments.

Now work with formal power series modulo p. For 1<=r<=p−1,

    (1+z)^(−r)=(1+z)^(p−r)/(1+z^p),

so the coefficient at degree mp+s, 1<=s<=p−1, is

    binom(−r,mp+s)=(-1)^m binom(p−r,s) mod p.          (8)

The equations G(−r)=0 are therefore

    sum_(s=1)^(p−1) c_s binom(p−r,s)=0 mod p.

Take r=p−1,p−2,...,1. This matrix is triangular with diagonal1, forcing
every c_s=0, a contradiction. Hence b<=a−1 and (4)–(6) follow. QED.

The extra excursion symmetry can thicken an existing prime-power
obstruction, but it cannot introduce a new prime and cannot amplify
squarefree raw content. The order-eight (h,g)=(9,27) example attains the
upper bound in (5), so that bound cannot be uniformly reduced.

Parity is essential. The integer polynomial
`F=x^4−6x^3+11x^2+3x=9 binom(x,1)+24 binom(x,4)` has raw content3
and scaled content9. It obeys factorial divisibility and has odd F(1)=9,
but lacks reflection parity, so it violates the strict square bound.
It is a formal hostile, not a tournament polynomial.

## 3. Exact eventual content of repeated ordinal joins

**Incoming attribution.** After this proof was completed, incoming
`origin/main` at `9138a603e2` supplied
[Session295, ARITHMETIC.md](../../session295_construction_atlas_20261113/ARITHMETIC.md).
It independently proves the stabilization law below and strengthens it to
the least onset, using the exact identity
`g(T)=gcd_(j>=1) (j−1)! F_T(j)`. It also gives an effective finite
lower-envelope calculus for arbitrary fixed lists of ordinal factors.
The proof here remains a direct Newton-coordinate derivation; the new
parity bound in Section2 additionally forces its intercept `b_p<=a_p−1`.
The power law is therefore credited as an independent overlap, not as an
exclusive new result of this packet.

Let T^[r] denote the ordinal join of r copies of T, r>=1. Its vertices
and arcs are literal: every vertex of copy i points to every vertex of
copy j when i<j. Define

    R(T)=max(1, max_(p|h(T)) b_p),

with the inner maximum zero when h=1.

**Theorem B.**

    h(T^[r])=h(T)^r                 for every r>=1,
    g(T^[r])=h(T)^r E(T)            for every r>=R(T). (9)

Thus the context-content amplification factor eventually becomes a
constant. The threshold is explicit and at most
`max(1,max_(p|h)(a_p−1))`; it is sufficient, not claimed minimal.

*Proof.* Ordinal multiplication gives F_(T^[r])=F_T^r. Taking the gcd
of its values shows that its raw fixed divisor is h^r exactly: for each
prime, `min_x v_p(F(x)^r)=r min_x v_p(F(x))`.

Fix p|h and abbreviate a,kappa,b. The first minimum Newton coordinate
is also the first minimum evaluation:

    v_p(F(x))>=a+1 for 0<=x<kappa,
    v_p(F(kappa))=a.                                  (10)

Indeed Newton evaluation and finite-difference inversion are triangular
integer transforms with diagonal1 at these consecutive arguments.
For F^r, all Newton coefficients below kappa consequently have valuation
at least r(a+1), while its kappa-th coefficient has valuation exactly ra.
All its Newton coefficients have valuation at least ra by the raw
fixed-divisor identity. Therefore:

* for k<kappa, the scaled coefficient has valuation at least ra+r;
* for k>=kappa, the factorial has valuation at least b, giving at least ra+b;
* at k=kappa, the valuation is exactly ra+b.

Once r>=b the minimum is exactly ra+b. The equality of prime support
from Theorem A excludes every prime not dividing h. Combining the finitely
many prime valuations proves (9). QED.

If h is squarefree, E=1 and R=1, so ordinal powers have exactly raw and
context content h^r at every exponent. Context-primitive T (h=g=1) can
never create nontrivial context content merely by taking ordinal powers.
Different factors can still interact: the inherited H5/H11 complementary
projector example has individually primitive factors but joined content5.

## 4. The asymmetric eight-vertex source and infinite exact families

Use the inherited literal source

    C=(56,117,129,22,100,76,141,59),
    (C_1,...,C_8)=(513,6102,32508,99432,187200,216000,141120,40320).

Its raw content is9, context content27, and automorphism group is trivial.
The first minimum at p=3 is kappa_3=4, so a_3=2,b_3=1,E(C)=3.
The exact polynomial is

    F_C(x)=x^8+34x^6+232x^4+246x^2.                   (11)

Theorem B now gives, with no further tournament census,

    |C^[r]|=8r,       H(C^[r])=513^r,
    h(C^[r])=9^r,     g(C^[r])=3*9^r       (r>=1).    (12)

These tournaments are asymmetric. Each copy of the strong asymmetric
C is a strong component, and the component order is fixed by all arcs
between components; an automorphism cannot permute the ordered copies.
Each restriction to a copy is then trivial.

For a direct valuation check, F_C(0),...,F_C(3) are all divisible by27,
whereas `v_3(F_C(4))=2`. The fourth Newton coefficient of F_C^r therefore
has valuation exactly2r. Its excursion factor3! raises this to2r+1;
all earlier and later weighted coordinates have at least that valuation.

Already at r=2,

    g(C→C)=243 < 27^2=729.                            (13)

Together with Session139's primitive-factor example whose join has
content5, this proves that g is neither submultiplicative nor
supermultiplicative under ordinal join. It is not merely “not exactly
multiplicative.” The source of the discrepancy is the common evaluation
argument in F_T F_U and the factorial weighting of its Newton coordinates.

## 5. What the first order-eight example does and does not explain

Incoming
[Session294, README.md](../../session294_context_27_20261113/README.md)
independently proves the order-seven exclusion and gives the stronger
complete order-eight classification: fourteen exceptional classes with
four distinct slot polynomials. Its labelled order-seven replay supplies
an additional independent route beyond our inherited small-bank control.
Our general prime-dependent lower bound and the formal polynomial hostile
below address the separate question of which arithmetic hypotheses alone
can explain that minimum.

If p divides g/h, Theorem A requires p^2|h. Since h divides C_n=n!,
an odd-prime amplification needs n>=2p. In fact n=2p is impossible.

At n=2p, a=v_p(h)<=2. Amplification therefore requires a=2 and
v_p(g)>=3. The Newton coefficients of G=F/p^2 modulo p vanish through
k=p, because the excursion factorial is a unit there. The only possible
support is k=p+1,...,2p, and its top coefficient is nonzero because F is
monic and C_(2p)=(2p)!. Hence G(p)=0 mod p, but

    (1+z)^(-p)=1/(1+z^p) mod p

shows `G(−p)=C_(2p)/p^2 !=0 mod p`: all the other allowed Newton indices
are not multiples of p. This contradicts parity. Consequently

    p | g/h  implies  n>=2p+1.                       (14)

For p=3 this arithmetic lower bound is7, **not8=3^2−1**. In fact even
monicity, positive even lower monomial coefficients, and a positive
palindromic descent vector do not force the order-eight boundary. The
formal polynomial

    F=x^7+4x^5+4x^3+18x

has ordered-cover Newton vector

    (27,270,2430,9360,17280,15120,5040),

with raw content9 and scaled content27. Its descent vector is

    (27,108,1485,1800,1485,108,27),

which is positive, palindromic, and sums to7!. It is nevertheless not
an actual tournament profile. For example its x^5 coefficient would
give exactly two directed triangles, while its x coefficient would
give nine spanning directed7-cycles by OCF. A spanning cycle makes the
tournament strong; every vertex of a strong tournament lies in a
directed triangle. Two triangles cannot cover seven vertices, so those
two requested statistics are incompatible. This is a physical
realizability obstruction beyond the arithmetic coefficient lattice.
For the used vertex fact, partition the other vertices into the
out-neighbours A and in-neighbours B of v. If no a in A points to a b
in B, the nonempty set A cannot reach v or B, contrary to strongness.
An arc a->b gives the triangle v->a->b->v.

The independent finite control uses the inherited complete bank of321
representatives up to reversal at orders3 through7, accounting for530
ordinary isomorphism classes; orders1 and2 are transitive. Reversal
preserves every path-cover count. Every one has h=g. Together with C,
this establishes first amplification order8 **within that complete
finite universe**, not from a universal “p^2−1 threshold” theorem.

## Reproduction and limits

Run `python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/content/check.py`.
The receipt names the exact bank revision and hash, records the whole
source profile, the formal hostile, and a bounded check of ordinal
powers. The finite powers are controls; Theorems A and B supply their
all-order quantifiers. The script does not enumerate large tournaments.

The context ideal is an ideal of **signed integer combinations** of
physical responses. Neither h=1 nor g=1 asserts that any prescribed
positive count is a physical context response. These results constrain
and organize constructors; they do not close the Hamiltonian spectrum.
