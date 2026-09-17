# The Collatz inverse-fibre braid and its missing global coordinate

**Status: PROVED elementary identities and obstructions; FINITE-EXACT controls.**
Collatz convergence, positive-cycle uniqueness, and completeness of the three
known negative cycles remain **OPEN** here. No novelty claim is made for the
standard inverse or exponent-word formulas.

## Inheritance and research posture

Anchor: the user's three Collatz rows and summand/doubling decomposition.
Niche: divisor supports and local prime sandwiches. Wildcard: normalized
Pythagorean geometry and quadratic period three. The six live concepts are
odd-core/height, paired summand incidence, prime-exponent support, oriented
neighbor sieve, circle/scale, and marked quadratic cycle/trace.

Closest proved mechanism: [THM-2422, operation fibres and twin-center ancestry](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md),
especially its swap-fixed diagonal and synchronous dyadic closure. The
historical [rapidity identity](../../07-reflections/collatz-rapidity-defect.md)
is independently recovered below, without inheriting its conjectural uniform
defect bound. The canonical hostile is a long exponent-one growth prefix.
The corrected near miss is the historical promotion of sampled defect bounds
to a theorem. The least-used relevant sidecar is the **ordered halving word
and its affine carry**, not just its mean or residue support.

Research cards used: "Separate unbounded local support from a height-bounded
modular cover", "Search the statement before the method", and "Type every
analogy and every implication" in [META-PATTERNS](../../00-navigation/META-PATTERNS.md).
The startup Collatz match in THM-3480 concerns a Collatz-Wielandt matrix
certificate, not the integer iteration; it supplies no dependency here.

## 1. Exact typing of the three rows

For positive odd n, set

    F(n)=(3n+1)/2=2^h u,       u odd, h>=0,
    T(n)=u,                   k=h+1=v_2(3n+1).

F is the single-halving map; T is the accelerated odd-to-odd map. For j>=0,

| Input | Single-halving image |
|---|---|
| 6j+1 | 9j+2 |
| 6j+3 | 9j+5 |
| 6j+5 | 9j+8 |

The powers of two in these three image rows have exponents respectively
1,5,3 modulo6, since the order of2 modulo9 is6. This characterizes **which
powers of two lie in a row**; it does not make every row member a power of two.

Every image odd core satisfies 3∤u. Conversely every positive odd u with
3∤u occurs, infinitely often, in **each** of the three source rows. Multiples
of3 can be starting nodes but can never be T-images. They are leaves in the
odd reverse graph; their own forward trajectories remain consequential.

For u≡1 mod6 the admissible h are odd; the lowest target is2u and its source
(4u-1)/3 is1 mod8. For u≡5 mod6 the admissible h are even; the lowest target
is u and its source (2u-1)/3 is3 mod4. These are exactly the user's entry rules.

## 2. The inverse-fibre theorem and the triadic odometer

**PROVED.** Fix positive odd u with3∤u. Put h0=1 for u≡1 mod3 and h0=0
for u≡2 mod3. All positive odd sources with T(n)=u are precisely

    h_j=h0+2j,
    n_j=(2^(h0+1) 4^j u-1)/3,              j>=0.                 (B1)

Indeed integrality is equivalent to2^(h+1)u≡1 mod3, fixing h's parity;
the numerator divided by3 is odd. Therefore

    n_(j+1)=R(n_j),       R(n)=4n+1,
    F(R(n))=4F(n),        T(R(n))=T(n).                         (B2)

Modulo6, R cycles1→5→3→1. Returning to the same source row takes three
applications: n→64n+21 and h→h+6. This is the exact two-unit height shift
and three-way interleaving; the repeated small odd cores are inverse fibres.

**Stronger PROVED refinement.** R is one cycle on all3^s odd residue
classes modulo2·3^s, for every s>=1. More precisely, for every integer n
and t>=1,

    R^t(n)-n=(4^t-1)(3n+1)/3,
    v_3(R^t(n)-n)=v_3(t).                                      (B3)

Proof: 3n+1 is a3-adic unit, and v_3(4^t-1)=1+v_3(t).
The latter follows by binomial expansion for multiplication of t by3,
and by reducing the geometric sum modulo3 when3∤t. R preserves oddness;
its exact period modulo2·3^s is3^s, exhausting all odd classes.
Thus j→n_j is an isometry in the3-adic metric and, after completion, a
copy of addition by1 on Z_3. This is an actual small-rule/all-scale law.

Its limit is crucial: the whole residue cycle above can lie over a **single
fixed target u**. Full local support says nothing about the later orbit of u.
[The summand companion](arithmetic_braids_20260917_summand.md) supplies a
stronger hostile:5n+1 has the analogous full5-adic inverse braid and a
nontrivial positive cycle13→33→83→13.

## 3. An exact paired-summand lift

For odd n>1 the shortcut edge n→F(n) has the distinct positive co-parent
(n+1)/2. The labelled summand triple is

    (a,b,z)=(n,(n+1)/2,(3n+1)/2),   a+b=z, a=2b-1.

R lifts to

    (a,b,z)→(4a+1,4b-1,4z).                                    (B4)

Opposite unit corrections in the two parents cancel, while the target
travels two steps up its doubling tower. Both summand identity and parent
relation survive. At n=1 the parent pair is(1,1), the excluded diagonal;
R sends it to the admissible pair(5,3). Thus the trivial cycle is a genuine
boundary exception, not a proof that distinctness survives inverse transport.

Connection contract: source = labelled Collatz summand triples; target =
odd-core-labelled dyadic towers; map = (B4); preserved = sum, parent affine
relation, target odd core; lost by quotient = height and both parent labels;
restoring sidecar = h and a (which reconstruct b); decisive test = n=1 versus5.

## 4. Corrected row and why four early cores do not give Goldbach

The first35 compressed pairs for the3 mod6 input row are

    (5,0),(7,1),(23,0),(1,5),(41,0),(25,1),(59,0),
    (17,2),(77,0),(43,1),(95,0),(13,3),(113,0),(61,1),
    (131,0),(35,2),(149,0),(79,1),(167,0),(11,4),(185,0),
    (97,1),(203,0),(53,2),(221,0),(115,1),(239,0),
    (31,3),(257,0),(133,1),(275,0),(71,2),(293,0),
    (19,4),(311,0).

The next pair is(5,6), since9·35+5=320=2^6·5. In particular the first
pure tower is(1,5), not(2,5);77 is odd,95 was omitted, and the middle entry
is79, not39. The complete three rows are in the JSON artifact.

The first four cores5,7,23,1 are explained by odd-part extraction of
5,14,23,32. This constructs no prime-sum representation:1 is not prime,
and later cores include25,35,77,95. A Goldbach connection would need a
map preserving primality and a specified target sum; none follows here.

## 5. Exact valuation laws, and the obstruction to a residue-only proof

Since9 is odd, j→9j+b permutes residues modulo2^H. On every full block
of2^H consecutive row indices, exactly2^(H-h-1) have valuation h<H,
and one has valuation at least H. Hence each row's natural density of
height h is2^(-h-1). Its T-output belongs to1 mod6 with density1/3 and
to5 mod6 with density2/3. These are counting densities, not deterministic
independence assertions along a specified orbit.

More generally fix a positive halving word(k1,...,kL), put K_i=sum_(j<=i)k_j,
K_0=0, and

    B=sum_(i=0)^(L-1) 3^(L-1-i) 2^K_i.

Its exact source condition is

    3^L n+B ≡ 2^K_L (mod 2^(K_L+1)).                           (B5)

This specifies one odd residue modulo2^(K_L+1). Sufficiency follows
backwards: the last numerator is2^kL times an odd value, and reduction
successively gives each prior exact valuation. Necessity is the iterated
identity2^K_L T^L(n)=3^L n+B. CRT then places infinitely many positive
sources for this same word in **each** row modulo6, with relative density
2^(-K_L) within that row. The residue condition is elementary; compare the
broader2-adic parity conjugacy of [Bernstein--Lagarias](https://websites.umich.edu/~lagarias/doc/bernstein.pdf).

For a concrete hostile, n0=2^(L+1)-1 has L consecutive k=1 steps and

    T^j(n0)=3^j 2^(L+1-j)-1  (0<=j<=L),
    T^L(n0)=2·3^L-1.

The expansion factor is unbounded with L. CRT gives such long prefixes
in any chosen source row. No fixed block length guarantees descent for
every integer; no residue-only positive weight w_(n mod6) can make
w_(T(n) mod6)T(n)<w_(n mod6)n at every n, because arbitrarily large
exponent-one sources have both source and target5 mod6.

**Stronger finite-congruence obstruction (independent audit).** For any
modulus M>=1, block length L>=1 and integer q>=1, choose

    n0=2^(L+1) Mq-1.

Then T^j(n0)=3^j 2^(L+1-j) Mq-1 for0<=j<=L. Every step has k=1,
every node has residue -1 modulo M, and T^L(n0)>n0. Hence for **any**
positive periodic weight w modulo M, the proposed Lyapunov quantity
V(n)=w(n mod M)n increases over this block. The same holds for
V(n)=log n+b(n mod M), with arbitrary real b. No fixed finite collection
of congruence coordinates and no fixed block length can make either
form decrease everywhere. This does not rule out adaptive stopping times,
unbounded moduli, or different nonperiodic Lyapunov functions.

This directly challenges a transfer from the sandwich CRT product: finite
prime-local data are exact on their finite universe, but adding any fixed
set of such data cannot fix this deterministic descent obstruction.

## 6. The global quantity that is still owed

For any finite positive odd orbit segment a0,...,aL,

    log(aL/a0)=L log3-K log2+D,
    D=sum_(i=0)^(L-1) log(1+1/(3a_i)).                         (B6)

This identity is unconditional for every finite segment. D is positive
when L>0. Descent is exactly K log2>L log3+D. For a segment terminating
at1 it gives the inherited rapidity formula; at n=1 with zero steps D=0.
It neither proves termination nor gives a uniform bound on D.

The historical reflection incorrectly treated D<0.2257 as proved and
identified a maximum-orbit envelope with Tao's theorem. [Tao, Theorem1.3](https://arxiv.org/pdf/1909.03562)
instead concerns orbit **minima**, for almost all starts in logarithmic
density and every diverging comparison function. No upper bound for the
whole orbit maximum is imported. A geometric iid model may motivate a
drift, but applying it to all deterministic trajectories is the missing step.

A cycle with ordered word k has the exact integrality gate

    n=B/(2^K-3^L).                                             (B7)

Both B and K matter; permuting the word generally changes B. The script
checks all106,761 ordered words with1<=L<=8 and L<=K<=18, reconstructs
actual intermediate integers/valuations, and finds the cycles

    (1), (-1), (-5,-7), (-17,-25,-37,-55,-41,-61,-91).

This is a FINITE-EXACT bounded-word statement. Negation conjugates signed
3n+1 to positive3n-1 and swaps1/5 mod6 while fixing3 mod6. Existence of
three negative cycles is exact; completeness is not proved.

## Evidence and next obligations

Reproduce with

    python3 04-computation/experiments/arithmetic_braids_20260917_collatz.py

[Source](../../04-computation/experiments/arithmetic_braids_20260917_collatz.py)
and [output](arithmetic_braids_20260917_collatz.json) freeze all universes,
filters, controls, and the source hash. The proofs use modular inversion,
binomial valuation, and exact identities; finite controls additionally use
forward division, inverse generation, CRT realization, and full residue sets.

The strongest next target is an operation-stable restriction on **ordered
halving words plus carry B and integer height** that survives passage to
the next odd core and forces eventual inequality(B6). The inverse braid
alone cannot supply this. No transfer to LRC(14) is claimed: its analogous
local-support/height warning transports as a research method, not a theorem.
