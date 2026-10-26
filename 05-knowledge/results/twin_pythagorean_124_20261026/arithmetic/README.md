# Twin primes as a prime inradius–exradius fibre

**PROVED elementary statements; FINITE-EXACT controls; independent audit pending.**
These are exact translations and boundary classifications, with no historical
novelty claim and no proof of twin-prime infinitude or abc. The central result
is a bijection: primitive right triangles whose inradius and exradius opposite
the even leg are both prime correspond exactly to twin-prime pairs.

## 1. Inheritance and scope

The source checkout is `52fde1af0df4` (26 October 2026). The bounded startup,
current routes, protocol, guardrails, primary-source ledger, and targeted
messages were read. All this lane's output was written to isolated staging.
Source hashes are in [SOURCES.json](SOURCES.json).

Closest mechanisms:

* [THM-3333 / Gaussian-square light cone](../../../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md),
  equations (12), (13), (18): the Euclid/Gaussian lift, radius, and parity content.
* [Perfect-core twin edges, §§1–3](../../perfect_core_twin_edges_20261019/PERFECT_CORE_TWIN_EDGES.md):
  the guarded map `6p -> 6(p+2)` for primes `p>=5`.
* [S69 aliquot report, §§2 and 5](../../session69_prime_ports_power_aliquot_20261026/session69-aliquot/REPORT.md)
  and [late S69 integration](../../aliquot_power_124_20261026/INCOMING_S69.md):
  `s(2p)=s((p+2)^2)=p+3`, inverse bounds, and distinct abc maps.
* [Unique-six theorem](../../aliquot_power_124_20261026/aliquot/README.md)
  and [weighted monoid carrier](../../session66_root27_aliquot_graph_interfaces_20261026/session66-root/DIVISOR_F1_BRIDGE.md).
* [THM-2422 / operation fibres](../../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md):
  factor-pair diagonals and the exceptional twin centre 4. Its additive-parent
  census is finite evidence only; no all-height ancestry claim is imported.

Canonical hostiles are the prime-square factor diagonal, the failed prime guard
at `6*3`, and unweighted divisor chains `4` versus `9`. The corrected near miss
is treating every twin centre as divisible by six: centre 4 is exceptional
(MISTAKE-268). The least-used sidecar here is **which exradius**, with the even
leg distinguished by parity. The board is radii / Euclid factors / aliquot
fibres / abc height and radical / Gaussian norm / weighted divisor support.
The existing “search the statement” and “controlled forgetting” methods apply.

## 2. Radius coordinates and the exact twin-prime theorem

An oriented primitive Pythagorean triple means positive integers `(A,B,C)`
with `A` the odd leg, `B` the even leg, `A^2+B^2=C^2`, and `gcd(A,B,C)=1`.
The standard Euclid parametrization has a unique pair

    m>n>=1, gcd(m,n)=1, m-n odd,
    A=m^2-n^2, B=2mn, C=m^2+n^2.                 (1)

This classical parametrization is recorded in
[Tripathi, introduction, p.331](https://www.fq.math.ca/Papers1/46_47-4/Tripathi.pdf).
For completeness, `(C+A)/2` and `(C-A)/2` are coprime positive integers
whose product is `(B/2)^2`; each is a square, giving `m,n`. The parity and
coprimality conditions follow from primitivity, and give uniqueness.

Let `h=(A+B+C)/2=m(m+n)` and let `Delta=AB/2`. The inradius `r=Delta/h`
and exradius `rho=Delta/(h-B)` opposite the **even leg** satisfy

    r=n(m-n),       rho=n(m+n),
    rho-r=2n^2=C-A, rho+r=2mn=B.                (2)

Thus the marked radius pair determines the triangle:

    n^2=(rho-r)/2, m=(rho+r)/(2n),
    A=r*rho/n^2,  B=rho+r,  C=A+rho-r.        (3)

Formula (3) is an inverse on the actual primitive-triangle image, not a claim
that arbitrary positive integer pairs are valid radius data.

**Theorem A (prime radii).** The following predicates on such a triple are
equivalent:

1. `r` and `rho` are both prime.
2. `n=1` and `m-1,m+1` are prime.
3. `C-A=2` and `Omega_arith(A)=2`, where multiplicity is counted.
4. For a unique odd twin-prime pair `(p,p+2)`,

       (A,B,C)=(p(p+2), 2(p+1), (p+1)^2+1).   (4)

**Proof.** If `rho=n(m+n)` is prime, `m+n>1` forces `n=1`; then (2)
gives `r=m-1`, `rho=m+1`. Opposite parity makes `m` even; primality of
`m-1` excludes `m=2`. Hence `p=m-1>=3` is odd. Conversely these conditions
give primitive (4) and its claimed prime radii. If `C-A=2`, (1) forces
`n=1`. For `m=2`, `A=3` has multiplicity one; otherwise the two factors
`m-1,m+1` exceed one. Their total multiplicity is two exactly when both
are prime. This proves all directions and uniqueness. QED.

The first triples in parity orientation are `(15,8,17)`, `(35,12,37)`, and
`(143,24,145)`. They have prime-radius pairs `(3,5)`, `(5,7)`, `(11,13)`.
The root `(3,4,5)` has `(r,rho)=(1,3)` and is excluded. For every pair after
`(3,5)`, the centre `m` is divisible by six and `B` is divisible by twelve.
The first centre `m=4` and leg `B=8` must stay separate.

### The full prime-inradius fibre has a second branch

The cardinality is a classical specialization of
[Robbins, Theorem2, p.369](https://www.fq.math.ca/Papers1/44-4/quartRobbins04_2006.pdf),
which counts all fixed-inradius fibres by distributing coprime prime powers.
We reproduce the prime case to retain its two explicit branches and parity
exception; the count itself is not a new result.

**Theorem B.** For an odd prime `p`, exactly two oriented primitive triples
have inradius `p`. Their Euclid parameters are

    (m,n)=(p+1,1) or (p+1,p),

and their side triples are respectively

    (p(p+2), 2(p+1), p^2+2p+2),
    (2p+1, 2p(p+1), 2p^2+2p+1).              (5)

For `p=2` the only one is `(5,12,13)`, with `(m,n)=(3,2)`.

**Proof.** The positive factorization `p=n(m-n)` permits only `(n,m-n)`
equal to `(1,p)` or `(p,1)`. Both have `m=p+1`; for odd p both obey the
coprimality and opposite-parity conditions and are distinct. For p2 the
first pair `(3,1)` has odd/odd content two and is not primitive, whereas
`(3,2)` is valid. These exhaust the factorizations. QED.

Only the first branch can have prime `rho`: the second has
`rho=p(2p+1)`. For example `(7,24,25)` has prime inradius3 but rho21;
`(63,16,65)` lies on the first branch with prime inradius7 but rho9.
Thus neither “prime inradius” nor “gap C-A=2” alone recognizes twins.

## 3. What transports to aliquot fibres and what does not

For prime radii `(r,rho)=(p,p+2)`, the inherited S69 formulas become

    s(2r)=1+2+r=rho+1=s(rho^2).               (6)

Conversely, for odd primes p,q, equality `s(2p)=s(q^2)` is exactly `q=p+2`.
Consequently ordered prime-radii triangles, twin pairs, and these **marked
two-witness aliquot collisions** are bijective. The target `M=p+3` determines
p, but the two marked witnesses need not exhaust the full fibre `s^-1(M)`.
The first failure of exhaustion among twin sources p<=97 is

    p=71, M=74,       s^-1(74)={70,142,5329}.  (7)

Here 142 and5329 are the twin witnesses; `s(70)=74` is extra. Equation (7)
is FINITE-EXACT with a proof of the inverse cutoff: for M>1, every source is
composite. If t is its least prime divisor and q=N/t>=t, the proper divisors
1 andq give `q<=M-1`, hence `N<=q^2<=(M-1)^2`. A separate additive sieve and
factor-product check certify all fibres M2..100, including the endpoint5329.

The witness `2p` is fixed exactly when `2p=p+3`, namely p3. Moreover the
inherited **global** unique-six theorem says every perfect target receiving
a prime-square predecessor is6. Its mechanism is retained: if `s(q^2)=N`
is perfect then `N=q+1`. The case q2 fails; otherwise N is even perfect,
`N=2^(a-1)(2^a-1)` with prime a. For odd a>=3, `N-1` is a multiple of3
greater than3, contradicting primality of q. The remaining a2 gives N6,q5.
This uses Euclid–Euler, cited and proved in the inherited note; it assumes
nothing about odd perfect numbers. The pair `(3,5)` and triangle `(15,8,17)`
are exactly the fixed-target boundary of (6).

The other inherited map `s(6p)=6(p+2)` is valid only for prime `p>=5`.
It gives an actual aliquot edge, while (6) gives two incoming arrows.
At p3 the first map is invalid: `s(18)=21`, not30. No equality of these
two diagrams, continuation of a prime guard, or conjugacy of aliquot orbits
with a Pythagorean tree is asserted.

## 4. Two abc maps, and their exact radical cost

Write `rad(t)` for the product of distinct primes dividing t. For twins
`p,p+2` let `m=p+1`. The primitive twin-centre abc triple is

    1 + p(p+2) = m^2.                          (8)

The three integers p,p+2,m are pairwise coprime, so its radical is

    R_0=p(p+2)rad(m)=(m^2-1)rad(m)>m^2.        (9)

Indeed m is even and at least4, so `R_0>=2(m^2-1)>m^2`. Thus this subfamily
already satisfies the stronger elementary `height<radical` bound. Applying
abc to (8) supplies no exclusion or prime-counting lower bound.

The Euclid square map sends (8) to another primitive abc triple,

    (m^2-1)^2 + 4m^2 = (m^2+1)^2.             (10)

Since the sides A,B,C of (4) are pairwise coprime, its radical is

    R_1=rad(ABC)=R_0 * rad(m^2+1).             (11)

Its height is `(m^2+1)^2`, compared to `m^2` in (8). Equations (9)–(11)
give the exact support and height change; equality of the two abc qualities
is not preserved. Assuming abc would give the correctly scoped bound

    rad(m^2+1) >= (m^2+1)^[2/(1+epsilon)]
                  / (K_epsilon^[1/(1+epsilon)] R_0),        (12)

with an unknown uniform `K_epsilon`; no effective cutoff is extracted.
The ordinary source map `(N,s(N),sigma(N))/gcd(N,s(N))` is yet another map:
the equal-target sources6 and25 give `(1,1,2)` and `(25,6,31)`, whereas
the twin-centre map at p3 is `(1,15,16)`. Map identity is an indispensable
sidecar even when some other constructions share outputs.

For provenance, [Granville–Tucker, pp.1226–1227](https://dms.umontreal.ca/~andrew/PDF/abc.pdf)
give the uniform-epsilon conjecture and explicitly exhibit the polynomial
Pythagorean identity as sharp for polynomial abc. That polynomial theorem
does not turn an integer specialization into a theorem on twin primes.

## 5. Gaussian norm loses the radius predicate

The inherited Gaussian map is

    (m+in)^2=A+iB,       N(m+in)=C.             (13)

On the twin branch n1, squaring `m+i` retains all data in A+iB; the norm
alone does not. At norm145, `(m,n)=(12,1)` gives the twin triangle
`(143,24,145)` with radii11,13, whereas `(9,8)` gives `(17,144,145)`
with radii8,136. Thus even a norm known to support a twin triangle does not
make every primitive lift a twin triangle. The restoration sidecar is the
chosen positive Euclid lift or its marked radius pair.

Norm primality is also independent of the displayed twin test: twins11,13
give composite norm145; centre10 gives prime norm101 but factors9,11,
which are not both prime. These two witnesses refute either implication.

## 6. The F1-compatible carrier makes the lost arithmetic explicit

Let the free multiplicative monoid on primes have monoid ring
`Z[X_l : l prime]`. For an integer N, its proper-divisor enumerator `D_N`
has one monomial per proper divisor. The inherited maps are arithmetic
evaluation `E(X_l)=l` and augmentation `epsilon(X_l)=1`:

    E(D_N)=s(N),          epsilon(D_N)=tau(N)-1.              (14)

For a twin pair p,q=p+2, the two incoming sources in (6) have

    D_(2p)=1+X_2+X_p,       D_(q^2)=1+X_q.

Their difference is the nonzero polynomial `X_2+X_p-X_q`. Evaluation E
annihilates it exactly when q=p+2, but augmentation gives1. Consequently

    E(D_(2p))=E(D_(q^2)),   epsilon(D_(2p))=3 != 2=epsilon(D_(q^2)). (15)

The twin collision is not a collision of unweighted divisor shapes; it
is created by the additive arithmetic evaluation. Conversely all prime
squares have augmentation2 but `s(9)=4 != 6=s(25)`. Thus neither equality
relation refines the other. Prime-labelled weights restore the arithmetic.
In particular augmentation cannot descend through E to a map on its
arithmetic image: a relation killed by E need not be killed by epsilon.

This is compatible with [Deitmar, §1 and Theorem1.1](https://arxiv.org/html/math/0404185v7):
in that model F1 objects are commutative monoids and base extension is the
monoid ring. It is not a new field structure on aliquot fibres. No F1
construction in this note determines additive prime gaps after the labels
and evaluation map have been forgotten.

## 7. Transfer ledger and exact reproduction

| Source -> target | Preserved predicate | Lost data / restoration | Cheapest hostile |
|---|---|---|---|
| Twin pair -> primitive triangle (4) | Both marked radii prime; bijection | None in parity orientation | Root `(3,4,5)` has radius1 |
| Triangle -> inradius only | That radius is prime | Which branch in (5); retain rho | r3 has `(15,8,17)` and `(7,24,25)` |
| Twin pair -> marked collision (6) | Same aliquot target; bijection of marked witnesses | Other preimages; retain full fibre | M74 has extra70 |
| Twin pair -> centre abc (8) | Primitive identity; invertible on the marked image | No loss there; dropping the twin domain permits composite neighbours | m8 gives1+63=64, not twins |
| Centre abc -> square abc (10) | Exact polynomial identity | Height and new norm support change | R1/R0=rad(m²+1) |
| Euclid lift -> Gaussian norm | Norm/hypotenuse | Lift and prime radii | Two lifts of145 |
| Divisor monoid -> augmentation | Number of proper divisors | Prime weights and collision | Twin values agree; counts3 and2 |

Run from this directory:

    python3 -I -B check_arithmetic.py > normal.json
    python3 -I -B -O check_arithmetic.py > optimized.json

The checker uses explicit exceptions, not Python assertions. Universes are
Euclid pairs m<=200, direct side enumeration C<=1000, twin lower primes
p<=10000, prime-inradius reconstruction p<=199, divisor polynomials N<=1000,
and complete aliquot targets2..100 with source bound9801. Side enumeration
is independent of the Euclid producer; a divisor sieve is independent of
factor-product sigma. Polynomial expansion controls (10) symbolically.
The receipts and manifest state actual counts and hashes. These finite
controls support the proofs and do not establish infinitude, abc, the open
twin-centre additive ancestry, or an aliquot/Pythagorean dynamical conjugacy.
