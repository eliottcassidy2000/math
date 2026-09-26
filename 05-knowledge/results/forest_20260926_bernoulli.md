# Bernoulli boundaries, constructible polygons, and a lossless carry graph

**Status:** CITED classical identities; PROVED elementary deductions and
obstructions; FINITE-EXACT controls. No Collatz convergence or Fermat-prime
finiteness claim. The carry-budget identity in section 6 was independently
proposed by the root lane and checked here. Classical consequences below
are not presented as new number-theory theorems.

## 1. Inheritance and the corrected seed

Closest proved mechanism: [the Bernoulli boundary note](bernoulli_boundary_20260925.md)
already makes the periodic B1 jump an exact divisibility detector, retaining
the singleton lost by a density calculation. The hostile is an arbitrarily
long expanding positive shadow of a negative rational Collatz cycle.
[THM-4507 / finite polynomial-valuation obstruction](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md)
shows why an arbitrary function of finitely many polynomial valuations
cannot repair an every-step logarithmic rank. The corrected near miss is
the idea that five initial examples establish a finite catalogue. The
least-used sidecar here is the **ordered carry path through every binary
position**, rather than its total or its finite prime support.

Anchor / niche / wildcard: exact source-preserving carry depth; Bernoulli
denominators and polygon construction; a three-state transducer and its
discounted orbit budget. The live concepts are prime support, dyadic depth,
endpoint convention, affine source residue, carry path, and excursion cost.

The first-five-primes statement concerns the **Fermat numbers**

    F_j = 2^(2^j)+1: 3,5,17,257,65537,...,

not Bernoulli numbers. In the usual convention the latter begin
`1,-1/2,1/6,0,-1/30,0,1/42`. In fact no Bernoulli number is a positive
integer prime: every positive even index has denominator divisible by six,
the odd indices after one vanish, and B0=1. Five Fermat primes are known;
it is not proved that there are no others. The current
[PrimePages Fermat-divisor catalogue](https://t5k.org/top20/page.php?sort=FermatDivisor)
was checked on 2026-09-26. The first failure of primality is exact:

    F_5 = 4294967297 = 641 * 6700417.

The repo's [Proth unification, THM-1355](../../01-canon/theorems/THM-1355-proth-unification-observer-hypotenuse-two-axes.md)
already places Fermat numbers inside `u*2^r+1`. Its
[arithmetic-seam audit](arithmetic_seams_20260921_primes.md) already separates
a new prime divisor from primality of the whole Fermat number. Neither
classification is changed here.

## 2. An exact Bernoulli-to-polygon bridge

Let `D_r=denominator(B_(2^r))`, in lowest terms, for `r>=1`.
Von Staudt--Clausen gives

    D_r = 2 * product_(j>=0, 2^j<=r, F_j prime) F_j.       (1)

Proof: a denominator prime satisfies `p-1 | 2^r`. Apart from p=2 it is
`p=2^d+1` with `d<=r`. If d has an odd factor bigger than one, the
difference-of-powers factorization makes p composite; hence `d=2^j`.
Conversely every prime of that form satisfies the divisor criterion.
The imported denominator theorem is
[DLMF 24.10.1](https://dlmf.nist.gov/24.10.E1).

Consequently `D_r/D_(r-1)` is a new Fermat prime precisely at an index
`r=2^j` for which F_j is prime (r>=2). The initial plateaux are

| r | D_r |
|---|---:|
| 1 | 6 |
| 2..3 | 30 |
| 4..7 | 510 |
| 8..15 | 131070 |
| 16..32 | 8589934590 |

The last displayed range uses only primality of F0..F4 and compositeness
of F5; it is not an assertion that all later plateaux persist. More
generally, eventual stabilization of D_r is equivalent to finiteness of
the Fermat primes. Stabilization at D16 is equivalent to absence of any
further Fermat prime.

By the classical regular-polygon criterion, for every integer n>=3,

    the regular n-gon is constructible
      iff oddpart(n) divides D_r/2 for some r.             (2)

The criterion is that n be a power of two times distinct Fermat primes;
see [Hochster's notes, chapter 12, pp. 95--96](https://sites.lsa.umich.edu/hochster/wp-content/uploads/sites/1337/2024/11/fib20B.pdf).
Equation (2) is an immediate combination with (1), not a novel
classification. The divisor condition retains the required squarefreeness.

**Loss and hostile.** D_r records prime support, not depth of an integer's
binary expansion or location of an orbit. Constructible sizes are not
invariant under the odd Collatz map: `17 -> 13`, whereas the regular
13-gon is not constructible. Even if only five Fermat primes exist,
regular `2^L`-gons and their increasing quadratic construction depth exist
for every L. A finite prime catalogue does not bound dyadic resolution.

The difference-table lane supplies one exact joint bridge: `x_n=2^n+1`
has `Delta^r x_i=2^i` for every r>=1, so every noninitial left edge is
one, even after replacing differences by absolute differences. The Fermat
numbers are the sparse subsequence `x_(2^j)`. This constant-edge property
therefore coexists with the composite F5. It cannot by itself express
the arithmetic restriction in Gilbreath's conjecture that the entire
initial row consists of consecutive primes. The formula follows directly
by taking a first difference and then noting that `Delta 2^i=2^i`.

## 3. The same cyclotomic factors give an exact source selector

For L>=1, N=2^L, an integer d, and a primitive Nth root zeta,

    1_(N|d) = (1/N) sum_(j=0)^(N-1) zeta^(jd)
             = (1/N) product_(j=0)^(L-1) (1+zeta^(2^j d)). (3)

This follows from the polynomial identity

    product_(j=0)^(L-1) (1+z^(2^j)) = sum_(j=0)^(2^L-1) z^j.

At z=2 those factors are the Fermat numbers, whether prime or composite.
At z=zeta^d they are an exact dyadic congruence detector. All required
regular N-gons are constructible independently of any odd Fermat prime.

With the endpoint convention `psi(x)={x}-1/2`, the identical detector is

    1/N + psi((d-1)/N)-psi(d/N) = 1_(N|d).                (4)

Here `psi(integer)=-1/2`; replacing it by the Fourier midpoint loses the
exact detector. This is the previous Bernoulli-boundary mechanism.

For a length-L shortcut Collatz word w, let a be its number of odd
steps and let C be its affine carry. Then

    T_w(n) = (3^a n+C)/2^L,

and the word is legal at the integer n iff `2^L | 3^a n+C`.
Put `d=3^a n+C` in (3) or (4). This transfers the polygon/Bernoulli
formula to the **exact fixed source**, not a distribution of sources.
Word construction uses `C_(j+1)=3^e C_j+e*2^j`.

What is preserved: the complete length-L parity obligation. What is
destroyed: n and n+2^L still share the selector. The needed sidecar is
the full coherent residue path, or an actual integer source plus its
quotient. Coherent residues define a 2-adic integer; they define a
nonnegative ordinary integer iff the least nonnegative residues
eventually stabilize. These are different quantifiers. Formula (3)
does not turn arbitrarily deep positive shadows into one positive
infinite orbit.

## 4. A three-state graph carrying the whole integer

For n>=0 and s>=1 define the scale-s carry

    c_s(n) = floor((3*(n mod 2^s)+1)/2^s),
    c_0(n) = 1.                                         (5)

Each c_s is 0, 1, or 2, and the path eventually stays at zero. If b_s is
binary digit s of n, direct division gives

    c_(s+1) = floor((c_s+3*b_s)/2).                       (6)

The graph is small, but its **ordered path has unbounded length**:

| present carry | next carry for bit 0 | next carry for bit 1 |
|---:|---:|---:|
| 0 | 0 | 1 |
| 1 | 0 | 2 |
| 2 | 1 | 2 |

The outgoing endpoints in each row differ, so the entire carry path
recovers each b_s uniquely. Conversely every path starting at c0=1
that eventually stays at zero encodes a unique nonnegative integer.
An arbitrary infinite path is a 2-adic source; eventual zero is precisely
the extra ordinary-integer condition. For an odd source the first edge
is `1 -> 2`.

The Bernoulli expression for every scale is also exact:

    c_s(n) = 3 psi(n/2^s)-psi((3n+1)/2^s)+1+2^(-s).       (7)

There is no convergence issue at a fixed integer: (5) vanishes for all
2^s>3n+1. Thus the root forest can carry the finite vector (5), its
source, and its location inside the affine block. This is a representation,
not a descent theorem; its losslessness means it is an exact recoding of
the integer rather than a finite-state shortcut to the Collatz problem.

## 5. A coordinate unbounded inside every fixed valuation fibre

Write S(n) for the binary digit sum and define the total carry

    K(n)=sum_(s>=1) c_s(n) = 1+3S(n)-S(3n+1).             (8)

Indeed `c_s(n)=floor((3n+1)/2^s)-3floor(n/2^s)`, and
`sum_(s>=1) floor(m/2^s)=m-S(m)` for every m>=0.

**PROVED fibre escape.** Fix finitely many primes and integer polynomials,
and an odd a>0 at which none of the polynomials vanishes. Choose an even
Q>a divisible, at each selected prime p, by a power strictly higher than
all selected `v_p(P_i(a))`. Put b=Q-a and

    n_r = Q*2^r-b.                                       (9)

For every r, `n_r=a mod Q`, so every selected polynomial valuation is
exactly the one at a. Once `2^r>3b-1`, binary complementation gives

    S(n_r) = r + S(Q-1)-S(b-1),
    K(n_r) = 2r + 1 + 3S(Q-1)-3S(b-1)
                      -S(3Q-1)+S(3b-2).                 (10)

To see the first line, write
`n_r=(Q-1)*2^r+(2^r-b)`, whose two summands occupy disjoint bits, and use
`S(2^r-b)=r-S(b-1)`. Apply the same calculation to
`3n_r+1=3Q*2^r-(3b-1)` for the second line.

This supplies the requested unbounded information inside a fixed finite
feature fibre, with no probability assumption and no convergence premise.
For the six polynomials and six primes in the script, a=27 yields

    Q=24913785600, n_r=Q*2^r-(Q-27),
    S(n_r)=r+4, K(n_r)=2r+10, r>=37,

while all 36 selected valuations remain fixed. These n_r are a family of
inputs, not successive points on one orbit.

**Immediate hostile to collapsing depth.** The growth edge `51 -> 77`
has identical `(S,K)=(4,9)` at its ends. There are arbitrarily high
instances:

    n_L=15*2^L+51,
    U(n_L)=45*2^(L-1)+77 > n_L,       L>=9,
    (S(n_L),K(n_L))=(S(U(n_L)),K(U(n_L)))=(8,17).          (11)

The bit blocks are separated, and `S(15)=S(45)=S(135)=4`.
Thus no `c log n+f(S(n),K(n))`, c>0 and arbitrary f, is nonincreasing
at every sufficiently large U edge. This statement does not by itself
exclude adding an arbitrary bounded correction: its one-edge logarithmic
gain is bounded. A stronger multi-edge argument would be needed for that.
The pair of totals loses the positions and order of the carries. The
three-state path retains them.

## 6. An unconditional orbit budget, with its limitation

For odd n let `n_j=U^j(n)`, `S_j=S(n_j)`, `K_j=K(n_j)`.
Division by powers of two does not change digit sum, so (8) becomes

    S_(j+1)=3S_j+1-K_j.

For every h>=1,

    S(n) = sum_(j=0)^(h-1) (K_j-1)/3^(j+1) + S_h/3^h.  (12)

Also K_j>=3: for odd inputs c1=2 and c2>=1. Since `U(x)<=2x` for x>=1,
`S_h<=floor(log2(n))+h+1`; the remainder tends to zero, independently of
whether Collatz terminates. Hence

    S(n) = sum_(j>=0) (K(U^j(n))-1)/3^(j+1).              (13)

This is a positive discounted carry budget for every fixed odd source,
including any hypothetical divergent one. Consequently its existence is
not a convergence certificate. The discount suppresses late high-scale
costs. The forest's remaining target is a constraint linking carry
**position, repeated scale crossings, and enclosing excursions**, strong
enough to control an undiscounted height debt. Neither (10) nor (13)
provides that inequality.

There is a full positional hierarchy behind (12). Let

    D_n(z)=sum_(s>=0) b_s(n) z^s,
    C_n(z)=sum_(s>=1) c_s(n) z^s.

The coefficient relation for multiplication by three and addition of one
is `b_s(3n+1)=3b_s(n)+c_s-2c_(s+1)`. Therefore, for odd n and
`a=v_2(3n+1)`,

    z^a D_(U(n))(z)=3D_n(z)+1-(2-z)C_n(z)/z.              (14)

The quotient is a polynomial because C_n has zero constant term.
At z=1 this gives the total-carry identity. At z=2 the carry coefficient
vanishes, leaving `2^a U(n)=3n+1`. At 1<z<2 it is an exact interpolation
between digit mass and integer value, with the actual halving count a
retained. It does not imply a sign for an orbit potential.

Differentiate at z=1. With `M(n)=sum s*b_s(n)` and
`H(n)=sum_(s>=1) s*c_s(n)`, the first position-sensitive identity is

    a S(U(n))+M(U(n)) = 3M(n)+2K(n)-H(n).                 (15)

This is a concrete replacement for the phrase "carry depth": H counts
where the carry is, and its exact update contains the true valuation a.
Higher derivatives keep further ordered-position moments. Any finite
truncation still requires its own hostile probe; (15) is a conservation
identity, not a new descent inequality.

## 7. Reproduction and audited scope

Run from the repository root:

    python -B 04-computation/experiments/forest_20260926_bernoulli.py
    python -B -O 04-computation/experiments/forest_20260926_bernoulli.py

[Script](../../04-computation/experiments/forest_20260926_bernoulli.py),
[retained output](forest_20260926_bernoulli.out). Normal and optimized
outputs agree. Universe and independent controls:

- Bernoulli numbers computed by exact rational recurrence through B256,
  compared with the denominator theorem at eight dyadic indices.
- Primality of the first five Fermat values by trial division; F5 by its
  displayed factorization. No test of all later Fermat numbers.
- Polygon catalogue versus independently computed phi(n), n=3..10000.
- Exact cyclotomic quotient-ring and B1 projectors through depth eight;
  all 2046 parity words of lengths 1..10, with a neighbouring wrong source.
- Carry path, B1 formula, inverse graph decoder and coefficientwise
  polynomial identity for n=0..20000; positional moment on 10000 odd inputs.
- 195 sources in three fixed fibres, each with 36 polynomial valuations;
  exact unbounded-family formulas, not empirical extrapolation.
- 1016 padded growth hostiles, L=9..1024.
- 1000 odd sources with forty-step discounted budgets and exact remainders.

Primary source scope: DLMF supplies the denominator identity; Hochster
supplies the polygon criterion; PrimePages supplies the checked catalogue
status. The arithmetic bridges, carry graph, fibre escape, and hostile
proofs are included here explicitly. No claim that the external statements
solve a Collatz proof obligation is made.
