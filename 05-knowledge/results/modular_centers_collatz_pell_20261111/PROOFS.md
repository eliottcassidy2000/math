# Modular centers, exact trajectory guards, and Pell birth ports

**11 November 2026. PROVED elementary statements below; CITED Dirichlet
corollary; FINITE-EXACT checks in verification.json. No priority claim.**
These are scoped arithmetic and graph-birth results. They do not prove
Collatz, Chan-Pak, infinitely many quadratic primes, or a classification
of all tournaments. [BOARD.md](BOARD.md) records inheritance and losses.

## P1. The central blocks and the four residue classes

Use the multiplication table on indices 1,...,m-1, with entries represented
by 0,...,m-1. For odd m>=3 the two middle indices are r=(m-1)/2 and
s=(m+1)/2. In Z/m, r=-1/2 and s=1/2, so the center is

    [ d  c ]       d = 4^(-1) mod m,   c = -4^(-1) mod m.
    [ c  d ]

Both representatives are positive and d+c=m. Thus the two copies of m
are **row sums** (also column sums), not additional entries of the 2x2
block. The user's numerical formulas are precisely

    m=4x+1: [ 3x+1   x   ]       m=4x-1: [   x    3x-1 ]
             [   x   3x+1 ]                 [ 3x-1    x   ].

This recovers the proved calculation in
[tiling_modular_atoms_20261003.md, sections 1–3](../tiling_modular_atoms_20261003.md).
For both signs sigma in {-1,+1}, put m=4x+sigma, x>=1, and b=3x+sigma.
The center is the unordered complementary pair {x,b}, each repeated twice;
the displayed positions retain the sign. Decode an edge by

    x -> b                  when x is odd,
    x -> (b-sigma)/6 = x/2  when x is even.

Consequently the user's four table classes encode exactly the raw map
C_sigma(x)=3x+sigma on odd inputs and x/2 on even inputs:

| m modulo 8 | sign | source | destination decoded from the other entry b |
|---:|---:|---|---|
| 5 | +1 | x odd | b=3x+1 |
| 1 | +1 | x even | (b-1)/6=x/2 |
| 3 | -1 | x odd | b=3x-1 |
| 7 | -1 | x even | (b+1)/6=x/2 |

The division by 6 is exact division of integers. It is not inversion of 6
in Z/m: 3 can divide m. "Half the edges" is correct as the partition by
input parity, with each class of inputs having natural density 1/2; it is
not a statement about the proportion of steps in a particular orbit.

There is no primality assumption. In fact both central row indices and
both central entries are units for every odd m. This immediately prevents
using central invertibility alone to distinguish prime from composite.

The earlier [paired-prime note, section 3](../paired_prime_collatz_20261003.md)
uses the *same pair* to encode two different shortcut-map inputs. That is
consistent with the present raw-edge decoder, but the input coordinate
must be kept: an expression 3x+1 is not by itself a choice of dynamics.

## P2. The table-modulus graph and its exact inverse gates

The bijection psi_sigma(x)=4x+sigma takes positive integers to
D_sigma={4x+sigma:x>=1}. Conjugating the raw map gives

    F_sigma(m) = 3m+2sigma       if m = 4+sigma mod 8,
                 (m+sigma)/2    if m = sigma mod 8.

These two cases preserve D_sigma. This proves, for every k>=0,

    F_sigma^k(4x+sigma) = 4 C_sigma^k(x)+sigma.

Every target m has the even-input parent 2m-sigma. The only possible
additional parent is (m-2sigma)/3, and it exists exactly when

    m=17 mod24   for sigma=+1;
    m= 7 mod24   for sigma=-1.

**Proof of the inverse condition.** In x coordinates the other parent
is y=(x-sigma)/3. It is a positive odd integer iff x=4 mod6 for sigma=+1
or x=2 mod6 for sigma=-1. Substitution m=4x+sigma gives the two residues.
These residues ensure positivity at their first positive representatives.
The even parent is always positive. The two parents, when both exist,
are distinct. This proves a complete incoming-edge classification.

Here 24 is the least period that simultaneously retains the modulo-4
coordinate and the modulo-6 odd-predecessor gate. It is a different
object from vertex 24 in a square-sum graph; their matching numeral does
not supply a map between them.

Controls: on the plus sheet 5->17->9->5 encodes 1->4->2->1. On the minus
sheet 3->7->3 encodes 1->2->1, and
19->55->27->79->39->19 encodes 5->14->7->20->10->5.

Changing sign while retaining x sends m to m-2sigma. This commutes with
the even branch, but on an odd input the two orders of operations differ
by 8sigma in modulus coordinates. Thus it is not a conjugacy of the two
positive dynamics. Negation *is* a conjugacy on the signed integer
domains: C_-(-x)=-C_+(x). Neither assertion proves termination.

## P3. A finite affine rule for every finite trajectory pattern

Here switch explicitly to the **shortcut** convention

    T_sigma(x)=x/2 if x even, (3x+sigma)/2 if x odd.

In modulus coordinates its step is

    Fhat_sigma(m) = 3^b (m+sigma)/2,  b=((m-sigma)/4) mod2.

Let w=(b_0,...,b_(L-1)) be any binary word, L>=1, with s ones. Start
A=1, B=0, D=1 and, reading bits in time order, update

    bit 0: (A,B,D) -> (A,B,2D),
    bit 1: (A,B,D) -> (3A,3B+D,2D).

At the end A=3^s, D=2^L. There is **exactly one** source residue modulo D
with this parity word, and it is

    r_sigma(w) = -sigma B A^(-1) mod D.                 (1)

For every positive x in that residue class,

    T_sigma^L(x) = (A x+sigma B)/D.

In table coordinates all such sources, and exactly those sources, satisfy

    m = 4r_sigma(w)+sigma mod 4D,
    Fhat_sigma^L(m) = (A m+sigma K)/D,
    K = 4B+D-A.                                       (2)

**Proof.** Composition gives the affine numerator by induction. For the
parity assertion, a length-L prefix occupies one residue modulo 2^L.
Its two lifts modulo 2^(L+1) have endpoints differing by the odd number
3^s, so exactly one has each next parity. The induction begins with
the even and odd classes modulo 2. Thus every word has one residue;
integrality of its endpoint and invertibility of odd A identify that
residue as (1). Substitution x=(m-sigma)/4 gives (2).

The affine numerator **without its source guard** is not an exact graph
edge description. The guard is the retained compatibility coordinate.
Composition consumes parity precision: this is not a deterministic map
on a single fixed finite residue ring. In particular x and x+2^L can
have endpoints of opposite parity after L steps.

This is an elementary finite form of the established parity-vector
machinery; see Bernstein–Lagarias, *The 3x+1 conjugacy map*,
[primary paper](https://websites.umich.edu/~lagarias/doc/bernstein.pdf),
and Laarhoven–de Weger, *The Collatz conjecture and De Bruijn graphs*,
[arXiv:1209.3495](https://arxiv.org/abs/1209.3495).
Their finite graph relation and 2-adic constructions do not settle
positive-integer convergence.

**CITED corollary (Dirichlet).** Every finite word w on either sheet has
infinitely many **prime table moduli** and infinitely many composite table
moduli. Indeed the residue in (2) is odd, hence coprime to 4D; apply
[Dirichlet's theorem, original paper in translation](https://arxiv.org/abs/0808.1408).
For composites choose any positive representative h>1 of the same
residue and use m=h(1+4Dt), t>=1. Both constructions respect the domain
after discarding any finite initial exceptions.

For the all-ones word, the source classes are x=-sigma mod2^L,
or m=-3sigma mod2^(L+2). Taking sufficiently large representatives
produces L consecutive strictly increasing shortcut steps, even with
prime table modulus. Thus a bounded parity window plus initial-table
primality cannot exclude long finite growth. Nothing here says that
every intermediate table modulus is prime or that one integer realizes
an arbitrarily prescribed infinite positive orbit.

## P4. What prime tables optimize, with a specified information model

Fix m>=2 and row multiplier a in Z/m. Put g=gcd(a,m), including
gcd(0,m)=m. The row map y->ay has image size m/g and each nonempty fibre
has size g. To prove this, write a=g a', m=g m', and use invertibility
of a' modulo m': solutions of ay=t exist exactly for g|t and form one
residue class modulo m', with g lifts modulo m.

Hence the nonzero rows are all permutations **iff m is prime**. The
number of reversible rows is phi(m). This is maximal reversible
coverage, not a minimization of the number of entries in a table.
The full two-input product map is many-to-one even over a prime field.

For an explicitly uniform column Y in Z/m and a known fixed multiplier,

    H(Y | aY) = log_2 gcd(a,m),
    H(aY) = log_2(m/gcd(a,m)).

If the multiplier is also selected uniformly, but revealed to the decoder,
the average lost information is

    (1/m) sum_a log_2 gcd(a,m)
      = sum_(p^e || m) [(1-p^(-e))/(p-1)] log_2 p.       (3)

For each prime p the expected exponent of p in gcd(a,m) is
sum_(j=1)^e Pr(p^j|a)=sum p^(-j), proving (3).
For a uniform **nonzero** known multiplier the mean loss is
(m times (3) - log_2 m)/(m-1); this vanishes exactly for primes.
No universal Shannon claim is made without these sampling and observation
choices. This refines the inherited row-map calculation in
[modular_fourier_zeckendorf_routes_20261003.md](../modular_fourier_zeckendorf_routes_20261003.md).

Two useful complete finite signatures are

    number of rows with gcd(a,m)=d: phi(m/d), d|m;
    number of ordered products ab=t:
      N_m(t)=sum_(d | gcd(t,m)) d phi(m/d).              (4)

In particular the zero-product count Z(m)=sum_(d|m) d phi(m/d), and the
number of ordered triples (a,b,c) with ab=ac is m Z(m). The latter
counts within-row collisions, including b=c; it is not the global
collision count of two independent row/column pairs.

For m=25 there are 20 reversible rows, four nonzero defective rows, and
the zero row. Row 5 is 5-to-1 onto {0,5,10,15,20}, while its central rows
12 and 13 remain permutations. For m=r^2, row r is r-to-1 onto r outputs
and r is a nonzero square-zero element when r>1. Z/25 is not the field
F_25: the latter has 24 reversible nonzero multipliers.

**Finite composition rule.** Chinese remaindering decomposes the actual
table into its prime-power tables. At p^e, multiplying row labels adds
their p-adic valuations and caps the sum at e. Valuation records suffice
for fibre sizes, but exact values also need the unit coordinates and
their residues. This classifies the row-map family for every modulus
without storing its m^2 entries; it does not supply free factorization.

## P5. How far from the center is the first singular row?

For an odd modulus and odd offsets u,v in
{-(2R-1),...,-1,1,...,2R-1}, the centered 2R by 2R window uses indices
(m+u)/2 and (m+v)/2, for 1<=R<=(m-1)/2. Its entry is

    uv * 4^(-1) mod m.

The row at offset u is reversible iff gcd(u,m)=1. Therefore if m is
composite with least prime divisor p, the first radius containing a
singular row is exactly

    R_first(m)=(p+1)/2.

If m is prime, there is no singular row in the nonzero-index table.
For 25 the first defect has radius 3 and offset 5.

For every fixed R there are infinitely many composites p^2 with p>2R-1
whose window has the same all-unit mask as a prime table. This is a
limitation of the **unit/nonunit mask**, not an information-theoretic
limitation of the numeric entries: the central row sum already recovers
m. Scanning until p<=sqrt(m) is trial division in different coordinates,
not a new faster primality algorithm.

For completeness the even m=2k>=4 middle 3x3 block, with indices
k-1,k,k+1, is

    [ d+1   e   d-1 ]
    [  e    d    e  ]  modulo m,
    [ d-1   e   d+1 ]

where d=0 if k is even and d=k if k is odd, and e=k-d. Expansion of
(k+i)(k+j) proves the formula. Thus every triple {4t,4t+1,4t+2}, t>=1,
has three explicit center templates, with different window sizes:

    m=24: [ 1 12 23 ]   m=25: [19  6]   m=26: [14  0 12]
          [12  0 12 ]         [ 6 19]         [ 0 13  0]
          [23 12  1 ]                        [12  0 14].

## P6. Pell plateaus retain explicit parent ports

Let Q_n be the simple undirected graph on {1,...,n}, with i adjacent to
j iff i+j is a square. Its new vertex n has parent set
P(n)={q^2-n: n<q^2<2n} and degree

    d(n)=floor sqrt(2n-1)-floor sqrt(n).

The inherited plateau classification
[Session 148, PELL_PLATEAUS.md](../session148_square_sum_20261103/pell_plateaus/PELL_PLATEAUS.md)
states that square centers R=r^2>=4 with d(R-1)=d(R) are exactly

    epsilon=2: r+t sqrt2=(3+2sqrt2)^j, j>=1; s=2t;
    epsilon=1: s+r sqrt2=(1+sqrt2)^(2j+1), j>=1;
    in both: s^2=2R-epsilon.

They also have d(R)=d(R+1). Here is the retained endpoint refinement:

    P(R-1) = {1} union {q^2-R+1: r<q<s},
    P(R)   = {q^2-R: r<q<s} union {R-epsilon},
    P(R+1) = {q^2-R-1: r<q<s} union {R-epsilon-1}.       (5)

To prove (5), list the possible square roots: at R-1 they run from r
through s-1; at R and R+1 from r+1 through s. The endpoint inequalities
follow directly from s^2=2R-epsilon, with epsilon=1 or2 and r>=3.
This proves the formulas without equating vertex identities just because
the three degrees agree. Each degree is s-r.

For R=25, the actual ports are

    P(24)={1,12}, P(25)={11,24}, P(26)={10,23}.

Equal birth counts therefore retain neither the endpoint labels nor the
Hamiltonian continuation state. This is the relevant link to the prior
[square-birth window packet](../square_birth_windows_20261103/README.md):
enumeration must carry these ports, not just the scalar 2.

## P7. Odd-square centers join the plateau graph to table dynamics

Every odd square R=(2h+1)^2 with h>=1 has

    R=4h(h+1)+1, x=h(h+1), b=3h(h+1)+1.

It lies in the 1 mod8 table class and decodes the halving edge

    h(h+1) -> h(h+1)/2,

whose destination is triangular. Furthermore an odd square is 1 or9
mod24, never17. By P2 it has **exactly one** incoming table modulus,
2R-1, and its outgoing modulus is (R+1)/2. Thus

    2R-1 --> R --> (R+1)/2                            (6)

is a completely classified local segment for every odd square R>=9.
The first arrow is the even raw-Collatz branch: the corresponding
integer sources are 2x -> x -> x/2.

For a Pell plateau, the same incoming modulus has a sharper type:

    epsilon=1: 2R-1=s^2       (always composite),
    epsilon=2: 2R-1=s^2+1     (quadratic prime candidate).

Examples: R=25 has 49->25->13, encoding 12->6->3. R=9 has
17->9->5, encoding 4->2->1; R=289 has 577->289->145.

Consequently the epsilon=2 branch supplies a specified subsequence of
**candidates** for [A002496](https://oeis.org/A002496), the primes k^2+1.
By the proved inverse-proper-divisor-sum criterion in
[infinite_birth_families_20261111/ARITHMETIC.md](../infinite_birth_families_20261111/ARITHMETIC.md),
it also has the exact equivalence

    2R-1 prime  iff  A135244(2R)=(2R-1)^2.              (7)

This is a typed bridge with an explicit map and preserved predicate.
It makes no claim that square-sum Hamiltonian paths correspond to
Collatz orbits, and no claim that Pell centers exhaust quadratic primes.

## P8. The incoming Pell trace: pruning primes, retaining infinite factors

Write r_j+t_j sqrt2=(3+2sqrt2)^j, j>=0, and
M_j=2r_j^2-1. Squaring the Pell unit gives M_j=r_(2j), so

    M_0=1, M_1=17, M_(j+2)=34M_(j+1)-M_j.

This is [A056771](https://oeis.org/A056771), beginning
1,17,577,19601,665857,22619537,768398401,... . The recurrence needs two
integers to generate this entire infinite branch. The underlying trace
divisibility mechanism was already proved in
[thirtysix_bridge_20260925.md, prime-index obstruction](../thirtysix_bridge_20260925.md).

If j=dq with q>1 odd, then M_j=T_q(M_d), where T_q is the integer
Chebyshev polynomial. Its odd parity implies x divides T_q(x). Since
M_d>1 for d>=1 and the traces strictly increase, this is a proper
factorization. Hence

    M_j prime => j is a power of 2,  for j>=1.          (8)

This is a necessary condition only. Let z_a=M_(2^a), a>=0. Then

    z_0=17, z_(a+1)=2z_a^2-1,
    z_1=577, z_2=665857,
    z_3=886731088897=257*1409*2448769.

This is a tail of [A001601](https://oeis.org/A001601); it is not a new
sequence. Also M_3=19601=17*1153 is an inherited control for (8).

**Infinite survivor.** The z_a are pairwise coprime: modulo z_a the next
term is -1, the following term is 1, and every later term stays 1.
All terms exceed 1 and strictly increase. Therefore this family has
infinitely many distinct prime divisors, although infinitely many
**prime terms** remains unproved here. In fact every prime divisor is
1 mod4, since each z_a=s^2+1 with s even, so s has multiplicative order
4 modulo that divisor. This supplies infinitely many prime factors of
the desired residue type without conflating factors and terms.

## Verification and limits

Run `python check.py` from this directory; `python -O check.py` produces
identical JSON and leaves every exception-based check active. The check
compares literal multiplication tables with the formulas, brute incoming
edges with the residue gates, independently generated parity words with
the affine guards, literal square-sum parent sets with (5), and Pell-pair
generation with the scalar trace recurrence. Exact factor certificates
are checked by trial division. No network or external library is needed.

The proofs establish the quantified results. The finite controls inspect
all odd moduli through10001 and even moduli through10000, 100000 raw
inputs per sign, 10000 incoming targets per sign, every word of lengths
1..12 on both sheets, every table of order2..201, centered singularity
profiles through1001, and 100 terms of each Pell family. Literal parent
sets are checked at the generated centers <=1000000; the trace
divisibility control runs through index128.

These are exact reductions of specified infinite families. The missing
global obligations are integer-orbit descent, prime values along thin
recurrences, and a compatible tournament realization of any proposed
path-count operation. None is replaced by a matching numeral or a finite
pattern census.
