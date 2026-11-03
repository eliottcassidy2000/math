# Prime arrival, triangular pairs, and complete unit-wheel run formulas

**PROVED elementary formulas; FINITE-EXACT prime searches.** Primality is
always distinguished from local coprimality. This extends the recent
[two-prime wheel theorem](../session146_carrier_primegame_20261103/prime290/WHEEL_GENERAL.md),
not its first discovery. Its source packet and the earlier
[prime-arrival audit](../positive_source_jacobian_prime_clocks_20261103/ARITHMETIC_CLOCKS.md)
already prove the 290/227 facts used here.

## P1. All even moduli and every shared-port pattern

Let m be even, P its distinct odd prime divisors, and

    L=m/(2 product_(p in P) p).

For a nonempty finite labelled offset set D, consider odd residues a mod m
for which every a+2d, d in D, is a unit. Let nu_p(D) be the number of
distinct residues of D modulo p. Then the exact number is

    A_m(D)=L product_(p in P) (p-nu_p(D)).             (P1)

If P is empty, this is m/2 for every D. The proof is CRT: oddness allows
2^(v_2(m)-1) choices, and modulo p^e exactly nu_p(D) residue classes are
forbidden, each with p^(e-1) lifts. Repeated named ports do not create new
forbidden residues. Appending a port changes only the factors where its
residue was previously unoccupied; P1 is an exact shared-port recurrence.

For arbitrary future offsets the number A_m(D) is not a closed carrier.
Modulo 10, D={0,1} and E={0,2} both have count three. Append offset 1:
the first still has count three and the second has count two. One must
retain the occupied residue sets (or an equivalence proved sufficient for
the allowed future ports), not merely their cardinalities. This is the
arithmetic counterpart of checking whether an operation descends through
the line quotient; it is not an identification of the two state spaces.

## P2. Exact length distribution of every +2 unit wheel

Make a directed graph on the odd units mod m with edge a->a+2 whenever
both are units. If P is empty it is one directed cycle of length m/2.
Otherwise put p0=min(P) and, for r>=1, define

    S_r=L product_(p in P)(p-r),  if r<p0;
    S_r=0,                      if r>=p0.             (P2)

S_r counts starts of r consecutive unit vertices, by P1 with offsets
{0,...,r-1}. The graph is a union of directed paths, each of length at
most p0-1 vertices. The number of its maximal paths with exactly ell
vertices is

    N_ell=S_ell-2S_(ell+1)+S_(ell+2),  ell>=1.         (P3)

**Proof.** An ell-window is a maximal run iff its predecessor and successor
are both holes. Inclusion-exclusion subtracts the two (ell+1)-windows
and restores the (ell+2)-window. A hole exists because P is nonempty, so
no cyclic component is missed. Equivalently

    S_r=sum_(ell>=r) (ell-r+1)N_ell,

and P3 is the inverse second difference. This gives the complete path
classification at every even modulus, including prime powers and repeated
wheel periods, with the power-of-two cycle exception explicitly separated.

When P has h primes, the interior histogram is the second difference of
a degree-h polynomial, hence has degree at most h-2 away from truncation.
For two primes p<q this recovers two paths of each length 1,...,p-2 and
q-p+1 paths of length p-1. For three primes p<q<r:

    N_ell=2(p+q+r)-6ell-6,  1<=ell<=p-2,
    N_(p-1)=(q-p+1)(r-p+1),                            (P4)

for squarefree m=2pqr; multiply by L for repeated prime powers or a
higher power of two. There is no empirical extrapolation in P4.

| m | Exact maximal run histogram, vertices per run : count |
|---:|---|
| 290=2*5*29 | 1:2, 2:2, 3:2, 4:25 |
| 1160 | 1:8, 2:8, 3:8, 4:100 |
| 5510=2*5*19*29 | 1:94, 2:88, 3:82, 4:375 |
| 2310 | 1:210, 2:135 |
| 30030 | 1:2790, 2:1485 |

The primorial examples connect the older 2,6,30,210,2310,30030 work by
an actual multiplicative sieve operation. They do not assert a result
about jumping champions: these are finite residue graphs, not prime gaps.

## P3. The precise triangular pictures and their clocks

There are two related but different triangular arrays.

1. **Odd pairs:** p=s-h, q=s+h with 0<=h<s and s,h of opposite parity.
   They represent positive odd integers (allowing p=q at h=0). The actual
   prime-pair predicate is primality of both entries. Its local wheel mask
   replaces that by coprimality to m; this is necessary for primes not
   dividing m, never sufficient. Fixed s is a sum row; fixed h is a gap row.
   The inherited machine chooses odd entries >=3, so its finite domain
   additionally excludes entries equal to one.
2. **First-arrival array:** for an odd prime p, the even moduli hit by its
   residue +/-1 are exactly the even divisors of p-1 or p+1. All entries
   p<m-1 are blank. Processing primes in increasing order and marking each
   modulus at its first divisor incidence is an exact algorithm for the
   first-arrival table. The divisor mask must keep that processing order.

The boundary pair (km-1,km+1) joins these pictures when BOTH entries are
prime: its center is km and its half-gap is one. Individual first hits
need not have the same k. For m=290 the exact values are:

    first +1 prime: 1451=5*290+1;
    first -1 prime: 2029=7*290-1;
    first crossing twin pair: (6089,6091)=(21*290-1,21*290+1).

The eligible rank of 1451 is 227: pi(1451)=230 minus the three prime
divisors 2,5,29. Its ten earlier progression candidates are composite,
as certified in the inherited packet. This deterministic first-arrival
clock is not a geometric random variable with mean phi(290)/2=56.

The prime 2 is exceptional for sums: the sum of two odd primes is even,
whereas 2+p=q prime means a twin-prime step. Keeping even sum nodes gives
the pair-incidence picture; keeping only prime sum nodes gives the much
thinner twin-step picture. There is a positive coverage statement: by
the first +/-1 hit, every prime needed for an even sum row <=m has
already appeared, since its largest possible odd summand is m-3. Thus
all such finite rows can be computed by that clock, including rows with
no representation. Incoming [hourglass/port work](../hourglass_quadratic_ports_20261103/ARITHMETIC.md)
already proves this saturation and the streaming convolution recurrence.
It does not assert nonemptiness of every row or identify residues with
actual prime summands. The parallel [odd-pair note](../quadratic_2916_odd_pair_triangles_n4_20261102.md)
records the concrete alias (5,47) modulo26: residues (5,21) sum to26,
but the actual primes sum to52 and21 is composite. That note was still
marked pre-audit when read; we retain only the directly checked alias here.

For finite prime flags x_i on an odd-number interval, the pair predicate
is x_i*x_j, while the layer lines observe x_i XOR x_j. XOR observations
identify complementary flag words. A known prime anchor, e.g. the entry
3, fixes that bit and restores every flag along a connected line graph.
The anchor restores information; it does not compute previously unknown
primality. Hourglass relations test consistency of that reconstruction.

The existing [finite triangle machines](../session146_carrier_primegame_20261103/machines/README.md)
already compile these flags into ordered FRACTRAN and exact TM/CA models,
and obtain small prime-flag inputs from both versions of PRIMEGAME.
This session extends the arithmetic run formulas and supplies the explicit
quadratic machine in [HOURGLASS.md](HOURGLASS.md); it does not relabel the
finite triangle compiler as a uniform Goldbach decision procedure.

## P4. Modulo 290 remains a hostile control

The inherited quadratic numerator map F(a)=(a^2-29)/4 descends to the
145 odd residues mod290. Its eight modular three-cycles are explained by
CRT mod5 and mod29, including the squaring map on F29. Only one cycle
lifts to a bounded odd-integer orbit. The integer orbit

    9 -> 13 -> 35 -> 299 -> ...

escapes although its residues run on 9->13->35->9. Thus neither the
eight modular cycles nor the finite wheel supplies the missing integer
bound. [DYNAMICS.md](DYNAMICS.md) supplies that bound independently.
Likewise the known 81/27 tournament pair still needs named-port data even
when cyclic triples and unmarked path-cover profiles agree. Numerical
coincidences among these different graphs are not identifications.
