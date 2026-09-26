# 223: Fermat sextics, finite-prime returns, and exact affine genealogy

**Status: PROVED elementary transfers; CITED Hasse-Weil/genus input plus
FINITE-EXACT exhaustive census gives a complete local classification;
CITED unit-equation bound gives a PROVED fixed-word specialization theorem.
Collatz termination and chronological paired-valuation independence remain
OPEN. No claim of literature priority.** Geometry lane, 2026-09-26.

Reproduce with
`python 04-computation/experiments/crossroads223_20260926_geometry.py --rank`.
The matching `.out` records the script SHA256. The program uses only Python's
standard library and exact integers/Fractions. No random sample is involved.

## Inheritance and concept board

Anchor: preserve the odd Collatz exponent word together with integer height.
Niche: multiplicative characters, finite-field curves, and prime revisits.
Wildcard: the genealogy of rational and shifted integer cycles.

* PROVED closest local mechanism: THM-4473,
  `01-canon/theorems/THM-4473-collatz-digit-chains-rotation-repunits.md`.
  Residue chains under Haar-distributed initial sources do not imply temporal
  randomness on a single deterministic integer orbit.
* PROVED integer-cycle source: THM-4484,
  `01-canon/theorems/THM-4484-free-and-sporadic-cycles.md`.
  An affine word with multiplier numerator A, denominator B and carry C has
  rational periodic point C/(B-A); scaling the additive shift can make it
  integral. The negative-sheet cycles remain hostile controls.
* PROVED geometric source: THM-3334,
  `01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md`.
  Here 223 is a finite record *parameter*: c_223=99905=5*13*29*53 gives
  eight primitive Gaussian allocations. The prime 223 is 3 modulo 4 and
  cannot divide any spine norm t^2+(t+1)^2. These two roles must not be merged.
* PROVED valuation source: THM-3131,
  `01-canon/theorems/THM-3131-prime-resonance-newton-slope-separation.md`.
  Its 223 is merely an endpoint of a finite prime control range. The useful
  inherited operation is local noncancellation, not numerical recurrence.
* FINITE-EXACT old sidecar:
  `05-knowledge/results/collatz_mod6_20260917_zsigmondy_triad.out` already has
  2^37-1=223*616318177. Recovering this factor is provenance, not new work.
* OPEN carrier: G2 in
  `05-knowledge/results/crossroads_20260926_geometry.md`, which keeps both
  coordinate slots of (-3m,3m+1). A generic set of eight distinct sources and
  distinct successors already has a multiplicative dependence; chronology
  is essential. A new prime at every rising step is also false (55,83,125).

The live concepts are finite characters, affine carry, integer scale,
chronological prime reuse, and rational-cycle normalization. The Hadamard
order 892 construction at q=223 supplies no map preserving a Collatz target;
its row sums alone are not a certificate. This lane does not force a
tournament onto these objects.

## 1. Complete classification of empty Fermat sextics

Let C_p be x^6+y^6+z^6=0 in P^2(F_p), and let C_p^* denote its part with
xyz nonzero. The following is a **local arithmetic-geometric theorem**:

```
C_p(F_p) is empty iff p in {7,31,67,79,139,223}.
C_p^*(F_p) is empty iff
 p in {2,5,7,13,31,61,67,79,97,139,157,223,277}.
```

Thus 223 is the largest prime with no projective point, while 277 is the
largest with no point in the torus. This distinction matters for Collatz
unit cosets: coordinate-zero points cannot be used there.

**CITED input.** For p>3 the partial derivatives show C_p is smooth, and
the smooth plane sextic has genus (6-1)(6-2)/2=10. The Hasse-Weil bound gives

```
N_p = #C_p(F_p) >= p+1-20 sqrt(p).                         (1)
```

The author's primary lecture notes provide both the Fermat genus statement
(Example 21.2.4) and the curve theorem with point-counting consequence
(Theorem 21.3.1 and Remark 21.3.2, take extension degree one):
[Piotr Achinger, Algebraic geometry, Lecture 21, pp. 3-4](https://achinger.impan.pl/ag2026/lec21.pdf).
The genus and Hasse-Weil theorem are cited inputs, not proved anew here.

**PROVED reduction to a finite census.** At p=401,
`(p+1)^2-400p=1204>0`; this quadratic increases for p>=401. Hence (1) is
positive for all p>=401. A coordinate-zero projective point has exactly
one zero coordinate, and on each of the three axes there are at most six
solutions. Thus there are at most 18 boundary points. At p=439,
`(p-17)^2-400p=2484>0`, increasing thereafter, so (1)>18 for p>=439.
There are no primes between 433 and 439. It suffices to check primes<=397
for projective emptiness and primes<=433 for torus emptiness.

**FINITE-EXACT complete census.** The script checks every prime<=433.
Put d=gcd(6,p-1), H_p={x^6:x in F_p^*}, and
`R_p=#{u in H_p : -1-u in H_p}`. Normalizing z=1 gives exactly

```
#C_p^* = d^2 R_p,
#C_p   = d^2 R_p + 3d [ -1 in H_p ].                     (2)
```

These formulas hold even at p=2,3: the map on F_p^* still has exactly d
preimages per image. The finite census yields the displayed complete lists.
A second implementation directly enumerates the disjoint charts z=1 and
z=0,y=1 for all six empty projective primes and eight positive/boundary
controls. At p=277 it finds exactly 18 projective points and no torus point.

The singular-characteristic controls are explicit. At p=2, [1:1:0] is a
projective point, but [1:1:1] is the only possible torus point and fails;
the exact counts are (3,0). At p=3, [1:1:1] is a torus point, with exact
counts (4,4). No smooth-curve bound is applied at either characteristic.

This proof's finite component is a transparent exhaustive certificate, not
a floating-point threshold or a finite cutoff being silently extrapolated.

## 2. What the sextic actually says about a Collatz step

Write U(n)=(3n+1)/2^k on positive odd integers, k=v_2(3n+1). Exact arithmetic
gives

```
ord_223(2)=37,  ord_223(3)=222,  3^6=2^21 mod 223.
H=<2>=(F_223^*)^6,   |H|=37,   -1 notin H,   -3 notin H.
```

Let C_a=3^a H with a modulo 6. Division by 2^k preserves each coset, so
the possible nonzero class transitions are determined by 3n+1 alone.
The exact 6-by-6 count matrix, recording each nonzero source residue once
and omitting the single zero-target residue, is

```
5 4 7 6 8 7
4 7 5 6 7 8
9 5 4 9 5 4
5 6 7 8 4 7
4 7 6 8 7 5
9 8 8 0 6 6
```

The sole forbidden nonzero edge is C_5 -> C_3. Here 3C_5=H and C_3=-H,
so that edge would say u+1+v=0 with u,v in H, exactly a torus point on
C_223. Conversely such a point would give the forbidden edge. Every
nonzero matrix entry has an integer-step realization by the Chinese
remainder theorem together with an exact dyadic valuation cylinder.

This is the claimed connection with all its coordinates: source is the
Fermat sextic, target is a deterministic forbidden odd-Collatz residue
transition, map is n -> (class(n),class(3n+1)), preserved predicate is
nonvanishing plus the affine equation, lost information is exponent size
and integer height, needed sidecar is the exact valuation word, and the
cheapest hostile probe is to retain the zero-target residue and compare
the minus sheet. Under n -> -n, class labels shift by 3; the minus-sheet
forbidden edge is C_2 -> C_0. Consequently this local obstruction does not
choose the allowed positive Collatz cycle over negative-cycle genealogy.

## 3. A finite-prime/real-height decoder

**PROVED.** If m=U(n) is a unit modulo 223, endpoint residues determine

```
2^k = (3n+1)/m mod 223,
```

and hence k modulo 37. If in addition 3n+1<2^37, then 1<=k<37 and the
exponent is recovered exactly. This combines a real bound and an odd-prime
observer without a probabilistic assumption.

The precision lift is also exact. The certificate checks
v_223(2^37-1)=1; the binomial valuation formula gives

```
ord_(223^t)(2)=37*223^(t-1), t>=1.                       (3)
```

If a=v_223(m)=v_223(3n+1), keep both endpoints modulo 223^(a+t), cancel
223^a, and use (3). A bound k<37*223^(t-1) again makes k unique.
Retaining only residues modulo 223 when m=0 modulo 223 instead produces
0/0 and cannot decode anything. This is the explicit missing precision
coordinate at a prime revisit. It supplies observation, not descent.

## 4. Prime-return budget: a proved gate and its sharp hostile

For a positive word w=(k_1,...,k_r), let S_i=k_1+...+k_i, S_0=0. Then

```
U_w(n)=(3^r n+C_w)/2^S_r,
C_w=sum_(i=0)^(r-1) 3^(r-1-i) 2^S_i.                    (4)
```

If both endpoints are divisible by 223, then C_w=0 modulo 223. There is
no one-step return, and no two-step return, since the latter would require
2^k_1=-3 modulo 223. A three-step return requires
`9+3*2^a+2^(a+b)=0 modulo 223`; the exact allowed (a,b) modulo 37 have
least positive representatives

```
(10,1),(15,5),(21,3),(19,11),(8,23),(24,12),(9,32).
```

**PROVED exhaustive bounded gate.** If an actual positive odd trajectory
goes from a 223-multiple to a 223-multiple using at most eight total
halvings, it strictly descends. The exhaustive set of positive
compositions with total S<=8 has 255 elements. Testing C_w modulo 223
leaves only (2,3,1,1) and (2,3,1,2), with respective affine maps

```
(81n+223)/128,   (81n+223)/256.
```

Both are below n for every positive multiple n of 223. This proves the
claim for all integer heights, because the finite certificate exhausts
words, not initial integers. The first return has a concrete realization
`46161 -> 34621 -> 12983 -> 19475 -> 29213`.

**REFUTED strengthening: every 223-return descends.** Already with nine
total halvings the word (1,1,1,1,2,1,1,1) has A=6561, B=512 and
C=6913=31*223. Its least positive 223-divisible realization is

```
14495 -> 21743 -> 32615 -> 48923 -> 73385 -> 55039
      -> 82559 -> 123839 -> 185759.
```

All intermediate terms are not divisible by 223, so this is a first
return as well as an expanding return. The failed implication is from
local additive cancellation to real contraction; the missing coordinate
is the length r relative to total valuation S. The strongest survivor is
the sharp budget-eight gate, not a uniform drift assertion.

For completeness every fixed word has infinitely many actual positive
223-divisible sources. Its exact source residue is
`n=(2^S-C_w)*(3^r)^(-1) modulo 2^(S+1)`, followed by CRT with n=0
modulo 223. Thus the expanding hostile cannot be dismissed as a residue
path that has no integer realization.

## 5. Why 223 appears as a carry: a rational-cycle genealogy

The minimum-cost return word w=(2,3,1,1) has

```
C_w=27+36+96+64=223,  B-A=128-81=47.
```

It therefore gives the rational plus cycle

```
223/47 -> 179/47 -> 73/47 -> 133/47 -> 223/47.
```

Scaling by 47 gives the integer 3n+47 cycle
`223 -> 179 -> 73 -> 133 -> 223`. Its half-step shape is (P,A)=(7,4)
and its gap 47 divides the additive shift, exactly the free-cycle
mechanism of THM-4484. This is not a positive integer 3n+1 cycle.
The earlier repository occurrence 223/45 in the transversality foundry
has a different eventually periodic word and denominator.

There is an elementary rigidity of the *exact carry*. For fixed r,

```
C=3^(r-1)+2^k_1 C_tail,    C_tail odd,
```

so k_1=v_2(C-3^(r-1)), after which the tail is determined recursively.
The final exponent never enters C. For C=223, r<=5; the decoder succeeds
only for r=4 and internal exponents (2,3,1). All words of carry exactly
223 are therefore (2,3,1,k), k>=1. Their gaps are 2^(6+k)-81. None divides
223, since positive divisors are 1 and 223, requiring powers of two 82
or 304. Hence no integer 3n+1 cycle has this unreduced carry exactly 223.
This deliberately says nothing about carries that are larger multiples
of 223 or about cycles with another numerator.

## 6. A fixed-word specialization theorem for paired rank

**PROVED application of a CITED unit-equation theorem.** Fix a positive
exponent word w of length N. Among the infinitely many positive integers
realizing this word, only finitely many have dependent first-slot values
`3m_0,...,3m_(N-1)`. In particular, only finitely many can violate G2's
paired multiplicative independence. The statement fixes the entire word
before quantifying over its integer realizations; no uniform cutoff in N
or exceptional-height bound is asserted.

Here is a quantitative exceptional-*count* version. Choose the least
positive exact cylinder representative a and write its sources as
`n=a+2^(S_N+1)t`, t>=0. Each node is an integer affine form

```
L_j(t)=a_j t+b_j,
a_j=3^j 2^(S_N+1-S_j),
b_j=(3^j a+C_j)/2^S_j.                                 (5)
```

For i<j the determinant `D_ij=a_j b_i-a_i b_j` is nonzero, because

```
C_j/3^j=sum_(ell<j) 2^S_ell / 3^(ell+1)
```

is strictly increasing. Thus the forms have distinct roots even over the
function field Q(t). This gives immediate first-slot independence over
Q(t); the remaining issue is arithmetic specialization.

Let S be the finite prime set dividing
`6 product_(0<=i<j<N) |D_ij|`, and put s=|S|. For N>=2 the number of
exceptional nonnegative integer t is at most

```
binom(N,2) * 2^(16s+8).                                 (6)
```

This bank is intrinsic to the subwords. If C_(i,j) is the affine carry of
the subword from node i to node j, then direct substitution in (5) gives

```
D_ij=-3^i 2^(S_N+1-S_j) C_(i,j).
```

Thus S is exactly {2,3} together with primes dividing the subword carries.
For t repetitions of the critical return block (2,3,1,1), its block carry
is `223*(128^t-81^t)/47`. This is a precise connection to the old
primitive-divisor lane, concerning carries rather than a guarantee of
fresh primes in actual orbit nodes.

Proof: suppose `product_j (3L_j(t))^c_j=1` for a nonzero integral vector c.
If p outside S divides L_i(t), it cannot divide any other L_j(t), since
a common divisor divides D_ij. Its valuation then forces c_i=0. Hence
every node in the nonzero support of c is an S-unit. This support has at
least two indices: a single integer 3L_i(t)>=3 cannot generate a relation.
Choose two supported indices i<j. The identity

```
x+y=1,
x=a_j L_i(t)/D_ij,   y=-a_i L_j(t)/D_ij                 (7)
```

places (x,y) in the square of the rational S-unit group, of torsion-free
rank 2s. The map t -> (x,y) is injective. Beukers--Schlickewei bounds the
number of such pairs by 2^(8(2s)+8); union over the binom(N,2) index pairs
gives (6). For N=1 independence is immediate.

The type boundaries deserve explicit checks. The affine forms can never
be proportional because their roots are distinct. Their slopes are also
distinct (a nontrivial power of 3 cannot equal a power of 2), so an actual
collision L_i(t)=L_j(t) occurs at at most one parameter for each pair.
Every such collision among the N source nodes is already covered by the
exceptional set: its common integer value divides D_ij, making both
values S-units. A collision involving only the final node m_N does not
necessarily destroy rank; the one-step word (2) at n=1 gives the single
pair (-3,4), which is independent despite its closing loop. Two laps
repeat the source and are dependent, again covered by the finite set.
Finally, a node equal to 1 is an S-unit and presents no exception to the
argument. The first-slot value is 3 rather than 1, which is why a
one-index relation is impossible. No injectivity assumption is needed
for the fixed-word finite-exception theorem itself.

The cited input is Theorem 1.1 of F. Beukers and H.P. Schlickewei,
*The equation x+y=1 in finitely generated groups*, Acta Arithmetica 78
(1996), 189--199: [primary first page and paper](https://matwbn.icm.edu.pl/ksiazki/aa/aa78/aa7826.pdf),
[publisher record](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/78/2/109078/the-equation-x-y-1-in-finitely-generated-groups).
The exact bound is verified from the primary indexed first page. No
novelty claim is made for this elementary specialization application.

For clarity about both slots and torsion, a dependence among the actual
transition pairs `(-3m_j,2^k_j m_(j+1))` requires

```
3^(sum c_j) product m_j^c_j = 1,
sum k_j c_j = 0,   product m_(j+1)^c_j = 1,
sum c_j = 0 modulo 2.                                  (8)
```

The last line is the first-slot sign constraint. Rational valuation rank
forgets this torsion, but any integral kernel vector can be doubled to
satisfy it. Consequently full rational valuation rank is equivalent to
multiplicative independence of the pairs; omitting torsion does not turn
a genuine kernel into independence.

The theorem isolates a useful necessary relation: a putative exceptional
source must make at least two affine nodes S-smooth, where S is determined
by the fixed subword carries/resultants. It also explains why high-height
lifts tend to have full rank. It does not control the exceptional source
at an increasing sequence of word lengths, and therefore does not prove
G2 or Collatz.

## 7. High-height chronological rank certificates

**FINITE-EXACT.** The optional `--rank` round tests sixty prescribed
valuation cylinders instead of more small starting integers. The words are
repetitions of (2,3,1,1) and (1,1,1,1,2,1,1,1), both returning to zero
modulo 223, together with all-one runs. Their lengths range from 16 to
222. For each word, let n_0 be its least positive 223-divisible source and
test

```
n=n_0+223*2^(S+1)*q,
q in {0,1,223,2^64+223,2^256+223}.
```

Every word is verified by exact integer division, every complete node
list is injective, and each prescribed 223-return boundary is checked.
The largest node has 617 bits. Every paired matrix, and in this selected
universe every first-slot matrix, has full rank.

To avoid unproved factorizations of large integers, the script constructs
an exact gcd-free basis of all 3m_j and 3m_j+1. It repeatedly splits bases
b,c at gcd(b,c), preserving their integer exponent rows, and merges equal
bases. The terminal bases are pairwise coprime. Each may be composite,
but every genuine prime-valuation row is a nonzero scalar multiple of its
basis row, so their rational row spaces coincide. The script reconstructs
every input integer exactly from the final bases, then computes rank
modulo the prime 1000000007 while retaining the two coordinate slots.
Full modular rank exhibits a nonzero integer minor and proves full
rational rank; a deficient modular result alone would not prove deficiency.
This minor belongs to the gcd-basis matrix: its rational equivalence to
the genuine prime-valuation matrix is what transfers the certificate.
Equality of those matrices' ranks modulo the auxiliary prime is not
needed and is not asserted.

Independent small controls use actual trial factorizations. The generic
eight-node hostile from the inherited G2 note has paired rank 7 and
first-slot rank 6. The injective orbit segment from 27 through its first
arrival at 1 has 41 transitions: paired rank 41, first-slot rank 36. Thus
the second slot does carry essential chronological information, even
though the high-height cylinder sample already succeeds in the first.

## Current consequence and next question

The successful transfer is precise: finite-field additive geometry can
produce exact forbidden Collatz coset edges and universal, finite-budget
return statements. The hostile return shows why this information alone
does not produce descent. The next chronological question is whether
prime-return spacing and exact valuations impose paired multiplicative
rank constraints stronger than generic unit-equation counting. That
question remains OPEN. The fixed-word finite-exception theorem is a
precise advance in the direction of G2, while exposing the missing
uniformity. A finite matrix statistic is not being presented as a
Collatz proof.
