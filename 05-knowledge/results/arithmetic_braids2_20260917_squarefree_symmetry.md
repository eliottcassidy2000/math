# Squarefree arithmetic, growing Collatz prefixes, and the labelled Fano bridge

**Status:** PROVED squarefree densities, simultaneous-squarefree growth theorem,
squareclass maps and sign-gauge kernel; FINITE-EXACT Fano tournament census;
CITED Clifford periodicity and the tournament odd-cycle formula. No Collatz
convergence or topological implication is asserted.

## Inheritance and the working board

This continuation inherits the exponent-box classification in
[the first divisor note](arithmetic_braids_20260917_divisors.md), the exact
word cylinder in [the Collatz note](arithmetic_braids_20260917_collatz.md),
and the floor identities in
[the floor-reciprocity companion](arithmetic_braids2_20260917_floor_reciprocity.md).
The closest proved mechanism is CRT applied to finitely many labelled prime
conditions. The hostile is an arbitrarily long exponent-one growth prefix.
The corrected near miss is inferring a joint realization from separate
existence theorems. The least-used sidecar is the unsieved cofactor, together
with the XOR product index when arithmetic support becomes a Fano diagram.

| Live concept | Operation and preserved predicate | Lost coordinate / hostile |
|---|---|---|
| Squarefree inputs | sieve away `p^2`; retain reducedness of `Z/nZ` | no orbit-height control |
| Finite Collatz word | one 2-adic cylinder, then CRT | future word and termination |
| Simultaneous squarefreeness | sieve all prefix linear forms | later iterates and stopping time |
| Prime squareclasses | reduce exponent vector modulo two | exponent heights; `2,8,32` collapse |
| Fano/octonion signs | orient labelled pairs using multiplication signs | XOR result label and basis gauge |
| Tournament H | count directed Hamiltonian paths | vertex/edge counts are different observables |

For the forbidden H-values use the current repaired
[THM-1370, H-spectrum omits 7 and 21](../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
with its explicit hypotheses and frozen small spectra. Its `H=63` witness
refutes the old `7*3^k` gap-tower reading. The old corollary in THM-200 must
not override this correction: padding an H=63 tournament with a source also
produces a nonstrong H=63 tournament. MISTAKE-502 records the separate role
of the carry atoms 63 and 343.

## 1. Density `6/pi^2` acquires a concrete role through the floor identities

The floor-reciprocity companion proves, with the user's confirmed floor
interpretation, that the cube and inverse-cube floor equalities each hold
exactly when `n>=2` is squarefree. Their exact defects are

```text
(R_3(n)-1)/2,   R_3(n)=prod_(p^a||n) p^(a-ceil(a/3)).
```

Thus these equalities hold on a set of natural density `6/pi^2`. This is
an actual transfer of a predicate: `R_3(n)=1` iff `Z/nZ` has no nonzero
nilpotent elements iff `n` is squarefree. The product-floor equality instead
detects a field, hence primality; `n=6` separates the two predicates.

For completeness, squarefreeness has indicator

```text
mu(n)^2 = sum_(d^2|n) mu(d).
```

Summing through `X` gives
`X sum_(d<=sqrt(X)) mu(d)/d^2 + O(sqrt(X))`; extending the sum to infinity
adds another `O(sqrt(X))`. The Euler product is
`prod_p(1-p^(-2))=1/zeta(2)=6/pi^2`.
The analogous arithmetic-progression main term is recorded in
[Nunes, *Squarefree numbers in arithmetic progressions*, equation (1.1)](https://arxiv.org/pdf/1402.0684).
Only fixed moduli are used here; none of that paper's difficult uniform
error-term results is needed.

**PROVED row-conditioned densities.** Relative to the indicated progression,

```text
P(squarefree | n=1 mod6) = 9/pi^2,
P(squarefree | n=5 mod6) = 9/pi^2,
P(squarefree | n=3 mod6) = 6/pi^2.                         (SF1)
```

These are limiting frequencies, not probabilities on a uniform infinite
integer. In all odd rows the prime-two square obstruction is absent.
In the two unit rows the prime-three obstruction is also absent, giving
`prod_(p>=5)(1-p^(-2))=9/pi^2`. In the third row, two of its three lifts
modulo nine avoid divisibility by nine; this contributes `2/3`, giving
`6/pi^2`. Signed intervals give the same statement for absolute
squarefreeness; negation exchanges rows 1 and 5 and preserves row 3.

The finite audit through 600,000 has exactly 100,000 inputs in each odd row.
The squarefree counts are respectively `91,169`, `60,789`, `91,211` in rows
`1,3,5`. These counts verify a finite sample only; the proof supplies (SF1).

## 2. Squarefree inputs realize every finite halving word

Write `T(n)=(3n+1)/2^{v_2(3n+1)}` on positive odd integers. Fix any positive
word `(k_1,...,k_L)` and put `K=sum k_i`. The first companion proves that
realizing this word is exactly one odd residue class modulo `2^(K+1)`.
Intersecting it with any row modulo six gives one compatible progression
modulo `3*2^(K+1)`.

Consequently, **for every fixed word**, the squarefree density within that
word's sources in row `r` is precisely the same constant `c_r` as (SF1).
Equivalently, the density of squarefree word sources within the whole row is

```text
c_r * 2^(-K).                                             (SF2)
```

The reason is structural: the word fixes only a power-of-two congruence,
and oddness has already removed the prime-two square obstruction. It imposes
no further squarefree obstruction at any odd prime. The fixed-modulus
squarefree sieve just used proves existence and the density, not merely
CRT compatibility for a finite collection of primes.

Thus even after restricting all starting numbers to those satisfying the
two floor equalities, no finite halving word is forbidden. This statement
does not yet say that its intermediate odd iterates are squarefree. The
following theorem establishes that stronger property for the growth words.

## 3. Arbitrarily long growth survives squarefreeness at every intermediate node

**PROVED.** For every integer `L>=1`, there are infinitely many positive
integers `q` such that all `L+1` integers

```text
n_j(q)=3^j 2^(L+1-j) q-1,             0<=j<=L,             (SF3)
```

are squarefree. Such q have positive natural density. Moreover there are
infinitely many such starts in each row `1,3,5 mod6`.
For every q, these integers satisfy

```text
T(n_j)=n_(j+1),   v_2(3n_j+1)=1,   n_(j+1)>n_j.
```

**Local density formula.** Let `nu_L(p)` count the residues `q mod p^2`
for which at least one integer in (SF3) is divisible by `p^2`. Then

```text
nu_L(2)=0,
nu_L(3)=1,
nu_L(p)=min(L+1, ord_(p^2)(3/2)),       p>=5.              (SF4)
```

Here `3/2` denotes `3*2^(-1)` in the unit group modulo `p^2`. The density is

```text
delta_L = prod_p (1-nu_L(p)/p^2) > 0.                      (SF5)
```

**Proof of the local formula.** All coefficients in (SF3) are even, so no
term is even. At three, only `j=0` has an invertible coefficient; every
later term is `-1 mod3`. At `p>=5`, all coefficients are invertible, and the
excluded roots are the first `L+1` terms of a geometric progression with
ratio `2/3` modulo `p^2`. Their count is exactly (SF4). Every excluded root
is nonzero. In fact it is a unit, so even if the progression exhausts its
cyclic subgroup, `nu_L(p)<=p(p-1)<p^2`. There is no local obstruction.

**Proof that the infinite sieve is valid.** For fixed `L`, impose the
conditions through a cutoff `z`. CRT gives density
`delta_(L,z)=prod_(p<=z)(1-nu_L(p)/p^2)` exactly. Put
`A_L=max_j 3^j2^(L+1-j)=2*3^L`. For `1<=q<=X`, an omitted square divisor
has prime `p<=sqrt(A_L X)`. For each of the at most `L+1` forms and each
prime `p>z`, there are at most `X/p^2+1` relevant q. Therefore the number
missed by the finite sieve is bounded by

```text
(L+1) X sum_(p>z) 1/p^2 + O_L(sqrt(X))
 <= (L+1) X/z + O_L(sqrt(X)).                             (SF6)
```

First let `X` tend to infinity with `z` fixed, then let `z` tend to
infinity. This proves (SF5). The product is positive because every factor
is positive and `sum_p nu_L(p)/p^2 <= (L+1)sum_p p^(-2)<infinity`.
All constants depend on the fixed word length; no uniform density in L
is claimed.

**Row refinement.** Conditioning q on the residue that puts `n_0` in a
chosen row does not alter any local factor for `p>=5`. At three, the
conditional factor is `1` in rows 1 and 5 and `2/3` in row 3. Thus the
three relative q-densities are `P_L,(2/3)P_L,P_L`, where
`P_L=prod_(p>=5)(1-nu_L(p)/p^2)>0`.

There is an exact separation formula behind the root collisions:

```text
2^(j-i)n_j-3^(j-i)n_i = 3^(j-i)-2^(j-i),
gcd(n_i,n_j) divides 3^(j-i)-2^(j-i).                      (SF7)
```

So repeated square obstructions are controlled by a geometric order,
not by an independence assumption about consecutive Collatz iterates.

**Finite controls.** Among `1<=q<=10,000`, the exact counts for all
intermediate nodes squarefree are:

| L | 1 | 2 | 3 | 6 | 10 |
|---|---:|---:|---:|---:|---:|
| count | 7,378 | 6,720 | 6,096 | 4,400 | 2,916 |

For example, at `L=10,q=1` every node in the growing path

```text
2047,3071,4607,6911,10367,15551,23327,34991,52487,78731,118097
```

is squarefree. The full three-row witnesses are in the JSON.

The code also supplies rigorous rational intervals for delta_L. For
`z>L+1`, the union bound on the omitted Euler factors gives

```text
delta_(L,z) * (1-(L+1)/z) <= delta_L <= delta_(L,z).         (SF8)
```

At `z=101,L=10`, the stored exact rational bounds have decimal displays
`0.2656554306...` and `0.2981244277...`. These are bounds for the limiting
density, not confidence intervals or replacements by a truncated product.

**Consequence for the other lanes.** Every node of these paths satisfies
both squarefree-detecting floor equalities while the orbit grows for an
arbitrarily prescribed number of steps. Those equalities and squarefree
status alone cannot force a fixed-length descent. Separately proving that
each word has a terminating completion does **not** establish a squarefree
terminating completion: two nonempty realization sets need not intersect.

## 4. The exact squareclass algebra and what it forgets

The positive rational squareclasses form an F2-vector space with prime
basis: the class of `prod p^a` is the finite vector `(a mod2)_p`.
Multiplication becomes vector addition. Each class has a unique positive
squarefree integer representative; multiplying two representatives means
ordinary multiplication followed by removal of square factors.

For three selected primes p,q,r, this gives an explicit F2^3 with eight
elements. Its seven nonidentity elements are

```text
p, q, r, pq, pr, qr, pqr.
```

The seven triples `{a,b,a*b modulo squares}` are the seven Fano lines;
each of the 21 pairs belongs to exactly one line. This is an exact finite
geometry of multiplication modulo squares.

This object occurs **literally inside the third divisor-balance family**:
for `N=p^2qr` the seven proper nontrivial squarefree divisors are exactly
the seven numbers displayed above. The observed equation `10=7+3` therefore
has a genuine Fano sidecar: its S-term counts the Fano points, and its
U-term counts a chosen prime basis of the three-dimensional binary space.
For squarefree `pqr`, the divisor `pqr` itself is excluded by properness,
so only six of those seven points remain. Squaring p restores that seventh
proper divisor. This describes an exact change to the divisor fibre.

The Fano sidecar does not by itself characterize equality: `p^3qr` has the
same seven squarefree divisors but defect four. The exponent-box count is
still necessary. The choice of three named primes identifies the binary
coordinates; it supplies neither a canonical octonion sign table nor an
orientation gauge. Those are additional declared choices in the next section.

Odd powers preserve squareclass and even powers erase it:
`[n^(2t+1)]=[n]`, `[n^(2t)]=[1]`. But squarefreeness itself is not a
squareclass invariant: `2` and `8` have the same class. Nor is the earlier
divisor equation: `2,8,32` share a class but only the first two satisfy it;
`p^2qr` has the same class as `qr`, which does not satisfy it. The exact
prime-exponent profile is the restoration sidecar.

The divisor cube also indexes a basis of an exterior algebra, giving the
same graded dimension `(1+t)^r`, but its multiplication differs. Repeated
support is **zero** under wedge product, is **cancelled** modulo squares,
and is **retained as exponent two** in ordinary multiplication. A Clifford
algebra instead uses anticommuting generators with prescribed nonzero
squares. These are useful related representations, not interchangeable
algebras. The inherited
[THM-3191, factorial-block exterior Clifford law](../../01-canon/theorems/THM-3191-factorial-block-exterior-clifford-law-and-global-carry-smith-profile.md)
has a specific operator and cubic relation; its use of the word Clifford
supplies no automatic arithmetic transfer here.

## 5. An actual octonion-to-tournament map, with an auditable sidecar

Use the Cayley-Dickson basis indexed by `a in F2^3`, with `e_0=1`. Fix the
convention

```text
(a,b)(c,d)=(ac-d*b*, a* d+c b),   (a,b)*=(a*,-b),
```

as in [Baez, *The Octonions*, section 2.2](https://math.ucr.edu/home/baez/octonions/node5.html).
For distinct nonzero binary labels, `e_a e_b=alpha(a,b)e_(a xor b)`, where
`alpha(a,b)` is `+1` or `-1` and `alpha(b,a)=-alpha(a,b)`.

Define a tournament on the seven imaginary basis labels by

```text
a -> b iff alpha(a,b)=+1.                                  (F1)
```

Each Fano line becomes a directed triangle. Each vertex belongs to three
lines and has exactly one outgoing edge in each, so the tournament is
regular of outdegree three. Here a genuine intrinsic binary relation exists;
unlike the earlier sandwich matrix, tournament analysis is appropriately
typed.

Changing basis signs `e_a -> sigma_a e_a` changes the edge sign by
`sigma_a sigma_b sigma_(a xor b)`. The gauges that change no edge are
exactly the characters `sigma:F2^3 -> {+1,-1}`. Indeed the unchanged-sign
condition is precisely the homomorphism law; its values on three basis
vectors are arbitrary and determine it everywhere. Hence there are eight
kernel gauges and `2^7/2^3=16` distinct orientations from 128 sign gauges.

**FINITE-EXACT exhaustive census.** The companion constructs the table from
the displayed recursion, then independently counts Hamiltonian paths by
dynamic programming and by all `7!=5040` vertex permutations:

| Universe | Distinct orientations | Directed Hamiltonian path counts H |
|---|---:|---|
| 128 imaginary-basis sign gauges | 16 | all `H=189` |
| Orient each of the seven Fano lines independently | 128 | 16 have `H=189`; 112 have `H=171` |

The sixteen H=189 orientations in the second row are exactly the sign-gauge
orbit in the first. The remaining line orientations are not being declared
octonion multiplication tables.

The full odd-cycle census supplies an additional hostile to premature
compression: the octonion tournament has **14** directed triangles, 42
directed five-cycles and 24 directed seven-cycles. There are seven disjoint
pairs of odd cycles, so its OCF count is `1+2*80+4*7=189`. The seven Fano
lines account for only half of its triangles. They pairwise intersect;
mistaking just those lines for the entire cycle conflict graph gives K7
and the incorrect value `I(K7,2)=15`. The missing non-Fano cycles are
load-bearing, even though the Fano triangles already orient every edge.

There is a further sidecar distinction. The 168 XOR-preserving relabelings,
`GL(3,2)`, split those sixteen orientations into two orbits of eight,
exchanged by reversing every arrow. Under unrestricted vertex permutations
all sixteen are isomorphic, in a labelled orbit of size 240, hence with
automorphism group of order 21. Thus forgetting the Fano/XOR sidecar merges
a real two-way distinction which the labelled arithmetic geometry retains.
These orbit counts are finite-exact, with the maps enumerated explicitly.

**Connection contract.** The source is a chosen labelled octonion basis;
the target is its sign tournament; (F1) preserves the sign of every product
of distinct basis elements. The tournament alone loses the product's output
label, the identity and the diagonal rule. Restoring XOR, `e_0=1` and
`e_a^2=-1` reconstructs the multiplication table. Changing the sign gauge
is declared explicitly. The hostile is forgetting XOR: it merges the two
GL-orbits above. A second hostile is nonassociativity:
`(e_1e_2)e_4=+e_7` while `e_1(e_2e_4)=-e_7` in this convention.

The seven vertices and 21 arcs are therefore connected to Fano arithmetic
by an actual map. They are **not the Hamiltonian counts** seven and 21.
The H-gaps proved in THM-1370 concern the latter observable, while this
specific bridge yields H=189. No proof of those gaps from octonions follows.

## 6. Oddness and period eight: precise connections and boundaries

There are several distinct parity mechanisms:

- Odd powers are equivariant under `x -> -x`; this is the involution used
  in the floor companion to pair residues and isolate zero residues.
- A permutation cycle of length ell has sign `(-1)^(ell-1)`, so an odd
  cycle has positive permutation sign. The tournament odd-cycle formula
  `H(T)=sum alpha_j 2^j` is an exact separate theorem; see
  [THM-002, OCF](../../01-canon/theorems/THM-002-ocf.md) and
  [Irving--Omar, *Revisiting the Redei-Berge symmetric functions*](https://arxiv.org/pdf/2412.10572).
- Multiplication of squareclasses adds prime-exponent parities; the sign
  `(-1)^omega(n)` on squarefree n is a degree grading. It does not determine
  the orientation or Hamiltonian paths of a tournament.

These admit a common language of involutions and gradings. Transporting a
theorem still requires a map between the actual objects and its preserved
predicate, as in the Fano construction above.

For real Clifford algebras with generators squaring to `-1`, the cited
algebraic periodicity is `Cl_(m+8) ~= Cl_m tensor Mat_16(R)`;
[Baez, section 2.3](https://math.ucr.edu/home/baez/octonions/node6.html)
states the convention and table. Octonions give a representation of `Cl_7`
by imaginary left-multiplication operators, but are not themselves a
Clifford algebra: Clifford multiplication is associative. The eight
squareclasses of three primes supply basis labels, not this periodicity
theorem or a Collatz invariant. The missing extra data are the quadratic
form, anticommutation law, and appropriate representation/topological
construction. No Bott implication is claimed.

## Reproduction and stopping boundary

```bash
python 04-computation/experiments/arithmetic_braids2_20260917_squarefree_symmetry.py
```

[Source](../../04-computation/experiments/arithmetic_braids2_20260917_squarefree_symmetry.py)
and [JSON](../../04-computation/experiments/arithmetic_braids2_20260917_squarefree_symmetry.json)
freeze every finite universe, exact rational density interval, multiplication
table and source hash. Controls include direct square testing through
10,000; independent trial-square tests on the first 300 q at each growth
length; direct residue-root enumeration through prime 19; and independent
brute-force Hamiltonian counts for all 128 Fano orientations. Ordinary and
optimized replay are checked separately. The universal density theorem uses
the two-stage limiting argument (SF6), not extrapolation of these counts.
Both runs passed and their LF-normalized JSON files were byte-identical:

```text
source SHA-256: 12f16b6daccc968925b6e8c842bec73e7aea42d1fff6d41da30bfd6d7f017a77
output SHA-256: ebda1beed54f983d4f1d483100eb0082ff8048ff58706622b59bb96985c9856b
```

The most useful new obstruction is now sharp: a complete finite growing
prefix can satisfy the squarefree floor predicates at every node. Further
progress must retain more than that predicate. The Fano construction gives
a separate faithful algebra/combinatorics bridge with a visible sidecar,
while Bott periodicity remains a cited neighboring structure with no
claimed dynamical transfer.
