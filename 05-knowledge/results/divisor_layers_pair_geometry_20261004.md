# Divisor layers, simplex pairs, and two different four-state Fourier maps

2026-10-04 (America/Denver).

**Status:** PROVED elementary complement obstruction, rank/character formulas,
and symmetry statements under the hypotheses below; FINITE-EXACT for the
declared audits. The divisor classification and labelled ten-letter bijection
are inherited. No Collatz transition, convergence, or prime-novelty theorem
is transferred by these maps.

## 1. Inheritance and the objects being compared

[DB1--DB3 in the divisor note](arithmetic_braids_20260917_divisors.md) define
`F(N)` as the number of divisors `1<d<N`, `S(N)` as the squarefree ones among
them, and `U(N)` as the prime ones among them. If
`N=prod_i p_i^(a_i)` has `r` distinct prime factors, then

```text
F = prod_i(a_i+1)-2,
S = 2^r-1-[N is squarefree],
U = r-[N is prime].                                      (1)
```

For `N>1`, the inherited classification is

```text
F=S+U  iff  N=p, p^3, or p^2*q*r,
```

where the displayed prime labels are distinct. The respective triples are
`(0,0,0)`, `(2,1,1)`, `(10,7,3)`. Prime divisors are already included in
the squarefree set: the equation is a count equality, not a partition into
the literal squarefree and prime sets. In particular, `63=3^2*7` gives
`(4,3,2)`, whereas `126=3^2*2*7` gives `(10,7,3)`; see
[sixth-clock branches, section 7](sixth_clock_branches_20261004.md).

The [earlier balance-family note, sections 3.1--3.2](collatz_mod6_20260917_divisor_balance_family.md)
already supplies the Boolean nonsquarefree face, its labelled prime-to-extra
map, and the k-free extension. The
[tournament/clock/divisor synthesis, section 4.3](procgen_tcpc_20261001_tournament_clock_prime_collatz.md)
already separates two different `3,4,3` partitions of the ten divisors:

* prime / squarefree composite / nonsquarefree proper;
* ranks `Omega(d)=1,2,3` in the exponent box.

They exchange `p^2` and `pqr`; only the rank layers are reversed by
`d -> N/d`. For `p^3`, the first partition is `1,0,1`, while the actual
nonempty rank layers have sizes `1,1`.

The [duck decoder, section 2](duck_decoder_20260925.md) already identifies
the ten divisor letters with `Sym^2_set(F_2^2)` and with `E(K_5)`. The
[tetrahedral pair refinement, sections 1--2](tetrahedral_pair_refinement_20261004.md)
identifies distinct pairs with the six tetrahedral edge midpoints and its
Walsh matrix with the character table of `F_2^2`. These are the closest
proved mechanisms. The canonical hostile is calling the transported
three-cycle a divisor-poset automorphism, refuted in
[the duck prime companion](duck_primes_20260925.md). The corrected near
miss is treating the two `3,4,3` partitions as equal. The least-used
sidecars here are the complement involution, the omitted bottom/top
divisors, and which arithmetic exponents a Fourier character measures.

The working board is divisor order, rank layers, pair histograms,
complement symmetry, arithmetic character charge, and endpoint deletion.

## 2. The ten-point map exists, but complement cannot be geometric

Assign the distinct primes `p,q,r` the nonzero labels `P,Q,R` in
`V=F_2^2`, so `P+Q+R=0`. The inherited map is

| Divisor | Pair in `Sym^2_set(V)` | Point in the tetrahedron |
|---|---|---|
| `p,q,r` | `{0,P},{0,Q},{0,R}` | three midpoints incident to corner `0` |
| `pq,pr,qr` | `{P,Q},{P,R},{Q,R}` | three opposite-face midpoints |
| `pqr` | `{0,0}` | corner `0` |
| `p^2,p^2q,p^2r` | `{P,P},{Q,Q},{R,R}` | the other three corners |

The second map sends `{u,v}` to `(e_u+e_v)/2`, including repeated pairs.
It is a bijection onto the four corners plus six edge midpoints of a
tetrahedron. Thus the seven squarefree divisors become one corner and
all six midpoints; the three nonsquarefree proper divisors become the
remaining corners. The numeric prime labels and the table recover the
original divisors. Divisibility and multiplication are not geometric
incidence operations under this map.

There is also a genuine `3,4,3` **parallel-height** partition of the ten
tetrahedral points. Partition its four corners into two pairs and let
height count the number of endpoints of a pair letter in the second
part. Heights `0,1,2` contain respectively the three letters within the
first part, four cross-part letters, and three within the second part.
Geometrically the end layers are edges with their midpoints; the middle
layer is a parallelogram. A fully specified graded-set bijection is to
sort exponent triples within divisor rank `j` and pair letters within
height `j-1`, then match their positions. The partition and the two
sort orders are additional choices, not arithmetic invariants.

This is a different map from the inherited table. In fact no affine
height function on that table can be constant on divisor rank layers
with distinct rank values: equality on the three rank-one midpoints
forces the heights of corners `P,Q,R` to agree; equality on rank-three
corners `0,Q,R` then forces the remaining corner height to agree as
well. All ten heights would be equal. Thus equal `3,4,3` counts do not
make the inherited squarefree map a geometric rank map. Even replacing
it by the graded-set map cannot fix the complement obstruction below.

**PROVED obstruction, independent of this particular table.** Let `N>1`
and let `v>=3`. Suppose the proper nontrivial divisors of `N` have
`v(v+1)/2` elements. No bijection from these divisors to the degree-two
point set of a simplex with `v` vertices intertwines divisor complement
with an affine symmetry of that point set.

To prove this, any affine set symmetry permutes the extreme corners and
therefore acts on pairs as the induced permutation in `Sym^2_set(V)`.
If it represents complement, that vertex permutation is an involution.
Write `f` for its number of fixed vertices and `t=(v-f)/2` for its
transpositions. The fixed pair letters are exactly the unordered pairs
of fixed vertices, with repetition, and the transposition pairs. Hence

```text
number of fixed pair letters = f(f+1)/2+t
                             = (v+f^2)/2 >= ceil(v/2).       (2)
```

Divisor complement fixes at most one divisor: `sqrt(N)`, when `N` is a
square. For `v>=3`, (2) is at least two. This is a contradiction.
For the ten divisors of `p^2qr`, complement has zero fixed points;
the possible tetrahedral involution counts are `10,4,2`.

The boundary is sharp. The three proper divisors of `p^4` have exponents
`1,2,3`; map them to `0,1/2,1` in a segment. Complement is its reflection,
with one fixed midpoint. The obstruction is to **affine symmetry of the
simplex point set**, not to an arbitrary permutation of its letters.
Under the inherited ten-point map, a concrete failure is
`p <-> pqr`: complement exchanges an edge midpoint and a corner, which
an affine symmetry cannot do.

## 3. A natural all-height layer model retains complement and order

Fix distinct primes `p,q,r` and an integer `k>=2`. Write a divisor of
`N=p^kqr` as `p^a q^b r^c`, with

```text
0<=a<=k, b,c in {0,1}; remove (0,0,0) and (k,1,1).
```

The unimodular coordinate change

```text
(a,b,c) -> (j,b,c) = (a+b+c,b,c),
(j,b,c) -> (j-b-c,b,c)                              (3)
```

turns rank into the first coordinate. Rank one has the three bit pairs
`00,01,10`; each rank `j=2,...,k` has all four bit pairs; rank `k+1`
has `01,10,11`. Thus the proper ranks are exactly

```text
3, [4 repeated k-1 times], 3; total 4k+2.            (4)
```

For each middle layer, `(b,c)` recovers its divisor as
`p^(j-b-c) q^b r^c`. This gives a family of actual Boolean four-point
layers, not a choice of ten unrelated labels. In these coordinates,

```text
complement: (j,b,c) -> (k+2-j,1-b,1-c),
divisibility: (j,b,c) <= (j',b',c') iff
  b<=b', c<=c', and j-b-c <= j'-b'-c'.              (5)
```

Complement is a central affine reflection about
`((k+2)/2,1/2,1/2)`. The replacement point configuration therefore
retains precisely the symmetry that the simplex-pair geometry cannot
retain. The full triple in (3) is essential: rank alone loses the prime
support and does not decide divisibility.

The order automorphism group of this proper divisor poset is exactly
`S_2`, exchanging `q,r`. To prove it, adjoin the unique bottom and top
to recover the full lattice. Its join-irreducible poset consists of the
chain `p,p^2,...,p^k` and the two isolated elements `q,r`. Since `k>=2`,
the chain is distinguished and is fixed pointwise; only the two isolated
elements can be exchanged. Conversely that exchange is an automorphism.
Composing with complement gives all order-reversing bijections, so the
group of order-preserving or order-reversing symmetries is `C_2 x C_2`.
In particular there is no order-three order automorphism.

This extends the rank geometry, not the original squarefree balance:
for all `k>=2`, `S=7,U=3,F=4k+2`, so `F=S+U` holds only at `k=2`.
The different k-free balance is the inherited result cited in section 1.

## 4. The rank-Fourier theorem and its arithmetic meaning

Define the actual exponent charge

```text
chi(d)=(v_q(d) mod2, v_r(d) mod2) in F_2^2.
```

It extends to a multiplicative-to-additive homomorphism from the positive
rational group generated by `p,q,r`. For proper divisors it is just
`(b,c)`. It forgets the `p` exponent; retaining `j` restores it by (3).
For a character label `(s,t) in F_2^2`, define

```text
W_st(z)=sum_(1<d<N, d|N) z^Omega(d) (-1)^(s*b+t*c).
```

**PROVED, all `k>=2`:**

```text
W_st(z)=(1+z+...+z^k)(1+(-1)^s*z)(1+(-1)^t*z)
        -1-(-1)^(s+t)*z^(k+2),

W_00(z)=3z+4(z^2+...+z^k)+3z^(k+1),
W_10(z)=W_01(z)=z-z^(k+1),
W_11(z)=-z-z^(k+1).                                (6)
```

The product enumerates the exponent box; the two subtractions remove
exactly bottom and top. Alternatively every middle layer in (4) is a
full Boolean square, so its three nontrivial character sums vanish.
The two boundary layers supply all nontrivial Fourier coefficients.
Complement gives the useful identity

```text
z^(k+2) W_st(1/z)=(-1)^(s+t) W_st(z).              (7)
```

Thus the individual `q` and `r` modes are antipalindromic, while the
trivial and joint modes are palindromic. At `z=1`, the charge counts
and unnormalized Fourier spectrum, in the order `00,01,10,11`, are

```text
counts:   (k,k+1,k+1,k),
spectrum: (4k+2,0,0,-2).                           (8)
```

The two Klein-four groups here also have different roles. The charge
group has XOR addition. The order/duality group in section 3 acts on
that charge by `(b,c)->(c,b)` and `(b,c)->(1-b,1-c)`; it is not the
four-element translation action. On Fourier labels, prime exchange
swaps `(s,t)`, while complement gives the sign and rank reversal in
(7). These formulas supply the connection without identifying the
actions merely because both groups have four elements.

This is **different from** the inherited unordered-pair XOR map. Its
counts are `(4,2,2,2)` and its spectrum is `(10,2,2,2)`. At `k=2`,
even the sorted count multisets differ, so no relabelling of the four
states can identify these two maps. They use the same character table
on different observables. Moreover the inherited pair XOR is not an
arithmetic homomorphism: take `N=60`, `p=2,q=3,r=5`. It assigns
`4 -> 0`, `3 -> Q`, but `12 -> 0`, although `4*3=12` remains a proper
divisor. The exponent charge does respect that multiplication.

One usable consequence is an exact word-product parity counter. For
length `ell>=1` words of proper divisors of `p^kqr`, whose total charge
is the `q,r` parity of their arithmetic product, the number with
specified charge `(b,c)` is

```text
((4k+2)^ell + (-2)^ell*(-1)^(b+c))/4.              (9)
```

Character diagonalization of XOR convolution proves (9). The arithmetic
product need not divide `N`; its parity charge remains defined. This
counter does not retain the full product or any Collatz trajectory.

The same mechanism extends without a new classification. For
`N=p^k q_1...q_m`, all primes distinct, `k>=1,m>=1`, retain the `m`
Boolean exponents and use `s in F_2^m`. Then

```text
W_s(z)=(1+...+z^k) prod_i(1+(-1)^(s_i)*z)
       -1-(-1)^|s|*z^(k+m),
z^(k+m)W_s(1/z)=(-1)^|s|W_s(z).
```

At `z=1`, the trivial coefficient is `(k+1)2^m-2`; every nontrivial
odd-weight character has coefficient zero, and every nontrivial
even-weight character has coefficient `-2`. These coefficients record
the two deleted antipodal corners. Flipping Boolean coordinates is an
action on the full box, but need not preserve the deleted set or divisor
order. Fourier analysis needs the characters; it does not assert such
flips are order automorphisms.

## 5. Exact verification and the connection boundary

The standalone [script](../../04-computation/experiments/divisor_layers_pair_geometry_20261004.py)
uses integer arithmetic and explicit checks that survive `-O`. It imports
no other experiment. The [saved output](divisor_layers_pair_geometry_20261004.out)
records:

* direct divisor counts and the inherited classification for `2<=N<=500`;
* 496 general rank-polynomial, complement, and integrated-spectrum cases,
  with `1<=k<=8`, `1<=m<=5`, and every character;
* the rank chart and inverse for `2<=k<=20`, including 49,324 pairwise
  order checks and all four rank-character formulas;
* an independent order-matrix backtracking enumeration giving exactly two
  proper-poset automorphisms for each `2<=k<=6`;
* all 1,115 vertex involutions for simplex sizes `1<=v<=8`, checking (2);
* the ten-point bijection, the two inequivalent charge histograms, the
  separate `2+2` geometric height partition and its graded-set bijection,
  the multiplicative hostile `4*3=12`, and exact word convolution for
  lengths one through five.

Reproduce from the repository root:

```text
python -X utf8 -B 04-computation/experiments/divisor_layers_pair_geometry_20261004.py
python -O -X utf8 -B 04-computation/experiments/divisor_layers_pair_geometry_20261004.py
```

Normal and optimized runs agree byte for byte. These finite tests audit
the implementations; the proofs supply the unbounded statements.

The useful bridge is now explicit. The ten-letter map retains a labelled
set, its chosen `7+3` split, and the inherited color permutation; its
sidecar is the divisor decoding table. The exponent/rank map retains
divisibility, complement, and prime exponents; its sidecar is the full
triple `(j,b,c)`. The Fourier projection retains the chosen product
parities; its missing coordinate is rank or exponent height. Equation
(2) is the decisive obstruction to identifying all these structures as
one geometry. No intrinsic tournament relation is introduced by their
matching cardinalities.
