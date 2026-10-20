# Raw Hamiltonian-path spectrum: the known conjecture and finite-core reductions

**Current answer: OPEN.** It is proved that no tournament has exactly 7 or
21 directed Hamiltonian paths. The assertion that these are the only missing
positive odd integers is a published conjecture. This packet does not claim
to settle it. It gives a fresh, small finite proof certificate for the two
exclusions, an all-order logarithmic core bound, exact rational series for
each fixed raw count, and a constructive one-apex recurrence with a decisive
hostile showing why that recurrence alone cannot prove completeness.

## 1. Types, inheritance, and primary provenance

Here `H(T)` is the **raw number of directed Hamiltonian vertex words** in T.
It is not `h(T)=H(T)/|Aut(T)|`, and is not the blue orbit count b(T) or
the black deficit h(T)-b(T). Likewise P_N below is a new raw-count series,
not any earlier strong-center D_beta or orbit-mass G_h.

The closest proved mechanism is ordered-SCC multiplication, followed by
the strong-tournament lower bound. The hostile to scalar completeness is
an insertion whose attachment data have been forgotten: even a source and
a sink marking of the same two-vertex tournament respond differently.
The corrected near miss is the historical assertion that finite odd-count
coverage proves all-order completeness. The useful additional data are the
entire one-apex cut word and the ordered list of nontrivial SCC types.

Primary sources checked on 2026-10-20:

* S. H. Chan and I. Pak, *Computational complexity of counting coincidences*,
  Theoretical Computer Science (2024),
  [DOI 10.1016/j.tcs.2024.114776](https://doi.org/10.1016/j.tcs.2024.114776),
  [author PDF, section 6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).
  Their Conjecture 6.2 is exactly the proposed raw spectrum. They credit
  bof and Gordon Royle's March 2016 question and computation. **CITED OPEN
  PROBLEM**, not a new conjecture of this repository.
* A. H. Busch, *A Note on the Number of Hamiltonian Paths in Strong
  Tournaments*, Electronic Journal of Combinatorics 13 (2006), N3,
  [full primary PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/).
  Theorem 1 and Corollary 1 give the exact strong minimum used below.
  The four-page proof was read, including its overlapping-strong-subtournament
  and order-preserving path-merge lemmas. **CITED THEOREM**.
* [Royle's original MathOverflow answer](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament),
  updated January 2022, reports that all positive odd values below 80557
  except 7 and 21 occur by order 11. This is **CITED COMPUTATIONAL EVIDENCE**;
  the underlying order-eleven census was not rerun here. It supplies no
  uniform construction for all larger counts.

The later primary article *About the number of directed Hamiltonian paths
in special tournaments*, [DOI 10.1016/j.dam.2026.07.017](https://doi.org/10.1016/j.dam.2026.07.017),
describes formulas for selected families. Its accessible abstract and
introduction do not claim spectrum completeness; they must not be promoted
to such a result. The current bounded literature search found no proof
of the conjecture. This is a search result, not a theorem that no proof
could exist elsewhere.

Relevant current canon is [THM-1370, all-order omissions](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
and [THM-4094, two-prime-lane reduction](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
Both explicitly leave completeness open. The stateful cut-response
obstruction in [THM-4118](../../../01-canon/theorems/THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md)
also prevents replacing actual attachment states by a response gcd.

### Historical proof audit

The status banner of [THM-115](../../../01-canon/theorems/THM-115-h21-permanent-gap.md)
records a later all-order cutoff, while its body retains older conjectural
language. Use the repaired theorem, not that historical body, for status.
Its pancyclic cutoff is weaker than Busch's and does not imply completeness.

Two older local proof sentences must not be inherited:

* [THM-029](../../../01-canon/theorems/THM-029-alpha1-gap.md) formerly claimed that three
  pairwise-intersecting triples must have a common vertex. As a set claim
  this is false: `{1,2,3}`, `{1,4,5}`, `{2,4,6}` are a counterexample. Additional
  tournament structure would have to do the work; the asserted pigeonhole
  implication does not. The revised local strong-union proof in
  [THM-343](../../../01-canon/theorems/THM-343-H7-impossible.md) avoids it.
* [THM-200](../../../01-canon/theorems/THM-200-h7-impossibility.md), Step 4,
  formerly called the conflict graph P3 and evaluated its independence polynomial
  at 2 as 15. In fact `I(P3,2)=11`; the described H=15 configuration is
  `K2 disjoint-union K1`, whose value is 15, as THM-343 later records.

These do not refute the exclusions. The independent proof below avoids
both erroneous arguments.

## 2. The two omissions need only a small certified base

For a strong tournament C of order n>=3, Busch gives

```text
H(C) >= b(n),
b(3q)   = 3*5^(q-1),
b(3q+1) = 5^q,
b(3q+2) = 9*5^(q-1).                                  (1)
```

In particular the minima at orders 3,4,5,6,7 are 3,5,9,15,25, and
`H(C)>=5^((n-1)/3)`.

Every Hamiltonian path of an ordered join concatenates one path in each
SCC, so `H(T)=product H(C_i)`. A value of 7 therefore requires a strong
component with raw count 7. Equation (1) restricts it to order at most
four, whose strong counts are 3 and 5. A value of 21 requires a strong
count 7 or 21. The former is excluded; the latter has order at most six.
The complete small control below shows that 21 is absent there. Thus
`H(T) notin {7,21}` at every order.

Only **1,100 standard-Hamiltonian-path presentations** of orders 1–6 are
needed, not all labelled tournaments. Force the consecutive arcs forward
and enumerate all remaining arcs. Every tournament has a Hamiltonian path,
so relabelling one such path to the standard order puts it in this universe.
Raw path counts are then recomputed by subset DP, with no quotient division.
For type counts, minimize the adjacency word over all actual Hamiltonian
words: this gives a canonical form inside the same universe.

| order | presentations | ordinary isomorphism classes | strong raw H values |
|---:|---:|---:|---|
|3|2|2|3|
|4|8|4|5|
|5|64|12|9,11,13,15|
|6|1024|56|15,17,19,23,25,27,29,31,33,37,41,43,45|

This is a **CITED + FINITE-EXACT proof** of the exclusions. The two infinite
attainability lanes in THM-4094 remain genuine open obligations; finite
prefixes, ordered products, and monotone padding do not settle them.

## 3. Every raw-count witness has a logarithmic core

**Theorem (PROVED relative to Busch).** Delete every singleton SCC from
a tournament T with `H(T)=N>1`, retaining the order and actual tournaments
of all other SCCs. The resulting core C satisfies

```text
H(C)=N,                 |C| <= floor(3 log_3 N).       (2)
```

The sharper individual bound for a strong component of raw count d is
`|C_i|<=1+3 log_5 d`, or the exact cutoff obtained by inverting (1).
Equality in `3^|C|<=N^3` occurs exactly when every component of C is C3.

**Proof.** For a strong component of size s=3, H=3. For s>=4, (1) gives
`H^3>=5^(s-1)>3^s`; the strict inequality holds at s=4 and strengthens
with s. Multiply `H(C_i)^3>=3^|C_i|` over the nontrivial SCCs. Singleton
components contribute one to H and are precisely the discarded data.
Equality forces all nontrivial components to have size three. Conversely
a chain of k triangles has core size 3k and count 3^k. QED.

This is an effective finite test for any proposed N: a witness, if one
exists, has at most the order in (2). For binary input length L, enumerating
that finite universe and counting paths takes `2^O(L^2)` time. No such
large enumeration is run here. Alternatively enumerate only strong types
at the sharper divisor-dependent cutoffs and use multiplicative closure
on the divisors of N. The reduction does not provide a uniform positive
construction for every N, and does not put verification in polynomial
time in L.

In Chan--Pak's Definition 5.7, a counting function is *concise* when each
attained value has a witness of polylogarithmic size. Thus H is concise:
the core has O(log N) vertices and O(log^2 N) adjacency bits, with N=1
handled separately. This is a direct consequence of the classical bound,
not a claimed new literature theorem. Witness compression is unconditional;
almost-completeness remains their Conjecture 6.2.

## 4. Every fixed raw-count class series is rational

Let P_d(z) count ordinary isomorphism classes of **nontrivial strong**
tournaments with raw count d, by vertex number. Equation (1) makes every
P_d a finite polynomial, whether or not d is attainable. Let

```text
F_N(z)=sum_[T]:H(T)=N z^|T|,
```

including the empty tournament only when N=1. Then `F_1(z)=1/(1-z)` and
for odd N>1,

```text
F_N(z) = (1/(1-z)) sum_(d|N, d>=3) P_d(z) F_(N/d)(z). (3)
```

Equivalently,

```text
F_N(z)=sum_(k>=1) (1-z)^(-k-1)
           sum_(d1*...*dk=N, di>=3) product_i P_di(z). (4)
```

**Proof.** There are k+1 independent singleton gaps around an ordered
word of k nontrivial strong components. Isomorphism preserves the unique
SCC order, so this decomposition neither loses nor duplicates ordinary
tournament classes. The raw count is the product of component counts.
At most Omega(N) factors are possible. This proves both formulas. QED.

Let kappa(N) be the greatest number of nontrivial SCCs in a witness, when
N is attainable. The pole at z=1 has **exact order kappa(N)+1**: all
highest-order terms have positive numerator at one, so none cancel.
Consequently the number of ordinary n-vertex types with raw H=N is
eventually a polynomial of degree kappa(N). Its leading coefficient is

```text
(1/kappa(N)!) sum_(d1*...*dk=N, k=kappa(N)) product_i P_di(1). (5)
```

These are ordinary classes, not converse-merged classes, and not the
earlier orbit-deficit series. A zero F_N decides that N is absent; knowing
each F_N is rational does not decide that all but two are nonzero.

The exact small strong polynomials, complete at all orders by (1), are

```text
P3=z^3, P5=z^4, P7=0, P9=2z^5,
P11=P13=z^5, P15=2z^5+2z^6,
P17=P19=2z^6, P21=0, P23=4z^6.                       (6)
```

For example,

```text
F9  = 2z^5/(1-z)^2 + z^6/(1-z)^3,
F15 = (2z^5+2z^6)/(1-z)^2 + 2z^7/(1-z)^3.            (7)
```

The square and cubic denominators here count singleton gaps around one
or two raw-count factors. They must not be confused with the Narayana
denominators in the self-converse orbit-deficit problem.

One explicit all-order consequence is useful. For fixed nonnegative a,b
with a+b>0, the number of ordinary n-vertex types with `H=3^a 5^b` obeys

```text
[z^n]F_(3^a 5^b)(z) ~ n^(a+b)/(a! b!) as n->infinity. (7a)
```

Indeed a+b is the maximum possible SCC-factor count, and attaining it
forces a triangles and b copies of the unique strong four-vertex type.
There are `binom(a+b,a)` orders for those factors. All other cores have
fewer gaps and give lower-degree terms. More generally, for an attained
prime-factor inventory whose primes all occur as strong counts, (5)
has leading coefficient `product_p P_p(1)^e_p/e_p!`. These are fixed-raw-H
class asymptotics and are not an evaluation of the earlier M_k notation.

## 5. A constructive recurrence, and its first missing attainable prime

Adjoin a marked vertex v to the transitive tournament `0->1->...->m-1`.
Let `s_i=1` mean `v->i` and `s_i=0` mean `i->v`. This is an actual
tournament construction. Its raw Hamiltonian count is

```text
H(s)=1+sum_(i<j, s_i=1,s_j=0) w(j-i),
w(1)=2,                  w(d)=2^(d-1) for d>=2.       (8)
```

For a proof, a Hamiltonian path has an increasing prefix, then v, then
an increasing suffix. If the prefix ends at a and the suffix begins at b,
the case a<b forces b=a+1; the case b<a leaves `2^(a-b-1)` choices for
the intermediate vertices. Paths with v at an endpoint and the adjacent
`0,1` transitions combine to `1+#adjacent(1,0)`. Adding the split paths
gives (8), including the extra unit for an adjacent `1,0` pair.

There is a constant-state integer recurrence. Start with `(H,A,l)=(1,0,0)`.
For each next bit s, apply

```text
s=0: (H,A,l) -> (H+A+l, 2A,   0),
s=1: (H,A,l) -> (H,     2A+1, 1).                    (9)
```

Here A is the binary integer encoded by the processed word and l its last
bit. Its parity equals l, so the H increment is even. Formula (9) follows
directly from (8): a new zero receives one contribution from every earlier
one, with an additional unit for an immediately preceding one.

This is not a recurrence in H alone. The one-bit words 0 and 1 both give
H=1, but appending zero gives H=1 and H=3 respectively. The retained
attachment word, or its state `(A,l)`, is essential. A useful explicit
family is

```text
H(1^a 0^b)=(2^a-1)(2^b-1)+2,             a,b>=1.      (10)
```

However **no one-apex/transitive-core tournament has H=13**. Delete leading
zeros and trailing ones from s: these are singleton SCCs, so H is unchanged.
If an inversion remains, the interval from the first one to the last zero,
together with v, is the unique nontrivial strong component. If its order
is at least six, (1) gives H>=15. At order at most five, (8) gives only
`1,3,5,9,11`, never 13. If no inversion remains, the tournament is transitive.

Yet the actual five-vertex tournament with outneighbor bit masks

```text
(10,4,9,18,7)
```

has H=13. Its single-vertex deletions are all nontransitive, and a
two-vertex deletion makes it transitive. Thus feedback-vertex number two
is necessary already at this attainable prime. Because 13 is prime,
even the ordered-SCC product closure of the entire one-apex family still
cannot realize it. This is a sharp stopping reason for that proposed
binary scalar construction, not an additional gap in the full spectrum.

## 6. Reproduction and surviving target

```text
python 05-knowledge/results/spectrum_basins_20261020e/spectrum_probe.py
python -O 05-knowledge/results/spectrum_basins_20261020e/spectrum_probe.py
```

The [checker](spectrum_probe.py) uses 1,100 standard-path presentations,
independent subset-DP counts, canonical Hamiltonian-word forms, and the
complete class-count controls through order six. It separately builds
all 511 one-apex cuts through transitive base order eight and compares
literal counts with both (8) and (9). It also checks the unique strong
interval, the H=13 hostile, 64 instances of (10), SCC multiplication, and
the integer form `3^core_order<=H^3` of the logarithmic bound.
It also compares (3) against the complete ordinary class census through
order six for every odd raw count at most 23 and retains its coefficients
through order twelve. Normal and optimized runs agree byte-for-byte in
their JSON receipt.

The proof-level next target remains a stateful construction supplying the
two infinite strong-count lanes from THM-4094, or a new certified missing
odd value. The present recurrence records the necessary attachment state
but cannot meet that target with only one feedback vertex. No bound on
further finite coverage would by itself discharge the all-order conjecture.
