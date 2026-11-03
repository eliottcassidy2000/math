# Forget internal allocation, retain endpoint pairing and path parity

Status: **PROVED closed rational carrier and all-order dimension formula**;
**FINITE-EXACT minimum571 at six, attained by actual words of at most five
births**. All-order minimality remains OPEN.

Priority/integration: concurrent [Session145](../session145_endpoint_fence_ramsey_20261103/carrier/POOLED_PARITY_CARRIER.md)
published the same parity quotient and minimum571 first, using exterior
coefficient contexts. This packet independently recovers them and adds
a small bank of actual S/B words, a literal full-row rank audit, the
uniform-birth variance obstruction, and a one-cut second-moment witness.

The previous endpoint-block carrier kept too much information. At order6
its661 coordinates can be reduced to571, and at order7 its3928 coordinates
can be reduced to2668, while preserving every fixed attachment mean.
The sharp dimensions7,28,126 at orders3,4,5 are unchanged.

## 1. Definition

Keep the same source universe and operations as
[FOREST_RESPONSE_RECURRENCE](../primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md):
all labelled Boolean tournaments; old edges retained; all new arcs fixed
independently of the old orientation. Fixed relabelling and ordinary
converse are also allowed. A label-dependent adaptive operation or native
path complement is NOT included by this statement.

Put z_ij=2a_ij-1=-z_ji. For a linear forest F, traverse each nonsingleton
path from its smaller endpoint to its larger endpoint and define sigma(F)
as the product of directed-edge signs relative to i<j. Its type kappa(F)
retains:

1. the named set S of isolated vertices;
2. the matching of named endpoint pairs;
3. for each pair, the parity of the number of internal vertices.

All remaining labels form the internal-vertex set. Their allocation among
the endpoint pairs is deliberately forgotten. Define the polynomial

    R_kappa(z)=sum_(F: kappa(F)=kappa) sigma(F) product_(e in F) z_e.  (1)

Types have disjoint, nonempty Fourier supports. Their polynomials are
therefore linearly independent on the Boolean cube. The constant type is
the empty forest. A path with b internal vertices has b+1 edges; reversing
it changes its signed factor by(-1)^(b+1). This is exactly why its parity
must be retained.

## 2. The coefficient transport theorem

For every fixed attachment C, its Hamiltonian count polynomial H(C(T))
lies in span_Q{R_kappa}. Take two old forests F,G of the same type. In
every Hamiltonian word contributing to the Fourier coefficient of F,
each F-path is a consecutive block. Keep the collapsed block order,
traversal directions and all external vertices fixed, but substitute the
corresponding G-path for each block.

This is a bijection. Old labels reassigned to different blocks remain
internal; no external edge is incident to them in such a word. Every
interface uses the same endpoint names. Every unused old adjacency still
contributes1/2. F and G have the same total number of forced edges, so
the normalization is unchanged. Reversal multiplies the two corresponding
blocks by the same sign because their internal parities agree. Hence
the two Fourier coefficients, corrected by sigma, are equal.

This proof improves the old replacement lemma: it permits different
internal vertex sets in corresponding blocks, provided their parities and
the total internal-label pool agree. It does not assert that the substituted
word is a directed path in a fixed Boolean source. It compares polynomial
coefficients before that source is evaluated.

## 3. Closure under actual fixed births

The vector q(T)=(R_kappa(T)) has rational linear updates for every fixed
attachment. It is enough to prove this for one appended vertex w. In each
term of a child R polynomial, remove w and substitute its fixed incident
edge signs. The old support is a linear forest. Fix any two old forests
of the same type and apply the same internal-allocation substitution.

Every extension of the first forest to a child forest transports uniquely:
w was either isolated, an endpoint grafted to an old endpoint, or internal
and joining endpoints of two distinct old paths (an isolated old vertex
can be a length-zero path in this description). These choices depend only
on named endpoints/isolates. The child's path parities depend only on the
old path parities and the number of added edges. The substitution retains
these, as well as the signs of every edge incident to w. All coefficient
ratios are therefore the old traversal signs, independent of allocation.
So the pullback belongs to the old R span. Iterate for any fixed attachment.

Ordinary converse multiplies each R by(-1)^(number of forest edges), which
is constant within its type. A fixed relabelling permutes types with the
canonical traversal sign, determined by the paired endpoints and parities.
This proves closure, not merely an upper bound on one particular readout.

For an actual source law, m=E[q] and R=E[qq^t] consequently update as
M m and M R M^t. These are possibly signed coordinate matrices; no
nonnegative construction rule or Boolean realization of an arbitrary PSD
matrix follows.

## 4. Exact dimension formula

Suppose there are j>=1 nonsingleton paths and m internal labels. Choose
2j endpoints, pair them, and select the m internal labels; the rest are
isolates. A parity vector of weight w is feasible exactly when w<=m and
m-w is even: allocate one internal label to each odd path, then allocate
the remaining even number to any one path. Therefore

    d_n = 1 + sum_(j=1)^floor(n/2) [binom(n,2j)(2j)!/(2^j j!)]
                * sum_(m=0)^(n-2j) binom(n-2j,m)
                    * sum_(0<=w<=min(j,m), w=m mod2) binom(j,w).    (2)

The leading1 is the all-singleton type. This gives

    d_0..d_10 = 1,1,2,7,28,126,571,2668,13028,65908,347626.

An equivalent formal exponential generating function is

    D(x)=e^x [ e^(x^2/2) sum_(w,r>=0)
                   x^(3w+2r)/(2^w w! (w+2r)!) - cosh(x)+1 ].       (3)

To derive it, sum the matching count over j>=w using
sum_(j>=w) (x^2/2)^j binom(j,w)/j! = e^(x^2/2)(x^2/2)^w/w!.
Subtract the spurious j=0,m>0 even cases: internal labels cannot exist
without a path. This explains the -cosh(x)+1 correction. All statements
are coefficientwise formal identities; no analytic convergence is needed.

At n=6 only the two-path/four-endpoint types newly collapse. For each of
the45 matchings of four chosen endpoints, the two remaining labels have
four allocations, but only two feasible parity patterns. This removes
45*2=90 coordinates, giving661-90=571. Thus the former candidate that the
endpoint-block bound might be sharp at every size is REFUTED at six.
No previously proved small-order minimum is retracted.

## 5. A strictly positive source-law hostile

On six labels take forests

    F1: (0,4,1) and (2,5,3),
    F2: (0,5,1) and (2,4,3).

They have the same endpoints and two odd internal parities. In the standard
i<j edge-bit order, their masks are10376 and5392, and both canonical
traversal signs are+1. Let chi_i=product_(e in Fi) z_e and define

    mu_i(T)=2^(-15)(1+chi_i(T)/2).                     (4)

These are different, normalized, strictly positive laws. For every fixed
attachment C, E_mu1 H(C(T))=E_mu2 H(C(T)), by Section2. But the finer old
endpoint-block polynomial for path blocks{0,1,4},{2,3,5} has expectation
1/32 under mu1 and0 under mu2 in the y=z/2 normalization. Thus the removed
allocation information is real and is invisible to the stated query class.

The one-vertex cut17 (newborn beats old0 and4) already separates their
second moments. Both means are1323/16, but the difference between their
second moments is405/64, by the independent literal Boolean counter.
Thus agreement for all separate mean queries still does not determine
shared-source second moments.

The failure is in law reconstruction, not deterministic source faithfulness.
The inherited all-size growth theorem still separates distinct individual
tournaments. The81/27 control still separates after a named-port attachment.

## 6. A minimum carrier recovered by short actual words

`carrier_rank.cpp` uses signed path-block contraction to search words in
S and B_0,...,B_5, always referring to the six original labels. It finds
571 independent actual query words of length at most5, recorded verbatim
in `rank_contexts.txt`. A contracted block has two endpoint states; reversal
has weight(-1)^(edge count). A free old-old join has weight1 and a fixed
join incident to a new vertex has weight2 or0. Summing signed block orders
therefore computes2^(N-1) times each canonical old Walsh coefficient for
a query on N vertices. The row scaling is nonzero and does not change rank.

The discovery rank is571 modulo both101 and103. More decisively,
`rank_literal.cpp` reads only these query words and evaluates their actual
Hamiltonian counts on768 distinct, explicitly recorded Boolean old sources.
It uses ordinary subset dynamic programming and imports no forest formula.
The resulting571 by768 INTEGER matrix has rank571 over both prime fields.
Thus its rational rank is at least571. The all-context upper bound proves
equality: the minimum rational linear dimension is exactly571, or570
varying coordinates for normalized source laws. This certificate uses
438528 actual graphs, with no exterior-coefficient interpolation required.

The full S/B word banks through lengths0,1,2,3,4 have ranks1,7,48,203,492
modulo both primes in the discovery run. Those are modular ranks; we do
not infer matching rational upper bounds from their agreement. Five births
are a proved sufficient query depth, not a proved minimum depth here.

## 7. Audit and next decisive test

`carrier.py` independently generates every labelled linear forest at orders
six and seven by cut permutations:1486 and12412 forests. It checks the
dimension formula and all170752 basis pullbacks (2668 child types times
all64 cuts) from seven to six, including zero coefficients and cancellations.
`carrier_literal.cpp` separately counts actual Boolean Hamiltonian paths
on all32768 six-vertex sources under89 specified fixed contexts, checking
the two Fourier coefficients directly. This avoids deriving the test from
the contraction argument itself.

Next: determine whether actual attachment queries span all2668 coordinates
at seven and prove or refute all-order sharpness of (2). Closure under native
path complement or source-dependent re-marking requires a separate audit.
The full spectrum theorem remains OPEN; shrinking an exact mean observer
does not construct a target tournament.
