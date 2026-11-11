# An extension-divisibility filtration that always terminates

**11 November 2026. PROVED P13 below; FINITE-EXACT census through seven
vertices.** This refines the qualifying-path compiler by a computable
arithmetic invariant. Gcd one is not a realization theorem for any
specified target. Chan-Pak remains OPEN. No priority claim.

## P13. Definition, composition, efficient computation, and terminal bound

For a tournament T on n>=1 vertices and k>=0, define

    C_k(T) = gcd { H(U): U extends T by k named new vertices }.

Every orientation of all old-new and new-new pairs is permitted. Naming
the new vertices duplicates some isomorphism classes but does not change
the gcd. Relabelling T or reversing every arc preserves C_k, so these are
ordinary isomorphism and converse invariants.

**Divisibility filtration.**

    C_0(T)=H(T), C_(k+1)(T) divides C_k(T), C_k(T) is odd.

Appending a universal sink to any k-vertex extension preserves its path
count, which proves the divisibility direction. Oddness is the classical
Hamiltonian-path parity theorem, also immediate from the current proved
[THM-002 / Odd-Cycle Collection Formula](../../../01-canon/theorems/THM-002-ocf.md).
More precisely there is the exact recursive composition law

    C_(k+1)(T) = gcd_{U a one-vertex extension of T} C_k(U).          (1)

Each larger extension is obtained by selecting its first named new
vertex and then the other k; all choices are allowed on both sides.

**Polynomially many samples for fixed k.** Let r=kn+binom(k,2) be the
number of free orientation bits, and let D=min(2k,r). The exact
qualifying-path polynomial F_(T,k) has degree at most2k, as proved in
[the extension compiler](../odd_c_seed_symmetry_20261110/EXTENSIONS.md).
Use any fixed orientation-bit convention. Then

    C_k(T) = gcd { F_(T,k)(1_S): S subset [r], |S|<=D }.             (2)

Thus at most sum_(j=0)^D binom(r,j) actual counts suffice, instead of
2^r. For fixed k this is polynomial in n; computing an individual
Hamiltonian count still has its own cost. This is an efficiency result
for the **common divisor**, not for the full support or histogram.

To prove (2), write the unique multilinear integer polynomial
F(x)=sum_A c_A product_(i in A) x_i, with c_A=0 for |A|>D. Finite
differences give c_A=sum_(B subset A)(-1)^(|A|-|B|)F(1_B). Conversely
F(1_S)=sum_(A subset S)c_A. The gcd of all values therefore equals the
gcd of all coefficients and already equals the gcd of the listed small
supports. This also shows independence from the chosen bit gauge.

**Universal termination.** For every n>=1,

    C_(n-1)(T)=1.                                                   (3)

The n=1 case has H=1. For n>=2 take k=n-1, label old vertices
v_1,...,v_n and new vertices w_1,...,w_(n-1), and consider the alternating
undirected path

    v_1, w_1, v_2, w_2, ..., w_(n-1), v_n.

Let E be its2n-2 old-new edges, viewed as free orientation variables.
The coefficient of the degree2n-2 monomial product_(e in E)x_e in F
is **+2 or-2**. Indeed any contributing qualifying permutation must
prescribe every edge in E; it has exactly2n-2 adjacencies total, so it
must be one of the two orders along this same path. Both are qualifying,
since there are no consecutive old vertices and no old arcs to check.
Reversing the path reverses all2n-2 demands. Their top-degree signs agree
because2n-2 is even, giving coefficient twice one sign.

By (2) the common divisor divides that coefficient2. It is odd, so it
equals1. Appending sinks gives C_k=1 for every k>=n-1. This is a universal
finite bound on a particular arithmetic obstruction, not an assertion
that all sufficiently large odd targets occur as extension counts.
Nor does gcd1 guarantee one count coprime to every prime dividing H(T)
simultaneously: a set such as{3,5} already illustrates that logical gap.
It does ensure, for each fixed prime, that some extension count escapes
divisibility by that prime.

Define the **divisibility lifetime**

    L(T)=min{k>=0:C_k(T)=1};  0<=L(T)<=n-1.

It is zero exactly for transitive tournaments, by THM-002: a nontransitive
tournament has a directed triangle and hence H>=3. This creates a
well-defined family of size sequences, for example
D_(n,k,d)=#{ordinary n-vertex T:C_k(T)=d} and
E_(n,l)=#{ordinary n-vertex T:L(T)=l}. Their sums over the last parameter
are A000568(n). This is a partition of that census, not a new formula
for the total census.

## Why the representation of coefficients matters

The gcd of coefficients in the **ordinary multilinear polynomial**
equals the gcd of values. The gcd of the positive counts grouped by
orientation cylinders can be smaller: overlapping cylinders introduce
integer cancellations when expanded. P11's endpoint-group theorem is
a sufficient divisibility certificate, not an iff classification.

The smallest nonunit C_1 example has six vertices and H=45. It is a
vertex deletion of Paley7. Its start/end and individual directed
two-path-cover counts have gcd1, yet its actual one-vertex extension
values have gcd3. Hence local cylinder divisibility alone misses this
first example. The degree-two finite-difference test detects it exactly.

At n=9, R_(3,2) has

    C_0=3159, C_1=3, C_2=1.

Add vertex9 pointing to{0,1,2,3}; the ten-vertex child is asymmetric and
has H15033. Add vertex10 pointing just to{0,1}, with all other old
vertices pointing to10; the resulting eleven-vertex tournament has
H51085. Since gcd(15033,51085)=1, and appending a sink to the first
child supplies15033 as a two-vertex count, these two extensions certify
C_2=1. Independent direct path traversal confirms51085.

## Complete finite census

`content_checks.py` constructs every ordinary class through seven
vertices by canonical one-vertex augmentation. Canonicalization sorts
outdegree cells and minimizes the adjacency code over all permutations
within those cells. Induction gives completeness. The independent labelled
mass sum n!/|Aut(T)| equals2^binom(n,2) at each order, and the class
counts reproduce A000568. No converse quotient is taken.

| n | ordinary classes | C_1=1 | C_1=3 | C_1 other |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 | 0 |
| 2 | 1 | 1 | 0 | 0 |
| 3 | 2 | 2 | 0 | 0 |
| 4 | 4 | 4 | 0 | 0 |
| 5 | 12 | 12 | 0 | 0 |
| 6 | 56 | 55 | 1 | 0 |
| 7 | 456 | 452 | 4 | 0 |

For every class, the gcd from masks of weight at most2 matches the gcd
of all2^n attachments. At order7 the four nonunit classes have H values
9,45,45,189. Adjacency rows and canonical codes are retained. All of
these classes reach C_2=1, as certified using (1) and actual first-child
counts. Thus lifetime distributions through seven have only levels0,1,2.

The alternating coefficient in the proof of (3) is also checked by
literal finite differences for every ordinary tournament of orders2..4.
The all-n terminal bound is a proof, not an extrapolation of this census.

