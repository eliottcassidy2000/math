# Complete trinomial rows: an independent path proof and exact two-row signs

**Status: PROVED analytically using the explicitly CITED classical inputs
below, independently audited by the wildcard lane; FINITE-EXACT and
independently audited for the two-row sign bank. The unbounded two-row
separation remains OPEN.**

The [independent audit](overnight3_20260906_moments_root_geometry_audit.md)
accepts the full-support path proof, the ordinary-product and simplicity
source chain, and the complete positive-mass consumer. Its separate checker
imports neither the producer nor SymPy: it reconstructs all 442 rows from
literal charge equations, verifies all 1,015 first-root intervals by exact
Sturm chains and all remainder enclosures, and counts all 2,242 second roots
as simple negative roots. The sign split is 552 compressed-negative values
at lower index 0 and 463 compressed-positive values at lower index -1;
every restored raw Laurent sign is negative.

While this proof was being written, incoming commit `50e2a436c` promoted
[THM-4436, complete factorial-row simple negative roots and trinomial phase
wall](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
Its [source report](nc2_channel_realrooted_overnight_hexagon_sep05.md), read
from `origin/main` at `d2f64b809`, proves the same all-row conclusion by
factorial splitting, Laguerre symbols and a finite hyperbolicity preserver.
The proof here is an independent lattice-path / finite Pólya-frequency /
ordinary Hadamard route to that concurrently proved conclusion. No duplicate
discovery or external-priority claim is made. Its additional content is the
exact prefix-truncation hostile and certified signs of the true doubled row
at all 1,015 first roots in a declared 221-support bank.

## 1. Inheritance and the live board

The closest proved cross-row mechanism is
[THM-4432, two first channels with both carries](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md),
whose [proof](trinomial_two_channel_empty_core_returns_sep06.md) was read.
It proves a negative second-row value at the unique first root, but does
not cover three or more first channels. The canonical hostile `(-13,1,8)`
has two first channels and five second channels. The corrected near miss is
[MISTAKE-544](../../01-canon/MISTAKES.md), with its
[semigroup repair](synthesis_20260905_moments_trinomial.md): first-row
channels do not generate all later channels. The least-used relevant
sidecar is the **complete finite support of each auxiliary binomial
sequence**, which is larger than the support of its eventual product.
Deleting that sidecar destroys the proof, as Section 5 shows exactly.

The concept board is: semigroup carries; complete binomial path sequences;
finite coefficientwise products; simple scalar roots; cross-row root
separation; phase on the coefficient torus. Carries change the complete
next row but not the path argument. Path positivity supplies individual
root location. The finite product retains every multinomial weight.
Simplicity sharpens the separation target but does not establish it.
The coefficient phase reduces the remaining cancellation locus to one
real ray. Incoming THM-4436 confirms the same conclusion by a different
factorization and supplies another audit route.

Targeted exact-statement, `Hadamard`, `Malo`, `real-rooted`, `interlacing`
and correction searches were performed; the incoming result above is
explicitly incorporated. `CORE-PAPERS.md` was consulted before the
primary sources below.

## 2. A complete binomial-line lemma

Let `p,q>0` and `U,V>=0` be integers. Define a two-sided finite sequence

```text
b_j = binom(U+V+(p-q)j, V+p*j)
      if U-q*j>=0 and V+p*j>=0, and 0 otherwise.
j_min=-floor(V/p),       j_max=floor(U/q).
```

**Lemma.** After translating its **entire** support to start at zero, the
coefficient sequence is Pólya frequency: every minor of its Toeplitz matrix
is nonnegative. Its generating polynomial has only strictly negative real
roots, with multiplicities allowed at this stage.

First prove the matrix assertion before translation. Its entry `b_(j-i)`
counts directed east/north unit lattice paths from

```text
S_i=(-q*i,p*i)       to       T_j=(U-q*j,V+p*j).
```

Every source has level `p*x+q*y=0`, and every sink has level
`p*U+q*V`. Each directed edge strictly increases this level. For any finite
minor all relevant paths lie in a finite region, so the determinant
expansion is a finite signed sum of path families. Swapping the tails of
the first intersecting pair gives the usual sign-reversing involution.
Only vertex-disjoint families survive.

Both source and sink x-coordinates strictly decrease with their indices.
If a pairing has an inversion, follow the two paths as continuous
functions of their common level. Their x-order reverses between the two
boundary levels, so they intersect. East/north unit lattice edges can
intersect only at an integer vertex or along a shared edge, which also
contains a shared vertex. Thus a surviving family uses the increasing
pairing and has positive sign. This proves every minor nonnegative.
The case `U=V=0` has only `b_0=1` and is immediate.

Translation is not an appeal to a truncation theorem. Put

```text
U'=U-q*j_min,        V'=V+p*j_min,       0<=V'<p.
```

The translated sequence is the same path construction with offsets
`U',V'`; its negative indices are all zero. The finite PF characterization
now gives real nonpositive roots, and its positive constant term excludes
zero. The characterization is stated in
[Brändén, *On linear transformations preserving the Pólya frequency
property*, Section 2](https://arxiv.org/pdf/math/0403364).
The parallel-boundary path construction is proved above; planar path
matrices as a method are also used in
[Yu, *Confirming Two Conjectures of Su and Wang*, Section 2](https://arxiv.org/pdf/0901.0385).
Yu's particular increasing-top binomial family is not being substituted
for the decreasing-top family used here.

## 3. The full shifted product, including simplicity

Let integers satisfy

```text
0<A<B, d=B-A, n>=0, N>=0, 0<=t<A,
h=floor(N/B), g=n+N+t.
P(z)=sum_(j=0)^h g! z^j/[(n+d*j)!(N-B*j)!(t+A*j)!].       (1)
```

No coprimality, charge positivity, or minimal-return hypothesis is needed
in this factorial statement. If `h=0`, `P` is a positive constant. Otherwise
it has exactly `h` distinct, strictly negative real roots.

Factor its coefficients exactly as

```text
K_j = alpha_j * beta_j,
alpha_j=binom(g,t+A*j),
beta_j =binom(n+N-A*j,n+d*j).                              (2)
```

The COMPLETE alpha support is `0<=j<=floor((g-t)/A)` and follows from
the lemma with `(p,q,U,V)=(A,A,g-t,t)`. The COMPLETE beta support is
`-L<=j<=h`, where `L=floor(n/d)`, and follows with
`(p,q,U,V)=(d,B,N,n)`. Define the polynomials

```text
F(z)=z^L * sum_(j>=0) alpha_j*z^j,
G(z)=sum_(j=-L)^h beta_j*z^(j+L).
```

Both are real-rooted. `G` has only strictly negative roots; `F` can also
have roots at the origin. Their degrees need not be equal. Their
**ordinary** coefficientwise product is precisely

```text
(F star G)(z)=sum_k [z^k]F * [z^k]G * z^k = z^L P(z).       (3)
```

Indeed the overlap consists of exactly `j=0,...,h`:
`floor((g-t)/A)>=h` follows from `g-t=n+N>=B*h>A*h-1`.
There is no binomial-degree normalization in (3).

For the exact imported strengthening, Brändén's Theorem 3.8 says
`F S G=sum_k k! [z^k]F [z^k]G z^k` is real-rooted when `F` is real-rooted
and `G` has roots of one sign. Section 7, printed page 18, explicitly
identifies ordinary Hadamard composition with `Gamma[F S G]`, where
`Gamma` divides coefficient `k` by `k!`, and states that `Gamma` gives
simple roots away from the origin. These hypotheses hold here, and
`F S G` is nonzero because its coefficient at `z^L` is
`L! alpha_0 beta_0>0`. Removing the exact factor `z^L` in (3) therefore
gives **simple** real roots of `P`. Positivity of all its coefficients
and of its constant term forces them to be strictly negative.
[Primary source, Theorem 3.8 and Section 7](https://arxiv.org/pdf/math/0403364).

Only finite polynomial composition is used. No general closure claim for
Hadamard products of infinite totally nonnegative Toeplitz matrices is
needed. The missing-support counterexample below also blocks replacing
the complete `G` by the desired prefix of beta.

## 4. Exact all-mass consumer and what the quotient loses

For arbitrary charges `a>=1`, `1<=b<c`, let

```text
f(u)=c0*u^(-a)+c1*u^b+c2*u^c,             c0*c1*c2!=0,
g0=gcd(a+b,a+c), A=(a+b)/g0, B=(a+c)/g0,
z=c0^(B-A)*c1^(-B)*c2^A.
```

At a positive mass `m` with a nonempty return fibre, choose its canonical
positive-charge counts `(N,t)` with `0<=t<A`. They satisfy
`A*N+B*t=a*m/g0`. All nonnegative solutions have counts
`(N-B*j,t+A*j)`, for `0<=j<=floor(N/B)`. Put `n=m-N-t`.
The balance identity `a*n=b*N+c*t>0` gives `n>0`, so the corresponding
negative-charge counts are precisely `n+(B-A)j`, with no additional
truncation. Thus

```text
CT(f^m)=c0^n*c1^N*c2^t * P(z),                          (4)
```

with the complete polynomial (1) and `g=m`. This proves simple negative
root location at **every nonempty mass**. The first such mass is not
always `g0`: for `(-1,1,2)`, `g0=1`, while the first return is at mass 2.
The checker reconstructs its first six scalar constant terms as
`[0,2,3,6,20,35]`. Singleton rows are positive constants and have no
cancellation parameters.

Source: a complete balanced multiplicity fibre with its multinomial
weights. Target: (1) together with its nonzero monomial prefactor.
Map: (2)--(4), with scalar invariant `z`. Preserved predicate: exact
constant-term vanishing, including multiplicity of scalar roots.
Lost information: how the root sets at different masses are positioned
relative to each other. Needed sidecar: cross-mass root ordering or a
signed remainder. Cheapest decisive test: evaluate the true doubled
Laurent row at every first-row root using rational isolation.

In particular, outside the strictly negative real `z` ray **every**
nonempty support moment is nonzero. The first nonzero moment is then the
first support return. The pair `(-a,b)` supplies one at mass
`(a+b)/gcd(a,b)<=a+b<a+c`. This phase wall and its smooth scalar-root
interpretation are already part of concurrent THM-4436; the independent
proof does not assert cross-row coprimality.

## 5. Exact obstruction to the tempting prefix argument

The valid primitive support `(-5,2,9)` has

```text
n=2,A=1,B=2,h=2,r=1,t=0,N=5,g=7,
beta_j=binom(7-j,2+j).
```

Its nonnegative-index prefix is

```text
beta_0+beta_1*z+beta_2*z^2=21+20*z+5*z^2,
discriminant=20^2-4*21*5=-20.
```

So even this genuine three-channel row does **not** have a real-rooted
beta prefix. The complete support is `j=-2,...,2`, giving

```text
G=1+8*z+21*z^2+20*z^3+5*z^4,
F=z^2*(1+z)^7,
F star G=z^2*(21+140*z+105*z^2).
```

The actual first row has discriminant `10780>0`. The failed implication
is preservation of real-rootedness after deleting the lower auxiliary
indices. The strongest survivor is the complete-support product lemma;
the missing coordinate is the shift `L` shared by both factors. This is
an explicit obstruction, not a claim of minimality among all supports.

## 6. The remaining two-row problem and finite evidence

For a collided first support return, write

```text
a=A*(B*h+r)+B*t, 0<=r<B, 0<=t<A, h>=1,
N=B*h+r, n=g0-N-t, d=B-A.
P(z)=sum_(j=0)^h g0! z^j/[(n+d*j)!(N-B*j)!(t+A*j)!].
```

The true doubled **Laurent** row is

```text
Q_raw(z)=sum_(j=-floor(2t/A))^(2h+floor(2r/B))
           (2g0)! z^j/[(2n+d*j)!(2N-B*j)!(2t+A*j)!].       (5)
```

Both endpoint carries remain in (5). At a root of `P`, its sign is
meaningful because every first root is negative. The unbounded statement
`Q_raw(rho)<0 for every P(rho)=0` remains **OPEN** for arbitrary `h>=2`.
It would imply two-rung noncancellation and the sharp trinomial width
bound in the collided stratum. Individual simple negative root location
does not prove it: even `1+z` and `(1+z)*(1+2z)` have simple negative roots
and a common root. Those are a logical hostile, not a trinomial example.

The exact bank fixes `(A,B)` in
`{(1,2),(2,3),(3,5),(4,7),(5,9),(5,12)}`, `h` in `{2,3,4,6,8}`,
and distinct corner residues `(r,t)` in
`{(0,0),(B-1,0),(0,A-1),(B-1,A-1)}`. It takes the first two positive
`n` satisfying `b=A*n-(B-A)t>0` and `gcd(a,g0)=1`, then adds
`(n,A,B,h,r,t)=(2,1,2,2,1,0)` and `(16,5,9,3,6,4)` if absent.
This defines 221 supports without a random sampler.

All first and doubled rows have simple negative roots. Their gcd is 1,
and the raw signs at **all 1,015 first roots are strictly negative**.
The first-degree counts are `44,45,44,44,44` for `2,3,4,6,8`.
The carry-signature counts `(floor(2t/A),floor(2r/B))` are
`(0,0):60`, `(0,1):60`, `(1,0):50`, `(1,1):51`.

The bank records primitive integer row coefficients, isolating rational
intervals, and rational interval-Horner enclosures for `Q mod P`.
Multiplication by `(-1)^floor(2t/A)` restores the raw Laurent sign.
Every first root is covered, so there is no untested aggregate cancellation
inside a sampled support. Very large coefficients can make midpoint
evaluation give a wrong sign even when a root approximation looks precise;
only certified interval bounds are used here.

Independent ordered-word dynamic programming reconstructs six raw rows
for `(-4,1,6)`, `(-5,2,9)` and `(-13,1,8)`, retaining the count of the
third monomial. This path uses charge addition, not the channel formula.
It agrees with all multinomial weights, including both carries in the
last support. The nonsemigroup example above is an additional scalar
normalization hostile.

The next exact target is a cross-row ordering or signed-remainder
mechanism that retains both carries. The network theorem gives each
row's complete root geometry; it provides no automatic coupling of the
two distinct networks. No larger endpoint census is needed to expose
that missing piece.

## 7. Reproduction artifacts

- [Producer](../../04-computation/overnight3_20260906_moments_root_geometry.py).
- [Normal output](../../04-computation/overnight3_20260906_moments_root_geometry.out).
- [Optimized output](../../04-computation/overnight3_20260906_moments_root_geometry_optimized.out).
- [Exact root/sign certificates](../../04-computation/overnight3_20260906_moments_root_geometry_certificates.json).

```text
python 04-computation/overnight3_20260906_moments_root_geometry.py
python -O 04-computation/overnight3_20260906_moments_root_geometry.py
```

Every verification uses an explicit exception-raising `need`, so optimized
Python cannot disable a gate. The program is a finite corroboration and
cross-row probe; the all-parameter proof is Sections 2--4 with its cited
inputs. Source SHA-256:

```text
61d7f81c02d56d0809dbc4081d7961325a43558a6073ebacc702a3ca553f676d
```

Normal and optimized runs pass all **23,980 active gates** and have
byte-identical LF output. The exact certificate file is also explicitly
written with LF line endings. SHA-256:

```text
7b37f94b0ccda0a60d238677c20531514c70440e712db7eb902bc9729472e0d6  normal and optimized output
3768e534d170bae5dba0ac3d9c817ddbef789b2ac1a43519362495a934b5d9eb  exact root/sign certificates
```
