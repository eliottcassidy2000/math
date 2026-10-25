# Finite-horizon moment transfer across exact clone excess

2026-10-25. **PROVED:** an exact linear aggregation for the signed
pair-free-core transfer, with rank min(S,h) on a fixed defect fibre with
remaining total degree S and h distinct cycle-count values. It removes
cycle count and final actual degree from the key, using a finite moment
vector as the weight. An adaptive sparse/moment implementation has no
more scalar weight coordinates than the old f-state transfer at every
stage, with strict measured reductions. **FINITE-EXACT:** four-path
coefficient controls through18. No bit-cost or runtime improvement is
claimed; the dense-only implementation remains a cost hostile.

## 1. Inheritance and the precise change of observable

The closest mechanism is the
[ordinary fixed-degree defect quotient](../../ports_regular_actions_hn_20261025/counting/README.md).
Its [excess extension](../../ports_regular_actions_hn_20261025/counting/EXCESS_EXTENSION.md)
proved that, without fixing actual degree, no scalar normalization can
merge distinct affordable f profiles while preserving every individual
core continuation. Fixing actual and excess degrees separately repairs
that scalar quotient. The
[clone-cut involution](../../ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md)
proves why complete signed inflation fibres cancel, with literal clone
mergers retained before quotienting their cut sets.

The new target is different from merging two individual histories: preserve
the aggregate future contribution of a weighted bank of histories. Signed
weights and linear combinations are retained. The live concepts are fixed
total degree / defect response / cycle-count moments / signed clone excess /
finite continuation horizon. Canonical hostile: a moment bank invisible
through horizon2 is visible at3. Corrected near miss: fewer vector keys
need not mean less arithmetic or memory.

## 2. A bank needs only S moments

Fix the target total degree N and process odd cycle lengths in decreasing
order. At a stage ell let the committed degree be u, including previously
assigned substitution excess, and put S=N-u. For a prefix let m be its
number of actual cycles and

    f(d)=sum_(past actual c) gcd(c,d),
    D(d)=d*m-f(d),       D(1)=0.

Retain affordable odd d<=min(ell,S). Group prefixes having the same S,D.
If their accumulated signed rational weights are w_H, define

    W_r=sum_H w_H * 2^(r*m_H),       1<=r<=S.         (M1)

For S=0 retain the single terminal value W_0=sum_H w_H. For S>0 no
W_0 channel is needed: a positive remaining substitution degree cannot
be completed with zero future actual cycles, because excess at each
length requires a positive actual multiplicity there.

Consider any decorated future consisting of actual odd cycles mu, with
actual degree a and excess degree S-a. Its interaction with H is

    2^[a*m_H-sum_(d in mu)D(d)].

All intrinsic cycle factors and signed excess coefficients depend only on
that future. Summing over H replaces the interaction by

    2^[-sum_(d in mu)D(d)] * W_a.                    (M2)

Consequently the vector (M1) determines every individual future's
aggregate contribution. It also determines their final signed sum. No
actual-degree target M, excess target E, or cycle count m belongs in the
key; the key is (S,D), and its weight is the moment vector.

This does not contradict the earlier scalar obstruction. A vector carries
more than one scalar and may annihilate a signed combination of histories
without identifying any two single histories as interchangeable.

## 3. The closed batch transfer

At odd length ell choose k actual cycles and e excess units. Define

    beta(k,e)=(-1)^e binom(k+e-1,e)  for k>0,
    beta(0,0)=1, beta(0,e)=0 for e>0,
    A_ell(k)=k*(ell-1)/2+ell*binom(k,2),
    S'=S-ell*(k+e).

On affordable smaller odd d put

    D'(d)=D(d)+k*(d-gcd(ell,d)).                     (M3)

For every destination channel r=1,...,S', or r=0 if S'=0, add

    W'_r += beta(k,e)/(ell^k*k!)
            * 2^[A_ell(k)-k*D(ell)+r*k]
            * W_(ell*k+r).                         (M4)

All choices have k,e>=0 and ell*(k+e)<=S, with the k=0 convention above.
If ell>S, only k=e=0 is possible and the unavailable D(ell) is not read.
The source channel in (M4) lies within1..S whenever S>0, since
ell*k+r<=ell*k+S'=S-ell*e. The terminal S=0 state simply carries W_0.

**Proof.** The ordinary source transition multiplies w by

    beta(k,e)/(ell^k*k!) * 2^[A_ell(k)+k*f(ell)],

and changes m to m+k. Multiplying by2^[r*(m+k)] and using
f(ell)=ell*m-D(ell) gives exactly the scalar factor and source moment
in (M4). Summing the equality over the bank proves the transition.
Equation (M3) follows directly from the definition of D. QED.

Initialize the empty prefix at S=N,D=0 with every W_r=1. Descend through
all odd lengths. At ell=1 force k+e=S, then collect terminal W_0. The
result is the exact number of pair-free tournament classes of order N,
because the inherited signed substitution and the clone-cut identity
remain unchanged. At ell=1 the entire remaining transfer can also be
written as the scalar boundary functional

    sum_(k=1)^S beta(k,S-k)*2^binom(k,2)/k! * W_k,   S>0.  (M5)

Weights in (M4) may be signed and rational, and its displayed exponent
may be negative. The reference implementation uses exact Fractions.
No unchecked integer division or positive-exponent assumption is made.

The ordinary unrestricted A000568 transfer is the no-excess restriction:
only the future actual degree r=S is observable, so its scalar D-gauge
is precisely the one-channel version. The new contribution is extending
that aggregation to exact clone excess without choosing separate final
actual-degree targets. It does not improve the already proved ordinary
one-channel state graph by itself.

## 4. Exact rank and strict preterminal compression

Suppose a fixed (S,D) bank has h distinct admissible cycle counts
m_1,...,m_h, with arbitrary aggregate weights at those counts. The map to
moments is the matrix

    V_(r,j)=(2^m_j)^r,       1<=r<=S.                (M6)

Its rank over Q is exactly min(S,h): multiplying column j by2^(-m_j)
leaves an ordinary Vandermonde matrix with distinct nonzero abscissae.
This is also the minimum dimension for a linear representation preserving
all individual future aggregates. For every r=1,...,S, the future of r
actual singletons and S-r singleton excess units exists and has nonzero
beta coefficient. Since D(1)=0, this future reads exactly W_r up to a
known nonzero scalar. Therefore none of the independent moment directions
can be removed for that observable.

There are genuine preterminal fibres with h>S. Take S>=ell>=3. At stage ell choose an
odd L>ell divisible by every future odd length1,3,...,ell. Fix j and
consider, at processed length L, k=1,...,j actual cycles together with
j-k excess units. These prefixes all have committed degree L*j and

    m=k,       D(d)=k*(d-gcd(L,d))=0.

They are legal source terms of the signed substitution, with nonzero
weights. Examples are ell=3,L=9; ell=5,L=15; ell=7,L=105. For a target
N=L*j+S they lie in the same (S,D) fibre before the remaining nontrivial
lengths are processed. Restricted to those source terms, the h=j scalar
response states have linear rank
min(S,j). For j>S this is a strict dimension reduction, and j/S is
unbounded while the future horizon is fixed.

The rank statement concerns the linear operator on weighted banks. It is
not a claim that arbitrary chosen bank coefficients are the particular
prefix weights produced by one unmodified counting run. The real family
above supplies the admissible histories; the linear map describes how
all weighted contributions on those histories propagate.

## 5. The horizon is load-bearing

At the cycle counts m=1,2,3, the formal signed bank with weights8,-6,1
has

    W_1=0, W_2=0, W_3=192.                          (M7)

It is invisible to every continuation of total degree at most2 within a
fixed D fibre. At total degree3 the three-singleton future detects it.
Thus a moment quotient cannot be reused at a larger remaining degree
without restoring additional moments. The transfer only moves to a
smaller horizon, which is why (M4) is safe.

Dropping the entire vector to its final scalar is likewise valid only at
the final boundary (M5), after all remaining options have been evaluated.
The existence of the clone-cut involution does not provide missing future
moments or make a low-order zero a global recurrence certificate.

## 6. Adaptive storage is dominated by the old scalar transfer

There is an additional sign invariant that makes adaptive compression
safe without reconstructing cycle counts from moments. At fixed committed
degree u and actual cycle count m, every prefix term has sign

    (-1)^(u-m).                                     (M8)

Indeed its sign is (-1)^(sum e_d), while
u-m=sum(d-1)*k_d+sum d*e_d is congruent to sum e_d modulo2 because all
lengths d are odd. All other factors are strictly positive. Therefore
terms with the same u,m never cancel, even when their cycle histories or
excess assignments differ. In particular the old nonterminal f-state
weights are nonzero: f(1)=m keeps the sign classes separate.

At each (S,D) key store the sparse map m->w_m while its support size h
is at most S. Once h>S, convert to the S moments (M1). A sparse
transition uses the original cycle factor; a dense transition uses (M4).
When merging a dense contribution with a sparse one, convert the latter
to moments and add. This gives one adaptive transfer, without choosing
a separate final actual degree.

**Storage theorem.** At every nonterminal stage, each adaptive key uses
exactly min(S,h) scalar weight coordinates, where h is the number of
old f-state cycle counts in that fibre. Its total number of scalar weight
coordinates is therefore no larger than the old f-state count.

**Proof.** Sparse storage has h coordinates and conversion is made only
when h>S. If a dense parent had h>S and takes a nonzero (k,e) transition,
its cycle counts shift injectively by m->m+k. None of their coefficients
vanishes. Contributions arriving at the same destination m have the same
sign by (M8), so merging cannot remove this shifted support. Hence the
destination has at least h>S>=S' cycle counts. It still needs the dense
form and no sparse reconstruction is required. Sparse-only mergers are
converted exactly at the threshold. At S=0 both algorithms keep the
same single summed terminal weight, discarding it if zero. QED.

The implementation retains even an all-zero **nonterminal moment vector**.
Its semantic cycle-count support can be nonempty, as (M7) demonstrates;
retaining it makes the support proof above literal. Discarding such a
vector would also be safe for its permitted continuation horizon, but
that is a separate exact-kernel pruning rule, not a claim that its
underlying prefix terms have disappeared.

This bound counts independent scalar weight coordinates. Key indices,
conversion work, rational bit lengths, and the reference implementation's
unused index-zero padding are outside the count. The theorem establishes
real arithmetic-state storage compression, not a wall-clock speedup or
an asymptotic polynomial bound.

## 7. Controls, costs and the remaining target

Run `python check.py` and `python -O check.py`. Both pass with explicit
exception gates. Every coefficient0..18 agrees by four paths: the dense
moment transfer, the adaptive transfer, the old scalar f-state transfer,
and direct unmerged
odd partitions with independent truncated substitution products. The
inherited pair-free prefix through10 is also checked. Exact rational row
reduction verifies all63 Vandermonde cases S=1..7,h=1..9; the realizable
same-defect families, horizon hostile and invalid inputs are retained.

The implementations give these summed stage sizes:

| N | Old scalar keys | New vector keys | Dense scalar cells | Adaptive scalar cells |
|---:|---:|---:|---:|---:|
| 6 | 7 | 7 | 24 | 7 |
| 10 | 22 | 19 | 88 | 20 |
| 14 | 49 | 38 | 222 | 45 |
| 18 | 94 | 67 | 463 | 85 |

These are exact structural work counts, not timings. The new key graph is
smaller from10 in this table, but the dense vectors cost more scalar cells.
Adaptive storage removes that regression and is checked against the old
scalar count separately at every stage, not just in the summed table.
The proved storage saving must not be advertised as a runtime gain.

The next implementation target is combining adaptive weights with the
inherited odd-factorial normalization, followed by fair conversion-cost,
bit-length and runtime measurements. The present result establishes the
exact linear quotient, a working adaptive transfer, and its sharp
information and scalar-storage boundaries.

The moment key forgets actual tournament arcs, module owners, and ports;
it preserves their exact signed counting total through the inherited
cycle-index identity. Literal clone-merger transport continues to use
the separate canonical-module cut-set carrier. Neither representation is
silently substituted for the other in a connectivity argument.
