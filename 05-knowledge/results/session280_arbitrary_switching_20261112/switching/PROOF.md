# Arbitrary switching: certified56- and151-letter languages and a three-block obstruction

**Session280. PROVED by the rational certificate below, independently reviewed.**
The theorem concerns observation of a physical rooted interface after a **known**
switching word. It does not assert that the language covers all allowed odd
Hamiltonian-path counts. Chan–Pak remains OPEN.

## 1. Inheritance and the missing coordinate

Use the inherited moving-root block construction and cone from
[Session279, physical Sections1 and8](../../session279_interval_tightening_20261112/physical/PROOF.md):

    C = {(A,W,d): 0<=d<=A<=W},
    r0=(1,1,0), r1=(1,1,1), r2=(0,1,0).

An active binary word w starts1 and ends0. With P its binary value,
Q=sum_j (1−w_j)2^j, t=2^|w|, h=g(P), its physical block matrix is

    M_w = [[h,Q+1,-1], [P+1,t,0], [1,0,0]],
    g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod2) (v>=1).

The physical operation appends the one-apex tournament and reverses the arc
between successive roots. These matrices and their physical realization are
inherited, not claimed new. They preserve C and are strictly positive in
coordinates (A,o,e)=(A,(W+d)/2,(W−d)/2).

**Closest proved mechanism:** one fixed Pisot block eventually reveals the
integer interface from its scalar output (Session279). **Hostile:** the
individually non-Pisot word10101010, and the new switched triple below.
**Corrected near miss:** common ratio convergence is not absolute uncertainty
contraction. **Useful sidecar:** the actual preceding letter, rather than just
the current letter's eigenvalues.

An unrelated source was [Session172, Collatz contextual memory Sections2–4](../../session172_memory_dynamics_20261103/collatz/NOTE.md).
Its preserved finite observations do not supply an autonomous update after a
missing next bit. The concrete transfer here is to carry the immediately
preceding block as an observer-state label. It chooses a quadratic error
metric; discarding that label makes the contraction proof fail. This is a
connection of proof mechanisms, not an equivalence with Collatz convergence.
The nonautonomous square-sum transfer in
[Session162](../../session162_mrb_square_20261103/square_transfer/NOTE.md)
also emphasizes retaining an absolute source together with its changing
response. Here scalar size is explicitly retained throughout normalization.

## 2. A transverse cocycle, not a product of fixed Perron eigenvectors

Let x be a reference ray and y another ray with positive A coordinates. Define
its *scalar-conditioned residual*

    z=(y_W−(y_A/x_A)x_W, y_d−(y_A/x_A)x_d).

Thus y=(y_A/x_A)x+(0,z_W,z_d). Write q=Q+1 and p=P+1.
After applying M=M_w to both, put r'=x'_W/x'_A and s'=x'_d/x'_A.
Direct subtraction at the NEW scalar coordinate gives

    z'=K_w(r',s')z,
    K_w(r',s') = [[t−q r', r'], [−q s', s']].             (1)

This identity is exact. No approximation to a Perron ray is made. K depends
on the current reference ray, so replacing it by the fixed-block limit would
lose the switching history.

Let u be the preceding block and v the next block. Since x is in M_u C,
its normalized next image is in the triangle with vertices

    normalize(M_v M_u r_i), i=0,1,2.                     (2)

The coefficients are positive scalar-coordinate barycentric coefficients,
not necessarily the original cone coefficients. Formula(1) is affine in the
two output ratios. If a symmetric positive definite Q_v satisfies

    K_i^T Q_v K_i <= tau Q_u, i=0,1,2,                 (3)

then the same inequality holds throughout the triangle. Indeed for
K=sum_i alpha_i K_i, alpha_i>=0, sum alpha_i=1,

    sum_i alpha_i K_i^T Q_v K_i − K^T Q_v K
      = sum_i alpha_i (K_i−K)^T Q_v (K_i−K) >=0.

This is the exact convexity step that turns finitely many rational tests into
an all-input statement. The familiar quadratic-norm method is not new; the
explicit physical language, its finite certificate, and its scalar-interface
consequence are the results here.

## 3. A 56-letter arbitrary-switching theorem

Let L56 consist of every active word of length2 through7 having at most three
bit transitions. Such words have two or four constant runs. Its size is

    sum_(n=2)^7 [binom(n−1,1)+binom(n−1,3)] = 56.

`metric56.json` supplies a rational symmetric2x2 matrix Q_u for each letter,
using the common denominator10^6. The standard-library verifier checks:

1. The list is exactly L56, with no inherited filter hidden.
2. Every Q_u−(1/2)I is positive definite.
3. For all56² ordered transitions u→v and all three rays in(2),
   tau Q_u−K_i^TQ_vK_i is positive definite, where tau=249/250.

There are9,408 rational matrix inequalities. A2x2 symmetric matrix is positive
definite exactly when its leading entry and determinant are positive (the
verifier additionally checks the second diagonal). Thus the certificate needs
no floating eigenvalues, semidefinite solver, or conjectured stability test.
The metrics were discovered numerically, then rounded and checked with a
strictly weaker rational contraction factor. Discovery does not enter the
validity proof.

**Theorem S56.** For every switching word w1,...,wk in L56, and each ordered
pair of the three initial cone rays, its residual after k blocks obeys

    z_k^T Q_wk z_k <= tau^(k−1) z_1^T Q_w1 z_1.          (4)

Proof: x_j is in M_wj C because the cone is invariant. Apply(3) at every
transition, then multiply the inequalities. There is no bounded-depth or
periodicity restriction. In particular, the actual choices can be arbitrary;
no assumption that products of individually Pisot matrices stay Pisot occurs.

## 4. Absolute target normalization and an effective integer decoder

To turn(4) into an absolute coordinate error, one more uniform factor is
needed. For every w in L56, the three images M_w r_i in positive coordinates
(A,o,e) satisfy pairwise coordinate ratios between1/2 and2. This is checked
exactly in the certificate. All later positive matrices preserve these
componentwise comparisons. Consequently, for every nonempty switching word,

    a_max/a_min <=2,                                    (5)

where a_i is the A coordinate of its ith image ray.

A cone source x has unique decomposition

    x=(A−d)r0+d r1+(W−A)r2,

whose nonnegative coefficients sum to W. If W<=B, its output scalar N obeys
N<=B a_max. For any two image rays x,y, their conditional output difference at
scalar N is (N/y_A)z. Therefore(5), Q_u>(1/2)I and(4) imply, for EACH coordinate,

    width² <= 8 B² E0 tau^(k−1),                        (6)
    E0 = max_(w,i!=j) z_1(i,j)^T Q_w z_1(i,j).

The verifier retains BOTH orders i,j because scalar-conditioned residuals
are asymmetric. It obtains

    E0 = 6112850971/845000000.

For B=96, exact integer/rational comparison in(6) gives width<1 for every
k>=3292. This threshold is conservative and not claimed optimal. It is a
uniform theorem over all56^k possible switching words, not a finite table of
successful instances. A tighter per-word ray interval can decode much earlier.
For any other fixed B, iterate the rational factor until(6)<1; termination
follows from tau<1.

For an exact target N and a supplied switching word, compute its three image
rays, then the two tight coordinate intervals

    N min_i(W_i/A_i) <= W <= N max_i(W_i/A_i),
    N min_i(d_i/A_i) <= d <= N max_i(d_i/A_i).

At width<1, a true output has a unique integer W and d. For arbitrary N,
reject empty intervals and retain ambiguity when not yet narrow enough. Then
invert each known letter in reverse chronological order:

    A_old=d,
    W_old=(W−(P+1)d)/t,
    d_old=h d+(Q+1)W_old−A.

Require exact divisibility by t, cone membership, and finally membership in a
supplied bank of actual physical seed packets. These last conditions are not
implied by real contraction. Conversely, a surviving actual seed followed by
the supplied legal block word constructs the promised tournament.

Thus known-word scalar observation is injective on bounded integer cone
packets at sufficiently great depth. It does NOT distinguish different
isomorphism classes with the same packet, recover an unknown switching word,
or prove that some word accepts every prescribed odd target.

## 5. Individual Pisot letters really can fail under switching

Take the chronological block word

    (10, 101010, 1000).

Each block individually passes the exact Session279 criterion f(1)<0 and
f(−1)<0, hence has a Pisot Perron root. Their product is

    P=M_1000 M_101010 M_10
     =[[6801,8157,−923],[7125,8545,−967],[221,265,−31]],
    f_P(z)=z³−15315z²−19568z−4096,
    f_P(−1)=156>0, f_P(−2)<0.

So P has a real conjugate in(−2,−1). Repeating this switching period expands
an absolute transverse direction. This is not a numerical near-unit verdict:
the endpoint signs certify it exactly. All1,023 active triples with total
word length<12 pass both signs; the displayed total12 is therefore the first
possible total length for a bad three-letter period. This minimality is scoped
to active triples, not all constructions or all period lengths.

Why this obstructs the blanket extension: the three initial cone rays span
R³. Their normalized image triangle under P^m has a nonzero component along
the expanding non-Perron eigendirection. Conditioning at a physical seed
scalar N_m asymptotic to a positive multiple of lambda^m cancels the Perron
normalization, leaving a width component of order |mu|^m. Thus the entire
source-cone slice does not shrink to subunit absolute width. This refutes
uniform cone-based absolute contraction for unrestricted switching of all
individually Pisot letters. It does **not** prove scalar collisions between
integer physical seeds or rule out arithmetic decoding by another method.

The bad letter101010 has five bit transitions and is outside L56. The safe
language excludes precisely such high-oscillation blocks at these lengths;
no claim that this restriction is necessary or globally optimal is made.

## 6. Connection to the coverage task and remaining experiments

The exact certificate supplies an infinite switched language with recoverable
root interfaces. It is not a coverage theorem. The companion coverage lane
separates branch growth from distinct scalar support: even many legal words
can coincide or miss an absolute target.

A natural larger language consists of the98 active words through length8
with at most three transitions. It has favorable branch-pressure evidence
(certified separately), but the current last-letter quadratic searches did
not certify its contraction, even after restricting to the six-ray invariant
projective hull of all first images. This is a failed certificate search,
not a counterexample. Random period tests are not promoted to stability.

The next missing sidecar may be a finer context label, a nonquadratic norm,
or a target-dependent switching policy. The switched triple shows why
individual spectra cannot silently replace that work.

## Reproduction

    python3 switching/verify.py
    python3 -O switching/verify.py

Both execute the same explicit checks, with no Python assertions removed by
optimization. `metric56.json` is the proof certificate; `results.json` records
its hash, test universe, all-depth constants and hostile control. The exact
finite inequalities plus the induction in Sections2–4 prove the infinite-word
claim. They do not enumerate or claim all odd tournament counts.

## 7. A compatible 151-letter extension, with supercritical presentation supply

The negative98-letter metric search was not the end of the construction.
Freeze a successful small-language metric, then add a proposed block only if
its incoming, outgoing, and self-transition inequalities admit a compatible
new metric. This is a *compatibility-first* enlargement, rather than adding
all individually attractive blocks and hoping their products behave.

The resulting explicit language L151 is every active word of length2 through9
with at most three bit transitions, EXCEPT these eleven words:

    10000010, 11100000, 11110000, 11111000, 11111110,
    100000010, 111000000, 111100000, 111110000,
    111111000, 111111110.

Thus it contains all56 prior letters,37 length8 letters, and58 length9 letters.
The exclusions are a discovery-dependent sufficient choice, not a theorem of
necessity. In particular, some excluded words are monotone and individually
excellent decoders. The issue is compatibility with the already chosen
metrics, not intrinsic badness.

`metric151.json` supplies all151 rational metrics with denominator10^8.
`verify_extension.py` checks exact equality with this declared language and
all68,403 transition/ray inequalities, now with

    tau=9991/10000,  Q_u>(1/4)I,
    E0=12225700863/1690000000.

The SAME proof in Sections2–4 applies. In particular all known switching words
over L151 have decaying absolute transverse uncertainty, and W<=96 source
packets are uniquely recovered at every length k>=15417 by the conservative
uniform bound. This weaker numerical threshold is the cost of certifying the
larger language with this metric, not a claim about optimal decoding depth.
The new threshold uses4B²E0 tau^(k−1)/(1/4)<1 exactly.

The companion [coverage certificate](../coverage/PROOF.md) establishes
supercritical **presentation** growth for this SAME list. This removes a real
compatibility mismatch: an ample library and a stable library need not be
different objects. It still does not make distinct words have distinct scalar
counts, fill all absolute-target gaps, or provide a target-dependent inverse
selector. The retained open obligation is precisely that support/ownership
step, not the existence of a switching-stable library with enough formal
branch supply.

Reproduce the extension independently of any optimizer:

    python3 switching/verify_extension.py
    python3 -O switching/verify_extension.py

The optional `discover.py` demonstrates the cutting-plane metric search for
the56-letter baseline. `greedy_discover.py` demonstrates compatibility-first
extension; both require NumPy/SciPy and produce only NUMERICAL CANDIDATES.
The two standard-library rational verifiers, not those discovery programs,
are the proof certificates.
