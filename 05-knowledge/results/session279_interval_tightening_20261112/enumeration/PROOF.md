# Adaptive symmetry backbones and a finite marked-moment compiler

**PROVED IN NOTE; finite implementation controls are separately stated.**
This counts ordinary unlabelled tournaments (A000568), not Hamiltonian-count
fibres. Chan–Pak remains open. No claim of a faster full exact counter.

## 1. Inheritance and the new obligation

Session278 `enumeration/PROOF.md` supplies the gcd-to-min upper envelope,
Bonferroni correction coefficients, two-sided signed remainders, and exact
resummation of every interaction involving a length-three cycle. Its final
algorithmic question is computing low correction coefficients without
unrestricted partition enumeration. Session277 `engine/PROOF.md` supplies
absolute-interval/residue reconstruction, which is unchanged here.

The new result is a closed finite **marked-moment hierarchy** computing every
fixed number of those coefficients in polynomially many arithmetic operations,
for a fixed finite backbone of exactly handled cycle lengths. Explicit first-
and second-order implementations accompany the general recurrence. More
backbone lengths also yield markedly tighter certified intervals.

Unexpected route: THM-685 `kronecker-transfer-and-relation-lattice-resummation`
warns that coherent relation families should be resummed, not discarded in a
low-order expansion. Here the coherent families are all interactions touching
specified cycle lengths. This is an analogy of operations, not a transfer of
an LRC theorem. THM-4845 `a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four`
already supplies cycle-source/Gaussian response. The finite moment closure below
uses a different, envelope-weighted correction hierarchy; it does not claim
that the full gcd response is finite-dimensional at every accuracy.

## 2. Any finite backbone, including non-prefix choices

Let B be a finite **set** of odd integers >=3. For odd lengths a,b define

    K_B(a,b) = gcd(a,b),  if a in B or b in B,
               min(a,b), otherwise.

For an odd cycle partition lambda, replace the cross-cycle gcds in the usual
Burnside exponent by K_B, leaving diagonal/single-cycle terms unchanged. Call
the resulting positive weight W_B(lambda). For each pair of distinct cycle
occurrences outside B put

    d = min(a,b)-gcd(a,b),     q=1-2^(-d).

Discard zero q. The exact Burnside weight is W_B product(1-q). Thus the inherited
signed inequalities apply verbatim to

    C_j(B,n)=sum_lambda W_B(lambda)e_j(q(lambda)).

Adding a length to B moves ALL its interactions into the baseline. It is not
keeping only small-support partitions: arbitrarily many cycles of those lengths
are treated exactly. Every B yields a valid interval. Intersecting intervals
from several B is safe; monotonicity of a particular signed truncation as B
changes is **not assumed**.

## 3. Backbone-conditioned size transfer

Condition on multiplicities t_a, a in B, with s=sum a t_a<=n and r=sum t_a.
Their exact internal weight is

    w_B = 2^[sum_a (t_a(a-1)/2+a choose(t_a,2))
                 +sum_(a<b) t_a t_b gcd(a,b)]
          / product_a a^t_a t_a!.

Use DP mass m only for cycles outside B, including singleton cycles. Initially

    V(f)=2^[choose(f,2)+fr]/f!,     0<=f<=n-s.

Process the remaining odd lengths b>=3 increasingly. Adding j occurrences has
mass increment bj and multiplier

    A_b(m,j)=2^[jm+j(b-1)/2+b choose(j,2)+jR_B(b)]/(b^j j!),
    R_B(b)=sum_(a in B) t_a gcd(a,b).                         (1)

All old off-backbone lengths are smaller than b, so their min-interaction with
each new occurrence is exactly old mass m. Interactions with the backbone use
R_B exactly, even if a backbone length exceeds b. This proves the size transfer
for C0. After the last length, multiply by w_B and sum all backbone choices.
The backbone cardinality, not the maximum length, controls its outer dimension.

## 4. First jet: C1 without unrestricted partitions

Alongside V retain C=sum W e1(q) and M_a=sum W m_a for every already processed
off-backbone non-singleton length a. Set

    Q=sum_a q_ab M_a,    q_ab=1-2^[gcd(a,b)-a].

For the transition (m,j) with multiplier A, add

    V' += A V,
    C' += A(C+jQ),
    M_a' += A M_a  (old a),      M_b' += A jV.               (2)

There are j m_a new distinct edges of weight q_ab, hence (2) is exact. Two old
partitions with the same mass need not have identical response; their weighted
first moments carry exactly the lost information needed at this order.
`backbone.py::envelope(marked=True)` implements (1)-(2).

For fixed |B|, there are O(n^|B|) backbone choices, O(n^3) size transitions per
choice (a deliberately loose bound), and O(n) moment work per transition.
This is O(n^(|B|+4)) arithmetic operations and O(n^2) coefficient storage per
backbone choice, rather than visiting every unrestricted odd partition.
The initial rational implementation is not faster at the tested sizes; the
integer compilation below changes that practical conclusion for this limited
C0/C1 workload, not for full enumeration.

## 4a. Compile denominators instead of paying for them at every transition

For a fixed target n let Q=2^n. Normalize each size-m state by m!; normalize
C1 by the additional factor Q. Because d<=n, every

    qhat_ab = Q q_ab = 2^n-2^(n-d_ab)

is an integer. The transition multiplier becomes

    Ahat_b(m,j)=2^[exponent in(1)] (m+bj)!/(m! b^j j!).

This is an integer: it counts choices and permutations of j disjoint b-cycles
on the new labels, times the displayed power of two. Apply(2) with Ahat and
qhat and the correspondingly scaled C1; all transition arithmetic is integral.
At the end use the integer factor

    n! 2^E_B/(z_B (n-s)!).

It is integral because z_B divides s! and n!/(s!(n-s)!) is binomial. Divide the
final sums by n! and Q n! respectively. There is no rounding or modular alias.

Factorial clearing is inherited from the older Burnside engines and Session277;
the new use is compiling the finite marked hierarchy and its dyadic correction
coordinate into integer transitions. For general ell, 2^(ell n)m! clears the
corresponding moment denominators, although our optimized code implements only
the first jet. `integer_jet.py` agrees with direct partition C0/C1 through80 and
with the rational jet through96. In a matched single-run Python comparison:

| n | integer first jet (s) | direct partition C0/C1 (s) |
|---:|---:|---:|
|32|0.00695|0.00609|
|48|0.0437|0.0684|
|64|0.1748|0.5343|
|80|0.5561|3.4949|

At96 the integer jet took1.767s versus16.470s for the rational jet. These are
measured local workloads, not an asymptotic benchmark or a comparison against
the repository's optimized full GMP counter. They do demonstrate a practical
improvement for the specified first-correction computation after the earlier
rational prototype was slower. The arbitrary fixed-order theorem is separate
from these implementation timings.

## 4b. Primitive binary words compress the correction response

Let P(d) be the number of binary words of length d with least period d. The
classical period decomposition gives

    2^g = sum_(d|g) P(d),     P(d)=sum_(e|d) mu(d/e)2^e.

This is the same Möbius/primitive-period mechanism appearing in the repository's
Witt and necklace work (for example THM-505
`the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect`), but
here it is applied to the **correction kernel**, not used to identify simple
cycles with primitive walks. For previously processed a<b,

    q_ab = 1-2^(-a) sum_(d|a,d|b) P(d).

Thus the first correction response can be retained as

    T=sum_a M_a,
    R_d=sum_(a:d|a) 2^(n-a) M_a,
    Qhat_b=2^n T-sum_(d|b) P(d)R_d.                         (6)

Only odd d<=n/3 are needed. For distinct positive odd a<b<=n, b/gcd(a,b)
is an odd integer >=3, so gcd(a,b)<=n/3. No larger common divisor can occur
in the future interaction being queried. This boundary uses distinct lengths;
within-equal-length interactions were already exact and have zero defect.

The integer transition is

    T' += Ahat(T+jV),
    R_d' += Ahat[R_d+1_(d|b) j 2^(n-b)V],
    C1' += Ahat(C1+jQhat_b).

This retains a divisor response, rather than each cycle length separately.
It is a legitimate quotient for the first jet because (6) factors every
future response. The backbone-dependent external field remains separately
retained. It is not a lossless quotient for every higher correction order
without corresponding mixed divisor moments.

`divisor_jet.py` agrees with the integer length-moment jet through128. It uses
floor((floor(n/3)+1)/2) divisor coordinates plus T instead of all previous odd
length coordinates. This is an exact structural compression. Local times were
not uniformly smaller: at64 about0.211s versus0.200s; at128 about7.68s versus8.77s.
No robust further speedup is claimed from those single runs. A useful outcome
is the new precise bridge: **primitive-period counts implement the response of
an enumeration-error certificate**. A numerical resemblance to necklace counts
has become an explicit map and sufficient-state identity.

## 5. All fixed orders close: the triangular precision budget

For processed off-backbone lengths and a nonnegative integer multi-index alpha,
retain

    J_(ell,alpha)(m) = sum W e_ell(q) product_a m_a^alpha_a,
                    ell+|alpha| <= K.                      (3)

K is any prescribed fixed correction order. Upon inserting j new b-cycles,
the new-edge correction polynomial is

    P_b(z;m,j)=product_(old a) (1+z q_ab)^(j m_a).

Write P_h=[z^h]P_b. For fixed j, P_h is a polynomial in old multiplicities of
total degree at most h: use binomial(j m_a,h_a), a degree-h_a polynomial,
and sum h_a=h. This is true even when j=0 and at empty states.

The exact update for a destination state with old multi-index alpha and new
b-exponent beta is

    A j^beta sum_(h=0)^ell sum_gamma [m^gamma]P_h
                              J_(ell-h,alpha+gamma).       (4)

Every source needed in (4) satisfies

    (ell-h)+|alpha+gamma| <= ell+|alpha| <= K-beta <= K.

Therefore the finite family (3) closes: no untracked higher moment is needed.
At the end, J_(j,0) gives C_j for every j<=K. This is a constructive exact
recurrence, not a moment-closure approximation. The convention j^0=1 includes
j=0; no negative exponents arise.

There are O(n^K) moment indices per mass and finitely many correction levels.
Naively expanding all polynomial coefficients in (4) gives the loose bound
O(n^(|B|+2K+3)) arithmetic operations for fixed B,K and polynomial storage.
Binomial-polynomial coefficients can be computed with signed Stirling numbers
or multiplication of the finitely many linear factors; neither requires
unrestricted partition enumeration. Rational numerators/denominators have
polynomial bit length for fixed K: weights have O(n^2) exponent bits, the
correction factors add at most O(Kn), and all factorial denominators divide
n!, up to the fixed-order coefficient arithmetic. Thus this is also a
polynomial-bit algorithm at each fixed correction order.

This does NOT imply polynomial-time exact A000568 enumeration: the K or B
needed for integer accuracy need not remain fixed as n grows (section8).
The general compiler is given by (3)-(4); production code currently implements
K=1 and K=2, not an optimized arbitrary-K backend.

## 6. Explicit second jet and why squaring the mean is wrong

For K=2 retain raw second moments S_ac=sum W m_a m_c and mixed moments
N_a=sum W e1(q)m_a. Define

    Q2=sum q_a^2 M_a,    QQ=sum_(a,c) q_a q_c S_ac,
    NQ=sum q_a N_a.

The second correction update is

    C2' += A [C2+jNQ+(j^2 QQ-jQ2)/2].                       (5)

The last expression is e2 of the jm_a new edge occurrences, not e2 of distinct
length types. In particular QQ is not Q^2/V. The remaining updates are

    N_a' += A[N_a+j sum_c q_c S_ac],   N_b' += A[jC1+j^2Q],
    S_ac' += A S_ac,                  S_ab' += A jM_a,
    S_bb' += A j^2V.

`second_jet.py` implements these formulas. This example makes the reusable
principle precise: spend one unit of retained multiplicity degree for each
unit of future signed-correction order.

## 7. Sharp finite support thresholds and actual interval gains

Take B={3,5,...,a-2}, where a>=5 is odd. A residual interaction uses lengths
at least a. Four residual edges require at least four cycle occurrences.
With four occurrences, fewer than two lengths exceeding a gives at most three
edges; hence support is at least 4a+4. Five or more occurrences cost at least
5a>=4a+4. The partition (a,a,a+2,a+2) attains four edges at support4a+4.
Consequently S3 is exact for n<4a+4 and is not exact at n=4a+4.

The **two-sided interval** using C4 is stronger than S3 alone. If M<=4
residual edges, its C4 lower and upper remainder bounds coincide (or vanish).
Five edges require support at least4a+6: for four occurrences, three distinct
lengths are necessary and cost at least a+a+(a+2)+(a+4); for five occurrences,
at least two exceed a and cost at least5a+4; six cost at least6a.
The partition (a,a,a+2,a+4) attains five nonzero edges at4a+6.
Thus the rational interval obtained from S0..S3 and their two-sided C_(k+1)
bounds is exact through n=4a+5 and first has positive rational width at4a+6.
At that sector every q is strictly between0 and1, so each relevant nontrivial
error bound is strict. This does not assert the first integer-rounded width
for arbitrary a without checking its magnitude.

In particular B={3} already gives an exact interval through25 although S3
ceases to be exact at24. The integer-rounded interval first has32 bits of
width at26 in the supplied experiment.

With the **same four signed partial sums**, sufficient low residue bits at64:

| Backbone B | Bits |
|---|---:|
|3|738|
|3,5|493|
|3,5,7|310|
|3,5,7,9|187|
|3,5,7,9,11|121|

These high-order table coefficients still use partition enumeration; the
polynomial jet is independently controlled at lower orders. Do not present
the121-bit interval as generated by the implemented second-order jet.
The exact reference integer is used only as a check, not to define bounds.

## 8. Fixed hierarchy depth still cannot settle every integer

Fix finite B,K. Choose incompatible odd a,b outside B and r=K+1 copies of a,
plus one b and f singleton cycles. There are r residual edges, all q>0.
Its signed S_K error is W_B q^r, which tends to infinity as f grows because
its numerator includes 2^[choose(f,2)+f(r+1)] and denominator f! times a
constant. All sector errors have the same sign. Hence fixed B,K cannot support
universal nearest-integer recovery from S_K alone. Higher two-sided bounds
can exactly repair this particular r-edge family, but choosing more residual
edges than the retained depth again leaves a positive, growing gap. Explicitly
use r=K+2; every bound through depth K has a strictly positive sector slack,
which scales by that same diverging singleton factor.

The useful algorithmic choices are adaptive B, adaptive K, exact residues, or
a compatible alternative certificate such as the sibling ground-state
resummation. Each preserves an absolute target; none turns ordinary class
mass into Chan–Pak existence. A physical tournament realization still needs
its full endpoint/interface data.
