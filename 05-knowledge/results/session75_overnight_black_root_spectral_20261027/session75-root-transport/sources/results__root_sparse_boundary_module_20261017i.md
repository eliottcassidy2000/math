# ROOT matching is a sparse boundary-forced transport

2026-10-17, clock_duality_oct15. **PROVED, scoped and self-audited:** the
all-source localization theorem below, its dyadic specialization, the
eventual clock-period theorem R3, the growing-depth matching window R4,
and the source-preserving frontier renewal laws R5--R6.
**FINITE-EXACT:** whole-support comparisons, including their omitted-mass
correction, and complete phase banks giving all-source dyadic error bounds
through depth5. **OPEN:** uniform energy at arbitrary depth and universal Collatz.
This is a reduction of the matching error, not a proof of the energy hypothesis.

## Inheritance and live board

Anchor: uniform ROOT energy. Niche: coupled golden contact circuits
([companion](golden_contact_circuit_20261017i.md)). Wildcard: finite extraction,
elementary inverses, and golden norm filters
([paper reading](finite_extraction_and_inverse_clocks_20261017i.md)).
Board: **literal source / omitted support / boundary forcing / integer clock /
balanced flow / contact circuit**.

Closest proved mechanism: R8 of
[root_braid_logistic_20261017f.md](root_braid_logistic_20261017f.md), and P2 of
[root_positive_transport_20261017g.md](root_positive_transport_20261017g.md).
Canonical hostile: transport expands some balanced vectors (THM-4667,
[energy ledger](../../01-canon/theorems/THM-4667-energy-ledger-of-the-signed-root-transport-partial-currents-merge-cross-terms-and-the-octave-trace.md),
statement7); sparse support alone permits an arbitrarily large ROOT atom.
Corrected near miss: the incoming matching census omitted ideal-only support.
Least-used sidecar: the exact lattice point next to a homogeneous band boundary.

Incoming work through6f078dca77 was read: ideal energy through depth12,
whole-source traces through2^24, the two matching scripts, THM-4668's
whole-linear-space Cartesian exclusion, and HYP-9289. No stronger all-depth
matching theorem was found in the targeted inheritance search. The finite
extraction manuscript suggests separating exact balance from approximation
cost; here that means retaining the entire error vector, including its zeros
on one side of the comparison.

## 1. Definitions on the unchanged source space

Let S be the positive odd integers other than5 mod8. Let B(s)=4s+1 and
P(s)=seed((3s+1)/2^a), where a=v2(3s+1) is1 or2 on S; seed removes B powers.
The signed pushforward is L delta_s=(-1)^a delta_P(s). In particular L fixes
delta_1. Use the original H_N and ideal band table K_d from the cited R8.

Define, on every t in S, not only on the actual support,

    c_(N,d)=L^d H_N,
    k_(N,d)(t)=K_d(t mod3^(d+1), b),
    b=floor(log2(M_d N/t)),  M_d=4(3/2)^d,
    k_(N,d)(t)=0 when b<0,
    e_(N,d)=c_(N,d)-k_(N,d).

All these vectors have finite support at finite N,d. The band inequalities
are closed on the right, as in R8. The ideal vector k is a lattice sampling;
its squared norm is not automatically the continuum ideal energy times N.

Write m=oddpart(N). For j>=0 define the explicit finite set

    J_j(m)={t in S: for some ell>=0,
              -1/3 < t - 3^j m/2^ell <= 1/2}.

Only 2^ell<=2*3^j*m can contribute. Thus

    |J_j(m)| <= floor(log2(2*3^j*m))+1.

Set D_0(m)=J_0(m) and D_(j+1)(m)=P(D_j(m)) union J_(j+1)(m).
These are literal integer sets with their genuine P edges.

## 2. Uniform support theorem

**R1 (PROVED).** For every positive integer N and every d>=0,

    supp(e_(N,d)) subset D_d(oddpart(N)).                         (1)

In particular,

    |D_d(m)| <= (d+1)(log2(m)+2)
                   + (log2(3)/2)d(d+1).                        (2)

For every dyadic source N=2^k, the same set D_d(1) works for every k>=0.
It has O(d^2) points and lies in1..3^d. This is uniform in both k and d.
It does not bound the coefficients at those points.

**Proof, initial boundary.** Replacing B^j(s)=4^j s+(4^j-1)/3 by4^j s
in H_N gives exactly the lattice vector k_(N,0); this is P2's homogeneous
source row. Each cutoff moves from theta to theta-(1-4^(-j))/3. A changed
primitive integer s must lie in(theta-1/3,theta], with theta an integer
power of2 times m. A positive exponent makes theta even, and no odd s can
be within this interval. Exponent0 or a negative exponent is included in
J_0(m). For m=1 the only possibility is s=1, agreeing with P2's six-periodic
ROOT correction.

**Proof, one transport.** An inverse branch of a primitive terminal t has

    s=(2^A t-c_A)/3,
    c_A=(3+2^a-2^A)/3,  A=a+2j, a in{1,2}, j>=0.              (3)

It is included precisely when s is integral; then s is positive, odd and
primitive, and its exact sign is(-1)^a. The residue recurrence for K uses
this same congruence and sign. The only discrepancy is its homogeneous
size test, which replaces s by2^A t/3.

For an input cutoff s<=M_d N/2^b, the corresponding terminal cutoffs are

    t<=theta+delta_A versus t<=theta,
    theta=M_(d+1)N/2^(b+A-1),
    delta_A=c_A/2^A=(3+2^a)/(3*2^A)-1/3,
    -1/3<delta_A<=1/2.                                     (4)

The difference of the two step functions is therefore supported in
(theta-1/3,theta+1/2]. Writing theta=3^(d+1)m times an integer power of2,
positive powers again give even integers and cannot meet this interval at
an odd integer. Every remaining point is in J_(d+1)(m). This argument
applies to each jump of the band function, including the outer cutoff.
For a fixed finite source envelope, only finitely many inverse branches
can contribute. Consequently

    b_(N,d+1):=L k_(N,d)-k_(N,d+1)
    has support contained in J_(d+1)(m).                     (5)

Now e_(N,d+1)=L e_(N,d)+b_(N,d+1), proving (1) by induction.
Counting at most one integer in each interval of width5/6 proves (2).
For m=1, each point of J_j(1) is at most3^j. Also
P(t)<=(3t+1)/2, which maps1..3^j into1..3^(j+1).
This proves the final assertion. No termination assumption was used.

The proof preserves the starting source, the terminal, every residue guard,
and the sign. The displayed d is primitive depth. If an implementation uses
canonical A letters, it must still retain the braid endpoint correction
physical clock=sum A_i+2j_0-2j_d. Localization does not change that clock.

## 3. What this buys, and what remains unpaid

The exact inhomogeneous recursion is

    e_(N,d)=L^d e_(N,0)+sum_(j=1)^d L^(d-j) b_(N,j).          (6)

Thus the matching problem is a sparse, boundary-forced system on moving
integer sets. It has not become an unmarked random walk. The candidate
states can be generated before computing any orbit amplitudes. The fixed
source N enters the forcing amplitudes and band phases, even when the set
D_d(1) is independent of the dyadic exponent k.

**R2 (PROVED, elementary).** Let I_(N,d)=sum_t k_(N,d)(t)^2 and
Q_(N,d)=sum_(t in D_d(m)) e_(N,d)(t)^2. Then

    |sqrt(E_d(N))-sqrt(I_(N,d))| <= sqrt(Q_(N,d)).             (7)

The same triangle inequality holds after summing over the original octave
of N. This separates three distinct bills: continuum ideal energy,
sampling that energy on the integer lattice, and the sparse forcing energy Q.
A single large coefficient can make Q large, so (1) alone pays none of them.

There is an automatic but weak stable norm. Put
W_d(e)=sum_t |e(t)|(t+1)/3^d. Since P(t)+1<=3(t+1)/2,

    W_(d+1)(Le)<=W_d(e)/2,
    W_d(e_d)<=2^(-d)W_0(e_0)+sum_(j=1)^d2^(j-d)W_j(b_j).     (8)

This is a true all-depth contraction preserving the sources. Its weakness
is explicit: a ROOT atom has weight2/3^d. Turning (8) into unweighted ROOT
energy costs the missing inner-boundary amplitude estimate. It is useful
for separating recent outer forcing from old inner accumulation, not for
claiming convergence.

Next precise targets: bound the octave Gram matrix of the vectors in (6);
look for cancellation among the injections b_j before taking absolute
values; and retain the six-periodic initial ROOT correction when testing
golden or ternary clock decompositions. An occupation certificate must
include these forcing rows, rather than only the homogeneous L edges.

## 4. Whole-support audit and correction

`root_full_support_20261017i.py` uses integer K recurrences and exact integer
band comparisons. It enumerates every primitive seed through the maximum
of the actual and homogeneous envelopes. Squared sums are checked to fit
int64. Independent scalar replay from each original odd source checks768
source/depth pairs, N=1..128 and d=0..5.
The same768 cases check the complete error support, and640 comparisons
check the one-step forcing support(5) independently.

The declared matching universe is90 cases: d=0..8 and
N in{4095,4096,4097,7910,16383,16384,16385,65535,65536,65537}.
Every error support satisfies (1). Examples of the omitted-mass correction:

| N | d | old support-only squared error | full squared error | omitted energy |
|---:|---:|---:|---:|---:|
|4096|5|18|34|16|
|4096|6|9|38|29|
|16384|7|180|184|4|
|65536|8|593|597|4|

At N=4096,d=6 the full error seeds are1,91,729. For N=7910,d=6 they are
1,11,89,353,5633: a non-dyadic source shows why the qualifier m=1 matters.
At the smallest depth, a minimal source-size witness is N=3,d=0: the actual
coefficient at seed3 is0, but K_0(0,2)=-1. N=1,2 have no omitted ideal-only
seed at depth0. The smallest dyadic depth0 witness is N=4, with actual
ROOT coefficient0 and ideal coefficient-1. These cases are included in the
independent small-source controls.
The old scripts' reported quantities remain reproducible on the restricted
support. They were not full l2 norms. The actual energy census and the
ideal energy ladder are unaffected.

Reproduce from the repo root:

```bash
python3 04-computation/experiments/root_full_support_20261017i.py --out /tmp/root-full.json
```

Receipt: `04-computation/experiments/energy_extraction_20261017i/root_full_support.json`.
No claim of an independent referee, universal energy control, or universal entry.

## 5. The error itself has an exact finite clock

The support theorem permits evaluating the forcing without enumerating the
whole source interval. A subsequent experiment suggested a periodic law;
the following proof establishes it independently of that observation.
Write e_(k,d)=e_(2^k,d), and p_d=2*3^(d+1).

**R3 (PROVED).** For every d>=0 and k>=d,

    e_(k+p_d,d)=e_(k,d) as vectors on the literal seeds.       (9)

Consequently the exact finite bank0<=k<d+p_d contains every possible
dyadic error vector at depth d. The maximum squared error over all k>=0
is a finite, effectively computable integer C_d. The effective construction
does not by itself bound C_d suitably as d grows.

**Proof.** The base case is P2's exact six-periodic ROOT atom. Consider the
forcing b_(2^k,d) for d>=1 at one of its possible terminals t in J_d(1).
Put q=3^d, so the preceding K table has residue modulus q and band period2q.
By (3), its inverse source at a legal exponent A is

    s_A=[2^A(3t+1)-(3+2^a)]/9.                              (10)

Legality is integrality; its sign is(-1)^a. Define

    C=2^(3-d)3^(d+1), D=3t+1,
    kappa=floor(log2(C/D)), g=D-C/2^(kappa+1)>0.

The gap is strictly positive because C/D cannot be a power of2: D is
prime to3 whereas C contains a positive power of3. Choose A0>=1 so that
2^A0*g>7. For A>=A0, the exact input band is

    floor(log2(M_(d-1)*2^k/s_A))=k-A+kappa.                 (11)

Indeed its nonintegral logarithm is k-A plus the logarithm of
C/(D-(3+2^a)/2^A); the chosen gap keeps this ratio in[2^kappa,2^(kappa+1)).

We may choose A0<=kappa+d+1. To see this, t<=3^d implies
kappa>=2-d. The positive number

    2^(kappa+d+1)*g
       =2^(kappa+d+1)D-8*3^(d+1)

is therefore a positive multiple of8, hence greater than7. All legal head
sources A<A0 satisfy

    s_A<2^(kappa+d)D/9<=8*3^(d-1)=M_(d-1)*2^d.

Thus, for k>=d, every head band is nonnegative and shifts by p_d when k
does. The target band is also nonnegative. R8's band periodicity leaves
these head and target terms unchanged.

For the tail A>=A0, (11) restricts it to A<=k+kappa. Changing k to k+p_d
keeps the common terms unchanged and adds exactly p_d=6q consecutive A
terms. Group them in triples A,A+2q,A+4q. Their signs, integer guards and
bands modulo2q agree. By the elementary3-adic valuation identity
v3(2^(2q)-1)=d+1, (10) shows that their source residues run through all
three lifts of one residue modulo q/3: consecutive differences have exact
valuation d-1 and agree modulo q. R8's last-digit balance annihilates the
sum of the three coefficients. Hence the added block vanishes.

This proves b_(2^(k+p_d),d)=b_(2^k,d) for k>=d at every possible forcing
terminal. Elsewhere both vectors vanish by R1. The error recurrence(6),
the induction hypothesis, and divisibility p_(d-1)|p_d finish (9).

The proof retains actual integer sources and the full native branch guard.
It uses an eventual period, not a numerical fit or a density statement.
Earlier phases k<d are kept explicitly in the bank.

## 6. Complete finite phase banks and their cost

`root_boundary_clock_20261017i.py` evaluates only the sparse forcing sites,
using exact inverse sources(10), integer band tests and R8's finite K tables.
It checks every k=0..1600 at depths0..5. This includes the complete banks
required by R3 and repeated phases beyond each bank. Eighteen cases agree
with the independent whole-support computation above, including every signed
error coefficient, and every tested period comparison succeeds.

Combining those exact finite banks with R3 gives the following **global
dyadic-source conclusions at the stated fixed depths**, not just samples:

| d | p_d | exact max_(k>=0) sum_t e_(k,d)(t)^2 | a maximizing k | error vector there |
|---:|---:|---:|---:|---|
|0|6|4|0|-2 delta_1|
|1|18|9|8|-3 delta_3|
|2|54|52|37|-4 delta_1-6 delta_9|
|3|162|484|123|22 delta_27|
|4|486|39593|212|-187 delta_1-68 delta_81|
|5|1458|4503509|756|delta_1+878 delta_15-1932 delta_243|

In particular, for every k>=0 and d=0..5,

    |sqrt(E_d(2^k))-sqrt(I_(2^k,d))|<=sqrt(C_d).              (12)

This strengthens a fixed-depth error with an unspecified log-source factor
to an exact source-uniform constant on the dyadic slice. The large phase
k=212 at depth4 also warns against extrapolating small-source amplitudes.
Its large absolute error is still tiny relative to N=2^212.

The clock p_d is exactly the order of2 modulo3^(d+2), the extra two ternary
digits arising from division by9 in(10). The earlier golden Q clock has
order4p_d at that same modulus. These are two observations of the same
integer exponent, not identical group elements; keep the exponent lift.
This is a concrete connection to golden clocks and the18-band first-depth
period, rather than a resemblance between unlabelled numbers.

The next sharper target is the **coupled phase/depth region** needed by the
energy argument. Bounding every one of the exponentially many clock phases
at each depth may be wasteful when the useful depth is tied to log2 N.
Neither the finite maxima nor their growth establish that coupled bound.

```bash
python3 04-computation/experiments/root_boundary_clock_20261017i.py --depth 5 --max-k 1600 --out /tmp/root-clock.json
```

Receipt: `04-computation/experiments/energy_extraction_20261017i/root_clock.json`.

## 7. A proved growing-depth window, with an explicit unpaid outer region

The finite-clock argument also supplies coarse bounds without constructing
its exponentially large phase bank. Define

    B_d=2^d * 3^(1+d(d+3)/2),
    R_d=2+6d(d+1)B_d,
    q_d=3^(d+1), L=floor(log2(M_d*2^k)).

Let G_d be the continuum ideal **squared** energy:

    G_d=sum_(b>=0) [M_d/2^(b+1)] [3/(8q_d)]
                          sum_(r modq_d) K_d(r,b)^2.

This is distinct from the l1 tariff gamma_d of R8.

**R4 (PROVED).** For all d,k>=0,

    max_(r,b>=0)|K_d(r,b)|<=B_d,
    ||e_(2^k,d)||_1<=R_d,
    |I_(2^k,d)-2^k G_d|<=B_d^2 [3q_d(L+1)+1].              (13)

Consequently, with H_(k,d)=R_d+B_d sqrt(3q_d(L+1)+1),

    |sqrt(E_d(2^k)/2^k)-sqrt(G_d)|<=2^(-k/2) H_(k,d).      (14)

This holds at every depth, but the bound becomes useful only in a specific
region. For every fixed c<1/sqrt(log2(3)), the right side tends uniformly
to zero for0<=d<=c sqrt(k) as k tends to infinity. Moreover, without assuming
ideal dissipation,

    sup_(0<=d<=c sqrt(k)) |E_d(2^k)/2^k-G_d| -> 0
             whenever c<1/sqrt(2log2(3)).                  (15)

These constants are about0.794 and0.562, respectively. They are sufficient
constants from a coarse bound, not claimed optimal thresholds. This proves
a growing-depth matching region; it does not bound G_d uniformly or the
arbitrarily deep ROOT coefficient.

**Proof.** The finite recursion for one K entry has at most b+1 terms.
Reduce b modulo p_d by R8. Thus its maximum is at most p_d times the
previous maximum. Starting from3 yields B_d exactly.

For forcing at depth j>=1, R3's proof reduces k to a representative
0<=k<j+p_j. Its input envelope U=M_(j-1)2^k is at most
4*3^(j-1)*2^p_j. By(10), any legal A with s_A<=U satisfies

    2^A <= (9U+7)/4 <= 3^(j+1)*2^p_j+7/4
                         < 2^(p_j+2j+3).

Hence A<=p_j+2j+2<=2p_j. There are at most2p_j summands,
so each forcing coefficient has absolute value at most
2p_j B_(j-1)+B_j=3B_j. The number of forcing sites is at most2j+2.
Signed pushforward is an l1 contraction, regardless of termination or
collisions. Since ||e_(2^k,0)||_1<=2, the recursion gives

    ||e_(2^k,d)||_1 <= 2+6 sum_(j=1)^d (j+1)B_j <= R_d.

For sampling, each primitive residue r modq_d consists of three progressions
modulo8q_d. Its count in a half-open band differs from3/(8q_d) times the
band length by at most3. Sum over q_d residues and the L+1 bands that
can contain a positive integer. The continuum tail occupies a total
interval shorter than1; its contribution is at most(3/8)B_d^2. This proves
the last part of(13), including the harmless larger constant1.
Triangle inequalities and |sqrt(x)-sqrt(y)|<=sqrt(|x-y|) give(14).

Finally, log2(B_d)=log2(3)d^2/2+O(d), and L=k+O(d+1).
The right side of(14) is at most
2^(-k/2+log2(3)d^2/2+O(d+log(k+2))), uniformly in that region.
For(15), use G_d<=(3/8)M_d B_d^2 and
|x^2-y^2|<=|x-y|(2y+|x-y|). The larger cross term is bounded by
2^(-k/2+log2(3)d^2+O(d+log(k+2))). Both claimed regions follow.

The unmatched region is now explicit: depths beyond a square-root window
in the binary source clock. Kathara's composition growth pays for a crude
finite recursion; three-lift balance removes the unbounded clock tail;
finite extraction then applies to the remaining phase bank. What is still
needed is cancellation that improves this quadratic-in-depth exponent,
plus control of the ideal energy and the late, inner-boundary ROOT mass.
The proof does not erase those obligations by choosing a larger source.


### Comparison with prior depth windows

THM-4493 ([gcd-height rank certificate](../../01-canon/theorems/THM-4493-collatz-gcd-height-rank-certificate.md))
already gives a square-root-logarithmic window with constant0.9600473 for
multiplicative independence along no-dip individual odd trajectories. That
larger constant does not dominate(15): the hypothesis, source and predicate
differ. It bounds pairwise gcd overlap in one orbit, whereas R4 bounds the
full signed coefficient vector from an interval of sources. Both quadratic
exponents come from multiplying depth-dependent arithmetic costs; this
identifies a common obstruction to improve, not a transfer theorem.

Likewise the d<=m/10 carry-cost window in
[outer packets, O3](outer_packet_golden_clocks_20261017h.md) is longer but
controls a different remainder. Its separate uniform o(N) result concerns
one signed outer packet. It does not estimate the full coefficient energy
or the inner forcing in(6). R4 therefore supplements those results; it does
not supersede them or reprove their stronger statements for their observables.


## 8. Incoming first-hit frontiers: an exact renewal interface

Incoming commits51674124d9,1dafe2fa3a anddd9b0d1f0f introduced
[root_frontier_penrose_20261017.md](root_frontier_penrose_20261017.md) and
[HYP-9290](../hypotheses/HYP-9290-frontier-energy-cone-of-the-stopped-root-transport.md).
Their first-hit frontier avoids transporting small coefficients onward to
ROOT. The enlarged census through N<2^28 is incoming evidence, not a new
reproduction here. The old C=2.4 suggestion at cutoff100 is superseded by
its E_100/N>2.62 record. At N=1, E_Y=10 for every Y>=3, so all-N constants
also need a finite-prefix payment. The conjecture with an unspecified C
remains open.

Let pi select the finite primitive frontier U={s<=Y}, and split the native
signed transport on low/high coordinates:

    L = [ A  B ],  A=pi L pi, B=pi L(1-pi),
        [ C  D ]   C=(1-pi)L pi, D=(1-pi)L(1-pi).

Here A,B,C,D name operators, not the previous exponent or coefficient bound.
The genuinely stopped operator is [I B;0 D]. For initial vector(a0,h0),
put g0=a0 and g_j=B D^(j-1)h0 for j>=1. These are signed **first-arrival
increments**; the stopped frontier vector at depth d is sum_(j<=d) g_j.
The moving remainder is D^d h0 and is not discarded.

**R5 (PROVED, formal renewal).** If a_d is the low projection of the
unstopped vector L^d(a0,h0), then

    a_d = g_d + sum_(ell=1)^d K_ell a_(d-ell),
    K_1=A, K_ell=B D^(ell-2) C for ell>=2.                  (16)

Equivalently, as formal power series only,

    a(z)=[I-zA-z^2 B(I-zD)^(-1)C]^(-1)
                        [a0+zB(I-zD)^(-1)h0].             (17)

**Proof.** Partition each native path ending at a frontier state according
to whether it has a previous frontier visit. No previous visit gives g_d.
Otherwise its last frontier visit is followed either by one low-to-low
edge or by an excursion through high states, giving exactly K_ell. Every
word is counted once with its original source weight and sign. The formal
geometric series has finite meaning at each coefficient; no analytic
convergence or eventual hitting assertion is required.

The same identity holds with full native edge words as labels. Retain each
edge's a and output braid j, hence canonical A=a+2j, and recover the literal
shortcut clock as sum A+2j_initial-2j_terminal. The variable z above counts
primitive stages, not physical elapsed time. Frozen waiting steps do not
advance a source's actual clock. Evaluating these labels prematurely would
lose precisely the marked-clock Gram information still needed for energy.

The frontier lemma survives without universal entry. Define f=0 off the
ROOT basin. Every initially ROOT-supported source eventually hits U, and
nonhitting sources have f=0. Since initial support is finite, the frozen
vector stabilizes after all hitting arrivals even if a moving complement
persists. Thus R(N)=<f,c_U> and R(N)^2<=|U| ||c_U||_2^2. At finite depth,
however, the energy identity is

    ||frozen_d||^2+||moving_d||^2
       =E_0+2(MM_hi,d+MM_lo,d+MF_d).                       (18)

The incoming statement that only above-frontier merges contribute was
incorrect: frozen arrivals MF remain, as its own decomposition measures.
The endpoint identity without moving energy is valid only when that vector
vanishes. These scope repairs are applied to the incoming note and hypothesis.

For a height frontier, C is usually nonzero. Thus neither the complete K
clock period nor R4's matching theorem automatically passes through the
first-hit restriction. The finite matrix in(17) contains infinite-memory
series; finite frontier size does not make those series finite polynomials.
This is the finite-extraction manuscript's degree obligation in a concrete
new form. The next sufficient bound is the source-wise two-time Gram sum
sum_(i,j)<g_i,g_j>, with the original integer clocks attached.

## 9. A forward-closed frontier removes the unknown return kernel

There is a more direct reroute than treating C as an uncontrolled error.
Enlarge the frontier until it is forward closed under the actual P map,
provided that enlargement has a finite certificate. This is a change of
stopping set, not of the original source or dynamics.

**R6 (PROVED, finite-set criterion).** Let V be a finite set containing1,
with P(V) subset V, and suppose a rank r:V->nonnegative integers has
r(1)=0 and r(P(s))<r(s) for s!=1. Stop on first entry into V. Then:

1. Every point of V is certified to reach ROOT, with its exact native word.
2. The outside vector of stopped and unstopped transport is identical at
   every finite depth: C=0 in the block matrix, so both are D^d h0.
3. If R=max r, the internal signed matrix satisfies A^R=Pi, where
   Pi=delta_1 f|_V^T, and A^j=Pi for j>=R. Therefore

       (I-zA)^(-1)=sum_(j=0)^(R-1) z^j A^j
                                      + z^R Pi/(1-z).      (19)

4. The frontier current inequality remains R(N)^2<=|V| E_V(N).

**Proof.** Strict rank descent proves1 and the stated power of A; the ROOT
self-edge has sign+1. Forward closure gives C=0, establishing2. Splitting
the formal geometric series at R gives3. The same finite-source argument
used after(17) gives4; universal entry into V is not assumed.

Exact closure witnesses, including every vertex and its descending rank:

| original height Y | closed frontier size | largest seed | maximum native rank |
|---:|---:|---:|---:|
|3|2|3|1|
|10|6|17|5|
|100|74|2051|42|
|1000|627|63665|64|

For example, replacing the375-point height1000 frontier by this627-point
certified set completely removes excursions from the frontier back into
its complement. This costs a known larger Cauchy-Schwarz factor; it does
not prove an energy comparison between the two stopping sets. The full
clock-labelled internal paths remain in the certificate. Formula(19) is
for primitive-depth bookkeeping and does not identify primitive and literal
shortcut times.

The external arrival series g(z) is still the open part. A finite controller
with a known transient polynomial and ROOT pole is now available; the
remaining task is to bound its source-specific influx rather than assume
that every source enters it. This is a concrete finite-extraction/renewal
bridge, and an example of changing the stopping set to eliminate a return
term instead of estimating it indefinitely.

`root_frontier_renewal_20261017i.py` checks(16),(18), and direct original-source
aggregation in3328 cases each (N=1..64, four frontiers, depths0..12), plus
87709 literal shortcut-clock and endpoint-braid comparisons. No test waits
for universal hitting. The four finite closure constructions have a10000-step
cap per new path, a cycle check, and independently checked descent ranks.

```bash
python3 04-computation/experiments/root_frontier_renewal_20261017i.py --out /tmp/frontier-renewal.json
```

Receipt: `04-computation/experiments/energy_extraction_20261017i/frontier_renewal.json`.
