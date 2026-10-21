# A conditional scale cone, and why the native lift does not supply it

2026-10-21. **PROVED:** six-axis compression of the route-free gate, a
conditional one-scale propagation cone, exact native source-lift density,
and the guard/clock translation cylinder. **FINITE-EXACT:** a partial-source
hostile refutes a tempting fixed-cut scale contraction. **OPEN:** uniform
153N coverage, actual satisfaction of the propagation cone, and Collatz
entry. This is separate from the frozen `routefree_suffix_gate_20261021`
packet; none of its receipts or claims are changed.

## 1. Exact finite-dimensional certificate state, not closed evolution

Assume all old sources have arrived at the selected cut h. Keep the base
packet's actual terminal vectors T_u=W_u+H_u, physical increment Delta0,
and unfinished new source charges. Opposite alphabet vectors satisfy
v_(i+6)=-v_i. If n_i is the pending multiplicity of angular type i, define
six antipodal-axis counts and supports

    c_i=n_i+n_(i+6),
    ell_i=max_(u in V)|beta(T_u,v_i)|,  0<=i<6.

The signed-sum zonotope depends exactly on c, since the summands for v and
-v are the same segment. Thus D(n)=D(c), with a representative multiplicity
vector (c_0,...,c_5,0,...,0). The linear support contribution is exactly
L=c·ell. Consequently the all-future sufficient inequality is

    Delta0+D(c)+c·ell <=153N.                       (1)

This compresses the **certificate value** to thirteen integers
(Delta0,c_0,...,c_5,ell_0,...,ell_5), plus the old-completion flag. It does not
compress the underlying dynamics faithfully. The literal source/guard/clock
ledger remains the authentication sidecar; the next section shows a lost
coordinate needed to update the state.

### Conditional one-scale propagation theorem

Let N'=64N, with selected cuts h,h' respectively, and suppose all old sources
have arrived at both cuts. Let the actual certificate states satisfy

    c'<=8c coordinatewise,
    ell'<=8ell coordinatewise,
    Delta0'<=64Delta0.                             (2)

Then

    Delta0'+D(c')+c'·ell'
       <=64[Delta0+D(c)+c·ell].                    (3)

Indeed Z(c') is contained in8Z(c), so Q-homogeneity gives D(c')<=64D(c).
Nonnegativity of counts/supports gives c'·ell'<=64c·ell. Adding the assumed
increment inequality proves(3). Thus a valid suffix certificate at N
propagates to N' under(2). If all earlier cuts through h' are also covered,
it gives all-horizon coverage at N'. This can be iterated wherever the
actual inequalities hold; it is not asserted that they hold natively.

The inequalities are sufficient, not necessary, and coordinatewise failures
do not refute coverage. They identify a precise missing scale obligation:
a mechanism controlling pending multiplicity and terminal coherence together.

### Exact nested-domain interface, and an entry-strength caveat

There is a genuine lossless part of the scale recursion which does not use
an affine source lift at all. At the same horizon h, the old source domain
at64N is exactly the union of the old and new source domains atN, since
4(64N)=256N. The relevant charges are identical: old scale64N uses A_(64N),
and A_(64N)=A_N on the previous old subdomain. Hence, terminal by terminal,

    W_(64N)(h)=W_N(h)+H_N(h),
    c_old,(64N)(h)=c_old,N(h)+c_new,N(h).             (4)

These are actual same-source identities, retaining terminal, depth, full
clock and both guards unchanged. For h'>=h, monotonic first arrival gives
c_old,(64N)(h')<=c_old,N(h)+c_new,N(h) coordinatewise. The outer new source
population and its arrivals/support remain missing; equation(4) alone does
not close the recursion. Squared zonotope energies do not add when the
pending multisets are joined, so replacing D(c_old+c_new) by D_old+D_new is
not justified.

There is also an important logical cost to the simpler old-complete cone.
For any base N0>=1, the statement

    for every k>=0 there is a finite cut h_k at which all old
    sources of scale64^k N0 have arrived in V

is **equivalent to universal native primitive entry into V**. In one
direction every primitive source lies in some old domain along the scale
tower. In the other, each old domain is finite, so universal entry supplies
a finite maximum first-arrival depth for that domain. This is not a weak
technical precondition for uniform iteration: it already contains the
unproved entry obligation. For example, unfinished new sources at one scale
become unfinished old sources at the next scale at the same horizon.
The base packet's general pending-old gate, or its mixed refinement, avoids
assuming this old-completion prerequisite. The conditional13-integer cone
is valid and useful where completion is certified, but cannot silently be
used as a convergence-free uniform induction.

## 2. A fixed-cut native scale contraction actually fails

The frozen base receipt has N=4096,h=105, pending nonzero new count386,
and D=91876. It has axis totals

    c=(88,38,95,32,108,25).

At N'=262144=64N,h'=111=h+6, the new checker exhibits4096 distinct unfinished
sources in the upper source band [128N',256N'). It does not enumerate that
whole annulus. It scans743154 primitive candidates until that many witnesses
are found, then verifies each witness by a literal111-step stopped prefix.
Every witness retains source, mod15, mod35, endpoint after111 steps, actual
accumulated clock, and original charge. All4096 have norm-one charge. Their
axis totals and exact signed-zonotope energy are

    c_subset=(1371,0,1344,0,1381,0),
    D_subset=7499929.

Adding further source segments cannot decrease the signed-zonotope maximum,
so the complete pending set satisfies

    D(N',111)>=7499929>5880064=64D(4096,105).

Likewise m(N',111)>=4096>3088=8m(4096,105). Already three axis counts of the
subset exceed8c. Hence the naive rules

    D(64N,h+6)<=64D(N,h),
    m(64N,h+6)<=8m(N,h)

are genuinely false in the intended N>=4096 range. This does not say that
N'=262144 violates153N', or that no later cut h' can satisfy the scale cone.
No complete energy or pending census at262144 is claimed. Six extra
**primitive steps** here are merely a generous test choice; the +6 in the
native lift below is a division-clock change, not an arrival-depth change.

## 3. Why the apparent scale lift loses almost all new sources

For a primitive source s, put a=2 when s=1 mod8 and a=1 otherwise, and let
P(s)=u have actual clock t=a+2j. The exact inverse formula is

    9s=2^t(3u+1)-(2^a+3).

Therefore

    L(s)=64s+7(2^a+3)

has P(L(s))=u and edge clock t+6. For source and lifted source both outside
V, the stopped primitive depth is the same, and the full arrival clock
increases by6. It is not a transport adding six primitive arrival steps.
If the original source itself belongs to V, stopped-depth equality need
not hold and the frontier membership sidecar must be retained.

Moreover the image of positive primitive sources has exactly the three
residue classes

    L(s) mod512 in {113,227,483}.

For s=1,3,7 mod8 these are respectively113,227,483. Conversely every positive
integer in one of these classes has the corresponding positive primitive
preimage under the displayed affine formula. Thus the image has density
3/512, while primitive sources have density3/8: the lift covers only1/64 of
the primitive sources. This is an exact periodic density statement; finite
annulus endpoint clipping is not silently ignored. One old source address
cannot represent the other63/64 of the enlarged source population.

The source guards transform as L(s) mod105=(64s+35) mod105 in branch1 or
(64s+49) mod105 in branch2. The full clock, stop membership, source-band
boundaries and missing addresses cannot be recovered from those guards
alone. Thus the familiar inverse-address lift is useful but does not prove
any of the population/coherence inequalities(2).

## 4. Guard, charge and actual edge clock do not close next arrival

At N=4096, consider two actual new sources

    x=32881, y=247921.

Both are outside V, both have charge(-1,-1), both have joint source guards
(mod15,mod35)=(1,16), and both have first-edge clock8. Nevertheless

    P(x)=385 in V,
    P(y)=2905 not in V.

They have different first-arrival status at the next primitive step. A
one-source subpacket made from either has the same cut-zero six-axis state,
zero terminal vectors, and identical listed guard/clock data, but different
next-cut states. This proves nonclosure for the projected native subpacket
update. It does not claim impossibility of every finite-dimensional model,
or pretend these are two different full-annulus universes at the same N.
The missing coordinate is actual magnitude and its grounded-frontier
membership, not another independent copy of the shared mod5 guard.

### General exact prefix cylinder

Let s=s_0->s_1->...->s_r be any finite native P path, with individual clocks
t_1,...,t_r and T=sum t_i. For any integer k>=0 define

    s_0'=s_0+840 k 2^T.

Then its first r native clocks are identical and, with T_j=sum_(i<=j)t_i,

    s_j'-s_j=840 k 3^j 2^(T-T_j),  0<=j<=r.        (5)

Proof: the inverse formula shows that a source difference divisible by2^t
induces a parent difference3 times as large divided by2^t. The displayed
shift is a multiple of8 at every endpoint, so the target remains primitive;
the same branch and exact number of B-strippings remain valid. Inducting
gives(5). Each shift is also divisible by105, so every corresponding prefix
state has the same joint mod15/mod35 registers. The clocks are unchanged as
integers, not merely modulo anything. Finite membership in V, however, need
not agree. Formula(5) concerns the unstopped native path; stopping is exactly
the extra predicate lost under the projection.

This explains the small hostile and supplies a structural reason why a
finite guard/clock description alone cannot decide grounded arrival. Keeping
source identities or magnitude-sensitive suffix certificates repairs the
loss; it does not by itself yield a uniform population bound.

## 5. Connection contract and reproducibility

Source -> target: the actual stopped ledger maps to six-axis pending
segments and terminal supports. Preserved predicate: the exact sufficient
suffix bound(1). Destroyed information: which actual source reaches which
frontier point at which next cut. Required sidecar: literal source values,
first-entry membership, clocks and prefix records. Cheapest decisive tests:
the two-source hostile, then the explicit4096-source subset disproving the
most direct scale contraction. This is an obstruction plus a repaired
conditional propagation statement, not merely a larger sweep.

The checker imports only the explicitly selected unchanged base checker,
authenticates its frontier and base receipt, checks the exact lift and
six-axis identities, and stores the4096 literal pending witnesses. Each
witness is independently walked through the claimed cutoff rather than
trusted from the memoized search. Normal and optimized immutable replays use
non-removable checks.

    python3 check.py --repo-root /path/to/math --base-dir /path/to/base-packet --base-receipt /path/to/base-packet/receipt.json --output receipt.json
    python3 -O check.py --repo-root /path/to/math --base-dir /path/to/base-packet --base-receipt /path/to/base-packet/receipt.json --verify receipt.json

Stop here: the coherent outcome is the explicit scale cone and a structural,
finite-exact obstruction to obtaining it from naive guard/clock scale
transport. A new proof must adapt h', exploit more than those projected
coordinates, or use a weaker aggregate inequality than coordinatewise(2).
