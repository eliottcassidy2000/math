# Triangular-valuation Newton coordinates for exact clone-excess residues

2026-10-26. **PROVED:** at fixed absolute binary precision K, the integral
moment-bank observation image has an explicit mixed-precision Newton
presentation. Its full polynomial image requires only

    t(S,K)=min(S,L(K)),
    L(K)=#{j>=0 : j*(j+1)/2 < K}

cyclic residue coordinates, with respective precisions
`K, K-1, K-3, K-6, ...` while positive. **IMPLEMENTED / FINITE-EXACT:** an
adaptive sparse/Newton transfer preserves every allowed signed clone-excess
batch and all state mergers. This is a genuine reduction beyond the inherited
`min(h_support,S,K)` constant-tail weight representation. It does not improve
the already scalar ordinary A000568 state, and does not provide normalized
binary digits after division by a power of two.

Dependencies are the integral gauge in
`05-knowledge/results/port_residual_integral_20261025/counting/README.md`,
the absolute-precision constant-tail quotient in that directory's `MODULAR.md`,
and the exact clone-cut/excess semantics in
`05-knowledge/results/ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md`.
No new literal clone merger, owner, or port quotient is introduced.

## 1. Inheritance and the smaller carrier

Closest mechanism: the proved integral moment transfer. Canonical hostile:
a full integer product lattice need not satisfy the moment congruences; even
the constant-tail module is larger than the actual polynomial evaluation
image. Corrected near miss: storing K channels each at K-bit precision does
not establish a minimum-size modular module. Least-used sidecar: after the
integral gauge, every cycle-count exponent is a **nonnegative integer**.
The concept board is integer bank / future evaluation / triangular valuation /
signed exact excess / semantic support / absolute versus normalized precision.

At a nonterminal `(S,D)` state put `b=b(D)` and write

    P(X)=sum_m x_m X^(m-b),     x_m=O(u)*w_m in Z,    m>=b.

Its observations are exactly the integral moments

    Y_r=P(2^r) mod 2^K,        1<=r<=S.

The coefficient lattice here is an integer polynomial lattice, not arbitrary
independent Y coordinates. The constant-tail theorem used some consequences
of that distinction; this result uses the full evaluation image.

## 2. Exact invariant factors of the observation image

Define

    q_0(X)=1,
    q_j(X)=product_(i=1)^j (X-2^i),
    tau_j=j*(j+1)/2.

The q_j are monic of degree j, so they are an integral unimodular basis of
Z[X]. Every integer polynomial has a unique finite expansion
`P=sum_j c_j q_j` with integers c_j. At a positive evaluation node,

    q_j(2^r)=0                        if r<=j,
    v2(q_j(2^r))=sum_(i=1)^j i=tau_j if r>j.       (N1)

Indeed `2^r-2^i=2^i*(2^(r-i)-1)`, with the last factor odd when r>i.
Consequently q_j is invisible on the horizon if j>=S or tau_j>=K.
For the remaining j, only

    c_j mod 2^(K-tau_j)                              (N2)

is needed.

These conditions describe the entire kernel, not merely sufficient
vanishing tests. Suppose P(2^r)=0 modulo 2^K for all r=1..S. Evaluating at
r=1 first gives `c_0=0 mod 2^K`. Inductively assume each earlier c_i is
divisible by `2^max(K-tau_i,0)`. By (N1) those earlier terms vanish at
every later node. At r=j+1 the new coefficient is multiplied by a number
of valuation exactly tau_j, so it must satisfy (N2) with residue zero.
This proves necessity; (N1) proves sufficiency.

Equivalently, divide column j of the evaluation matrix by `2^tau_j`.
Its first j entries are zero and its next entry is odd, leaving a lower
triangular matrix with unit diagonal entries on the first
min(degree slots,S) rows. Thus the image,
as a `Z/(2^K)` module, is exactly

    direct_sum_(0<=j<S, tau_j<K) Z/(2^(K-tau_j)).    (N3)

For polynomials of degree less than h_deg, replace j<S by
`j<min(h_deg,S)`. The exact minimum number of module generators is
`min(h_deg,S,L(K))`; the exact information size in bits is

    sum_(0<=j<min(h_deg,S)) max(K-tau_j,0).          (N4)

This is a theorem about the specified full integer coefficient lattice
`Z^h_deg`. It is not a claim that all such coefficients are realized by one
particular counting run, or a minimum bit encoding for every reachable fibre.

For S>=64 and K=64, there are 11 coordinates of precisions
`64,63,61,58,54,49,43,36,28,19,9`, totalling 484 residue bits. The old
constant-tail presentation used 64 coordinates of 64 bits, totalling 4096.
These are mathematical residue payload sizes, not Python-object memory.

### Degree slots and arbitrary sparse support must not be confused

Let h_support be the number of distinct semantic m keys in a sparse bank.
They need not be consecutive. Adaptive storage below uses at most
`min(h_support,S,L(K))` scalar coordinates. This is an upper bound, not an
exact invariant-factor formula for every sparse exponent set. The hostile
`P(X)=x*X^K` has h_support=1, but all its positive-node observations vanish
modulo 2^K. Sparse exponent gaps and additional actual-prefix divisibility
can therefore reduce the observation image further. Formula (N3), with its
h_deg refinement, is the exact result on the explicitly stated lattice.

## 3. All legal integral transitions descend

For a batch of k actual ell-cycles and e excess units, set

    a=ell*k,  S'=S-ell*(k+e),
    D'(d)=D(d)+k*(d-gcd(ell,d)),
    delta=k+b-b' >= 0,
    J0=A_ell(k)-v2(k!)+k*(ell*b-D(ell)) >= 0,
    T=beta(k,e)*O(u')/[O(u)*ell^k*O(k)] in Z.

These inequalities and integrality are inherited from the proved integral
gauge, including coordinate deletion, skips and excess. The polynomial
transfer is simply

    P'(X)=T*2^J0*X^delta*P(2^a X).                 (N5)

At a destination node `X=2^r`, it reads the source node `2^(a+r)`.
For S'>0, `1<=a+r<=S`. Therefore any source kernel polynomial maps to a
destination kernel polynomial: its source observation is zero modulo 2^K,
and the only further operations multiply by integers and nonnegative powers
of two. This proves closure under the complete kernel from section 2,
not only under an upper-bound representation.

When S'=0<S, a completing batch necessarily has k>=1. Its terminal value
is `T*2^J0*P(2^a)`, again an available positive source node. The factor
`X^delta` becomes 1. An already terminal state is a separate scalar and
carries that scalar unchanged. No unavailable value P(1) is reconstructed
from the source quotient.

The sign in beta is retained exactly. Source and destination b are functions
of their exact keys, and the map is linear, so arbitrary additions, sparse
conversions and incoming state mergers commute. No canonical-owner or clone
identity is guessed from the polynomial: the inherited signed cycle-index
and clone-cut interpretation remains the semantic justification of the count.

### No division-by-two loophole

All operations in (N5) have J0,delta>=0. No moment normalization with a negative
power of two is performed after reduction. Dividing by two does not preserve
the kernel: at K=2 the constant P=4 is invisible, while P/2=2 is visible.
Likewise P=X^2 is invisible but P/X=X is visible. The latter is a hostile to
a hypothetical negative delta. These operations are excluded, not repaired
by the Newton presentation.

This theorem computes absolute residues modulo 2^K. It cannot by itself
compute `a(N)/2^B mod 2^K`; the reverse potential in the separate ordinary
implementation changes coordinates before truncation and addresses a
different question.

## 4. Auditable arithmetic and adaptive support

The recurrence can operate entirely in mixed Newton residues without any
division. Use

    X*q_j=q_(j+1)+2^(j+1)*q_j.                     (N6)

To form q_j(2^a X), repeatedly multiply by `(2^a X-2^i)` in the target
Newton module. Combine these basis polynomials with the source c_j, multiply
by X^delta using (N6), then by T*2^J0. Discard exactly the invisible basis
coordinates and reduce retained coefficients at their own precisions.
Every operation is integer multiplication, addition, nonnegative shift, or
an odd-unit inverse for T. Representatives of the smaller coefficient
residue classes are safe because the source kernel descends by section 3.
Python temporaries can be wider than their stored residues; no machine-word
width bound is asserted.

Keep sparse m->x_m while h_support<=t(S,K), then convert to Newton form.
Retain each semantic sparse m key even when its coefficient is zero modulo
2^K, and retain zero nonterminal dense vectors. The inherited exact sign
law `(-1)^(u-m)` implies that same-m integer contributions cannot cancel.
Every batch shifts support injectively by m->m+k. Since t(S',K)<=t(S,K),
once dense, a successor cannot need sparse reconstruction. Thus the adaptive
implementation uses `min(h_support,S,L(K))` scalar coordinates per
nonterminal defect fibre, with mixed precision in its dense form. Semantic
keys, graph metadata and representation tags are additional storage.

Strict compression occurs on genuinely admissible prefix families, not just
invented exponent slots. Before stage ell=3, process four total 9-cycle units:
k=1,2,3,4 actual cycles and e=4-k excess units. All have u=36, D(3)=0 and
m=k. At target N=40, S=4. The old K=4 constant-tail bank has four available
weight slots, while L(4)=3 suffices in the Newton quotient. This family
retains the exact excess coefficients; it does not discard the clone merger.

## 5. Reproduction, exact controls and cost boundaries

    python newton_modular.py > newton-normal.json
    python -O newton_modular.py > newton-optimized.json
    cmp newton-normal.json newton-optimized.json

Supply `--repo /path/to/repository` if not running in the research repository
or the session workspace. All gates use explicit exceptions.

The checker compares every pair-free/core coefficient N=0..44 at every
absolute precision K=1..9 with both the inherited exact integer engine and
the previous constant-tail modular engine. Every intermediate sparse/Newton
bank is also compared with the independent rational f-response recurrence,
scaled by its exact odd factorial, not with another Newton transfer. Direct
literal pair-free class enumeration checks N=0..5. The final receipts record
all counts and stagewise payload comparisons.

Additional controls cover arbitrary signed integer polynomials, every legal
decreasing source/destination evaluation horizon in the declared finite
universe, exact kernel generators, exhaustive mixed-residue injectivity,
exhaustive bounded-degree image cardinalities, semantic zero support, the
sparse-gap rank hostile, both invalid divisions, and invalid inputs.

The ordinary one-channel recurrence is unchanged. The `(S,D)` key graph is
also unchanged in core mode: the gain is in weight-coordinate storage and
stored residue payload. The reference Newton transform may perform more
operations than the old direct evaluation-vector transfer; it includes
basis conversion work. Neither elapsed runtime, peak memory nor end-to-end
bit complexity is established by the coordinate count. Sparse sets may
admit still smaller images, but deriving their exact adaptive invariant
factors is separate work.
