# Native ROOT boundary modes: domain classification, cutoff obstruction, and the marked C3 lattice

Session 75, 2026-10-27. **Phase 1 frozen packet.**
**PROVED scoped deductions** and **FINITE-EXACT controls**.
No all-source ROOT arrival, Collatz convergence, uniform source-debt estimate,
or all-scale spatial ROOT-current bound is claimed.

## Findings

1. A homogeneous first-hit boundary mode is a genuine formal solution after the
   ROOT self-edge has been removed. Its root-basin values are uniquely forced:
   `h(s)=epsilon(s) 4*2^T/(3^d*(3s+1))`, for `h(1)=1`. It is ROOT-conditioned;
   setting it to zero outside the ROOT basin is compatible with every interior
   equation and does not establish membership of any omitted integer.
2. Its membership in a specified all-source weighted Hilbert space is exactly a
   summability question, stated below. The recovered native source documents do
   not supply one fixed all-source Hilbert weight for this mode. In counting
   `ell2`, even ROOT's direct inverse family rules out a nonzero boundary mode.
3. For the actual finite-cutoff weights `mu_N(s)=Q(A_N(s))`, the mode has finite
   seminorm automatically. **For every N>=2, however, native pullback and its
   forcing difference do not descend to the quotient `ell2(mu_N)`.** An explicit
   source of positive mass maps outside the mass support. This is a domain
   obstruction, not a divergent finite sum.
4. Zero interior forcing is not free transport of the prescribed source. To
   produce its unchanged terminal charge, the initial boundary energy must be
   inflated by the exact factor `B_s^(-2)`, where `B_s=product b` on its path.
   Keeping original initial energy instead attenuates the terminal charge.
5. The initial boundary invoice has a new universal **polynomial-in-depth**
   bound on first-hit words: `B_s^(-2)<=(1+3d)^(2/9)`. On an h-edge
   nonrepeating exterior prefix of the certified 627-frontier, the sharper
   factor is `(1+3h/754)^(2/9)` and paid forcing grows at most logarithmically.
   A visit-multiplicity sidecar covers arbitrary prefixes without assuming
   that unknown positive cycles do not exist.
6. The actual three ROOT rows embed in an **index-three sublattice** of the
   integer centered plane. The directed C3 skew map is exactly multiplication
   by `1+2 omega`, with square `-3`. It preserves source labels and dyadic pair
   ownership, but multiplies both energy and debt by three: it does not pay debt.
   The first row remains a marked observer, not a rotational average.

## 0. Inheritance and boundaries

All recovered sources are pinned at commit
`96e35642110104cc547516a63f37eebcf930f707` in `SOURCE_PINS.json`.
The closest mechanism is `paid_root_cycle_budget_20261026/README.md`, R1--R3:
positive native affine carry, its two-port identity, and dyadic source ownership.
The source/row definitions are E1--E2 in
`root_eisenstein_first_arrival_20261018j.md`. The sparse-boundary note R5--R6
supplies the distinction between first hit, moving sources, and ROOT support.

The finite-norm boundary-mode sentence is inherited transcript evidence; no
additional source defining its fixed Hilbert measure/domain has been located.
The statements below therefore separate every possible reading. In particular,
**the S73 fully odd map U and its probability lambda/occupation W are not used**.

Concept board: native edge / boundary domain / source norm / marked observer /
pair debt. Canonical hostile: a positive-mass source whose target has zero mass.
Corrected near miss: finite norm on a cutoff is not a closed operator domain.
Least-used sidecars: direct inverse families and the index-three row-lattice
condition. Anchor: boundary classification; niche: quotient-domain obstruction;
wildcard: C3's exact integral observer rather than its spectrum alone.

## 1. Exact native definitions

Let `S={s>0: s odd, s !=5 mod8}` and `B(x)=4x+1`. For s in S, put

```
a(s)=2 if s=1 mod8, else 1,
(3s+1)/2^a = B^j(y),  y in S,
P(s)=y,  t(s)=a+2j,  epsilon_s=(-1)^t,
X_s=3s+1,
b_s=3X_s/(3X_s+2^a),  d_s=1-b_s^2,  c_s=epsilon_s b_s.
```

The B-chain is stripped while its odd member is 5 mod8; this is the unique
primitive representative. Direct substitution gives the corrected identity

```
2^t X_(P(s)) = 3X_s+2^a.
```

The correction by three units is real: at s=3, y=1, a=1, t=3, it reads
`32=30+2`, not `32=30+5`. All `0<b_s<1` and `d_s>0`.

For any real or Eisenstein-vector values e,F, with
`e'=c_s e+F` and `r=d_s e-c_s F`, polarization gives

```
Q(e')+Q(r)/d_s = Q(e)+Q(F)/d_s.                    (1)
```

This signed-gain version is exactly the same scalar contraction dilation as
R1, with the actual edge sign retained. It is a path identity, not by itself a
spatial operator on a prescribed measure.

## 2. First-hit homogeneous boundary mode

Let R be the actual ROOT basin: the seeds whose forward P orbit first reaches
1 in finitely many steps. The definition is mathematical and makes no claim
that `R=S`. Write `d(s),T(s),epsilon(s)` for the first-hit native depth, complete
integer clock, and sign, only when `s in R`. ROOT has zero first-hit depth and
clock. Its outgoing self-edge is not executed.

For a first-hit word `s=s_0,...,s_d=1`, the affine identity telescopes:

```
B_s := product_(i<d) b_(s_i) = 3^d X_s/(4*2^T),
product_(i<d) c_(s_i) = epsilon(s) B_s.            (2)
```

### Proposition B1 — precise rooted solution

For a prescribed scalar or Eisenstein boundary value z, the unique function on
R satisfying

```
h(1)=z,
h(P(s))=c_s h(s) for s in R\{1}
```

is

```
h(s)=epsilon(s) B_s^(-1) z
    =epsilon(s) 4*2^T/(3^d X_s) z.                (3)
```

**Proof.** Multiply the nonzero gains along the unique first-hit word and solve
for h(s). Conversely, removing its first edge from that word checks the local
equation. No limit or arrival assertion about another source is used. QED.

Extending h by zero on `S\R` satisfies every remaining interior equation:
`P(s) in R` iff `s in R`. Thus a nonzero boundary solution supported on R cannot
prove a positive value at an arbitrary primitive integer.

If the equation is also imposed at the ROOT self-edge, then `P(1)=1`,
`a=t=2`, `c_1=b_1=3/4`; hence `h(1)=3h(1)/4` forces `h(1)=0`.
**A nonzero boundary mode must explicitly remove this edge or pay its forcing.**
For h(1)=z the omitted self-edge would have `F=z/4`, costing `Q(z)/7` in (1).
This is boundary forcing, not interior forcing.

### Proposition B2 — exact Hilbert criterion and counting no-go

For any supplied strictly positive weights w(s), the ROOT-supported solution
(3) belongs to `ell2(S,w)` iff

```
Q(z) * K_w < infinity,
K_w = sum_(s in R) w(s) B_s^(-2).                 (4)
```

This is an iff, not a claimed bound on the unprovided weights. In particular
there is **no nonzero counting-ell2 mode**, even after ROOT's self-edge is
removed.

**Proof of the no-go.** For all k>=0 the actual primitive seeds

```
r_k=(32*64^k-5)/9 = 3,227,14563,932067,...
P(r_k)=1,  a=1,  t=3+6k,
b_(r_k)=1-1/(16*64^k).
```

These assertions follow either by direct substitution into P or by the native
inverse compiler. They are distinct and positive. Equation (3) therefore gives
`h(r_k)=-z/b_(r_k)`, so `Q(h(r_k))>=Q(z)`. The counting sum diverges. The same
argument excludes every w with `sum_k w(r_k)=infinity`. QED.

**General classification, if a fixed w is supplied.** Delete ROOT's outgoing
edge and use the literal equations on the resulting functional graph. Its
ROOT component contributes the one-dimensional mode (3) iff K_w is finite.
Every other component containing a finite directed cycle contributes zero,
because the modulus of the product of its gains is strictly less than one.
On a component with an injective infinite forward ray, choose a base o. If
`P^i(s)=P^j(o)`, define

```
H_C(s)=product_(k<j)c_(P^k(o)) / product_(k<i)c_(P^k(s)).
```

This is independent of the chosen meeting point: extending both paths
multiplies numerator and denominator by the same gains. Any two vertices of a
weak component have coalescing forward paths, so H_C is defined throughout it.
Every formal solution is `z_C H_C`. It is admissible in ell2(w) exactly when
`sum_C Q(z_C) sum_(s in C) w(s)|H_C(s)|^2<infinity`; components with divergent
inner sum must have z_C=0. This completely classifies homogeneous solutions
once the graph, first-hit boundary domain, and positive weights are fixed. It
neither asserts the existence of nonroot components nor decides ROOT arrival.

## 3. The original finite-cutoff measure is not an operator domain

For integer N>=1 the genuine source charge is

```
A_N(s)=sum_(j:B^j(s)<4N) omega^(B^j(s)) *
       [(omega-omega^2) if B^j(s)<2N, else omega],
mu_N(s)=Q(A_N(s)) in {0,1,3}.
```

In particular `mu_N(s)=0` when `s>=4N`. The expression
`||h||_N^2=sum_s mu_N(s)Q(h(s))` is a seminorm on literal functions and a norm
only after zero-mass coordinates are identified. Every ROOT-supported mode
(3) has finite seminorm here, simply because the support of mu_N is finite.
That finiteness is automatic; it supplies no all-scale or all-source bound.

### Proposition H1 — all-cutoff quotient obstruction

For **every N>=2**, the maps `Ch(s)=h(P(s))` and
`Fh(s)=h(P(s))-c_s h(s)` do not descend to the quotient `ell2(mu_N)`.

**Proof.** Choose

```
s_N=4N-1 if N !=1 mod4,
s_N=4N-5 if N=1 mod4 (so N>=5).
```

In both cases s_N is positive, primitive, and in `[2N,4N)`. Its B-chain has
only that point below 4N, so `A_N(s_N)=omega^(s_N+1)` and `mu_N(s_N)=1`.
The native branch has a=1,j=0 and

```
y_N=P(s_N)=6N-1 or 6N-7, respectively.
```

For the first case y_N is not5 mod8 precisely because `N !=1 mod4`; for the
second case y_N is7 mod8. In both cases `y_N>4N`, so `mu_N(y_N)=0`.
The function `delta_(y_N)` represents zero in ell2(mu_N), but
`C delta_(y_N)(s_N)=F delta_(y_N)(s_N)=1`. Hence neither map is well-defined
on equivalence classes. QED.

The minimal example is `N=2,s=7,P(s)=11`: `mu_2(7)=1`, `mu_2(11)=0`.
The boundary case N=1 is different: the positive support is `{1,3}`, which is
forward closed. Thus the stated threshold N>=2 is exact.

More generally a composition operator on a discrete weighted quotient needs
`w(P(s))>0` whenever `w(s)>0`. Its boundedness, once this is satisfied, is
exactly the finiteness of
`sup_(y:w(y)>0) [sum_(P(s)=y) w(s)]/w(y)`. This follows by regrouping the
squared norm, and shows why finite support alone is not enough.

Repairs are legitimate but must be declared: a source-labelled path space,
a forward-closed finite authenticated carrier with explicit extra weights,
or another supplied all-source measure. Such repairs change the domain or
norm. Zero-extending representatives merely selects a section; it does not
make the failed quotient operator intrinsic. Nothing here licenses replacing
mu_N by S73's lambda or W.

## 4. Zero interior forcing moves the invoice to the initial boundary

Fix one authenticated first-hit path and its prescribed initial source charge
A, with `mu=Q(A)`. The actual signed charge transport has terminal
`epsilon(s) A`.

If F=0 and the initial value remains A, (2) gives instead
`e_terminal=epsilon(s) B_s A`. Unless A=0 or the path is empty, `B_s<1`, so this
is a different current. To obtain the actual terminal value with F=0 one must
set `e_initial=A/B_s`. Telescoping (1) then gives the exact invoice

```
initial energy = mu/B_s^2,
terminal energy = mu,
dissipative second-port energy = mu*(B_s^(-2)-1).         (5)
```

Thus zero interior forcing does not mean zero debt. The new initial condition
owes `sum_s mu_s(B_s^(-2)-1)` in addition to the original source norm. Neither
its uniform bound nor coherent terminal cross terms are paid by existence of
this mode. At source 3 already `B_3=15/16`, so unit terminal charge costs
`256/225`, not 1, at the initial boundary. Truncating the word gives the same
identity with the exact prefix product and its actual endpoint; no ROOT
boundary value may be silently imposed at an unarrived endpoint.

The original R1 forcing invoice and (5) are alternative representations of
fixed data, not permission to alter the prescribed source family. Choosing an
input that solves a homogeneous equation is not a bound for the actual input.


## 4a. A positive bound: polynomial boundary inflation and logarithmic forcing

### Proposition B3 — sorted native losses and the no-repeat improvement

At every nonroot primitive seed define `q_s=1/lambda_s`, where
`lambda_s=(1-b_s)/(1+b_s)`. The full multiset of these denominators, sorted
in increasing order, is exactly

```
q_k=18k+10+3(-1)^(k+1),  k>=1:
31,43,67,79,103,115,... .
```

To verify this without enumeration, the class `s=8j+1` gives `q=36j+7`;
excluding ROOT leaves `43,79,115,...`. The other two primitive classes
`s=8j+3,8j+7` together give `31,67,103,...`. These two disjoint progressions
interlace as stated. Therefore, for any h distinct nonroot native edge sources,

```
product b_s^(-1) <= product_(k=1)^h (q_k+1)/(q_k-1)
                <= (1+3h)^(1/9).                        (9)
```

**Proof.** Reciprocal gains decrease with q, so the first inequality is the
finite rearrangement maximum. Since `q_k-1>=18k+6`,

```
log product (q_k+1)/(q_k-1)
 <= sum_(k=1)^h 2/(q_k-1)
 <= sum_(k=1)^h 2/(18k+6)
 <= integral_0^h 2/(18x+6) dx = (1/9) log(1+3h).
```

Here the logarithmic bound is the elementary `log(1+x)<=x`; the last bound
uses a decreasing integrand. A first-hit ROOT word has no repeated seed,
because a repetition in a deterministic orbit would prevent later first
arrival. Hence its initial energy inflation obeys

```
B_s^(-2)<=(1+3d(s))^(2/9).                               (10)
```

The ranked product is optimal over unconstrained h-element selections of
nonroot sites; no optimality over *legal native path words* is asserted.
This is a depth-dependent bound, not uniform in depth and not an arrival proof.

### The actual 627-frontier improves the constants

Let V be the unchanged forward-closed frontier in `sources/frontier.txt`.
Its finite descent ranks are independently authenticated by the checker.
The minimum exterior primitive seed is 1009. Consequently every exterior
reciprocal loss denominator is at least 4543. The full denominator multiset
above, starting at 4543, is

```
4543,4567,4579,4603,...,
```

whose kth value is at least `4543+18(k-1)`. Exterior sites form a subset of
this tail; the tail includes some interior sites, which only weakens the bound.
For any h-edge exterior prefix with no repeated seed, exactly the same proof
and a direct sum of reciprocal denominators yield

```
B_ext^(-2) <= (1+3h/754)^(2/9),                          (11)
sum_edges lambda_s <= (1/18) log(1+18h/4525).             (12)
```

For (11), replace `q_k-1` by `4524+18k` in the integral. For (12), use
`q_k>=4525+18k` and integrate `1/(4525+18x)` from 0 to h.
No transcendental approximation is used in the proof.

For any finite collection of original charged sources whose retained exterior
prefixes have no repetitions and native depth at most h, (12) gives

```
F_N <= (S_N/18) log(1+18h/4525).                         (13)
```

Every successful first-entry word has no repeats. A prefix ending before
entry need not have this property globally, and this hypothesis is not dropped.
For an arbitrary finite exterior prefix, let m be its number of distinct sites
and M its largest visit multiplicity. Group visits by site before applying the
same inequalities to obtain the fully unconditional finite-prefix repair:

```
sum_edges lambda_s <= (M/18) log(1+18m/4525),
B_ext^(-2) <= (1+3m/754)^(2M/9).                         (14)
```

Source weights stay attached to their original source when (14) is summed.
Thus moving sources and potential repetitions have not been erased. The M
sidecar is exactly why this is not a uniform all-depth forcing payment.

Finally the finite internal controller supplies an exact boundary coefficient.
The maximum internal `B_u^(-2)`, for u in V, is

```
C_V = 332306998946228968225951765070086144
      /211753373595316860199691410086465225
    < 1.57,
```

attained at u=993, native rank31; there
`B_u=460166680231540515/576460752303423488`. Forward closure and strict descent
certify all 627 internal words. If an actual source first enters V after h
exterior edges, (11) gives the full ROOT boundary invoice

```
B_s^(-2) <= C_V (1+3h/754)^(2/9).                        (15)
```

This finite constant plus a mild explicit depth tariff is a proved survivor
of the zero-forcing idea. It is neither the original unchanged source energy
nor a bound for coherent source merging.

## 5. Directed C3, the exact integral lattice, and its marked observer

Let `omega^2+omega+1=0`, with positive imaginary part, and
`Q(a+b omega)=a^2-ab+b^2`. Define

```
E(a+b omega)=(2a-b,-a-b,-a+2b),
R(x0,x1,x2)=(x1,x2,x0),
K=R-R^2 = [[0,1,-1],[-1,0,1],[1,-1,0]].
```

The map E is the actual source-row map, not an arbitrary isospectral change:
`E(A_N(s))=(c_N^(0)(s),c_N^(1)(s),c_N^(2)(s))`. Direct calculation proves

```
E(omega z)=R E(z),
E((1+2omega)z)=K E(z),
(1+2omega)^2=-3,  K^2=-3I on x0+x1+x2=0,
||E(z)||_2^2=6Q(z).
```

The integral image is precisely

```
L={x in Z^3: x0+x1+x2=0, x0=x1=x2 mod3}.
```

Indeed the inverse is `a=(x0-x1)/3`, `b=-(x0+2x1)/3`. Relative to the usual
basis `(1,-1,0),(0,1,-1)` of A2, E(1),E(omega) have determinant -3, so L has
index three. For example `(1,-1,0)` is centered but is not in L. Calling E an
integral identification with all A2 would lose this congruence sidecar.

The original physical observer is

```
ell_0(a+b omega)=2a-b,   ell_0(z)^2+3b^2=4Q(z).            (6)
```

It is sharp at b=0. Rotating z but keeping ell_0 fixed changes the target:
`ell_0(1)=2`, `ell_0(omega)=-1`. After K the conjugated observer must be
`ell_0 composed with K^(-1)`; on the plane `K^(-1)=-K/3`.
No unmarked spectrum identifies these readouts.

### Exact transfer to source-debt ledgers, and its limit

For `beta(z,w)=Q(z+w)-Q(z)-Q(w)`,

```
E(z) dot E(w)=3 beta(z,w),
beta((1+2omega)z,(1+2omega)w)=3 beta(z,w).                (7)
```

Therefore applying K to **every named charge** commutes with native signs,
terminal aggregation, and source partition refinement. Every E_q, every
signed level debt Delta_q, source norm S, and prescribed positive-after-level-
sum debt D^+ scales by exactly 3. Ratios such as D^+/S are unchanged. The map
preserves all pair owners but supplies no cancellation estimate.

There is a useful anisotropic identity instead. Write every aggregate
`W=a+b omega`, and let `B_q=sum b(W)^2` and `M_q=sum ell_0(W)^2` over exactly
the same bins as E_q. Equation (6) yields

```
M_q=4 E_q-3 B_q,
Delta_q^marked=4 Delta_q^Q-3 Delta_q^b.                  (8)
```

Thus a marked-source estimate can exploit a certified missing-leg covariance;
it cannot erase that covariance. An **authenticated native hostile** at N=512 has source 799 with
`A_N(799)=-2-omega` (already at terminal 799, clock 0) and source 1065 with
`A_N(1065)=omega`, `P(1065)=799`, clock 2. Source 799 belongs to V and 1065
does not. Thus their signed charges are unchanged and share the actual first-
arrival terminal 799, yet `beta_Q=0` while `2ell_0(z)ell_0(w)=6`. Their odd-
source coordinates 399 and 532 differ by 133, so the dyadic pair owner is q=1.
This exact witness was found by the dyadic-ledger worker and independently
verified here. It distinguishes an exact geometric transfer from a new
source-specific debt theorem; other pairs may still cancel at the same level.

## 6. Verification and stopping boundary

Run inside this packet directory with Python 3 and no third-party packages:

```
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json
```

The checker has explicit failure checks, not assertions disabled by `-O`.
It independently compares source charges with the original three K rows for
all N=1..128; tests every cutoff hostile N=2..8192; checks 16 exact direct-ROOT
inverse witnesses; authenticates and checks the full first-hit carry/product/
boundary-invoice identities for every positive primitive seed <=20000;
authenticates all 627 frontier ranks and its exact internal maximum; checks
the polynomial depth/exterior bounds and actual marked pair; and
checks the C3 lattice, observer, polarization, and debt scaling on declared
integer grids. The finite first-hit universe is not a universal arrival claim.

The bounded question is resolved as follows: formal first-hit boundary modes
are classified, finite-cutoff seminorm existence is automatic but does not
provide an operator domain, and counting-ell2 admissibility is decisively
excluded. The depth-dependent invoice (15) is the positive quantitative result.
A different fixed native all-source weight, if supplied, must be judged by (4),
not guessed. Primitive-source debt, arbitrary-source arrival,
and the actual spatial source-covariance estimate remain open.
