# Sunflower supports, a common horizon space, and source ownership

20 October 2026. **PROVED:** exact overlap and finite-horizon observer
identities. **FINITE-EXACT:** native source receipts and interpolation controls.
**CITED INPUT:** the supplied sunflower preprint's stated common-space lemma
and interpolation mechanism. The paper's main theorem is not re-audited here.
The unbounded native 153N inequality and universal entry remain **OPEN**.

Artifacts: [checker](verify.py), [exact receipt](receipt.json).

## 1. Inheritance and the useful paper mechanism

The source is Bucić–He–Lee, *Proof of the Erdős–Szemerédi sunflower
conjecture*, arXiv:2610.26686v1, dated 19 October 2026:
[primary record](https://arxiv.org/abs/2610.26686) and
[the supplied PDF](C:/Users/Eliott/Downloads/2610.26686v1.pdf).
The relevant complete pages were read; pages4,7,9 were also rendered and
visually inspected.

The precise interfaces are:

- pp3–4: freedom to choose nonzero values on rainbow tuples while keeping
  the allowed equality-pattern support;
- Lemma3, pp4–5: a tensor diagonal on V, decomposed with all one-variable
  factors from one common space P, forces |V|≤dim P;
- Lemma6, pp6–8: successive interpolation puts degree≤h−1 in the h-th
  selected variable, leaving an error supported exactly on rainbow prefixes;
- Lemma7, pp8–9: the coordinate product and a pigeonhole count produce
  a common monomial space. This is an algebraic counting mechanism;
  it does not assume statistically independent samples.

No numerical sunflower bound is needed below. The transferable operation
is to identify **one exact common space for all requested observers**.
The paper's ability to tune nonzero amplitudes is different from permission
to alter physical signed source weights.

The nearest native mechanism is
[the exact arrival-cut ledger](../session18_native_annulus/README.md),
especially its equations(A)–(E). It retains original source labels,
old/new status, actual first terminal, P-depth and integer clock.
Its canonical hostile is an opposite pair with unequal arrival depths.
Its corrected near miss is deleting a mixed old/new zero packet without
its old-baseline correction. The least-used sidecar is the unique birth
time max(d_s,d_t) of an unordered cross term.

Our live board is **support / owner multiplicity / actual clock / common
feature space / signed Gram energy / unarrived collateral**. The anchor is
the finite native annulus bill; the niche is common-space compression; the
wildcard is a literal native sunflower with overlapping cancellation cores.

## 2. One common observer for every cut

Fix a finite authenticated source family and a finite horizon D.
A source s has old/new tag, original charge A_s, signed charge
v_s=(-1)^(integer clock) A_s, and—if it has arrived by D—actual terminal u_s
and arrival depth d_s. No terminal is required for a source still moving
at D. Write Q(a,b)=a²−ab+b².

Partition arrived sources into bins with the **same tag, terminal and depth**.
For an old bin use the feature column v_s∈Q²; for a new bin use
(v_s,Q(v_s))∈Q³. Give all new sources still moving at D one further
scalar feature Q(v_s). Old sources still moving contribute no arrived
observer through D.

For a rational coefficient vector c on the original labelled sources, let
Φ_D(c) be the sum of these columns in each bin, retaining the moving scalar.
The native application uses c_s=1 throughout. Arbitrary rational c below is
a linear-algebra test domain, not a changed native source measure.

**Common-space theorem (PROVED).** Φ_D determines simultaneously, for
every 0≤d≤D:

1. the old and new signed charge in each actual terminal;
2. the new original source norm already arrived in each terminal;
3. the total new norm still moving;
4. consequently the exact energy increment and the formal residual
   Δ(c)−Σ_arrived new c_s Q(v_s). At the native coefficients c_s=1,
   this is the sum of the two owned cross ledgers.

If r_b is the column rank of bin b and ε is1 when the moving scalar is
nonzero on at least one source, otherwise0, then

    dim P_D = Σ_old bins r_b + Σ_new bins r_b + ε.       (1)

This is the minimum dimension of a linear compression which recovers this
specified **terminal-wise charge-and-norm observer**, for every rational
c. Here the norm bill is the linear observer Σc_sQ(v_s), not the norm
of the scaled atoms. This is a fixed-horizon statement; it is not asserted
minimal for the coarser scalar energy alone, for one fixed native
coefficient vector, or uniformly over all future horizons.

**Proof.** Terminal prefix sums recover the first three observers from bin
sums. Their quadratic expressions give the fourth. Conversely, consecutive
differences of every terminal prefix recover each arrived bin sum;
the horizon-D moving total supplies the last coordinate. Thus the family
of all prefix observers has precisely the same kernel as Φ_D.
The bins have disjoint source supports, so their ranks add. Within each
bin an exact column basis gives coordinates and reconstructs its sum.

For minimality, any linear map L with a decoder recovering those observers
must satisfy ker L⊆ker Φ_D: two coefficient vectors in one L-fibre must have
the same observers. Rank-nullity yields rank L≥rank Φ_D. This argument
allows a nonlinear decoder and assumes no independence of sources.
The checker performs rational basis selection and exact reconstruction.

For arbitrary c, the residual in item4 includes the additional diagonal
term Σ_new(c_s²−c_s)Q(v_s), beyond the cross terms of the scaled charges.
This term vanishes for native coefficients and 0/1 source subsets.
Two identical columns with coefficients (1,−1) and (0,0) have identical
linear features, but different scaled-atom diagonal bills. Thus the
arbitrary-coefficient test domain does not grant preservation of physical
pair ownership after reweighting; its purpose is the rank statement.

For the full annulus, with old charge W_u and new charge H_u at a cut,

    Δ(d)=Σ_u[Q(W_u+H_u)−Q(W_u)],
    C(d)=Δ(d)−M_new(d),
    153N−Δ(d)=R_new(d)−C(d).                           (2)

These are the inherited source-owned identities. Observer compression
preserves them exactly. It does not make their right-hand side nonnegative.

### The diagonal obstruction from the paper

Lemma3 applies to P_D as a space of functions on the original source IDs.
If a k-diagonal tensor on V has the paper's specified decomposition with
all singleton factors in P_D, then |V|≤dim P_D. In particular one cannot
distinguish more original sources with such a tensor than this observer
has rank, even by tuning nonzero rainbow amplitudes.

This is a useful boundary, not a promise that the native cancellation
predicate admits that decomposition. Equal observer columns permit exact
ledger compression but do not authenticate equal source histories.
The original source table and route evaluator are retained.

## 3. The exact sunflower ownership defect

Let S_1,...,S_r be distinct sets of **labelled source atoms**, with
S_i∩S_j=C for i≠j. The charge and actual route data live in an immutable
sidecar keyed by those source IDs. No claim is made that the unlabelled
sets of states along their trajectories form a sunflower.

For any additive cut observer Z_d, including the separate old/new terminal
charge vectors,

    Z_d(union_i S_i) = Σ_i Z_d(S_i) − (r−1)Z_d(C).      (3)

This is simply the incidence identity: a petal atom occurs once and a
core atom r times in the sum on the right. It holds at every actual cut,
not only when the packet has completed.

In one actual terminal, if all packets are all-new and have completed with
Z(S_i)=0, then their distinct union has vector −(r−1)Z(C). For a one-source
core of charge v and r one-source petals of charge −v,

    union energy = (r−1)² Q(v),
    original union source norm = (r+1) Q(v).           (4)

At r=3 these are equal. At r=4 they are9Q(v) and5Q(v). Thus the existence
of arbitrarily many zero-packet sunflower petals cannot itself prove a
source-norm upper bound.

There is a positive ownership rule. If nonnegative rational packet
multiplicities λ_i are used to remove zero packets from a unit-mass
source ledger, the source-capacity guard is
Σ_(i:s∈S_i)λ_i≤1 for every source s. For a nonempty sunflower core this
includes Σ_iλ_i≤1. Choosing λ_i=1/r consumes the core once and only
1/r of each petal; the residual is retained. At a partially completed cut,
even that removal still needs the exposed packet-boundary vector from
the inherited equation(D). Fractional removal is an algebraic bookkeeping
operation, not an authorization to replace the native integer source family.

This isolates the first false implication:
“each selected support cancels” does **not** mean “their union cancels.”
The missing coordinate is source multiplicity, followed by its actual
arrival clock. Merely retaining intersection sizes does not restore either.

## 4. A literal native four-petal example

At base N=4096, use enlarged-scale charges A_(64N). All five sources below
belong to the new annulus and first arrive in the actual frontier terminal955.
Every original charge is (0,1).

| Role | Source | Integer clock | P-depth | Signed charge |
|---|---:|---:|---:|---|
| common core |33435|40|22|(0,1)|
| petal1 |32955|59|34|(0,−1)|
| petal2 |33891|21|10|(0,−1)|
| petal3 |33897|21|10|(0,−1)|
| petal4 |33921|21|10|(0,−1)|

The four source supports {33435,s_i} are a set sunflower. Each is an
authenticated opposite pair after both members arrive. The checker stores
every full native route, independently recomputes the charge from shell
coordinates, and verifies strict first arrival with the inherited627-point
frontier and its decreasing rank.

At cut10, three negative petals have arrived: energy9, arrived norm3,
moving norm2, and cross debt6. At cut34 the whole union has energy9,
arrived norm5, moving norm0, and cross debt4. These are **subpacket**
ledgers. They do not refute the full153N inequality at N4096; other native
sources supply signed cross terms which this subpacket does not contain.

The five-source common observer has dimension **three**, one direction
at each depth10,22,34. The three depth10 petals have the same column.
The exact compiler reconstructs every cut from these three coordinates,
yet keeps all five source IDs and routes for authentication. In particular,
the three depth10 sources cannot be treated as a single source of unit
multiplicity.

The pair {33435,32955} alone has energy1 at cut22 and energy0 at cut34.
This gives a native clock hostile in the same packet, without inventing
alternate trajectories.

The separate [same-parent classification](../sibling_annulus_guards_20261020i/README.md)
shows that two new-annulus sources with the same immediate primitive parent
have unequal original charge norms. Such siblings therefore cannot be
these exact opposite pairs. The example above retains distinct histories
and a later common terminal; it does not replace the immediate-parent
constraint by eventual coalescence.

For a separate complete control, the checker enumerates all672 native
sources at enlarged scale64·7. Its horizon-complete common rank is662;
small observer width is not automatic. At cut8 it reproduces the inherited
counterexample

    Δ=1084, arrived-new norm1054, moving-new norm17,
    cross debt30, and 153·7−Δ=−13.

Again this N7 example lies outside the intended N≥4096 target.

## 5. Why amplitude tuning is not signed-norm preservation

For k=3, the paper's interpolation construction gives, on distinct
a_0,a_1,a_2,

    q(a_0,a_1,a_2) = −(a_2−a_1)/(a_0−a_1).

Thus q(1,2,3)=1 while q(1,3,2)=−1/2. Both are nonzero, as the support
test requires. But their signs and magnitudes differ. The construction
does not require permutation symmetry, positivity, or preservation of a
pre-existing physical Gram coefficient.

The checker evaluates the interpolation decomposition on all64 tuples
for (k,D)=(3,4) and all625 for (4,5). It checks monochromatic value1,
rainbow nonvanishing, vanishing on the other equality types, and each
specified low-degree factor. These are method controls, not a new audit
of the main sunflower theorem.

A native pair coefficient is already fixed:

    beta(v_s,v_t) times
    1[actual terminals agree] times
    1[max(d_s,d_t)≤d].

Replacing it by a tuned nonzero value preserves at most its zero pattern.
It need not preserve Δ, its sign, or the153N budget. Any useful tensor
application must retain the physical coefficients or provide a separate
inequality comparing the tuned tensor with that exact owned ledger.

## 6. The surviving target and executable contract

The positive transfer is exact, finite and source-preserving:
**compile all horizon observers once, with one common space, and keep
overlap multiplicity and exposed boundaries as sidecars.**
This can simplify repeated cut evaluation or identify the exact kernel
of a proposed statistic. It assumes no independence and erases no
moving-source collateral.

The next mathematical obligation remains a bound on the authenticated
native C(d) relative to R_new(d), uniformly over the intended scales.
A sunflower route would need an additional owner-capacity and
clock-boundary inequality. The paper's cardinality theorem alone supplies
neither. A dimension bound on an actually constructed source-diagonal
tester would be usable, but existence of that tester is a separate task.

The source/target/map contract is therefore:

- source: the original finite native rows, scale, old/new marks and routes;
- target: all terminal-wise finite-cut charge/norm observers;
- map: exact bin basis coordinates followed by prefix summation;
- preserved: every stated finite-cut energy and owned cross bill;
- lost unless sidecar retained: individual source identity, within-bin
  history, and future arrival information;
- decisive tests: the native rank-three/five-source flower and its
  energy9/bank5 ownership defect.

Run from the repository root:

    python -B -X utf8 05-knowledge/results/sunflower_owned_support_20261020i/verify.py --write
    python -B -O -X utf8 05-knowledge/results/sunflower_owned_support_20261020i/verify.py

The compiler authenticates source scale, exact types, charge, every native
edge, integer clock and first-arrival boundary before compression.
A route cap is a finite computation boundary: failure at that cap is not
a nonarrival theorem. Complete small-source controls and the five supplied
large-source routes are the only dynamical universes claimed here.
All checks are explicit and remain active under optimization.

Normal and optimized runs reproduce **88,302 checks** and identical saved
JSON. LF SHA256:
`c036b4b149fe5423dc7ab8763ed6bbf542981244c6f8dfd3e27c18d9e77dc575`.
The final coefficient-domain control separates the linear norm bill from
the quadratic diagonal of reweighted atoms; no physical reweighting is
used in the native examples.

Independent peer audit **PASS**: full proof and production API review,
plus a separate optimized replay matching the saved receipt. The review
checked terminal-wise rank minimality, literal source and first-arrival
binding, overlap capacity, partial-clock boundaries, and the distinction
between the linear bill and the quadratic diagonal of reweighted atoms.

