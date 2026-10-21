# Exact clock-aware cancellation beyond pairs and three-type circuits

2026-10-21. **PROVED:** an exact finite-inventory cancellation compiler,
its independently checkable optimality certificate, and a temporal
non-nesting obstruction. **FINITE-EXACT:** the declared native source
populations below. Uniform long-horizon `R'` control and the native
`153N` target remain **OPEN**.

## 1. Inheritance and the new object

The source is the original primitive integer, with first arrival in the
fixed627-point frontier and its actual integer **P braid clock**. This is
not the fully odd Collatz clock. We inherit the charge and route definitions
from [native annulus](../session18_native_annulus/README.md), chronological
antipodal matching from [owned pairs](../owned_pair_blocks_20261020i/README.md),
and the bank of six opposite pairs and forty three-type circuits from
[circuit completion](../eisenstein_circuit_completion_20261018g.md).
The current [cross-scale credits](../cross_scale_frontier_credits_20261021/README.md)
pay the sibling bank but leave `R'=R-C`, with negative `C` removed only once.
The [recovery overview](../session21_recovery_20261021/README.md) explicitly
distinguishes recovered evidence from reconstructed historical prose;
no missing historical receipt is a dependency here.

The canonical hostile is actual `N=7, d=8`: `I=30 > M=17`. Repacking
sources cannot repair that physical deficit. The corrected near miss is
promoting a larger cancellation count into a smaller physical interaction.
The least-used sidecar in this packet is a four-node divergence integer,
which is lost by ordinary cycle-flow compression.

Concept board: original source and clock / twelve charge types / integer
flow defect / actual arrival cut / mature source ownership / unchanged
quadratic remainder. We optimize cancellation **at a supplied cut**, not
the arrival time, terminal, clock, or physical energy.

## 2. A four-point difference model with an indispensable defect

Write `Q(a,b)=a²-ab+b²` and let

    p0=(0,0), p1=(0,1), p2=(-1,-1), p3=(1,0).

The twelve directed differences `p_j-p_i`, `i!=j`, are exactly the six
norm-one and six norm-three signed native charges, with no repetition.
An inventory `c_ij` records the number of **already arrived new sources**
of that type at one actual terminal. A selected inventory satisfies
`0<=x_ij<=c_ij`, integrally. Put

    b_i = incoming selected count - outgoing selected count.

The selected charge sum is `sum_i b_i p_i`, and `sum_i b_i=0`. Since
`p1+p2+p3=0` and the only linear dependence of these three vectors is this
one, the selected sources sum to zero exactly when

    b = t*(-3,1,1,1),  t an integer.                    (1)

Their owned diagonal norm is `sum Q(p_j-p_i)*x_ij`.
This gives an iff, not merely a sufficient encoding. Ordinary circulation
retains only `t=0` and misses the three outward unit spokes: their total
charge is zero, norm3, and divergence `(-3,1,1,1)`.

There is a second exact boundary. The four charges

    (-1,-1), (1,-1), (1,2), (-1,0)                     (2)

sum to zero and have total source norm8, but no proper nonempty submultiset
sums to zero. Every type occurs once. Thus neither an opposite pair nor
any primitive three-type circuit can be extracted. All46 inherited short
circuits leave (2) unchanged; the new compiler cancels it completely.
This is a real extension of that extraction language, without contradicting
its bounded-residual theorem.

The map retains type capacities, the integer `t`, and then original source
IDs. Without those capacities or IDs, the flow would not be a native
source certificate. Mod15 and mod35 are retained as compatible projections
of each actual source; they are never treated as independent populations.

## 3. Complete optimization and a small proof certificate

For each nonempty subset `S` of `{1,2,3}`, let `in(S)` and `out(S)` be the
total capacities entering and leaving it. The exact feasible integer
interval for (1) is

    max_S ceil(-out(S)/|S|) <= t <= min_S floor(in(S)/|S|).    (3)

Necessity follows by summing divergence on `S`. Sufficiency is the
source/sink flow cut criterion: subsets containing0 are complements of
the seven displayed subsets, so these are all required cuts. With integral
capacities and prescribed integer divergence, augmenting paths produce
an integral feasible flow. Zero is always feasible.

For fixed `t`, maximize the selected source norm. The producer first finds
a feasible flow and then augments positive-cost residual cycles. Costs are
the original norms1 or3, with their negatives on reversal. Each augmentation
strictly increases an integer objective bounded by the inventory norm, so
this finite algorithm terminates. It is not advertised as a new complexity
bound or as the fastest flow algorithm.

A fixed-`t` receipt consists of the selected integral counts and four
integer potentials `pi_i` satisfying, on **every** residual arc `i->j`,

    residual_gain <= pi_j-pi_i.                         (4)

Any competing flow with the same divergence differs by residual cycles;
summing (4) on them bounds its gain by zero. The consumer checks capacities,
divergence, exact zero charge, norm, and every residual inequality. It
does not have to trust the producer's cycle choices.

Let `f(t)` be the fractional optimum for real feasible `t`. Mixing two
feasible flows proves concavity of `f`. At integer `t`, a linear objective
has an integral network-flow optimum: if a vertex of the feasible polytope
had fractional edges, integer divergence would make their undirected
support contain a cycle, along which a small two-sided perturbation
contradicts extremality. Thus the integer optima are the samples of this
concave function. Consequently a fixed-`t` optimum plus certified optima
at its feasible immediate neighbours proves the global integer maximum.
Endpoint cases use one neighbour; a singleton interval uses none.

The producer uses discrete concavity to locate this maximum. The verifier
checks (3), (4), the complete neighbour set and the neighbour inequalities.
The resulting value `Z*(c)` is the **maximum owned source norm in any zero
submultiset of this supplied finite inventory**, and `A-Z*` is its minimum
residual source norm. This includes all packet sizes and repeated types.

## 4. Arrival times prevent an irreversible optimal packet bank

Suppose the first three source charges are

    U0=(0,1), U1=(-1,-1), -U0=(0,-1).

The unique maximum zero subset is `{U0,-U0}`, of norm2. When `U2=(1,0)`
arrives, the unique maximum becomes `{U0,U1,U2}`, of norm3. The first chosen
subset is not contained in the second. This four-source example proves
that all-alphabet maximum cancellation cannot always be realized by one
irrevocably increasing packet bank. The stronger simultaneous optimality
of inherited rank-pair matching remains valid **inside each antipodal
class**, and is not a theorem for this enlarged optimization.

Our certificates are therefore cut-indexed. At each cut every selected
source has already arrived, every ID is used at most once, and all selected
sources have the same terminal. Reassignment at another cut is allowed;
credits from distinct cuts are never added. A completed zero packet has
zero cross interaction with every external vector at that cut. Before its
last member arrives, its partial charge remains fully charged.

There is also an integrality hostile: the inventory

    (0,1), (-1,-1), (1,-1), (1,2)

has norm8, vector `(1,1)` of norm1, and no nonempty zero submultiset. One
cannot silently replace integer packet selection by a fractional cone or
assume every zero-free small inventory lies in a strict half-plane.

## 5. Native first-arrival controls and larger mature blocks

The checker reconstructs all primitive sources `1<=s<256N` at `N=7,4096`,
retaining each source's mod15/mod35, terminal, integer clock, depth, unsigned
charge and norm. It uses the inherited route cap10000 and authenticated
frontier; an unresolved source stops the run rather than being discarded.
Independent palette coordinates check every charge. The route DAG is
replayed, and source-row hashes and all257 physical horizon rows are saved.
All sources in these finite populations arrive by depths22 and149,
respectively, so their ledger is constant after those authenticated cuts.

The exact optimizer is applied at all cuts0..256 to terminals1,87,711 at
N4096; N7 terminals1,3,7 are small boundary controls. These are selected
terminal controls, **not** a claim to optimize every terminal at N4096.

| terminal | first cut maximizing improvement over antipodal matching | arrived norm | antipodal cancellation | exact cancellation | residual norm |
|---:|---:|---:|---:|---:|---:|
|1|11|210|142|199|11|
|87|35|207|164|203|4|
|711|47|6378|5930|6317|61|

At terminal711 this is one mature zero block with3747 distinct source IDs
and norm6317. The deterministic lexical46-circuit greedy comparison
cancels6277 there; its40-unit gap is **not** a claim that every possible
short-circuit ordering has that gap. The unconditional language separation
is the irreducible quartet (2).

At each displayed best cut, the receipt contains the full selected source
ID list, all terminal source rows, the optimum flow and neighbour/potential
certificates. Every selected source is already present at that cut. The
source compiler preserves the original clock sign rather than recomputing
it from a residue or a newly chosen path.

### Inherited all-scale lift of the mature711 block

The selected3747 new sources give one exact infinite family of **lifted
subsets**, using the inherited inverse law, not another full-annulus census.
For each source let `a=2` when `s=1 mod8`, and `a=1` otherwise. For `k>=0`
put `m=64^(3k)`, `N_k=4096m`, and

    s_k = m*s + (2^a+3)*(m-1)/9.                        (7)

The [native phase-lift identity](../session18_native_phase_lifts/README.md)
gives the same literal P parent and an edge clock increased by18k. The
first-arrival terminal711 and depth are consequently unchanged provided
the lifted source is outsideV. At `k=0` this is authenticated: all selected
sources are new and their terminal711 is old, so none is initially inV.
The least selected source is33183; for `k>=1`, `s_k>=262144*33183>63665=maxV`.
Thus the guard holds at every lift.

The additive term lies in `[0,m)`, so every even-ended integer annular band
is preserved for odd sources. Moreover `s_k=s mod3`: `(m-1)/9` is divisible
by3 for this three-scale lift. The six-band charge formula and even clock
increment therefore preserve every signed charge. Different source IDs
remain distinct: equal branches give affine injectivity, and the two
branch offsets differ by less thanm, so cannot cancel a nonzero multiple
ofm. The lifted block is mature by depth47, has zero vector and owned
norm6317 for everyk. It remains zero at all later finite cuts.

For a requested integer modulusq, compute the numerator
`9*s*m+(2^a+3)*(m-1)` modulo `9q` and only then divide by9. This preserves
the exact mod15/mod35 registers without expanding a huge source. Bounded
supplementary checks replayed all3747 selected sources at `k=1,2,7`
(11241 lifts), including the actual parent/clock, band, charge and full-
denominator residues modulo15,35,105. The all-k conclusion is the proof
above, not extrapolation from those probes. No optimality in the *complete*
larger inventory, proportional budget growth, or total-energy gain follows.

## 6. Exactly what this does, and does not, bound

At one terminal and one cut let `W,H` be the actual old/new signed vectors,
`A` the arrived-new source norm, `Z` a certified mature zero-subset norm,
and `T=A-Z` the residual source norm. The residual source vector is still
`H`. Define

    Gamma = Q(H)+beta(W,H)-T.

Then the source-owned interaction is exactly

    I_terminal = Gamma-Z.                              (5)

Using `Z*` maximizes mature cancellation and minimizes the retained source
norm, but increases `Gamma` by exactly the same difference. Thus no
physical energy saving has been manufactured. The packet makes the
strongest zero-subset compression at that cut; the coherent residual still
owes its quadratic bound.

After summing terminals, the exact sufficient/necessary guard is

    sum Gamma <= M + sum Z.                            (6)

If external/sibling interaction owners have already been extracted,
`R'=I-J'` is unchanged as well. The new source compression is not an extra
negative pair credit to be added to `J'`. All original cross terms and the
transferred `R'` remain owned. The finite N7 deficit is unchanged: at depth8
`A=1054, M=17, I=30, DeltaE=1084`; its full finite population later has
excess26. At N4096 the largest observed full-ledger excess is `-201674` at
depth149. These inherited finite facts are freshly replayed here, not
promoted to a uniform theorem.

The exact remaining structural question is whether the coherent residual
in (6), or an independently certified stopped-prefix tail, can be bounded
uniformly at long horizons. Neither a finite optimal flow nor eventual
arrival of the tested sources answers that question.

## 7. Reproduction and hostile controls

    python -B 05-knowledge/results/long_horizon_pairs_20261021b/check.py --write
    python -B -O 05-knowledge/results/long_horizon_pairs_20261021b/check.py

The standard-library verifier checks all4096 binary twelve-type inventories
against independent exhaustive subset enumeration, plus all5940 inventories
with three nonzero types and multiplicities1..3. It retains exact-int guards,
the circulation/star hostile, irreducible quartet, non-nesting example,
zero-free integer hostile, forged objective/potential and missing-neighbour
rejections. Checks are explicit exceptions and remain active under `-O`.

Native proof inputs are pinned by the frontier LF hash; the finite source
population and every selected row are reconstructed. The source-ownership
helper deliberately assumes its rows have already been authenticated by
that reconstruction; it is not an independent route-discovery oracle.

Normal producer:740191 checks. Saved receipt LF SHA256:
`a365ec43470e2c34907eb5239c1df7df1d1a9338673cb5013a632dde92363a68`.
Independent normal/optimized producer comparison matches; optimized immutable
receipt verification reports740192 checks, including its final equality gate.
Independent proof/API audit and optimized saved comparison: **PASS**.
The reviewer checked the exact divergence kernel, seven cut inequalities,
residual-potential optimality, integer-flow/concavity neighbour certificate,
mature source ownership and both language hostiles. The inherited all-scale
corollary also passed an independent proof check. No correction remains;
the three source/receipt files are frozen at the hash above.
