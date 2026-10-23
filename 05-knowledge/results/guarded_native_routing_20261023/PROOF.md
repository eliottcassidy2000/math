# Moving repair intervals give legal native routes with polynomial congestion

23 October 2026. **PROVED, scoped:** the pinned contiguous cyclic-cell
families below, and authenticated count-residue guards. **CONDITIONAL:**
the finite-state extension requires an explicit connector and native lift.
Arbitrary original electrical guards and all native metagraph moves remain
outside this all-size theorem. The independent computational paths in this
session are same-agent checks, not an independent-subagent audit.

## 1. Inheritance and the exact improvement

The incoming [coupled-move theorem](../native_guard_moves_farey_20261023/moves/COUPLED_PHASE_MOVES.md)
already proves connectivity, sharp minimum phase-support budgets, and a
quadratic sorting route. Its original electrical example also authenticates
a two-cell move inside a three-lift automorphism orbit. The incoming
[guard theorem](../native_guard_moves_farey_20261023/guards/README.md) closes
the actual electrical guard on complete source-automorphism orbits.

The [return theorem](../hamiltonian_observer_restrictions_20261023/trace/TRACE_WITH_CLOCKS.md)
and [source-tail bound](../shared_frame_restrictions_20261023/RETURN_COST.md)
leave a different obligation: convert accepted return transitions into
legal moves and pay their length and congestion. Connectivity or diameter
alone does not supply that estimate. Sending all paths through one sorted
state can concentrate an exponential population there.

We replace that global sorting hub by a moving interval that repairs the
shared phase. Complementary splices encode the discarded endpoint data,
bounding how many routes can use one edge. Local sorting is used only in a
bounded-size repair interval. This gives an all-size polynomial congestion
bound using the ALREADY DECLARED swaps and coupled seam flips.

The general path-comparison mechanism is classical: Levin--Peres, with
Wilmer, [Markov Chains and Mixing Times, second edition, section 13.4,
Theorem 13.20 and Corollary 13.21](https://pages.uoregon.edu/dlevin/MARKOV/mcmt2e.pdf).
Our carrier, encoding, native expansions and explicit constants are proved
below. No novelty claim is made for canonical paths or fractional routing.

## 2. Actual source, pins, guards, weights and moves

Use the incoming regular host with q owners per cell and t=r+1 cells,
where q>=3 and t>=3 are odd. Its internal and outer tournaments are regular
circulants; on an outer arc i->j, a cross arc (i,u)->(j,v) is present iff
v-u belongs to {1,2} modulo q. The source has qt actual vertices and stays
fixed throughout this proof.

An observer concatenates the cyclic-cut path in each cell. Its actual phase
tuple a_0,...,a_r has differences d_i in {0,1}; keep the unique binary lifts
of these differences. The full actual owner at clock iq+u is

    (i, a_i+u mod q),   0<=u<q.

Fix actual phases at the two endpoint cells. This fixes the actual first
and last path owners, and indeed every owner clock in those two cells.
An additional guard may require the INTEGER weight k=sum d_i to equal a
specified residue modulo h. This is an explicit guard assumption; it is
not a substitute for the original electrical potential.

Put m=lcm(q,h), with h=1 in the unfiltered case. Compatibility of the two
residue prescriptions is tested first. If compatible they specify k=b
modulo m. Translate the initial phase to zero in formulas while retaining
its actual label alpha. Define

    X = { binary words d of length r : sum d_i=b (mod m) }.

Assume X nonempty. Every word determines exactly one actual path with the
fixed alpha. We retain all words, including those with different integer
weights in the same residue. Let N=|X|, with uniform weight 1/N initially.
The augmented phase c_j=sum_(i<j)d_i mod m records the guard as well as the
native phase alpha+c_j mod q. It is computable from the actual word; it
does not introduce extra multiplicities or discard an integer weight.

The legal primitive moves are:

1. Swap adjacent differences 01 and 10. This changes one internal cell
   phase and keeps the integer weight exactly.
2. Replace a consecutive 0^m by 1^m or conversely. This changes weight by
   exactly +/-m and changes exactly (q-1)m/q actual internal cell phases.

Both keep the fixed source, all owners, endpoint pins, residue guard and
Hamiltonian junctions. The second move is an ATOMIC coupled block, as in
the incoming move definition; intermediate partial rotations are not
asserted legal. It is not automatically an elementary blue/black metagraph
edge. The exact integer weight before and after each move is retained.

For r>=m there are M=(r-1)+(r-m+1)=2r-m proposal locations. With probability
1/2 do nothing; otherwise choose one of these locations uniformly and apply
its involution if enabled, doing nothing if disabled. Call this chain P.
Every proposal preserves X and is an involution, so P is reversible for the
uniform law and is PSD by laziness. Every nontrivial directed edge has
capacity mu(x)P(x,y)>=1/(2NM). Self-loops still consume one integer attempt.

## 3. A fixed-length connector and the moving splice

Set L=m-1. For ANY two augmented phases u,v there is a length-L binary
connector: let delta=v-u mod m in {0,...,m-1}, and use

    1^delta 0^(L-delta).

Fix this connector once. For two accepted phase paths x,y with the same
endpoints, and 0<=i<=r-L, define z_i to have:

- y's prefix through clock i;
- the fixed connector from y_i to x_(i+L);
- x's suffix from clock i+L onward.

Every z_i is accepted. In particular it retains the augmented terminal
phase, not merely its native residue. The route

    x, z_0, z_1, ..., z_(r-L), y                         (1)

has at most R=r-m+3 transitions, counting possible idle transitions.
Consecutive states change at most L consecutive internal phase positions.
The first and last transitions change at most L-1. If r<=m, use one whole
internal block directly instead; the bounds below deliberately overestimate
this small case.

This is already a linear-length route by bounded repair blocks. It remains
to control multiplicity and expand the blocks into the inherited primitives.

## 4. Complementary splices retain the lost information

For a fixed intermediate stage i and a fixed directed edge z_i->z_(i+1),
form the complementary accepted path w with x's prefix through i, the
fixed connector from x_i to y_(i+L), and y's suffix after i+L. Record also
the two ORIGINAL length-L binary blocks of x and y at that interval.

Given z_i, i, w and these two blocks, reconstruct x from w's prefix, its
recorded block and z_i's suffix. Reconstruct y from z_i's prefix, its
recorded block and w's suffix. Thus the encoding is injective. There are
at most N*4^L possible records for a fixed edge and stage. The initial and
terminal transitions need at most N*2^L records each.

A nontrivial edge changing phase positions J can be a moving-splice edge
only if J is contained in {i+1,...,i+L}. There are at most L possible i.
Including the initial and terminal transitions, its total ordered-pair
multiplicity, counted WITH occurrences, is therefore at most

    (L+2) N 4^L.                                       (2)

This is the decisive replacement for the global hub: the auxiliary record
w remains in the same accepted carrier and absorbs almost all endpoint
information. The bounded repair words carry the remaining information.
The estimate uses actual labelled words, not quotient representatives.

## 5. Expansion and an explicit polynomial congestion bound

Inside one repair block there are m binary differences. Its endpoint
augmented phases determine their weight modulo m. If that residue is
nonzero, the weight is uniquely determined, and any two completions are
joined by adjacent swaps, sorting locally and unsorting. If the residue
is zero, the only completions are 0^m and 1^m, joined by one seam flip.
The maximum number of completions is

    C_m = max(2, binom(m,floor(m/2))).

Every block edge has a primitive expansion of at most

    J_m = 2 floor(m^2/4)

events. Sort ones before zeros and count zero-before-one pairs: this count
strictly decreases during sorting, and unsorting has the reverse direction.
A fixed oriented primitive edge occurs
at most once in an expansion. Every expanded route in (1) has length at
most D=R J_m and phase-work at most D*(q-1)m/q.

Fix a primitive directed edge. It lies in at most m length-m difference
windows. Within any such window its exterior is fixed, and at most C_m^2
ordered pairs of block completions can have an expansion using the edge.
Using (2), at most

    m C_m^2 (m+1) N 4^(m-1)

ordered endpoint pairs/occurrences can use this primitive edge. Give every
pair the uniform demand 1/N^2 and majorize every path length by D. Dividing
the resulting load by the capacity 1/(2NM) gives

    A(m,r) = 2(2r-m)(r-m+3) m(m+1) J_m C_m^2 4^(m-1).   (3)

Consequently

    Var_mu(f) <= A(m,r) E_P(f),
    gap(P) >= 1/A(m,r).                                (4)

For completeness, (4) follows by writing f(x)-f(y) as the sum of its
increments along the route, applying Cauchy--Schwarz with its integer
length, and interchanging the sums over ordered pairs and directed edges.
Loops contribute zero energy; repeated edge occurrences are counted.
This proves the consequence, not merely a diameter or acceptance statistic.

For q=m=3 and the regular-host lengths r>=4,

    A(3,r)=13824 r(2r-3),
    D=4r,     each primitive changes at most two cell phases.   (5)

The constants are conservative; the useful fact is that congestion is
quadratic in r at fixed guard modulus. It is not a depth-uniform constant.
No optimal spectral-order claim is made.

If 2<=r<m, the endpoint residue determines one exact weight. Sorting and
unsorting the entire word uses J_r=2 floor(r^2/4) swaps. The crude bound

    A_small=2(r-1)2^r J_r

follows by routing at most N^2 pairs through any edge and using N<=2^r.
Here r<m is bounded at fixed m. A singleton has no nonconstant mixing
claim, and r=0,1 families are handled as deterministic boundaries.

## 6. Realize any stationary return demand without concentrating it

Let R0 be ANY Markov kernel on this same X with stationary law mu. In
particular it may be a guarded trace of an ambient uniform chain. For
each requested transition x->y, independently choose z uniformly in X and
route by the canonical path x->z followed by the reverse of y->z.

The first-leg endpoint demand, summed over y, is mu(x)/N. The second-leg
demand, summed over x, is mu(y)/N by stationarity. Thus both are uniform
all-pair demands. Each complete route has at most 2D events. Using the
constant-length majorant in section 5, its total length-weighted congestion
is at most 4A(m,r). Every intermediate state is a legal accepted native
observer. No rejected ambient excursion state is executed.

This is a fractional routing theorem, not an assertion that the primitive
chain has exactly the return kernel R0. The endpoint transition law is
retained by the compiled random route; its microstep law and duration differ.

Uniform z can be sampled without listing all words. Backward transfer
counts satisfy

    N_i(c)=N_(i+1)(c)+N_(i+1)(c+1),
    N_r(c)=1_(c=b).

They give exact branch probabilities or an integer rank/unrank decoder in
O(rm) table entries. The integers have O(r) bits. This is where the shared
phase transfer becomes an actual routing tool. It keeps common-phase
labels and all kernel lifts required by the input parametrization.

## 7. Nonuniform source weights and honest clocks

Suppose the actual paths have supplied positive weights a<=w(x)<=b and
pi(x)=w(x)/(N Z), Z=E_mu w. Apply Metropolis acceptance to the symmetric
primitive proposal P. Its conductances obey

    pi(x) P_w(x,y) = min(w(x),w(y))/Z * mu(x)P(x,y)
                  >= (a/Z) mu(x)P(x,y).

Also Var_pi(f)<= (b/Z) Var_mu(f). Therefore gap(P_w)>=a/[b A(m,r)].
This keeps the actual weights; it does not assume b/a is depth-independent.

For any pi-stationary return kernel, still choose the intermediate z
uniformly. The two leg demands are bounded by b/Z times the uniform pair
demands, while native conductance is at least a/Z times its uniform value.
The explicit routing congestion is at most 4(b/a)A(m,r). This is an upper
bound on the original weighted demand, not a replacement of pi by uniform
measure.

Retain an integer ledger with: original ambient excursion length tau;
the sequence of actual owner clocks; integer weights along the route;
native event length ell; and total changed-cell work W. The construction
gives ell<=2D and W<=2D(q-1)m/q. Since tau>=1, the same constants bound
ell/tau and W/tau. It does NOT assert ell=tau. If a literal same-duration
simulation is required, padding by lazy native steps works when tau>=ell;
shorter tau can fail. Even a one-step refresh can request two phase paths
that require more than one legal bounded-support move.

Every nonzero cyclic change in one cell moves all q actual owner clocks,
so the corresponding changed-owner work is exactly qW. This is retained
as a separate cost, rather than identifying an owner change with a phase
change or with one event.

Thus a preserved integer clock is an explicitly charged clock, not a
renamed return-event count. If the consumer requires preserving the exact
integer weight k as an invariant, forbid seam flips and stay in one layer;
the residue theorem does not give permission to erase that stronger mark.

## 8. Multiple pins, kernels, converse and the guard boundary

Additional fixed cell phases divide the word into independently pinned
intervals. If the authenticated guard factors into interval count-residue
conditions, the accepted set is the product of the corresponding X_j.
Choose interval j with probability A_j/sum A_j and run its local chain.
Tensor variance decomposition gives gap at least 1/sum A_j. Shared pins
never move. Correlated guards need their own joint compatibility proof;
the product statement must not be applied to them.

When the initial common phase is unpinned, retain the q actual lifts.
They are connected by the actual whole-source common rotation, with its
cost of r+1 changed cell phases recorded. A refresh of that separate
coordinate plus the word chain gives a product chain; choosing one lift
without retaining its multiplicity would change the source measure.
Additional indexed pins may remove or correlate this freedom and must be
substituted before invoking a product description.

Converse transports the actual path to its reversal in the converse host,
and clock j to qt-1-j. Conjugate each full move and all pins by this map.
Support budgets and event counts are unchanged. A printed residue is not
silently substituted for the transported guard record.

Finally, a constant-size guard automaton alone is insufficient. Accept
only the words 0^r and 1^r, with r a positive multiple of 2q. They satisfy
the same endpoint pins and are recognizable with a fixed small automaton.
Any move connecting the two accepted states changes exactly r(q-1)/q
internal native phases. Hence no uniformly bounded support can connect
these guarded families. The failed implication is small transfer width
implies bounded connector length. The repair is a PROVED connector property
of the complete accepted state, including pending owners and guard memory.

The original electrical guard is not known to have that property on the
unbounded cyclic-cell family. The finite tests in ELECTRICAL.md retain its
actual character and target ownership, and do not promote them to an
all-size factorization. Our new bound closes the routing/congestion step
for the declared count-guard families; authenticated electrical connectors
and a native Collatz identification remain open.
