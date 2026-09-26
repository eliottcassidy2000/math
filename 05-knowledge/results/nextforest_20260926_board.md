# From an exact carry forest to certified first descent

**PROVED scoped results / INDEPENDENTLY AUDITED / FINITE-EXACT certificates /
Collatz OPEN, 2026-09-26.** This continues the
[ordered-carry forest](forest_20260926_board.md). No literature-priority claim.

## 1. The strongest new result

For odd n let U(n) be the odd part of3n+1. Write its binary digits as b_s,
set c_0=1 and c_(s+1)=floor((3b_s+c_s)/2), and define

    C_n(2)=sum_(s>=1)c_s 2^s,       E(n)=8n-C_n(2).

The [boundary lane](nextforest_20260926_boundary.md) proves

    2n+6 <= C_n(2) <= 8n,
    E(n) >= 2[(3n+1)-2^floor(log2(3n+1))].

Equality in the second line holds exactly when n's finite binary word
contains no00. The defect is the nonnegative charge of two transitions
in the three-state carry graph, with the charge weighted by its bit position.

**All-height certified theorem:** every odd n>1 with E(n)<=4096 has a
smaller odd iterate within104 applications of U. There is no upper bound
on n. In fact the certificate covers the larger set whose displayed
dyadic remainder is at most2048. It checks683 positive offsets, the
separate zero-offset family, and10,977 finite exceptions; the remaining
infinitely many inputs follow by an exact affine lift. The maximum
large-input horizon is66; the finite exceptions raise the uniform bound
to104. A hand proof gives four steps for E<=12.

This is a first-descent theorem. A later iterate can have larger defect,
so convergence of every member does not follow just from class membership.
The defect itself is not a Lyapunov function: even a descending return
can create unbounded defect.

## 2. Beyond27: an exact family, its frequency, and the transient

The family

    n_k=(4^k+17)/3, k>=3,

has constant defect36 and recurrence `n_(k+1)=4n_k-17`.

| k | source | first smaller odd iterate: number of U steps |
|---:|---:|---:|
|3|27|37|
|4|91|28|
|5|347|6|
|6|1371|6|
|7|5467|6|
|every k>=5| (4^k+17)/3 |6|

For k>=6 the sixth iterate is `243*2^(2k-10)+5`; direct comparison proves
it is smaller, while each of the first five is larger. The k=5 collision
has an extra division and is checked separately. Thus this recursion
preserves a precise arithmetic defect, but does not preserve27's long
first excursion. It supplies a useful hostile to identifying all
self-similar families with increasing stopping times.

There are exactly `max(0,floor(log_4(3X-17))-2)` such members below X
for X>=27, hence logarithmic growth and density zero. Its dyadic limit
is17/3, since `n_k-17/3=4^k/3`. Extending the arithmetic to odd rationals
with odd denominator sends17/3 to9 in one odd step. These are positive
rational-basin shadows, contrasting with the expanding rational-cycle
shadows used by THM-4507. The integer source varies with k; no limit
argument replaces an orbit of one fixed integer.

The exact core is9. The available binary precision is2k-1. At27 this
precision expires before a contracting core prefix can be lifted. At
larger k that same prefix becomes legal and gives the six-step descent.
The [precision lane](nextforest_20260926_precision.md) retains what happens
at the failed boundary, including the occurrence of233 in91's orbit.

## 3. Inheritance and six live concepts

Anchor: selected first descent for a fixed integer. Niche: dormant GMC
observer/certificate mechanisms and Rule-30 residual-state methods.
Wildcard: source-preserving precision exchange at a valuation collision.

Closest mechanism: the exact forest polynomial and source-labelled affine
word map. Canonical hostile: expanding rational shadows. Corrected near
miss: finite feature count does not itself imply information loss; one
exact rational can encode every digit. Least-used sidecar: the precision
spent at each core step, including its creation at a new representation.

| Live concept | New finding | Consequence / remaining test |
|---|---|---|
| Ordered carry charge | E>=2d; all-height bounded-defect certificate | Use absolute defect and exact source, not normalized proximity |
| Return selection | bounded-defect families descend with a finite certificate | Control the restart when the core budget exceeds available precision |
| Additive digit ranks | six-edge positive current kills all signed fixed weights | Introduce real interactions or selected returns; more linear coordinates cannot help |
| Observer precision | Prouhet twin transcripts collide; rational evaluations are injective | Audit the observer's exact arithmetic domain and precision |
| Residual automata | complete ternary unit tower; exact3^a macro-carry states | Keep numerator and denominator separately; real contraction is insufficient |
| Precision exchange | exact split/swap/collision update | Preserve source and clock; subtraction of exponents is not a decreasing global budget |

Method cards used: *recover niches by more than citation counts*; *one item
left requires a typed residual*; *separate local support from bounded-height
coverage*. No new universal method card is needed: this session adds
distinct evidence and sharper counterindications to those existing cards.

## 4. Recovered mechanisms that actually transferred

**GMC positive circuits -> an all-height additive-rank obstruction.**
[THM-3258, affine Farkas clutch](../../01-canon/theorems/THM-3258-depth-two-affine-farkas-clutch-and-complete-reset-distance-gauge-no-go.md)
motivated the [six-edge certificate](nextforest_20260926_circuit.md).
For Delta_n=V(Tn)-V(n), where V(n)=sum b_s w_s,

    2Delta_3+Delta_4+Delta_5+2Delta_6+Delta_8+Delta_9=0.

All six inequalities force the first four weights to vanish. Opposite
power/Mersenne currents propagate this to every position. Padding proves
the same conclusion if monotonicity is required only above any finite
threshold. Arbitrary signed weights are covered. This is a positive
linear certificate, not a Collatz cycle: the actual integer boundary is
`{6,6,9}->{2,5,14}`. Nonlinear bit interactions need not cancel.

**GMC nonlinear recovery and Hermite jets -> a precise observer distinction.**
[THM-2631](../../01-canon/theorems/THM-2631-homogeneous-wick-channel-linear-decoder-and-private-row-no-go.md),
[THM-2639](../../01-canon/theorems/THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate.md),
and [THM-4443](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md)
suggested auditing linear recovery, nonlinear recovery and degree bounds
separately. The [moments lane](nextforest_20260926_moments.md) constructs
actual Prouhet source pairs sharing any fixed finite transcript of digit
and carry jets plus a specified valuation bank, while their ratio tends
to infinity. Their common moment values vary: exact first-moment fibres
are finite. Conversely one exact rational evaluation D_n(p/q), p/q>0
and p/q!=1, already determines n. Counting scalar coordinates is therefore
not an information-loss argument.

Adaptive Newton identities recover the complete digit support from its
sparsity and sufficiently many moments. The nonlinear support polynomials

    Q_n(X)=product_(b_s=1)(X-s),
    R_n(X)=product_(s>=1)(X-s)^(c_s)

have the exact transport law, for m=U(n), a=v2(3n+1),

    Q_n(X)^3 X R_n(X)=Q_m(X-a) R_n(X+1)^2.

This retains positions, multiplicities and shifts. A fixed log|Q_n(x)|
is still an additive digit functional and falls under the circuit theorem.
An actual new rank would need interactions, a source-dependent observer,
joint carry information, or a selected-return inequality.

**Rule-30 sections -> an exact complexity and reset calculation.**
[THM-3471](../../01-canon/theorems/THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum.md),
[THM-4204](../../01-canon/theorems/THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md),
and [THM-4210](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md)
transfer a resource question, not their evolution rules. The
[automata lane](nextforest_20260926_automata.md) classifies residual states
as `(0,0)` or `(a,b)` with a>=1,0<b<3^a,3 not dividing b. There are
exactly3^A states through level A, identified with the finite ternary grid.
At each level, the selected output-even transitions cycle through all
units modulo3^a. They are not a fixed integer's ordinary even orbit.

A macrostep q->3^a q+b needs exactly3^a synchronous binary carry states.
Its reset-word count at length h is `max(0,2^h-3^a+1)`. After a fixed
integer's finite input bits end, the exact normalized state beta=b/3^a
contracts by at least a factor2/3 at each step, unconditionally. Its
reduced numerator is the original orbit value and need not decrease.
Indeed `beta*|beta|_3=b`: forgetting the denominator growth loses the
target even though beta is itself an exact, injective state encoding.

## 5. What the new certificate does not settle

The smaller core u=oddpart(d) is always below n. Nevertheless, knowing
its convergence by induction does not ensure enough input precision to
shadow its useful prefix. At27 the available budget5 fails at cumulative
core budget6, before the contracting prefix of budget9. Increasing the
defect cutoff leaves a finite exceptional set that may contain the very
integer under consideration. This blocks a circular strong-induction proof.

Normalized defect cannot supply a uniform horizon either. The exact family

    n_(L,k)=(2^(L+2k)+2^(L+1)-3)/3, L,k>=1,

has E=2^(L+2)-4 and L initial growing odd steps, while E/n tends to0
as k grows. A bound depending on absolute defect must grow with L;
relative nearness to the zero-defect family is insufficient.

The strongest next target is an amortized inequality across the exact
precision split/swap/collision events. It must control the new core and
coefficient after a swap and keep indefinitely open forest intervals.
The three-case update is now exact and tested; a well-founded quantity
for the resulting coupled process remains OPEN. A different positive
target is a genuinely nonlinear support-polynomial inequality that
survives the translated carry factors above. Neither is presently proved.

## 6. Incoming work and a boundary on analogies

During this session, incoming
[THM-4510, sum graphs as reflection orbits](../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md)
connected the square-sum chain at15 to the primitive triple(3,4,5): target
squares9,16,25 compose reflections into a finite rotation. This is a real
answer to the earlier8/9 endpoint observation. It does not identify those
additive rotations with the ternary-unit multiplication cycles above.
The operations and target predicates differ, so no Collatz transfer is claimed.
The other incoming commit repaired the oscillation/Gilbreath finite census;
its whole-orbit shell estimate remains open.

## 7. Audit and reproduction

The four main lanes have self-contained proofs, explicit finite universes,
normal/-O replays, and independent proof audits. The precision follow-up
has its own exact triple and fixed-source checks. The
[independent audit script](../../04-computation/experiments/nextforest_20260926_audit.py)
imports no producer code. It computes carries by residue floors and
independently regenerates every offset certificate and finite exception,
reproducing the104 bound and the27-family times. Its
[output](nextforest_20260926_audit.out) pins the generated transcript.

Run the matching `04-computation/experiments/nextforest_20260926_*.py`
files; the [audit record](nextforest_20260926_audit.md) lists the exact
scope and evidence. The session's theorem is bounded-defect first descent,
not Collatz convergence, an optimal104 bound, or a density-one statement.
