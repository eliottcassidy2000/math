# Lonely-runner wall movies: a two-lap bound, a mirror cancellation, and invisible boundary atoms

**Status.** Exact constructions and proofs, with rational and literal-tournament certificates. The main synthesis sharpens S168's four-lap signed-unit bound to **two laps for threshold-return lifts**, and to **one lap for a reflection-compatible unshifted runner movie**. A shifted two-runner example attains two laps. An exact simultaneous-wall example shows how an asymmetric scheduler can create holonomy even while every boundary atom is retained.

The threshold movie is an externally timed, marked system. It is not an intrinsic rule on unmarked tournaments, is not identified with the blue/black complement-line metagraph, and proves no new case of the general lonely-runner conjecture. “Lap” below means one unit-time lap of this specified periodic input, not the least period of the coarser threshold word. No literature-priority claim is made.

## 1. The old interfaces being joined

The repository already isolates the relevant information losses:

- [LRC Carrier Pullback Index, CPI-009/010/017/071](https://github.com/eliottcassidy2000/math/blob/e92fc100d5e174b839e432deb31ccd1c86202f17/00-navigation/LRC-CARRIER-PULLBACK-INDEX.md): observer/source status, threshold colors, phase-swept class paths, lonely masks, and relative-worldline frame data must travel together.
- [LRC Tournament Technique Index, LTT-010/013/056/060](https://github.com/eliottcassidy2000/math/blob/e92fc100d5e174b839e432deb31ccd1c86202f17/00-navigation/LRC-TOURNAMENT-TECHNIQUE-INDEX.md): a Hamilton count is not a gap meter; safe/danger/boundary need distinct states; arbitrary tournaments need not be realizable clock images; observer-relative speeds and time reversal matter.
- [S168 literal binary compiler](https://github.com/eliottcassidy2000/math/blob/e92fc100d5e174b839e432deb31ccd1c86202f17/05-knowledge/results/session168_computational_structures_20261103/turing/NOTE.md): two strong five-vertex seeds give literal one-arc changes inside one complete path-cover-profile fiber.
- [S168 four-lap clock theorem](https://github.com/eliottcassidy2000/math/blob/e92fc100d5e174b839e432deb31ccd1c86202f17/05-knowledge/results/session168_computational_structures_20261103/clocks/NOTE.md): a signed-unit transport has the form `(s,x) -> (s+f(x),x+a)` and therefore order dividing four.

This note supplies a concrete movie retaining the first payload, compiles it through the second interface, and proves that the special geometry removes some of the third interface's possible monodromy.

## 2. An exact threshold movie, including isolated lonely times

Fix nonzero integer observer-relative velocities `v_1,...,v_m`, rational phases `phi_i` modulo one, and a rational threshold `0<delta<1/2`. On the unit time circle set

    d_i(t) = ||phi_i + v_i t||,
    L_> = {t : d_i(t)>delta for every i},
    L_>= = {t : d_i(t)>=delta for every i}.

Here `||x||` is distance to the nearest integer. A zero relative velocity is a constant constraint and can be retained separately; it is excluded from the wall-count formula below. Choose a cut at `t=0` which is not a wall. This condition holds automatically in the unshifted case. The ordinary anchored LRC input has `phi_i=0` and the appropriate threshold `1/(m+1)`.

The walls are exactly the rational times

    t = (k + sigma delta - phi_i)/v_i,    sigma in {-1,+1}.

Keep those in `(0,1)`, merge equal times without discarding their runner labels, and keep all intervening open time cells. Also retain the **cut atom `t=0=1`**, which is not a wall but is omitted by this open-cell presentation. There are exactly `2 sum_i |v_i|` wall **incidences**, although simultaneous walls reduce the number of distinct times; the cut atom adds no wall incidence or transition.

Give every runner the trit

    danger = 0,    equality = 1,    strict-safe = 2.

It is constant on each open time cell and is evaluated exactly at each wall atom and at the cut atom. Then `L_>` is the union of cells whose trit word is all twos, **together with the cut atom if its word is all twos**. The set `L_>=` consists of those cells and all retained atoms whose trits are at least one. Since the cut is not a wall, its strict and closed membership agree. Equivalently, the first and last safe cells can be joined through a safe cut into one circular open cell. Keep the rational endpoints: a class sequence without time labels cannot recover the lengths of its safe windows. In the certificate, `strict_intervals` are open intervals in the cut chart; the separate `cut_atom.strict_safe` field restores the possible missing point. The chart endpoints `0` and `1` denote the same circle point.

**Proposition 1 (exact finite receiver).** The ordered rational cell endpoints, wall-owner sets, trit words, and cut-atom membership determine both lonely sets exactly, including all isolated points of `L_>=`, and determine their common Lebesgue measure by summing the safe-cell lengths.

**Proof.** A continuous `d_i(t)-delta` can change sign only at its listed walls. Its zeros are transverse because `delta` lies strictly between zero and one half. On each complementary cell its sign is constant. Evaluating all signs at one rational midpoint therefore determines the whole cell. Evaluation at the wall atoms and the cut atom decides every omitted point. The atoms have measure zero. ∎

**Cut regression.** For `v=(1)`, `phi=(1/2)`, and `delta=1/3`, the strict set in `[0,1)` is `[0,1/6) union (5/6,1)`. Its circular open component crosses the cut, so replacing it by just `(0,1/6) union (5/6,1)` would incorrectly delete the strict lonely time zero. The certificate records the two open chart intervals and a strict-safe cut atom; its mass is `1/3`.

This is an exact computation for each specified input, not a uniform theorem asserting that a lonely cell or atom must exist.

## 3. A literal constant-profile tournament movie

Encode the trits by the two-bit path

    danger 00  --  equality 01  --  strict-safe 11.       (1)

At a simultaneous wall, first move every affected coordinate from its left-cell code to `01`; this completes the actual wall macrostate. Then move the affected coordinates from `01` to their right-cell codes. The order within each phase is part of the scheduler. Intermediate words are controller microstates, not additional physical times.

Use S168's seeds `B_0=8`, `B_1=10`, in lexicographic pair-bit convention with a one meaning that the smaller vertex points to the larger. They differ only on arc `{0,2}` and have

    P_B(z) = 9z+33z^2+31z^3+10z^4+z^5,
    F_B(r) = f(r) = r^5+6r^3+2r.

Here `P` counts unordered path covers, and `F(r)=sum_j p_j(r)_j` counts labeled path slots, including empty slots. For a binary word `w` of length `K`, take the ordinal sum

    E(w)=B_(w_1) ⊕ ... ⊕ B_(w_K).

The ordered strong components recover the word, so this map is injective on tournament isomorphism classes. A single bit change reverses exactly arc `{5j,5j+2}` for its zero-based coordinate `j`, and

    F_(E(w))(r)=f(r)^K,       H(E(w))=9^K.              (2)

These interface facts are inherited from the linked S168 proof. The verifier independently recomputes both seed profiles in two ways.

**Proposition 2 (threshold compiler).** The threshold movie on `m` constraints has a literal single-arc realization on `10m` vertices, in one complete profile fiber `f^(2m)`. One input lap makes exactly `4 sum_i |v_i|` threshold arc changes. The fixed grouping into consecutive bit pairs decodes the danger/equality/safe trits at every designated macrostate.

**Proof.** Every wall incidence crosses two edges of (1), and every edge changes one bit. Equation (2) gives the actual arcs and common profile. The retained macrostate flags exclude partially processed simultaneous events from the lonely-time predicate. ∎

The encoding is injective on state words, **not on times**: distinct time cells may have the same word. The event cursor, rational times, fixed runner order, and microstep phase are explicit sidecars.

## 4. Closed threshold returns allow at most two transport laps

At Cayley–Dickson level `n`, the signed basis units are

    U_n = {(-1)^s e_x : s in F_2, x in F_2^n}.

Fix a basis index `a_i` for each threshold coordinate. When its first bit in (1) changes `0->1`, apply the signed-unit permutation `L_(a_i)`. When it changes `1->0`, apply the inverse permutation. Changes of the second bit carry identity. Thus the forward/reverse convention is genuinely inverse-edge transport, not repeated left multiplication regardless of crossing direction.

S168's signed-XOR lemma gives, for every composite and inverse,

    g(s,x) = (s+q(x), x+a).                            (3)

**Theorem 3 (two-lap threshold ceiling).** Over any closed threshold path with the fixed coordinate assignments above, its one-lap return has

    g(s,x) = (s+q(x),x),        g^2=identity.           (4)

Thus an initial unit returns after one or two specified base laps; four laps cannot occur in this model.

**Proof.** The first bit of each coordinate returns to its initial value. Hence it changes an even number of times. Each application of `L_(a_i)` or its inverse contributes the same index displacement `a_i` to (3). The total XOR displacement is therefore zero, proving (4). Squaring adds the same Boolean sign twice. ∎

The sign `q(x)` can depend on the initial basis index. The theorem does not silently replace an arbitrary higher-level return by a global sign. Nor does it apply when the transport label depends on the edge position rather than the fixed threshold coordinate, when a crossing uses the same non-involutory map in both directions, or when a separate winding clock is transported. A genuine torus cover driven by integer coordinate windings can still realize the four-lap case of S168.

Append `n+1` bits for the signed-unit state to (1). Equation (2) now gives order

    5(2m+n+1),       F(r)=f(r)^(2m+n+1).             (5)

Each register update is serialized by flipping precisely its changed index/sign bits. Since `L_(a_i)` and its inverse change `wt(a_i)` index bits and at most one sign bit, a complete input lap uses at most

    2 sum_i |v_i| (wt(a_i)+3)                        (6)

literal arc reversals, including the threshold track. Intermediate register bits are again latched microstates. This is a bound on changed arcs, not on the cost of discovering walls or computing signs.

## 5. An unshifted runner movie can be made to retrace exactly

For zero phases,

    ||v_i(1-t)|| = ||v_i t||,

so the entire trit word is invariant under `t->1-t`. Also `t=1/2` is never a wall: its distances are zero or one half. Pair the walls at `t` and `1-t`. Pick any serialization in the first half, and use its literal reverse at the partner in the second half. In particular, reverse the coordinate order within simultaneous-wall phases.

**Theorem 4 (mirror cancellation).** The resulting full threshold path is literally `P` followed by `P^{-1}`. Every inverse-edge transport on it, not only a Cayley–Dickson transport, has identity return on every sheet.

**Proof.** Reflection identifies the second half's successive macrostate words with the first half's words in reverse order. By construction it also identifies each serialized microedge with its inverse. Every edge transport therefore cancels, in reverse order, without a commutativity assumption. ∎

The assertion concerns this reflection-compatible lift. A simultaneous geometric event does not by itself choose an order for noncommuting transports. If register microsteps are to retrace too, compile the inverse register update in the reverse bit order of the forward one.

This gives a useful hierarchy for clock design:

    arbitrary S168 signed-unit clock: up to four laps;
    fixed-coordinate closed threshold clock: up to two laps;
    mirror-compatible unshifted threshold clock: one lap.

The distinctions are structural, not numerical experiments.

## 6. A sharp shifted example, and its invisible lonely interval

Take two observer-relative velocities, phases, and threshold

    v=(1,2),       phi=(0,1/4),       delta=1/3.

Its six distinct walls occur at

    1/24, 5/24, 1/3, 13/24, 2/3, 17/24.

The exact lonely sets are

    L_>  = (13/24,2/3),
    L_>= = [13/24,2/3],        measure = 1/8.          (7)

At quaternion level `n=2`, assign the first coordinate `i=e_1` and the second `j=e_2`. The chronological transport word is

    j, j^{-1}, i, j, i^{-1}, j^{-1}.

Its first pair cancels and the remaining commutator acts as multiplication by `-1` on every signed quaternion unit. All eight sheets therefore need exactly two base laps.

This is a whole phase chamber, not an isolated numerical example. Keep `phi=(0,theta)` with `0<theta<1/3`. The six wall times, in order, are

    1/6-theta/2, 1/3-theta/2, 1/3,
    2/3-theta/2, 2/3, 5/6-theta/2.

The chronological transport word is unchanged, while

    L_> = (2/3-theta/2,2/3),      measure = theta/2.   (9)

The open inequalities on `theta` verify the displayed order directly; intersecting the two runners' safe intervals proves (9). Every rational parameter gives the exact finite compiler above. The geometric formula also holds for real parameters. Its measure derivative is `1/2`, although the complete compiled profile is the constant polynomial `f^7`. This is a concrete observation-kernel family, not a claim about the polynomial Jacobian conjecture. At `theta=0`, multiple walls and isolated lonely atoms appear, so a generic chamber's event order must not be extrapolated through that boundary.

For comparison, the same speeds with `phi=(0,0)` give

    L_> = empty,       L_>= = {1/3,2/3}.              (8)

Their atom-respecting reflected transport word reduces to identity. Both systems have the same `m=2,n=2` tournament register size: **35 vertices, profile `f^7`, and Hamilton count `9^7=4,782,969` throughout**. The checker constructs every actual adjacency matrix and verifies all literal reversals in one lap: 18 for (8), 19 for (7). The second lifted lap in (7) returns its unit as promised by Theorem 3; a one-lap state comparison alone does not claim closure.

Thus the complete profile cannot distinguish a strict lonely interval from an isolated lonely boundary time, cannot recover the safe measure of the retained timed trajectory, and cannot detect its one-versus-two-lap unit return. The marked state and timed path do retain all three pieces of information. The construction is not evidence that the profile of a natural runner tournament has these properties; it is an explicit compiler into a profile fiber.

## 7. Exact endpoints do not remove the scheduling ambiguity

Keep the unshifted input `v=(1,4)`, `delta=1/3`, with the same quaternion assignments. There are simultaneous entries at `t=1/3` and simultaneous exits at `t=2/3`. Both of the following schedules pass through the exact `01,01` wall macrostates:

1. Process coordinate order `1,2` in the first half and `2,1` in the second. The return is identity.
2. Process coordinate order `1,2` in both halves. After canceling the other second-coordinate excursions, the transport word is `i,j,i^{-1},j^{-1}`, whose return is `-identity`.

The strict lonely measure is `1/6` for both; their geometric trit macrostates and time labels agree. Only the noncommuting resolution of simultaneous events changed. This is a scheduler-created square loop, not an intrinsic monodromy of the unshifted runner geometry.

More generally, exchanging two adjacent different-coordinate wall updates changes the lifted word by a conjugate of their commutator. A proposed holonomy statistic must therefore either specify the wall resolution, prove invariance under all permitted exchanges, or work with a genuinely commuting cover. Retaining endpoint owners is necessary for the lonely predicate but is not sufficient to choose a nonabelian transport through a multiple wall.

## 8. Reproduction, certificates, and next uses

Run with the Python standard library:

    python3 check.py > certificate.json
    python3 -O check.py > certificate.optimized.json
    cmp certificate.json certificate.optimized.json
    python3 independent.py > independent.json
    python3 -O independent.py > independent.optimized.json
    cmp independent.json independent.optimized.json

All validation checks use explicit exceptions, so optimization cannot delete them. The main verifier:

- recomputes both seed profiles by subset/path-partition counting and permutation/cut counting;
- compares high-digit Cayley–Dickson recursion with the independent low-digit transducer on every basis pair through level seven;
- checks all `168` unshifted triples from speeds `1..8` at thresholds `1/4,1/3,2/5` for literal reflected paths and identity all-sheet transport;
- checks `400` deterministically sampled signed-speed/rational-phase triples, with varying levels and assigned basis indices, for zero XOR return and the two-lap ceiling;
- compares the wall-cell safe sets with an independently built intersection of rational closed safe intervals, including wall and cut-atom membership, and includes the strict-safe-cut regression above;
- constructs every 35-vertex tournament in the two main examples and verifies each changed literal arc;
- checks both endpoint-retaining resolutions of the simultaneous-wall hostile.

The separate `independent.py` imports no main-verifier code. It uses the exact `1/24` grid containing every wall of the two main examples, and an explicit quaternion multiplication table, to recover their lonely measure, isolated atoms, and full return permutation. It also checks the strict-safe-cut regression, reconstructing strict membership from open chart intervals plus the cut atom at all `1/24` grid points, including zero, one, and both equality endpoints.

The all-input bounds are proved above; the finite sample is a control, not evidence of a general LRC solution. Useful next directions are: exact wall-resolution invariants, transport labels compatible with wall commutation, and coupling the timed lonely barcode to a separate genuine winding cover without conflating their return mechanisms. The central retained payload is now explicit: **metric threshold + endpoint atom + rational duration + chosen wall resolution + initial sheet**.
