# Uniform strict scalar descent for unrestricted rooted grafts

Status: PROVED candidate, awaiting independent audit. This is a reduction, not a Chan–Pak completeness proof. No new unrestricted coverage range is claimed.

## Inheritance and scope

The literal graft and its scalar row are inherited from `session225_graft_recovery_20261108/PROOF.md`. Physical four-state representative bounds are inherited from `session228_recursive_constructions_20261108/physical_state/PROOF.md`. The nonsink endpoint inequality is THM-4876, `01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md`, statement 7, independently audited there. Reachable-root strongness is Session226 `strong_grafts/PROOF.md`. The new synthesis below applies to ALL nonempty gadgets, without fixing a menu, requiring an invertible matrix, or a prescribed output state. It strengthens the elementary strict-descent observation by combining the endpoint bound with the physical membership oracle. It should not be confused with a proof that a successful branch exists for every target.

Closest mechanism: the Session228 finite cyclic inverse. Hostile: exact cone and lattice membership do not imply physical membership. Corrected near miss: a decreasing exhaustive inverse can still have no successful branch. Least-used mechanism: SCC removal before the root, coupled to scalar prime targets. Board: graft positivity, endpoint bound, physical state oracle, source SCC normalization, prime targets.

## 1. A uniform two-input descent theorem

For a nonempty rooted tournament (A,a) let x=(D,E,U,d) have the inherited meanings. A proper graft here means that a is a nonsink; it need not reach all vertices. Let B be any nonempty tournament, H=H(B), and Q its ordered two-bin path-cover count, allowing empty bins. Form G by B→a and A−a→B. Let N=H(G), V=U−d, f=U−D−d+E, J=D−E+U. Then

    N = H J + (Q−2H)V
      = Hd + Q(D−E) + (Q−H)f.                 (1)

These are rearrangements of the inherited graft row. Physicality gives
D≥d≥1, d≥E≥0, f≥0, U≥D, Q≥2H, and V≥1. Also D−E≥1 for a nonsink. Hence

    H+Q ≤ N,     H ≤ N/3,
    U ≤ floor(N/H)−1.                          (2)

Indeed the second form in (1) gives N≥H+Q. Subtracting HU in the first gives N−HU=H(D−E)+(Q−2H)V≥H.

The inherited nonsink endpoint bound E≤(D−1)/2 now gives

    N ≥ H(3D+1)/2 + (Q−2H),
    D ≤ floor((2N−2Q+3H)/(3H)) < N.            (3)

To check the strict inequality at the small boundary: when H=1,Q=2 it gives D≤(2N−1)/3<N for N≥3; larger H,Q only strengthen the bound. Thus BOTH input Hamilton counts H and D strictly decrease, even for singular gadgets and nonsimple source graphs.

Every assignment of vertices of B to two ordered bins contributes at least one cover, by Rédei existence on each induced part. Thus Q≥2^{|B|}, and (2) gives

    |B| ≤ floor(log2(N−1)).                    (4)

This is a bound on all possible gadget orders at a fixed scalar target, not merely a chosen finite menu.

The nonsink hypothesis is necessary: for a sink source, x=(d,d,d,d) and the output is N=Hd. Singleton sources can therefore reproduce every gadget's scalar unchanged. Including this boundary would make a supposed universal inverse tautological rather than decreasing.

## 2. Complete finite scalar inverse over all gadgets

Given an odd N≥3:

1. Enumerate every labelled nonempty B up to (4), retaining its pair (H,Q) and one physical representative; discard H+Q>N.
2. Enumerate integers D,d,E,U subject to D,d positive odd, d≤D, 0≤E≤min(d,(D−1)/2), d<U≤floor(N/H)−1, U≡E modulo2, f≥0, and (3).
3. Retain only tuples satisfying (1). Apply the Session228 physical four-state oracle. Its finite representative bound here is

       R = floor(log2(U−d))+2             if d=1,
       R = d+floor(log2(U−d))+2           if d>1.

   Exhaust every rooted labelled tournament of order at most R and test its full state. Return actual grafts from successful representatives.

This algorithm terminates and accepts exactly the counts having a proper graft representation. Necessity is (1)–(4) and the representative theorem. Sufficiency is literal construction from a physically verified source and gadget. No graph-image assumption or rational matrix inverse is used. In particular it covers rank-deficient gadgets, rather than discarding them when an inverse matrix fails.

For an explicitly strong output, filter source representatives by a reaching every source vertex. This remains complete: the representative compression in the inherited proof only removes initial inactive singleton components. Such components do not exist when a reaches all vertices, so the order bound already holds for the original source. The inherited strongness theorem supplies exactly this filter.

The source state oracle is potentially enormous. This is not an efficient algorithm. It is not enough to recurse on scalar D alone: replacing a source by an arbitrary graph with the same D can change E,U,d and destroy the target equation. The retained full state is indispensable. A scalar success certificate must return the actual source graph or an equivalent full-state representative.

## 3. Prime targets have reachable-root predecessors

Let C be the ordered prefix of all SCCs of A strictly before the SCC containing a, and let R be the remaining rooted suffix. Then

    x(A,a)=H(C)x(R,a),
    H(G_B(A,a))=H(C)H(G_B(R,a)).                (5)

The first identity follows by concatenating Hamilton paths; the rooted-start bin cannot contain vertices of C, so the other bin contains all of C in every U-cover. The second follows by linearity of (1), or directly by ordered decomposition of the graft. The source root reaches all of R. If a is a nonsink of A it is still a nonsink of R, since a cannot beat an earlier SCC. The latter graft has count at least3.

If N is prime, (5) forces H(C)=1. Removing C preserves the entire state, and hence gives a reachable-root predecessor and a strong output with the same target. Thus the strong filtered scalar inverse in Section2 is already complete for prime proper-graft targets. This is not a statement that every prime has such a predecessor.

Moreover gcd(H,Q) divides N by (1). Since H≤N/3, any proper-graft representation of prime N must use a gadget with gcd(H,Q)=1. The latter divisibility is inherited in essence from Session227's scalar gcd obstruction; the point here is its combination with the strictly decreasing, unrestricted inverse and reachable-root normalization.

## 4. Exact low-target compiler and checks

`verify.py` implements a small complete reachable-root source compiler under D≤20,U≤31. The source begins with the strong SCC containing its root, then an arbitrary ordered list of strong components. A singleton is allowed as the initial SCC. If a nonempty SCC T with (h,q) is appended after the rooted graph, its state becomes

    (D h, 0, U q, d h).                         (6)

Both bins may independently use a final path in T, proving the U factor q; the other three entries are direct concatenation. Rooted initial SCCs are exhaustively taken from all labelled tournaments through order6. Completeness of that order cap uses the CITED Busch strong floor f(7)=25 and monotonicity, recorded in `camion_busch_gaps_polyhedra_collatz_20261002.md`, Section2. A strong initial SCC with D≤20 cannot have seven or more vertices. Later SCCs have q≤31 and therefore order≤4. Formula (6) makes source-state merging safe for every later suffix. Each append doubles U at least, so the compiler terminates even with arbitrary singleton tails.

The compiler finds 77 distinct states. Testing all gadget pairs through order4 is complete for target≤31 by (4). Its exact reachable-root proper-graft scalar image through31 is

    3,5,9,11,13,15,17,19,23,25,27,29,31.

For any source producing N≤31, (3) implies D≤20 and (2) U≤30, so the slightly larger U≤31 bank contains every required source. All witnesses are checked by literal subset path DP on the actual graft, independently of (1). This is a deliberately small validation of the ALL-GADGET inverse and strict decrease, not a new unrestricted spectrum range; the repository already has much larger attained ranges.

The absence of7,21 in this filtered image is a scoped exact computation. Their known unrestricted impossibility is not proved anew here. The presence of every other odd target through31 does not establish future branch coverage.

The file `inherited.py` is copied unchanged from Session225's portable verifier. `verify.py` supplies the new SCC compiler and literal controls; `CHECK.json` is its deterministic output. Run `python3 -I -B verify.py > replay.json` and `cmp replay.json CHECK.json`, and repeat with `-O`.

## 5. What constructor is still missing

This closes a specific logical gap: unrestricted proper graft inversion CAN be made finite, physically sound, and quantitatively decreasing in BOTH input scalar counts. It does not close the existential gap. A successful induction still needs a full-state-preserving family whose inverse candidates hit a physical state for every admissible N. The fixed-menu obstructions do not prove failure of this unrestricted family, and the present inverse does not prove its surjectivity. A minimal missing target can now be stated without an unspecified feasibility oracle, but no such target beyond the known exceptions is produced here.
