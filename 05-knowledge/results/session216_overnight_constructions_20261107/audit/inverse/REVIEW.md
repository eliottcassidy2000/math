# Independent audit: exact apex-mask inverse

Status: APPROVED mathematical results in the explicit scopes below. Two presentation clarifications are requested; neither changes the two-apex certificate or its proof. Auditor: Session216 independent audit worker. No eventual-spectrum theorem is approved.

## Approved scope

For a fixed k-apex tournament, the positive gap transfer counts Hamilton paths of every nonuniform incidence-mask word exactly. Capping all coordinates at h+1 yields an exact finite target inverse. The reachable state bound is (h+2)^(k! 2 3^k); breadth-first predecessor traces supply shortest witnesses, and exhausted search supplies class-relative absence. This includes k=0 if H(empty)=1 is explicitly adopted. The provided executable target compiler is specialized to k=2, not a generic executable for arbitrary k.

The two-apex executable exhausts 708 states for target7 and 4626 for target21 with output-monotonicity pruning. These are entire-class exclusions, not new ambient exclusions. All 4,998 other positive odd counts through10,000 have supplied literal chain-plus-two-apex witnesses on at most12 vertices. This does not extend the inherited unrestricted 32-bit finite range.

## Mathematical checks

1. A path uniquely determines the apex permutation and the chain vertices assigned to each of its intervening gaps. Within a gap, the chain order is forced. First and last boundary tests are exactly the information retained. Empty interior gaps must additionally satisfy the apex-to-apex edge. States where the final gap has value1 are permanently unreachable; hence the sharper coordinate count k! 2 3^k.
2. min(x,h+1) is a homomorphism for nonnegative capped addition and multiplication. Applying it after each transfer preserves the cap of every exact coordinate. Nonnegative final summation detects equality to h exactly. A capped coordinate may later vanish under a zero coefficient; no rejection based merely on that coordinate is permitted or implemented.
3. Actual tournament Hamilton-path counts are induced-subtournament monotone by deterministic insertion of the added vertex. Restricting an inserted word recovers its parent. Therefore current output above h safely prunes every extension, regardless of hidden coordinates.
4. A shortest successful automaton word cannot repeat a saturated state, since deleting its loop preserves the remaining deterministic transitions and final acceptance. This proves the word bound and BFS correctness; it does not prove target occupancy.
5. At k=0 the empty permutation and one binary gap give two coordinates. The initial readout is1. Appending chain vertices preserves H=1. This is correct under H(empty)=1, whereas the producer's auxiliary hp([],k=0) returns0. That helper is used only at k=2 in the certified packet. Clarify the convention and repair the generic helper if it is advertised at k=0.
6. Explicit budget termination returns INCOMPLETE rather than ABSENT. A crash or external termination supplies no absence certificate. The astronomical theoretical bound is not a practical guarantee or a bit-polynomial algorithm.

## Independent verification

`independent.py` rebuilds the machine using only potentially reachable coordinates, independently reconstructs transitions, and compares against literal permutation enumeration and the producer's transfer. It runs596 controls at k=0,...,4, including every apex orientation through k=3, sampled arbitrary four-apex orientations, empty words, and cap thresholds0,1,3,7,21. The Python auditor uses explicit exceptions, not assertions; normal and optimized receipts are identical.

The independent compressed inverse reproduces all supplied positive/negative target controls. For target7, an additional run WITHOUT output pruning exhausts13,579 saturated states and still proves absence. This specifically tests that coordinate saturation remains sound even when hidden large entries are retained. Budget1 on targets7,21,23 returns INCOMPLETE.

Both witness formats are checked for exact equality, mask-domain validity, maximum length10, unique target coverage, and exact target set {1,3,...,9999} minus {7,21}. `count_all.cpp` uses a separately written pull-style subset recurrence and independently constructs the literal tournament. Every witness passes in optimized and debug builds. Overflow is impossible because each12-vertex path count is bounded by12!, well below2^64.

A separate replay copy regenerates witnesses.json byte-for-byte. Census, exit-budget, suffix-digit, and prefix/suffix/duality-digit receipts agree after dropping timing fields. The finite digit claims remain exactly bounded:18 and96 first-step candidates respectively; neither has a depth-two survivor. The exit-budget statements are restricted to chain lengths at most10.

## Provenance and presentation

The packet correctly credits Session215's multi-exit transfer and Session203's saturated-profile-membership. No generic automata/saturation priority is asserted. This audit confirms a specialization and a class-relative finite certificate, not a new unrestricted coverage range.

Requested presentation repairs:
- Remove README's initial instruction to run probe.py immediately before warning that its main exploratory driver is memory-intensive. Start with the certificate commands instead.
- State H(empty)=1 for the k=0 boundary and identify inverse.py as k=2-specific. If the auxiliary hp function is offered generically, return1 for its empty graph.

## Frozen sources and reproduction

`frozen/` contains the exact reviewed producer inputs. `FROZEN_SHA256SUMS` lists all reviewed original file hashes; only files listed there belong to the freeze (Python may create an irrelevant __pycache__). `REPLAY.json`, `INDEPENDENT.json`, `INDEPENDENT_OPTIMIZED.json`, `ALL_WITNESSES.json`, and `ALL_WITNESSES_DEBUG.json` contain receipts.

From this audit directory:

    python independent.py > independent-replay.json
    cmp independent-replay.json INDEPENDENT.json
    python -O independent.py > independent-optimized-replay.json
    cmp independent-optimized-replay.json INDEPENDENT_OPTIMIZED.json
    g++ -O3 -std=c++17 count_all.cpp -o /tmp/session216-count-all
    /tmp/session216-count-all < frozen/witnesses.txt
    cd frozen && sha256sum -c ../FROZEN_SHA256SUMS

No inherited source files were edited. No git publishing was performed. Parent is the sole publisher. Any later producer revision requires a recorded diff or separate hash review before being represented as these exact frozen inputs.
