# Independent audit of logarithmic growth and stationary inverses

**PASS.** Session 195, 2026-11-05. This audit concerns structural bounds for the specified rooted grammar, not an unbounded spectrum-coverage theorem.

Frozen source proof SHA-256:

    1ac36e4a683bd9e976f1357b96b2ef4fb80c03dd7c6791ebe6a7d4db0de982dc

The verifier and three inventory snapshots are also pinned in `receipt.json`.

## Proof findings

The invariant R>=2d holds for every reachable coordinate state other than the initial singleton state. The cases for S and a non-singleton gadget use the reset d'=h and d'=dG*h respectively; the latter also needs B>=2dG, correctly supplied by the two distinct empty-slot covers of the nonempty root-deleted gadget. K preserves the invariant, D strengthens it, and the only initial K orbit is a fixed coordinate state.

Every entry of the two-operation table is correct. In particular, after a non-singleton gadget the new deletion count must be compared to its new R; this follows from the invariant and yields the claimed N-then-K lower bound. The initial-state qualification is essential and is correctly handled by deleting leading K loops. Combining the first-step bound with the pair bounds yields exactly R_j>=2^ceil(j/2).

The added non-D contraction is valid: deleting a terminal root gives e<=d, so at a noninitial reachable state e<=R/2. Both singleton-root and stationary-root births satisfy h'>=h+R/2. For a gadget, L>=g>=s implies L-s/2>=1/2, giving the same inequality. Consequently 2h'>=2h+R>=3h. The initial singleton exception is necessary: S and some transitive gadget branches preserve H=1 there. The verifier correctly restricts its added abstract-cone test to e<=d.

The logarithmic witness proof handles terminal overshoot correctly. It bounds the **preterminal** R by the final H, not the final R. A last gadget is allowed to produce R>H. Leading K loops can be removed without changing the entire state; terminal D operations can be removed without changing H. Applying the exponential R bound to the t-1 preterminal steps gives t<=2 floor(log2 H)+1. Multiplying by the maximum number of vertices added per step and including the initial vertex proves the order bound. These are existence bounds for a normalized witness, not bounds on every redundant presentation.

The K/D inverse normal form is valid for every physical stationary seed because d is odd. Parity identifies the final operation, and each inverse strictly decreases R whenever R>d. At R=d, physicality forces (h,e,R,d)=(d,d,d,d), and K really is a fixed state. Decoding the R word is not by itself a full-state membership test: the source correctly requires replay on all four seed coordinates. The singleton f(q) recurrence and the even-q bound q<=2H-2 follow directly from the physical operations.

The root-moving inverse formulas recover precisely the coordinates claimed. The old deletion count is genuinely forgotten, and an s=0 gadget also leaves the old endpoint coordinate unresolved. No arbitrary feasible predecessor is asserted reachable.

## Correction and hostile controls

One precision correction was requested and applied. Two omitted small values do not establish noncofiniteness; the source now says the finite criterion does not establish eventual coverage. No mathematical formula was corrected.

The audit explicitly covers:

- arbitrarily repeatable initial K loops at the singleton;
- long terminal D tails with H=1 and unbounded R;
- singleton and non-singleton first steps preserving H=1;
- final non-D outputs with R>H;
- physical stationary seeds with R=d, including d>1;
- nonempty words whose leading K letters are erased only in the R=d case;
- targets with the correct R but wrong remaining coordinates, rejected by replay;
- gadgets with s=0 and the singleton correction.

Thus the proof's qualifications are necessary, and the tested hostile cases remain inside their correctly stated boundaries.

## Independent checks

The independent verifier uses the actual previously audited bank of 67 physical matrices. It tests every reduced word through length three over the whole bank plus K,D, and longer words through length seven over a smaller alphabet containing the slow singleton branches, stationary moves, a transitive edge and both six-vertex pieces. It checks the invariant, every-two-step doubling, strict non-D contraction, normalized length and order bounds, and every applicable inverse formula.

The stationary normal form is separately checked on 273 different physical four-coordinate seeds from the complete order-at-most-five and order-six inventories, and on every K/D word through length ten. Seven seeds have R=d. Each decoded word is compared to the original canonical word, replayed on the full state, and tested against a target with altered H. The singleton response is independently checked on 32,768 binary inputs by the explicit pair expansion, not by importing the source recurrence.

Exact counts and snapshot hashes are in receipt.json. Normal and optimized Python runs agree byte-for-byte. The source's own abstract-cone verifier is also replayed in both modes from a copied, hash-matched script, with explicit input paths.

## Reproduction

    python3 verify.py --source PATH_TO_LOGARITHMIC_PACKET --extension PATH_TO_GADGET_PACKET
    python3 -O verify.py --source PATH_TO_LOGARITHMIC_PACKET --extension PATH_TO_GADGET_PACKET

The checked source proof, verifier and inventory hashes are recorded separately. The checker rejects changes during a run. No statement here proves that every odd target other than 7 and 21 belongs to the grammar.
