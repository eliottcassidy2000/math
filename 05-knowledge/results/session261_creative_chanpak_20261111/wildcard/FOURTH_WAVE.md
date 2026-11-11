# Fourth wave: reversal, two-ended contexts, and finite-control certificates

**PROVED reductions and physical operations; FINITE-EXACT probes; OPEN digit coverage.** This wave follows the all-size single-letter obstruction, rather than increasing that scan's length.

## 1. Random recovery and changed question

Seed261108 sampled10 further files from the same11,638-file historical pool; its choices are recorded in `RANDOM_WAVES.json`. Read:

* `05-knowledge/results/duality_20261019_multiplication_addition_golden_circle_gilbreath.md`: the useful mechanism is changing charts while preserving an actual orbit, not numerically matching an affine formula.
* `05-knowledge/results/presentation_counts_single_view_classes_rotation_inequality_20261020.md`: rotation/converse presentations retain the object or a count, while a preferred monotone H-descent can fail.
* `05-knowledge/results/session29_framed_profiles_20261022/session29-collatz/phase2/DENSITY_ONE_DESCENT.md`: endpoint success does not preserve the family needed for another induction step.

No historical claim from these files is imported as a new proved dependency. The new question is whether a count-preserving presentation change can alternate the interface used by successive binary digits. Appending both digits in one fixed chart is now known impossible for single-letter two-apex digits.

## 2. A literal count-preserving involution

Keep a_1->a_2. Reverse every arc of the tournament, then relabel a_1 and a_2 by swapping them and reverse the spine order. This again has the same core and forward spine.

On bitmask letters the operation is

    phi(0)=3, phi(1)=1, phi(2)=2, phi(3)=0.

Thus define w*=phi(reverse(w)). Reversing a Hamilton path is a bijection between paths of a tournament and its converse, so

    H(T_A(w*))=H(T_A(w)),       (w*)*=w.

This is an unconditional physical return, available at every word. It is not a same-q marked-seed theorem and does not assert preservation of additional marks. It changes append operations into prepend operations: (wb)*=phi(b) w*.

## 3. Why the old vector is insufficient, and the exact replacement

The positive vector v(w)=M_w v_0 retains all future right extensions, but need not determine a left extension M_w M_b v_0. A count-preserving reversal therefore cannot simply be installed as an arbitrary map on v(w).

Retain instead the context pair

    X=M_w,       Y=M_(w*).

This has at most2D² nonnegative integer coordinates; D=36 for k=2. The identity word starts at (I,I). The exact updates are

    append b:  (X,Y) -> (M_b X, Y M_(phi(b))),
    prepend b: (X,Y) -> (X M_b, M_(phi(b)) Y),
    converse:  (X,Y) -> (Y,X).

All three are nonnegative integer linear operations on the flattened pair. The output is H=L X v_0. The second component supplies the needed sidecar for reversal; ordinary scalar H or the right-response vector alone does not.

The carrier is deliberately not minimal. It is finite-dimensional, explicit, and preserves every two-sided context. A smaller representation may be obtained by taking the reachable/observable quotient, but a signed quotient would lose the simple saturation argument unless positivity is restored.

Since these are positive linear operations, capping every coordinate at N+1 again yields an exact finite automaton for target N over arbitrary sequences of append, prepend, and converse operations. This does not enlarge the ultimate two-apex tournament class; it enlarges the available physical induction steps and makes their legal domains explicit.

## 4. Finite-control binary induction criterion

Let Q be a finite set of control states. Each q has two specified physical operation blocks F_(q,+), F_(q,-) on the context pair and specified destination states tau(q,+),tau(q,-). The operation blocks may include append, prepend, and converse moves. Augment the carrier by a constant1 and denote its output row by ell.

Initialize the row space R_q with the two digit-defect rows

    ell F_(q,+) - 2ell - e_const,
    ell F_(q,-) - 2ell + e_const.

Repeatedly add to R_q every row of R_(tau(q,e)) pulled back through F_(q,e). This finite simultaneous closure stabilizes after at most |Q|*(2D²+1) total dimension increases.

**Theorem.** A literal seed context in control state q annihilates R_q iff every binary digit sequence, executed according to this finite controller, obeys h->2h±1 at every step. Proof: the closure is precisely the span of all future digit defects pulled back to the initial state.

This is a genuine extension of the stationary-cylinder criterion. A return may change both the physical presentation and the control state; it need not return to one common affine line or one marked port. An explicit seed bank for odd counts in [B,2B-1], each satisfying its controller's row kernel, again gives every odd count>=B by the inverse-digit induction in `PROOF.md`.

No such complete controller is supplied. The criterion reduces a guessed finite-control grammar to exact linear algebra plus literal seed reachability, not to a scalar numerical coincidence.

## 5. Connection to the other lane's marked reselection

The parent reports a physical five-vertex tournament with row masks(2,4,9,19,7) and a mark-set switch {0,3}<->{1,4} preserving its full marked-cover tensor. That is a promising different presentation change. It is an inherited lane input here, not independently audited or assumed to preserve this spine carrier.

To combine it with the finite-control criterion, retain the marked-cover tensor used by the actual digit constructors, express the reselection as its exact coordinate map, and use distinct control states for the two port menus. Equality of scalar q alone would not justify the map. If full tensor equality is indeed the constructors' sufficient interface, the reselection is a zero-cost control change and can occur inside a digit block. The remaining obligation is a closed finite menu of physical transitions with the two required count slopes.

## 6. Focused controls and status

`probe_two_ended.py` examines all two-apex seed words of length<=6, with the eight possible single-letter append/prepend operations. It tests immediate plus/minus pairs and then four iterations of the same chosen pair. The receipt is `two_ended_controls.jsonl`. This is a finite control, not a proof of absence and not a search over state-dependent controllers.

`probe_block_returns.py` separately checks append-only digit blocks of lengths1,2,3 on seeds of length<=6. It finds234 immediate candidate pairs in total, but none survives four iterations. This identifies a real family-return failure: satisfying the first digit equations is common enough, but the interface is not restored. The all-size theorem in `EXTENSIONS.md` applies only to single-letter append digits.

The new positive result in this wave is the literal converse return and its complete positive two-sided carrier. The new reduction is the finite-control rowspace certificate. Coverage and a successful controller remain OPEN.

## 7. Transfer ledger

Source: converse/presentation changes and chart returns. Target: two-apex binary count induction. Map: w->phi(reverse(w)), together with the paired context matrices. Preserved predicate: H exactly and every future two-sided query. Destroyed information: the word after its matrix pair is retained. Sidecar: the converse context Y and a literal seed word. Cheapest decisive test: compare H(w),H(w*) and each matrix update with direct word construction.

Source: framed finite descent with a lost return register. Target: finite-control digit grammar. Map: attach a separate row-kernel domain to each control state and pull destination constraints back through each digit. Preserved predicate: every future transition is valid. Destroyed information: presentation history after the carrier and controller are known. Sidecar: physical reachability and the actual port-changing map. Cheapest decisive test: a two-step failure of an otherwise valid immediate ± pair; the supplied probes show this is a serious obstruction, not bookkeeping decoration.
