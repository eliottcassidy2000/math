# Independent review: adaptive active-depth inverse bounds

**PASS, no correction.** The weighted local lower bound, adaptive product cutoff, actual internal-order bound and weighted triangle budget follow with the stated hypotheses. Internal passive padding is retained, and no finite enumeration of unrestricted context templates is inferred.

## Frozen source binding

The reviewed packet is under `session197-odd-cycle-path-polynomials/`:

- `ACTIVE_DEPTH_INVERSE_BOUNDS.md`: `3c5ac22cb8a730b9fe01926154012fd98475f84213248ced717d7976c46c32f8`.
- `verify_active_depth_inverse_bounds.py`: `efa065fb084794c001d0b02d53f6e0ad56a036804664d841f50eed3b200a40f3`.
- Normal and optimized receipt: `e01dc56b02432f8b1c8ce56b406bfcb3145cbfb46cac97f5a69fae86b6757dd4`.

The exact named universal-context/profile proof and code dependencies were copied into this audit's `source/` directory. Both source replay modes reproduce the frozen receipt byte-for-byte; no original file was changed.

## All-order proof audit

The selected skeletons form disjoint classes. Skeletons with no crossing cycle give multiplication by B_(Q-r). For a specified root triangle, retaining one crossing cycle with its two named exterior vertices and any disjoint exterior cycle packing gives the operator X B_(Q-triangle)(X) delta. The cyclic skeleton is fixed by that named pair, so different root triangles are not counted twice. The factor2 in the odd-cycle identity cancels the parity-projector half; no additional factor2 is missing. Multiplication acts after delta, in the order used in the proof.

The omitted skeleton terms are nonnegative in the extended central-factorial basis by the explicitly cited all-basis universal-context theorem. This is a cone inequality, not an unsupported ordinary-monomial coefficient comparison. Actual input bin profiles lie in that positive cone. Evaluation at1 is its positive coordinate-zero readout, and delta B_T evaluated there is B_T(2)-B_T(0)=B_T(2) for a nonempty input. These facts give exactly the weighted local scalar inequality, including exterior empty graphs with H=1.

The first two-bin bound is a genuine injection of a Hamilton word and one of n+1 labelled cuts. Concatenation and the first bin's length recover both data. Empty slots do not create duplicates. The second bound retains the empty odd-cycle collection, whose contribution is2^n. Consequently every active step gives both H_out>=(n+2)h and H_out>=h+2^n.

A rooted context contains at least its root, so it cannot decrease order. An active one contains a root triangle and adds at least two vertices. The local inequality also makes every step nondecreasing in H, including passive steps with nontrivial exterior components. Thus before active step j the actual order is at least n0+2(j-1), yielding the displayed product bound without any fixed-word or fixed-template assumption. Adaptive choices do not alter this pointwise induction.

The integer cutoff is the largest active-step count not ruled out by that necessary product bound; it is not a claim that the maximal count is attained. When M<h0 there is no bounded construction, and when M<=1 there is no active node. The internal-order inequality uses the actual graph entering the active node. For an earlier partial graph feeding a later active node, monotonicity of both order and count gives the same necessary bound. It does not permit removal or movement of internal padding.

The weighted triangle budget follows by subtracting d_Q h and dividing by positive b. Bounding this weight also bounds the number of root triangles, since each exterior H is at least one. It does not bound exterior vertices lacking those incidences. The proof explicitly preserves this distinction, so it does not silently turn unrestricted rooted-template enumeration into a finite algorithm.

The examples and periodic consequence have the correct scope. Repeated active steps force factorial-type, hence superexponential-in-period, growth; this does not rely on parity or Hankel positivity. Passive terminal padding can leave H fixed while increasing final order, so the internal bound is not a final-order bound.

## Independent physical and arithmetic controls

`independent.py` imports no source implementation. It constructs actual rooted substitutions with explicit block incidences and computes subset start-path counts directly. Its two-bin count is the sum H(S)H(V\S) over all labelled first-bin subsets. Every literal graph is checked to be a tournament.

Normal and optimized runs agree on:

- 285 rooted labelled contexts through order4 and3,420 actual weighted comparisons, including nontransitive inputs. Of these,1,056 are strict.
- A separate six-vertex context with one root triangle and an exterior directed triangle, giving d_Q=w_Q=3 and the exact weighted equality H_out=3(h+b).
- All1,099 labelled inputs through order5, with47,096 literal Hamilton-word/cut pairs and both two-bin lower bounds.
- 192 nonuniform three-step context sequences, containing288 active steps, checking the adaptive product, actual internal order and weighted budget.
- 24,024 integer depth-boundary cases, including targets below the initial count.
- Transitive inputs of orders1–9, giving the exact wrapper counts3,5,9,17,33,65,129,257,513 and preserving the internal-padding boundary.

These are finite controls of the physical interpretation and implementation. The all-order result rests on the retained-skeleton inequality, the cut injection and the adaptive induction reviewed above.

Independent checker SHA256: `58a6f7306208c06b8e066364c5e90bf481f36d087d8398c4c273e916c54fc87c`.
Independent normal/optimized receipt SHA256: `45e64e731c6dbe915b4f88e03c21a738c576d65e8cbef27cb5814e86ce94a9ce`.

Run `python3 independent.py` or its `-O` form. The copied original replay is `python3 source/verify_active_depth_inverse_bounds.py`. All code uses the standard library only.
