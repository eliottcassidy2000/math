# Independent audit: leaf cofactors and unit-free maximum-endpoint entry

**Verdict: PASS.** The analytic arguments, all six actual-entry controls,
complete mixed-support scope, explicit integer cutoffs, and strict-versus-weak
boundary distinctions in
[second_20260906_entry.md](second_20260906_entry.md) are accepted without
mathematical repair. The tree bound is unconditional; the LRC conclusions
remain relative to the named inherited suppliers and their actual-entry
domain. Automatic balanced closure and general unit-free entry remain OPEN.

## 1. Full analytic review

I read the entire candidate and its primary source. The inherited actual
crossing/weighted-kernel, finite-box internal-height, and phase-grid arguments
in **THM-3818**, `01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
were already read during this branch's entry work. The completed audited
[overnight12 signed-box result](overnight12_20260906_lrc_gcd_semigroup.md),
[overnight15 larger-unit result](overnight15_20260906_lrc_larger_unit.md),
and the current joint-shadow caps are the explicit inputs. Their reserved
namespaces, if any, are not treated as proved suppliers.

For the tree inequality, each side of an edge cut contains at least one
original leaf: following the finite tree away from the cut eventually
reaches a leaf. Thus the leaf-cut fraction really lies in
`[1/ell,1-1/ell]`. Averaging the logarithms of the leaf-rooted cofactors
makes the edge contribution `a^alpha b^(1-alpha)`, up to interchange.
Weighted AM-GM gives the claimed maximum under `a+b<=C`. The derivative
of `log phi(alpha)` is `log(alpha/(1-alpha))`; symmetry and monotonicity
therefore put its maximum at the interval endpoints. Since a tree with
`n>=3` has at most `n-1` leaves, the uniform constant follows. The
leaf-dependent expression and its fractional exponent are correctly typed.

The star equality construction checks exactly. Its sharpness belongs to
positive-real edge coefficients and their cofactors. Integer common-content
division can reduce the primitive minimum to one, and repeated leaf labels
are allowed in that equality construction. No sharpness claim is therefore
transferred to distinct primitive integer speed shapes. The inequality
still transfers in the correct direction by dividing positive integer
cofactor content. The separate two-vertex integer bound `177` correctly
uses positive distinct primitive labels and sum at most 356.

For entry, all coefficients and scales retain their physical meanings.
The internal pair height `B<=Q` permits the exact signed-box lemma. When
`c=gD/gcd(gD,tv)<=Q`, a point in the central interval produces a genuine
support-at-most-three mixed relation. Its distinguished coefficient is
positive, so its partial sum on that component is nonzero. Actual equality
entry excludes it and forces `tv/delta>R>QL/D`. The direction is necessary
under entry, then sufficient after comparison with the full larger maximum.
No collective gcd is used as a pair gcd, and no selected pair maximum
replaces the full phase Lipschitz constant.

The native criterion forces `t>7(b+1)L/a`. Coprimality of the two physical
scales turns the preserved smaller-component lifts into a complete translated
grid in the larger clock. Its spacing is strictly less than the protected
open interval length, yielding strict clearance for both components. The
integer refinement correctly keeps `delta|tv`: writing
`delta=z*d,v=z*e`, with `gcd(d,e)=1`, gives `t=d*h` and
`h>=floor(R/e)+1`. The refinement is the exact lower bound under those
retained arithmetic conditions, not a claimed classification of the
original gcd identity at every candidate scale. A closed interval handles
the stated weak equality boundary.

In a hypothetical strict failure, the larger physical component's gcd is
`g<=C_b` by the current subset caps. Hence `D<=Q/C_b` pays the cleared
coefficient gate. The second uniform condition gives the required strict
scale separation. For `g>C_b`, the inherited conclusion is only weak
safety, and the candidate correctly avoids extending strictness to that
branch. All five displayed cutoffs and the balanced zero output agree
with independent exact arithmetic. The `2+11` cutoff is explicitly credited
to the inherited endpoint theorem.

## 2. Independent controls and consequences

The [independent verifier](../../04-computation/second_20260906_entry_audit.py)
imports no producer. It enumerates labelled trees by taking connected
`(n-1)`-edge subsets of the complete graph, rather than by Prufer codes.
It computes each rooted cofactor as an actual maximal minor using integer
Bareiss elimination, rather than multiplying cut-selected edge coefficients.
The finite universe contains all labelled trees on three, four, and five
vertices with all positive integer edge coefficients of sum at most six,
five, and three, respectively: **26,800 weighted tree cases**. Every case
checks the kernel relation, uniform bound, leaf-count bound, and primitive
normalization direction.

Separate exact star controls at `n=3,...,10` attain the real/cofactor bound
and show its primitive-content loss. A five-vertex star refutes the tempting
`(C/2)^(n-1)` uniform minimum bound. A three-vertex asymmetric star shows
why the result cannot be relabeled as a maximum-cofactor bound.

The verifier independently regenerates allowed decoder sums by a prime
sieve and products of inert prime powers. For each of the six supplied
nonunit larger shapes, it reconstructs the actual graph and the physical
box. It verifies all **1,001 mixed supports** across the six rows using
a complement-semigroup membership routine, independently controlled on
**3,708 literal signed-box points**. Dominance excludes the full mixed
relation span, while connected internal graphs provide rank eleven. This
is actual equality entry, not a rank calculation on a selected family.
All larger shapes lack one, and every pair with their full maximum has gcd
greater than one.

Each native condition is evaluated on its actual physical row. An independent
nearest-grid-point selector then constructs a full thirteen-speed phase;
its clearance is checked against every physical label and is strictly above
`1/14`. The selector differs from the producer's first-interior-grid-point
construction. A further **4,320** exact arithmetic cases verify the
integer divisibility refinement and its first possible scale.

The separate
[balanced hostile](second_20260906_decoder.md) remains consistent with this
stronger tree theorem: its minimum exceeds the balanced entry cutoff but
lies below the improved six-vertex cofactor bound. That hostile forbids
inferring the failed sufficient balanced conditions from actual entry and
the retained gcd profiles; it does not challenge this native criterion.

## 3. Reproduction

    python3 -B 04-computation/second_20260906_entry_audit.py
    python3 -B -O 04-computation/second_20260906_entry_audit.py

Both runs pass **207,160 always-active gates**, with byte-identical LF
[frozen output](second_20260906_entry_audit.out). Source SHA256:
`68bb61eb1f4f7ae8065075ea9ab26c6e61f97a57c57be2b77daa6052a32de1cb`.
Output SHA256:
`8f36195a7d3c6fb8a39be07c522e1c2dd0e72cfd8140d13f6b84d837a7dc1ed0`.

## Incoming cross-divisor integration — independent acceptance

**PASS.** I read the complete independently audited incoming
[open_frontier_sep06_decoder.md](open_frontier_sep06_decoder.md) at
`origin/main` commit `8e560f2142` through `git show`, including its native
and cross-divisor proofs. Its native criterion overlaps the present entry
mechanism; the leaf-cofactor bound is the additional input here. Taking
`v=min V`, `d=gcd(D,v)`, and `H=D/d`, one has `d|delta` and exactly
`lcm(D,v)=Hv`. Every numerical cutoff in the present table is therefore
sufficient for **H**, under the same actual-entry hypotheses: if `M_a` is
the proved bound on `min V`, the table ensures `C_b H<=Q` and
`H v<=H M_a<=aQ/[7(b+1)]`, which are precisely the incoming cross-divisor
conditions. The original cutoffs on `D` remain valid weaker subclasses
because `H<=D`. This combination preserves the inherited weak conclusion
when `g>C_b` and the strict construction in the remaining branch. The two
balanced controls have `d=1` for all six endpoint partners, so this incoming
cancellation does not change their recorded failures. This is an analytic
integration check; the frozen independent computation and its stated finite
universe above are unchanged.
