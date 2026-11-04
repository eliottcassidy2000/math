# Independent audit of uniform exposed-anchor compression

**PASS.** Reviewed `../returns_next/UNIFORM.md`, sections U1–U4. The stronger query-level claim is correct: an exposed-anchor separator can be chosen individually uniform on every proper selected module. No mathematical corrections were required. This is not a proof for fully covered laminar families.

## Support-level involution

The product `s_u G_wv` enumerates ordered pairs of vertex-disjoint paths, so the starts `u,w` are different physical names. Since the endpoint `v` is exterior to the module, the second path must leave its initial module run. Swapping the two initial runs retains all support, both path slots, and the second endpoint. Each join remains valid by module homogeneity. The first remainder may be empty. Later reentry is harmless: the initial run still stops at the first exterior vertex, so applying the swap again restores the original pair. Nonempty disjoint initial runs cannot be equal. This proves fixed-point-freeness coefficientwise, before applying the support functional.

The audit separately checks a source with actual later reentry; the absence of reentry in a small exhaustive corpus would not itself verify that boundary case.

## Boolean aliasing and battery

For a singleton exterior anchor, identifying all first-port variables within a module can contribute to the mixed quadratic only through the original mixed quadratics or cubics with two first-port labels in that module. The latter sum vanishes by the involution. The shared-anchor cubic is the sum of the old shared-anchor cubics. No second-port label in another block aliases to the singleton anchor. The incoming whole-converse polynomial has degree at most three; the proof does not truncate a general ambient internal-reversal contrast.

Grouping E1 therefore preserves its right side of one. The compressed block count and outdegree must be used in the battery. The proof does so. In particular, compressed order two is valid because E1 was applied before compression to an original tournament of order at least three. The even/odd battery sizes, the endpoint correction, and distinctness of the listed masks all check.

## Lifting and scope

Uniformity makes each restored proper submodule a genuine module of the actual two-port extension, so its converse is exactly invisible by ordered path-cover reversal. Laminar reversals preserve the required modules. Other disjoint maximal selected modules are also uniformly seen by the lifted ports. The surviving whole-core contrast lifts with its nonzero residue through the ambient quotient. An exposed singleton is essential to this argument; the known fully covered cyclic blow-up remains outside its scope.

## Independent finite checks

`verify_uniform.py` imports only this audit packet's independently written tournament primitives, not producer code.

- 534 module/anchor cases across 111 source presentations, including all labelled sources of orders three and four, a declared 32-source order-five sample, larger substituted sources, and an explicit reentry control.
- 6,652 ordered disjoint path pairs checked for valid switching, support preservation, fixed-point-freeness, and involution. There are 3,290 cases with an empty first remainder and 12 later-reentry paths.
- 570 full-attachment batteries, including 80 compressed-order-two batteries and 2,506 actual query checks. The small corpus includes every disjoint proper-module packing of every labelled order-three/four source; the larger controls reach order eight.
- Every restored-module query agrees over the integers, and each battery plus endpoint control contains a genuine residue-two separator.
- Normal and `python -O` receipts are byte-identical. Tests use explicit exceptions, not assertions.

The receipt pins the reviewed proof and both independent checker files. These finite controls supplement the all-order proof and are not presented as exhaustive at larger orders.
