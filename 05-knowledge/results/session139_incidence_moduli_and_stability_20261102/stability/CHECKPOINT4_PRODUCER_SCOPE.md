# Final stability-lane selector for Session139 checkpoint4

## Accepted results

1. **All-order absolute margin for matching uncertainty.** With r disjoint
   uncertain pairs, q_i=p_i(1-p_i), and arbitrary binary cross orientations,

       det(I+A)-abs(det A) >= product_i(1-q_i)       when n>2r,
       det(I+A)-abs(det A) >= product_i(1-q_i)-product_i q_i
                                                        when n=2r.

2. **All-order nonlinear matching gap, with the binary baseline retained.**
   For independent binary rounding T and B=E Delta(T),

       Delta(A)>=B+1-product_i(1-q_i)-1_{n=2r}product_i q_i.

   The primary proof is `../tensor_matching/PROOF.md`. The independent review
   and arithmetic second proof are selected from this directory. The arithmetic
   proof is specifically an integer-parity/Eisenstein/F4 argument; it does not
   establish an arbitrary-modulus analogue or a Riemann-hypothesis result.

3. **Sharp stability for every n>=140.** `PROOF_N140.md` combines the new
   nonlinear matching input with the accepted shared wedge budget and dual
   certificate. The exact constant is

       kappa_(2R)=(4/R)[1-(3/4)^R-(1/4)^R],
       kappa_(2R+1)=(4/R)[1-(3/4)^R],    R=floor(n/2).

   Positive-variance equality is exactly ordinal fair pairs, with one singleton
   in odd order. Relative to the initially unresolved10..1499 range, this closes
   1360 orders,140..1499, and leaves10..139 (130 orders) unresolved. It is NOT
   an all-order sharp-stability theorem for arbitrary uncertainty graphs.

4. **Two exact local-factor hostiles.** The conditional factor1-2q fails with
   fractional cross edges at n5; the conditional factor1-q fails even with
   matching uncertainty at n8. These remain valid counterexamples. Their
   historical note left the absolute matching floor open; result1 above now
   resolves that formerly open absolute floor affirmatively. Conditional and
   absolute statements must not be conflated.

## Audit and freeze state

The primary tensor theorem and arithmetic second proof are independently
accepted in `../tensor_matching_audit/{AUDIT.md,ARITHMETIC_AUDIT.md}`. The
historical phrase in the independent review saying exact controls were being
checked is superseded by these completed audits. The n140 proof is accepted
in `../independent_threshold140/AUDIT.md`, with a final status-only review
receipt if supplied by that audit. The local-factor controls were independently
recomputed in `../independent/audit_local_factors.py` and its matching normal
and optimized outputs.

`CURATED_CHECKPOINT4.json` is the producer selector: it pins only the selected
new producer documents, verifiers and receipts, and separately identifies the
independent evidence. Earlier n>=748 and n>=526 packets remain untouched and
belong to earlier checkpoints. Float probes and research-only files are
explicitly excluded. `verify_curated_checkpoint4.py` verifies the selector's
hashes and lengths without depending on a particular absolute workspace path.
