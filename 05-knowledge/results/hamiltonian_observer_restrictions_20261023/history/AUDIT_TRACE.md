# Independent review: guarded return sampling and actual clocks

Verdict: **PASS** for the stated finite reversible-chain theorem and auxiliary
tournament lift. No mathematical correction requested. The native
port-preserving move realization remains explicitly **OPEN**.

Reviewed `../trace/TRACE_WITH_CLOCKS.md` and `../trace/check_trace.py`.
This review does not reprove the incoming full-cube spectral theorem; its
use as an inherited premise is explicit. No Git or routing changes were made.

## Proof and implementation checks

1. The weighted trace formula and harmonic extension have the correct order:
   `B=(I-Q_CC)^-1`, `R=Q_SS+Q_SC B Q_CS`, `h_C=B Q_CS f`.
   Irreducibility plus nonempty acceptance ensures the inverse exists.
   The identity `E_Q(h)=pi(S) E_R(f)` has the required normalization.
   The variance bound follows from conditioning on `S`; no uniform-law
   assumption is hidden in it.
2. The trace gap argument proves the asserted lower bound for the Poincare
   gap. Separately, weighted symmetrization makes the direct principal block
   PSD, and the excursion term is `A B A^T` with positive definite `B`.
   Together these give nonconstant spectrum in `[0,1-gamma]`. On a singleton
   the statement is correctly treated as vacuous.
3. The duration generating matrix counts positive returns: direct transitions
   carry `z`, excursions carry `z^2(I-zQ_CC)^-1`. Differentiation at one
   gives the displayed mean-clock vector. Kac's conditional stationary mean
   is `1/pi(S)`, not one ambient step. Stationary conditional start is
   explicitly required for the averaged formula.
4. The rational PSD routine is a valid unpivoted Schur-complement test.
   A positive pivot permits exact elimination. A zero diagonal with any
   nonzero remaining off-diagonal entry is impossible for a PSD matrix and
   is correctly rejected. All comparison matrices are symmetric. The gap
   matrix `(1-gamma)D + gamma pi pi^T - D R` is the correct full-space
   expression of the mean-zero eigenvalue upper bound.
5. The compressed order-four owner computation keeps actual owner orders.
   If `D_yy=1-k_y/24`, rejected intermediate owner choices multiply on the
   right of `Q`, so `F=(I-QD)^-1 Q` is the sum of all rejected-word
   excursions. The final particular accepted owner has probability `1/24`.
   The clock equation is `t=1+QD t`, hence `t=(I-QD)^-1 1`.
   Neither formula loses a factor of 24 or replaces endpoint acceptance by
   acceptance of a whole auxiliary word.
6. Literal tournament reconstruction and converse transport retain actual
   owner labels. All nonempty accepted subsets are tested in the four
   isomorphism representatives, with Hamiltonian-list counts `1,3,3,5`.
   Thus the subset count is `1+7+7+31=46`. Actual relabeling transports
   these fibres bijectively without changing the auxiliary word kernel.
   Indexed pins and any exact guard define some tested subset; no particular
   unproved native mod-8 formula is inserted.

## Reproduction and independent controls

Both normal and optimized executions of the source passed:

```text
python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/trace/check_trace.py
python -O -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/trace/check_trace.py
```

Each run reports 255 nonempty subsets of the eight-word auxiliary carrier,
seven nonuniform subsets, and 46 actual owner subsets. Harmonic energy,
return clocks, rational PSD, sharp ambient gap, rejection-versus-return,
source reconstruction and converse checks all passed. The author script
reproduces `trace_receipt.json` when run.

Two additional independent controls were evaluated:

- For `Q=(I+J_pi)/2` with `pi=(1/2,1/3,1/6)`, conditioning on the first
  refresh gives the closed form `R=(I+J_(pi|S))/2` for every nonempty `S`.
  Comparison with all seven computed traces passed. This argument is
  independent of Schur elimination and detects an incorrect uniform
  stationary law. Its return-clock mean is `1/pi(S)`.
- Direct enumeration of all 64 labelled four-vertex tournaments and their
  24 literal owner orders gave Hamiltonian-count histogram
  `{1:24, 3:16, 5:24}`, totaling 192 Hamiltonian presentations. This agrees
  with the lift universe `8*24` and the representative fibre sizes.

## Scope boundary

The paper does **not** assert that an ambient excursion preserves the source,
pins, ports, or native guard at intermediate states. It expressly separates
endpoint acceptance from a native move-lifting theorem and retains elapsed
ambient time. The three objects `R`, `R(z)`, and `B Q_CS` have distinct roles:
return targets, duration, and harmonic extension. The comparison with the
earlier electrical odd-denominator hostile is methodological, not an
identification of the two operators. The final obligation correctly requires
legal native paths, congestion, and integer lengths before transferring a
native spectral or wiring conclusion.
