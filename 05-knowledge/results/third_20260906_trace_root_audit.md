# Independent audit of the all-height carried trace and norm jets

**ACCEPTED: PROVED exact all-height boundary jets and the specified positive
trace-semiring obstruction; FINITE-EXACT controls.** This audits the
[trace/norm note](third_20260906_trace.md).

The logarithmic trace formula follows from the generic root power sums,
including the inverse-root sum and its sign. Polynomial carry cancellation
extends the identity across the singular parameter; inversion at a zero
root is never asserted. The order of every omitted term in the derivative
calculation was checked, including the top channel whose coefficient does
not vanish and whose logarithmic coefficient instead has a double zero.
The h=1 exception is necessary and explicitly retained.

The norm proof uses h distinct leading Puiseux roots after y=epsilon^h.
Their product and the characteristic-constant sign cancel as displayed.
Polynomiality then rounds the next order to an integer power of y. The
same check after deleting the carry gives a different exact order, with
positive leading coefficient. Dividing the full forced factor contributes
precisely the triangular-number sign. Every proposed semiring generator
has a nonnegative value at zero; evaluation at zero therefore already
excludes the negative norm constants. All of this takes place at a
negative-parameter continuation, outside the actual support domain.

The [independent audit source](../../04-computation/third_20260906_trace_root_audit.py)
imports no repository producer. It reconstructs quotient multiplication
traces at heights1,...,5, and norm polynomials by resultants at heights1,...,4,
with and without the carry. It checks exact leading coefficients and every
lower coefficient, full determinant deflation, the four-periodic sign,
and a literal positive-trace/negative-norm hostile. A separately constructed
rational recurrence matrix has nonzero determinant and residue29 modulo101,
so the asserted five-parameter recurrence fails over the reals as well as
over the rationals. These controls do not extend the finite root-location
bank to an all-height root theorem.

[The frozen output](third_20260906_trace_root_audit.out) records45 always-active
gates. Normal and optimized outputs agree byte for byte.

```bash
python3 -B 04-computation/third_20260906_trace_root_audit.py
python3 -B -O 04-computation/third_20260906_trace_root_audit.py
```

Raw LF SHA256:

```text
source d6ec79de0872243fff3646dac133868bfed099c59464e688689bc37d93df6d1a
output ae4203dca9d9dda3354ef88895867c93a7f88e71f0a77be2e6aea3ee215591b3
```
