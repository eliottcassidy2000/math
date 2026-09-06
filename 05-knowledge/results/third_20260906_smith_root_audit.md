# Independent full-matrix audit of every mixed (m,2,1) Smith form

**ACCEPTED: PROVED in the stated all-m affine and projective scope;
FINITE-EXACT independent full-matrix controls.** This audits the
[mixed-bank theorem](third_20260906_smith.md). No producer code is imported.

The analytic audit checked both ideal containments: each of U,V,W is an
actual second minor, and the other six minors are integral combinations.
This retains the entire cokernel after clearing the complete heavy bank.
The valuation cases exhaust A<=B and A>B, including primes dividing m+1
or m+2, zero linear numerators, and the resonant boundary d=v_p(m).
The cap min(e,md) retains both competing minors. Plucker's identity
preserves the normalized unit ratio modulo p^e, exactly the digits used
by the cap. The residue ladder, its dyadic missing zero state, and the
precision-kernel interval follow with the stated equality boundaries.

Projective transport is used only for complete Hasse banks. At two,
three occupied residue classes make every pair bracket a unit, so the
exceptional chart case is already a unit observer. Multiplicity labels
cannot be permuted as if equal. The incoming cardinal-column primitive
content theorem justifies value-datum attainment; the independent inverse
controls below test the claimed source perturbation, including unchanged
heavy jets. No result for residual rank four is inferred.

The [audit source](../../04-computation/third_20260906_smith_root_audit.py)
constructs literal homogeneous Hasse matrices and calls SymPy's integer
Smith normal form on the whole matrix, independently of the producer's
residual DVR elimination. Its [frozen output](third_20260906_smith_root_audit.out)
records **120,392 always-active gates**, including 1,750 full integer
matrices, 7,000 full prime partitions, 105 projective controls, eight
literal inverse ownership controls, and 104,106 adjacent-kernel precision
checks. The affine bank has m=2,...,12 and distinct nonzero a,b in [-6,6],
augmented by deep resonant cases for (m,p)=(3,3),(4,2),(5,5),(6,2),
(8,2),(9,3),(10,5),(12,3),(16,2). Projective controls include infinity,
unimodular coordinate changes, altered tangents and all three dyadic
residue classes. Finite controls challenge rather than replace the proof.

```bash
python3 -B 04-computation/third_20260906_smith_root_audit.py
python3 -B -O 04-computation/third_20260906_smith_root_audit.py
```

Normal and optimized output agree byte for byte. Raw LF SHA256:

```text
source f5a8b8336002428a57d2a751bb1ffc8447d5332c93ed359a21c2e69c27c10254
output edcd37f22165277cafe768653f7a01666e922de16cb6ff2430997f405133da11
```
