# Independent audit: Session139 threshold n >= 748

**Conclusion: PASS.** The proof establishes `K_n = kappa_n` for every integer
`n >= 748`, including the inherited equality classification, using the
accepted Session138 inputs. It requires neither a new graph census nor the
invalid termwise Pfaffian unit-coefficient claim. It does **not** establish
orders 10 through 747.

This audit independently reconstructed all new inequalities before reviewing
the producer's final proof. All stated comparisons below use exact integer
or rational arithmetic. Uniform extension is mathematical, not a scan.

## Reviewed sources

Producer files, relative to the Session139 research root:

- `stability/PROOF_N748.md`, SHA-256
  `b34365ee7b3748969432be6ec0d430d1b4461385cddbfd9d864df69cc9c317f7`.
- `stability/verify_n748.py`, SHA-256
  `9fd557a3ca6a8c3967cd161b3fffda3d73bf81733752bce6b02880bdb3fa38d1`.

Accepted inputs were read from the Session138 source tree:

- `threshold_refinement/PROOF_N1500.md`, SHA-256
  `fec495e0e7845b723b8d964eeb242bc537119d1baccbb683b2dd2b1000227636`.
- `asymptotic_stability/PROOF.md`, SHA-256
  `af7a86c2390b5e7507c1c689fc4cd2c0d3a191c2f6019207eed3b7c9e472dae0`.

Both lie below
`05-knowledge/results/session138_stability_and_actual_contexts_20261102/`.
This audit accepts their previously reviewed matching, rounding, binary,
and equality theorems; it does not claim to reprove all of Session138.

## Mathematical review

1. **Bootstrap is comfortably strict.** For positive V, g(V)/V decreases:
   its derivative has the sign of `(4V+1)exp(-4V)-1 < 0`, by
   `exp(4V)>1+4V`. The partial exponential sum `16319/315 > 50`
   proves `g(1)>49/200>8/747`. Thus candidate configurations have V>1,
   then `V>49(n-1)/1600`. At n=748 this is `36603/1600>4`.
   Consequently
   `g(V)>(1-50^(-4))/4 = 249999/1000000 + 3/3125000`.
   No small-variance branch is lost.

2. **The retained B is used correctly.** From
   `kappa V >= B+g(V) > D^2/(2n)+c`, one obtains
   `V-D > c/kappa + (D-n kappa)^2/(2n kappa)-n kappa/2`.
   All divisions are by positive quantities. This genuinely improves the
   old separate bounds; it does not infer an upper bound on B from a lower
   bound on B.

3. **A nonempty heavy matching is established before use.** Its mass exceeds
   `c(n-1)/8 - 4n/(n-1) - n/110`. After subtracting n/60, the endpoint
   surplus at 748 is `428091991/5976000000>0`. The derivative is
   `1449967/264000000 + 4/(n-1)^2>0`. Thus its size exceeds n/15.
   Under the temporary assumption L>=1, `r_h <= floor(n/2)-L` also proves
   `n>2r_h`, so the coefficient a_(r_h) is applicable for both parities.

4. **The exponential tail is uniform.** Raising its endpoint comparison to
   the fifteenth power gives the exact true integer inequality
   `(2000*748)^15 * 3^748 < 4^748`.
   Since `log(4/3)>1/4` (integrate 1/t on [1,4/3]), the logarithmic
   derivative of `n(3/4)^(n/15)` is less than `1/n-1/60<0` throughout.
   Hence the real exponent n/15 and nonintegral tail cause no gap.

5. **The light-count estimate is parity-safe and strict.** Multiplication
   by the positive integer `R-L` is legitimate. Using `kappa R<=4`,
   `V<=n/7`, and the preceding tail gives `4 epsilon V<=1/3500`.
   Because L>=1 and `kappa V>c`, the result is strictly
   `L<CD`, where `C=308000000/13627923`. In particular D>0 before any
   division by D occurs. No parity error is hidden in replacing R by n/2.

6. **The combined wedge bound does not double-count.** Distinct greedy
   groups count disjoint unordered edge pairs. Heavy contributions can
   therefore be added to the Cauchy-Schwarz light contributions:
   `W >= eta D_h + D_l^2/(4L) >= eta D_h + D_l^2/(4CD)`.
   With `x=D_h/D`, completing a square gives
   `eta x+(1-x)^2/(4C) = h + (x-(1-2C eta))^2/(4C)`.
   Here `h=8027923/749535765`,
   `1/2-C eta=2427923/27255846>0`, and
   `h-8/747=190787/186634405485>0`.
   Thus `W>=hD>kappa D>=B>=W` is indeed impossible.

7. **Boundary and equality cases are covered.** V=0 is irrelevant to the
   defining infimum, whose domain is V>0. Under L>=1, D=0 is excluded by
   `L<CD`; D_l=0 is valid in the weak combined wedge inequality; and B=0
   cannot escape the strict contradiction. When L=0, eta>kappa and any
   D>0 gives `W>=sum w_i S_i>kappa D`; hence D=0 and B=0. Positive V
   ensures the greedy matching was not initially empty. All supported
   binary completions then have zero nonnegative gap and are transitive
   by the accepted n>=4 theorem. The accepted matching/H=1 minimization
   supplies precisely the stated fair-pair equality classification.

## Executable controls

`audit_n748.py` is an independent standard-library certificate with explicit
exceptions, never removable `assert` statements. It checks 22 conditions,
including coefficient-wise square identities and both parity bounds.

- `normal.json` and `optimized.json`: PASS; byte-identical.
- `negative.json` and `negative_optimized.json`: a deliberately false
  extension to n=747 is REJECTED, exit status 1, by
  `h-8/746 = -3727781/279576840345 < 0`; byte-identical.
  This is a control of this certificate, not a disproof of the theorem
  at order 747.
- `producer_normal.json` and `producer_optimized.json`: independently rerun
  producer certificate; PASS and byte-identical.
- `producer_mutation.json`: replacing h by zero is rejected by the producer
  certificate in both normal and optimized Python, exit status 1.

An initial verification-only defect was found: the producer used `assert`,
so Python -O disabled its comparisons. The producer replaced these with
explicit RuntimeError checks; the corrected script is the pinned source
above. No mathematical defect was found, and no Git or remote mutation
was performed by this audit.
