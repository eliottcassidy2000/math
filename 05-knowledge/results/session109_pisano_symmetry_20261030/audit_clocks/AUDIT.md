# Session109 independent audit: observer clocks and affine carries

Status: mathematical claims verified within the scopes stated below. This is a proof review and a finite independent cross-check, not a novelty claim, general modular-data saturation proof, or Collatz convergence argument. No publication was performed by this audit.

## Findings and corrections

1. The first clocks draft inferred a trace-loss factor in `{1,5}` merely from injection into a five-element set. Injection alone proves at most five. The producer corrected the argument: first state the size bound, then obtain the exact alternatives from prime-power and CRT calculations. Independently, the trace radical here has square zero, so the invisible-shift map is additive and supplies another valid repair.
2. The first clocks checker used Python `assert`; these disappear under `python -O`. The producer converted them to explicit exceptions. This independent checker uses explicit `check` calls throughout; normal and optimized outputs are identical.
3. Minimal sidecar cardinality means **one present-time sidecar value added to the entire trace history**, not the entire time history of the sidecar. For the latter interpretation a binary indicator of one point in an orbit already has the full orbit clock, so no five-symbol lower bound holds. Resolved: both five-valued and general p-valued statements now explicitly restrict this lower bound to one contemporaneous sidecar value.

## Proof review

### Fibonacci observer clocks

The algebra has basis `(I,B)` for both `B=Q,Q²`, since the off-diagonal entry is one. The adjacent-trace matrices have determinant five and Smith factors `(1,5)`. Thus their radicals have size `gcd(N,5)` and adjacent traces determine the complete recurrence history. The marked-vector proof uses a unit determinant, correctly avoiding the false implication that nonzero determinant is enough over a ring.

The local lifting proof at five is valid: `Q^20=I+5D` with `D mod5` nonzero and in the trace radical, and fifth powering raises exact five-adic divisibility by one. The initial `Q^4-I mod5` is nonzero square-zero; hence the initial order is twenty. The exact trace clocks follow from an exhibited order-five invisible subgroup and the kernel-size bound. CRT then gives the stated formulas and the masking criterion `v_5(ord_u B) >= e`.

The scalar/signed distinction is essential and correct: `A^10=11I mod15`, whereas `-I` never occurs. Scalar observations retain separate local scalars; a global sign requires compatible local signs.

For the sidecar, `5b=2U-tau*T` determines b modulo `N/gcd(N,5)`. The top base-five digit at the five-primary component resolves the remaining ambiguity; the stated a-reconstruction formulas are correct. This is injective on the whole algebra. The pigeonhole lower bound holds exactly in the present-time sense above, and only on cyclic histories actually having five-point fibers; CRT-masked histories need none.

### Ramified metallic extension

Under `p>=5` and `v_p(a²+4)=1`, `delta=2mu-a` satisfies the Eisenstein relation `delta²=a²+4`; `O=Z_p[delta]` is a ramified quadratic DVR. The trace radical is exactly `delta^(2k-1)O/delta^(2k)O`: direct trace containment and cardinality prove equality without an appeal to an unspecified different theorem.

The identity `mu^4-1=a*mu²*delta` is exact. In the pth-power binomial expansion the first term has valuation `s+2`, every middle term has strictly larger valuation, and the last has valuation `ps>s+2` for `p>=5`. Consequently `v_delta(g^d-1)=1+2v_p(d)`. Solving full-return and trace-return inequalities gives `p^floor(r/2)` and `p^floor((r-1)/2)` respectively. The coefficient of mu in g is `a(a²+2)=-2a modp`, a unit, so two g-traces really span the full trace pairing. Zero is correctly separated. The exclusion of p=3 and ramification valuation greater than one is substantive.

### Affine obstruction and orders

All three principal affine theorems pass review. For odd stride the entry gcd divides two by determinant reduction, and modulo-two order three settles the exact gcd. At even stride the two factorization identities reduce Smith form either to a scalar unimodular matrix or to `F_m*(2Q-I)`, with Smith factors `(1,5)` for the latter fixed factor.

Translation conjugation and the image criterion have consistent signs. Smith reduction counts cokernel classes and fixed points, but correctly does not purport to identify a particular carry class without a basis transformation or direct image test. The global affine order law follows from the first linear return, after which all iterates are translations. It concerns the permutation of the whole board, not every point orbit. CRT gives an lcm, not a product.

Both hostile examples are correct: stride five modulo eleven permits orders two and twenty-two with the same linear part; stride four modulo five has nonremovable carries but unchanged order five, since the return-carry map is zero. Therefore neither obstruction-group size nor the presence of a nontrivial obstruction alone determines order growth.

## Independent computation

Run `python3 check_independent.py` and `python3 -O check_independent.py` from this directory. No producer modules are imported.

- Q/A clocks: moduli 2..100 using Fibonacci-pair recurrence and KMP word-period computation, distinct from producer matrix-power enumeration and exhaustive shift tests.
- Whole-algebra trace-plus-sidecar injectivity and radical square-zero: every element modulo 2..40 for Q and A.
- Affine orders: every carry, strides 1..10, moduli 2..11, **5,050 cases**. Orders come from explicit permutations of all board points and lcm of cycle lengths, not affine-map iteration.
- Closed Smith factors: strides 1..300 via independently generated Fibonacci entries.
- Ramified extension: every state at `(p,a)=(5,1),(5,4),(13,3),(17,8)`, exponents k=1,2. Valuation comes from the integer norm `c²+a*c*b-b²`, not producer ideal enumeration or selected layer representatives; full cycles and KMP trace periods check both clocks.
- Total: **151,516 explicit checks**, all pass in ordinary and optimized Python, with byte-identical JSON output.

`results.json`, `run_normal.json`, and `run_optimized.json` preserve the exact tested domains. These checks support but do not replace the unrestricted proofs.

Independent audit reran all three producer checkers (`clocks/check_clocks.py`, `clocks/check_metallic.py`, `affine/check.py`) under `python3 -O`; all exited zero.

## Final synthesis review

`SYNTHESIS.md` §4 is correct: if full clock m=5t and trace clock t, the t-step action cyclically permutes each five-point fibre. Cyclic word `00011` gives the five distinct length-three windows `000,001,011,110,100`; two binary readings have at most four outcomes. The conclusion is explicitly restricted to a designed predicate and fixed-length serial observations, not a free implementation or a snapshot lower bound. The checker correctly realizes this on all twenty Q-powers modulo five, partitioned by two traces into four fibres.

The prime table gives `(p,pi(p),stride,bits)=(3,8,8,4),(5,20,4,4),(11,10,2,12),(19,18,2,20)`, agreeing with independent arithmetic. For the no-reversal claim at eleven and nineteen, the implicit reason is that any odd power alleged to reverse has an odd power equal to the half-order involution. That involution has determinant minus one and therefore a nonzero fixed vector, impossible for a tournament reversal fixing zero; even powers have odd order and cannot reverse. Scalar minus identity is not required for arbitrary reversing maps, but this argument validates the stated specific cyclic examples.

The odd-prime ramified hypothesis in §2 already implies p≥5: three does not divide a²+4 for any integer a. Thus it does not accidentally extend the DVR/valuation theorem to an excluded case. The composition and Collatz scope limits remain explicit; the tournament attachment theorem itself was audited separately, not independently re-proved here.

`synthesis/check.py` passes both ordinary and optimized Python with byte-identical outputs, preserved as `synthesis_normal.json` and `synthesis_optimized.json`. All known clock/carry audit items are resolved.
