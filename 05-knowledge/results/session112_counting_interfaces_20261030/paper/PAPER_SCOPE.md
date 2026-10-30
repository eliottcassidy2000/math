# Factoring inspiration: exact scope and reusable mechanisms

## Source inspected

Isaac Rudich and Louis-Martin Rousseau, *Deterministic Factoring in Time N^{o(1)}*, [hexagon:2610.00600v1](https://hexagonmath.org/pdf/2610.00600v1), six pages. The source's [Palomar v1 record](https://palomar-registry.org/entry?id=PALOMAR-2026-10-28-000015&version=1) pins [commit 5f5ff8f3823fa8fef8ee02a6f29049db3c76a0ce](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/5f5ff8f3823fa8fef8ee02a6f29049db3c76a0ce). Its Challenge.lean and Solution.lean were inspected directly. The registry reports successful Comparator checking and proof checking by Lean, nanoda and con-ron. Those checks were not independently rerun here. Registry review is AI-mediated, not human peer review or a novelty certificate.

The stated theorem quantifies **one machine before every accuracy parameter**: there exists M such that for every positive integer k some C_k bounds its runtime by C_k 2^{floor(n/k)}. The machine computes Nat.minFac with binary input/output and finite stack alphabets. The stronger explicit rate exp(O(log N/sqrt(log log N))) is a proposition named as a conjecture, not a proved theorem. The Challenge's placeholder is deliberate: Solution supplies Aux.subexponentialRate. A placeholder in the challenge is not evidence of an admitted proof in the solution.

The paper reduces minimum-factor search to prefix divisor counts, realizes multiplication as uniquely extended local constraints, and uses a re-proved network-counting engine. It computes a scaled integer count through modular observations. These are useful mechanisms independently of adopting its full complexity theorem.

## Asymptotic precision

The phrase that N^{o(1)} is slower than every L_N[alpha], alpha<1, is not a consequence of this asymptotic class alone. Both poly(log N) and exp(sqrt(log N)) belong to N^{o(1)}. The class includes fast and slow functions. An explicit bound exp(C log N/sqrt(log log N)), for fixed positive C, does grow faster than each exp(c(log N)^alpha(log log N)^{1-alpha}) for fixed 0<alpha<1,c>0: the ratio of exponents is (C/c)(log N)^{1-alpha}/(log log N)^{3/2-alpha}, which diverges. An upper-bound comparison alone is not a runtime lower bound for a particular machine.

Uniform subexponential bounds do absorb polynomial overhead, but the parameter should be selected carefully. Suppose T(n)<=C_k 2^{floor(n/k)} for every k and runtime is at most D n^a (T(n)+1)^b. For a desired target r, choose k>=2br and absorb n^a into a constant times 2^{n/(2r)}. The result is O_r(2^{n/r}); increasing the input parameter k offsets the slowdown. The theorem does not supply an efficient practical crossover.

## Exact older bridge: signed octonion multiplication

The earlier repository note [Paley/Fano/octonion design](../../collatz_mod6_20260922_paley_fano_octonion_design.md) constructs octonion multiplication from the Paley tournament on seven points plus a chosen Fano line system. The multiplication signs and third-point outputs together matter. The factoring paper's use of an eight-square identity therefore meets an actual old algebraic object, not merely the number eight.

There is an elementary obstruction to dropping signs. In the two-dimensional subalgebra with basis 1,i and i^2=-1, set x=1+i and y=1-i. Then xy=2 and ||xy||^2=4=||x||^2||y||^2. Replace the multiplication signs by positive support, so i*i=1. The same x,y now multiply to zero. Thus unsigned multiplication support does not preserve nonzeroness, the norm identity, or the ability to use a signed sum as an existence test.

Conversely, coherent basis-sign changes preserve norm composition. For any orthogonal diagonal sign matrix D with D(1)=1, transport multiplication by x*y=D^{-1}((Dx)(Dy)). Then ||x*y||^2=||Dx Dy||^2=||Dx||^2||Dy||^2=||x||^2||y||^2. This is a gauge transformation, not arbitrary sign deletion. It explains why some changed orientations preserve the algebra while unsigned quotients can lose crucial information.

No 11/12 exponent saving follows merely from the old Paley construction. That saving needs the precise crossing decomposition, tensor interfaces, grouping overhead and recursion from the counting engine. Nor does any of this establish a new general factoring algorithm or a Collatz theorem. The usable session principle is narrower: retain the signs, shared interfaces, normalization and positivity conditions required by each observation-to-search step.

## Incoming formal-engine scale caveat

The later [source audit](../../scaled_count_fft_carriers_20261030/SOURCE_AUDIT.md)
and [theory C2](../../scaled_count_fft_carriers_20261030/THEORY.md), inspected
at commit 05b452976f450d06a705c11b30718f35173e11a8, report that the generic
Lean counting interface has multiplier `(t+1)^A * 2^(t*R)`, not a purely
dyadic scale. Our odd-modulus corollary applies to its explicitly promised
`2^E` normalization; it is not automatically a speedup of that full engine.
The additional factor must be retained in the general effective-modulus
analysis. The incoming packet also reports a local Lean build at a later
source revision; that is distinct from our inspected revision above, and
we did not independently rerun its build.
