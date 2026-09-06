# Independent analytic audit of the all-R completed Hamiltonian theorem

**Status: ACCEPTED ANALYTICALLY.** This audits the universal polynomial
`R(p,y)` statements in
[the completed Hamiltonian report](planar_jc_long_20260906_hamiltonian.md).
The accepted objects are the full actual depth modules, their closures
in the specified `D`-adic completion, and the separately stated logarithmic
completion. No polynomial specialization or termination is inferred.

The generator formulas follow directly from the inherited Poisson brackets
and the product rule for `S=c+p²D R`. I independently checked both signs
in `delta p=-D B_R`, `delta y=D A_R`, and the cancellation giving
`delta D=-D²(4yR+2pyR_p+3p³R_y)`. The expressions for `delta x` and
`delta u` are polynomial in `p,y`; rewriting `p³=D(1+u)` and `y²=Du`
also makes each divisible by `D` in the actual ring `B_2`, with the
quotient of depth at most one. These two different divisibility statements
must be distinguished, as they are in the report.

The key depth-lowering step is sound: `Dx=yp` and `Du=y²` absorb a factor
`D` into any positive-depth factor. Differentiating an `x` or `u` factor
lowers its exponent count, while differentiating the coefficient in
`K[p,y]` supplies exactly such a factor `D`. Consequently
`delta P_d` is contained in `P_(d-1)` for `d>=1` and `delta P_0` in `D P_0`.
The same product rule gives the separate uniform estimate
`delta(D^N P_d)` contained in `D^(N+1)P_d`. Iterating proves both the
claimed depth/cusp estimate and the uniform topological nilpotence on the
entire `B_2/D^N B_2` quotient. The nonzero leading action for `R=1` on
`x^d` verifies that the one-step depth bound cannot be improved uniformly.

The quotient exponential is a genuine algebra automorphism. Its sum is
uniformly finite on each quotient, and the derivation's Leibniz formula
gives multiplicativity. The inverse is the exponential with the opposite
parameter, and the binomial formula gives the additive group law.
Compatibility with reduction therefore gives inverse continuous
automorphisms of the full inverse limit, not just a map defined on
polynomial inputs. Since `B_2` is Noetherian, the reduction kernels agree
with the extended `D`-adic ideals. The displayed geometric-series unit
also independently confirms preservation of those ideals. Taking closures
then proves the filtration and affine-source-family assertions without
assigning a finite depth to every element of the completion.

The logarithmic derivation raises every integer `tau` valuation by at
least one. This separately defines its exponential on `K[s]((tau))`;
it also extends to the Laurent field `K(s)((tau))` when needed. The two
constructions agree on finite-depth polynomial inputs because both use
the same iterates and their estimates apply there. The sequence
`D^N x^(2N)` tends to zero in the first topology but has logarithmic
valuation `-N`. It decisively rules out a claimed continuous embedding
of the entire `D`-adic completion into that Laurent completion. The report
correctly keeps the two completed rings separate.

The added comparison to the literal source completion is also valid.
Since `D=t³(1+x²t)²`, every `D^N B_2` maps into `t^(3N)K[x,t]`.
This supplies a continuous map from the `D`-adic completion to
`K[x][[t]]`, agreeing with the exponential on polynomial inputs by
termwise continuity. The report does not need, or claim, its injectivity.

The regular two-form input was checked against
[THM-3973 / exact-volume-simple-cubic-determinantal-affine-plane-completion](../../01-canon/theorems/THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion.md):
the form extends regularly to the smooth source surface, with a zero
along its boundary rather than a pole. The function `D=p³-y²` is not
that boundary divisor. The constructed derivation is regular, and
`i_delta omega=dS` on the dense source chart implies the global identity
because the differential module is locally free. Thus its Lie derivative
vanishes. Differentiation loses at most one `D`-adic order, so the uniform
function estimates justify taking the completed pullback and substituting
an arbitrary scalar time. The literal logarithmic density is `1/tau`;
it cannot be replaced by a unit ordinary coordinate Jacobian. The local
identities and smoothness descend to every characteristic-zero base field.

Finally, non-local-nilpotence only excludes a polynomial additive-group
action with the given infinitesimal generator. It does not by itself
exclude a rational or polynomial value at an isolated nonzero time. The
report retains that qualification. The subsequent
[nonrationality audit](planar_jc_long_20260906_nonrational_audit.md) supplies
a separate generic-curve argument for its explicitly smaller family.

This is an analytic all-`R` audit. It does not claim that the finite
controls for the smaller invariant family enumerate all polynomial `R`.
