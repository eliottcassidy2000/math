# Independent root review of the direction and walk transfers

**Status: independent analytic audit ACCEPTED; exact replays documented below.**
Date 2026-09-06. This review is separate from the two producing agents.

## Primitive-direction Hermite precision

The [direction theorem](creative_20260906_smith_bridge.md) was read in full,
including its inherited [THM-4439 / terminal-cluster input](../../01-canon/theorems/THM-4439-all-node-twojet-metric-precision-by-terminal-clusters.md).
The two cardinal forms do have the specified value and derivative; their
degree is `2n-1` and the other nodes are killed to order two. The coefficient
content argument is exact: every bracket factor is primitive, the local
linear coordinates are unimodular, and `gcd(d_i,c_i)` divides `d_i`.
Consequently the proposed denominators are attained, and their lcm is the
largest Smith factor, rather than merely an upper bound on it.

I independently checked the tangent-shift and determinant-minus-one actions,
the cancellation of powers in the affine determinant, and the residue-class
multiplication map. The last map is unimodular because its determinant has
valuation zero after comparing the full and block determinant formulas.
This avoids requiring a nonexistent single integral affine chart when all
`p+1` projective residue classes occur. Integer approximation to precision
greater than the determinant valuation transfers every determinantal divisor
to local integer nodes, justifying the stated use of THM-4439. Singleton,
infinity, depth-zero and primitive-content boundaries are retained. The
general full-partition metric claim remains false, as required.

A separate [SymPy verifier](../../04-computation/creative_20260906_smith_independent.py)
imports no producer code. It constructs the jets by literal symbolic
substitution and coefficient extraction, uses a full integer Smith algorithm,
and compares both the largest factor and exact matrix inverse denominators
with the proposed formula. It also verifies the bracket determinant.
All **60 configurations** pass, with one through six nodes, all projective
residue classes at primes two and three, signed directions, infinity and
deep infinity clusters. The [matching output](creative_20260906_smith_independent.out)
is frozen. The computation uses installed SymPy 1.9; the producer's verifier
uses only the Python standard library, giving a distinct arithmetic path.

    python3 -B 04-computation/creative_20260906_smith_independent.py

No external library's success is treated as a proof of an unstated
mathematical consequence. The universal result is the audited argument.

## Divisor-component walk compiler

The [walk theorem](creative_20260906_lrc_bridge.md) and its inherited
[exact cancellation identity](overnight14_20260906_lrc_endpoint_walk.md)
were read in full. For any visited set of gcd `h`, its induced divisible
component has gcd both divisible by and dividing `h`; this proves the
necessary equality. A spanning-tree walk proves the converse and attains
the maximum support. The endpoint path is used once and other tree edges
twice, giving `2s-2-ell<=2n-3` with distinct endpoints. Enumerating gcds
over subsets of the `n-2` other vertices is complete without factoring.

The budget version correctly omits exact gcd equality: a component's actual
factor divides the nominated budget divisor. The LRC application keeps the
physical core scale, maximum endpoint, actual graph, finite box and equality
entry. Visiting *at least* r vertices permits the inherited r-subset cap
because the full visited gcd divides every r-subset gcd. Hence its final
inequality really supplies the old endpoint-gcd criterion. No new safe
class or universal existence of an eligible walk is claimed.

The producer's exact three-way comparison covers all 100,544 declared
graph/label/endpoint cases with **559,903 always-active gates**. A root
replay is compared byte for byte to the frozen output. The two hostiles
respect their different scopes: `6--2--3` is an actual decoder path, while
the repeated-walk star is a general labelled-graph example only.

    python3 -B 04-computation/creative_20260906_lrc_bridge.py

These two results are accepted in their stated scopes. They do not turn
the coefficient lattice or a decoder walk into a physical LRC solution
without the explicit consumer hypotheses.
