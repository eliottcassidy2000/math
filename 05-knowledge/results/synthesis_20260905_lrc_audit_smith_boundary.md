# Independent audit of the complete-residue two-jet Smith boundary

**Status: ACCEPT, analytical audit.** I independently checked the proof in
[synthesis_20260905_wildcard_smith_boundary.md](synthesis_20260905_wildcard_smith_boundary.md)
against the scope of
[THM-4080](../../01-canon/theorems/THM-4080-confluent-two-jet-single-scale-smith-partition.md).
The new `s=p` formula and consecutive-node corollary through `n=p^2` are
valid as stated. This audit does not claim an independent full execution
of the producer's script; root separately replayed it, and the following
proof checks address the two delicate universal steps.

## Saturated derivative rows really extend to an h-minor

For `p<h<2p`, use the original `p x h` derivative matrix `R` on degrees
`0,...,h-1`. Its rank modulo `p` is exactly `p-1`; its columns `1,...,p`
have a determinant of exact valuation one. Perform row operations over
`Z_(p)` lifting ordinary elimination modulo `p`. The resulting first
`p-1` rows are independent modulo `p`, and its last row is divisible by `p`.
Divide that row by `p` to obtain `R*`.

The exact-valuation-one minor forces `R*` to have rank `p` modulo `p`.
Otherwise every maximal minor of `R*` would be divisible by `p`, and every
maximal minor of the original `R` by `p^2`, a contradiction. Thus this
step can be performed using row operations alone; a Smith column change
does not have to be applied to the value rows.

The complete value-and-derivative matrix has rank `h` modulo `p`:
a polynomial of degree less than `h<2p` in its kernel would be divisible
by `(X^p-X)^2`. This double-root criterion holds also in characteristic
two because the first Hasse derivative equals the ordinary derivative.
The old derivative span is contained in the saturated derivative span.
Consequently the latter `p` independent rows can be extended to a basis
using `h-p` rows from the actual value-row list. Undoing the single division
by `p` and the row-unimodular operations gives an original residual
`h`-minor of exact valuation one. The selected value rows are unchanged;
the determinant includes all `p` derivative rows, so the row transform
acts by its unit determinant on this exact minor.

This validates attainment of every middle determinantal-divisor bound,
not only full determinant valuation or modular rank. The degree-cost lower
bound remains valid for `e=1`, when a nonminimal degree choice can tie the
same valuation: existence of additional tied witnesses does not affect it.

## CRT and the endpoint n=p^2

For `n` consecutive integer nodes, each residue class modulo `p` contains
a consecutive progression `b,b+p,...,b+(s-1)p`. When `n<=p^2`, one has
`s<=p`, so all nonzero differences among the rescaled nodes
`0,...,s-1` are p-adic units. Every class is therefore exactly in either
the inherited `s<p` single-scale theorem or the new `s=p,e=1` theorem.

Let the monic factor for a class be the product of `(X-x)^2` over its
nodes. Factors for different classes have unit resultants over `Z_(p)`;
their polynomial quotient modules split by CRT. The CRT map between the
monomial bases of total degree below `2n` and the block remainder bases is
p-unimodular. Evaluation of values and first derivatives factors through
these remainders because their difference vanishes to order two at every
node in that class. Hence the observer cokernel is exactly the direct sum
of class cokernels, and the p-primary Smith multiset is their union.

The boundary `n=p^2` is included: every class has size exactly `p`, with
one and only one internal p-adic scale. At `n=p^2+1`, one class has size
`p+1`; its endpoints differ by `p^2`, so the proof correctly stops.

## Independent exact lowest-prime control

For `p=2` with two nodes translated to `0,d`, `v_2(d)=e`, two unit pivots
reduce the observer by integral row/column operations to

```text
diag(I_2, [[d^2,d^3],[2d,3d^2]]).
```

The remaining first determinantal divisor has valuation
`min(2e,3e,e+1,2e)=e+1`, and the determinant has valuation `4e`.
The resulting full profile is exactly `(0,0,e+1,3e-1)`, including the
degenerate tie `(0,0,2,2)` at `e=1`. This independently confirms both
the `+1` correction and its independence from the geometric scale `e`.

No scope, saturation, quantifier, or characteristic-two error was found.
