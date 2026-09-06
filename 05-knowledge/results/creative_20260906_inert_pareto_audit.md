# Independent audit of the inert-pair five-profile envelope

**Verdict: PASS.** The all-height envelope, unique scalar maximizers,
inclusive mass/component disjunction, and improved original-body gates in
[creative_20260906_inert_pareto.md](creative_20260906_inert_pareto.md) are
accepted without mathematical repair. Their scope is sufficient entry under
the stated arithmetic and body hypotheses. No universal body lower bound
or LRC(14) closure follows.

## Analytic audit

I read the complete candidate proof and source, plus the relevant current
proved supplier arguments in:

* **THM-4153**, `01-canon/theorems/THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer.md`:
  exact Bernoulli overlap formula, `1/(2pq)` tail, physical two-half-lift
  carrier, compact-to-open mass strictness, and literal width containment;
* **THM-4453**, `01-canon/theorems/THM-4453-lrc14-inert-sum-five-ray-disjointness-and-dyadic-entry-closure.md`:
  arithmetic filter, actual labelled decoder-component scope, and the
  original absorption gates;
* **THM-4450**, `01-canon/theorems/THM-4450-lrc14-absorbed-label-overlap-hierarchy-and-component-address-decoder.md`:
  the ten-distinct-label absorption hierarchy and the doubled-body half-mass
  identity.

The exact gap is `20/469-2/49=6/3283`. Hence `pq>=274` is strictly below
the claimed mass maximum. The complete 21-element smaller-product universe
has unique maximizer `(1,67)`. The all-height tooth bound `b<=2/(7q)` makes
`q>=15` strictly smaller than `1/49`; the unique admissible smaller-height
pair is `(7,13)`, whose two literal open intervals each have width `1/49`.
These are valid finite reductions with explicit infinite-tail proofs.

Every `q>67` has both profile coordinates strictly smaller than those of
`(1,67)`. The remaining 46 profiles have exactly the claimed five maximal
points, all attained and mutually incomparable. The five ratios also lie in
the bounded actual ratio atlas. That latter fact does not supply actual
entry for an arbitrary proposed component split.

The carrier is a proper open subset of the circle, and neither zero nor one
belongs to it. Within the fixed chart, opposite-parity tooth intersections
retain all strict endpoints and are separate components; the source does
not merge across excluded walls. For an odd common scale `g`, the true
failure carrier is its inverse image under multiplication by `g`, with the
same mass and component widths divided by `g`. Thus no clock factor is lost.

For a nonempty compact safe set contained in this proper open carrier,
the mass inequality is strict: equal measures would make the open set
difference empty, forcing a nontrivial clopen set. Each positive-length
closed body component is likewise strictly shorter than its containing
open component. Therefore the source's `>=` thresholds, including equality,
are sound. Choosing a table profile dominating the actual profile makes
the rowwise mass-or-width disjunction sufficient. Zero-mass or empty safe
sets cannot pass its positive thresholds.

The absorption arithmetic is exact:

    2*(20/469)=40/469 < 8/91,
    20/469+8/63=716/4221 < 20/117.

Here a ten-body means ten distinct positive labels. THM-4450 states its
clock-two overlap argument with the absorbed label outside the body; if it
is already present, the claimed lower bound is immediate because adding
it changes no safe set. No extra exception or approval is needed. The
half-mass estimate for the clock-four form follows directly from the two
physical lifts. These facts justify both inclusive improved gates.

## Independent exact computation

The [independent source](../../04-computation/creative_20260906_inert_pareto_audit.py)
imports no candidate program. It generates admissible sums by a prime sieve
and bounded products of inert prime powers, rather than factoring each
candidate sum. Its geometric routine constructs the walls directly from
all equations `||s(y+epsilon)/2||=1/14`, with `s=p,q` and `epsilon=0,1`.
On every rational cell it evaluates the literal two-half-lift Boolean
predicate. It merges adjacent cells only when the intervening wall itself
belongs to the carrier, and separately verifies that every final component
endpoint is excluded. This is independent of the producer's parity-labelled
two-pointer tooth intersection routine.

The independent engine checks the complete 21-pair product head, complete
46-pair `q<=67` head, all 582 primitive odd pairs with `q<=75`, and all 548
odd-3-unit ratios of the inherited `p+q<=356` atlas. The 582-pair bank compares
literal physical geometry against the Bernoulli formula, including excluded
arithmetic families. All three hypothesis hostiles, the empty carrier,
exact width-extremal intervals, joint-gate example, and absorption constants
are independently reproduced. The five-profile frontier is recomputed from
all 46 profiles, not supplied as a pruning assumption.

    python3 -B 04-computation/creative_20260906_inert_pareto_audit.py
    python3 -B -O 04-computation/creative_20260906_inert_pareto_audit.py

Both runs pass **9,505 always-active gates**, with byte-identical LF output.
The matching [frozen output](creative_20260906_inert_pareto_audit.out)
is retained. Source SHA256:
`c45772a01f2b36f86252231854747b87f4848cc466a225c61449bed0ea4d650c`.
Output SHA256:
`214bf0c73c508a557299764860004c392302b295e6df8f08e4433e32d212156d`.
