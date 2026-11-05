# Independent review: capped moving-port automata

**PASS.** Frozen source snapshot in `../session197-orbit-wildcard/capped-ports/`:

- `PROOF.md`: `77be8657a64155c261aebb00e06b3e3db7114ad764d6b417a51c119ff85e4ad1`
- `verify.py`: `6486bfe3df04a1d6a680dd374167ab772c5a251f68c2761923958caebf2244bc`
- `receipt.json`: `4c86a2bd36ded8f0108c418ef85f60170cbe688709e99535b8adbea51728da44`.

A copied source verifier run under normal Python and `-O` exactly reproduces the frozen receipt. No source files were changed.

## All-order proof review

A consecutive run of unexposed singleton SCCs is a transitive module, including its incidences to exposed vertices elsewhere. Permanent non-exposure is the crucial invariant: every future birth has the same default incidence to all members, and forgetting other ports cannot invalidate that. The same alphabet word therefore produces the same quotient graph and corresponding ordered ports after any positive clone length is substituted.

The quotient-SCC dichotomy is exhaustive. If the quotient vertex remains a singleton SCC, the expanded clones remain individual singleton SCCs and contribute the factor one. If it enters a nontrivial SCC, routes to and from exterior vertices lift through every clone because the exterior incidences are uniform. No partial clone merge is possible, including after later port changes. This justifies exact continuation equivalence, not merely equality on the tested words.

With b=max{t:5^t<=N^3} and M=b+1, a merged compressed block belongs to a strong component of order at least b+2, so Busch gives H^3>=5^(b+1)>N^3. The threshold is strict and safe. Both original and compressed graphs then have capped output N+1. The argument only applies Busch to an actual nontrivial strong component; at M=1, such a component automatically has order at least three.

For bounded H, the product of SCC counts bounds their number by floor(log3 N), and the product of Busch inequalities bounds the sum of their orders minus their number by b. Exposed singleton SCCs separate at most k additional runs. This yields the stated V. There is no mistaken assumption that intermediate H divides the cap N.

The deterministic machine proof maintains equivalence for every suffix, rather than claiming that a compressed graph is physically reachable from the original seed. That distinction is correctly explicit. Monotonicity under induced extension makes the rejecting sink sound. The finite labelled state count includes every possible ordered port list. The edge cases N=1 and k=0 are covered. The theorem establishes regular capped languages, not finite unrestricted linear rank, unrestricted spectrum coverage, or permission to re-expose old vertices.

No mathematical or scope correction is required.

## Independent finite controls

`capped-ports-audit/independent.py` builds the complete cap-one machines using **only transitive vertex orders**. A new vertex preserves transitivity exactly when its win/loss pattern along that order is a suffix; otherwise it rejects. Compression keeps the first vertex of each unexposed order interval, then performs the specified old-label relabelling. This calculation uses no SCC algorithm, path counter, or source compression routine. It independently obtains:

- k=0: one reachable state, two transitions, maximum stored order one;
- k=1: twelve states, ninety-six transitions, maximum order three;
- k=2: 188 states, 4,512 transitions, maximum order five.

Every resulting transition is also compared with the source's graph-level machine, including the reachable sink.

A separate literal path-word DFS stops only after N+1 complete Hamilton paths have been found. It compares original long modules and shortened modules on 720 continuation prefixes of length at most four, at caps 1,3,5, with two moving ports and all letter types sampled deterministically. A separate reachability-search SCC implementation checks the singleton/full-merge dichotomy on the same 720 physical graphs. All checks pass. These are independent controls, not a substitute for the all-continuation proof.
