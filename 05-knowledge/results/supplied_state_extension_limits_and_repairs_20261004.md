# Completing supplied Collatz states and repairing the extension method

2026-10-04. **FINITE-EXACT:** all 239 frozen unresolved requests now have
checked first-hit certificates, starting from ROOT1 alone. **PROVED:** fair
conditional completion, an all-length obstruction to the fixed half-child
decoder, and a child-suffix operation that improves learned families.
**OPEN:** completion for every positive odd integer. The finite result and
the universal assertion have different quantifiers.

## What now succeeds

The [fair extension program](../../04-computation/experiments/fair_frontier_extension_20261004.py)
retains each original request, its actual forward cursor, exact guarded
word arcs, auxiliary jobs, and the component grounded at ROOT1. A smaller
child is work to do, never an automatic proof seed. Checked actual words
transport home certificates in either direction: prepend the word, or cut
it from an existing first-hit certificate. An unsupported cycle supplies
no certificate.

Every round reserves forward progress for an original request before
performing a finite auxiliary action. If a supplied source has first-hit
odd rank T and the batch has N requests, it is grounded within NT rounds,
even if some other request or symbolic search never finishes. This proves
that optional symbolic work cannot prevent completion of a finite route;
it does not prove that all routes are finite. Saved states resume exactly.

The inherited 239-source list is frozen by its hash in the
[full report](fair_frontier_extension_20261004.md). No nonroot certificate
is supplied to this experiment.

| Policy | All 239 completed at round | Fresh literal queries | Arithmetic guard tests | Post-search arc letters replayed |
|---|---:|---:|---:|---:|
| Literal only | 4101 | 3482 | 0 | 3482 |
| Guarded extension | 3708 | 3467 | 2588 | 3637 |

Both policies construct all 239 canonical first-hit certificates, checked
against an independent literal encoder. The guarded policy uses fifteen
fewer fresh literal observations and additional symbolic and validation
work. This is not a runtime speedup claim. A separate 16-source family-lift
control also completes, with 551 fresh literal queries under either policy.
An auxiliary job that yields forever cannot prevent 27 and 703 completing.
Guard checks themselves can compute valuations and endpoints; the fifteen
saved fallback observations are not fifteen fewer total arithmetic evaluations.

The [certificate bundle](fair_frontier_extension_20261004.json) retains
the 239 individual first-hit valuation words, their source identities and
the universe hash. Its words are outputs of the grounded search, not inputs
to it. The reproduction command is in the full report.

## Why increasing the old decoder cannot solve every request

For U(n)=oddpart(3n+1), a valuation word w has affine map

\[
 F_w(X)=\frac{3^{r_w}X+B_w}{2^{A_w}}.
\]

For positive odd n>=3 with n=3 mod4, so its half-child is positive odd,
the particular identity

\[
 F_w(X)=F_v((X-1)/2)
\]

forces equal word lengths and halving-cost difference one. For supplied
home routes define I(n)=A(n)-2r(n). Even allowing formal root loops, this
identity can exist at a source n precisely when

\[
 I(n)-I((n-1)/2)=1.
\]

This is the [proved all-length criterion](half_child_extension_obstruction_20261004.md),
not an extrapolation from a bounded word search. At 7 and its half-child 3,
both charges equal one. Consequently no longer word enumeration can make
that exact representation work at 7. This says nothing against 7 reaching
the root; its complete route has five odd steps.

There is also a finite source-aware rejection test. If the supplied child's
first-hit rank is R, at most R source steps decide the strict half-child
ansatz. Preserve those observed steps and their frontier when rejecting
the ansatz. Never reject the original request along with the failed model.

Changing the child map immediately recovers a guarded infinite family:

\[
 7+256t\xrightarrow{(1,1,2,3)}5+162t
       \xleftarrow{(1)}3+108t,\qquad t\ge0.
\]

Its smaller child is (27n+3)/64. It requires a certificate for that child
at the supplied parameter; the seed certificate alone does not prove the
whole family converges.

## Reusing a child suffix can enlarge a learned family

Suppose checked words give n to J from one side and m to J from the other,
with m<n. Let their lengths be r,s and halving costs A,D. The lift slope is
lambda=2^(D-A)3^(r-s). Cut p edges with total halving cost C from the
supplied child's prefix, leaving a new child x<n. The
[child-cut theorem](supplied_state_join_extensions_20261004.md) gives

\[
 \lambda'=\lambda\frac{3^p}{2^C},\qquad
 \frac{P}{P'}=3^{\min(p,\max(s-r,0))},
\]

where P and P' are the old and new least source periods. If x<m, positive
carry makes the slope strictly smaller. The new source progression contains
the old one, strictly when P'<P. All-height decrease still requires lambda'<=1; merely
finding a smaller x is insufficient.

The expanding pair 233/231 illustrates the difference. Cutting the child's
route at 31 changes slope 2^27/3^17>1 into 8192/59049<1 and yields

\[
 233+118098t\ \longleftarrow\ 31+16384t.
\]

The new source progression contains the old one and is 2187 times denser
within their shared arithmetic setting. In contrast, simply moving the
same pair's common endpoint later leaves its slope unchanged. The 75/73
example, cut at 71, is the hostile control: the child decreases but the
slope remains above one.

## Research boundary and next operation

Inheritance: grounded observation union is the closest proved mechanism;
7 and the expanding 233/231 pair are the decisive hostiles; a successful
short decoder was the near miss; retained first-hit rank, halving cost and
child provenance are the necessary sidecars. The incoming
[inverse shields](collatz_join_shields_and_lifts_20261004.md) explain why
deeper inverse search at a fixed early endpoint can saturate. The
[binary and ternary bank](collatz_binary_ternary_guard_fusion_20261004.md)
has a positive-density necessary remainder, so its present one-stage
coverage cannot be promoted to a universal selector.

The live board is original source / fixed versus variable child / exact
ordered word / rooted closure / first-hit rank and cost / fair scheduling.
The anchor is completion of supplied states; the niche is a finite
decision for a failed representation; the wildcard is broadening learned
families by retaining alternate child suffixes.

| Connection | Preserved predicate | Information that must remain attached | Cheapest decisive test |
|---|---|---|---|
| Actual word to undirected proof edge | Reaches ROOT1 iff its endpoint does | Directed word, exact source guard, grounded provenance | Disconnected component must remain uncertified |
| Completed route to charge I | Invariant under root padding | Separate first-hit rank and halving count | 7 versus half-child 3 |
| Child prefix cut to broader family | Exact common future | Child identity, suffix certificate, periods and slope | 233/231 repair and 75/73 partial repair |
| Fair search to finite-batch completion | Each finite route eventually receives enough work | Original cursor and finite action boundary | An endlessly yielding decoder beside 27 and 703 |

The next policy experiment should prioritize guarded arcs by which unresolved
request components they can connect, keeping the same guaranteed original
turns and separate cost counters. A bounded or shielded decoder should
retire without discarding its observed path. The sharper mathematical
question is whether an expanding join must admit a useful cut on its
particular child's prefix; 2914 bounded examples support it, but no proof
is supplied. Of those examples, 2822 already have an endpoint below the
source; only 92 need an earlier child cut. Even such a theorem would not
prove that every source finds a grounded join.

All three component packages have explicit hostile controls and identical
normal and optimized Python outputs. Independent audits check the
half-child proof, child-cut formulas, rooted closure and certificate
construction. The full proof obligation that remains is universal
grounded coverage, not a larger successful finite census.
