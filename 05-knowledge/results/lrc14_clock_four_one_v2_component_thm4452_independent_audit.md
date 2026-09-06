# Clean-room referee audit: q=4 one-v2 strict components

**Verdict: verified**, assuming—as in the LRC speed setting—that `r,a,b` are
positive odd integers and `a != b`.  No flaw or endpoint correction was found.
The audit is self-contained and does not import the author's implementation.

## Pointwise decomposition

An odd tail kills at most one of the four points `x+j/4`, since multiplication
by an odd number permutes the four quarter phases and a danger arc has phase
width `1/7<1/4`.  The even tail `2r` kills sheets `0,2` together when
`x in D_(2r)`, sheets `1,3` together when `x in D_(2r)-1/4`, or no sheet.
The two alternatives are disjoint.  The remaining pair must be killed by the
two odd tails, giving pointwise—including at strict endpoints—

```text
P = [D_(2r) intersect (P_ab-1/4)]
    union [(D_(2r)-1/4) intersect P_ab].
```

Here `P_ab` is the two-sheet failure set for `a,b`.

## Universal width bound

Expanding `P_ab` gives the disjoint open sets

```text
(D_a intersect tau(D_b)) union (D_b intersect tau(D_a)).
```

Every component therefore has width at most
`min(1/(7a),1/(7b))=1/(7 max(a,b))`.  Components of `D_(2r)` have width
`1/(14r)`, and its quarter translate is disjoint from it.  Applying the exact
decomposition componentwise proves

```text
lambda(P) <= min(1/(14r),1/(7 max(a,b))).
```

This reduces possible all-odd equality or violation at `1/98` to 84 rows:
`r in {1,3,5,7}` and unordered distinct `a,b` in
`{1,3,5,7,9,11,13}`.  The odd-3-unit target `1/110` reduces to 30 rows:
`r in {1,5,7}` and unordered distinct `a,b` in `{1,5,7,11,13}`.

## Exact strict results

The original four-sheet predicate was evaluated at every rational defining
wall and every complementary cell.  It gives:

- all odd: maximum `1/98`, attained exactly by
  `(r,{a,b})=(1,{7,9}),(3,{7,13}),(5,{3,7})`; runner-up `1/105`;
- odd 3-unit: maximum `1/110`, attained exactly by `(5,{1,11})`;
  runner-up `17/2002`.

Each equality row has eight physical components of the equality width.  Direct
quotient recomputation gives two components of four times that width: `2/49`
or `2/55`.  Generally `P` is invariant under translation by `1/4`; the bound
above makes every physical component shorter than `1/4`, so `x -> 4x`
identifies each four-component orbit and multiplies width exactly by four.

## Reproduction

```powershell
python -B 04-computation/lrc14_clock_four_one_v2_component_thm4452_independent.py
python -O -B 04-computation/lrc14_clock_four_one_v2_component_thm4452_independent.py
```

The runs are line-identical and end in `PASS`.
