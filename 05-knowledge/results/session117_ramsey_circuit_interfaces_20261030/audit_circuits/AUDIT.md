# Independent circuit and local-potential audit

Session 117. Scope: `../circuits/POINT_DEFECT_CIRCUITS.md` and the local-weight
claim in `../challenge_audit/AUDIT.md`. No novelty or unrestricted lower-bound
record is asserted. The full basis throughout is all 16 binary Boolean gates,
with unrestricted sharing/fanout and a single output. This is not a U2 audit.

## Verdict: exact point-defect classification passes

For n>=3 and k=wt(r), the claimed costs n-1 (k<=1) and n+k-2 (k>=2)
follow from the displayed construction, essential-variable bound, branching
lemma and independently reproduced three-variable finite base.

Detailed checks:

1. Input translations can be absorbed separately into every receiving gate's
   truth table, even with shared fanout. Output complementation is absorbed
   into the output gate. Since these functions have at least two essential
   inputs, an output gate exists. No free NOT gate is being assumed.
2. In a minimum circuit, unary and constant gates can be eliminated by
   substituting them into their successors. Repeated input pins at a gate
   likewise make that gate unary. Gates outside the output ancestry are
   unnecessary. All inputs of the target family are essential.
3. If every input has exactly one receiving gate, the earliest gate uses two
   distinct inputs, and the rest of the DAG only sees those two inputs through
   its one output bit. Shared descendants do not invalidate this factorization.
   The parity and parity-plus-product pair cofactors have three joint responses,
   even after either is independently complemented. They cannot both factor
   through a single bit. This proves a branching input for pure-parity defects.
4. Setting a branching input to one makes every one of its distinct receivers
   unary or constant. These are distinct original gates, so removing two is
   legitimate even if one receiver lies downstream of the other. Subsequent
   simplification cannot restore a removed gate. With at least two essential
   residual inputs, a terminal negation can be incorporated into a remaining
   gate; there is no one-gate normalization penalty.
5. Therefore C(F_m)>=C(F_(m-1))+2 for m>=4. The three-variable base, not the
   exceptional two-variable OR case, correctly starts the induction. Removing
   inactive coordinates one at a time costs at least one gate; for k=2 the
   process correctly stops at three variables, rather than claiming the false
   two-variable base.

## Independent computational certificate

`independent_checks.py` imports neither producer implementation. Its raw DAG
search uses explicit pointwise evaluation of all 16 binary truth tables, both
ordered pins, repeated pins, and every previous signal. It retains ordered gate
histories without complement quotienting. Every circuit with at most three
non-input gates has such a topological description; therefore the search covers
sharing, not merely formulas. A circuit whose output is an earlier gate is
already covered at a smaller depth.

Results:

- 38 one-gate histories and 1,744 two-gate histories;
- 232 distinct raw three-input truth tables reachable through three gates;
- costs 2,2,3,4 for canonical support weights 0,1,2,3;
- the same costs checked for every r, every defect v and both output signs;
- actual wire-indexed gate lists (not just a scalar identity) match the formula
  and gate count on 20,464 evaluations for n=3,...,10;
- all four additive-constant variants of the three-response obstruction pass.

Absence through three gates supplies the lower bound four, while the separately
validated construction supplies the upper bound four. This is a finite exhaustive
certificate, not a proof-assistant or independently checked SAT proof log.

## Verdict: potential bounds pass, with strictly local scope

The archived source's Section 7 row table implies, after g=1 normalization,

    i <= 2+p,
    t <= 2+2i-p,
    t <= 3+2p,
    t <= 3+i+p.

For i,p>=0 these are exactly the nontrivial constraints: the other zero-loss
rows are automatically nonincreasing, and GP with k>1 is weaker than k=1.
For c=t-i, the first two give c<=2+i-p<=4. AP2 and GP give
p>=(t-3)/2 and i>=(t+p-2)/2>=(3t-7)/4, hence t>=4c-7.
These deductions do not require any claim that the rows describe every real
circuit state.

The proposed segment i=3c-7, p=2c-5, t=4c-7 for 5/2<=c<=4 is feasible for
all rows, with GP and AP2 binding. P2 gives c<=4; AP1 gives c>=5/2.
Thus coefficient four forces penalty at least nine, attained at (i,p,t)=(5,3,9).
For a circuit telescope the usual initial/final-state invariants and t>=0 also
matter; all weights on this reported segment have t>=3. Establishing local
algebra cannot establish those invariants or the claimed universal circuit
elimination theorem.

Independent tests check 251 exact rational points against every row type
(including GP k=1,...,9), plus 619 feasible integer triples in the tested box.
The analytic monotonicity of GP in k covers larger k. These tests corroborate,
not replace, the exact inequalities above. The audit does not certify the
provisional four-coefficient challenge submission, its disperser construction,
or its full-family assembly.

## Reproduction

    python3 session117_research/audit_circuits/independent_checks.py
    python3 -O session117_research/audit_circuits/independent_checks.py

Both executions produced byte-identical JSON, archived as
`independent_checks.json`. Failure checks use explicit exceptions and survive
optimization. No producer code, publication metadata or ledger was edited.

## Addendum: complete rowwise frontier independently reviewed

The subsequent `../challenge_audit/POTENTIAL_FRONTIER.md` extends the segment
above to the complete nonnegative-coefficient classification. Its algebra passes
independent review. Setting i=b-c and t=b transforms the constraints exactly into

    b >= c,
    max(0, b-c-2, (b-3)/2, c-3) <= p <= b-2c+2.

Comparing each lower endpoint with the upper endpoint gives, respectively,
b>=2c-2, c<=4, b>=4c-7 and b>=3c-5. The last is redundant: when c<=3,
2c-2>=3c-5; when c>=3, 4c-7>=3c-5. Together with the specified c>=0, this proves
both necessity and sufficiency of

    0 <= c <= 4,  b >= max(c, 2c-2, 4c-7).

The first boundary segment (0<=c<=2) uses (i,p,t)=(0,0,c), and its binding
restriction is t<=2. The second (2<=c<=5/2) uses
(i,p,t)=(c-2,0,2c-2); GP binds and AP2 supplies c<=5/2. The third is already
checked above. Adjacent witnesses agree at both joins. At c=4,b=9 both interval
endpoints are 3, proving the claimed uniqueness of normalized optimal weights.

This is a complete classification of the stated local row-certificate cone for
c>=0, not of all potentials, all circuit arguments, or actual circuit complexity.
The candidate's invariant, coverage, termination and explicit-family claims
remain outside this independent algebra certification.
