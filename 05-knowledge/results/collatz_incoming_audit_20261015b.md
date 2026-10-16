# Incoming two-octave work: an assumption repair and a joint-branch hostile

2026-10-15. **PROVED:** the repaired density implication, the explicit
mean-zero branch-correlation hostile, and the annular source-row formula.
**FINITE-EXACT:** 10,923 two-octave identities,14 integer energy-seam checks,
and the stated finite-cylinder calculation. **OPEN:** any all-scale
Collatz energy budget. This is a scoped audit of incoming commit
`378a9b6390`, not an audit of its full 2^34 census or Gaussian surrogate.

[Script](../../04-computation/experiments/collatz_incoming_audit_20261015b.py),
[output](collatz_incoming_audit_20261015b.json),
[main Gram/glue packet](collatz_gram_slack_20261015b.md).

## 1. Repair scalar-to-Haar promotion before using it

[THM-4658, two-octave-memory-energy-descent-critical-lines](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md)
4(a) originally inferred that every residue class has its Haar share of
a leading term from the total-sum asymptotic S(N)=cN^rho+o(N^sigma).
That inference omitted the excess coordinate. The exact identity is

    C(2 mod3;M)=(S(M)+H(M))/3.

Thus the stated scalar consequence D(rho)=0 is proved if one also assumes
H(M)=o(M^sigma), or the stronger Haar-leading class asymptotics. Without
this hypothesis the class-sum recursion leaves a leading H contribution.
This is a proof-gap repair, not a counterexample to the same implication
for the actual Collatz sequence. The native mean-zero term cannot simply
be discarded. The theorem statement and proof now retain that hypothesis;
its exact Haar eigenvalue calculation and excess-distribution formula4(b)
are unaffected.

The incoming note also inferred nonexistence of any Pell-type closure
from variation of a measured cross-correlation. Variation does not exclude
a richer deterministic state. Its wording now says the proposed scalar
closure is unestablished. Our changing-generator identity explains why
that distinction is productive: exact Pell algebra can appear in a larger
state without supplying native Collatz energy control.

## 2. A mean-zero finite-cylinder hostile to branchwise contraction

The incoming operator has branches, for real rho=sigma,

    U f(y)=f(2y),
    V f(y)=1_(y=2 mod3) f((2y-1)/3),
    P* f=-a U f-b V f, a=2^(-sigma), b=(3/2)^sigma.

Their individual norm identities give
phi(sigma)=a^2+b^2/3. This is less than1 for1/2<sigma<1.
It does **not** follow that their sum contracts, even on mean-zero inputs.

Take the function constant on classes modulo9 with values, in order0..8,

    f=(-1,-2,2,-1,-2,2,-1,1,2).

Exact integration on the27 output cylinders gives

    integral f=0,
    ||f||^2=||Uf||^2=8/3,
    ||Vf||^2=8/9,   <Uf,Vf>=2/3.

Consequently

    ||P*f||^2/||f||^2=phi(sigma)+ab/2.

At sigma=3/4 this equals

    (1+sqrt(3))/(2sqrt(2))+(1/2)(3/4)^(3/4) > 1,

although the first term alone is less than1. These are exact strict
inequalities, not floating-point tests. For example the first term is
between5/6 and1, and the second exceeds3/8.

A construction explains the witness. Start with the first-level wavelet
h(y)=3[y=2 mod3]-1 and the refinement operator
W h(y)=1_(y=1 mod3) h((y-1)/3). Then f=h+Wh. The two scales are
orthogonal, but W links them under the next transfer, creating the
positive cross term. This is a small recursive gadget that detects what
branchwise squared norms forget.

The zero-mean rank-one correction in THM-4658 4(b) has dual term
proportional to integral f, so it vanishes on this witness. This does not
prove a slow eigenmode exists or refute any discrepancy budget. It rules
out obtaining a full contraction from the displayed branch square sum
alone. The joint-terminal Gram construction retains exactly such cross
terms before taking absolute values.

For complex rho, the displayed formula is the bilinear transfer dual;
the Hilbert-space adjoint conjugates its scalar coefficients. Our witness
uses real rho, so this convention does not affect the calculation.

## 3. A new source row suggested by the exact memory identity

The incoming identity is

    H(4N)=H(N)+H_(2,2)(3N)-H_(2,8)(6N).

Its increment has the particularly simple source representation

    Delta(N)=H(4N)-H(N)=sum_(N<n<=4N) w(n)t(n).           (1)

The old prefix through N cancels exactly. For each N in[X,2X), push this
annular weight row by any common authenticated finite controller, keeping
its actual source positions, and call the resulting row b_N. Then

    D_X=sum_N Delta(N)^2=z^T G_D z,
    G_D=sum_N b_N b_N^T,
    D_X<=||G_D||_1.                                      (2)

The same signed pair-fibre slack and observed-graph closure apply to G_D.
This is a precise new interface to THM-4658's 4-adic energy descent. The
row has no incoming-prefix weights because it is a *difference* of H;
that cancellation does not justify deleting the incoming current from
our original coarse J rows.

It is worth testing both (2) and the full joint rows for H(N),Delta(N).
Their cross term can be strongly negative or positive. Bounding the two
marginal energies separately may throw away the information needed to
make the next octave inexpensive. An all-scale bound for this annular
Gram matrix is not established here.

## 4. Two rigorously conditional growing-budget routes

Incoming [HYP-9283](../hypotheses/HYP-9283-collatz-energy-budgets-constant-fail-log-suffice.md)
provides an interval-arithmetic route under energy budgets
E_X<=X^2(64+16j), X=2^(34+j). Its proposed coarse target15+4j is
consistent with the fine-scale allowance. We read the tail script's
Cauchy--Schwarz sum and geometric tail bound; we have not rerun its Arb
calculation or promoted the Gaussian model to an all-scale fact.

Independently, G2 in the Gram/glue packet proves the sufficient coarse
budget k^2/100 by analytic comparison with the **already verified** B3
logarithmic pointwise tail allowance. Thus it does not depend on the new
Arb run. Neither coarse budget dominates the other at every scale. Both
retain explicit constants and thresholds, and both remain unproved
premises for actual native controllers.

The incoming surrogate suggests that constant energy budgets may fail;
this is HEURISTIC. It does not refute the constant4/9 Gram target. The
productive rerouting is to retain both certified implication routes and
focus the arithmetic search on weighted graph reentry and cross terms.

## 5. Checks

    python3 -B 04-computation/experiments/collatz_incoming_audit_20261015b.py

An independent finite first-descent table through65536 verifies the
memory identity for every0<=N<=10922 (the6N cutoff is retained) and the
exact remainder inequality on14 dyadic windows. The27-cylinder witness
is integrated using rational arithmetic. The native complementarity
control t(3)-t(8), with selectable edges3->5 and5->8, gives Gram norms
4,4,4,0 for the four subsets. It confirms the non-submodular payment
obstruction in the main packet. These are controls on the proved
identities, not a universal entry or convergence argument.
