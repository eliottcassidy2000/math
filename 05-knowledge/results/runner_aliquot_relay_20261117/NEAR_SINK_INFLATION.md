# Which inflatable near-sink cores can be exact selectors?

This follows [the reduction to two rooted equations](NEAR_SINK_PEEL.md). Write `h=H(B)` and `R=P₂(B)−2h` for a nonempty module B, so R counts ordered covers by **two nonempty** directed paths. Splitting each Hamiltonian path at each internal gap yields distinct covers and proves `R≥(|B|−1)h`. Every tournament of count one is transitive, so for h=1 we have `R=2^|B|−2`.

## One inflatable vertex is not enough

Replace an unmarked indegree-one vertex of a fixed template by B, copying the exterior contacts uniformly. Each Hamiltonian path or rooted ordered two-cover enters B at most twice: once through the sole exterior predecessor, and once at a free path start. (The first path starts at the root outside B.) Collapse each maximal B-run to a symbol. Every collapsed pattern has either one run, contributing h, or two, contributing R. Consequently H and E+O are fixed linear forms in h,R.

For E alone there is no need to track the parity of internal runs. The [marked endpoint cancellation](NEAR_SINK_PEEL.md#endpoint-cancellation) gives 2E=E+O−B_end, where B_end counts Hamiltonian paths ending at the root. Those paths have the same at-most-two-runs bound. Thus the one-sink equation or twice the direct equation becomes

    a h + b R = −d,       d=28 or 56 respectively.

If b is nonzero, `|b|R/h≤|a|+d`, so |B| is bounded independently of the target. If b vanishes, h is an odd divisor of d; the only possibilities are 1 and 7, and 7 is impossible. Then B is transitive and H is affine in a power of two, possibly constant. A finite bank of such templates therefore gives at most O(log X) exact counts up to X across these two inverse branches, even if the target chooses B.

## Two inflatable vertices: a mixed-term test

Two unmarked indegree-one vertices in a strong tournament form a chain x→y: the predecessor of y is x and the predecessor of x is exterior. There cannot be three such vertices in a strong tournament of order greater than three, since their three internal arcs exhaust their incoming arcs and they would beat the entire exterior.

Replace x and y simultaneously by X and Y. Even in this coupled situation a Hamiltonian path or rooted two-cover has at most two visits to each slot. X has at most one arrival from outside plus one unrooted start; Y can follow X at most once per visit to X. When the unrooted start is in Y instead, X has at most one visit. Collapsing both slots consequently gives a bilinear count in the pairs (h_X,R_X),(h_Y,R_Y). This applies to the one-sink defect and, by endpoint cancellation, to twice the direct defect.

Normalize by h_X h_Y and put u=R_X/h_X, v=R_Y/h_Y. A defect equation becomes

    a + b u + c v + e u v = −d/(h_X h_Y).

If e≠0, at least one of u,v is bounded: for both sufficiently large `|e|uv` dominates `|a|+|b|u+|c|v+d`. One module therefore has bounded order. Fix it and apply the one-variable argument (bilinearity is sufficient; its frozen replacement need not have indegree one). The same bounded-one-slot conclusion holds if e=0 and b,c are not of opposite sign, provided one of b,c is nonzero. If b=c=e=0, both modules must have count one and hence be transitive; the count is a nonnegative combination of products of their dyadic split counts and has at most O(log² X) distinct values up to X. Thus with these two slots, **a nontrivial candidate for escaping a finite or polylogarithmic set of exact returns must erase the mixed defect and leave opposite-signed single-slot defects**. This is a necessary test, not a sufficient construction.

### Neither branch can pass the test

In fact the coefficient of h_X R_Y equals that of R_X R_Y, separately in H, E+O, and B_end. Here is the pattern bijection. A pattern with two Y-runs and one X-run must have one Y-run at the unrooted path start and its other Y-run immediately after X. Prepend a new X-run to that unrooted path. This is legal since X→Y. Conversely, if both slots have two runs, the exterior predecessor and the free start must both have been used to enter X; each Y-run must immediately follow an X-run. Remove the initial X-run of the unrooted path. These maps are inverse and count **patterns**, not literal internal path covers.

For B_end, the same move changes the start but leaves the marked endpoint fixed. Consequently e=c for both E+O−3H and 2E−6H. Either c≠0 and the mixed term bounds one slot, or c=e=0 and there is only one normalized defect variable left. In particular no fixed strong template with two such slots can supply more than O(log² X) distinct exact returns below X, across both branches. This conclusion does not require a census. Endpoint cancellation is essential here: adding a second X-run could change the parity of the portion of X assigned to the rooted first path, so the pattern bijection need not preserve E term by term.

## Both new examples are rigid

In the p=31 order-seven core from the preceding note the eligible vertices are 2 and 5, with 2→5. For its one-sink equation the simultaneous defect is

    E+O−3H = −28 h₂h₅ −27 R₂h₅ −13 h₂R₅ −13 R₂R₅.

In particular changing either slot, or both together, cannot keep even one exact return beyond the original: equality at −28 forces h₂=h₅=1 and R₂=R₅=0, hence both modules are singletons. The [replay](NEAR_SINK_INFLATION_CHECK.py) verifies the four evaluations that determine this bilinear identity, taking each slot as a singleton or a transitive pair (both h=1, with R=0 or 2).

The p=47 core has just one eligible vertex, 1. Its defect is −28h+9R, so equality forces `9R=28(h−1)`. This yields h≡1 mod 9 and, by the split bound, |B|≤4. For orders up to four the maximum Hamiltonian path count is 5 (also checked by the replay), hence h=1, R=0, and B is again a singleton.

The obstruction is not a lack of addresses within the successful examples. It is the sign and interaction of their exact two-cover defects. The prospective next construction is an interface with different external contacts, so that the coefficients are not forced to coincide, or one permitting a genuinely independent third run.
