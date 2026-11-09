# Bounded feedback-arc spectra are polylogarithmically sparse

Status: PROVED and independently audited ([audit](../ribbon_audit/AUDIT.md)). Chan–Pak remains OPEN.

## Inheritance and scope
Read current AGENTS, protocol and guardrails, THM-4748 (two reversed arcs), THM-4940 (packing counts in every base), MISTAKE-739, and Session229/230 syntheses. The strongest hostile control is the three-back-arc realization of every Mersenne count: an infinite arithmetic family does not provide a dense spectrum. THM-4878 section 6 already proves the qualitative bounded-back-arc obstruction, with its detailed counting proof stated for connected cutwidth-2 words. The present contribution is an explicit arbitrary-cutwidth vertex kernel and description bound, not the qualitative sparsity conclusion. Its fixed-k exponent is 2k rather than the earlier gap-count exponent 2k+1. The moving-k bound here is uniform without any implicit k-dependent constant; it is weaker than the inherited fixed-cutwidth result, as expected because arbitrary chord patterns are allowed.

## Theorem
Let X>=3 and k>=1 be integers. Put L=floor(log_2(X-1)) and N=k(L+k+1). The number of distinct Hamilton-path counts at most X realized by tournaments with feedback arc number at most k is at most

    sum_(j=0)^k binom(N(N-1)/2,j).

In particular this count is O_k((log X)^(2k)). No fixed feedback-arc budget realizes all admissible odd counts, or even a positive-density subset of the integers.

## Proof
Choose a vertex order with r<=k backward arcs. Fix a backward arc y->x, with x<y and span s=y-x. For each other backward arc delete one endpoint other than x,y. This is possible because no distinct arc has endpoint set {x,y}. At most r-1 vertices are deleted. Restrict to the original interval [x,y]. In the induced remaining tournament the only backward arc is y->x, and at least s-r interior vertices survive. If s>=r+1, the classical single-back-arc formula gives H>=1+2^(s-r). The odd-cycle packing formula gives vertex-deletion monotonicity. Thus H(T)<=X implies s<=L+r. If s<=r the same upper bound is automatic.

Delete every vertex outside the union of the closed backward-arc intervals. Such a vertex is a singleton strong component: no backward arc crosses it, so every earlier vertex beats it and it beats every later vertex. Removing it does not change H. At most r(L+r+1)<=N vertices survive, since every interval has at most L+r+1 vertices. Append sinks until there are N vertices, preserving H and the <=k back-arc description. Each resulting tournament is specified by at most k reversed pairs out of N(N-1)/2 pairs. Counting descriptions, allowing collisions, proves the bound.

The case r=0 has H=1 and can be represented by the transitive N-tournament.

## Quantitative moving-budget consequence
For each fixed 0<epsilon<1/3, let

    k(X)=floor((1/3-epsilon) log X / log log X).

The count represented with at most k(X) back arcs is at most X^(1-3epsilon+o(1))=o(X). Indeed N=O(k log X), M=N(N-1)/2=O(k^2 log^2 X), and sum_(j<=k) binom(M,j)<= (k+1)(eM/k)^k. Its logarithm is at most k(3 log log X+O(1)). Thus for density one of target integers, any realization (if one exists) requires more than k(X) back arcs. This is a target-complexity obstruction, not an assertion that the targets have realizations.

## Implication for the constructive program
A fixed finite menu of sparse reversed-arc templates cannot close Chan–Pak if the total number of reversals stays bounded. Growing-budget constructions, fixed cutwidth with unboundedly many backward arcs, apex words with unrestricted losses, and proper substitution are not ruled out. The next positive target is a bounded-interface reset that admits unboundedly many arcs while implementing two arithmetic digit choices; fixed back-arc constructions are the wrong invariant to preserve.

## Base-b version
For b>=2 and X>=b+1, the same proof gives an order kernel

    N_b = k(k+2+floor(log_2((X-1)/b))).

A single backward arc with t>=1 surviving interior vertices has exactly 2^(t-1) odd directed cycles, all through that arc, hence H_b=1+b*2^(t-1). The target arc of span s leaves t>=s-r interior vertices after removing endpoints of the other arcs. Thus s<=r+1+floor(log_2((X-1)/b)), with the small-span cases automatic. Outside all backward intervals, singleton strong components contribute factor one in every base. The same binomial-description bound holds using N_b. This is not a spectrum-completeness claim in any base.
