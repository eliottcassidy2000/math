# All-order multi-star rounding stability

Status: PROVED, independently audited. The final matching-uncertainty
corollary follows directly from the same audited Pfaffian counting mechanism.

## Statement

For a weighted tournament A of order n>=3, let T be its independent binary
rounding, with each edge oriented i->j with probability A_ij. Write
Delta(A)=H(A)-det(I+A)-|det(A)| and V=sum_{i<j}A_ij A_ji. For every integer
1<=k<=n,

    Delta(A) >= E Delta(T) + gamma_(n,k) V,

    gamma_(n,k)
      = 4[n-k+2-(n+k+2)/2^k] / [n(n-1)].             (1)

The k=1 case is the star coefficient 2/n. Optimizing k near log_2 n yields
max_k gamma_(n,k) ~ 4/n. This improves the inherited ~24/n^2 variance
coefficient in the large-n regime. The sharp H=1 constructions show the
unknown all-order optimal coefficient is at most ~8/n, so (1) establishes
the correct order with an asymptotic factor-two gap.

## One-edge interpolation

For e={i,j}, p=A_ij, q=p(1-p), R=A-{i,j}, and G_epsilon=H-u+epsilon v,
direct determinant expansion gives

    G_epsilon(A)=(1-p)G_epsilon(A^0)+pG_epsilon(A^1)
                  +q(u_R-epsilon v_R).

Here A^0,A^1 are the endpoint orientations of that edge. Taking the minimum
of the two signed margins gives

    Delta(A)>=(1-p)Delta(A^0)+pDelta(A^1)
                   +q(u_R-|v_R|).                 (2)

The endpoint choice and determinant sign conventions must be consistent,
but the minimum form (2) is invariant under reversing the edge gauge.
The remainder is nonempty for n>=3. Its final coefficient is strictly
positive, so any preliminary edge rounding may safely discard that
nonnegative correction. Every rounding order produces the SAME product
law for the terminal binary tournament T.

## Binary-core determinant lemma

**Lemma.** Let B have order r>=1. Suppose all its nonbinary edges have a
vertex cover C of size c. Then

    det(I+B)-|det(B)| >= 2^(-c).                    (3)

Put K=B-B^T, and adjoin an observer with skew entries 1 to all old vertices.
The inherited augmented-skew identities express u+v and u-v as, respectively,
2^(1-r) times the sum of principal Pfaffian squares not containing the
observer, and the sum containing the observer. The binary core has size
m=r-c. Every even principal Pfaffian of a binary skew tournament is an odd
integer: its expansion has (2l-1)!! odd many summands, each +/-1. Hence its
square is at least one, including the empty Pfaffian.

If m>=1, retain only even subsets of the binary core for u+v, and the
observer together with odd subsets of the core for u-v. There are
2^(m-1) such subsets in each case. Every retained augmented Pfaffian is
also an odd integer because observer entries are +/-1. Both margins are
therefore at least 2^(1-r)2^(m-1)=2^(-c).

If m=0, c=r, use the inherited uniform bounds u+v>=2^(1-r) and
u-v>=r 2^(1-r), both at least 2^(-r). This covers the degenerate core
without pretending an empty set has 2^(-1) parity subsets.

## Ordered-center rounding

Choose distinct ordered vertices x1,...,xk. First round every edge having
neither endpoint among these centers, discarding its nonnegative correction
in (2). Then round the still-unrounded edges incident to x1, then x2, and
so on. Each edge touching the center set is assigned to its earliest
center endpoint. If assigned at position t, set s=k-t.

- If its other endpoint is outside the center set, the pair-deletion
  remainder has uncertain-edge vertex cover consisting of the s later
  centers. Its correction in (2) is at least q_e 2^(-s).
- If its other endpoint is a later center, that endpoint also disappears
  in the pair deletion, so the cover has at most s-1 vertices. Its
  correction is at least q_e 2^(1-s).

Previously rounded edge values are binary. Other unrounded star edges
incident to the deleted current center disappear, so no uncertain edge
outside these covers has been overlooked. These bounds hold on every
conditioning branch; q_e remains its original variance before that edge
is rounded.

Now average the ordered centers uniformly without replacement. For a fixed
edge, the probability of assignment at position t with its other endpoint
outside the center set is 2(n-k)/[n(n-1)], and the probability that the
other endpoint is a later center is 2(k-t)/[n(n-1)]. Hence its expected
correction coefficient is

    2/[n(n-1)] sum_(s=0)^(k-1) [(n-k)+2s] 2^(-s),

which evaluates to (1). The terminal binary expectation is independent of
the chosen center set and order, so it is retained once, not multiplied.

The s=0 second case has zero probability; no negative-size cover is used.
For k=n all n-k outside cases vanish; the formula remains valid though
this k is usually far from optimal.

## Scope and consequences

Order two is excluded: Delta=0 for every weighted 2-tournament while V can
be positive. The empty pair-deletion remainder has u-|v|=0 and is precisely
the failed coefficient. The binary-core lemma itself has r>=1.

Any nonnegative binary lower bound Delta(T)>=F(T) can be transferred as
Delta(A)>=E F(T)+gamma V. For example a sum of directed-cycle indicators
multiaffinely extends to the corresponding weighted cycle sum. This
transfer is only as strong as the separately proved binary F; it does not
assume a weighted even-cycle recurrence that has not been proved.

For n>=4 one explicit simplified coefficient follows by k=ceil(log_2 n):

    gamma >= 4[n-ceil(log_2 n)]/[n(n-1)],

because n+k+2<=2n<=2^(k+1). The exact optimized (1) is preferable at small
orders. Consecutive numerators differ by

    -1+(n+k+1)/2^(k+1),

which locates the maximizing integer k without an unbounded search.

## Matching-uncertainty corollary of the graph-rounding framework

The companion root proof `../rounding_graph/PROOF.md` supplies the general
independence-polynomial remainder bound

    u-|v| >= 2^(1-r) min(I_even(G),I_odd(G)), r>=1,

where G is the graph of uncertain edges in the remainder. We specialize
that framework here rather than duplicate its general proof.

For a matching of c uncertain edges and t binary isolated vertices in G,
the independence polynomial is (1+2z)^c(1+z)^t. Formula (4) becomes

    u-|v| >= (3/4)^c,                         t>=1;
    u-|v| >= (3/4)^c-(1/4)^c,                t=0.

These bounds are sharp, attained by ordinal sums of fair pairs and
singletons. More generally, suppose A's entire uncertainty graph is a
matching of r>=1 edges (all other edges may form arbitrary binary cross
orientations). Randomly order and round those r edges. For n>2r, every
pair-deletion remainder has at least one isolated vertex in its uncertainty
graph, so an edge at position j gains q_e(3/4)^(r-j). For n=2r, only the
first remainder has no such isolate, subtracting q_e(1/4)^(r-1).
Averaging the uniformly random order therefore proves

    Delta(A) >= E Delta(T) + a_r V,          n>2r;
    Delta(A) >= E Delta(T) + b_r V,          n=2r>=4,

where

    a_r=(4/r)(1-(3/4)^r),
    b_r=(4/r)(1-(3/4)^r-(1/4)^r).

Thus the conjectured dimension-sharp constants are proved on the full
matching-uncertainty class, not merely on its H=1 subclass. The binary
cross-edge orientations need not be transitive. This is a scoped theorem;
it does not prove the same coefficient for arbitrary uncertainty graphs.
