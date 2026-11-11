# From census error to predictions for whole metric distributions

**PROVED IN NOTE.** Elementary orbit counting and positive response transfer.
No external novelty claim. The exact labelled response mean is inherited from
Session270's carrier proof; the mixture and error propagation make it usable
for ordinary-class prediction without suppressing automorphism weights.

## 1. The missing automorphism weight is an exact positive measure

Let C_n be the ordinary tournament isomorphism classes, a_n=|C_n|,
N=binom(n,2), b_n=2^N/n!, and

    delta_n = a_n/b_n - 1 >= 0.

Let U_n be uniform on C_n. Let L_n be the distribution on C_n induced by a
uniform labelled tournament. If g(T)=|Aut(T)|, orbit-stabilizer gives

    U_n(T)=1/a_n,       L_n(T)=1/(b_n g(T)).

Define s_n=sum_T(1-1/g(T))=a_n-b_n. When s_n>0 define

    D_n(T)=(1-1/g(T))/s_n.

It is supported on symmetric classes. Direct substitution proves

    U_n = L_n/(1+delta_n) + delta_n D_n/(1+delta_n).          (1)

If delta_n=0, U_n=L_n and the second term is absent. In particular,

    TV(U_n,L_n) <= delta_n/(1+delta_n).                       (2)

For every event E, |U_n(E)-L_n(E)| is bounded by (2), with no hypothesis
on E: primeness, strongness, regularity, a prescribed H interval, a spectrum,
or a property of the modular tree. For every bounded real f,

    |E_U f-E_L f| <= delta_n/(1+delta_n) * (sup f-inf f).     (3)

These are absolute probability/expectation errors. They do not give relative
accuracy for events rarer than the symmetry correction. In particular,
regular tournaments are already an exceptionally thin event; conditioning on
them can magnify the bound by the reciprocal of their probability.

### A direct interpretation of the Burnside correction

Under L_n,

    E_L |Aut(T)| = a_n/b_n = 1+delta_n.

Thus delta_n is exactly the mean number of nonidentity automorphisms, not just
an analytic remainder. Markov gives L_n(nonrigid)<=delta_n. Because tournament
automorphism groups have odd order, a nontrivial group has at least three
elements, yielding the stronger L_n(nonrigid)<=delta_n/2. This uses the
standard elementary observation that an involution would reverse an arc,
and Cauchy's theorem. The mixture explains why symmetric classes are more
visible under U_n, while their absolute share is still exponentially tiny.

## 2. Predict every path-cover coordinate, not merely H

For n>=1 let q_r(T) count ordered covers by r nonempty directed paths,
1<=r<=n. Choosing a vertex permutation and r-1 cut positions shows

    M_(n,r) = n! binom(n-1,r-1),
    0 <= q_r(T) <= M_(n,r),
    mu_(n,r) := E_L q_r = M_(n,r) / 2^(n-r).                 (4)

The n-r required arcs are distinct; each has probability 1/2. In polynomial
form, with Q_T(t)=sum_(r=1)^n q_r(T)t^(r-1),

    E_L Q_T(t) = n! 2^(1-n) (1+2t)^(n-1).                  (5)

Combining (1) and (4),

    1/(1+delta_n) <= E_U q_r / mu_(n,r)
                   <= (1+delta_n 2^(n-r))/(1+delta_n).      (6)

Thus the downward relative error is at most delta_n/(1+delta_n), and the
upward relative error at most delta_n(2^(n-r)-1)/(1+delta_n).
A rigorous upper bound B_n>=delta_n may replace delta_n in either error.
For q_n and q_(n-1), use their universal exact values rather than this loose
range bound. Every q_r has a huge magnitude; the relative statement is the
useful prediction, not exact integer rounding.

In particular, the ordinary-class mean Hamiltonian count is predicted by

    E_U H ~= n!/2^(n-1),

with upper relative error at most B_n(2^(n-1)-1)/(1+B_n). A Burnside error
of order n^3 2^(-2n) becomes a metric-mean error of order n^3 2^(-n), still
exponentially small. The absolute mean error can nevertheless be enormous.
This repairs the obstacle noted at the end of Session270's carrier proof:
the missing weights are controlled, not declared to vanish.

The two scales matter: bounded event probabilities transport with order
n^3 2^(-2n), whereas this coarse bound for H means loses a factor 2^(n-1).
Large moments can lose much more. No concentration theorem for H follows
from its first moment.

## 3. Positive recursive contexts transport forecast intervals

Fix a quotient tournament Q with d marked ports. Draw its children T_i
independently, with fixed orders n_i. The inherited substitution law has form

    H(Q[T_1,...,T_d]) = F_Q(q(T_1),...,q(T_d)),              (7)

where F_Q has nonnegative coefficients and is multilinear in the d response
vectors: every monomial uses exactly one coordinate of each nonempty child.
Deleting all intermodule arcs from a Hamiltonian path proves this directly:
each module contributes one ordered internal path cover; quotient words
supply the nonnegative integer gluing coefficient.

Independence therefore implies E H=F_Q(E q(T_1),...,E q(T_d)). Let P be the
result of inserting the exact labelled means (4). Let delta_i be the census
correction at n_i and epsilon_i=delta_i(2^(n_i-1)-1)/(1+delta_i). Then for
independent uniform ordinary-class children,

    P product_i (1+delta_i)^(-1)
       <= E H <= P product_i (1+epsilon_i).                 (8)

Sharper bounds use coordinate-dependent factors from (6). Monotonicity and
multilinearity prove (8) term by term. The same argument applies to each
parent path-cover coordinate, so intervals propagate through a tree whose
primitive leaves are independently sampled. At a large structured hierarchy,
errors accumulate by leaf multiplicity, not merely by depth.

**Important scope:** this samples a marked assembly with independently drawn
children. It is not uniform over resulting unmarked parent classes. Equal
child types and quotient automorphisms require the response-weighted cycle
index from Session268. Canonical uniqueness does not make a product sampling
law uniform. A nonuniform seam may have a signed formula; cancellation can
destroy a relative error estimate. There one must propagate absolute boundary
state intervals or use a positive endpoint-state realization.

## 4. A useful large-scale forecast dictionary

- A000568 magnitude: log b_n plus the tiny log(1+delta_n) correction.
- Rigidity: mean nonidentity automorphisms is delta_n; nonrigid probability
  under labelled sampling is at most delta_n/2.
- Primeness: the structural lane's exact module union bound transports to
  ordinary classes with the additive error in (2).
- Full path-cover response: polynomial (5), with all-coordinate errors (6).
- Recursive random modules: positive contraction (8), retaining marking.
- Two-core/ordered-spine words: a different ensemble; the Chan–Pak lane's
  word-growth law applies there. Uniform words, labelled tournaments and
  ordinary tournament classes must never be interchanged.

## Sources and boundary

Inherited sources, accepted in their stated scope:
`../session270_chanpak_enumeration_20261111/carrier/PROOF.md` (relative to
this session's parent results directory) for (4) and the marked carrier;
Session268 `modular/PROOF.md` for response contraction and weighted orbit
counting; its `spectrum/PROOF.md` for contextual completeness of q. Burnside
and orbit-stabilizer are classical. This note's bounds are derivations, not
an external priority claim and not an exact A000568 validation campaign.
