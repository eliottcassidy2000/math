# Quantitative run complexity for fixed-apex spectrum constructions

**PROVED counting/encoding consequences of cited strong bounds; finite tests
are controls, not an existence proof. Universal spectrum completeness OPEN.**
The fixed-run density-zero mechanism is inherited from the feedback-word
packet's `run_density.md`. This note makes the resource grow with the target,
using the exact induced core developed in Session133. No priority claim.

## 1. Actual presentations and monotone resource

A k-apex presentation consists of an ordered marked set A, |A|<=k, such that
T-A is transitive. In its unique chain order, encode each vertex by the bits
specifying which marked vertices beat it. If |A|=j this is a word on q=2^j
letters; the marked tournament itself has 2^binom(j,2) choices. Define the
change count R(w)=#{i:w_i!=w_(i+1)}, zero for empty/singleton words.
This is a presentation resource, not an invariant of an unmarked graph unless
minimized over all such presentations.

Delete every chain vertex lying in a singleton SCC, but retain all marked
apices, including singleton ones. H is unchanged. The remaining chain word
is a subsequence of the original, and R cannot increase: the inequality
1[a!=c]<=1[a!=b]+1[b!=c] controls deleting an interior letter. If one also
discards singleton apices, coordinate projection can only identify letters.

Every remaining nonsingleton SCC C_i contains a_i>=1 marked apices: a
subtournament of the transitive chain alone cannot be strongly connected
and nontrivial. If n_i is its order, the surviving spine length is

    m=sum_i(n_i-a_i)<=sum_i(n_i-1).

Busch's cited bound H(C_i)^3>=5^(n_i-1), multiplied over SCCs, therefore gives

    H(T)^3>=5^m,  hence m<=floor(log_5(H(T)^3)).

This is stronger for marked-spine purposes than using the whole-core bound
K(X)<=3 ln(X)/ln3. It applies to every target, not only prime counts, and
retains exactly the original marked apex tournament. The empty spine is
allowed; if H=1 then m=0. The fixed marked graph itself must still be
compatible with that target. This compression preserves scalar H and the
run budget, not the complete future-response state. Long trailing sink
words can have the same H and different continuation probes; the logarithmic
bound concerns existence of a reduced representative, not every encoding.

This refinement was obtained jointly with the sibling two-apex decoder
lane after checkpoint1. The earlier K(X) bound and prime-only stronger bound
remain correct; this is an improvement, not a retraction. Counting every
presentation rather than quotienting by isomorphism only enlarges an upper
bound on the number of distinct H values.

## 2. Exact finite encoding bound

For integers N,R>=0 set

    E_k(N,R)=sum_(j=0)^k 2^binom(j,2) *
      [1+2^j sum_(r=0)^min(R,N-1) binom(N,r+1)(2^j-1)^r].

The sum is empty when N=0. There are exactly
q binom(m-1,r)(q-1)^r words of length m with r changes: choose the first
letter, change positions, and a different letter at every change. Summing
m=1,...,N uses the hockey-stick identity and adding the empty word gives the
bracket above. Thus the number of positive H<=X with a <=k-apex presentation
of change count <=R is at most E_k(floor(log_5(X^3)),R). The logarithmic floor is evaluated
with integer powers, avoiding numerical ambiguity. The fixed apex graph
choices are deliberately overcounted; its number of vertices does not
consume any of the spine-length bound.

It is NOT a lower bound on the number of attainable values. Different words
or different apex markings may have the same H. No identification of these
fibres is required for an upper bound.

## 3. Growing resources, not merely fixed templates

Fix k. If R(X)=o(ln X), the preceding count is X^o(1). Indeed N=O(ln X), and
for R<N/2 the elementary bound binom(N,s)<=(eN/s)^s gives a logarithm of
order R log(N/R)+R+log N=o(ln X). The fixed number of marked-graph choices
and letter choices do not change that conclusion. Bounded or sublogarithmic
change counts therefore describe only a subpolynomial population of targets.

More sharply, let q=2^k, C>0 be the relevant order coefficient, and
R=c ln X+o(ln X), with 0<c/C<1-1/q. Standard binomial entropy gives

    E_k(C ln X+O(1),R) <= X^(beta+o(1)),
    beta=C[-rho ln rho-(1-rho)ln(1-rho)+rho ln(q-1)],
    rho=c/C.

The formula follows equally from factorial bounds in the exact finite sum;
it is a counting upper bound, not a probabilistic assumption on actual words.
For every target take C=3/ln5 using the spine-specific bound. Whenever beta<1,
these low-change presentations cover a zero-density subset even of the primes
(by the classical [prime number theorem, DLMF27.2.3](https://dlmf.nist.gov/27.2.E3)). This does not prove the remaining
prime values are attained: they may be unattained, or require more changes.

## 4. Explicit two-apex consequences without numerical root fitting

For k=2 and all targets, use C<2 and R=floor((ln X)/5). The binomial bound
already gives exponent at most

    (ln30+1)/5 < 9/10.

The strict inequality follows from exp(7/2)>30, certified by a finite rational
partial sum in the checker. Therefore only X^(9/10+o(1)) values up to X can
have a two-apex presentation with <=(ln h)/5 changes, where h is that value. In particular, for
almost all primes, either no two-apex presentation exists or every such
presentation has more than (ln p)/5 changes. This is a necessary resource
condition, NOT an assertion of prime realization.

The weaker checkpoint1 bound C<3 and R=floor((ln X)/6) give exponent less than5/6 since
ln54<4. Again exp4>54 has a rational partial-sum certificate. Thus the earlier valid bound was
X^(5/6+o(1)) target values admitting that restricted resource budget.

A universal two-apex construction, if one exists, must use profiles with
logarithmically many changes on a large target population. This does not
rule out two apices: the word length itself may be logarithmic, allowing
that many changes. It rules out replacing the remaining problem by fixed
or sublogarithmic run patterns.

## 5. Checks and precise scope

`check.py` verifies exact run histograms, the hockey-stick count, deletion and
coordinate-projection monotonicity, and all marked <=2-apex presentations
of tournaments through order five against exact SCC deletion and path DP.
Normal and optimized receipts agree. The all-order proof is the encoding
and counting argument above, not an extrapolation of finite examples.

The arc-budget entropy theorem in the sibling interval lane is complementary:
a bounded number of apices may have an unbounded number of backward arcs.
Neither barrier establishes a missing global H value beyond7 and21, nor
closes the prime or7p realization obligations.
