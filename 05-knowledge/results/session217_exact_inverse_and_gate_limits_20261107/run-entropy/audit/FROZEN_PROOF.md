# Reversal product and the logarithmic run defect

Status: PROVED low-run zero-density obstruction and elementary identities; no all-target occupancy theorem.
The three-core inverse conjecture remains OPEN. No literature novelty claim.

## Inheritance and scope

Input: Session216 three-core-portable/PROOF.md formulas (1)--(4) and its
80 core-cut states; Session216 numeration-portable/PROOF.md supplies the exact
inverse reductions and the scoped uniform radix-rewrite obstruction. The
closest mechanism is the cut-coordinate returning compiler. Hostiles are the
permanent grammar omissions 47 and 63 and the balanced-band misses 265 and 287.
The least-used sidecar here is the reversed cut coordinate G, rather than
only the output H or intermediate A. Concept board: exact occupancy; binary
reversal products; run defects; target-dependent inverse division. These
identities reorganize the same grammar and do not add constructions.

## 1. Exact defect of the binary quadratic

For a nonnegative integer n written with a binary digits (leading zeroes
allowed), define v(n)=rev_a(n)/2^a and

    delta(n)=n v(n)-f(n).

Both quantities are independent of leading zeroes. For n=0 use v=delta=0.
Writing epsilon=n mod 2, direct use of the inherited f recurrence gives

    v(2n)   =v(n)/2,       delta(2n)=delta(n),
    v(2n+1) =(1+v(n))/2,
    delta(2n+1)=delta(n)+(1+v(n))/2-epsilon.       (1)

Partition the binary expansion, read from the most significant end, into
maximal runs of ones, of lengths l_1,...,l_r. Let u_j be v of the prefix
strictly preceding run j, including the immediately preceding zero if present.
Thus 0<=u_j<1/2. Iterating (1) through a run proves

    delta(n)=sum_j [u_j+(1-u_j)/2^l_j].           (2)

To see the run increment, its first bit contributes (1+u)/2. The subsequent
l-1 bits contribute minus (1-u)/2^2 through minus (1-u)/2^l; their geometric
sum leaves u+(1-u)/2^l. Zeros leave delta unchanged. Therefore

    0<delta(n)<=(3/4)r <= (3/4)ceil(a/2)   (n>0),
    f(n)=n v(n)-delta(n).                         (3)

The first inequality is strict for n>0, and the 3r/4 bound is deliberately
not claimed sharp. The number of one-runs in an a-bit word is at most ceil(a/2).
For the one-run word n=2^a-1 the defect is exactly 2^-a; it has no uniform
positive lower bound. For alternating words n_m=(4^m-1)/3, with 2m-1 digits,

    delta(n_m)=2m/3-2/9+2/(9*4^m).                (4)

Indeed before run j the reversed prefix is (1-4^{-(j-1)})/3, so (2) sums
an elementary geometric progression. Thus logarithmic growth genuinely occurs:
there is no constant-error replacement f(n)=n v(n)+O(1).

## 2. Full cut-coordinate product, with only logarithmic error

Fix one core state (h,q,A,G,T,S,E). Write prefix P in a digits and suffix R
in b digits, with t=rev_a(P), s=rev_b(R). Empty padding has length and reversed
word zero. Set

    alpha=2^b(qP+A)+hR,
    gamma=2^a(qs+G)+ht,
    Q=q 2^(a+b).

These are exactly the final body's cut coordinates A_L,G_L,q_L from the
inherited ordinal composition formula. Its two-boundary coordinate is

    tau=P G+h f(P)+T+(P mod2)S
          +(qP+A)s+h f(R)+E first_b(R).

Define D=alpha*gamma/Q-tau. Expanding the product and applying (3) proves

    D=h[delta(P)+delta(R)] + B,
    B=AG/q + hA v(P)/q + hG R/(q 2^b)
          +h^2 R v(P)/(q 2^b)-T-(P mod2)S-E first_b(R).  (5)

Here R/2^b=0 when b=0. Consequently the exact strong output count is

    N=alpha(2-gamma/Q)+h+D.                       (6)

The core-dependent correction B is bounded uniformly in both padding lengths:

    -T-S-E <= B <= (AG+hA+hG+h^2)/q.

Thus the difference from the explicit reversal-product main term is at most
linear in a+b, rather than the exponential size of the individual terms:

    -T-S-E <= D
       <= (3h/4)[ceil(a/2)+ceil(b/2)]
                          +(AG+hA+hG+h^2)/q.      (7)

All claims include empty padding and empty core cuts algebraically; strongness
must still be imposed exactly as in the inherited grammar. They do not infer
strongness from a numerical count. The formula works for any nonempty fixed
cut-labeled core with the stated coordinates, not just the bank of 80.

## 3. What this does and does not reduce

Equation (6) is an exact reformulation, with D retained. It makes an analytic
route precise: first control a multiplication-like image of a cut encoding
and its reversal, then resolve its logarithmic run correction. Simply rounding
the main term is invalid: (4) shows unbounded error even for one padding side.
Nor does a product near the target imply exact representation. The same cut
word determines both alpha and gamma, so treating them as independent factors
would discard the principal constraint. No uniform selection, interval
coverage, new all-length exclusion, or finite-state occupancy certificate is
claimed here. This is a structural sidecar, not a decisive reduction to a solved
number-theory problem.

## 4. An all-length obstruction to low-run covering

For a legal strong recipe let r be the total number of maximal one-runs in
its two padding words; the chosen core cut is not charged. Fix 0<=rho<1/4.
Let B_rho(X) count distinct targets N<=X which have at least one legal recipe
with r<=rho log_2 N. Then

    B_rho(X) <= X^{H_2(2rho)+o(1)} = o(X),        (8)

where H_2(t)=-t log_2 t-(1-t)log_2(1-t) and H_2(0)=0.
Thus no grammar restricted to that run budget covers all sufficiently large
odd targets. More explicitly, if L=floor(log_2 X) and R is a nonnegative
integer, the number of targets <=X admitting a recipe with at most R runs is
at most

    40(L+1)(L+2) sum_{j=0}^R binom(L+2,2j),       (9)

with binomial coefficients zero outside their usual range. In particular a
fixed run budget R gives only O((log X)^{2R+2}) different targets.

Proof. The inherited positivity bound gives

    N >= 2^b(qP+A)+h.

If a>0, legality implies P>=2^{a-1}, and q>=32 for this bank, so
2^{a+b}<= (N-h)/16 < X. If a=0, the nonempty core cut gives A>=1,
so again 2^{a+b}<=N-h<X. Hence k=a+b<=L for every legal recipe
with target at most X. Encode its padding by concatenating P, a separator
zero, and R as binary words. The result has k+1 digits and exactly r one-runs.
For each fixed split a,b this is injective. The number of n-digit binary
words with exactly j one-runs is binom(n+1,2j), by choosing the 2j endpoints
of the runs among the n+1 gaps (intervals are nonempty and separated).
Thus for each split the count is at most sum_{j<=R}binom(k+2,2j), which
is bounded by the same sum with k+2 replaced by L+2. Sum over the k+1
splits and k=0,...,L and multiply by 80 core states to obtain (9).
This counts recipes, so collisions can only improve the target upper bound.

For (8), put R=floor(rho log_2 X). For all large X, 2R/(L+2)<1/2.
The binomial entropy estimate
sum_{j=0}^{m}binom(n,j)<=2^{n H_2(m/n)} for m<=n/2 follows, for example,
by weighting the binomial expansion at p=m/n: each weight p^j(1-p)^{n-j}
with j<=m is at least p^m(1-p)^{n-m}. Applying it with n=L+2,m=2R
upper-bounds the even-index sum and proves (8); rho=0 follows directly.

The threshold 1/4 is the entropy threshold of this counting argument, not
a proved sharp boundary for actual output support. Nothing here proves that
run budget >=(1/4)log_2 N is sufficient, or that every individual N needs
that many runs. The quantifier is that, for every fixed rho<1/4, all but
o(X) integers up to X have no low-run representation. This remains true
without assuming the full three-core inverse conjecture. If that conjecture
is true, almost every odd target must use high-run-complexity padding.

The same theorem holds for any fixed finite bank of nonempty strong cores
and cut states, with bank size B replacing 80 and prefactor
B(L+1)(L+2)/2 replacing 40(L+1)(L+2). For a>0, q>=2 implies
2^{a+b}<=2^b qP<=N-h. For a=0, a nonempty cut gives A>=1
by splitting a Hamilton path immediately before a selected vertex. No other
specific constants of the three-core bank are used in the counting proof.

This obstructs bounded-run automata as an all-target shortcut, not general
finite-state languages: a finite automaton can allow linearly many runs.
Complementary choices can preserve high run complexity, so (8) does not
rule out target-dependent complementary recoding.

## 5. Reproduction and novelty check

Run `python3 verify.py` in this directory. It independently evaluates f by
selected bit pairs and tests all a-bit words for a<=13, the run formula, the
alternating formula through m=40, and all 80 core states for a+b<=8. Formula
(5) is checked with rational arithmetic against the inherited exact tau, and
(6) against the inherited closed expression, not a fresh tournament DP.
This is a verification of algebraic identities; no spectrum coverage scan.

The local Session216 three-core, numeration, arithmetic, and density proofs
were checked. The density packet routes to prior canon THM-4885-one-apex-join-representation-mass-pi2-over-12.md and
THM-4890-one-apex-collision-balance.md: the former bounds a one-apex
reversal-product error by the number of 10 transitions; the latter (L0)
gives an exact local run defect. Consequently sections 1--2 are presented as
an elementary three-core specialization/extension of that existing mechanism,
not as a new general binary-product or logarithmic-error discovery. Targeted
run/entropy searches in the local Session216 sources found no matching
three-core low-run obstruction, but no unrestricted novelty or priority claim
is made. No changes were made to the dirty math tree.
