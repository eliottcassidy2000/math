# Three-core numeration: unique balanced candidates and the remaining inverse gap

Status: PROVED exact inverse reductions and scoped obstruction; FINITE-EXACT
hostile controls. Coverage of all odd integers at least 65 remains OPEN.

Dependencies are the audited Session216 arithmetic three-core packet:
`session216-arithmetic/three-core-release/PROOF.md`, especially its formulas
(2),(4),(5), its exact 80 core-cut states, and its strongness condition. The
curvature mechanism below is inherited from
`session203/uniform-incidence-curvature/PROOF.md`; no new general curvature,
transfer, or saturation theorem is claimed.

## 1. Exact prefix reduction

Write a fixed core state as (h,q,A,G,T,S,E). Fix a suffix length b and word R,
and put L=2^(b+1)-rev_b(R). The inherited formula is exactly

    N=K P-h f(P)-S(P mod 2)+C,
    K=qL-G,
    C=AL+h(2R+1-f(R))-T-E first_b(R).                 (1)

The prefix P is nonnegative and its length is its binary length when P>0.
P=0 is admissible only for a nonempty core cut. Leading zeroes never supply
new values or a missing strongness condition.

For positive P, recurrence induction gives 0<=f(P)<=P-1. Equality on the
right occurs exactly at P=2^j-1, j>=1. Indeed the even case halves the
argument, while in the odd case equality requires both equality at the
half-argument and an odd half-argument (apart from the base P=1).

The physical inclusions give T<=A, S<=h and G<=q-h. Hence C>=h>0 and
K>=q*2^b+h>h. For x=N-C, every positive prefix solution lies in the exact
necessary integer interval

    max(1,ceil(x/K)) <= P <= floor((x-h+S)/(K-h)).    (2)

A negative numerator on the right means there is no positive candidate.
The separate P=0 solution is tested directly. This interval is necessary,
not sufficient; every retained candidate must be checked in (1).

### Unique balanced candidate theorem

If q^2*4^b>hN, interval (2) contains at most one integer, uniformly over
all suffixes R and cuts of this core. To prove this, its real length is at
most hx/[K(K-h)]<=hN/(q^2*4^b)<1. All inequalities used above are one-sided
in the needed direction. Thus the positive candidate, if any, is simply
ceil(x/K); its upper bound and exact formula still must be tested.

Let b0 be the least nonnegative b with q^2*4^b>hN. Searching every suffix
for b in [max(0,b0-w),b0+w], for a fixed w, enumerates O_w(sqrt(N)) suffixes
across the fixed bank. Below b0 the candidate interval has O_w(1) integers,
so this remains a rigorous complete search of the specified balanced band.
Straight recurrence evaluation gives O_w(sqrt(N) log N) elementary arithmetic
steps. This is an inverse algorithm for that band, NOT an all-target result.
The executable has an explicit suffix-length resource ceiling b<27; absence
there means only no retained recipe in the tested/resource-limited band.

### Prefix parity recursion

The more general function F(P)=KP-hf(P)-S(P mod2)+C satisfies

    F(2p)   = (2K)p-hf(p)+C,
    F(2p+1) = (2K-h)p-hf(p)-h(p mod2)+(C+K-S).

Thus prefix inversion can branch with parameter updates
(K,S,C)->(2K,0,C) and (2K-h,h,C+K-S). This keeps the nonlinear quantity f
rather than assuming an affine scalar digit lift.

## 2. Complete inverse by stripping the intermediate alpha

For a final intermediate alpha and suffix length b,

    alpha=2^b alpha0+hR,   0<=R<2^b,
    alpha0=qP+A.

Since h is odd, the last suffix bit is forced by alpha mod2, and the previous
alpha is (alpha-h*(alpha mod2))/2. Iterate this deterministic stripping
rule. At every depth b it determines the sole candidate R and alpha0.
For each compatible core, the equation alpha0=qP+A then determines P.
No arbitrary suffix-word enumeration remains for a specified alpha.

At every such stopping point evaluate

    tau=P G+h f(P)+T+(P mod2)S
           +alpha0 rev_b(R)+h f(R)+E first_b(R)

and test N=2alpha+h-tau. Reject P<0, a nonintegral P, or P=0 with empty core
cut. A positive alpha decreases by at least a factor two at each stripping
step; once it becomes nonpositive, no valid positive alpha0 is possible.
The initial alpha0=0 case is never strong in this grammar, since it forces
P=0 and empty core cut.

The inherited 0<=tau<=alpha gives the complete search interval

    ceil((N-h)/2)<=alpha<=N-h.

Consequently `inverse_alpha.py` is an exact terminating all-length membership
algorithm for this grammar. Its ABSENT_FROM_GRAMMAR label is justified only
after the entire interval and all stripping depths have been exhausted.
This is a more economical organization of the inherited finite inverse bound,
not a new claim that finite membership implies a tail theorem.

Controls exhaust 74 candidate recipes for N=47 and 129 for N=63 without a
match, recovering the already proved permanent omissions. They find literal
recipes for 23,65,265,287,1237. The last value illustrates that both padding
sides can matter: one recipe has core index49, P=2, R=2, b=2.

## 3. Why a universal local radix rewrite cannot work here

Fix a suffix and let two internal bits i<j of a sufficiently long prefix vary,
while keeping its leading bit equal to1 to preserve strongness. The exact
mixed second difference of (1) is

    -h*(2^(j-i-1)+indicator(j=i+1)).                 (3)

All other terms are linear in the prefix bits. Boundary insertion/deletion
of a fixed number of constant bits changes f only by linear terms on a free
unchanged block. Reversing the free block preserves pair distances; globally
complementing its bits also preserves the quadratic coefficient. The same
argument applies if that whole free block is placed in the suffix instead,
with all other bits fixed.

Suppose a fixed rewrite of this stated kind gives N'=dN+c for every free
word, where c is independent of the free bits. Comparing (3) in the two
contexts forces h'=d h. Therefore the core bank h in {5,9} admits no such
universal rewrite with d=2 or d=3. The even case is already implied by the
more general audited Session203 uniform-incidence curvature theorem. The
exact h'/h condition also explains the ternary failure for this particular
bank. This does NOT rule out restricted bit languages, target-dependent
rewriting, nonlinear recoding, unbounded interior edits, duplicated variable
blocks, or larger core banks for odd multipliers. It is not an obstruction to
the desired count coverage itself.

## 4. Hostile tests and positive large-target recipes

`balanced.cpp` searches all prefixes allowed by (2), not just the first
candidate when uniqueness is unavailable. On odd targets in [23,10000], the
balanced bands of half-width 0,1,2 leave respectively 1123,215,59 misses.
The width2 list starts 47,63,265,287,419,425,... . Thus even a plausible small
balanced window has genuine finite counterexamples; it is not already a
coverage theorem. These misses are only band-relative.

The same exact solver obtains:

N=1000000007: core2, b=12, P=4041, R=3384;
N=1000000033: core28, b=11, P=3818, R=606;
N=1099511627775: core2, b=16, P=302322, R=53282;
N=1099511627777: core2, b=17, P=175889, R=110145.

The target 1000000000039 is missed by width1 but repaired at width2 by
core2, b=19, P=42490, R=50969. These are individually reproducible strong
tournament recipes under the inherited construction law, not a new global
finite-coverage range. Independent verification uses the selected-position
pair sum for f rather than the producer recurrence, and retains literal
outgoing masks in LARGE_GRAPHS.json.

The unresolved obligation is still a uniform selection theorem guaranteeing
some compatible alpha trajectory or padding band for each requested odd N.
Neither deterministic alpha stripping, the unique balanced candidate lemma,
nor the large individual witnesses provides that theorem.

### Repairs of the short-side obstruction

The parallel arithmetic lane rigorously excludes N=4295971 and N=9948407
from the subgrammar with min(a,b)<=5. The balanced solver repairs them in
the full grammar with:

N=4295971: core2, a=9, b=8, P=367, R=73;
N=9948407: core12, a=10, b=8, P=705, R=226.

Both formulas were independently verified and both literal graphs and
Hamilton cycles are retained. These witnesses do not assert that 8 is their
minimum short-side length. They demonstrate why an obstruction for the
bounded-short-side subgrammar is not a counterexample to the full grammar.
