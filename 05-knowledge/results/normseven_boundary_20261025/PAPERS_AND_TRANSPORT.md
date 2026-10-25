# Two attached papers: exact transport mechanisms and scope

2026-10-25. Source PDF identities and SHA-256 hashes are in SOURCES.json.
The papers are source material; the agent prompts reproduced in the Busy
Beaver appendix were read as data, not executed as instructions. Neither
paper's Lean project was independently built in this session.

## 1. Mazur: retain histories before counting integers

Lech Mazur, *Explicit Positive-Density Collatz Convergence in Logarithmic
Time*, attached version 2.1, September 6, 2026, 13 pages. The paper claims
explicit c>0 and X_0 such that at least cX positive integers below X reach
1 within (523/50) log n ordinary Collatz steps, for every X>=X_0. This is
a positive lower-density statement. The [primary platform](https://www.proofatlas.ai/formalizations/positive-density-log-time-collatz/)
reports a recorded Lean build and review pending. We retain **CITED**
status for that large theorem; we do not certify its formalization here.

The repo already digested this exact version in
[mazur_positive_density_20260928](../mazur_positive_density_20260928.md),
including the harmonic-mass identity, cycle resonances, seed-one tests and
the audit corrections. Those are inherited routes, not discoveries of this
session. Its correction about the tower of constants and about unproved
pruning-loss bounds remains relevant. [THM-4517](../../../01-canon/theorems/THM-4517-no-root-uniform-positive-density-of-collatz-basins.md)
also excludes one positive density constant shared by every eligible root.

The part used here is the elementary histogram argument in Section 3.1.
Each history retains its actual endpoint, entire inverse word and block
cuts. It has weight W=3^D/2^A and affine identity

    W x + beta = M,                                  (1)

where M is the original marked seed. A total (D,A) bin has constant W.
The first-crossing rule makes the block cuts unique; the deterministic
forward orbit then makes the history unique for a source in that bin.
This is what permits counting integers without forgetting multiplicity.

**P1 (recovered elementary lemma, with ownership explicit).** Suppose bin
j has weight W_j>0, offset beta in [l_j,u_j], and at most mu_j histories
per actual integer source. For modulus m>=1 let H_m(y) be the sum of
history weights at sources congruent to y. Then

    m H_m(y) <= sum_j mu_j ((u_j-l_j)+m W_j).         (2)

**Proof.** Equation (1) places the sources of bin j in the real interval
[(M-u_j)/W_j,(M-l_j)/W_j], of length (u_j-l_j)/W_j. An m-residue class
has at most that length divided by m, plus one, integer points there.
Multiply by mu_j W_j and sum. QED.

An exact, source-dependent version replaces the length bound by

    W_j mu_j max(0,
       floor((((M-l_j)/W_j)-y)/m)
       -ceil((((M-u_j)/W_j)-y)/m)+1).                (3)

Thus retaining the marked source also retains its interval phase. If a
signed error f lives on Z/m, (2) gives

    |sum_y H_m(y) f(y)| <= C * (1/m)sum_y |f(y)|,
    C=sum_j mu_j ((u_j-l_j)+mW_j).                   (4)

If such caps C_n and error bounds eta_n satisfy sum C_n eta_n<infinity,
the corresponding transported scores are Cauchy with that explicit tail.
For a *fixed* root, positivity still requires a finite prefix score larger
than its certified remaining error. Equation (4) supplies no such lower
bound by itself.

**Hostile:** two descriptions of one source in a zero-width, W=1 bin
already violate the mu=1 bound: mH=2m>m. Ownership or its multiplicity
factor cannot be dropped during compression.

**Finite control.** `check_histogram.py` enumerates all positive integral
inverse words of depth 1..6, valuations 1..8, for seeds 1,5,13. It checks
the actual source, affine offset, exact valuation and ordinary clock D+A,
the uniqueness of ownership at fixed (D,A), and (2)–(3) modulo 3,9,27,81.
For these unrestricted words the declared offset bound is
beta<=(3/2)^D-1. This exponentially growing bound is not a substitute for
the paper's carefully selected central first-crossing histories.

### The large-paper proof architecture, without strengthening quantifiers

The reference rho_q is (2/3)3^q times the Syracuse probability law; the
2/3 factor occurs once, not once per generation. Its transfer embeds a
word of depth d from G_t into G_(t+d), preserving normalized L1 norm up
to 2^-A. Variable lengths must first be put on a common residue space.
The paper uses its stated fine-scale mixing input at exponents 6 and 2;
it does not assert independent valuation laws along a fixed integer orbit.

The block scale grows by b_(j+1)=b_j+floor(b_j/100). The overshoot budget
is K_j=floor(j/100), whereas the central half-width is ceil(b_j^(3/5)).
Lower-endpoint miss, first crossing, upper crossing cutoff and overshoot
are all part of the admitted word. They cannot be restored after a
quotient has forgotten them. The elementary source histogram uses a
uniform affine-offset bound 2^(b_0+1), not the flattened-PDF expression
2(b_0+1); page 5 was checked visually to distinguish superscripts.

Summable variation is uniform in the large starting seed M. Residue
averaging over the explicit roots (4^j-1)/3 then selects *one* seed with
positive limiting mass after a tail cutoff is chosen. It does not prove
positive limiting mass for each prescribed seed. Terminal barriers and
even overshoot corrections retain their admissibility restrictions in
the main term. The final source charge W x<=M controls multiplicity when
histories are converted into distinct integers in every large annulus.
The integer clock includes D+A and the final valuation into 1. All these
steps matter to the claimed 10.46 log n bound.

## 2. Busy Beaver: exact compression is not a termination certificate

Dos Santos et al., *At the Limits of Computability: Steps Towards
Autonomous AI Discovery on the Busy Beaver Frontier*, [arXiv:2610.22603v2](https://arxiv.org/abs/2610.22603v2), 24 pages. The source reports eleven
new machine accelerations and fifteen rediscoveries/reformulations. Its
new Cryptid classification is provisional; the associated halting problems
remain unresolved. An acceleration is not a proof that its trajectory
terminates. We did not replay all machine proofs or their Lean builds.

The useful source mechanism is the explicit preservation of the machine
and starting configuration while refining the macro-rule grammar. The
paper's Section 3 allows correcting rules and adding guards, but requires
native transition correctness and exhaustive closure. Appendix D preserves
the head marker while compressing the two sides of the tape separately;
it uses a one-block lookahead to improve a compression boundary. Expanding
the chosen blocks recovers the original tape. Abstracting away the counts
does not, and must be treated as a different map.

The most relevant concrete merger is Appendix A.11:

    C(a,b,c,d)=0^inf 1^(a+1) C> 00 1^(b+4) (01)^c 1^d 011 0^inf,
    d in {3,6},
    C(a,b,0,6)=C(a,b+3,0,3).                         (5)

**PROVED literal identity:** when c=0 the adjacent runs of ones merge;
both right tapes have a run of b+10 ones. The head and left tape are
unchanged. Choose the second representation as canonical when c=0,d=6.
This is exact source ownership, not an arbitrary quotient by a statistic.

`check_native.py` independently parses the raw six-state machine

    1RB1LA_1LC0RE_1LF1LD_1RE0LC_1LB1RA_---0LA

and checks all 576 parameter sources a=0..5,b=0..7,c=0..5,d in {3,6}.
No blank-orbit reachability filter is used. Every stated simplified
macro branch reaches the literal target tape or halts; the largest
checked native clock is 210. The blank tape reaches C(4,0,0,3) after
exactly 162 steps with head displacement -4. The report states its halt
clock convention explicitly. This is **FINITE-EXACT**, not a proof of
all parameters or a resolution of the machine's halting problem.

**Printed-boundary discrepancy.** Page 5 labels the Antihydra halt-family
arrow by 3b^2+9b-1 native steps and declares b>=0. At b=0 this is -1,
so that literal universal clock statement is false. The rendered page
was inspected; this is not an OCR superscript mistake. The source's
later discussion explicitly permits repairing bounds during formalization.
We do not infer that its separate Lean files contain the same error.
The strongest supported conclusion here is that this displayed formula
needs a boundary exception or corrected clock; we have not established
its all-b corrected formula. It is not used in our native test of A.11.

## 3. A precise common construction, and an important difference

Incoming [valuation work](../unprotected_ports_cancellation_20261025/valuation/README.md)
retains the integral Frobenius carry

    C_p(x,y)=((x+y)^p-x^p-y^p)/p,

and its associative cocycle identity before reducing coefficients. Its
cyclic and elementary abelian orbit polynomials have different modular
root depths even when their ordinary derivatives agree. Hasse information
is an indispensable extra coordinate there.

Our root-algebra multiplication retains the ordinary digit carry

    alpha_m(a,b)=floor((a+b)/m)

for representatives 0,...,m-1. Associativity gives

    alpha(a,b)+alpha((a+b) mod m,c)
       =alpha(b,c)+alpha(a,(b+c) mod m).

Both arise by dividing an exactly divisible additive defect, and both
retain information lost after quotienting. Their coefficient modules
and arguments differ; no isomorphism between the two cocycles is claimed.
For the root algebra, the cocycle determines the missing monomial ideal
and yields the all-size conductor theorem C1. For the incoming tournament
polynomial, the retained integral carry determines exceptional coefficient
characteristics. That is a checked common method with distinct targets.

| Source -> target | Map and preserved predicate | Lost data / required record |
| --- | --- | --- |
| Root cover -> two quotients | Character restrictions; full product algebra | At r vanishing roots, r quadratic sheet classes and conductor valuations |
| Universal root algebra -> exponent monoid | Record every nonnegative monomial and carry | Additive line specialization and unit twists must be supplied separately |
| Tape parameters -> native configuration | Expand literal blocks with the marked head | Empty-block descriptions may duplicate the same source; canonical ownership repairs counts |
| Collatz history -> affine bin | Retain M,D,A,beta and exact word | Owner multiplicity, admissibility and integer clock are needed for physical source counts |
| Integral orbit polynomial -> derivative mod p | Ordinary differentiation | Positive p-multiple-degree coefficients; Hasse data recover the missing coordinates |

The root normalization separates *distinct* geometric sources identified
by a quotient. Tape canonicalization identifies *duplicate descriptions*
of the same source. These operations point in opposite directions. What
they share is the need to compute the exact fibre before transporting an
invariant or a weight. Calling both “compression” without this distinction
would hide the relevant mathematical obligation.

## 4. Productive next obligations

The norm-seven chart problem now has an explicit integral answer. Its
arithmetic extension is to compute the residual squareclass map for fixed
original twists, rather than average over twists or assume that a rational
product point lifts. Over a normal non-UFD base, determine which conductor
ideal classes replace the principal monomial in C1.

For the Collatz lane, the remaining target is a fixed-source positive
prefix score exceeding a proved tail budget while maintaining canonical
history ownership. The histogram bound and a finite state grammar alone
do not give that inequality. The native-machine example reinforces this:
even an exact closed acceleration may leave the reachability question open.
