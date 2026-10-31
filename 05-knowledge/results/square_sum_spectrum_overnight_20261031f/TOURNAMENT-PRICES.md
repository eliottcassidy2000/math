# Fixed boundary decorations and the sharp price of a repeated binary word

**Status: PROVED in the stated native grammars, with exact controls. These are structural boundaries and explicit families, not universal spectrum completeness.**

The inherited mechanisms are packet e's translated-cluster identity, THM-4748 / `01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`, and the marked-block action in packet b. The earlier fixed-feedback/bounded-run density-zero result is in `05-knowledge/results/feedback_word_hamilton_spectrum_20261020/run_density.md`. A repeated mixed word has unbounded run count, so the result below is a different obstruction.

## 1. Fixed boundary decorations cannot repair the midpoint deficit gap

Let L be a fixed finite subset of{2,3,...}; let S be a nonempty finite subset of nonnegative integers with minimum0; let R be a fixed finite set of nonnegative integers containing0. Fix an integer D>=0. For all sufficiently large a, set K=2a+D and take the one-apex tournament whose chain1,...,K is transitive and whose apex loses exactly to

    L union (a+S) union (K-R).

Require the three regions to be separated by gaps at least2. This eventually holds for any such fixed data. The apex beats1 and losesK, so the tournament is strong. Put

    U=sum_(l in L)2^(-l),          V=sum_(r in R)2^(-r),
    A=sum_(s in S)2^s,             B=sum_(s in S)2^(-s),
    E(X)=sum_(x<y in X)2^max(y-x-1,1).

Its raw Hamiltonian-path count is exactly

    H(a)=C_2 2^(2a)+C_1 2^a+C_0,                       (1)
    C_2=2^(D-1)V(1-U),
    C_1=(1-U)A/2-2^(D-1)VB,
    C_0=1+sum_(l in L)2^(l-1)-E(L)-E(S)-E(R).

**Proof.** Apply the known one-apex singleton-minus-pair formula. Right singleton terms and left-to-right pairs contribute C_2. Middle singletons and left-to-middle/middle-to-right pairs give C_1. Remaining within-region pairs and left singletons give C_0. Cross gaps are at least2, so only within-region pairs need the special adjacent-gap coefficient. `three_region_fans.py` verifies2560 actual adjacency counts against this exact rational formula.

**Theorem 1.** If H(a)=2^k(a)-c for a fixed integer c and infinitely many unbounded a, then L is empty, R={0}, D=max S, S=D-S, and c=E(S)-1. Thus every such positive fixed deficit already belongs to packet e's symmetric midpoint grammar.

**Proof.** Since H(a)/2^(2a) tends to C_2>0, the powers2^(k(a)-2a) tend to C_2 along the infinite subsequence. A convergent sequence of powers of two with positive limit is eventually constant. Therefore C_2 itself is a power of two.

Now0<=U<1/2 and1<=V<2, so1/2<V(1-U)<2. If C_2 is a power of two, V(1-U)=1. Both factors are positive dyadic rationals. A product of two dyadic rationals is1 only when both reduced odd numerators are1; equivalently both factors are powers of two. Their displayed ranges force V=1 and1-U=1. Hence R={0}, L empty, and C_2=2^(D-1).

The remaining equation C_1 2^a+C_0=-c at infinitely many a forces C_1=0. Thus A=2^D B. Uniqueness of binary expansion gives S=D-S. Since min S=0, its maximum is D. Equation(1) now gives c=E(S)-1, as claimed.

In particular adding arbitrary fixed boundary clusters cannot produce the missing deficit9 family. This does not exclude individual counts2^k-9, nor moving/overlapping clusters or a second feedback vertex.

### Near equality and its exact cost

There is a useful near miss. If L={2} and R={0,2,...,2r}, then

    V(1-U)=1-2^(-2r-2).                                (2)

The factors approach exact cancellation through15/16,63/64,255/256,... but no finite r gives1. Allowing r to grow with a might turn the leading error into a bounded amount. It also makes the right region spread toward the other regions; the separation assumption must then be rechecked, and the within-region energy changes. Equation(2) is not permission to substitute r(a) into a fixed-region conclusion while dropping those terms. The periodic-word calculation below exposes one such surviving price.

### Stronger closure: any fixed finite collection of affine clusters

The same rigidity extends beyond three regions. Let the losing chain positions be a fixed finite union of sets

    j t + S_j,

where the slopes j are nonnegative integers, the offset sets S_j are fixed finite integer sets, and t ranges over all sufficiently large integers. Require the largest slope M to be positive; constant presentations are excluded. Normalize to a strong one-apex tournament: **the apex beats1 and loses to the last chain vertex K=M t+D**. Thus K belongs to the losing set, M is the largest slope and D=max S_M. At slope0 all positions, if any, lie in{2,3,...}.

**Affine-cluster rigidity theorem.** If the raw count equals2^k(t)-c for one fixed c at infinitely many unbounded t, then the top cluster consists only of K, the slope0 cluster is empty, and either there is no other cluster or there is exactly one cluster at slope M/2. In the latter case its offset set is symmetric under b->D-b. Hence these families reduce to the previous midpoint fans; adding finitely many other affine clusters cannot produce a new fixed deficit.

**Proof.** The count is a polynomial in2^t. Its leading coefficient is2^(D-1)V(1-U), where U is the slope0 dyadic sum and V=sum_(b in S_M)2^(b-D). The preceding dyadic-unit argument forces U=0 and V=1. Thus the top cluster is the singleton{D}.

For each remaining positive slope j write A_j=sum_(b in S_j)2^b and B_j=sum_(b in S_j)2^(-b), with zero for absent clusters. Every lower coefficient must vanish, so

    A_j = sum_(v-u=j) A_v B_u,       0<j<M.          (3a)

All summands are positive. The set J of occupied positive slopes is consequently closed under positive differences. If g=min J, repeated subtraction shows that every member is a multiple of g, and J={g,2g,...,M}.

Suppose M>=3g. For j=M-g the only pair is(g,M), giving A_(M-g)=A_M B_g. Since A_M=2^D, uniqueness of binary expansion forces S_(M-g)=D-S_g, and therefore B_(M-g)=A_g/A_M. But in the equation for j=g, the pair(M-g,M) already contributes exactly A_g. The distinct pair(g,2g) contributes a further strictly positive amount. This contradicts(3a). Thus J is either{M} or{M/2,M}. In the latter case(3a) directly forces the offset reflection. The constant term is1-E(S_(M/2)), proving the claimed reduction.

The coefficient census `affine_cluster_fans.py` checks every offset profile in its declared finite universe, then verifies the accepted profiles by literal adjacency counts. The minimal three-slope hostile, losing at t,2t,3t, has

    H=2^(3t-1)-2^(t-1)+1.

The extra negative exponential cannot be treated as a constant deficit. This theorem concerns fixed affine parametrizations, not arbitrary one-apex words. Rational slopes on a fixed arithmetic progression reduce to integer slopes by rescaling; unbounded offset sets, non-affine locations and a second feedback vertex remain available.

## 2. Exact repetition formula for a native singleton word

The marked state is(F,q,g,h). A singleton0 or1 acts by

    S0(F,q,g,h)=(F+q+g,2q,0,h),
    S1(F,q,g,h)=(F,2q+h,h,h).

Here F is the actual raw path count, h>0 is the internal block path count, and q,g retain the source boundary data. Let w=b_0...b_(L-1) be a fixed nonempty binary word, B_0=2^L. After one period, at a compatible boundary g=h b_(L-1), its action has the form

    q'=B_0 q+beta h,
    F'=F+alpha q+gamma h,
    beta=the binary integer b_0...b_(L-1),
    alpha=sum_(b_j=0)2^j.                              (3)

Define

    delta(w)=gamma-alpha*beta/(B_0-1).                  (4)

After r periods,

    q_r=B_0^r(q+h beta/(B_0-1))-h beta/(B_0-1),
    F_r=F + alpha(q+h beta/(B_0-1))(B_0^r-1)/(B_0-1)
            + h delta(w) r.                           (5)

An arbitrary initial g is handled by applying one actual period first; the compatible boundary then holds. Formula(5) is just summation of the geometric recurrence, with no fictitious negative state treated as a tournament.

## 3. Renewal proof of strict positivity

For a mixed word, let u_j be the rational number with binary expansion

    u_j=0.overline(b_j b_(j+1)...b_(j+L-1)),

with indices cyclic. It satisfies u_(j+1)=2u_j-b_j. Substituting the formal stationary value q_j=-h u_j into the increment formula computes the coefficient delta; this is an algebraic evaluation, not an admissible native source. It gives

    delta = sum_(b_j=0)(b_(j-1)-u_j).                 (6)

Split the zeros into cyclic runs. A zero-run of length r contributes one transition1->0. If u is the phase at the following first1, the sum of phases along that zero-run is u(1-2^(-r)). Its contribution to(6) is

    1-u(1-2^(-r))>0.                                 (7)

Thus delta>0 for every mixed word. Constant words have delta=0. This is a renewal price attached to the actual integer period.

## 4. Sharp minimum at every period length

Let R_0 be the number of cyclic zero-runs, equivalently1->0 transitions. Define

    C_t=#{j:b_j=0 and b_(j+t)=1},      1<=t<=L-1.

Then C_t=C_(L-t). Every0->1 displacement by t positions crosses some0->1 boundary, and each such boundary has at most t possible starting positions. Hence

    C_t<=R_0 min(t,L-t).                              (8)

Expanding the binary phases in(6) gives the exact energy

    delta=R_0 - sum_(t=1)^(L-1) C_t 2^(L-1-t)/(2^L-1).

Use(8) and sum the two geometric pieces. For every mixed word of length L>=2,

    delta(w)>=R_0 mu_L>=mu_L,
    mu_L=(2^floor(L/2)+2^ceil(L/2)-2)/(2^L-1).        (9)

**Equality classification.** Equality delta=mu_L forces R_0=1. Such a circular word has one zero-run of length r and one one-run of length s=L-r, and

    C_t=min(t,L-t,r,s).

It attains equality in every positive-weight term of(8) exactly when min(r,s)=floor(L/2). Therefore the minimizers are precisely rotations of0^r1^s with |r-s|<=1. This proves the sharp minimum and its equality cases at all lengths. For even L=2r, mu_L=2/(2^r+1).

The minimum tends to0. Fixed-period positivity therefore supplies no uniform positive price when the period itself may change. This is the explicit escape left for an adaptive controller.

**Density is insufficient.** The primitive length-six words000111 and001011 have the same density1/2, but prices2/9 and68/63. The full cyclic correlation row, not just the number of ones, determines the price. Here “run lengths as balanced as possible” means two contiguous runs; it does not mean a balanced Christoffel word with evenly dispersed symbols.

The rational circle coordinate u_j also appears in prior binary-angle/necklace work, including THM-4714 / `01-canon/theorems/THM-4714-collatz-parity-numbers-are-mandelbrot-external-angles-and-minkowski-farey-addresses.md`. The only shared map used here is the explicitly proved cyclic binary shift u->2u-b. No Collatz convergence, angle-landing, or arithmetic-orbit assertion from that file is a dependency. THM-4485 / `01-canon/theorems/THM-4485-periodic-edit-price-feedback-sets.md` concerns a different price, for modified parity dynamics; its density-threshold warnings motivate retaining the full word, but its numerical bounds are not imported.

## 5. What fixed native contexts can and cannot do

Fix an actual prefix state with h>0, repeat a mixed w, and append any fixed word of the singleton/D grammar. The suffix changes the output F by a positive multiple of the incoming F plus terms depending only on q,g,h. The multiplier is the product of its internal block H values. At period boundaries q is a pure exponential plus a constant, and g,h are constant. Consequently the final count has the form

    A (2^L)^r + C r + D,       A>0, C>0.             (10)

**Corollary.** Such a fixed-context family can meet2^k-c for a fixed c at only finitely many repetition counts r. Otherwise division by2^(Lr) makes2^(k-Lr) tend to A. Discreteness forces A=2^b and k=Lr+b eventually along the subsequence. The remaining condition Cr+D=-c has at most one solution, a contradiction.

This extends the earlier fixed pure-zero-tail obstruction to every mixed periodic singleton controller. It does not exclude nonperiodic words, changing periods, or controllers that change their actual source.

### An exact family displaying the price

The actual block D with mask11 gives(F,q,g,h)=(9,18,4,5). Append S0 to get(31,36,0,5), then repeat1100. Formula(5) gives

    H_r=32*16^r+2r-1=2^(4r+5)+2r-1.                (11)

The number of vertices is7+4r in the retained presentation, including its named padding vertex. Literal adjacency counts verify r=0,...,8; the all-r identity follows from the action. The extra2r is the exact unpaid term that a pure-exponential fit would omit.

## Reproducibility and next obligation

`periodic_drift.py` checks every one of32766 binary words of lengths1..14, exact renewal/correlation identities, sharp inequality and equality cases,2540 state repetitions,90 literal singleton-context tournaments and the nine D-block examples. `three_region_fans.py` checks2560 actual decorated-fan counts and64 leading-factor controls. Normal and optimized receipts agree.

`deficit_nine_queries.py` uses the inherited exact source inverse to construct changing-word witnesses for every2^k-9 with9<=k<=36 and checks each against actual tournament adjacency. This is a finite test of the surviving route, not an eventual theorem. It deliberately keeps the complete words and vertex clocks instead of fitting a fixed-period formula.

The next useful tournament direction is a genuinely changing-word source controller that both tunes the exponential coefficient and pays its boundary/renewal terms. The sharp value mu_L quantifies how cheap a long period can become, but does not construct that payment or certify universal entry. The seven-parent inverse bank remains the stronger numerical coverage instrument; universal raw-H completeness remains OPEN.
