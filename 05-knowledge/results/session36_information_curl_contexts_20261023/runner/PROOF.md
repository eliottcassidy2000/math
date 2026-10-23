# Runner-comb observations: an exact blind ladder and a sharp interval sidecar

2026-10-23, Session36. **PROOF-COMPLETE SCOPED THEOREM + FINITE-EXACT CHECKS.** This is an observer theorem for circle carriers, not a new LRC row exclusion or a claim of historical priority. All interval measures are normalized Haar measure on T=R/Z. No Git mutation was performed.

## 1. Recovery, corrected near miss, and delta

The read-only recovery began with `00-navigation/CURRENT-FRONTIER.md` and `01-canon/ACTIVE-GUARDRAILS.md`, especially LRC controls 12, 23, 26, 27d, 27i and tournament control 41.

Closest proved mechanisms:

- THM-2072, `fixed-owner-clock-bank-no-go-and-half-shift-certificate`: a common-multiple body makes every prescribed rational packet empty even when a continuous safe pair is explicit. Sensor blindness is not geometric difficulty.
- THM-2161, `fixed-modulus-banks-miss-primitive-defect-seven-lrc-rows`: labelled residues can agree at every prescribed finite modulus while an easy adaptive phase separates their actual gaps.
- THM-4349, `adaptive-owner-clock-completeness-and-twelve-height-tail-cutoff`, and THM-4371, `minimal-wall-lcm-adaptive-owner-clock`: retain an actual safe arc, strict tail danger, and physical walls; adaptive height-dependent clocks are complete on the dyadic seam. These do not assert that every seam is safe.
- THM-4451, `lrc14-dyadic-three-tail-strict-component-width-caps`: an a.e. interval merge at the removed endpoint 1/14 for (1,7,13) gives a false physical component. Measure and topology must stay separate.
- THM-2408, `endpoint-prony-resultant-clock-separation-and-shared-node-boundary`: common-gauge jumps produce finite exponential spectra; shared endpoint nodes can cancel. A Fourier lift is not a physical terminal clock.
- THM-2545, `word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile`, and THM-3731, `valuation-owner-equivariance-and-semantic-packet-boundary`: owner labels, complete margins, and all colours do not provide the joint common-ancestry diagonal.
- `05-knowledge/results/tournament_future_boundary_20261020/README.md`, Proposition 1: the correct compressed state depends exactly on which future operations are still allowed to query it.

The new result in this recovery is the explicit kernel of a finite-depth runner-comb response language and the sharp component-complexity threshold that kills it. A bounded targeted search did not find this exact statement in canon. This is not a literature-priority claim.

## 2. Exactly which observer is being studied

For an integer b>=2 and positive integer n, let

    D_n^(b) = {t in T : ||nt|| < 1/(2b)}.

For f in L1(T), define the phase-resolved overlap response

    O_n(f)(s) = integral_T f(t) 1_(D_n^(b))(t-s) dt.

The LRC(14) danger radius is 1/14, so b=7. The probe is a translated full n-tooth comb, not a single chosen tooth. Its endpoints are strict, but endpoint choices do not change this Haar integral.

Fix J>=1 and q=b^J. The depth-J observer is

    R_J(f) = (O_1(f), O_b(f), ..., O_(b^(J-1))(f)).

Every coordinate records the whole phase-response function of s. Thus this is a finite number of probe speeds, not a finite number of scalar phase samples. No claim is made that these independent translated probes are already permitted on a current LRC owner/word packet.

## 3. Exact blind-ladder theorem

For f,g in L1(T), put h=f-g. Then

    R_J(f)=R_J(g)

if and only if

    integral h=0, and h(t+1/q)=h(t) almost everywhere.       (1)

Equivalently the observation kernel consists exactly of zero-mean functions pulled back by the q-fold covering t -> qt.

Proof. With Fourier convention fhat(k)=integral f(t) exp(-2pi i k t) dt,

    (1_D_n)hat(0)=1/b;
    (1_D_n)hat(k)=0                         if n does not divide k;
    (1_D_n)hat(nm)=sin(pi m/b)/(pi m)        for m!=0.        (2)

To see (2), average the n translates of one interval of length 1/(bn). The root sum vanishes unless n divides k; in that case the n identical terms give the displayed sinc. Thus a nonzero frequency k is detected by probe n=b^j exactly when b^j divides k but b^(j+1) does not divide k. Across j=0,...,J-1 this is precisely the set q does not divide k. The mean is detected by every probe. All responses vanish iff hhat(0)=0 and hhat(k)=0 for every k not divisible by q.

These Fourier conditions are equivalent to 1/q-periodicity: translation by 1/q multiplies the kth coefficient by exp(2pi i k/q), and uniqueness of Fourier coefficients in L1 finishes both directions. The convolution responses are continuous, so a.e. equality of responses is also equality at every phase. QED.

The result holds for composite b as well: the divisibility bands b^j Z minus b^(j+1) Z still partition Z minus qZ. No prime-adic valuation theorem is being smuggled in.

### For prime b, valuation depth is the entire probe language

Suppose b=p is prime. The ladder R_J has exactly the same distinguishing power as the possibly infinite bank of ALL probe speeds n with v_p(n)<J. Indeed every such n kills the blind space in (1): if q divides k and n divides k, then p divides k/n, so the last line of (2) vanishes. The reverse implication holds because the bank includes n=1,p,...,p^(J-1).

Thus adding arbitrarily many or arbitrarily large speeds at already retained p-adic levels gives no new distinctions once the ladder is present. Any finite speed bank is blind on some q-periodic hostile: take J=1+max_n v_p(n), q=p^J, and use (5). For LRC radius 1/14, the relevant prime is 7. This statement concerns phase-overlap probes, not residue-clock owner relations; the two observer maps must not be identified.

For composite b, the prefix theorem remains valid, but this valuation-only extension is not asserted: divisibility by a composite base does not have the required additive prime valuation law.

## 4. A sharp component-complexity repair

A regular-closed interval carrier is the closure of a finite union of nondegenerate open arcs. The empty set and whole circle are allowed separately. Let c(K) be the number of its proper positive-length connected components (set c(empty)=0 and c(T)=1).

**Theorem.** If K is such a carrier with c(K)<q, and 0<=g<=1 is measurable, then

    R_J(1_K)=R_J(g)  implies  g=1_K almost everywhere.       (3)

Consequently, for two regular-closed interval carriers, if either has fewer than q components, equal responses imply exact set equality. Conversely there are distinct disjoint regular-closed carriers K_q,L_q, each with exactly q components, having identical responses. Thus on the class c<=C, depth J is injective exactly when b^J>C.

Proof of uniqueness. Write the c proper arcs of K, after choosing a cut in a gap, as [a_i,b_i], with all 2c endpoints distinct and cyclically ordered. Define

    P(t) = - product_(i=1)^c sin(pi(t-a_i)) sin(pi(t-b_i)).  (4)

It is a real, 1-periodic trigonometric polynomial of degree c. It is positive on the interiors of K, negative on the complement of K, and zero just at the endpoints. The sign assertions hold because each paired factor is negative precisely between its two endpoints, and the arcs are disjoint.

By (1), every Fourier coefficient of h=1_K-g of absolute index <=c is zero, since c<q; this includes the mean. Hence integral hP=0. But 0<=g<=1 makes hP nonnegative everywhere, and strictly positive wherever h is nonzero away from the finite endpoint set. Therefore h=0 a.e. Empty and full K follow just from the mean and 0<=g<=1. If K and L are regular closed, a.e. equality implies equality because each is the closure of its interior. QED.

This also gives the independently useful finite-moment statement: the mean and first c complex Fourier coefficients of 1_K determine it among all functions in [0,1]. These are 2c+1 real quantities. This is a sufficient summary, not a claim that its coordinate count is minimal.

Proof of sharpness. Define, with endpoints understood modulo 1,

    K_q = union_(r=0)^(q-1) [r/q, (r+1/4)/q],
    L_q = K_q + 1/(2q).                                  (5)

Both have measure 1/4 and exactly q positive-length components; they are disjoint. Their difference is 1/q-periodic and has mean zero, so (1) applies. Indeed every retained response is identically 1/(4b). The next probe n=q distinguishes them at phase zero:

    O_q(1_Kq)(0)=1/(2b),   O_q(1_Lq)(0)=0.               (6)

For b=2, the first value is 1/4 and the second remains zero because only endpoints touch. Formula (6) respects strict endpoints. QED.

## 5. What the finite spectral summary uses

When c(K)<=C<q, for every 1<=k<=C there is a unique j<J with b^j dividing k and b^(j+1) not dividing k. Write k=b^j m. Then

    1_Khat(k) = O_(b^j)(1_K)hat(k) / [sin(pi m/b)/(pi m)],
    measure(K) = b * integral_T O_1(1_K)(s) ds.           (7)

The denominator is nonzero. The Fourier coefficients on the right are integrals across the whole phase-response function, not values at finitely many arbitrary phases. Thus a finite spectral summary can be sufficient even though a finite phase-sample bank is not.

### Quantitative bounded-component stability

Let K have c>=1 proper arcs, let 0<=g<=1, and set

    epsilon = ||1_K-g||_1,
    eta = max_(|k|<=c) |(1_K-g)hat(k)|.

Then

    eta >= epsilon^(2c+1) / [2(4c)^(2c)].                  (7a)

Thus bounded component complexity supplies an explicit inverse Holder estimate, not merely injectivity. No optimality of its exponent or constant is claimed.

Proof. The polynomial P in (4) has Fourier-coefficient l1 norm at most 1: each sine factor has two half-integer-frequency coefficients of modulus 1/2, and after 2c factors the frequencies are integers. Consequently integral hP<=eta. Put r=epsilon/(8c). The union of radius-r neighborhoods of the 2c endpoints has measure at most 4cr=epsilon/2. Since |h|<=1, at least epsilon/2 of the discrepancy mass remains outside. Each sine factor there has absolute value at least 2r, because |sin(pi x)|>=2 dist(x,Z). The sign match hP=|h||P| yields integral hP>=(epsilon/2)(2r)^(2c), proving (7a). The case epsilon=0 is immediate. Empty and full K are controlled directly by the mean.

A sharper optimization of the same argument gives

    ||1_K-g||_1 <= min(1, (2c+1) eta^(1/(2c+1))).         (7b)

For any r>0, the endpoint neighborhoods contribute at most 4cr discrepancy, and on their complement the sign-polynomial bound gives at most eta/(2r)^(2c). Minimize

    4cr + eta/(2r)^(2c)

at r=eta^(1/(2c+1))/2. Here eta<=1, so r<=1/2. If eta=0, let r decrease to 0. This proves (7b); it improves the constant in (7a), not its exponent.

For a uniform error bound on the entire retained response functions,

    rho = max_(0<=j<J) ||O_(b^j)(1_K)-O_(b^j)(g)||_infinity,
    c<q,

formula (7) gives

    eta <= A_(b,c) rho,
    A_(b,c)=max(b, pi*c/sin(pi/b)).                       (7c)

Indeed every Fourier coefficient of a response error is at most rho. For k=b^j m with b not dividing m, |sin(pi m/b)|>=sin(pi/b) and |m|<=c, while the mean is multiplied by 1/b. Therefore

    ||1_K-g||_1 <= min(1,(2c+1)(A_(b,c)rho)^(1/(2c+1))).  (7d)

All integer bases b>=2 are allowed in (7c)-(7d). These are L1 stability bounds, not endpoint/Hausdorff stability, null-point detection, or owner/arrival transport. Their constants retain the measurement bandwidth and component budget. No optimal uniform-response noise exponent is claimed.

## 6. Finite scalar phase samples still fail with one component

Fix ANY finite list of translated combs (n_i,s_i), with arbitrary positive integer speeds and phases, and observe only

    integral 1_K(t) 1_(s_i+D_ni)(t) dt.

Their finitely many walls partition the circle into finitely many open cells. Choose one nonempty cell U. Every probe indicator is constant on U. Choose two distinct disjoint closed subintervals K,L of U with the same positive length. Then all the observations agree: each is either that common length or zero. K and L have one component each and different locations. Their means agree too.

Adding any fixed finite point-membership samples does not help: include those points in the wall partition and choose the same cell. This theorem concerns a prescribed finite observation bank. It does not rule out adaptive algorithms with extra feature-size assumptions, a known endpoint grid, or instance-dependent stopping time.

### Forgetting phase loses a one-interval carrier even with every speed

If L=K+a, then O_n(1_L)(s)=O_n(1_K)(s-a). Thus every Fourier magnitude |O_nhat(k)| and every translation-invariant norm of the response is identical for K and L, at ALL probe speeds. Two distinct translates of one short closed arc already witness the loss. The component theorem requires the common labelled phase s; a magnitude spectrum or a phase-free multiset of response values does not supply it.

## 7. Haar recovery stops at the boundary skeleton

Every overlap integral ignores null sets. For a closed finite union F of arcs and isolated points, let

    F_reg = closure(interior F),
    I(F) = F minus F_reg.

The set I(F) is finite, F=F_reg disjoint union I(F), and all responses of F and F_reg agree. The component theorem can recover F_reg under its complexity hypothesis, but I(F) is an independent necessary sidecar for exact closed-set reconstruction. Conversely, giving F_reg and I(F) recovers F exactly.

This is physically relevant rather than a merely artificial endpoint objection. For the tight AP V={1,...,13}, the weak-safe set at threshold 1/14 is exactly

    {a/14 : 1<=a<=13 and gcd(a,14)=1}
    = {1,3,5,9,11,13}/14.                                (8)

Proof. If all ||jt||>=1/14, then the fourteen points 0,t,...,13t have every pairwise circular separation at least 1/14: their index difference is one of 1,...,13. Fourteen circular gaps therefore all equal 1/14. Hence t=a/14 with gcd(a,14)=1. Conversely such a t permutes the fourteen grid points and every required distance is at least 1/14. QED.

Thus an actual nonempty LRC weak-safe carrier has every Haar overlap zero. Comparing its observer with the empty set does NOT assert that an empty safe carrier is physically realized by some 13-speed row. It demonstrates why a mass-only observer cannot establish the existential boundary case. Open or mixed-endpoint carriers require additional endpoint membership bits as well; THM-4451's puncture is exactly such a warning.

## 8. Owner/arrival and tournament connection, with the loss ledger

Source: a circle carrier with its actual phase coordinate. Target: its responses to the specified translated-comb language. Map: convolution by those comb indicators. Preserved: all listed future overlap experiments. Lost: the exact zero-mean q-periodic difference, and all null-set topology. Restoration: either one carrier's component bound c<q plus its isolated-point sidecar, or more probe speeds together with boundary data. Cheapest hostiles: (5), two tiny equal arcs inside one finite probe cell, and the six-point AP carrier (8).

The precise connection to tournament continuation work is equality under a specified test language, not a tournament structure on runner speeds. In the tournament result, adding one permitted cycle witnesses inequality of a retained response coordinate; dropping that permitted cycle makes the coordinate unobservable. Here adding the next comb witnesses the periodic hostile by (6); retaining only the previous probes makes its q-periodic difference unobservable. Both quotient maps are exact relative to their operation language.

There is no map carrying LRC cover rows to tournament orientations, and no transfer of Hamiltonian-path positivity, Burnside weights, or tournament enumeration follows. In particular geometric set recovery does not identify the actual owner, terminal word, common ancestry, or a later target-active role. THM-2545's aligned/swap tables can retain all one-point geometric observations while changing the common-base arrival diagonal. Those sidecars remain separate.

## 9. Status and next decisive test

The abstract theorem applies only after a genuine carrier and lawful phase-response bank have been constructed. It suggests a concrete experiment on an already extracted interval packet: count its positive-length components, retain isolated points, and ask whether the available response language includes every needed b-adic band. It does not prove a uniform component bound or that independent translations are lawful for the current owner/word packet. It neither closes the dyadic seam nor any of the 165 valuation rows; LRC(14) remains OPEN.

The checker `check.py`, machine-readable `receipt.json`, byte-matched normal/optimized transcripts, `VALIDATION.md`, and exact `SOURCES.json` accompany this proof. The separate `NOISE.md` gives the optimized quantitative corollary. The 92,927 always-active exact checks include 152,064 binary carriers and 65,024 stability pairs. No theorem ID or canonical status has been changed.
