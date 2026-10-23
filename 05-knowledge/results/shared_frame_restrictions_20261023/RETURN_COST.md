# A guarded return needs a tail bound as well as an endpoint gap

23 October 2026. **PROVED scoped**, finite reversible positive-semidefinite
chains. This applies to the cut-preserving block coupling component by
component. It is not an independent estimate of native move congestion,
Collatz charge, or acceptance density.

## 1. Inherited return theorem and the extension

The incoming [guarded-return proof](../hamiltonian_observer_restrictions_20261023/trace/TRACE_WITH_CLOCKS.md)
shows that for a finite irreducible reversible PSD chain Q of gap at least
gamma>0, watching successive returns to any nonempty set S preserves a
gap at least gamma in RETURN-EVENT time. It retains the exact duration
generating matrix and stationary mean 1/pi(S). This is the standard trace
construction, with its original clock explicitly retained.

We add a source-indexed exponential excursion estimate. Write
alpha=pi(S), C=X minus S and B=Q_CC. For proper S,

    ||B||_(L2(pi|C)) <= rho := 1-gamma*alpha < 1.          (1)

The measure in this operator norm may be left unnormalized or normalized
on C; the norm is unchanged. For every actual starting state s in S and
integer t>=1,

    P_s(tau_S^+ > t)
       <= min(1, sqrt((1-alpha)/pi(s))*rho^(t-1)).         (2)

Consequently one sufficient source-specific expectation bound is

    E_s tau_S^+ <= 1 + sqrt((1-alpha)/pi(s))/(gamma*alpha). (3)

For S=X, the return is one step, with no killed operator. The tail is zero.
The deliberately simple bound (3) need not be sharp; summing the truncated
minimum in (2) improves it substantially. Let

    A=sqrt((1-alpha)/pi(s)),
    L=max(0, ceil(log(A)/[-log(rho)])).

For proper S, 0<rho<1, and

    E_s tau_S^+ <= 1+L+1/(gamma*alpha)
               <= 2+[1+log_+(A)]/(gamma*alpha).          (4)

Here log_+(A)=max(0,log(A)). Thus the cost of a small source weight is
logarithmic in this sufficient bound; the inverse acceptance-mass cost
remains explicit. No independence of source and accepted set is assumed.

## 2. Proof, retaining the source weight

Extend f:C->R by zero on S. Cauchy--Schwarz gives
(E_pi f)^2 <= (1-alpha)||f||_pi^2, hence

    Var_pi(f) >= alpha ||f||_pi^2.

The ambient gap inequality implies

    <f,Bf>_pi <= (1-gamma*alpha)||f||_pi^2.

The principal block B is PSD, so this upper Rayleigh bound is also its
operator norm bound. That proves (1). In fact the killed-kernel estimate
extends beyond PSD: for a nonnegative reversible killed matrix, Perron
theory bounds every eigenvalue's absolute value by its top eigenvalue.
PSD is retained here to match the inherited positive trace-chain theorem.

For t>=1 the exact source tail is Q_sC B^(t-1) 1_C. Regard the row Q_sC
as density h_s(y)=Q(s,y)/pi(y) on C. Reversibility gives

    ||h_s||_pi^2
       = sum_(y in C) Q(s,y)^2/pi(y)
       = [1/pi(s)] sum_(y in C) Q(s,y) Q(y,s)
       <= 1/pi(s).

Since ||1_C||_pi=sqrt(1-alpha), Cauchy--Schwarz and (1) prove (2).
Summing its geometric upper bound proves (3). This uses the weight of
the ACTUAL initial state; averaging it away would change the claim.
To prove (4), sum the first L terms using the bound one and the remaining
terms using A*rho^j. Since A*rho^L<=1, the latter sum is at most
1/(1-rho). Finally -log(rho)>=gamma*alpha and ceil(x)<=x+1
give the displayed logarithmic estimate. The case A<=1 uses L=0.

For comparison, under a stationary conditional start the mean is exactly
1/alpha, as proved in the incoming note. The sourcewise bound does not
replace that identity. For the complete-refresh chain Q=J_pi with gap 1,
B has top eigenvalue 1-alpha, so (1)'s dependence on alpha is sharp.
No positive uniform excursion bound follows as alpha tends to zero.

## 3. Combining this with cut-preserving restrictions

Lift the operator of RESTRICTIONS.md to all words as the sandwich
Q=P_G^H P_T^H P_G^H. Its nonzero spectrum is unchanged; the extra modes
are zero. For contiguous blocks, each nontrivial support component has
gap at least 8/9. Condition its stationary measure on that component
BEFORE computing alpha and pi(s).

Any actual source/port/clock predicate defines an accepted subset S of
that component. If nonempty, the return theorem and (1)–(3) give:

    return-event gap >= 8/9,
    killed excursion contraction <= 1-(8/9)*alpha.

This remains true if the predicate is not constant on graph classes,
because it is applied on the full word lift. It does not join different
common-information components. If the block choice gives a singleton
component, the mixing statement has no nonconstant content; choosing all
singleton blocks is not a proof of exploration or coverage.

The result is a precise combined design rule: preserve the observer cuts
that certify contraction, then account separately for acceptance density,
initial source weight, and integer duration. A meaningful optimization
must also retain required reachable endpoints; reducing the state space
until there is nothing to explore would make the certificate vacuous.

## 4. Remaining native and Collatz obligations

The intermediate ambient states of a return excursion can violate source
or port predicates. To obtain a native theorem one must realize the return
by allowed moves and control its integer length and congestion. Nor does
an abstract alpha lower bound prove universal Collatz entry: S must first
be an authenticated nonempty source-specific accepted set with the actual
charge/clock law. The missing independent estimate has become explicit,
not solved by renaming event time as physical time.

`return_cost.py` checks the PSD killed-kernel inequality and squared version
of every source-tail bound through time 12 for every nonempty subset of
the literal eight-word order-three coupling. It also checks all subsets
of a nonuniform three-state reversible chain, and the complete-refresh
sharpness control. It also solves the killed linear system independently
for every source expectation and tests both the square-root and capped-tail
bounds. The integer L is found by rational squaring, without approximating
logs. These are exact rational certificates, not numerical eigenvalue fits.
