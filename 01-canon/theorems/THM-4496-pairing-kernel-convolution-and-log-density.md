---
id: THM-4496
title: "Finite cutoff convolution exactly determines the attained minimum logarithmic pairing density"
status: "PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT controls"
source: "crossroads10-20260926"
depends_on:
  - 01-canon/theorems/THM-4491-pairing-two-step-tree-extra-density.md
  - 01-canon/theorems/THM-4492-pairing-two-cutoff-density-separation.md
proofs:
  - 05-knowledge/results/crossroads10_20260926_flow.md
scripts:
  - 04-computation/experiments/crossroads10_20260926_flow.py
  - 04-computation/experiments/crossroads10_20260926_audit.py
---

# THM-4496: finite kernels and the global logarithmic cost

**PROVED + independently audited.** Reserved as an empty stub in pushed
checkpoint 3130458ce, then promoted after root proof audit and a separate
geometry-lane audit. This concerns the global TWO-STEP PAIRING model from
THM-4491, not convergence of the unmodified Collatz map.

Let F be a globally legal flip set, A_F(X) its count through X, r=2/3,
and d_F(X)=A_F(floor X)/X. For a finite probability vector p on the
nonnegative integers, put a_k=p_k r^(-k). THM-4492 gives the limit

    B(p)=lim_(X->infinity) min_F sum_k p_k d_F(r^k X).

Rational p suffices initially; a uniform total-variation estimate extends
the definition to real p. Write u_m for the uniform vector on 0,...,m-1
and B_*=sup_p B(p) over all finite probability vectors.

## Conclusions and mechanism

1. B is concave, invariant under nonnegative integer shifts, and
   |B(p)-B(q)|<=TV(p,q). The finite objectives are minima of linear
   functionals taking values in [0,1]; shifts change only bounded floor
   errors before the limit.
2. B(p*q)>=max(B(p),B(q)). Convolution is a convex combination of shifts;
   use concavity and shift invariance. In the original unnormalized
   variables, J_(a*b)(X)>=sum_j b_j J_a(floor(r^j X)).
3. Uniform adjacent windows are complete:

       B(u_m)>=B(p)-sum_k p_k min(k/m,1),
       lim_m B(u_m)=B_*.

   Compare p*u_m with u_m in total variation, then fix p before sending
   m to infinity. Also m B(u_m) is superadditive. This does not assert
   monotonicity at every successive m. Global realization requires the
   separate phase-repair argument below.
4. Every global pairing satisfies

       lower_log_density(F)>=B_*
           >=821510388809/2677850419968
           =0.30677978974599224... .

   Integrate the fixed-kernel inequality in log X; its finitely many
   shifts contribute bounded endpoints. Abel summation identifies the
   average with (sum_(i<=X,i in F)1/i)/log X, up to a vanishing term.
   The exact rational certificate is inherited from THM-4492 and replayed.
5. Define H(X)=min_F sum_(2<=i<=X,i in F)1/i. Then

       lim_(X->infinity) H(X)/log X=B_*,

   and one global legal pairing has an actual LOGARITHMIC density B_*.
   Consequently B_* is the attained minimum of lower and upper
   logarithmic density, and of logarithmic density where that limit
   exists. A legal prefix through Y>=2 can be extended through any X>Y
   at cost at most H(X)+log Y+3. Natural-density attainment is not claimed.

For assertion 5, finite subtree root-bit gaps obey |Delta_i|<=3/(i-1).
The even recurrence has one child; the odd recurrence has one nonpositive
and one nonnegative child term, so their absolute sum is bounded by the
maximum child bound. Every child has j-1>=3(i-1)/2. Fixing the prefix
only restricts roots Y<j<=floor((3Y+1)/2), whose total gap is at most9/4.
The prefix itself costs at most log Y. This first proves lower-logarithmic
attainment at h=liminf H(X)/log X; the next argument identifies h and
upgrades attainment to an actual logarithmic limit.

## The phase repair that proves sharpness

For weights 0<=w_i<=C/i above a cutoff Y>=2, deleting edges crossing Q
log_(3/2)-phase cells changes the finite optimum by at most18CQ/(Y-1),
uniformly in the final cutoff and any legal fixed prefix. Restore edges
bottom-up, reoptimizing the entire child component at cost at most its
root gap3C/(j-1). A crossing is within1/2 of a geometric phase boundary;
at most two edges charge each boundary, whose reciprocal sum is geometric.

Let D_m(X) be the free harmonic optimum on (floor(X/R^m),floor X],
R=3/2 and integer m>=2. Two such shells shifted by a factor R^t,
0<=t<=1, have same-weight optima differing by at most4C for sufficiently
large X: an extra band costs at most C and its interface repair at most3C.
The uniform-kernel interior weight, after multiplication by m, is

    c_t(i)/i - 2/(R^t X),
    c_t(i)=R^(1-frac(log_R(X/i)+t))/(R-1),  2<=c_t(i)<=3.

Dropping the nonnegative bottom cap and paying at most2 for the displayed
correction, then moving the shell, bounds its free c_t-weighted optimum
by m f_m(R^tX)+14, where f_m is the normalized cutoff objective converging
to B(u_m). This retains the actual phase weight.

Use Q equally spaced phase cells and shifts t=j/Q. The average over shifts
of the infimum of c_t on any one cell is exactly

    a_Q=1/[Q(R^(1/Q)-1)] -> 1/log R.

Writing D_q for the phase-cut harmonic minima, weighted optima dominate
their respective sums of cell infima times D_q. Phase repair gives
D_m(X)<=sum_q D_q+18Q/(floor(X/R^m)-1). Average the weighted inequalities;
send X to infinity with m,Q fixed, then Q to infinity. This yields

    limsup_X D_m(X)/(m log R)<=B(u_m)+14/m.

Concatenate shells, paying at most3 at every interface. Thus

    limsup_X H(X)/log X<=B(u_m)+14/m+3/(m log R).

Finally m tends to infinity. The already proved lower bound gives the
claimed full limit B_*. The order of these three limits is essential.

For attainment, the free suffix optimum E(Y,X) obeys
E(Y,X)<=H(X)-H(Y), while any fixed prefix extends at cost at most its old
cost+E(Y,X)+3. Choose conditional optima at X_n=2^(n^2). Their excess
over H(X_n) is at most3n+O(1). Monotonicity of harmonic cost and
log X_(n+1)/log X_n->1 squeeze the resulting single global pairing to
logarithmic density B_*.

## Failure boundaries and exact controls

Finite-depth truncations are not shift invariant: B_H(S^s p)=B_(H-s)(p),
with zero for negative depth. Their valid convolution inequality is
B_H(p*q)>=sum_j q_j B_(H-j)(p). A delayed kernel can have zero retained
certificate while its exact value is unchanged.

At cutoffs486,216,96 the two separate costs210 and94 become306 when one
common assignment is required. In the ten-vertex cutoff11 control, the
apparent gain from4 to5 is instead entirely a nested-floor discrepancy.

The compact phase-band family from the proof note has every finite
grid-kernel certificate sqrt(6)-2, but every member has logarithmic
density1/2. Thus the generic inference from certificate completeness to
global sharpness is false. The actual pairing tree instead supplies the
quantitative phase-repair property used above. The numerical value of B_*
beyond the stated lower bound, natural-density attainment, longer-horizon
pairing legality, and Collatz convergence remain **OPEN** here.

Full proofs, scope labels and controls are in the linked proof note.
