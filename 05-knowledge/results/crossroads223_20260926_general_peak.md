# Every fixed odd multiplier: elementary sharp peak asymptotics

**Status: PROVED, with the new connector lemma independently audited by
the root and flow lanes.** This integrates the automata lane's rational
connectors with a hard-wall upper bound. It refines THM-4480's cited
small-deviation argument to a logarithmic error. Private pairing price
is handled only for q=3 in the companion theorem; no general-q private
or globally consistent pairing claim is made.

## 1. Exact statement

Fix an odd integer q>=3 and set

    c=log_q 2, d=1-c, lambda=d/c,
    D=2^(-(1-H_2(c))), A=pi^2*cd/2,
    z=q/max(1,lambda)>1,
    kappa=(3/2)(pi^2*cd)^(1/3)(log z)^(2/3).

For length-L binary words let S_j=e_j-cj, Bad={S_j>0 for 1<=j<=L},
and P=max_(0<=j<L)S_j. Then

    rho_peak,q(L)=2^(-L)sum_Bad q^(-P)
                 =D^L E_c[lambda^(S_L)q^(-P);Bad].       (1)

For constants depending on the fixed q,

    log rho_peak,q(L)=L log D-kappa L^(1/3)+O_q(log L).  (2)

The arbitrary-edit price epsilon_L(q) of
[THM-4480](../../01-canon/theorems/THM-4480-peak-discounted-provability-price.md)
has the same formula, by its polynomial capacity sandwich. The constants
are not claimed uniform as q grows with L. The endpoint shift at q=5
is retained: when lambda>1 the lower-bound paths must finish near the top.

## 2. A finite hard-wall upper bound for every c in (0,1)

For R>=1 kill the Bernoulli-c walk on arrival outside (0,R+1), allowing
state zero only at time zero. Put H=R+3, t=pi/H and

    beta=log[d sin(ct)/(c sin(dt))],
    f(x)=exp(beta*x)sin(t(x+1)),
    mu=d exp(-c beta)sin(t)/sin(dt).

The unrestricted kernel obeys B_full f=mu*f. Every killed child reachable
from (0,R+1), or from the initial state zero, lies in (-1,R+2), where f
is nonnegative. Therefore the killed kernel satisfies Bf<=mu*f.

Using 1-u<=sin(sqrt(6u))/sqrt(6u)<=1 and
-log(1-u)<=2u for 0<=u<=1/2 gives |beta|<=t^2/3. Consequently
(R+1)|beta|<1. On 0<=x<R+1 both distances of x+1 to the endpoints of
(0,H) are at least one, so f(x)>=exp(-1)f(0). Positivity then gives

    Q_c(S_j in (0,R+1), 1<=j<=L) <= exp(1)mu^L
       <= exp(1)exp[-A L/(R+3)^2].                    (3)

The last inequality follows from the same sine product as the q=3 proof:

    log mu=-sum_(k>=1) a_k[1-c^(2k+1)-d^(2k+1)]t^(2k)
           <=-cd*t^2/2.

All brackets are positive for every 0<c<1. This is an unreflected
hard-wall operator; it does not use the reflected top mass 2d.

Partition P into bins m<=P<m+1. Since S_L<=P+d,

    lambda^(S_L)q^(-P) <= C_0 z^(-m),
    C_0=max(1,lambda^d).

Such a word survives in (0,m+2), including its last state. Apply (3)
with R=m+1. There are at most L+1 bins, as P< L, and

    m log z + AL/(m+4)^2
      = (m+4)log z+AL/(m+4)^2-4log z
      >= kappa L^(1/3)-4log z.

Thus the following explicit upper bound holds for L>=1:

    rho_peak,q(L) <= exp(1)C_0 z^4 (L+1)
                         D^L exp[-kappa L^(1/3)].      (4)

## 3. Polynomial-cost endpoint connectors supply the lower bound

The full proof is in
[the independent audit note, section 5](crossroads223_20260926_bridge_audit.md).
Here are its exact objects and information costs.

For M of order L^(1/3), take h=floor(M/8), K=floor(M^2/log M).
A rationally centered K-step bridge in a strip of width h has probability
at least mu_(p,h)^K/[K(h-1)]. Rotating at a minimum or a maximum makes it
one-sided. To connect two heights a distance O(M) apart, change its mean
from c to p=c+O(M/K). Under Q_c the probability loses precisely
exp[-K KL(p||c)], at worst exp[-O_q(M^2/K)]=M^(-O_q(1)). Confinement
costs exp[-O(K/h^2)], and the rotation costs a polynomial factor. Each
half-connector therefore costs only M^(-O_q(1)).

Two such halves through a point near M/2 connect arbitrary guarded
endpoints in (0,M), respecting their dyadic-time phase. Reverse time
and reflect space for the second half; the increments bit-c and their
weights are preserved. A connector starting at zero uses positive drift
and is strictly positive at every later time.

Use two 2K-step connectors around a long N=L-4K trace bridge of width
W=M-4. Set p=ceil(cN)/N. All central diagonal paths have the same number
of ones, so the change-of-measure cost is exp[-N KL(p||c)]=exp[-O_q(1/N)].
The first connector's ending height determines its trace index; hence
full word sets for different indices are disjoint. The trace lower bound
can be summed without an extra unproved mixing or entrywise claim.

The resulting path probability is at least

    Q_c(0<S_j<M for 1<=j<=L; S_L=r)
       >= M^(-C_q)exp[-A L/M^2-O_q(1)],                (5)

with r in (M-2,M-1) if lambda>1, and r in (1,2) otherwise. Such endpoints
exist in the correct phase because cL is irrational for odd q. Changing
width M-4 to M costs O_q(L/M^3)=O_q(1); the connectors' KL loss is only
O_q(log M). The exact integer phase equations, strict inequalities,
and product measure argument are proved in the linked lemma.

Insert (5) into (1). On this event the weight is at least a fixed constant
times z^(-M). Choose M=ceil((2A/log z)^(1/3)L^(1/3)). This gives

    rho_peak,q(L)>=D^L exp[-kappa L^(1/3)-O_q(log L)],   (6)

which with (4) proves (2). THM-4480's
rho_peak,q/M_L^*<=epsilon_L(q)<=rho_peak,q, M_L^*=O_q(L^3), transfers it
to the arbitrary-edit price.

## 4. Why the new operation matters

A forced suffix of length O(M) costs exp[-O(M)] and changes the sharp
constant. The successful connector uses MORE time, K of order M^2/log M,
but makes a smaller change in the mean and loses only a polynomial in M.
This is a precise scale tradeoff, not a temporal-randomness assumption
about one integer orbit. Initial-height phase, endpoint reward, killed
boundaries, and word multiplicities are retained throughout.

Classical sine identity: [DLMF 4.22.1](https://dlmf.nist.gov/4.22.E1).
The transfer, connector and finite hard-wall proofs above replace the
Mogul'skii input for this specific two-valued walk. No claim is made to
replace the general small-deviation theorem, or to establish convergence
for q=3 or divergence for q>=5.
