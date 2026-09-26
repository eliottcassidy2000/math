# Rational bridges, positive transfer matrices, and the private Collatz price

**Status: PROVED + INDEPENDENTLY AUDITED by the flow and automata lanes.** This strengthens
the companion [Robin estimate](crossroads223_20260926_robin.md). Collatz,
the globally consistent pairing price, and the uniform one-unit Robin
comparison remain OPEN. Canon: THM-4488,
`01-canon/theorems/THM-4488-private-pairing-peak-price-rational-bridge.md`.

## 1. Objects and theorem

Use c=log_3(2), d=1-c, D=2^(-(1-H_2(c))), A=pi^2*c*d/2, and

    kappa=(3/2)(pi^2*c*d)^(1/3)(log 3)^(2/3).

For binary words let S_j=e_j-c*j, Bad_L={S_j>0 for 1<=j<=L},
and A_M(L)=#{w in Bad_L: S_j<M for 0<=j<L}. Define

    rho_peak(L)=2^(-L) sum_(w in Bad_L) 3^(-max_(0<=j<L) S_j).

N_m is the reflected-barrier count in the companion note. The private
price pi_L and the arbitrary-edit price epsilon_L have precisely the
definitions in the [pairpeak note](procgen_pairpeak_20260926_pairing_peak_price.md)
and [THM-4480](../../01-canon/theorems/THM-4480-peak-discounted-provability-price.md).
In particular, private certificates may differ for different sources.

The following hold:

1. For L>=3 and m>=2,

       N_m(L) <= 4 exp(1) (m+3)L A_(m+5)(L).                  (1)

2. For L>=3,

       2 rho_peak(L) <= pi_L
          <= [4 exp(1) 3^5 L(L+3)(27L+18)+4L 2^(-L)] rho_peak(L).
                                                                    (2)

3. For L>=8, with C=2(1+exp(1)3^6),

       D^L exp(-kappa L^(1/3)) / [108 exp(1) L^(4/3)]
          <= rho_peak(L)
          <= (C/2)(L+1)^2 D^L exp(-kappa L^(1/3)).            (3)

Consequently log rho_peak(L), log pi_L, and log epsilon_L all equal

       L log D - kappa L^(1/3) + O(log L).                   (4)

The q=3 conclusion is elementary: the matching lower bound no longer
requires Mogul'skii's small-deviation theorem. The polynomial comparison
in (2) is stronger than equality of stretched-exponential constants.

## 2. A finite positive-matrix lemma

If a square nonnegative matrix P has all 2-by-2 minors nonnegative
(TP2), then

       spectral_radius(P) <= trace(P).                      (5)

Proof: On any ordered subset of indices, uncrossing an inversion in a
permutation cannot decrease the product of its selected entries. Thus
the product around every simple directed cycle is at most the product
of the diagonal entries at its vertices. This uses multiplication only,
so zeros are harmless. Decompose a closed walk into simple cycles:
its edge product is at most the product of its visited diagonal entries,
with multiplicities. Summing over all n-tuples of visited vertices gives
trace(P^n)<=trace(P)^n. For a nonnegative finite matrix,
limsup_n trace(P^n)^(1/n)=spectral_radius(P), by the Perron decomposition
into irreducible blocks (take multiples of a block's period). This proves
(5). A generic nonnegative matrix need not satisfy it: [[0,1],[1,0]]
has trace zero and spectral radius one, and a negative 2-by-2 minor.

## 3. A rational bridge that retains the exact Collatz threshold

Set e=ceil(cL), p=e/L, and b=1-p. For L>=3, c<p<=4/5.
For L=3,4,5 this is direct; for L>=6 use c<19/30, which follows
from 2^30<3^19. Hence 0<p<1. Fix an integer M>=2.

At time j the states are (Z-pj) intersect (0,M). A step goes up by
1-p with probability p, or down by p with probability b, and is killed
outside (0,M). Since pL=e is an integer, the L-step matrix P starts and
ends on exactly the states 1,...,M-1.

Every two-by-two minor of P is nonnegative. Indeed, two paths joining
ordered initial states to reversed ordered endpoints must meet at an
integer time: their height difference is integer and changes by at most
one per step. Swap their tails at the first meeting. This injects the
reversed-endpoint pairs into the ordered-endpoint pairs, preserving
survival and total weight. This is the needed TP2 property; no claim
about infinite-dimensional eigenvectors is involved.

Put t=pi/M and

    beta_p=log[b sin(pt)/(p sin(bt))],
    f(x)=exp(beta_p*x) sin(tx),
    mu_p=b exp(-p beta_p) sin(t)/sin(bt).

The unrestricted one-step operator sends f to mu_p*f. Any killed
child lies in (-1,0] or [M,M+1), where f<=0, since M>=2. Killing
therefore INCREASES this value: B_j f>=mu_p*f on the surviving states.
Iteration and the return of the phase give P f>=mu_p^L*f. Positivity,
or the elementary Perron lower bound, yields r(P)>=mu_p^L. By (5),

    trace(P)>=mu_p^L.                                      (6)

Each trace path has exactly e up-steps and common probability
w=p^e*b^(L-e). A binary word is counted for at most M-1 starting
states. Rotate each such word at a minimum of its p-centered partial
sums. Its new p-centered sums lie in [0,M); its c-centered sums are
strictly positive at every positive time and less than M+1, because
0<(p-c)L<1. Each output word has at most L possible preimages. Thus

    A_(M+1)(L)/2^L >= mu_p^L / [(M-1)L 2^L w],             (7)
    rho_peak(L) >= 3^(-M-1) mu_p^L / [(M-1)L 2^L w].      (8)

All strict inequalities needed for Bad_L are preserved; equality in the
rational bridge is permitted only before adding the positive drift p-c.

## 4. Rounding entropy costs only a constant

Let h(p)=-p log p-(1-p)log(1-p). On [c,4/5], h'(p)>=-log4.
Because L(p-c)<1,

    2^(-L)/w=exp[-L(log2-h(p))] >= D^L/4.                 (9)

Euler's sine product gives the absolutely convergent identity

    log mu_p = -sum_(j>=1) a_j [1-p^(2j+1)-(1-p)^(2j+1)] t^(2j),
    a_j=zeta(2j)/(j*pi^(2j))>0.                          (10)

For p>=c>1/2 each bracket decreases with p. Thus mu_p>=mu_c
at the SAME t. Use M=m+4 in (7). The Robin supersolution proves
N_m(L)/2^L<=exp(1)D^L*mu_c^L at t=pi/(m+4). Equations (7),(9)
then prove (1).

The pairpeak note's Theorem C allows its K_0 to depend on L. Taking
K=5 and K_0=4 exp(1)L(L+3) proves (2). This is a shifted polynomial
comparison, not HYP-9142's uniform one-unit comparison.

## 5. The sharp constant with a logarithmic error

For t<=pi/2, the terms j>=2 in (10), after discarding the brackets,
sum to at most

    t^4/[180(1-t^2/pi^2)] <= t^4/135.

Here zeta(2j)<=zeta(4)=pi^4/90 and 1/j<=1/2. Also p(1-p)<=cd,
so log mu_p>=-cd*t^2/2-t^4/135.

Choose theta=(pi^2*cd/log3)^(1/3), M=ceil(theta L^(1/3)). The elementary
bounds 1<theta<2 ensure M>=3 for L>=8. Then

    (M+1)log3 + AL/M^2 <= kappa L^(1/3)+2log3,
    pi^4 L/(135M^4)<1,
    (M-1)L <= 3L^(4/3).

Substitute these and (9) into (8) to obtain the lower bound in (3).
The upper bound follows from the Robin note's finite private-price
upper bound and pi_L>=2rho_peak(L). This proves (3). Equation (2)
gives (4) for pi_L. THM-4480's proved sandwich
rho_peak(L)/M_L^*<=epsilon_L<=rho_peak(L), M_L^*=O(L^3), gives (4)
for epsilon_L as well. These are positive finite-horizon densities;
they do not establish descent for every positive integer.

## 6. Information preserved and lost

Source: a killed Bernoulli bridge with a rational slope. Target: exact
Collatz bad words with a hard upper wall. Map: round the endpoint slope,
take a closed transfer path, forget its initial height, and rotate at a
minimum. Preserved predicates: strict Collatz prefix growth and a
bounded peak. Lost data: the original source residue, rotation, initial
height, and the ownership or consistency of any pairing edits. Sidecars:
the multiplicities (M-1)L, the constant entropy loss, and the finite
phase-return condition pL in Z.

For q>=5, c=log_q2<1/2 and the tilted endpoint weight exceeds one.
The trace bridge ends near the bottom and does not supply the top-end
reward needed for the inherited sharp general-q constant. A bottom-to-top
killed-kernel lower bound is an additional requirement. The present
closed-bridge proof is restricted to q=3. The subsequent
[polynomial-cost connector construction](crossroads223_20260926_general_peak.md)
supplies that requirement and proves the general-q peak asymptotic;
the private pairing comparison remains restricted to q=3.

Classical identity: [DLMF 4.22.1](https://dlmf.nist.gov/4.22.E1).
The remaining matrix and path arguments are proved above. No literature
priority claim is made. The companion script records exact finite controls
separately from numerical probes of transcendental inequalities.

Reproduction: `python 04-computation/experiments/crossroads223_20260926_bridge.py`.
Normal and `-O` runs agree. Its exact controls include 304 transfer matrices,
60,648 two-by-two minors, 912 trace-power inequalities, and independent
word enumeration through L=14. Numerical spectral controls continue through
L=1024. A separate reviewer checked all 147 binary 3-by-3 TP2 matrices,
1,176 trace-power gates, and 108 independently built bridge universes.
