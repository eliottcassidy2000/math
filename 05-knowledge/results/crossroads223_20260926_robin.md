# Sine supersolutions for the reflected Collatz barrier

**Status: PROVED + INDEPENDENTLY AUDITED by the flow and automata lanes.** No canon ID reserved.
The proposed finite bound is elementary apart from the classical sine
product. The matching private-price constant uses the inherited CITED
Mogul'skii lower bound, through THM-4480. Collatz and the globally
consistent pairing-price question remain OPEN.

## 1. Inheritance and the actual bridge

The source of inspiration is the old LRC (2,2,3) lift work: a boundary
operator must retain the guarded return state. The mathematical input is
the reflected barrier and multi-scale private-price inequality in
[the pairpeak note](procgen_pairpeak_20260926_pairing_peak_price.md),
sections 2--5. Its HYP-9142 asks for a uniform comparison with a hard wall.
We attempt a weaker direct spectral estimate sufficient for the sharp
private-price constant. No LRC theorem is a proof dependency.

Put

    c=log_3 2, d=1-c, lambda=d/c,
    D=2^(-(1-H_2(c))), A=pi^2 c d/2,
    kappa=min_(x>0)(x log3+A/x^2)
         =(3/2)(pi^2 c d)^(1/3)(log3)^(2/3).

The reflected process starts at x=0, dies at x<=0 after positive time,
and has two original binary symbols. Below x=m-1 they give increments
d and -c. At x>=m-1 both give -c; one is a flip. It stays below m-c.
Let N_m(L) be the number of surviving binary words of length L.

After the exact Bernoulli tilt from the pairpeak note, the sub-Markov
kernel B on positive states has weights c upward and d downward below
the zone, and weight 2d downward in the zone. With g(x)=lambda^x,

    N_m(L)/2^L = D^L (B^L g)(0).                         (1)

The initial state 0 is allowed only at time zero; transitions arriving
at nonpositive states are killed. Since c>1/2, 2d<1: the top is partially
absorbing. This factor is essential.

## 2. Finite spectral theorem

For all integers m>=2 and L>=0,

    N_m(L)/2^L <= exp(1) D^L exp(-A L/(m+4)^2).         (2)

Define H=m+4, t=pi/H and

    beta=log[d sin(ct)/(c sin(dt))],
    f(x)=exp(beta x) sin(t(x+1)),
    mu=d exp(-c beta) sin(t)/sin(dt).                   (3)

The exact interior identity is

    c f(x+d)+d f(x-c)=mu f(x).                          (4)

Indeed c exp(beta d) sin(dt)=d exp(-beta c) sin(ct),
so the cosine-phase coefficient vanishes. The remaining coefficient is
the displayed mu by the angle-addition identity. Below zero, killed
values of f would be positive, since x-c+1>=d>0. Thus killing changes
the equality to an upper bound for Bf.

In the zone x in [m-1,m-c), put r=H-(x+1). Then r>=3+c. The inequality
`2d f(x-c)<=mu f(x)` reduces to

    2 sin(t(r+c))/sin(tr) <= sin(t)/sin(dt).            (5)

The left ratio decreases with r where its numerator and denominator are
positive, so it suffices to use r=3+c. Since sin(z)/z decreases on (0,pi),
the left side is at most `2(1+c/(3+c))`. On the right,
`sin(t)/sin(dt)>=(1-t^2/6)/d`. The exact coarse bounds
`5/8<c<2/3`, `t<=pi/6`, and pi^2<10 imply

    2d(1+c/(3+c)) <= 39/44 < 1-pi^2/216 <= 1-t^2/6.

Consequently `Bf<=mu f` everywhere, including the initial state.

The same sinc monotonicity gives beta<0. From
`sin(ct)/(ct)>=1-c^2 t^2/6` and `sin(dt)/(dt)<=1`,

    0<=-beta<=c^2 t^2/3,
    m|beta| <= (4/27) pi^2 m/(m+4)^2 < 1.

On 0<=x<m-c, both distances of x+1 to the endpoints 0,H are at least one,
so `f(x)>=exp(-1) sin(t)=exp(-1)f(0)`. As 0<g(x)<=1,
positivity and iteration give

    (B^L g)(0)<=exp(1) mu^L.                           (6)

## 3. The sine product preserves the sharp coefficient

Using Euler's sine product, or its convergent logarithmic power series,

    log mu = log sinc(t)-c log sinc(ct)-d log sinc(dt)
           = -sum_(j>=1) a_j [1-c^(2j+1)-d^(2j+1)] t^(2j),
    a_j=zeta(2j)/(j pi^(2j))>0.

Every bracket is positive; the j=1 term is `-cd t^2/2`. Therefore

    log mu <= -cd t^2/2 = -A/(m+4)^2.                 (7)

Equations (1), (6), (7) prove (2). The identity is classical:
[DLMF 4.22.1](https://dlmf.nist.gov/4.22.E1). Its positive coefficients
can equivalently be written in terms of the even Bernoulli numbers.
This is an exact coefficient-sign argument, not an appeal to an analogy
between the small and large indices of the Bernoulli sequence.

## 4. Sharp private-price consequence

Let pi_L denote the private pairing price, as defined in the pairpeak
note; each source is allowed its own pairing modification. Its proved
multi-scale inequality is

    pi_L <= 2L [3^(1-M)rho_L
         +sum_(m=2)^(M-1)3^(1-m) N_(m+1)(L)/2^L]
         +2 N_2(L)/2^L.                               (8)

Take M=L (L>=8). For x=m+5, (2) and elementary minimization give

    3^(1-m) N_(m+1)(L)/2^L
      <= exp(1) 3^6 D^L exp[-x log3-A L/x^2]
      <= exp(1) 3^6 D^L exp[-kappa L^(1/3)].           (9)

The last term of (8) obeys the same bound with x=6. Also rho_L<=D^L and
`(L-1)log3>=kappa L^(1/3)` for L>=8. Thus a simple uniform version is

    pi_L <= C (L+1)^2 D^L exp[-kappa L^(1/3)]  (L>=8),
    C=2(1+exp(1)3^6).                                 (10)

The known lower bound pi_L>=2rho_L^peak and THM-4480's sharp asymptotic
then give

    log pi_L = L log D-kappa L^(1/3)+o(L^(1/3)).        (11)

Only this matching lower bound invokes Mogul'skii, in the form of
[Gantert--Hu--Shi 2011, Lemma 2.1](https://www.numdam.org/item/10.1214/10-AIHP362.pdf),
already used and audited in THM-4480. The present finite upper bound
does not require a uniform Robin/Dirichlet ratio.

## 5. Boundaries and decisive tests

- (10) alone does not prove pi_L<=poly(L)rho_L^peak: the inherited
  asymptotic lower bound has an o(L^(1/3)) remainder, not O(log L).
- No inference about one globally consistent pairing follows from private
  certificates. The companion flow lane finds explicit incompatibilities.
- This does not prove HYP-9142's uniform K=1 Robin comparison.
- The absorption 2d<1 is load-bearing; at c=1/2 the boundary becomes
  reflecting and the Dirichlet leading constant need not survive.
- The sine comparisons are proved for all real states, so the irrational
  changing fractional part of ct does not require a finite periodic matrix.

The companion independent script compares exact reflected counts to
57,330 direct binary word cases, replays 72 positive integer realizations,
and numerically probes 158,772 real-state inequalities and 172,032 exact
count/transcendental-bound comparisons through L=4096. The latter are
VERIFIED numerical probes, not FINITE-EXACT proofs of transcendental
inequalities. Both independent reviewers rederived the kernel, boundary
comparisons, sine coefficient, and multi-scale consequence. The flow
reviewer additionally tested 398,199 real-state samples and 2,280 exact
DP count comparisons without importing the root script.
