---
id: THM-4488
title: "Private pairing price is polynomially comparable to peak-discounted bad paths; rational bridges give the sharp q=3 constant with O(log L) error"
status: "PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT controls and separately labelled VERIFIED numerical probes"
source: "crossroads223-20260926; root sine/bridge proofs, independent flow and automata audits"
depends_on:
  - 05-knowledge/results/procgen_pairpeak_20260926_pairing_peak_price.md
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md
proofs:
  - 05-knowledge/results/crossroads223_20260926_robin.md
  - 05-knowledge/results/crossroads223_20260926_bridge.md
audit: 05-knowledge/results/crossroads223_20260926_bridge_audit.md
scripts:
  - 04-computation/experiments/crossroads223_20260926_robin.py
  - 04-computation/experiments/crossroads223_20260926_bridge.py
  - 04-computation/experiments/crossroads223_20260926_bridge_audit.py
---

# THM-4488: private pairing price and rational bridges

**PROVED + INDEPENDENTLY AUDITED.** Reserved in checkpoint `5a0338e38`,
then promoted after two independent proof audits and separate exact controls.
This is a finite-horizon theorem for private pairing certificates and
arbitrary edits. Collatz and the globally consistent pairing comparison
[HYP-9140](../../05-knowledge/hypotheses/HYP-9140-pairing-price-is-peak-discounted.md)
remain OPEN.

## Definitions and statement

Put c=log_3 2, d=1-c, D=2^(-(1-H_2(c))), A=pi^2*c*d/2 and

    kappa=(3/2)(pi^2*c*d)^(1/3)(log3)^(2/3)=2.10758042868... .

For a binary length-L word, write S_j=e_j-cj and call it bad when
S_j>0 for all 1<=j<=L. Let A_M(L) count bad words with S_j<M for
0<=j<L, and put

    rho_peak(L)=2^(-L) sum_bad 3^(-max_(0<=j<L) S_j).

The reflected count N_m(L) uses the walk with increments d,-c below
m-1 and two down symbols in the zone x>=m-1, killed on arrival at x<=0;
the initial state zero is allowed. The private price pi_L is the limsup
mean, over sources n>=3, of the least certificate cost sum_(flipped i) n/i.
Each source may use a different pairing member. Epsilon_L is THM-4480's
least density of arbitrary edits forcing every n>=2 to descend within L.

For L>=3 and m>=2,

    N_m(L) <= 4 exp(1)(m+3)L A_(m+5)(L).                  (1)
    2 rho_peak(L) <= pi_L
      <= [4 exp(1)3^5 L(L+3)(27L+18)+4L2^(-L)]rho_peak(L). (2)

For L>=8, with C=2(1+exp(1)3^6),

    D^L exp(-kappa L^(1/3))/(108 exp(1)L^(4/3))
      <= rho_peak(L)
      <= (C/2)(L+1)^2 D^L exp(-kappa L^(1/3)).           (3)

Consequently, for f equal to rho_peak, pi, or epsilon,

    log f_L = L log D-kappa L^(1/3)+O(log L).            (4)

## Proof mechanism

The [Robin proof](../../05-knowledge/results/crossroads223_20260926_robin.md)
constructs f(x)=exp(beta*x)sin(pi*(x+1)/(m+4)) for the tilted kernel.
Its top weight is 2d<1. An exact trigonometric identity and boundary
inequalities give Bf<=mu*f and N_m/2^L<=exp(1)D^L*mu^L. The convergent
sine product gives log(mu)<=-A/(m+4)^2. Inserting this into the inherited
multi-scale private-price inequality gives the upper side of (3).

The [rational-bridge proof](../../05-knowledge/results/crossroads223_20260926_bridge.md)
sets p=ceil(cL)/L and kills a Bernoulli-p walk outside (0,M). Its L-step
transfer matrix P has nonnegative 2-by-2 minors by tail switching. A
finite matrix lemma proves trace(P)>=r(P). The unshifted sine function
is a subsolution, so trace(P)>=mu_p^L. Forgetting the initial height and
rotating at a minimum gives at least trace(P)/[(M-1)L] weighted words
with strict Collatz prefix growth and height less than M+1.

Rounding costs at most a factor four in entropy. Also mu_p>=mu_c for
p>=c>1/2, by the coefficient signs of the sine product. With M=m+4 this
proves (1); the pairpeak note's Theorem C gives (2). Optimizing
M=ceil((pi^2*cd/log3)^(1/3)L^(1/3)) and bounding the fourth-order
remainder proves the lower side of (3). The arbitrary-edit sandwich
rho_peak/M_L^*<=epsilon_L<=rho_peak with M_L^*=O(L^3) proves (4).
All constants, boundary checks, and multiplicities are explicit in the
linked complete proofs.

The only classical analytic identity used here is Euler's sine product,
[DLMF 4.22.1](https://dlmf.nist.gov/4.22.E1). Its logarithmic coefficient
is a_j=2^(2j-1)|B_(2j)|/[j(2j)!]>0. The even Bernoulli numbers therefore
provide an actual sign-controlled estimate. Mogul'skii is no longer
needed for this q=3 refinement of THM-4480.

## Scope and checks

- (1) has a five-unit shift and polynomial loss. The uniform one-unit
  Robin comparison [HYP-9142](../../05-knowledge/hypotheses/HYP-9142-robin-inequality-barrier-price.md)
  remains OPEN.
- A trace bridge controls paths, not the simultaneous ownership of flips.
  No globally consistent delta_L upper bound follows.
- For q>=5 the endpoint tilt has the opposite sign. The closed-bridge
  argument above alone does not give the sharp general-q constant.
- All exact scripts pass normally and with `-O`. Root checks include
  60,648 transfer minors and direct enumeration through L=14; an independent
  audit exhausts 147 binary TP2 matrices and 108 bridge universes.
  Transcendental numerical probes are not presented as exact proof.
