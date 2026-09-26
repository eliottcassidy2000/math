---
id: THM-4478
title: "Critical growth bands and affine integer capacity prove the sharp provability-price exponent"
status: >
  PROVED + INDEPENDENTLY AUDITED + FINITE-EXACT controls. Every deterministic
  modification of half-step Collatz which makes every positive n>=2 descend
  below itself within L steps must edit a set of lower natural density at
  least 2^(-(1-H_2(log_3(2)))L-O(sqrt(L log L))). The arbitrary-edit optimum
  has upper bound rho_L; the pairing-family optimum of THM-4475 has upper
  bound 2 rho_L. Consequently both have exponent 1-H_2(log_3(2)) and
  HYP-9137 is proved. This is a finite-horizon modification theorem, not
  convergence of the unmodified Collatz map.
source: collatz-crossroads-20260926 session
depends_on:
  - THM-4475-price-of-provable-descent-tends-to-zero
related:
  - HYP-9137-sharp-provability-price-exponent
  - THM-4477-cauchy-schwarz-price-bound-and-am-fair-criticality
  - THM-2163-radix-relation-carry-descent
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
scripts:
  - 04-computation/experiments/crossroads_20260926_ballot.py
  - 04-computation/experiments/crossroads_20260926_flow.py
  - 04-computation/experiments/crossroads_20260926_geometry.py
audit: >
  Independent geometry and automata lanes rederived the affine offset,
  integer interval count, allowed odd-count interval, first-hit injection,
  pair endpoint normalization, lower-density limit, cyclic-minimum count
  including repeated words, and entropy estimate. The direct-count proof
  avoids Haar averages, unproved independence, and weighted preimage
  multiplicity. Independent exact tests cover cyclic blocks through b=16,
  actual sources through 4*2^12 with L=6,8,10,12, and K=3,9,27; the geometry
  path independently checks sources<=5000, L=4,8,12, K=2,4,16.
---

# THM-4478 -- the sharp cost of forcing finite-time descent

**PROVED + INDEPENDENTLY AUDITED.** The reservation in `b54ea1ffe` was an
unproved namespace stub; this promotion supplies the proof and audited scope.

## 1. Setting and statement

Let T:N_(>0)->N_(>0) be the half-step Collatz map, T(n)=n/2 for even n and
T(n)=(3n+1)/2 for odd n. Put

    theta=log_3(2), h=H_2(theta), eta=1-h=0.050044... .

An L-step descent modification is any function G:N_(>0)->N_(>0) such that
for every n>=2 there is 1<=j<=L with G^j(n)<n. Let

    E(G)={v:G(v)!=T(v)}.

No computability, periodicity or density-existence assumption is imposed on
G. For a set E, write lower_d(E)=liminf_(X->infinity)|E intersect[1,X]|/X.

**Theorem A (general finite cut).** For every K>=1 let B_L(K) be the binary
words of length L whose prefix slopes w_j=3^(e_j)/2^j obey 1<=w_j<=K for
every 0<=j<=L. Put

    rho_L(K)=|B_L(K)|/2^L,
    R_L=sum_(k=0)^(L-1)(floor(k/3)+1),
    M_L(K)=(floor(log_3 K)+1) R_L.

Every L-step descent modification satisfies

    lower_d(E(G)) >= rho_L(K)/(K M_L(K)).              (1)

The same lower bound holds for the lower natural density of flipped *pair
indices* in the pairing family of THM-4475 (whose descent condition is only
n>=3). Finite exceptions do not affect either assertion.

**Theorem B (sharp exponent).** Let epsilon_L be the infimum of upper
natural edit densities over arbitrary L-step descent modifications. Let
delta_L be the analogous pairing-family infimum of THM-4475. For L>=8,

    2^(-eta L-O(sqrt(L log L))) <= epsilon_L <= rho_L,
    2^(-eta L-O(sqrt(L log L))) <= delta_L <= 2 rho_L,
    rho_L=2^(-eta L+O(log L)).                         (2)

Thus log_2(epsilon_L)/L and log_2(delta_L)/L both tend to -eta. The pairing
conclusion is exactly the exponent assertion of HYP-9137; it does not
assert a bounded or polynomial factor between delta_L and rho_L.

## 2. The affine interval retains the missing integer coordinate

Fix a parity word with all prefix slopes >=1. Composing branches gives

    T^k(n)=w_k(n+h_k),
    h_k=sum_(0<=i<k, bit_i=1) 1/(3w_i).

An even step changes neither h nor the normalized constant. An odd step
adds 1/(3w_i). Therefore

    0<=h_k<=k/3.                                      (3)

For a fixed endpoint v, depth k and number e of odd letters, all actual
integer ancestors following selected words satisfy

    v/(3^e/2^k)-k/3 <= n <= v/(3^e/2^k).              (4)

There are at most floor(k/3)+1 such integers. A starting integer has a unique
trajectory, so distinct word labels cannot create additional actual sources.
For 1<=3^e/2^k<=K, the allowable e lie in an interval of length log_3 K,
hence there are at most floor(log_3 K)+1 possibilities. Summing over k<L
shows that one endpoint v serves at most M_L(K) distinct selected sources.

This is an actual integer-spacing estimate. Counting symbolic inverse
branches or replacing incidence by a hub-weight distribution would lose it.

## 3. The first-hit cut and density normalization

The parity-word map is a bijection on residues modulo 2^L. For completeness,
each first branch reduces modulus by a factor 2; multiplication by 3 is a
unit on the odd branch. Induction proves bijectivity. Thus the selected
ordinary sources form a union S of residue classes of density rho_L(K).

For n in S every prefix has T^j(n)=w_j(n+h_j)>=n. If G descends by time L,
the two trajectories must first disagree at some state T^k(n) with k<L.
Assign n to that first edited endpoint. For n<=X, (3) gives

    T^k(n)<=K(X+L/3).

By the endpoint capacity (4),

    |S intersect[1,X]|-O(1)
      <= M_L(K) |E(G) intersect[1,K(X+L/3)]|.          (5)

Replacing X by Y/K-L/3 and taking the lower limit gives (1), even if E has
no natural density.

For pairing modifications each flipped pair {2i-1,2i} serves at most
2M_L(K) sources. Every hit pair has

    i <= (K(X+L/3)+1)/2.

The factor 2 in pair capacity cancels the factor 1/2 in this cutoff when
taking density of pair indices. Hence the same bound (1) holds.

## 4. Exponentially many sources inside a subexponential growth band

Choose a block length b, 1<=b<=L, and k_b=ceil(theta b). Every binary word
with k_b ones has total logarithmic slope in (0,log3). Rotating after a
minimum cumulative logarithmic slope makes every prefix nonnegative. A
rotated output has at most b input rotations, including when periodic.
There are therefore at least binom(b,k_b)/b admissible blocks.

Write L=tb+r, 0<=r<b. Concatenate t blocks and append r ones. Each block
ends with slope <3 and has internal slopes <=3^b. All the concatenated
prefix slopes lie in [1,K] for K=3^(t+b+r). Consequently

    rho_L(K) >= 2^-L (binom(b,k_b)/b)^t.               (6)

This gives an explicit finite bound on the right side of (1):

    2^-L (binom(b,k_b)/b)^t
       / [3^(t+b+r)(t+b+r+1) R_L].                   (7)

Stirling's elementary factorial bounds give

    log_2 binom(b,k_b)=b h-O(log b).

Thus the negative logarithm of (7) is at most

    eta L + O((L/b) log b + b + log L).

Taking b of order sqrt(L log L) proves the lower bounds in (2). Taking
b=floor(sqrt L) already proves the sharp exponent, with the weaker
O(sqrt L log L) error. No probabilistic meander approximation is needed.

## 5. Upper bounds and exact consequence

Let Bad_L^actual={n>=2:T^j(n)>=n for every 1<=j<=L}. On each fixed parity
cylinder, either every slope is >1, or some slope is <1 (there is no equality
3^e=2^j for j>0). In the first case every positive lift is bad; in the second
only finitely many lifts can be bad, by the affine formula. Hence
Bad_L^actual differs by finitely many points from the usual bad-word
residue union, of density rho_L.

Define G(n)=1 on Bad_L^actual and G(n)=T(n) otherwise. An originally bad
source descends immediately. An originally good source either follows T to
its descent within L, or meets an edited state first and then jumps to 1
within that horizon. This proves epsilon_L<=rho_L.

For the pairing-family upper bound use precisely Theorem A of
[THM-4475](THM-4475-price-of-provable-descent-tends-to-zero.md), which builds
a member of P_L with flip density at most 2rho_L for L>=8. Its partner
certificates are not replaced by the arbitrary-edit construction above.

Finally, rho_L is bounded above by the binomial tail with e_L>=theta L;
Chernoff's elementary entropy bound gives rho_L<=2^-eta L. The cyclic
minimum argument at block length L gives the matching lower bound up to a
polynomial factor. This proves (2).

## 6. Failure boundary and relation to moments

- The theorem concerns a modification with a *uniform fixed* horizon L.
  Collatz asks about one unmodified map with an n-dependent, unbounded
  convergence time. Nothing permits interchanging these quantifiers.
- Small-density or density-zero edits can change convergence. The known
  planted-divergence example remains a hostile to passing to a limit of G_L.
- The untrimmed conditioned second moment does not give the sharp exponent:
  its diagonal tilts odd letters to frequency 3/4 and stays bounded below.
  This proof selects a growth band, then uses exact integer incidence; it
  does not contradict THM-4477's distribution-only barrier.
- Both signs have the same word entropy. The mechanism is not a criterion
  selecting the positive plus root cycle over the minus cycles.
- For q>=5 the critical endpoint population is not the full bad population;
  the displayed q=3 sharp-rate conclusion does not transfer automatically.
- There is no Lean formalization claim and no literature priority claim.

Full connection ledgers, failed hypotheses and controls:
[ballot construction](../../05-knowledge/results/crossroads_20260926_ballot.md),
[flow and capacity](../../05-knowledge/results/crossroads_20260926_flow.md),
[independent geometry audit](../../05-knowledge/results/crossroads_20260926_geometry.md).
