---
id: THM-4492
title: "Finite cutoff kernels separate lower and ordinary two-step pairing density, force oscillation, and give upper density at least 0.3067797897"
status: "PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT rational certificates"
source: "crossroads233-20260926 flow lane; root proof audit and independent full-cost computation"
depends_on:
  - 01-canon/theorems/THM-4491-pairing-two-step-tree-extra-density.md
proofs:
  - 05-knowledge/results/crossroads233_20260926_flow.md
scripts:
  - 04-computation/experiments/crossroads233_20260926_flow.py
  - 04-computation/experiments/crossroads233_20260926_scale_audit.py
---

# THM-4492: shared pairings at several cutoffs

**PROVED + independently audited.** Reserved in checkpoint bb630ddb8,
then promoted after the proof and exact controls were audited. This is a
theorem about the globally consistent TWO-STEP pairing family. It neither
proves Collatz nor settles longer-horizon global peak-price gluing.

## Object and main consequences

Pair i={2i-1,2i} has a bit epsilon_i. Bit zero sends the odd/even members
to 3i-1 and i; bit one sends them to i-1 and 3i. Let F={i:epsilon_i=1}.
Global P_2 means every n>=3 descends below itself within two steps.
By [THM-4491](THM-4491-pairing-two-step-tree-extra-density.md), its exact
clauses on pair indices at least two are

    epsilon_(2k)+epsilon_(3k)=1,
    epsilon_(3k+1)<=epsilon_(2k+1)<=epsilon_(3k+2).

They form a rooted tree with every finite feasible prefix extendible.
Write A_F(X)=|F intersect[1,X]| and let ell,u be its lower/upper natural
densities. Let alpha be the attained minimum LOWER density from THM-4491.

Every global P_2 member satisfies

    u >= 821510388809/2677850419968
       = 0.30677978974599224... .                         (1)

Any existing natural density satisfies the same bound. Moreover, with

    b=4538724347/15109399071=0.30039079156439336...,

one has

    4u+9ell >=13b,       9u+4ell >=13b.                  (2)

Since THM-4491 proves alpha<=U=3089623223/10460353203<0.295366,
no P_2 member has natural density alpha. Every member attaining lower
density alpha necessarily has

    u >=13417603/43046721=0.3116986076593383...,
    u-ell >=170854306/10460353203=0.016333512137142698....

The previous lower-density attainment theorem remains valid. Its
construction must oscillate. Minimum ordinary and upper densities,
and optimality of the bounds here, remain OPEN.

## General finite cutoff functional

Put r=2/3. For nonnegative rational a_0,...,a_K, not all zero, define

    J_a(X)=min_(global P_2 F) sum_k a_k A_F(floor(r^k X)),
    W=sum_k a_k,       Lambda=sum_k a_k r^k.

The minimum retains ONE common assignment across all cutoffs. Its tree
weights are sum_k a_k 1{i<=floor(r^k X)}.

For a full depth-h subtree put w_h=sum_(k<=min(h,K))a_k, and give each
child the corresponding depth-(h-1) pattern. Let Delta_h=F_h^1-F_h^0
and tau_h=min_b F_h^b minus the sum of its child minima. Initially
Delta_0=a_0,tau_0=0. For child differences x,y,

    even i: Delta_h=w_h-x,
            tau_h=min(max(x,0),w_h);
    odd i:  Delta_h=w_h+min(x,0)+max(y,0),
            tau_h=min(max(-x,0),w_h+max(y,0)).

These values depend on i modulo2^h and satisfy
|Delta_h|<=W(h+1), 0<=tau_h<=Wh. Set

    B_h(a)=sum_(i mod2^h)tau_h(i),
    gamma_a=sum_(h>=1)B_h(a)/3^(h+1).

Then

    J_a(X)/X -> gamma_a,       u >= gamma_a/Lambda.      (3)

Every finite truncation of the nonnegative series is a certified lower
bound in (3). Its tail is at most W(H+3)(2/3)^(H+1).

**Proof of the scale limit.** A depth-j descendant of i lies between
(3/2)^j(i-1)+1 and (3/2)^j(i+1)-1. Thus roots in
(r^(h+1)X,r^h X] have full depth h apart from bounded endpoint errors,
and a depth-j descendant is below cutoff floor(r^k X) exactly when
j+k<=h, again away from those errors. The weight pattern is w_(h-j).
The band's length times its residue average gives
(1/3)r^h X * B_h/2^h = X B_h/3^(h+1).
Local tolls telescope J_a(X). For a general truncated subtree the same
depth bound holds with weights in[0,W], so indices i<=epsilon X contribute
O(W epsilon X(1+log(1/epsilon)))+O(W log X). This proves the series by
first summing finitely many bands, then sending epsilon to zero.
Finally every common assignment has weighted cost at least J_a(X),
while its limsup divided by X is at most Lambda u. This proves (3).

## Two scales and forced oscillation

Take a=(1,0,1), so Lambda=13/9. Then J_a(X) is the optimum of
A_F(X)+A_F(floor(4X/9)) and let beta=gamma_a/Lambda.
Taking the larger cutoff along a liminf sequence gives
ell+(4/9)u>=gamma_a. Taking the smaller cutoff along a liminf sequence,
with the larger cutoff rounded by9/4, gives u+(4/9)ell>=gamma_a.
Rounding changes only O(1) entries. These are (2), since beta>=b.

The twenty exact toll sums for this kernel are

    1,2,6,12,30,60,124,250,498,982,1996,3958,
    8008,16062,32318,64712,129682,259672,519980,1041548.

They give beta>=b as displayed. Already twelve terms give
beta>=684037/2302911>U, proving the strict separation from alpha.
Substitution ell=alpha<=U in the first inequality gives the stated
oscillation fractions.

The smallest finite witness in this cutoff family is X=5,N=2:
C(2)=0,C(5)=1 but min_F[A_F(2)+A_F(5)]=2. The clauses
epsilon_2+epsilon_3=1 and epsilon_3<=epsilon_5 expose the incompatibility.
This is a conflict between shared assignments, not a directed comparison
that would naturally define a tournament.

## Eight scales and exact validation

For (1), use

    a=[128,192,288,432,648,972,1458,2187],
    W=6305, Lambda=1024.

Each a_k r^k equals128. Exact recurrence gives

    sum_(h=1)^20 B_h(a)3^(20-h)=3286041555236.

Divide by1024*3^21 to obtain (1). Independently, direct arrays for the
two fixed-root costs give sum_(i mod2^20)min_b F_20^b(i)=3286041555236.
The equivalence follows from S_h=B_h+3S_(h-1), because each child
residue occurs three times across all parent residues. This independent
computation avoids the difference/toll recurrence.

The [full proof](../../05-knowledge/results/crossroads233_20260926_flow.md)
records the bounded positive-kernel comparison, all coefficients, and the
normalization and pruning audits. Controls include exhaustive binary
assignments, independent fixed-root tree costs, and direct prefix optima.
Both scripts pass normally and with Python assertions disabled.

The inherited P_4 example violates a P_2 clause. Therefore the exact
tree and these density bounds cannot be transferred to longer horizons
without rebuilding their legal path constraints. HYP-9140 remains OPEN.
