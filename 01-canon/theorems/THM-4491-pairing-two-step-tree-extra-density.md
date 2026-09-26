---
id: THM-4491
title: "The exact minimum lower density of two-step pairings is a convergent tree-cost series; disjoint clauses give 29/108 and twenty coefficients give 0.2907539"
status: "PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT coefficients and hostile controls"
source: "crossroads223-20260926 flow lane; independent root and geometry proof audits"
depends_on: []
proofs:
  - 05-knowledge/results/crossroads223_20260926_flow_global.md
  - 05-knowledge/results/crossroads223_20260926_flow_series.md
scripts:
  - 04-computation/experiments/crossroads223_20260926_flow_global.py
  - 04-computation/experiments/crossroads223_20260926_flow_series.py
---

# THM-4491: actual pair ownership in two-step descent

**PROVED + independently audited.** Reserved in checkpoint `455ba57a2`,
promoted after the root and geometry audits. This is a theorem about
modified two-step pairings. It does not prove Collatz or longer-horizon
price asymptotics. Minimum LOWER density is distinguished from natural
or upper density throughout.

## Map and exact tree

Pair i is {2i-1,2i}. For a bit epsilon_i,

    0: 2i-1 -> 3i-1, 2i -> i;
    1: 2i-1 -> i-1,  2i -> 3i.

Write F={i:epsilon_i=1}, and P_L for the property that every n>=3 has
an iterate strictly below n within L steps. Then P_2=P_3, and P_2 is
equivalent to the clauses, for source-pair indices i>=2,

    i=2k:   epsilon_(2k)+epsilon_(3k)=1;
    i=2k+1: epsilon_(3k+1)<=epsilon_(2k+1)<=epsilon_(3k+2).

Indeed first descent within two steps means D or UD; UUD never descends,
where U(n)=ceil(3n/2), D(n)=floor(n/2). The clauses form a tree on i>=2:
every j>=3 has unique smaller parent 2j/3, (2j+1)/3, or (2j-1)/3
according as j is 0,1,2 modulo3. Every feasible prefix extends globally.
Thus C(X)=min_(P_2)|F intersect[1,X]| is an exact finite tree optimum.

## Exact limiting series and attained lower density

Let F_h^b(i) be the full depth-h subtree optimum with root bit b. Set
Delta_h=F_h^1-F_h^0 and tau_h=min_b F_h^b minus the sum of its child
subtree optima. Initially Delta_0=1, tau_0=0. For h>=1, with x,y the
depth-(h-1) child differences,

    i even: Delta_h=1-x, tau_h=1{x>=1};
    i odd:  Delta_h=1+min(x,0)+max(y,0),
            tau_h=min(max(-x,0),1+max(y,0)).

One has |Delta_h|<=h+1 and 0<=tau_h<=h. The residue i mod2^h determines
tau_h. Put A_h=sum_(i mod2^h)tau_h(i). Then

    lim_(X->infinity) C(X)/X = alpha=sum_(h>=1)A_h/3^(h+1),
    0<=alpha-sum_(h=1)^H A_h/3^(h+1)<=(H+3)(2/3)^(H+1).

Every global P_2 member has lower density at least alpha, and one member
attains lower density alpha. Subsequent
[THM-4492](THM-4492-pairing-two-cutoff-density-separation.md) proves that
natural-density attainment of alpha is impossible: every member attaining
this lower density has upper density at least0.3116986076.
The first twenty exact coefficients give

    3041388727/10460353203 <= alpha <= 3089623223/10460353203,
    0.2907539227382626... <= alpha <= 0.2953650955221956... .

**Proof mechanism.** The local increments telescope C(X). All descendants
at depth h lie between (3/2)^h(i-1)+1 and (3/2)^h(i+1)-1. Thus, away
from bounded endpoint errors, roots in the band
((2/3)^(h+1)X,(2/3)^h X] have full subtree depth h. Their residue average
contributes A_h/3^(h+1). Roots i<=epsilon X contribute at most
O(epsilon X(1+log(1/epsilon)))+o(X), by the depth bound. This proves
the series and its geometric tail. Forcing a previously chosen prefix
through Y costs at most O(Y(1+log X)) at a later cutoff X: at most 2Y
disjoint child subtrees change their root bit. Conditional optima at
successive-square cutoffs therefore yield a consistent infinite member
whose density along those cutoffs tends to alpha. The universal lower
bound proves that its liminf is exactly alpha.

The [complete series proof](../../05-knowledge/results/crossroads223_20260926_flow_series.md)
retains the truncated-child cases, periodicity, tail and attainment
quantifiers. Different prefix optima alone would not supply attainment.

## A transparent finite cut and its horizon boundary

The canonical disjoint matching {2r,3r}, v_3(r) even, costs X/4+O(log X)
flips through X. The clauses additionally imply, for every s>=0,

    epsilon_(36s+17)+epsilon_(54s+23)>=1.

Use the chain at indices i=16s+7, j=24s+10, h=24s+11, l=36s+15,
a=36s+17, b=54s+23:
epsilon_j+epsilon_l=1, epsilon_j<=epsilon_i<=epsilon_h<=epsilon_a,
epsilon_l<=epsilon_b. The two endpoint progressions are disjoint and
lie outside the matching. Hence every P_2 member has lower density
at least 1/4+1/54=29/108. This direct proof is weaker numerically than
the twenty-term series but exposes actual once-per-pair ownership.

The cut fails for P_4. Start with the globally valid recursion
b1=b2=0, b_(3k)=1-b_(2k), b_(3k+1)=0, b_(3k+2)=b_(2k+1), and change
only bit23 from one to zero. Bits17 and23 are now both zero. The only
source needing more than two steps is

    30 -> 45 -> 68 -> 34 -> 17.

Only map values45,46 changed; sources above46 encounter them only after
descending. Exact verification of3..46 completes this global P_4 proof.
It refutes the longer-horizon CLAUSE, not the numerical density bound.

The [global proof and controls](../../05-knowledge/results/crossroads223_20260926_flow_global.md)
include 65,534 brute Boolean assignments and all64 local six-bit cases.
The series script independently checks all254 residue types through
depth7 and computes twenty exact coefficients. Normal and `-O` runs agree.
