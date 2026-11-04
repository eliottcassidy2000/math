# Session 178: a finite valuation bridge from LCM closure to cut obstructions

## Status and contribution

This is a proved refinement of Session 177's Euler-aware maximal sigma graph, not a proof of nonexistence of odd perfect numbers or an improvement on known lower bounds. It gives an explicit finite exponent box from a fixed prime support, sound capacity propagation, LCM conjunction pruning, and exact cut constraints. LTE, factor-chain methods, and fixed-support finiteness are standard ingredients; historical novelty is **not established**. The contribution here is the self-contained synthesis and executable certificates. No outside literature claim is needed for the proof.

Let S contain s>=2 odd primes and let q in S satisfy q=1 mod 4. Suppose N is perfect with precisely this support and Euler prime q. Write N=product p^{a_p}, n_p=a_p+1, with a_q=1 mod 4 and all other a_p positive even. Put b_q=2 and b_p=1 otherwise. Let E be Session 177's maximal allowed edge set. For p!=r let d_pr=ord_r(p); an edge is allowed if d_pr is odd for p!=q, or v_2(d_pr)<=1 for p=q. Define

- k_pr=r if d_pr=1, and k_pr=d_pr otherwise;
- ell_pr=lcm(b_p,k_pr);
- c_pr=0 if d_pr=1, and c_pr=v_r(p^{d_pr}-1) otherwise.

No self-edge exists since sigma(p^a)=1 mod p.

## 1. Exact valuation law

For an allowed edge,

v_r(sigma(p^{a_p})) = 0 if k_pr does not divide n_p;

otherwise it equals v_r(n_p) if d_pr=1, and c_pr+v_r(n_p) if d_pr>1.

Indeed, for order greater than one the denominator p-1 is prime to r; odd-prime LTE applied to (p^{d_pr})^{n_p/d_pr}-1 gives the formula, since r does not divide d_pr. For order one, LTE subtracts v_r(p-1) and leaves v_r(n_p). Disallowed edges have valuation zero for all exponents of the designated role. Thus perfectness forces the exact column balances

(**V**) a_r = sum_{p!=r} v_r(sigma(p^{a_p})).

These are balances at each odd prime, not merely graph incidence constraints.

## 2. Effective exponent-box theorem

Put C_r=sum_{p:(p,r) in E} c_pr and C=max_r C_r. Let B be any positive integer such that

(**B**) C+(s-1)(floor(log_3(B+1))+1) < B.

Such B can be found by successively doubling from 2. Then **every a_p < B**.

Proof. Let A=max a_p and choose r attaining the maximum. The exact law gives

A <= C_r + (s-1) log_r(A+1) <= C+(s-1)log_3(A+1).

Set m=s-1. Condition (B) implies B>m. The function f(x)=x-C-m log_3(x+1) is strictly increasing for x>=B, since f'(x)=1-m/((x+1)ln 3)>0, and (B) implies f(B)>0. Therefore A>=B contradicts the displayed necessary inequality. QED.

All constants are computable from the support alone. The theorem ignores primes outside S when deriving the bound; this weakens the test but does not invalidate it. It does not claim a uniform exponent bound independent of S.

## 3. Sound capacity and LCM propagation

Start with U_p=B-1. Delete any edge whose ell_pr exceeds U_p+1. Given a surviving edge, define its upper valuation capacity

M_pr(U_p) = c_pr+floor(log_r(floor((U_p+1)/ell_pr))) if d_pr>1;

M_pr(U_p) = floor(log_r(floor((U_p+1)/b_p))) if d_pr=1.

For deleted edges use zero. These capacities are upper bounds, not necessarily attained by an exponent of the required parity: this is enough. For d>1, n_p is divisible by ell_pr and r is prime to ell_pr, hence its r-adic valuation is bounded as stated. For d=1, write n_p=b_p t and use r odd. The deletion criterion already handles absence of any divisible n_p.

Replace U_r by the smaller of U_r and the largest exponent of r's role not exceeding sum_p M_pr(U_p). If this drops below 1 for q or below 2 for another prime, reject. Recompute edges and repeat. Every step preserves every hypothetical solution. If the surviving graph loses strong connectivity, reject by Session 177's theorem. These operations terminate because integral nonnegative upper bounds only decrease (or a contradiction is reached).

For a branch requiring several outgoing edges p->r, let

L_p=lcm(b_p,{k_pr: required targets r}).

If L_p>U_p+1, the branch is impossible. Otherwise every edge with k_pt dividing L_p is forced, exactly as in Session 177's LCM closure theorem. Thus **LCM closure and valuation bounds feed each other**: new required edges increase L_p; forced edges add incoming valuation lower bounds; capacity bounds shrink U; shrinking U can forbid a conjunction even when every edge separately survives. For example, at source 3 with bound U_3=10, requirements 3->11 and 3->13 are separately possible (labels 5 and 3), but together require n_3>=15, so are impossible. If allowed, their conjunction also forces 3->4561, an external-prime obstruction for supports omitting 4561.

One may strengthen M by enumerating n in its finite role-compatible domain and taking the actual maximum. The closed formulas are cheap safe bounds.

## 4. Exact cut balance and its strict direction

For any nonempty proper T subset S put A_T=product_{p in T}p^{a_p}. With actual valuations w_pr define

O_T=product_{r outside T}r^{sum_{p in T}w_pr},

I_T=product_{r in T}r^{sum_{p outside T}w_pr}.

The column balances (V), absence of outside-support odd factors, and v_2(sigma(N))=1 imply

(**Cut**) sigma(A_T)/A_T = 2^{[q in T]} O_T/I_T.

The factor 2 comes only from the Euler source. In fact v_2(sigma(q^{a_q}))=1, by q=1 mod4 and a_q=1 mod4; all other sigma factors are odd.

Every proper nontrivial divisor A_T of a perfect number satisfies

1 < sigma(A_T)/A_T < 2,

since sigma(N)/N is the product of the abundancies of A_T and its coprime nontrivial complement. Consequently:

- if q is outside T, then I_T < O_T < 2 I_T;
- if q is inside T, then O_T < I_T < 2 O_T.

This determines the multiplicative direction across every cut, beyond the existence of crossing edges. In logarithms, the difference log O_T-log I_T lies strictly between 0 and log 2 for Euler-free cuts, and strictly between -log 2 and 0 for Euler-containing cuts.

**Recovery of the graph obstruction.** If T has no outgoing edges, O_T=1. For q outside T the cut interval demands I_T<1, impossible. For q inside T it demands 1<I_T<2, impossible because I_T is an integer. Every proper nonempty T must therefore have an outgoing edge, which is equivalent to strong connectivity of the actual finite graph. Thus these exact cut intervals strictly retain the earlier graph obstruction while recording the sizes of crossing factors.

**Conditional cut certificates.** Suppose a branch yields valid bounds O_low<=O_T<=O_high and I_low<=I_T<=I_high. If q is outside T, reject when O_high<=I_low or O_low>=2I_high. If q is inside T, reject when I_high<=O_low or I_low>=2O_high. These implications use non-strict bounds to contradict the strict cut interval, so endpoints matter. The bounds must describe the same branch; independently attainable extremes are valid bounds but need not be jointly attainable. Rational/integer comparisons avoid floating-point ambiguity.

## 5. Concrete support rejected beyond the maximal-graph sieve

Take S={3,5,7,11,13}. Its abundancy ceiling is 1001/384>2. The only Euler candidates are 5 and 13. For each, the Session 177 maximal allowed graph is strongly connected.

**q=5.** Here C=3 and B=16 meets (B): 3+4*3=15<16. Thus every exponent is <=15. The only permitted source into 13 is 3. Its order label is 3 and c_3,13=1. The capacity is 1+floor(log_13(floor(16/3)))=1. Thus a_13<=1, contradicting its positive even role. This is an exponent-independent rejection of the whole q=5 branch.

**q=13.** Here C=4 and B=32 meets (B): 4+4*4=20<32. One capacity pass, rounded down by role, yields

(a_3,a_5,a_7,a_11,a_13) <= (4,2,4,4,1).

The only original outgoing edge of 5 is 5->11, with label 5. Its new source bound requires n_5<=3, so this edge vanishes. A vertex of outdegree zero contradicts necessary strong connectivity. Thus this branch also fails.

These five-prime examples were already excluded by known general prime-count bounds. Their purpose is to prove that this implementable arithmetic refinement detects something the preceding maximal-graph plus abundancy test does not. They are not evidence for a new numerical OPN bound.

## 6. Complete finite decision problem, and the unresolved issue

After obtaining the exponent box, enumerate each source's role-compatible exponent a<=U_p. Discard it if sigma(p^a) has any odd prime factor outside S. Store the exponent and its complete valuation row. Choose one row per source, imposing (V) in every column. This is an **exact finite decision procedure** for an OPN with support S and Euler prime q: necessity is immediate; sufficiency follows because all odd valuations of sigma(N) then equal those of N, its 2-adic valuation is exactly one, and there are no remaining prime factors. LCM requirements and cut constraints are sound filters for this finite problem, not substitutes for its final exact balances.

The unresolved global problem is to control the support itself. No argument here bounds the largest support prime or the number of primes independently of N, so no global nonexistence conclusion follows. The concrete cut filters above are proved, but no large-support computational gain or new benchmark record is claimed.

## Reproduction

Run `python check_primeflow.py` in this directory. The standard-library checker independently compares the LTE laws and capacity formulas with directly constructed sigma values, checks all LCM patterns for a finite range, and produces the two support certificates. Numerical testing corroborates the implementation; the proofs establish the unbounded claims.
