# Independent audit: Session30 arrival-frontier contraction

Verdict: PASS within the explicitly scoped claims. No uniform 153N theorem, universal entry theorem, or recovery of the user-reported composed Collatz bounds is established.

## Independence and reproducibility

The independent script `audit.py` imports no producer modules. It reconstructs the native P map, charges, source families, stopped routes, quadratic algebra, and exact residue counts. SymPy is used for exact polynomial/rational identities; literal checks use integers. All checks explicitly raise on failure, so optimization cannot remove them. Run `python3 audit.py` and `python3 -O audit.py` from this directory. Their receipts agree byte for byte. The producer's original SHA256SUMS also checks successfully. This audit reads the native frontier JSON and compares its size and maximum with the declared constants.

Final finite scope: 811,022 explicit checks; 32,563 clusters and 40,711 child sources, across every N from 1 through 512 and N=1440,4096,7958,7959,8192,63666,65536,100000. At N=63666, 6,631 independent child arrival computations agree with the parent-derived terminal, depth, and actual full clock within 500 native steps. Eleven exact symbolic identities independently verify the algebra and native inverse-chain formulas. Actual source IDs, clocks, and residues modulo both 15 and 35 contribute to the receipt digest. This is a finite control, not a universal arrival assertion.

## Algebra and exterior audit

Write Q(a,b)=a²−ab+b². For completed signed coefficients x,y, old charge is W=x(1,2)+y(1,−1) and incoming charge is H=x(−1,−1)−y(1,−1). Direct expansion gives

Q(W+H)−Q(W)=−2x²+3xy−3y²
=−2(x−3y/4)²−15y²/8.

The completed sum is (0,x), while W1+W2=3x. Thus S(W)=(0,(W1+W2)/3) and

4Q(W)−9Q(SW)=3(W1−W2)².

The 4/9 constant is sharp even on the bank lattice: W=(3,3), x=2,y=1. This sharp algebraic witness is not a claim that a native arrival configuration realizes those signed counts. The index-three bank lattice consists exactly of integer W with W1+W2 divisible by 3.

With current-layer old parents a(1,2)+b(1,−1), the exact increment additionally has −3ax+3ay−6by. Maximizing over real x,y gives 6/5(a²−ab+4b²). The quadratic matrix [[2,−3/2],[−3/2,3]] has determinant 15/4 and is positive definite. The extension to the exterior preserves its complete block polarization: its constant includes δ(O,Hext)+β(aL0+bL2,Hext), and its two linear coefficients include every bank/exterior cross term. The Schur contribution is (3l1²+3l1l2+2l2²)/15. An exact polynomial check verifies the full identity before the inequality is taken.

Hostile control: one completed B pair by itself contributes −3. Adding exterior old vector (−10,0) and no exterior new vector changes the total increment to +27. Therefore dropping a completed pair or reusing its isolated negative credit is invalid. The proposed gate does not drop that interaction. The 4/9 bound applies only to completed-bank energy, not full native energy.

## Native source proof and guards

For the four child formulas, the first division produces respectively B(u), B²(u), B²(u), B³(u), where B(u)=4u+1. The initial divisions use exponents 2,2,1,1. Primitive u ends the stripping process; every positive B-iterate before it is 5 modulo 8. The actual full clocks are therefore 4,6,5,7, without clock reduction modulo any register.

The A and C parent classes are distinct upper-half residue classes; D has another disjoint upper-half residue class; B occupies the lower interval. All parents are below 2N; every child is above 8N and below 256N. Children associated with distinct parents cannot coincide because P is a function. C's two children are distinct explicitly. Consequently all original source endpoints across all four families are disjoint; no later shared terminal is counted as new source ownership.

Since max(V)=63665, every child is outside V when N≥7959. Below this threshold one must discard any whole cluster with a child already in V before using the stopped-arrival statement. This removes at most 627 clusters because original child sources are disjoint. The audit of literal source formulas at small N does not waive this stopping guard.

When the parent has finite first arrival (t,k,d), a child with edge clock c has first arrival (t,k+c,d+1). Thus C's two children arrive at the same depth and terminal but have different actual clocks. Their signed charge sum is L1+U2=U1. If the parent has no finite first arrival, neither child has one; no universal entry is needed for this conditional statement.

Exact residue-floor counts match literal family counts at every audited N. For 144 dividing N each family has N/48 parents. A+C and B+D each supply N/24 clusters; all four together supply N/12 clusters with 3N/16 distinct original endpoints. This endpoint cardinality is not itself a norm budget or timely-arrival estimate. The doubling is of relation clusters, not a doubling of certified global energy credit.

## Required scope and editorial caution

The exterior-plus-boundary inequality at every horizon remains the open obligation. Congruence density alone supplies no bound on current-layer arrival concentrations or the exterior coupling. The user-reported bounds 520738,1039055,3848613 are not recovered or used here.

Minor wording change recommended in CONTRACTION_EXTENSION.md, verification paragraph: replace “proves the arithmetic identities by exact checks” with “corroborates the proved arithmetic identities by finite exact checks.” The all-scale proof rests on the displayed affine inverse formulas and congruence/interval reasoning, not on a finite test range.
