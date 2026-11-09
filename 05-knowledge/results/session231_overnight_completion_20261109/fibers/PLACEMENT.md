# Exact placement law for the twelve-bit intervals

Status: PROVED in the stated scope, independently audited in ../fiber_audit/MINIMALITY_PLACEMENT_AUDIT.md. No unrestricted completeness or global minimum claim.

## Context stability

For every k, all choices j in the short interval theorem share h,q,A,G,S,E and differ in T by2j. Consequently any common ordinal left/right SINGLETON context preserves their exact consecutive-odd window. With non-singleton contexts the step is multiplied by the product of their Hamilton counts, so the singleton restriction is substantive. For strongness use a selected first prefix vertex, or leave the prefix empty so that the original selected leading vertex remains first.

Let P be the prefix cut integer (a bits, first bit1 when a>0), R the suffix cut integer (b bits), and rev_b(R) its b-bit reversal. At fixed P,R,a,b, the baseline N_k has the form

    N_k=alpha_context·4096^k+(5870/91)k+gamma_context,

where alpha_context,gamma_context are rational constants determined by that context. This follows because q,A,G are affine in4096^k, h,S,E are constant, and the only linear k term enters through−T, whose coefficient is5870/91 in the explicit bare baseline. Therefore

    N_(k+2)=4097 N_(k+1)−4096 N_k−264150.                (2)

The interval is {N_k−2j:0≤j≤k}. Formula(2) is all-order for each FIXED context; allowing context to depend onk introduces additional choices and is precisely the remaining selector problem.

## Modular and quantitative boundaries

Modulo45, the common short-core interface has A=1,G=13,q=S=E=h=0. The inherited arbitrary binary context formula then gives

    N_k ≡ 2^(b+1)−rev_b(R)−13P+7+5k (mod45).

In particular D_k≡9+5k modulo45 for empty contexts. Once k≥44, each individual window meets every residue modulo45, but this says nothing about its position on the number line. It is not an occupancy theorem.

All these strong constructions have the first selected deletion vertex p and the final closure vertex v as fixed flanks. For every ordered path cover(L,R) of the intervening deletion graph, the sequence L,v,r,p,R is a directed Hamilton path. This injects all those covers into the counted paths. The middle cover count is

    5760·4096^k·2^(a+b),

so EVERY count in every displayed interval satisfies this lower bound. Thus a target N can use only

    k ≤ floor(log_4096(N/5760)),

when N≥5760. The repair width is at most logarithmic in the target size. For the bare sequence the exponential baseline growth dwarfs width2k and adjacent k-windows do not overlap. Arbitrary two-sided contexts could still provide enough moving baselines; their occupancy has not been proved or refuted here.

The exact remaining question is: for each sufficiently large odd N, do some legal P,R,a,b,k place N in[N_k−2k,N_k]? A proof needs actual baseline placement, not just interval length, gcd-one, residue coverage, or count entropy.

## Additional arithmetic controls

For the earlier unsaturated h95 fixed-step190 family, appending its Y word acts linearly on(h,q,A,G,T) modulo190 because Y starts0. The exact three-step action returns the physical initial five-coordinate vector modulo190. Thus the baseline residues are149,169,119 cyclically in k, and each is coprime to190. All displayed progression values are consequently coprime to190. This is not a prime-value theorem.

For any fixed moving carrier, d=gcd(h,q,A_0,G_0) is odd and remains a divisor of all reachable A,G, because dividing by2 is invertible modulo odd d. Every edge cocycle is divisible byd; every cyclic edge is even. Hence every physical loop gain is divisible by2d, in addition to2h/gcd(h,q). Their least common multiple is a valid stronger obstruction. Content-one saturated interfaces are therefore the useful targets for a unit-gain search; saturation alone is insufficient.

`verify_contexts.py` checks (2), the modular law, the lower bound, and all positions for every legal prefix/suffix of lengths≤4 andk≤4. The all-order statements follow from the displayed algebra, not that finite test.
