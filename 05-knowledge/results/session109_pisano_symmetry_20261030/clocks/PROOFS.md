# Exact observer clocks for Fibonacci and cat-map cyclic states

Session 109. All statements below are proved elementary deductions. They are a synthesis of classical recurrence/finite-group facts, **not a priority claim** and **not a proof of HYP-9372**. Set Q=((1,1),(1,0)), A=Q²=((2,1),(1,1)); all rings are R=Z/N, N≥2. Let B denote either Q or A, m=ord_N(B), and H=<B>.

## 1. Observation is a quotient of a cyclic history, not necessarily a group representation

For any function f:H→Y define its **whole-history clock** p_f as the least positive t such that f(B^{j+t})=f(B^j) for every integer j. The shifts satisfying this condition form a subgroup K_f of Z/m. Consequently p_f|m and |K_f|=m/p_f. Two states B^i,B^j give identical future histories iff i−j∈K_f, so every history fiber has size m/p_f. A coincidence f(B^t)=f(I) alone does not establish a period.

If f is an equivariant observation of a group action, this specializes to the stabilizer quotient. If f is a representation followed by projectivization, p_f=m/|H∩ker(projective representation)|. No claim is made that this kernel equals the matrix center: central matrices may act nonscalarly, and noncentral matrices may lie in the representation's projective kernel. This retains exactly the correction in HYP-9372.

For several simultaneous observations, the joint clock is the lcm of their clocks: its shift subgroup is their intersection. For pairwise coprime moduli, any observation retaining its entire tuple of local observations has clock equal to the lcm of local clocks. This follows directly from CRT and is not true for an arbitrary lossy postprocessing of the tuple.

## 2. Marked, scalar, and signed clocks

For both Q and A, v=e₁ is a cyclic marked vector: det(v,Bv)=1. If B^t v=v, commutation also gives B^t Bv=Bv, hence B^t=I. Thus the marked orbit B^j e₁ has full matrix clock m. More generally, det(v,Bv) a unit is sufficient for faithfulness of a marked vector over R. Being merely nonzero is not sufficient over a ring.

Let α(N)=min{k>0:F_k=0 modN}. The identity Q^k=F_k Q+F_{k−1}I gives the exact scalar clock for Q: s_Q=α(N). Its scalar multiplier λ=F_{α−1} is a unit, and π(N)=α(N)ord(λ). Since λ²=(−1)^α, ord λ divides4. For N≥3, odd α implies λ²=−1≠1 and ordλ=4; even α gives order1 or2. Scalar powers of A occur when α|2t, so s_A=α/gcd(α,2). The scalar loss for A is only1 or2; the scalar loss for Q is1,2,or4. These are classical restricted-period facts.

Under CRT, being an arbitrary unit scalar is a local condition, so s_B(N)=lcm_i s_B(N_i). Requiring one **global sign** ±I is more restrictive: local signs must synchronize. Example: ord_15(A)=20, first arbitrary scalar power A^10=11I mod15, and −I never occurs. Indeed local first scalar times are2 mod3 and5 mod5; the local −I conditions t≡2 mod4 and t≡5 mod10 are incompatible. Thus scalar clock10 versus signed clock20. Quotienting by all scalars must not be silently replaced by quotienting by {±I}.

## 3. Exact trace-history theorem and CRT masking

Write τ=tr B (1 for Q, 3 for A). The algebra R[B] is free on {I,B}, since B has off-diagonal entry1. For X=aI+bB, its two adjacent trace observations are

- Q: (trX,trQX)=(2a+b,a+3b);
- A: (trX,trAX)=(2a+3b,3a+7b).

Both coefficient matrices have determinant5 and an entry1 (Q) or relatively prime entries2,3 (A), hence Smith normal form diag(1,5). Their kernels over R have cardinality gcd(N,5). Cayley–Hamilton shows that these two traces determine the whole trace history tr(B^kX), k≥0: all later terms obey the same order-two recurrence. The trace clock t_B therefore equals m unless 5|N; in all cases m/t_B≤gcd(N,5), because the cyclic invisible-shift subgroup injects via B^t↦B^t−I into this kernel. The exact local and CRT calculation below strengthens this bound to m/t_B∈{1,5}. Thus the only possible trace-clock loss is a factor5 at the discriminant prime5, **not** the factor2 or4 of a scalar/projective clock.

The exact local formulas, for e≥1, are

| B | matrix clock at 5^e | trace clock at 5^e |
|---|---:|---:|
| Q | 4·5^e | 4·5^{e−1} |
| A | 2·5^e | 2·5^{e−1} |

Proof of the first row without any prime-power conjecture: Q mod5 has order20, checked from Q^4=I+D₀ where D₀ mod5=[[4,3],[3,1]] is nonzero square-zero, and its eigenvalue is3 of order4. Q^20=I+5D₁, D₁=[[2189,1353],[1353,836]], and D₁ mod5=D₀ mod5. Binomial expansion inductively gives
Q^{4·5^k}=I+5^kD_k for k≥1, D_k≡D₁ mod5.
Indeed (I+5^kD_k)^5≡I+5^{k+1}D_k mod5^{k+2} for k≥1. Thus Q^20 has exact order5^{e−1} modulo5^e, and Q has order20·5^{e−1}. Both trD₀ and tr(QD₀) vanish mod5. Therefore Q^{4·5^{e−1}}−I lies in the trace-pairing kernel (also when e=1), and this power has order5. The kernel has size5, proving the claimed trace clock. For A, its trace sequence is the even subsequence of Q's: 2·5^{e−1} is a trace period and its matrix clock is half Q's, while the same kernel bound makes this period exact.

Write N=5^e u with gcd(u,5)=1. Set ord_1(B)=1. CRT gives

    t_B(N) = ord_u(B)                        if e=0,
    t_Q(N) = lcm(ord_u(Q),4·5^{e−1})          if e≥1,
    t_A(N) = lcm(ord_u(A),2·5^{e−1})          if e≥1.

The local defect is masked globally iff v₅(ord_u(B))≥e. For instance at N=55, A has local matrix orders10 at5 and5 at11, so global matrix and trace clocks are both10. The factor5 loss at5 alone disappears when the 11-component is retained. No assertion about non-5 prime-power lifting or Wall–Sun–Sun primes is needed.

## 4. A minimal five-valued sidecar recovers the missing phase

For a cyclic state X=B^j, put T=trX and U=tr(BX). Then in either case

    5b = 2U−τT modN.

Thus the trace history determines b modulo N/gcd(N,5). If 5|N, write N=5^e u. Choose the canonical representative b₅ of b modulo5^e and add the **one-digit sidecar**

    d(X) = floor(b₅ / 5^{e−1}) ∈ {0,1,2,3,4}.

The history determines b₅ modulo5^{e−1} and b modulo u; d recovers b₅, hence b by CRT. Then recover a=U−3b for Q, or a=U−T−4b for A. This proves injectivity of the combined observation (T,U,d) on the entire algebra R[B], not just its cyclic subgroup. If trace history actually has five-state fibers on H, any ONE contemporaneous snapshot sidecar that distinguishes each fiber must have at least five possible outputs; this sidecar attains that snapshot lower bound. This is not a lower bound on the alphabet of a sensor whose temporal history is retained: a binary sensor can recover phase over several samples. See `../SYNTHESIS.md` §4 and `../synthesis/check.py` for the five-phase cyclic code 00011, which distinguishes phases with three binary samples at times 0,t,2t; two binary samples cannot distinguish five phases. If CRT masking leaves one-state history fibers on H, no sidecar is necessary for cyclic states (though the algebraic trace pairing still has a five-element kernel).

The sidecar is a choice of residue representative, not a canonical group character; the injectivity statement requires no homomorphism or cocycle claim. At e=1 it is simply the off-diagonal entry mod5. A full marked vector also restores the clock, but carries more information than this minimal sidecar.

## 5. Hostile examples and scope of traces

1. Over F₅, A and −I have the same determinant and the same traces of **every** power: tr(A^j)=2(−1)^j. Their matrix orders are10 and2. The loss is invisible even to the full unmarked trace history. A marked vector distinguishes them immediately.
2. Over Z/N, I and U=[[1,1],[0,1]] have trace2 at every power but orders1 and N. Over Z, U has infinite order despite that same constant trace history. Therefore no statement that trace histories generally determine matrix order is valid. The positive theorem above is specifically about shifts within the fixed cyclic algebra/orbit of Q or A, with its nondegenerate trace pairing away from5.
3. Even inside <A> over F₅, the five states A^{j+2r}, r=0,...,4, have identical future trace histories. Their off-diagonal residues distinguish them. Distinct scalar and trace clocks are genuinely incomparable: at N=5, A has matrix clock10, scalar clock5, trace clock2.
4. THM-4763 concerns characteristic-zero logarithmic representations and Lucas-weighted traces; the present finite-ring theorem is not a deduction that such representations have finite order. Reduction modulo N changes the observer and coefficient ring. Neither theorem proves the general modular-data saturation hypothesis HYP-9372.

## 6. Provenance and reproducibility

Read: THM-4762, THM-4763, and the corrected HYP-9372. The elementary group/kernel framework and sidecar formulation are the session's scoped synthesis. Fibonacci matrix/restricted period/CRT facts and Lucas periodicity are classical; no broad novelty is asserted.

Primary sources browsed 2026-10-30:
- D. D. Wall, *Fibonacci Series Modulo m*, American Mathematical Monthly 67 (1960), 525–532: https://www.tandfonline.com/doi/abs/10.1080/00029890.1960.11989541 (bibliographic landing page; not used as substitute for the explicit proofs above).
- D. W. Robinson, *The Fibonacci Matrix Modulo m*, Fibonacci Quarterly 1(2) (1963), 29–36: https://www.mathstat.dal.ca/FQ/Scanned/1-2/robinson.pdf . The primary paper explicitly develops period, restricted period, scalar multiplier, their relation, and lcm identities (pp.29–30). These are prior art, not new results.
- M. Fiebig, aBa Mbirika, J. Spilker, *Period patterns, entry points, and orders in the Lucas sequences: theory and applications*, accepted manuscript v2 (2024): https://arxiv.org/abs/2408.14632 . Its abstract and metadata identify prior study of both U and V sequence periods; cited as scope context, not as verification of an uninspected theorem.

Run `python3 check_clocks.py`. It uses exact integer modular multiplication, direct orbit enumeration, and full-word period tests, independent of the proof formulas. `results.json` records the successful finite domain. These computations verify finite cases only; the proofs carry the general assertions.

## 7. Unified ramified metallic trace/orbit theorem (p≥5)

Let p≥5, a∈Z_p, D=a²+4 with v_p(D)=1, and O=Z_p[μ] with μ²=aμ+1. Put δ=2μ−a, so δ²=D. Since2 is a unit, O=Z_p[δ] is a ramified quadratic DVR, normalized by v_δ(δ)=1 and v_δ(p)=2. Work in O/p^kO=O/δ^{2k}O. The involution sends μ to a−μ and δ to−δ; let Tr denote its trace. The trace-pairing matrix on (1,μ) is

    [[2,a],[a,a²+2]],  determinant D.

Its Smith form over Z_p is diag(1,D). Therefore the history radical has size p and is exactly

    J_k = δ^{2k−1}O / δ^{2k}O.

For a direct verification of containment, δ^{2k−1}=D^{k−1}δ, Trδ=0 and Tr(μδ)=D. These traces vanish modp^k; the ideal has cardinality p, equal to the kernel size. Thus no unproved use of trace-dual/different machinery is needed (though J_k is precisely the inverse-different explanation).

Let g=μ⁴. The exact identity g−1=aμ²δ gives v_δ(g−1)=1 because a and μ are units. For every positive integer d,

    v_δ(g^d−1) = 1+2v_p(d).

Proof: a prime-to-p power preserves positive valuation by the first binomial term. For z with s=v_δ(z)≥1, the term pz in (1+z)^p−1 has valuation s+2; the terms binom(p,i)z^i, 2≤i<p, have valuation ≥2+2s>s+2, and z^p has valuation ps>s+2 since p≥5. Hence raising to p increases the valuation exactly by2. Iteration proves the formula. This step is why p=3 is excluded without a separate argument.

Take a nonzero residue x and define its remaining depth r=2k−v_δ(x), 1≤r≤2k. Its full g-orbit clock and scalar-trace-history clock are respectively

    p^floor(r/2),       p^floor((r−1)/2).

Indeed full return means x(g^d−1)∈δ^{2k}O, while trace-history return means x(g^d−1)∈J_k, equivalently valuation≥2k−1. Substitution of the valuation formula gives these two least powers of p. For scalar trace along g, this really is the entire history radical: g=(a²+1)+a(a²+2)μ, and a(a²+2)≡−2a modp is a unit. Hence (1,g) and (1,μ) are equivalent bases of O; the two adjacent g-traces suffice just as above. The zero residue has both clocks1, separately.

**Parity is the exact criterion.** Even remaining depth r loses one factor p; odd depth loses none. In particular every unit (r=2k) has full g-clock p^k and trace clock p^{k−1}. On each fixed orbit with loss, there are p states per trace-history fiber, so a supplementary ONE contemporaneous snapshot value must have at least p possible outputs to distinguish that fiber. This minimality claim does not apply to the alphabet of a supplementary temporal sensor history; multiple binary samples can encode more than two phases. A constructive global sidecar is the top base-p digit of the μ-coefficient: the trace pair determines that coefficient modp^{k−1}, because D times it is determined and v_p(D)=1; that digit restores the coefficient, then Tr determines the constant coefficient since2 is a unit.

For Fibonacci a=1,p=5, g=Q⁴. This is the common ramified trace mechanism behind the missing factor5 above, now refined by valuation layer. It does **not** equate a trace-observation kernel with an affine-conjugacy cokernel for arbitrary strides; the latter has different stride-dependent Smith factors. The identity μ⁴−1=aμ²δ explains their meeting at this particular ramified stride, not a universal equality.

Additional independent check: `python3 -O check_metallic.py` verifies 240 deterministic nonzero valuation-layer representatives for (p,a)=(5,1),(13,3),(17,8),(29,5),(41,18), at k=1,2,3, including return and minimality for both clocks. It checks the ramification hypotheses before testing; `metallic_results.json` records the domain. This is supplementary finite evidence, not exhaustive testing at these larger moduli.
