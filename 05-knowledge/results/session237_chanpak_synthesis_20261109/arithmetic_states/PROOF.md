# Exact arithmetic memory and simultaneous modular projections of fixed-port tails

Status: PROVED (proof below), with independent FINITE-EXACT checks. Independent audit PASS; THM-5091. Chan–Pak, HYP-9480(C), and coverage by the full matrix orbit remain OPEN.

## Inheritance, scope, and live concepts

The closest mechanism is the O/Z tail law

    H_s(m) = h + R m - d Q(m) - e (m mod 2),
    Q(0)=0, Q(2a)=Q(a), Q(2a+1)=Q(a)+a+(a mod 2).

This is **KNOWN**, from repository-root `session205/finite-oz-tail-bank/ARITHMETIC_REDUCTION.md`, restated in **THM-4888**, `01-canon/theorems/THM-4888-finite-oz-banks-mod3-law-involution-rigidity-residue-confinement.md`. The modular pumping mechanism is **KNOWN**, `session206/exact-modular-image/DYADIC_SURJECTIVITY.md`. Older references incorrectly prepend `05-knowledge/results/` to these two source paths; the actual tracked paths are repository-root paths. THM-5085, `01-canon/theorems/THM-5085-fixed-port-tails-in-every-base-and-the-ternary-bank-to-ten-million.md`, recovers these identities at x=2 and extends them to general base x. No binary priority is claimed for that extension.

The canonical hostile examples are impossible tournament counts 7 and 21. The corrected near miss is modular surjectivity mistaken for actual integer coverage. The underused sidecar is the joint residue of the parameter and its digit-pair statistic, not the scalar output alone.

Concept board: genuine tournament counts; digit-pair arithmetic; finite automata; simultaneous modular lattices; finite-bank coverage. Anchor: exact target-preserving arithmetic of physical tails. Niche: minimal state complexity. Wildcard: Smith-normal-form style simultaneous correlations. Session236's periodic drift and all-word tower rigidity are compatible background, not dependencies.

All results below concern **binary genuine Hamilton-path counts**. They assert nothing about ternary cycle-packing evaluations or the even-graph count problem.

## 1. The retained arithmetic state is fully reachable

Put q(m)=Q(m)/2. This is an integer: induction uses the fact that a+(a mod 2) is even. Fix any positive integer M. Every pair

    (m mod 2M, q(m) mod M)

occurs for some nonnegative integer m. In fact every pair occurs for arbitrarily large m.

### Proof, with the pumping step exposed

Write 2M=2^k u with u odd and k>=1. Choose b>=max(k,2) with 2^b=1 mod u, and write x=2^b-1 and c=2^b-2. For any a>=0, choose A>=k+bit_length(a), and set

    m(a,t) = a + x sum_{j=0}^{t-1} 2^[A+j(b+k)]  (t>=1).

Each high block consists of b consecutive one-bits, with k zero-bits separating blocks. Its contribution to Q is c. The lower a and a high block have cross contribution

    x sum_{i:bit_i(a)=1} 2^(B-i-1)

for a block beginning at B; it is divisible by 2M. Two high blocks separated in starting positions by D contribute x^2 2^(D-b), also divisible by 2M. There is no cross-block adjacent pair. These statements follow from the digit identity

    Q(m)=sum_i bit_i(m) floor(m/2^(i+1))
         + sum_i bit_i(m) bit_(i+1)(m).

Consequently m(a,t)=a mod 2M and Q(m(a,t))=Q(a)+tc mod 2M. Now c/2=2^(b-1)-1 is odd and coprime to u: modulo u it is -1/2. Thus gcd(c/2,M)=1. As t varies, q(m(a,t)) runs through all residues modulo M, keeping the first coordinate fixed. Choose a in each residue modulo 2M. Increasing t by multiples of M gives arbitrarily large witnesses. This is the old Session206 pumping construction with both coordinates retained; the simultaneous formulation is the present corollary.

## 2. Exact minimal finite-state arithmetic for the root

Let t=(h,R,e,d)=(1,2,0,1), the rooted directed edge with its source as port. Thus H_t(m)=1+2m-Q(m). Interpret the binary expansion of m **most significant bit first**, with optional leading zeros. The function

    f_M(m) = (H_t(m)-1)/2 mod M = m-q(m) mod M

has a deterministic finite automaton with output and initial state (0,0), whose states are all pairs (r,z) in Z/(2M) x Z/M, transitions

    0: (r,z) -> (2r mod 2M, z),
    1: (r,z) -> (2r+1 mod 2M, z+ceil(r/2) mod M),

and output (r-z) mod M. The formula ceil(r/2) mod M is well-defined on r mod 2M, since changing r by 2M changes the ceiling by M.

**Theorem.** For every M>=2 this automaton is minimal and has exactly **2M^2 states**. At M=1 the function is constant and its minimal automaton has one state.

### Proof

The recurrences of Q prove the transitions and output by induction. Section 1 proves that every state is reachable. It remains to prove pairwise distinction by future output; the proof gives a suffix of length at most one.

Take distinct states (r,z), (s,w), with 0<=r,s<2M. If their current outputs differ, the empty suffix distinguishes them. Otherwise r-z=s-w mod M. After a zero the difference of outputs is r-s mod M. Thus a zero distinguishes them unless r=s mod M. If r=s as integers, equality of current outputs forces z=w, contrary to distinctness. The only remaining case is |r-s|=M. After a one the output is current_output+floor(r/2)+1. The difference is therefore, up to sign, floor((a+M)/2)-floor(a/2), where a=min(r,s). For M>=2 this is either floor(M/2) or ceil(M/2), lying strictly between 0 and M. So a one distinguishes the pair. Every pair is distinguishable and every state reachable, proving minimality.

This is an exact Myhill–Nerode argument, not an inference from the finite tests. The arithmetic state is not merely one convenient oversized carrier: no smaller deterministic MSB-first output automaton computes the same full residue function.

## 3. Exact simultaneous residue image of any finite seed list

Let s_i=(h_i,R_i,e_i,d_i), i=1,...,l, be fixed physical seeds covered by the binary tail law. Their h_i are odd and R_i-e_i are even. Use one shared tail parameter m. Set

    Y_i(m)=(H_(s_i)(m)-1)/2 mod M,
    A_i(b)=(h_i-1+(R_i-e_i)b)/2, b in {0,1}.

Then the **exact joint image** is

    union_{b=0,1} { A(b) + j R - z d mod M : j,z in Z/M }.

Every displayed vector occurs for arbitrarily large parameters m. The proof writes m=2j+b and applies Section 1 to choose j and q(m)=z independently modulo M. Necessity is direct substitution; sufficiency requires that reachability theorem.

Thus simultaneous residues are two translates of a two-generator module, independently of the number of seeds. Standard integer row/column operations (Smith normal form) can describe this finite module, but no external normal-form theorem is needed for the stated formula.

Two examples using physical seeds already recorded in THM-4888/HYP-9480:

* Root t and s=(5,6,2,3): Y_s=2+3Y_t-b mod M. The joint image is exactly these two affine lines, even though the root alone has full modular image.
* Root t and s=(9,11,3,5): for fixed b, the coefficient matrix of (j,z) is [[2,-1],[11,-5]], with determinant 1. Therefore every pair of residues occurs already within either fixed parity of m, for every M.

These are **synchronized-parameter correlations**. They are not correlations between the image sets with independently varying parameters, so they neither prove nor refute HYP-9480(B)'s image-union heuristic. That distinction is load-bearing.

## 4. What the connection preserves and loses

Source: actual rooted tournaments followed by O/Z letters. Target: modular finite automata and finite affine modules. Map: retain (m mod 2M,Q(m)/2 mod M), then read the exact count residues. Preserved predicate: Hamilton-path-count congruence, including simultaneous congruences for a shared tail. Lost information: integer magnitude, equality to a specified integer, collisions between unrelated parameters, and the actual labelled tournament. Sidecar for construction: the input digit word and the rooted seed.

The automaton reads arithmetic digits MSB first, but the physical tail word appends letters in increasing bit position: bit i=0 means a sink at step i+1 and bit i=1 means a bridge back to the fixed port. Confusing those two orders is a hostile failure. The checker constructs literal tournaments in the physical order.

Consequences and boundaries:

* Fixed-modulus residue checking is decidable by a fully explicit automaton; its exact root memory grows quadratically with M.
* This is not a lower bound on arbitrary algorithms, proofs, or unrestricted tournament encodings. It is a lower bound only for the specified deterministic MSB-first residue-output model.
* Full modular projection does not imply actual coverage. Counts 7 and 21 remain absent despite surjectivity at every fixed modulus.
* No positive-density miss theorem for arbitrary finite banks follows. The scalar projections were already characterized in Session206; the present contribution is the exact minimal memory and joint module formulation.
* Generalizations to x=3 would require retaining the base-x port state and its integrality conditions; they are not asserted here.

## 5. Reproduction and hostile controls

Run from the repository root:

    python3 04-computation/experiments/session237_arithmetic_states/verify.py

The script has no repository imports and uses explicit exceptions, so checks survive python -O. It independently explores all modular transitions for M=2,...,40 (exactly 2M^2 reachable states each), verifies distinguishing output triples, checks the digit recurrence for all m<65536 in six moduli, and constructs all 256 root tails m<256 as literal tournaments. Directed Hamilton paths are counted by subset endpoint DP, independently of the tail formula. It checks the simultaneous image and both examples for M=2,...,20. It also constructs 2990 literal block-pumped integers for M=2,...,16, verifying both preserved parameter residues and complete q residues independently of automaton BFS. Normal and optimized-Python outputs are byte-identical.

Hostile controls reject (a) retaining m only modulo M, (b) retaining only the current scalar output, and (c) reading modular surjectivity as actual attainment of 7 or 21. The concrete collapsed states (r,z)=(0,0),(4,0) at M=4 have identical current outputs and identical r modulo M but different outputs after a one.

Results: `04-computation/experiments/session237_arithmetic_states/results.json`. No large-count coverage census, asymptotic density claim, or conjecture promotion is made.
