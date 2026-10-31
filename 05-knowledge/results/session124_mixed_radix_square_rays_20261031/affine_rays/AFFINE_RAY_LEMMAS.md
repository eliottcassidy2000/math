# Uniform affine square-sum rays: arithmetic, parity, and Pell families

**Scope.** These are arithmetic identities and simplicity theorems for explicit infinite vertex sequences. No Hamiltonian-cycle extension follows without a separately verified compatible-pair supply and complete transition/consumer tables. In particular the four-prefix scale25 ray below is simple even though an attempted extension architecture can fail.

## 1. Complete two-prefix family in the bounded-offset regime

Let q>=3 be odd, A=q², and r>0 even. Set

`c=A+3-r²`, with `|c| <= (A-1)/2`.

For t>=0 put `Q_t=A^t`, `C_t=c*(Q_t-1)/(A-1)`, and

`H_t=(Q_t-C_t, 3Q_t+C_t)`.

Then concatenating H0,H1,... gives a positive, strictly increasing, simple square-sum ray. Its successive square roots are

`2*q^t, r*q^t` for t=0,1,... .

The finite prefixes `P_l=H0,...,H_l` have2l+2 vertices and satisfy

`P_0=(1,3)`,

`P_(l+1)=(1,3) concatenated with (A*p_j+(-1)^j*c)_j`,

where j is one-based.

### Proof

C_t is an integer geometric sum and `C_(t+1)=A*C_t+c`. The within-block sum is4Q_t. The boundary sum is

`(3Q_t+C_t)+(A*Q_t-C_(t+1))=(A+3-c)*Q_t=r²*Q_t`.

The bounded offset gives `|C_t|<=(Q_t-1)/2`, hence

`(Q_t+1)/2 <= H_t[1] <= (3Q_t-1)/2`,

`(5Q_t+1)/2 <= H_t[2] <= (7Q_t-1)/2`.

Both entries are positive and the second exceeds the first. The next block's first entry exceeds the current second by at least

`((A-7)*Q_t+2)/2 > 0`.

Thus the entire ray is strictly increasing, including negative-offset examples. The affine recurrence follows from C_(t+1)=A*C_t+c; the prefix length2 is even, so signs align after prepending.

### Exact parameter classification

For a fixed odd q, this bounded-offset family consists exactly of the positive even integers r satisfying

`(q²+7)/2 <= r² <= (3q²+5)/2`.

There are no further positivity or simplicity restrictions. In this even-r family, c is actually divisible by4 (because q²=1 modulo8 and an even square is0 or4 modulo8). The square-sum identities themselves do not require r even; that hypothesis selects the even-offset subfamily relevant to parity-compatible lifting architectures.

For comparison, if c is allowed to be any integer while A>=9 remains fixed, the displayed two-coordinate formula is positive for every t exactly when

`-3*(A-1) <= c <= A-1`.

It is strictly increasing across the entire concatenated sequence exactly when

`-(A-1) <= c <= A-5`.

To see this, write beta=c/(A-1). The coordinate growth coefficients are1-beta and3+beta. The within-block difference is `2*(1+beta)*Q_t-2*beta`; the boundary difference is `((A-3)-(A+1)*beta)*Q_t+2*beta`. Requiring these to remain positive at arbitrarily large t gives beta>=-1 and beta<=(A-3)/(A+1), with endpoint cases positive. For integer c, the upper bound is A-5 since `(A-1)*(A-3)/(A+1)=A-5+8/(A+1)`. These wider criteria do not classify every possible nonmonotone simple sequence outside the bounded-offset family.

## 2. Extreme offsets are two Pell-type families

At the positive extreme `c=(q²-1)/2`, the connector condition is

`q²-2r²=-7`,

and the blocks simplify uniformly to

`H_t=((Q_t+1)/2,(7Q_t-1)/2)`.

At the negative extreme `c=-(q²-1)/2`, it becomes

`3q²-2r²=-5`,

and

`H_t=((3Q_t-1)/2,(5Q_t+1)/2)`.

Examples are(q,r,c)=(5,4,12) and(3,4,-4), respectively. The first yields1,3,13,87,...; the second yields1,3,13,23,... .

There is a complete elementary orbit classification for positive q odd and r even:

- `q²-2r²=-7`: start from(1,2) or(5,4), and iterate
  `(q,r) -> (3q+4r,2q+3r)`.
- `3q²-2r²=-5`: start from(1,2) or(3,4), and iterate
  `(q,r) -> (5q+4r,6q+5r)`.

For the ray theorem, discard the initial q=1 solutions. These matrices correspond to the positive units3+2sqrt(2) and5+2sqrt(6), respectively (in the second case use x=3q in x²-6r²=-15).

**Completeness proof.** In the first equation, any solution with q>=9 is sent by the inverse matrix to positive integers `(3q-4r,3r-2q)` with smaller q; positivity of its first coordinate follows from q²>56. Parity and the equation are preserved. Among positive odd q<9, only(1,2),(5,4) occur. In the second equation, any solution with q>=7 descends to positive `(5q-4r,5r-6q)` with smaller q; the first coordinate is positive because q²>40. The only smaller odd-q solutions are(1,2),(3,4). In both cases deterministic descent proves completeness and separates the two orbits.

Thus the half-integer formulas are not isolated coincidences: their permitted dilation scales lie in explicitly classified Pell families.

## 3. Four-prefix identities: simplicity is an additional obligation

Take q odd, s>=4, r>0, with s and r of opposite parity. Put

`A=q²`, `c=A+s²-6-r²`, `H=(1,3,6,s²-6)`.

For Q_t=A^t and C_t=c*(Q_t-1)/(A-1), define

`H_t=(Q_t-C_t,3Q_t+C_t,6Q_t-C_t,(s²-6)*Q_t+C_t)`.

The concatenation is always an integer square-sum **walk**: its successive root multipliers are2,3,s,r, scaled by q^t. Here c is even, and the even prefix length4 preserves the alternating-coordinate signs. Positivity and distinctness are not automatic for arbitrary parameters.

A useful sufficient simplicity criterion is:

1. `|c|<=(A-1)/2`;
2. the four positive literals H are distinct;
3. no literal in H equals `A*m-c` or `A*m+c` for a positive integer m.

Indeed the affine lift is injective on any simple positive input sequence: opposite-sign collisions would require `A*(x-y)=±2c`, impossible under the bound unless c=0, in which case injectivity follows directly. The third condition keeps the new prefix disjoint from every lifted older vertex. Induction then proves simplicity of every nested prefix and hence of the infinite ray. This criterion is sufficient, not asserted necessary.

For a more general collision test, put d=A-1, epsilon_i=(-1)^i, and `a_i=d*H[i]+epsilon_i*c`. Equality between coordinates at levels t>=u is exactly

`(a_i*A^(t-u)-a_j)*A^u = c*(epsilon_i-epsilon_j)`.

Same-sign collisions require a coefficient ratio to be a power of A. Opposite-sign collisions require A^u to divide2c; when gcd(A,2c)=1 and c!=0, at least one level must be0. This gives an exact Diophantine test without assuming monotonicity.

### The proven scale49 arithmetic ray

For(q,s,r,c)=(7,5,8,4),

`H_t=((11Q+1)/12,(37Q-1)/12,(71Q+1)/12,(229Q-1)/12)`, `Q=49^t`.

All low literals1,3,6,19 lie below the least possible positive lifted value49-4=45. The sufficient criterion proves simplicity. The entries are strictly increasing within each block and across blocks. The root multipliers are2,3,5,8. Its Hamiltonian-extension theorem is a separate certificate result, not part of this arithmetic proof.

### The simple but nonmonotone scale25 arithmetic ray

For(q,s,r,c)=(5,6,7,6),

`H_t=((3Q+1)/4,(13Q-1)/4,(23Q+1)/4,(121Q-1)/4)`, `Q=25^t`.

The root multipliers are2,3,6,7. The first two blocks are

`(1,3,6,30), (19,81,144,756)`.

The potentially troublesome literal30 lies outside the low singleton band but is not a lifted positive vertex: lift residues are19 and6 modulo25, whereas30 has residue5. The literals1 and3 are also excluded by residue. Literal6 could only have the positive-sign preimage0, which is not allowed. Thus the sufficient criterion proves this ray simple for all depths.

Equivalently, after multiplying coordinates by4, same-sign collisions compare coefficients3,23 or13,121 times powers of25; none of those coefficients is divisible by5, so equal values force equal levels and coefficients. Opposite-sign collisions at two positive levels are impossible modulo25. The remaining level0 checks reduce to the disjoint first two blocks, with all later blocks already larger than the level0 maximum.

This ray is **not globally increasing**. It increases within each four-vertex block but has a descent at every boundary:

`H_(t+1)[1]-H_t[4]=(1-23*25^t)/2<0`.

Hence30→19 is the first member of a systematic pattern, not an exceptional numerical accident. There is no vertex collision despite those descents.

### A concrete warning against a generic simplicity claim

Take(q,s,r,c)=(5,9,10,0). All parity and bounded-offset hypotheses hold, but

`H0=(1,3,6,75)`, `H1=(25,75,150,1875)`.

The vertex75 repeats. The square-sum identities alone therefore yield only a walk, not a simple ray or an embeddable prescribed path.

## 4. Exact controls and scope

`test_identities.py` uses only integer and rational arithmetic. Normal and optimized runs are byte-identical, with **81,143 checks**:

- all2,643 admissible two-prefix parameter pairs for odd3<=q<=201, through eight blocks;
- the complete extreme-offset orbit lists through q=1,000,000;
- both specific four-prefix rays through21 blocks, including exact descents;
- residue-based prefix/lift disjointness and the generic75-collision counterexample.

These tests support, but do not replace, the universal proofs above. No computation here supplies a compatible Hamiltonian pair, validates a transition catalog, or proves a Hamiltonian-cycle extension for a new scale or a four-prefix candidate.

## 5. Arithmetic check for the later offset8 scale25 candidate

The same sufficient criterion also covers `(q,s,r,c)=(5,5,6,8)`, whose prefix is `(1,3,6,19)`. Its block formula is

`((2Q+1)/3,(10Q-1)/3,(17Q+1)/3,(58Q-1)/3)`, `Q=25^t`.

The root multipliers are2,3,5,6. Lifted positive vertices have residues17 or8 modulo25; none of the prefix literals has either residue. Thus the nested ray is simple. Each block increases internally, but every boundary descends by `(2-8Q)/3`, beginning19→17. This arithmetic fact is independent of whether any proposed relaxed-head transition catalog is complete or correct; no cycle-extension claim is made here.
