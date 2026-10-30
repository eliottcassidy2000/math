# Independent audit: sequence conservation and the strengthened weight frontier

**Verdict: PASS for the stated algebraic/transition claims.** This audit does not
certify the complete external elimination theorem, its affine-disperser assembly,
or an improved unrestricted circuit lower-bound record. The four-deletion AP2
lemma is a separately checked structural input. No producer implementation was
imported into the independent checker.

## 1. Independent conservation derivation

Use G=GP, B=AP2, A=AP1, P=P2+PE, Z=P0, H=PC, D=elementary deletions. Put
J=n-I_0, U=total usage decrease beyond the upper bounds, K=sum_GP(k-1),
and E=extra charged-object deletions beyond the guaranteed row counts. Usage
loss and gate loss are different slack variables. Every physical deletion must
be allocated once; the algebra cannot repair double counting of shared gates.

The row dimensions imply h=a_T=n-G-A-B. Pair births and removals give
q_T=G-2B-A-P-Z-H. Usage gives I_T=n-J-2G-A+P-U. Substitution yields exactly

    h=A+B+q_T+Z+H+U+I_T+J.

All summands are nonnegative, so A+B<=h and G=n-h-A-B>=n-2h. Also
P=n-h-3B-2A-q_T-Z-H. The final subtraction beyond h is at most 3h by the
same budget identity. Thus P>=n-4h. If h<=D_0 these imply A+B<=D_0,
G>=n-2D_0 and P>=n-4D_0. No independence of local worst cases is assumed.
These statements apply to every completed prefix with its own h, not just the
last state. When h is small relative to n, AP1/AP2 are sparse but GP/cleanup
steps are forced to be numerous.

## 2. Gate-slack identity verified independently by potential changes

For strengthened W=s+4I+2q, the baseline potential losses are:

- G: 8+(k-1), A: 9, B: 8;
- P: 0, Z: 2, H: 3, D: 1.

Additional usage loss U costs 4U and additional deletions cost E. Therefore

    W_0-W_T=8(G+A+B)+A+2Z+3H+D+K+4U+E.

Using W_0=s_0+4(n-J), W_T=s_T+4I_T+2q_T, and G+A+B=n-h gives

    s_0=4n-8h+s_T+4(I_T+J+U)+2q_T
         +A+2Z+3H+D+K+E.

This matches the producer exactly. The older AP2=3 table similarly gives
4n-9h with coefficients 5 on I_T+J+U, 3 on q_T, 2 on A, 3 on Z and 4 on H.
Every remainder coefficient is strictly positive, so the producer's aggregate
zero-surplus conditions are necessary and sufficient. Eliminating the remaining
counts then yields B=h,G=n-2h,P=n-4h. The supplied prefix-feasible equality
words show sharpness only of scalar inventories. They do not prove realizability
by a reduced DAG, an eligible function or a minimum circuit.

## 3. Independent frontier elimination

For g=1, let c=t-i and b=t. The inequalities become i=b-c>=0 and

    max(0,b-c-2,(b-4)/2,c-3) <= p <= b-2c+2.

Comparing the four lower endpoints with the upper endpoint yields b>=2c-2,
c<=4, b>=4c-8 and b>=3c-5. The last is implied by b>=2c-2 below c=3 and
by b>=4c-8 above c=3. Along with b>=c and c>=0 this proves the exact
projected feasible region claimed in the producer, including sufficiency.
Its lower envelope changes at c=2 and c=3. At c=4,b=8, both endpoints for
p equal 2, forcing i=4,t=8 and hence the unique normalized optimum (1,4,2,8).
The coefficient-four ceiling still follows from P and one-child GP alone.

The conditional telescope in FRONTIER.md has the correct sign for the initial
unused-input term +4(n-I_0), and retains terminal s_T,I_T,q_T rather than
silently setting them to zero. It is a valid implication under the named
sequence invariants/coverage, not a proof of those hypotheses.

## 4. Independent computational checks

Run `python3 check_independent.py` and `python3 -O check_independent.py` from
this folder. Both produce the archived `receipt.json` byte for byte. Assertions
are explicit exceptions and survive optimization.

The checker independently generates 6,000 abstract transition histories, 420,000
prefixes in total, with nonzero initial unused inputs and random extra usage/gate
losses. It checks the three conserved inventories, endpoint budget, frequency
bounds and both expanded gate identities at every prefix. These are scalar
regression histories, not claims of actual circuit realizations.

It also enumerates 109,265 rational weight triples directly in weight coordinates
and constructs 1,105 feasible coefficient/penalty pairs from the projected region.
The mathematical statements over all integers and all real nonnegative weights
follow from the preceding exact proofs, not the bounded enumeration.
