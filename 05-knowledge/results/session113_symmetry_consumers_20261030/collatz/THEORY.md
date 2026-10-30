# A consumer-specific quotient for sharp Collatz content bounds

Status: exact certificate theorem at E=6,7,8 (independent receipts); the all-E sharp content inequality remains OPEN, while the universal epoch-envelope budget is REFUTED below. The quotient is an upper-envelope observer, NOT an exact future-return observer, and not a source reset.

Inherited source: [fifth-layer content and epoch budget](../../layer_lifts_cyclic_intersections_20261030/COLLATZ.md). [Independent audit](../audit_collatz/AUDIT.md) verifies the proof bridge and finite certificates; [first-return refinement](../uniform_content/THEORY.md) records the surviving next target.

## 1. Inherited objects and closed quotient
For a positive valuation word, P=3^d,Q=2^A,B are its exact affine coefficients. Put z=(B+Q-P)/Q and epsilon=0 if m=P/Q is a square modulo5, otherwise1. Retain the original prefix clock A; no actual starting integer is replaced. At precision q=5^E, the pair (epsilon,z mod q) has 2q states, with

    z'=(3z+2^a-2)/2^a mod q,
    epsilon'=epsilon+a+1 mod2.

The identity starts at (0,0). Letters1 leave that state unchanged, so any initial run of ones can be separated without remembering its length in the state (its clock contribution is still paid explicitly). For arbitrary affine states this quotient is closed; its fibres have size phi(q)/2=2*5^(E-1).

It preserves the defect height truncated at E and square phase. It does not preserve m=1 and therefore cannot answer exact future identity queries. This is why the inherited 4*5^(2E-1) lower bound does not apply to this upper-bound consumer.

## 2. Epoch upper envelope
Write h(s)=min(v5(z),E), with h=E at z=0 modq. Let k(a)=min(v5(2^a-2),E), for a>1. Set

    D=8q, N=5q-5, lambda=N/D=(5/8)(1-5^-E).

An epoch consists of a letter a>1 followed by L ones. Its starting state is t. If epsilon(t)=1, the epoch has no identity returns even modulo5, so no content charge. If epsilon(t)=0, the content charge is at most

    h(t) + lambda_h(t) L,
    lambda_h=(5/8)(1-5^-h).

Indeed at level e<=h the multiplier runs through powers of3/2 of order 2*5^(e-1). There are at most floor(L/[2*5^(e-1)])+1 returns among the L+1 positions. Sum over levels. No higher level can contribute. The bound is deliberately an envelope, not an equality; it retains the common actual quotient state across all levels.

The initial ones contribute at most lambda times their length. Define integral transition costs after the first non-one letter:

    c(s,1) = -(5q/5^h(s)-5) if epsilon(s)=0, else -N;
    c(s,a) = D*max(h(t),k(a))*[epsilon(t)=0] - Na, a>1,
    t=T_a(s).

Because max(h,k)>=h, summing these costs from the first non-one letter gives an upper envelope C satisfying

    D S_E(word) - N A <= C.

The stronger max(h,k) expression exactly corresponds to the inherited active jump debt accounting; using h alone would also be a valid, tighter envelope. All costs are integers. Initial ones have been paid by N times their original clock, not discarded from A.

## 3. Finite certificate and infinite letter tail
A vector V on the 2q states certifies the claim if

    V(T_a(0,0)) >= c((0,0),a) for every a>1,
    V(T_a(s)) >= V(s)+c(s,a) for every state and a>=1,
    max V <=0.

Induction gives C<=V(final)<=0, including words ending in the middle of an epoch. All-ones words were treated separately. Therefore S_E<=lambda A for every finite positive valuation word.

For letters a>cap>=1, c(s,a)<=DE-Na, since both h and k are truncated at E. Thus it is enough to verify

    maxV+DE-N(cap+1)<=minV,
    DE-N(cap+1)<=minV.

These inequalities cover all omitted positive letters, without assuming a period or bounding actual valuations. The second is the entry condition. The first handles every possible preceding quotient state.

The independently checked fixed-precision results are:

| E | sharp coefficient | quotient states | full exact-return states |
|---|---:|---:|---:|
|6|1953/3125|31,250|195,312,500|
|7|19531/31250|156,250|4,882,812,500|
|8|48828/78125|781,250|122,070,312,500|

The C++ producer uses integer max-plus relaxation and exports V. The Python/NumPy auditor independently computes transitions using modular inverses, verifies every exported edge for a=1,...,64, all entry inequalities, and the two tails. Receipts give byte hashes and exact extrema. The vector need not be trusted as optimally relaxed: its checked inequalities alone prove the theorem.

## 4. Sharpness and source scope
For t=2*5^(E-1) consecutive ones, the content is 3^i-2^i at prefix i; v5 is0 for odd i and1+v5(i/2) for even i. The truncated sum is (5^E-1)/4, exactly lambda*t. The positive source 2^(t+1)-1 realizes the finite word. Thus the coefficient is sharp even among source-realizable words. No infinite all-ones positive orbit is claimed.

The theorem is about fixed-prime content and every finite word. It neither bounds all-prime content nor proves descent or convergence of a fixed positive Collatz source. The all-E sharp inequality remains open. Section6 refutes all-E feasibility of this particular epoch envelope; finite receipts are not induction.

## 5. Reproduction

    g++ -O2 budget.cpp -o /tmp/collatz_budget
    /tmp/collatz_budget 6 64 /tmp/collatz_V6.bin
    python audit.py 6 /tmp/collatz_V6.bin --output audit6.json

Replace6 by7 or8 for the other fixed-precision certificates. The producer's input is intended for these certified precisions; it is not a general arbitrary-size integer implementation.

## 6. Exact boundary: the epoch budget fails already at layer9

Let V=(3,1,3,1^6,3,1,3) and U=V 1^231 V. This is a finite positive valuation word of length255 and original clock271. All its epochs are active; their (letter,height,following ones) triples are

    (3,0,1),(3,2,6),(3,0,1),(3,5,231),
    (3,0,1),(3,2,6),(3,0,1),(3,9,0).

Each non-one letter has kappa=0. Hence the original untruncated budget J−B is exactly

    sum h - (5/8)*24 - (5/8)*(4 + 12/25 + 231/3125)
      = 769/5000 > 0.

This refutes the inherited universal J<=B hypothesis, including the tighter h-only epoch envelope. At E9 the sharp-baseline envelope excess is60112/390625>0. Therefore no vector V satisfying Section3's edge/entry/nonpositive inequalities can exist at E9: telescoping along this one word contradicts it. This is an algebraic obstruction to this envelope, not merely unsuccessful optimization.

Crucially actual untruncated content charge is154, below the E9 target132324151/781250 (about169.375), and below5A/8. Thus neither the sharp content conjecture nor Collatz is refuted. The loss is the independent +h allowance at each epoch: square phase alone does not retain the multiplier's precise first-return offset. A future refinement must price that offset, or replace this per-epoch upper envelope; changing the potential alone cannot repair its positive path.

`envelope_obstruction.py` independently iterates full integer P,Q,B, verifies the listed heights and actual gcd charges, checks both rational excesses, and constructs and verifies a positive integer source realizing the entire word. Its JSON receipt retains that source and endpoint. This closes the inherited budget hypothesis with a counterexample, while retaining the valid sharp E6–8 certificates.
