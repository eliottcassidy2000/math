# Carry read heads and the positive-state cost of exact contextual recovery

**Status: PROVED elementary interface theorems; FINITE-EXACT reproducible controls.** No external priority claim. These results specialize and combine previously established carry and observer mechanisms. They are not a Collatz convergence argument, a physical implementation of a new Collatz operation, or a new tournament Hamilton-count realization.

## 1. Source attribution and scope

Read sources:

- [source](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md), repository snapshot `1b91742ce55a202b73c80d56561ba56d00b503a6`: the affine word matrices, carry cocycle, faithful word representation, and blindness of all triangular word traces are inherited. The source itself attributes the underlying affine representation to classical Collatz literature; no independent literature audit is claimed here.
- [source](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/information_observers_collatz_20261014a.md), read from the same Git snapshot: the exact `(1,3)` / `(3,1)` source-guard example and distinction between statistical and verification sufficiency are inherited.
- [source](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/session33_continuation_states_20261023/metagraph-memory/NOTE.md), section 2a, in the same Git snapshot: the positive-state support-extremality proof for the two-time involution observer is inherited. It gives an actual order-five tournament example with 38 positive states versus 19 signed-linear dimensions.
- [S171 framework](../../session171_knots_information_20261103/framework/NOTE.md) and [S171 synthesis](../../session171_knots_information_20261103/synthesis/NOTE.md): contextual equivalence and continuation-stable quotients are inherited.
- [S171 information](../../session171_knots_information_20261103/information/NOTE.md): the constant-profile Reed–Muller tournament storage family and fixed 17/19 read head are inherited, including its independent finite verification. This packet does not repeat tournament enumeration.

New scoped deductions here are the exact lower-left-entry criterion for an external carry probe, the direct arbitrary-context positive-state lower bound, and their application to the specified Collatz and tournament examples. The general algebra and convexity proofs are elementary.

## 2. Exact triangular-context blindness and its smallest repair

Let k be a field, fix a,d in k, and consider the marked matrix family

    G_b = [[a,b],[0,d]],  b in k.

The matrix basis and its marking are part of the experiment. A trace context is a finite product of fixed upper-triangular matrices and copies of G_b, followed by trace. If a,d are nonzero, copies of G_b^(-1) may also be used. Fixed upper-triangular inverses are allowed whenever defined. Finite tuples of these traces and arbitrary subsequent processing are allowed.

**Theorem 1 (stable blindness and one-probe repair).**

1. Every such triangular trace context is independent of b. This includes arbitrary two-sided triangular continuations and arbitrarily many appearances of the unknown matrix.
2. For any fixed matrix P=(p_ij), the external context R_P(b)=tr(PG_b) is

       R_P(b)=p11 a+p21 b+p22 d.

   It separates every pair of distinct b if and only if p21 is nonzero, provided the source family has at least two b-values. When p21 is nonzero,

       b=(R_P(b)-p11 a-p22 d)/p21.

3. In particular J=[[0,-1],[1,0]] has determinant one and tr(JG_b)=b. It reads the previously invisible carry in one exact scalar response.

**Proof.** Sending an upper-triangular matrix to its ordered diagonal is an algebra homomorphism to k×k; it preserves products and defined inverses. The diagonal of G_b is (a,d), independent of b. Trace is the sum of the two diagonal entries, proving (1). Direct multiplication gives (2), and multiplication by J gives (3). The field condition ensures that a nonzero p21 is invertible. ∎

This is a classification within the declared probe family, not a universal minimum-complexity theorem. Over a ring with zero divisors, nonzero p21 alone need not imply injectivity. If a,d are not supplied, the single scalar need not identify the whole matrix. The measurement must access G_b or a retained equivalent: it cannot reconstruct b from discarded triangular trace data alone.

The distinction from Session171 is precise. Its one-vertex tournament attachment is an actual combinatorial operation with an independently counted Hamilton response. J is only an external matrix experiment here. It is not an odd Collatz step, a legal inverse step, or a construction of an integer orbit. The common mechanism is an enlarged observation alphabet that reaches a retained coordinate absent from the old statistic.

## 3. A literal Collatz instance and its source guard

For a positive integer exponent e define

    G_e=[[3,1],[0,2^e]],   U_e(n)=(3n+1)/2^e.

A chronological word w=(e1,...,el) has G_w=G_el ... G_e1 and

    G_w=[[3^l,B_w],[0,2^A]],   A=sum e_i.

This follows by multiplication. In particular,

    G_(1,3)=[[9,5],[0,16]],
    G_(3,1)=[[9,11],[0,16]].

Every old triangular trace experiment gives the same result on these two matrices, but J reads 5 versus 11. Their ordinary trace is 25 in both cases. The inherited carrier representation is faithful: retaining l,A,B_w determines the word; this packet only needs, and independently verifies, the displayed two-word instance.

For these words to describe the *accelerated odd* map at an integer source, the ordered valuations must be exact, not merely formal exponents. Direct congruence calculation gives

    (1,3): n=19 mod32,     endpoint (9n+5)/16;
    (3,1): n=29 mod32,     endpoint (9n+11)/16.

For example, 19 -> 29 -> 11 has exact valuations 1 then 3. At source 19 the word (3,1) is invalid already at its first step. At source 29, the legal word (3,1) gives 29 -> 11 -> 17. The two formal matrices have identical diagonal data but different authentic source cells.

Thus the J-read plus known diagonal data distinguishes these source guards; a carry response alone, without its source, is not a termination certificate. Where a first-hit ROOT convention is required, actual intermediate states must additionally be checked for earlier ROOT visits. The argument neither asserts a decreasing potential nor gives a uniform bound on orbit length.

## 4. Exact positive hidden-state lower bound

Let X be a nonempty finite source set. For each allowed context c let r_c:X→Y_c be a deterministic response, with Y_c finite. Let x~x' mean equality of every response. Let q be the number of equivalence classes.

An h-state *positive contextual realization* consists of a probability vector p_x in R^h for every x and, for every c, a column-stochastic matrix O_c with rows indexed by Y_c, satisfying

    O_c p_x = e_(r_c(x))  for every x,c.

There is no requirement that x→p_x be linear or injective. No dynamics is assumed in this definition; a stochastic dynamical implementation is a special case when its contextual output channels have this form.

**Theorem 2 (positive state cardinality).** Such a realization exists if and only if h≥q. In particular randomizing the encoding cannot reduce the exact number of positive hidden states below the number of contextual classes.

**Proof.** Fix a hidden state j with p_x(j)>0. Since O_c p_x is a vertex of a probability simplex and is a convex combination of the columns of O_c, every column with positive coefficient must be that vertex. Therefore O_c e_j=e_(r_c(x)) for every c. If p_x(j)>0 and p_x'(j)>0, all responses of x and x' must coincide. The supports of encodings from distinct classes are disjoint and nonempty; choosing one source from each class proves h≥q. Conversely, use one hidden state per class, encode x as its class point mass, and set each channel column to its class response. Extra unused states may be appended. ∎

An infinite context family causes no difficulty: the support argument applies separately to each context and X still has finitely many classes. The theorem requires nonnegative channels and probabilities, deterministic exact outputs, and an encoding for every pure source. Approximate outputs, signed models, or a restriction to mixtures of sources are different problems.

If X also has allowed maps f_i:X→X, then the q-state realization implements them deterministically exactly when ~ is preserved by every f_i. If contexts include every permitted future observation, this congruence property holds: every response after f_i is another declared context. Without that closure, the theorem's channel realization alone does not supply an update law. A dynamical realization still obeys the lower bound even if no q-state update exists.

## 5. Applications and resource accounting

For the two Collatz words at fixed diagonal data, triangular contexts yield q=1. Adding J yields q=2. No positive stochastic one-state model can reproduce both exact J responses. The matrices here are abstract source records, not an autonomous stochastic model of Collatz.

For Session171's storage family h in RM(d-r,d), a complete bit-address bank separates all 2^K sources, where

    K=sum_(j=0)^(d-r) binom(d,j).

Therefore every positive contextual realization reproducing all pure-source addressed responses exactly has at least 2^K hidden states. This lower bound is attained abstractly by class point masses. It does not assert that every such hidden state is a new physical tournament.

If the declared bank E consists only of individual binary linear responses and has rank a, there are 2^a contextual classes and its exact positive minimum is 2^a. Include richer subset-weight queries only after recomputing their actual fibers: they are not one binary linear measurement.

For d=4,r=2, K=11. The constant bulk profile has one class. The five addresses of weight at most one have rank five and give 32 classes; all addresses, or the eleven addresses of weight at most two, give 2048 classes. Those class counts require respectively 0,5,11 bits to name, not 1,32,2048 bits. The theorem concerns the number of states in a positive simplex model. It does not conflict with an eleven-bit representation or with a smaller signed-linear observer dimension.

No new Hamilton-count spectrum value is obtained here. S171's physical code/read-head formulas are inherited, while the lower bound describes their exact abstract positive realization cost.

## 6. Reproduction

Run `python3 check.py > observed.json` and `python3 -O check.py > optimized.json`, then compare both with `results.json`. Checks use explicit exceptions and exact integer or rational arithmetic. They cover triangular contexts and probe formulas, 5,460 nonempty words over exponents 1..4 through length six, literal guards for both example words, the exact finite response partitions of RM(2,4), and positive class-state channel realizations. Finite controls validate transcription; proofs above establish the general statements.
