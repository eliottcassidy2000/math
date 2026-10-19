# Independent review: finite absolute-difference strips

2026-10-19. Status: **PROVED proof review**, with separately bounded exact checks.
Reviewed root draft: `BINARY_ENTRY_AND_EXACT_MIXING.md`.
No checkout or Git modifications were made. This review makes no external novelty claim.

## Verdict and corrections

The maximum-rigidity proof, constant-one entry theorem, rational hitting-probability argument, exact fair F-step binary mixing theorem, no-synchronization distinction, and full-chain total-variation bound are correct as stated in the reviewed draft. Two stronger statements below improve the deterministic and stochastic conclusions. Sharpness of the constant-one word must not be described as optimality over all binary forcing words: the parent reports shorter nonconstant universal-entry words in some boxes.

Primary precedent checked directly: Chase–Hunter–Tao, *Gilbreath's conjecture: a Cramér random model and a deterministic analysis*, arXiv:2607.08712v1, Lemma 3.8 (Parentage), printed page 14, https://arxiv.org/pdf/2607.08712. Under a maximum cap equal to d, the shifted pair and the height-at-least-2d alternatives are excluded. Iterating the remaining parentage alternative proves integer maximum rigidity. The real-valued proof below is elementary. The paper cites Chase's 2024 *A random analogue of Gilbreath's conjecture* as an earlier source of related block lemmas; no priority beyond this verified precedent is asserted.

## 1. Rigidity, with all directions explicit

Let A_n(a_0,...,a_n) be the apex after n absolute-difference rows. For M>0 and all a_i in [0,M],

    A_n(a)=M implies a_i belongs to {0,M} for every i.

Induct on n. The two child apices must be 0 and M. Reflect if necessary so that the first n seeds produce M. The induction hypothesis makes those seeds endpoints, and every cell of their triangle an endpoint. The remaining child apex, as a function of the final seed z, is a composition of maps z -> |z-c| with c in {0,M}. On [0,M] each map is z or M-z, so the composition is z or M-z. Since that child apex is zero, z is an endpoint. This also covers n=1; n=0 is immediate.

Endpoint seeds alone do not imply apex M: their normalized binary Pascal parity must be one. For example (M,M) has apex zero. Integrality is unnecessary here, but is essential for the later decrement-by-one and binary-entry conclusions.

## 2. Exact characterization of universal F-step maximum-drop words

Let F>=1 and m>=2 be integers. Evolve x=(x_0,...,x_(F-1)) in {0,...,m}^F using binary boundary inputs w_0,...,w_(F-1):

    x_i(t+1)=|x_i(t)-x_(i+1)(t)| for i<F-1,
    x_(F-1)(t+1)=|x_(F-1)(t)-w_t|.

The following are equivalent:

1. Every initial state in {0,...,m}^F has maximum at most m-1 after these F steps.
2. The first forcing bit w_0 is one.

Proof of sufficiency. Every binary word w of length F is the first F time values of the leftmost cell of some binary spatial tail b_0,...,b_(F-1). Indeed,

    w_t = sum_(j=0)^t binom(t,j)b_j modulo 2,

and this is a unit-triangular linear system. Thus b_0=w_0=1. Append that binary tail to the initial state. Its unforced absolute-difference evolution realizes the prescribed boundary word because absolute difference equals XOR on bits. Every width-F strip cell at time F has a seed cone containing the spatial location F, whose value is b_0=1. All seeds lie in [0,m], and 1 is interior to that interval. Maximum rigidity excludes output m. Integer-valuedness gives the bound m-1.

Proof of necessity. If w_0=0, start from (m,0,...,0). At time F the leftmost cone has initial row (m,0,...,0,w_0), so its apex remains m. Later forcing bits cannot enter this cone. Thus the required universal decrease fails.

At all times the capped height H(x)=max(1,max_i x_i) cannot increase under binary forcing. Every binary forcing word of length F whose first bit is one decreases H by at least one unless H=1 already. Consequently the constant-one word of length (M-1)F maps the entire box {0,...,M}^F into the binary invariant region.

## 3. Sharpness of the constant-one bound at dyadic widths

For every F that is a power of two and every integer M>=2, the starting state (M,0,...,0) under constant-one forcing first enters the binary region at exactly

    tau=(M-1)F.

Proof. All coordinates except the leftmost are binary forever. Their evolution is the impulse response of the binary tail (1,0,0,...) appended at location F. The right neighbor of the leftmost cell at time t is

    b(t)=binom(t,F-1) modulo 2.

For F=1 this denotes the constant boundary itself. While the leftmost coordinate is greater than one, its next value is its current value minus b(t). When F=2^e, Frobenius gives

    (1+z)^t = product_(j>=0) (1+z^(2^j))^(t_j) modulo 2.

Thus binom(t,F-1) is odd exactly when all e low bits of t are one, or t is congruent to F-1 modulo F. The leftmost value therefore drops at update times F,2F,..., and reaches one exactly at (M-1)F. All other coordinates were already binary.

For arbitrary widths the same sparse initial state satisfies, up to entry,

    x_0(t)=M - #{0<=s<t : (s & (F-1))=F-1}.

This is an explicit formula for this starting state, not a worst-case formula for every width. Even at dyadic widths the constant-one word need not be the shortest universal-entry word among all binary forcing words.

## 4. Improved iid input bounds

Let boundary bits be independent Bernoulli(q), with 0<q<=1. Let tau be first entry to the binary region, and let r=max(0,M-1). The case r=0 is tau=0.

### Disjoint-block bound

For r>=1, partition time into blocks of F updates. A block succeeds whenever its first bit is one; these events are independent Bernoulli(q). Each success reduces capped height by one, and unsuccessful blocks cannot increase it. Hence, for every integer k>=0,

    P(tau>kF) <= P(Bin(k,q)<=r-1),
    E[tau] <= Fr/q.

This replaces the much weaker all-one-block probability q^F by q.

### Renewal bound

Start a phase at time s. Wait until the first one, at update time s+G-1, where G has the positive geometric distribution P(G=j)=(1-q)^(j-1)q. Then run the F updates beginning with that one. The phase ends at time s+G+F-1 and decreases capped height by at least one. Restart there. Strong independence of future iid inputs makes successive G values independent. After r phases the state is binary, so

    tau <=st r(F-1)+sum_(j=1)^r G_j,
    E[tau] <= r(F-1+1/q).

More explicitly, for every integer n>=0,

    P(tau>r(F-1)+n) <= P(Bin(n,q)<=r-1).

The same construction proves pathwise domination by the phase-end time on the original forcing sequence, not merely an expectation estimate. For q=1 it recovers the deterministic bound. For fair input the expectation is at most (M-1)(F+1).

## 5. Mixing, Markov rationality, and limitations

Inside the binary region, the update is x -> A x+e_F u over F_2, with A=I+S and S the nilpotent right-neighbor shift. Since A is invertible, a fixed forcing word is always a permutation of binary states, so no forcing word synchronizes all of them to one state. The vectors e_F,Ae_F,...,A^(F-1)e_F form a basis by the unit-triangular Pascal expansion. Consequently F independent fair inputs send any fixed binary state to the uniform law exactly. F is necessary for a fixed starting point because fewer than F fair bits produce fewer than 2^F possible states.

Conditioning on the stopping time tau validates the parent's full-chain bound

    TV(Law(x(t)),Uniform(binary states)) <= P(tau>t-F), t>=F.

Every event {tau=s} depends only on inputs before s; the subsequent F fair inputs are fresh. Uniformity is preserved thereafter. The finite chain has no nonbinary recurrent class and has a unique stationary law, uniform on binary states.

For rational q>0, designate any failure set before binary entry and stop at the first failure or binary entry. The transient submatrix Q has rational entries and spectral radius below one by the explicit absorption bound. Therefore (I-Q)^(-1) exists over Q, and all hitting probabilities and expected stopping times from rational initial distributions are rational. No irrationally weighted initial mixture should be silently included in the rationality claim.

The following boundaries are essential:

- At q=0, (M,0,...,0) is a stationary nonbinary state.
- For noninteger seeds, binary entry can fail: F=1 and x=2.3 under input one eventually alternate between 0.3 and 0.7.
- A fair iid spatial binary sea gives a fair iid boundary trace by the invertible triangular Pascal transform. A biased spatial sea generally does not give an iid biased temporal trace. The Bernoulli(q) theorem is a theorem about independent temporal forcing.
- A terminal binary row does not erase failures that occurred earlier. Failure events must be recorded before entry.
- Neither bounded finite strips nor independent boundary hypotheses describe the entire deterministic prime-gap sequence. Gilbreath's conjecture remains outside this theorem.

## 6. Reproduction

Run from this directory:

    python3 independent_reset_review.py > independent_reset_review.json
    python3 -O independent_reset_review.py > independent_reset_review_optimized.json

The independent checker includes rational-grid maximum rigidity, triangle/strip embeddings, every binary forcing word through width five, all small box states, exact rational stochastic tail checks, sparse-defect formulas through width 65 and magnitude 25, fixed-word binary bijectivity and exact fair mixing through width seven, and hostile zero-input/noninteger controls. Exact universes and counts appear in the source and JSON. The all-M dyadic claim is justified by its proof, not by this finite verification.
