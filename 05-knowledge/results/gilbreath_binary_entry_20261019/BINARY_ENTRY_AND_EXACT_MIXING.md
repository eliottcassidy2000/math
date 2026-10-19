# Finite Gilbreath strips: linear-time binary entry and exact mixing

Status: PROVED elementary consequences; independently reviewed.
No claim about the prime initial sequence or historical priority.

## Inheritance and question

The repository's THM-4511 (`01-canon/theorems/THM-4511-gilbreath-size-four-wall-theorem.md`)
settles size-four walls. The exact Markov reduction in
`05-knowledge/results/gilbreath_size6_extinction_20260926.md`, Theorem A,
encodes a finite defect region driven by the trace of a fair binary sea.
Its rational absorption calculations motivate a uniform absorption certificate
for every finite bounded region, rather than only individual small chains.

Chase–Hunter–Tao, *Gilbreath's conjecture: a Cramér random model and a
deterministic analysis*, arXiv:2607.08712v1, Lemma 3.8 (Parentage), already
contains the maximum-rigidity mechanism used below. We give a short independent
proof for readability and credit that source. The new use here is the explicit
finite-strip entry bound, its stochastic certificate, and the exact binary
mixing/no-synchronization distinction. External novelty has not been established.
Primary source: https://arxiv.org/pdf/2607.08712

## 1. Maximum rigidity (cited mechanism, elementary proof)

Let a length n+1 row have entries in [0,M], M>0. If its n-fold
absolute-difference apex equals M, every seed entry belongs to {0,M}.

Proof by induction. For n=0 this is immediate. For n>=1, the two child
apices belong to [0,M], so their difference has magnitude M only when they
are 0 and M. Reflect the row if needed so the left child is M. Induction
makes the first n seed entries endpoints 0 or M; every cell of their
triangle is therefore an endpoint. With the last seed denoted t, the other
child apex is a composition of maps t -> |t-c| with c in {0,M}. On [0,M]
these maps are respectively t and M-t. Its value is 0 only when t is an
endpoint. This completes induction.

This is also the bounded-maximum case of CHT Lemma 3.8: its shifted-endpoint
and height-at-least-2M alternatives are unavailable.

## 2. A universal binary-entry word and a single-hit maximum drop

Fix integers F>=1, M>=1. Let x(t) in {0,...,M}^F evolve by

    x_i(t+1) = |x_i(t)-x_(i+1)(t)|,  1<=i<F,
    x_F(t+1) = |x_F(t)-u(t)|.

For any binary input word of length F whose first bit is 1, the maximum
falls from at most m to at most m-1, for every integer m>=2. No condition
on the other F-1 input bits is needed.

Proof. The binary right tail can realize any prescribed boundary trace:
its time-t first cell is sum_j binom(t,j)b_j modulo two. This triangular
Pascal transform is invertible, and its first seed is b_0=u(0)=1. Extend
the resulting length-F seed arbitrarily by binary cells. Every one of the
F cell cones after F steps contains the interior seed b_0=1. All seeds
are in [0,m], so maximum rigidity prevents an output m. Integrality gives
the claimed drop. Binary forcing never increases a maximum already >=1.

This characterizes the universal F-step maximum-drop words exactly: a
word beginning 0 fails on the initial state (m,0,...,0), because the first
cell's F-step cone consists of m followed only by zeros, including b_0=0.
Its apex is still m; later input bits lie outside this cone.

In particular, constant input 1 throughout L=(M-1)F steps sends every
state into {0,1}^F. More generally it suffices to have M-1 input ones
separated by at least F time steps, and wait F steps after the last.
For M=1 the state is already binary. The binary region is invariant.

For the constant-one word, this bound is sharp at every dyadic width
F=2^a and every M>=2: start at (M,0,...,0). The right block is binary and
the neighbor driving the first coordinate is binom(t,F-1) modulo two.
Lucas's theorem makes it 1 exactly when t is F-1 modulo F. The first
coordinate therefore falls by one exactly every F steps until it reaches
1, at time (M-1)F. This concerns the constant-one word, not the shortest
possible forcing word; other words can be faster uniformly.

## 3. Negative-binomial transient tail and rational hitting probabilities

Let the inputs be independent Bernoulli(q), 0<q<=1, and let tau be first
entry to {0,1}^F. Put r=M-1. Starting from each current time, wait for
the next input 1, then allow F steps for the maximum drop; repeat r times.
The waiting times G_1,...,G_r, counted including the successful input,
are independent geometric(q) variables on {1,2,...}. Consequently

    tau <= r(F-1)+G_1+...+G_r  in stochastic order,
    P(tau>r(F-1)+n) <= P(Bin(n,q)<r),  n>=0,
    E[tau] <= r(F-1+1/q).

The bound is pathwise for the greedy selected input-one times; independence
only enters the displayed distributional evaluation. For a fair sea the
mean bound is (M-1)(F+1), polynomial in both width and height. For q=1
the dyadic-width examples attain the deterministic bound rF. For M=1,
tau=0 and the r=0 conventions are harmless.

For rational q and a fixed starting state (or rational initial mixture), probabilities of ever hitting any designated failure set
before binary entry are rational: the finite-state absorbing equations have
rational coefficients, and the above uniform absorption bound makes the
transient matrix I-Q invertible. This removes the possible closed-transient-
class loophole in an unsupported invocation of a finite Markov chain.
For q=1 the dynamics are deterministic and the same statements apply.

For the original unscaled Gilbreath region with even entries <=2M and a
right tail of independent fair {0,2} cells, divide by two. Its boundary
trace is independent fair Bernoulli: the time-t trace is the Pascal transform
of the first t+1 initial tail bits, with the newest bit occurring with
coefficient one. This triangular transform preserves the uniform product
measure. Thus q=1/2 applies exactly to that random-sea model. A biased
initial sea does not in general give independent biased boundary inputs;
we make no such identification.

## 4. Exact finite-time mixing inside the sea

Over F_2, write A=I+S, where S shifts each cell's right neighbor into that
cell and sends the last coordinate to zero. The binary evolution is

    x(t+1)=A x(t)+e_F u(t).

A is invertible. Therefore every fixed input word is a bijection of the
2^F binary states. In particular there is no synchronizing word sending
all states to a single state. Section 2 is entry into an invariant set,
not single-state synchronization.

Nevertheless, with fair independent inputs, after exactly F steps the
state is uniformly distributed on all binary states, independently of
its starting binary state. Indeed the control matrix with columns

    e_F, A e_F, ..., A^(F-1)e_F

is invertible: A^j e_F=sum_(r=0)^j binom(j,r) S^r e_F, a triangular
Pascal change of basis with diagonal ones. Thus each desired final state
has exactly one length-F forcing word. F is minimal for a point starting
state, since fewer than F random bits cannot support all 2^F outputs.

Consequently the full bounded fair-input chain has unique stationary law
uniform on the binary sea, all other states transient, and

    ||Law(x(t))-Uniform({0,1}^F)||_TV <= P(tau>t-F), t>=F.

Condition on each stopping-time value tau=s<=t-F; future independent
inputs give exact uniformity after F more steps. The remaining mass is
at most the displayed tail. Uniformity remains invariant under later steps.

## 5. Connections and exact limits

Source -> target: finite even defect region plus a fair binary right sea
-> finite Markov strip with independent boundary forcing.
Preserved: every finite space-time value and first failure event.
Lost if using only binary reduction: magnitude and earlier failures.
Restoration: retain the bounded integer state until entry and separately
record every pre-entry edge value. A safe terminal row alone does not
certify earlier Gilbreath rows.

The control matrix is Pascal; its invertibility explains both retention
of all binary-state distinctions for a fixed forcing word and loss of all
initial-state information after averaging F fair inputs. This is a precise
connection to the session's Hasse-deletion/divided-power lane, not an
isomorphism with Collatz or a transfer of prime-input randomness.

The real prime gaps are a deterministic sequence with infinitely many
unbounded entries. No finite F,M or independent boundary law covers that
entire sequence. These results do not prove Gilbreath, a prime-gap bound,
Collatz, or a Rule30 prize.

## Reproduction

    python3 verify_binary_entry.py > verification.json
    python3 -O verify_binary_entry.py > verification_optimized.json

Checks cover maximum rigidity, every state in explicitly enumerated small
boxes, exact F-step output multiplicities, and the hostile no-synchronization
control. Bounds and universes are encoded explicitly in the source.

## 6. Exact control problem and a new open question

Define ell(M,F) as the shortest binary input word sending the entire box
{0,...,M}^F into the binary sea. Complete breadth-first enumeration of the
subset automaton gives, for widths F=1,2,3,4:

- M=2: 1,2,3,4.
- M=3: 2,4,6,8.
- M=4: 3,6,7,9.
- M=5: 4,8,9,12.

These are FINITE-EXACT minima in the stated boxes, not asymptotic claims.
For example ell(4,4)=9, witnessed by 010111000, improving on the constant-
one length12. The BFS explores image sets in increasing word length and
retains each image set once, an exact equivalence because future evolution
depends only on that set; it uses no unproved pruning. A separate replay
checks every original state for each returned word.

OPEN HYPOTHESIS: ell(3,F)=2F for every F>=1. The upper bound is proved
above; only a lower bound is missing beyond the enumerated range. The
M=4 data already refute the blanket extrapolation ell(M,F)=(M-1)F.
A future lower-bound route should retain the full image set, since fixed-
word invertibility on the binary sea prevents single-state reset arguments.

Reproduce with `python3 shortest_entry_words.py`; its optimized replay
uses the same explicit checks and produces identical JSON.
