# A strong two-port counter carrier whose path spectrum stays sparse

**PROVED elementary construction and observer obstruction; FINITE-EXACT
controls.** This executes the declared two-counter instruction language.
No external universality theorem is used here, and this is not a proof of
Hamilton-count coverage.

The closest mechanism is the exact triangle wrapper from
[Session159, positive profile operator](../../session159_constructions_20261103/profile_operator/NOTE.md).
The hostile is a computational state whose scalar output does not change
when its registers change. The corrected near miss is treating an integer
code as realization of that integer as H. The least-used sidecar is the
ordered SCC word exposed after deleting the two wrapper ports. The live
board is counter contents / named ports / full cover profiles / native
transport / exact output image.

## 1. Literal strong sources and two birth instructions

Let TT_a be the transitive tournament on a vertices, including TT_0 empty.
Let C be a directed triangle. For a,b>=0 form

    X(a,b) = TT_a -> C -> TT_b,
    W(a,b) = X(a,b) -> u -> v -> X(a,b).

The named ports u,v persist throughout. All W(a,b) are strongly connected:
the three nonempty blocks X,{u},{v} form a directed triangle. Removing u,v
recovers X; its unique nonsingleton SCC is C. Thus the ordered singleton
SCCs before and after C recover the two counters a,b from the actual graph.

The following are literal operations, with no relabelling of old vertices:

* Increment a: add z losing only to v. Within X, z is a new universal source.
* Increment b: add z beating only u. Within X, z is a new universal sink.
* If a>0, decrement a by deleting the first singleton SCC of X.
* If b>0, decrement b by deleting the last singleton SCC of X.
* Test a=0 or b=0 by inspecting the corresponding side of the C component.

Every endpoint is again the displayed strong source with the same u,v.
The two increments are exactly the two-polarity single-port births already
supported by the eighteen-polynomial carrier in
[GRADED_OBSERVER G1/G2](../../source_selection_closure_20261103/GRADED_OBSERVER.md).
Both have old-to-new mask {v}; their homogeneous polarity is opposite.
Decrements are guarded deletions of permanent counter vertices. They are
an enlargement of the previous birth-only language, not temporary resets
silently reinterpreted as counter instructions.

Consequently every deterministic program with instructions
INC(i,next) and DECJZ(i,nonzero,next_zero) has an exact simulator on these
actual tournaments, with its finite control state retained. A DECJZ step
takes its zero branch without deleting when the register is zero and
otherwise performs the specified deletion and nonzero branch. Induction
on execution length proves equality of register contents and tournament
SCC sizes at every step. This statement does not require a halting claim.

There is already a literal three-fraction test of the connection. Encode
the two data registers in primes 2,3 and use control primes 5,7,11. The
ordered FRACTRAN list

    [7/10, 11/5, 15/7],       input 5*2^a*3^b,        (0)

executes DECJZ(A,increment_B,halt) followed by INC(B,loop). For a>0 the
first fraction gives 7*2^(a-1)*3^b and the third gives
5*2^(a-1)*3^(b+1). For a=0 the second gives 11*3^b, where no fraction
applies. Thus it halts after exactly 2a+1 steps at 11*3^(a+b).
The tournament program has exactly the same instruction-boundary states
and clock. Reversing the first two fractions makes input 10 halt at 22
instead of 33: priority is an actual guard, not unordered fraction data.
This tiny arithmetic proof is self-contained; it makes no universality
claim about the three-fraction program.

## 2. Exact image and a full-observer collision

Let F_r(Y) count ordered r-slot spanning path covers, allowing empty slots.
The inherited ordinal law is F_r(A->B)=F_r(A)F_r(B), and F_r(TT_a)=r^a.
The triangle has

    F_r(C)=r^3+2r.

For example its unordered cover polynomial is 3z+3z^2+z^3; evaluation in
the falling-factorial basis gives this identity. Hence, with s=a+b,

    F_r(X(a,b)) = r^(s+3)+2r^(s+1),
    H(W(a,b)) = H(X)+F_2(X) = 3+12*2^s.               (1)

In particular every value of this family is at least 15 and divisible by
3. It hits no prime and no target 7p with p an odd prime other than 3.
For p=3 the target is the inherited forbidden 21 and is also absent from
(1). Thus a complete counter computation can live in a family contributing
none of the target atoms. This is an exact all-parameter obstruction, not
an empirical density claim.

**Stronger observer theorem.** The full eighteen endpoint-cover polynomials
of W(a,b), marked at (u,v), depend only on s=a+b.

To prove this, start with the unmarked cover polynomial of X. First append
u as a sink, retaining its endpoint name; then append v losing to u and
beating X, retaining both u,v. Deleting each newborn from a path cover
expresses its marked coefficients using only the unmarked old cover
coefficients. This is precisely the endpoint joining bijection behind
GRADED_OBSERVER G2, applied with zero, then one old marked port. Since the
entire unmarked polynomial of X is determined by all F_r(X), the displayed
identity proves equality of every marked coefficient for a fixed s.

The finite-difference inversion is legitimate: the polynomials (r)_k are a
basis, so equality for every nonnegative integer r gives equality of all
unordered cover coefficients. No source graph is inferred from those
coefficients.

For the minimal counter collision s=1, W(0,1) and W(1,0) therefore have
identical entire eighteen-polynomial vectors, although the a=0 test is true
only for the first. They are distinct as (u,v)-marked sources, since their
deleted ordered SCC words differ. All future two-polarity birth responses
at these ports agree by the exact forward carrier. The legal deletion and
zero-test guards do not factor through that carrier.

This failure starts when the operation alphabet adds tests and deletion;
it does not contradict the proved sufficiency for forward births. Retain
the ordered SCC split, or its two counter sizes together with actual named
source data. Extra polynomial degree alone does not restore this guard.

The loss is unbounded even at a fixed output: for fixed s, the s+1 states
W(a,s-a), 0<=a<=s, have identical marked cover vectors and the same
initial control state in (0), but their halting clocks are respectively
1,3,...,2s+1. At every completed nonzero loop the cover vector returns
to itself. Its apparent recurrence forgets how much executable work
remains. This is an actual source family, rather than an abstract product
countermodel.

## 3. A native reset moves information invisible to the observer

If a>0, let x be the first vertex of TT_a. Make it the last vertex of TT_b
while fixing every other old name and every old arc not incident to x.
This reverses exactly s+2 arcs: x's incidences with X minus itself. Its
incidences with u,v do not change. The endpoint is W(a-1,b+1), with the
same H and the same eighteen-polynomial vector by (1) and the theorem.

The inherited [one-temporary-vertex arc compiler, METAGRAPH section 5](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md)
implements each reversal with one temporary birth, two legal native moves,
and deletion. Thus this information-changing reset has an explicit upper
budget of 2(s+2) native moves, s+2 temporary births/deletions, and peak
order s+6. All original vertex names, including u,v, survive. The
intermediate graphs need not preserve H or strongness, and temporary
exception ports may move outside {u,v}; the full graph and legal path words
remain required. This is not a closed native operation on the eighteen
polynomials.

## 4. Integer encoding and the target that remains

The map n -> W(n,0) encodes every n>=0 injectively even in scalar H:

    n = log_2((H-3)/12).

Its construction uses n indegree-one births from W(0,0), hence n+5 vertices.
It is unary in n, not a compact binary code. Reading 2n or 2n+1 then
encodes every nonnegative even or positive odd integer as a decoded
observable. It does not make any even number, 0,7,21, or arbitrary odd
number into a Hamilton count. More economical binary and square-sum codes
are separate constructions; their existence alone has the same typing
boundary.

The useful connection is now exact: cyclic gluing supplies strong sources;
opposite birth polarities supply counter increments; native moves transport
counter allocation within one entire forward-observer fibre; and zero
tests expose missing source information. A prime/7p induction still needs
an output selector that leaves this sparse image while maintaining the
guards of its next operation.

## 5. Reproduction

Run `python3 -I -B counter_check.py` and the same command with `-O`.
The script uses independent Hamilton subset DP and endpoint-labelled
component enumeration. It checks all 16 pairs 0<=a,b<=3 for (1), strongness,
and exact increment/decrement transitions; all pairs with a+b<=3 for the
complete marked observer collision; a transfer-counter program on all 16
inputs; and the s+2 old-arc difference for those applicable inputs.
It also checks the coverage-address map from AXIOM.md and small inherited
core-budget values. These are finite controls of the proofs, not a spectrum
scan or an all-order machine-universality computation.

The gluing lane independently audited the full-polynomial collision,
strongness, increment masks, zero-test separation, and inherited native
compiler budget, with no corrections requested.
