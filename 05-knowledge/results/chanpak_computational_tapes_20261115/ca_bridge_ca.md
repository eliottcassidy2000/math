# A left-permutive tournament trace, with exact modular memory

**Scope.** This is a genuinely physical trace along a fixed tournament
operation. It is an all-order observation/decoding result, not a claim
that arbitrary traces, decoded profiles, or terminal odd counts occur.

## What transfers from the cellular-automaton work

The archive's [THM-3456](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary.md)
recovers a free boundary from a temporal trace by permutivity; it does
not decide whether the boundary is the physically distinguished
single seed. [THM-4204](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md)
further distinguishes spatial inverse automata from the temporal
physical trace: an unpointed quotient eventually ignores the center.
[THM-4210](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md)
keeps a transverse channel even when the current one-scale center
does not read it, because later nonlinear generation does. Here an
exact tournament operation supplies the moving reader, and factorial
weights identify precisely which transverse grades can return at each
arithmetic resolution.

## Physical update

Let Q_j(D) count covers by j ordered optional directed-path slots for
nonempty D, so Q_0(D)=0 and Q_1(D)=H(D). Let W(D) be the tournament
obtained by placing D, a singleton a, and a singleton b cyclically:
D→a→b→D. For every j≥1,

    Q_j(WD) = j (Q_(j+1)(D) + j Q_j(D) − Q_(j−1)(D)).

To see this directly, choose the slot of b (j choices). If a is in a
different slot (j−1 choices), it is at that slot's end and b is at
its own slot's start, so the D assignment is arbitrary in j slots.
If a and b are in the same slot, the possibilities with both in
order a,b split that slot into two independent D pieces, adding
Q_(j+1). In the remaining possibility b precedes a and requires a
nonempty intervening D path, adding Q_j−Q_(j−1).
These possibilities are disjoint and exhaustive. This update is also
the singleton specialization of the prior higher-rank Gram formula
(see the Bott/Cayley session); the construction below applies it to an
entire trace rather than just one reading.

Make the integer change of observable e_j=(j−1)!Q_j for j≥1, e_0=0.
The effective dynamics is now the radius-one, right-permutive,
nonuniform linear cellular automaton

    e'_j = e_(j+1) + j² e_j − j(j−1)e_(j−1).

The observed temporal trace is h_t=e_1(W^t D). A depth-t cone from
the observation reaches at most position t+1. The only walk to that
position uses the unit right-neighbor coefficient every time. Thus

    h_t = e_(t+1) + ∑_(j=1)^t a_(t,j)e_j      (a_(t,j) ∈ Z).

Equivalently the finite trace map is unitriangular over Z, with no
division needed to recover the factorial-weighted state. In the
unweighted Q_j it has diagonal (0!,1!,2!,…).

## All-order consequences

1. **Exact history replacement.** The first r readings uniquely
   determine Q_1,…,Q_r over Z. All wrapped H readings agree for D,E
   if and only if their full slot polynomials agree. If there is a
   first different grade j, they first diverge at time j−1 by exactly
   (j−1)! [Q_j(D)−Q_j(E)]. In particular the two physical Cayley
   modules with H=495, Q_2=43200 and Q_3 difference 4320
   have identical readings 495 and 43695, then separate by 8640.
   That separation is exact on actual strong tournaments after the first lift because W
   always yields a strong tournament for nonempty input. It is
   invisible at first order, but it is not an invariant of the tempo.

2. **Exact modulus budget.** Two modules have the same first r
   readings modulo M if and only if, for each j≤r,

       Q_j(D) ≡ Q_j(E) mod M/gcd(M,(j−1)!).

   In particular if K is the least index such that M divides K!,
   the first K readings determine *every* future reading modulo M.
   The finite dynamical system stores e_1,…,e_K mod M and has zero
   boundary e_(K+1)=0. Indeed all higher e_j are multiples of M for
   every physical tournament, and the weighted recurrence respects
   this fixed zero boundary. For M=p^a the lost precision grade by
   grade is exactly the factorial tariff ∑floor((j−1)/p^i), capped
   at a, in the stated equivalence. This is a state reduction, not
   the claim that all formal states occur.

3. **Prime refinement with known order parity.** Modulo odd p,
   K=p. Since Q_p≡Q_0=0 and
   Q_(p−j)≡(−1)^|D|Q_j by parity of the slot polynomial,
   only d=(p−1)/2 coordinates are needed. W preserves order parity.
   The first d readings recover those d residues since their
   factorial diagonals are units modulo p. Thus for sources of the
   same order parity the first d readings determine every future
   reading modulo p, even if orders differ. As elementary examples,
   modulo 3 the update is h↦2h for even-order sources and h↦0 for
   odd-order sources. Modulo 5, using coordinates (Q_1,Q_2), the two
   matrices are [[1,1],[-2,6]] in even order and [[1,1],[-2,2]] in
   odd order. These are exact phase laws, not independence claims.

## Inverse-use proposal and boundary

A prescribed target together with several proposed continuation
counts now specifies exact linear congruences and divisibility
constraints on the *whole necessary input*, not just its terminal
value. This gives a reverse compiler: specify a desired small
adjustment with unchanged past, then recover the least grade that
must move and its exact signed magnitude. For example a return with
two frozen initial readings and terminal shift 2 requires a
one-unit Q_3 change, which is impossible between same-order profiles
(the coefficients of the difference of slot polynomials are even).
More generally all nonzero differences under r frozen initial
readings are multiples of r! in the next reading, before any
additional same-order divisibility. So hiding a change at very high
grade until a final wrap necessarily makes the first reveal coarse;
the useful target for a genuine small repair is a *changed* reader
that brings low weighted grades into play earlier while restoring
the needed continuation data by a second route. This is the literal
analogue of changing the observer rather than claiming a free CA
boundary is an actual seed.

The map on unrestricted integer profiles is freely invertible one
coordinate at a time, but actual input profiles are slot polynomials
of actual tournaments. They must obey parity, positivity in the
odd-cycle-packing representation, and all other physical constraints.
The exact decoder moves those constraints to the correct side of the
problem; it does not provide their sufficiency, construct a changing
reader that avoids the tariff, or prove Chan–Pak.

The included checker uses literal subset convolution of ordinary
Hamiltonian path counts to obtain optional slot-cover counts. It
compares 225 physical source/wrap profiles and 3,600 modular
observations including nonprime and prime-power moduli, with identical
normal and optimized Python results. The proof does not depend on
these tests.

## Adaptive side reservoirs within the one-singleton architecture

Replacing a by an arbitrary nonempty physical reservoir A gives
`W_A(D)=C3[D,A,1]`. The same slot argument, now including the two
optional pieces of A, gives the prior Gram formula in nearest-neighbor
form. Put u_j=Q_j(A), d_j=u_(j+1)−u_j. In factorial state it is

    e'_j = d_j e_(j+1)
           + j (u_j−d_j+d_(j−1)) e_j
           − j(j−1)d_(j−1)e_(j−1).

Proof of the same-slot rule: with b in a designated slot, before b
one visits D then A and after b one visits D then A. A nonempty D
piece before b requires nonempty A before b. A nonempty A piece
after b requires nonempty D after b. Inclusion-exclusion on these
two violations gives
`d_j Q_(j+1) + (u_j−d_j+d_(j−1)) Q_j − d_(j−1)Q_(j−1)`,
multiplied by j choices of the designated slot. All other slots
just visit D then A with either piece optional. For A a singleton
all d_j=1, recovering the unit decoder.

This yields an exact *delayed-write tariff* even when each step uses
a different reservoir, provided both runs use the same sequence.
If starting profiles first differ at grade r, then their readings
agree through time r−2 and at time r−1 differ by

    (r−1)! ΔQ_r · ∏_(s=1)^(r−1) d_(r−s)(A_s),

where A_1 is the first reservoir applied. The unique earliest
causal path explains the product. The increments d_j are positive:
add a previously unused slot and, for instance, put a Hamiltonian
path of A in it. Thus there is an immediate difference exactly
at its predicted time, and the lower factorial bound does not
go away merely by changing a shared reservoir schedule. If the
schedule is chosen deterministically from the shared observed
prefix, both runs necessarily use the same choices up to and
including the first distinguishing step. To bypass the tariff
one must use more than that scalar prefix to select the reader,
retain a genuinely different source record, or expand beyond the
one-singleton architecture. This is a restriction on delayed
discrimination, not on producing a small step between arbitrary
unrelated profiles or differing schedules.

### Closed-form stronger tariff and the unique delayed +2

For equal-order modules with a first differing grade r, the
difference of the full slot polynomials has even coefficients,
the same parity of powers, and zeros at 0, ±1,…,±(r−1). Therefore
it is divisible in Z[t] by
`2 t^epsilon ∏_(i=1)^(r−1)(t²−i²)`,
with epsilon=1 (odd order), 2 (even order). Evaluate at r:
the first distinguishing H shift is divisible by the stronger
number

    2 r^(epsilon−1) (2r−1)! (r−1)!
      · ∏_(s=1)^(r−1) d_(r−s)(A_s).

Even without fixed order, a first distinction of magnitude 2 after
at least one equal reading is exceptionally rigid. For r≥4 the
factorial alone rules it out. For r=3, Q₃ is odd (all nonleading
coefficients of the slot polynomial are even), so 2!ΔQ₃ is
divisible by 4. For r=2, Q₂ is a multiple of 4 for every order
at least 2: its leading term contributes 2^n, nonlinear lower
terms contribute multiples of 4, and a linear term has even
coefficient. Thus one source must be the singleton; equal H=1
makes the other a transitive tournament, since any directed
triangle adds an odd-cycle-packing term. Its Q₂=2^n, so the only
possible source pair is TT1, TT2. The reservoir A must itself be a
singleton, as d₁≥2^|A|−1 by the leading monomial and the
nonnegative odd-packing coefficients. Conversely those two
inputs do realize 3 and 5 on the shared wrap, so this boundary
is sharp. This classifies delayed exact +2 in a broad physical
branch architecture; unrelated targets and independently changing
readers remain outside the obstruction.

### A genuine controllable mod-3 register (and its absorbing state)

Write h=H(D), a=H(A), and let s_D=(-1)^|D|, s_A=(-1)^|A|.
The j=1 formula, using Q₂≡sH mod 3, simplifies to

    H(W_A D) ≡ a h (s_A s_D − s_A − s_D) mod 3.

The parenthesized factor vanishes modulo 3 exactly when both
orders are odd, and otherwise equals -1 modulo 3. Thus 0 is
an absorbing register state for every allowed reservoir. From a
nonzero register state the physical controls below can write any
next residue, until the first zero; no reservoir within this
architecture restarts a zero register.

* Reset: choose a directed triangle, so a=3.
* Negate: choose TT2, so a=1 and the reservoir order is even.
* Keep: choose an order-four tournament with five Hamilton paths,
  so a≡−1 and the reservoir order is even. For example use arcs
  1→0, 2→0, 0→3, 3→1, 2→1, 3→2.

Both nonzero controls toggle source order parity and work on
both parities. This is more than a congruence enumeration:
one can program **any nonzero ternary residue trace** of a given
length on actual strong extensions, or reset once. But it does
not program exact H, and cannot generate a later nonzero output
after a modular reset. The rule says which desired computations
are impossible in this particular physical architecture.

The old [October 19 automata/arithmetic synthesis](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/automata_arithmetic_geometry_session_20261019/SESSION_MAP.md)
gives two complementary reasons to keep these distinctions:
Gilbreath binary boundary forcing mixes exactly after entry
without being synchronization to a fixed state; and the tournament
Hasse deletion theory identifies factorials as dyadic/prime-adic
costs of ordinary iteration, not carriers that can be erased by
claiming a divided-power inverse modulo p. The wrap trace
presentation gives an explicit new instance: its unscaled
integer observation matrix has Smith factors 0!,1!,…, while its
factorially weighted observation matrix is unimodular.

### Modulo 5: the hidden slot really restores a zero reading

With coordinates (Q₁,Q₂) and reservoir TT1 or TT3 there are
two parity-preserving physical updates over F₅. For even
source order TT1 acts by

    U_even = [[1,1],[3,1]].

Its characteristic polynomial is (x²-2x+3), irreducible
modulo 5. An eigenvalue λ in F₂₅ has norm λ⁶=3, of order four.
Thus its order is 8 or 24; λ⁴=λ+2 rules out 8 (which would
give λ⁴=−1). Therefore TT1 alone cycles through all 24 nonzero
two-slot states from any such initial state.

For odd source order, TT1 and TT3 act as

    U_odd = [[1,1],[3,2]],   V_odd = [[0,2],[1,4]].

The first has irreducible characteristic polynomial
(x²-3x+4); its eigenvalue has norm 4 and order 12:
order 4 would make its trace zero. Hence its nonzero states
split into two 12-cycles. The composite (U_odd V_odd) has
characteristic polynomial (x²+2), irreducible modulo 5,
and eigenvalue order eight since λ²=3. Neither 12-element
orbit can be invariant under an order-eight update all of
whose nonzero cycles have size eight. Thus each orbit can
exit to the other, and the two **positive** wrapping moves
reach all 24 nonzero states from any nonzero starting state.

In particular one can reach H≡0 but Q₂≠0 and then expose
Q₂ as the next H residue using TT1. This exactly contrasts with
the absorbing mod-3 scalar residue. If both coordinates are
zero, no reservoir in this architecture can restore them;
the update is homogeneous and the prime reflection rule gives
Q₃=0 also. So the useful abstraction is not whether zero is
ever observed but whether its *whole state* is zero. This is
a physical modulus-5 reachability result, not exact placement.

### Exact independence of the mod-3 and mod-5 registers on one host

Keep D even-order, start at TT2, and apply only odd-order
reservoirs, preserving the source parity. With TT1, the mod-5
state rotates by U=U_even while the mod-3 value negates. There
is a second control with the *same* 24-cycle but a different
mod-3 effect. Take the order-five source code 4 (lexicographic
upper-pair bits, 1 is forward); it has H=5, Q₂=64. The
mod-5 state rotates by

    V = [[1,4],[2,1]] = 3 U^(-1) = U^5,

because U^6=3I. The mod-3 value stays unchanged. Therefore one
V followed by nineteen U moves is a genuine return on the full
mod-5 state and flips the mod-3 value. It permits either mod-3
nonzero phase at every nonzero two-slot mod-5 state, without
changing the even order parity. Neither move purports to
return the absolute count.

Each register can also be reset independently:

* A triangle reservoir has H=3,Q₂=12. It resets mod 3 but
  acts by [[2,4],[2,2]] modulo 5, an invertible matrix.
* Take the order-eleven ordinal composition of code 4 at
  order five, code 12 at order five, and a final singleton.
  Code 12 has H=11,Q₂=100. Thus the composite reservoir
  has H=55,Q₂=12800, both divisible by 5, resetting the
  full mod-5 state while preserving a nonzero mod-3 state.

Using these moves and the nonzero controls reaches the full
75-element product F₃ × F₅² with actual tournaments. The
physical seed is fixed and all positive words yield strong
tournaments; no unconstrained abstract seed is substituted.
This is a working model of a target choosing which register
to alter without corrupting the other one. Its limitation is
sharp for the intended application: the return preserves the
other register only at the specified modulus, and executes
many substantial positive lifts. The exact-delay obstruction
does not contradict this modular independence and remains the
reason this machine does not solve exact target placement.
