# Atomic refinement, integer countdowns, and the finite-phase boundary

**PROVED elementary extensions and scoped obstruction. OPEN:** universal
Collatz convergence, escape from every infinite growing-return history, and
payment of growth before escape. Same-agent proof review and alternate exact
checks; no independent-agent or proof-assistant audit is claimed.

## 1. Inherit the actual source measure

For a positive odd source `n=2m-1`, retain the October 5 Elias weight

    w(m)=1/(2*4^floor(log2 m)).

The weights sum to one and every integer has positive mass. The inherited
[atomic-cylinder theorem](../algorithmic_collatz_measure_20261005.md), §4,
gives, for a source-index residue `m=r mod 2^L`, `0<=r<2^L`,

    mass(r,L)=4^-L + (w(r) if r>0, else 0).              (1)

This is a measure on labelled integers pushed through their actual words,
not a fresh Haar law at each return. A cylinder's finite representative may
be zero, but that does not identify its positive integer members with zero.

### Exact refinement converts tail mass into a named atom

The two children have residues `r` and `r+2^L` at level `L+1`. The old atom,
if any, stays in the first child. The new atom in the second child has weight

    w(r+2^L)=1/(2*4^L).

Each child has tail `4^-(L+1)`. Consequently

    old tail = two new tails + newly exposed atom.       (2)

The newly exposed atom accounts for half the old tail. This proves exact
mass conservation including the `r=0` boundary. Which child receives which
next dynamical bit depends on the native affine carry; the displayed
children are source congruences, not assumed dynamical zero/one labels.
Dropping the new atom would manufacture a factor-two contraction.

## 2. Apply this to the exact five-return grammar

Use [INTEGER_RETURNS](../ternary_integer_returns_20261029/INTEGER_RETURNS.md)
and the stronger incoming seven-prefix proof cited there. For
`a=v2(3n+1)`, the queue starts at one and updates by `-1+3[3|a]` until zero.
Its termination is proved on every positive odd integer. Its five strictly
growing return words have `(J,D,C)`:

    (1,1,1), (4,6,179),
    (7,11,16177), (7,11,12289), (7,11,9697).

For a chronological concatenation `u v`, retain

    J=J_u+J_v, D=D_u+D_v,
    C=3^J_v C_u + 2^D_u C_v,
    endpoint=(3^J n+C)/2^D.                            (3)

Its least odd residue and source-index residue are

    r_n=(2^D-C)3^-J mod 2^(D+1),
    r=((r_n+1)/2) mod 2^D.                             (4)

Let E_k be the sources for which the first k such returns all grow. The
five full return words are prefix-free, so their k-fold source cylinders
are disjoint. Applying (1) to every word gives the exact all-k identity

    mu(E_k)=q^k+A_k,
    q=1/4+1/4^6+3/4^11=1049603/4194304,
    A_k=sum_(w in B^k, r_w>0) w(r_w).                  (5)

The symbol `w(r_w)` in the last line is the Elias source weight, not a
word length. The first term is a simple recurrence; the second depends on
the exact original-source residues (4). It cannot be recovered just from
division clocks, slopes, or the number of words. A_k need not be monotone;
the full mu(E_k) is monotone because the E_k are nested.

The checker constructs every word through six returns (15,625 at the last
level), verifies (5), and independently executes the native queue on source
indices 1..1023. The omitted source mass is exactly 1/1024, giving an
independent enclosure of each exact cylinder sum.

Because q<1, the following are equivalent:

1. no positive integer stays forever in the five growing-return branches;
2. mu(E_k) tends to zero;
3. A_k tends to zero;
4. the least surviving positive source mu_k tends to infinity.

The equivalence with the least source is inherited from
[LEAST_SOURCES](../ternary_integer_returns_20261029/LEAST_SOURCES.md).
Equation (5) isolates a specific arithmetic term for the remaining task;
it does not prove that term decays. Even proving these equivalent statements
would still owe the growth debt before escape from this particular language.

## 3. What a decreasing natural-number certificate would mean

For a deterministic map T with an absorbing certified target, termination
from every state is equivalent to existence of a rank R into N with
`R(Tx)<R(x)` outside the target. In one direction there is no infinite
strictly descending natural sequence. In the other, the first hitting time
itself is such a rank. The converse is an equivalence, not an independent
construction of a proof: calculating or bounding that hitting time is the
original difficulty.

A source-preserving version allows an independently authenticated finite
budget B(n) and state `(n,current,t)`. Until acceptance, use `B(n)-t`; native
execution changes current and t, but not the original n. To use this as a
proof, establish that acceptance occurs before the budget expires. One
sufficient route is an independently proved atomic survivor bound

    mu{tau>t} <= b_t, with a computable b_t -> 0.

Choose B(n) with `b_B(n)<w((n+1)/2)`. The source cannot remain a survivor
at that time, because its one atom already exceeds the entire upper bound.
For example, the hypothetical bound `b_t=2^-t` would give budget
`2 floor(log2 m)+2`. But it is FALSE for native odd Collatz root-hitting
time: the finite source bank already gives survivor mass
`410285/2097152 > 1/8` at odd clock three. Source 27 alone retains weight
1/128 through odd clock forty and reaches 1 at odd clock 41 (division
clock 70). The positivity principle cannot be used to invent a rate.

Block horizons, odd steps and divisions are different clocks. A countdown
in completed returns yields native termination only when every block is
itself total; the inherited cheaper queue supplies that part. It does not
turn a locally decreasing return into descent below the original source.

## 4. No finite-phase residue-polynomial rank at every odd step

**Theorem.** Fix a modulus M, a finite phase set S and an initial phase.
Allow the phase update to depend on the phase, current residue mod M and
actual valuation a. Suppose a proposed rank on the reachable states is,
for each phase and residue, a rational polynomial in the current positive
integer and is nonnegative integer-valued there. Such a rank cannot strictly
decrease at every accelerated `3n+1` step before reaching 1.

**Proof.** Set L=|S| and consider the unbounded family

    n_0=2^(L+1) M t-1,   t=1,2,... .

Its first L valuations all equal one, and

    n_j=3^j 2^(L+1-j) M t-1,  0<=j<=L.                (6)

All these values exceed one and are -1 mod M. The phase sequence is
therefore the same deterministic sequence for every t. Among its L+1
states some phase repeats at indices i<j. Strict decrease would imply
`P(n_j)<P(n_i)` for the same polynomial P, for every t. But n_j>n_i.
Nonnegativity on the unbounded reachable progression forces P either to
be constant or to have positive leading coefficient. It is then eventually
nondecreasing, contradicting the strict inequality. QED.

The same proof excludes any finite-phase ranks eventually nondecreasing
on those residue progressions. It applies also at every cheaper-return
boundary: valuation one is already a complete growing return. It does NOT
exclude nonmonotone arithmetic ranks, source-labelled budgets, unbounded
history, or macros that wait for a separately certified reference descent.
The theorem is an obstruction to a declared class of proofs, not to Collatz.

The extra binary bit in `2^(L+1)` is essential. A scratch checker using
`2^L` failed on the last valuation; the explicit native test caught and
repaired that indexing before any result was published.

### A successful unbounded coordinate, with the unpaid restart visible

The inherited INTEGER_RETURNS proof already supplies a useful positive
control: let `f(3)=32`, and `f(n)=n^3` otherwise. While the queue is active,
the integer potential `f(n)4^Q` contracts by at most `125/128` per native
odd step, for every positive odd n and every positive integer Q. At terminal
states it may be set to zero. This proves termination of one whole return.
The queue is unbounded across starting sources, so it escapes the finite-
phase obstruction above. Restarting a completed queue raises its unabsorbed
potential by a factor four; the shifted potential in the incoming proof
has a different restart cost, `3/2`. Neither restart is free. Thus a concrete
next construction is an unbounded source account that pays these resets,
rather than another fixed collection of residue-polynomial phases.

## 5. The next concrete questions

* Bound the residue-dependent atom sum A_k in (5) without assuming the
  sources have already escaped. A harmonic/Haar tail bound is insufficient.
* Find a return grammar whose atomic escape certificate also pays its
  accumulated growth below the fixed original source.
* Search ranks carrying an unbounded arithmetic coordinate, then prove a
  local update inequality. Finite phases alone cannot defeat (6).

The `5x+1` control must remain visible: `13->33->83->13` contributes the
permanent source mass 69/2048 under these same weights. No source-preserving
re-encoding removes that obstruction.
