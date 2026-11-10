# A two-index gap anchor, its complete inverse, and a finite-bank density barrier

**PROVED:** all-tail normal form, a physical two-index even-gap family, a
complete target-dependent inverse, and a positive-density omission theorem
for its bounded input banks. The general tail formula and finite-bank
parameter-capacity principle are INHERITED; this is their explicit application
to the new pair. **INHERITED:** joint modular universality of
`(N,g(N))` and the four-block construction language. **FINITE-EXACT:** one
complete response fiber obstructs a proposed induction. The unbounded
family's coverage and asymptotic density, the four-block lemma, and Chan–Pak
remain **OPEN**.

This round resumes interrupted work and records the exact formulas before
trying to use them as an anchor. It does not turn a large finite coverage
prefix into an induction.

## 1. Inheritance and the live proof obligation

The source laws and the canonical source `(g(N), N mod 2, N, 1)` for `N>=2`
are inherited from `session199/hidden-structure/PROOF.md`. At `N=1`, use the
literal singleton state `(1,1,1,1)` and the auxiliary convention `g(1)=1`.
Every construction below starts with O, so the singleton boundary immediately
enters the ordinary nontrivial source class.

The existing four-block lemma, Section 5 of that source, asserts that every
odd `h>27781` has a source with at most three J letters and balanced readout
`G=D+2U=h`. It is OPEN, although the interval through 1,000,001 was exhaustively
verified in Session199. Its chronological caterpillar interpretation is also
inherited, not a new Session244 result.

The modular theorem is stronger than merely having no residue obstruction:
`session199/congruence-intervals/PROOF.md`, Theorem 1, gives every compatible
pair `(N mod 2m,g(N) mod 2m)` in every sufficiently long aligned binary
interval, with a common correcting-suffix length. Theorem 3 there gives
congruence-exact lower approximations with error `O_m(sqrt(h))`. Neither is
exact target placement. That same inherited zero-J image has a positive
density of omissions.

The general O/Z scalar formula is already proved in
`session205/finite-oz-tail-bank/ARITHMETIC_REDUCTION.md`. The stronger exact
parameter-capacity theorem is inherited from
`session206/parameter-capacity/PARAMETER_CAPACITY.md`, and modularly weighted
versions are in `session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md`;
see also THM-4888
(`01-canon/theorems/THM-4888-finite-oz-banks-mod3-law-involution-rigidity-residue-confinement.md`).
The g(M) form below makes the two-index substitution transparent; no priority
claim for the underlying tail image or finite-bank obstruction is made.

Concept board: pure-D synchronized pairs; exact last-block bit collapse;
valuation-density budgets; same-target response fibers; and literal graph
rebasing. Closest positive mechanism: the sibling composition packet's pure-D
pair. Hostile: 173 is an actual count (ROUND2), but not in the two-index family
below. Corrected near miss: modular surjectivity is inherited and does not
supply bounded-height exact selection. Least-used sidecar: the valuation of
the independent tail parameter, rather than only the available gap.

## 2. An exact scalar collapse for every retained-root binary tail

Let a chronological O/Z word have bits `epsilon_i`, with Z represented by 1,
and write

    A = sum_i epsilon_i 2^i,
    C = sum_(j<i) epsilon_j epsilon_i 2^(i-j-1)
        + sum_(i>=1) epsilon_(i-1) epsilon_i.

Trailing O letters do not change A or C and do not change the full count.
The inherited binary-block compiler gives

    D' = D - (A mod 2) E + A U - C d.                    (1)

There is an exact reduction of its quadratic term:

    C + (A mod 2) = A+1-g(A+1).                         (2)

**Proof.** Both sides are zero at A=0. Write the first chronological bit as
`a=A mod 2` and the remaining integer as `B=floor(A/2)`. If a=0, deleting
the initial zero leaves C unchanged. If a=1, its interactions with later
bits contribute `B+(B mod 2)`, hence `C(A)=C(B)+B+(B mod 2)`.
Substituting these two identities into the recurrence for g proves (2) by
induction: for `A=2B`, use the odd rule at `2B+1`; for `A=2B+1`, use
`g(2B+2)=g(B+1)`.

Consequently, with `M=A+1>=1` and `a=(M-1) mod 2`, the entire full-count
image of all finite O/Z tails from any fixed input state is exactly

    D_M = D+(M-1)U-(M-g(M)-a)d-aE
        = (U-d)M+d g(M)+(D-U)+a(d-E).                  (3)

For the empty tail set M=1, and (3) still gives D. Conversely every M is
realized by reading the binary digits of M-1 least-significant first,
using O/Z accordingly. Thus (3) is an exact image statement, not a relaxation
of the word constraints. It does not retain the final U, which also depends
on the discarded trailing O length; it is used only for the scalar count.

## 3. Physical two-index synchronized pairs

On every formal input `(D,E,U,d)`, the following actual graph words have
exact transfer matrices

    P = OJJJJOZ:
        (17D+40U, 3D+6U, 21D+50U, 3D+6U),

    Q = OJOZJOZ:
        (17D+42U, 3D+6U, 21D+50U, 3D+6U).             (4)

They differ by `(2U,0,0,0)`. Any common retained-root O/Z tail preserves
this pure-D difference. The proof of (4) is direct multiplication of the
physical source transfers; the verifier additionally constructs and recounts
literal graphs, including the singleton boundary. The words themselves and
their synchronization were supplied by the sibling composition lane.

Apply (4) to the canonical source with response N, then apply the tail
indexed by M. The lower count is exactly

    L(M,N) = (18M+3g(M)-4) g(N)
             +(44M+6g(M)-10) N,       M,N>=1,          (5)

and the upper count is `L(M,N)+2N`. Indeed the common endpoint and deletion
coordinate in (4) eliminate the parity term in (3). This is an independent
absolute-count parameter M beside the selected gap 2N, but the dependence
on M is still constrained by g(M).

**Complete inverse.** Every endpoint at most h satisfies

    MN <= floor((h-17)/40).                            (6)

To prove it, use `g(M),g(N)>=1` in (5):

    L(M,N) >= (18M-1)+(44M-4)N
            =40MN+17+(M-1)(4N+18)
            >=40MN+17.

Thus enumerate all positive M,N in (6) and test (5) and `L+2N=h`. This is
necessary and sufficient for membership in this family, with no graph-order
or tail-length cutoff. The smallest pair is 57,59 at M=N=1. In particular
61 and 173 are excluded at every order in this family, although both are
attained by other tournaments. For 173, (6) gives MN<=3, so its exclusion is
a five-pair arithmetic check, not a truncated graph search.

## 4. There is no congruence obstruction, at any modulus

Already the M=1 lower subfamily

    L(1,N)=17g(N)+40N                                  (7)

meets every odd residue modulo every even modulus 2m. This is a consequence
of the inherited joint-residue theorem, not a new modular universality
result. Write `g=2y+1`. The residue to solve is

    17g+40N =17+2(17y+20N).

Since `gcd(17,20)=1`, any odd target residue has a solution in y,N modulo m.
The inherited theorem supplies an actual N with those simultaneous residues
of N and g(N); N can be made arbitrarily large. For example, the identity
`17*(-7)+20*6=1` supplies explicit residue choices. Moduli that are odd are
handled by passing to twice that modulus.

This refutes any proposed fixed congruence obstruction for the full family.
It does **not** bound the N needed to hit a residue in terms of the desired
absolute count. In fact the inherited modular theorem's scale depends on
its modulus, so choosing a modulus larger than a target does not solve (6).
The all-tail family has both complete local residue access and the exact
small-target failures above.

### Fixed-input congruences and the prime lane

The union is locally universal, but individual input gaps can be highly
constrained. For fixed N put `q_N=gcd(N,3g(N))`. Both endpoints satisfy

    h = -g(N) (mod 2q_N).                              (8a)

There is no further fixed-modulus obstruction: modulo any even modulus B,
the attainable residues are exactly those congruent to `-g(N)` modulo
`gcd(2q_N,B)`. This is a modular-image statement, not equality of the
absolute integer image with that arithmetic progression.

To prove it, rewrite (5) as

    h = A M + b g(M) + c,
    A=18g(N)+44N, b=3g(N)+6N, c=-4g(N)-10N.

The inherited joint-residue theorem makes `(M mod B,g(M) mod B)` range
over all pairs whose second coordinate is odd. Thus the output image is
the coset `b+c` of the subgroup generated by A and 2b. Now

    gcd(A,2b)=2 gcd(9g(N)+22N,3g(N)+6N)
             =2 gcd(N,3g(N))=2q_N,
    b+c = -g(N)-4N = -g(N) (mod 2q_N).

The upper endpoint adds 2N, a multiple of 2q_N, and has the same modular
image. This is the exact specialization of the inherited general modular
image/capacity framework to the present pair.

For an actual prime target, an immediate stronger divisibility filter is
`gcd(N,g(N))=1`: otherwise that nontrivial common divisor divides both
endpoints, while (6)'s lower estimate makes the endpoints larger than it.
In particular `N=2^k+1`, k>=1, has `g(N)=N`, so every endpoint in those
input families is composite. Their gap universality is still genuine, but
these infinitely many inputs cannot anchor a prime target. N=1 prevents
this per-input filter from becoming an obstruction for the union.

## 5. A density budget for fixed seed banks

Equation (3) gives a convenient elementary valuation bound for some finite
banks, including a large bank of the new gap pairs. The stronger inherited
Session206 parameter-capacity coefficient is
`log(U/(U-d))/d`, with parameter count equal to that coefficient times X
plus `O_seed(sqrt(X))`. We give the following self-contained, weaker bound
so that the explicit application uses only rational estimates.

Fix a physical seed `(D,E,U,d)` with `U>d>0`. Let `S_seed` be its entire
retained-root O/Z-tail count image. Then

    limsup_(X->infinity) |S_seed intersect [1,X]|/X
      <= beta(U,d)
       := sum_(j>=0) 1/[2^(j+1)(U-d)+d].              (8)

**Proof.** Write `M=2^j m` with m odd. Since `g(M)=g(m)` and
`g(m)>=(m+1)/2`, including m=1, formula (3) implies

    D_M >= [2^j(U-d)+d/2]m + O_seed(1).

Thus at most `X/[2^(j+1)(U-d)+d]+O_seed(1)` positive odd m can contribute
at a fixed j. Formula (3) also bounds M by `O_seed(X+1)`, so only
`O_seed(log(X+2))` valuations contribute. Sum and use the convergent series.
Counting parameter values overestimates distinct scalar values, making the
bound valid even in the presence of collisions. The same proof yields an
`O_seed(log X)` error. For any fixed finite bank, sum (8) over its seeds.

For either seed in (4) at response N, its coordinates give

    U-d =18g(N)+44N,
    d   = 3g(N)+ 6N.

Consequently its density budget is

    beta_N = sum_(j>=0)
       1/[2^(j+1)(18g(N)+44N)+3g(N)+6N].              (9)

Let `A_K` consist of BOTH endpoint families in (5), with **arbitrary M>=1**
and `1<=N<=K`. Then

    upper density(A_K) <= 2 sum_(N=1)^K beta_N.         (10)

The count differences 2N do not affect the asymptotic bound. A useful coarse
bound, discarding positive g(N) terms, is

    beta_N < (1/N) [1/94 + sum_(j>=1) 1/(44*2^(j+1))]
           = (91/4136)/N.                             (11)

For every `K<=32768`, the harmonic number `H_K<11`. An elementary exact
certificate is `H_32768<=H_32+10 log 2 < H_32+347/50 <11`.
The verifier checks the rational last inequality and `log 2<347/500` using
its convergent `2*atanh(1/3)` series with a rational tail bound.
It follows from (10)-(11) that

    upper density(A_32768) <1001/2068 =91/188 <1/2.            (12)

Therefore the lower limiting fraction of positive odd targets omitted by
this complete bank is at least

    1-2*(1001/2068) = 33/1034 =3/94 >0.0319.                 (13)

This is an all-order theorem with unlimited retained-root tail length and
all 32,768 chosen gap sizes. It is **not** a finite experiment and not a
nonexistence theorem for tournaments. It also does not extend to unbounded
N: the estimates in (11) are harmonic, so summing them over all N diverges.
No positive-density omission or cofinite coverage is proved for the full
two-index family (5).

As with the stronger inherited capacity criterion, a proposed finite input
bank whose summed beta budget is
less than 1/2 cannot give cofinite odd coverage by retained-root O/Z tails.
This is a necessary test for such a proof strategy, not a sufficient test
for success when the sum is at least 1/2.

## 6. A complete hostile fiber for the simplest balanced induction

A possible induction would rebase two sources to the same balanced target
`G=D+2U`, then append O to obtain adjacent outputs. For two equal-G sources,

    G(s.O)=D+4U=G+2U.

Their O-outputs differ by2 if and only if the two original U values differ
by1. The complete at-most-three-J fiber at G=27783 contains 33 distinct
profiles, but only these 20 distinct response values:

    9261, 9288, 9441, 9478, 9801, 9849, 10041, 10115,
    11160, 11200, 11350, 11502, 11526, 11534, 11538,
    11570, 12348, 12360, 12376, 12564.

No pair is consecutive. The verifier exhausts all 216,030 source profiles
with at most three J letters and G<=27783; strict growth of G under every
letter proves completeness at all word lengths. It stores the entire
33-profile fiber and minimum J costs in ROUND4.json.

Thus the exact proposed invariant, “every covered target in the proposed
cofinite range has consecutive U values in its same-target <=3-J fiber,”
is false already at its first target. This does not refute the four-block
coverage lemma: that lemma asks for one witness per target, not this extra
induction invariant. Other rebases or different suffixes are not excluded.

## 7. Consecutive root changes admit a literal matching presentation

There is a useful actual graph rebase, rather than only a scalar equality.
Start from ANY rooted tournament `(T,r)`, and append k consecutive J births
`v_1,...,v_k`. Keep the old graph fixed and reorder only the newborns as

    v_2,v_1,v_4,v_3,...,

leaving `v_k` last when k is odd. Relative to the ordering with all old
vertices first, the new backward edges are exactly

    {r,v_1}, {v_2,v_3}, {v_4,v_5}, ...                 (14)

with endpoints present. They form a matching. Among old vertices retain
exactly T's own edges; this statement does not require T itself to be
transitive or to have no pre-existing reversals.

**Proof.** Before reordering, the new backward edges are the path
`r-v_1-v_2-...-v_k`. Each swapped adjacent pair reverses its order, turning
its path edge forward. No other newborn pair changes relative order.
All old vertices keep their positions, so old edges and old-new incidences
are unchanged. The remaining path edges are exactly (14), which are
pairwise disjoint. This is an edge-by-edge isomorphism under relabeling.

The marked root must also be transported. For positive even k, `v_k` is
second-last in the new order, not last; for odd k it is last. Resetting the
mark to the last position would not preserve the rooted state. The verifier
checks 360 such literal rooted isomorphisms across all source prefixes of
length at most3 and k<=8.

This removes unnecessary reversal depth from a consecutive J run. It does
not flatten arbitrary separated J blocks: moving a spine vertex across its
intervening Z leaves changes additional edge orders. No general forest or
depth-two image equivalence is claimed from the earlier small experiments.

## 8. Reproduction, surviving invariant, and stopping point

From this directory:

    python3 -I -B verify_round4.py > ROUND4.replay.json
    cmp ROUND4.replay.json ROUND4.json
    python3 -I -B -O verify_round4.py > ROUND4.optimized.json
    cmp ROUND4.optimized.json ROUND4.json

The checks include every O/Z block through length11 on four seed profiles,
16,384 two-index identities and inverse bounds, 32 independent literal graph
profiles, modular closure controls through modulus64 and 1,024 fixed-input modular
images, the exact density
budget, a complete target fiber, and 360 actual matching rebases. The
infinite statements follow from the displayed proofs, not these cutoffs.

The sufficient exact anchor condition for this new family is now explicit:
for each target h, one must find positive M,N satisfying (5) or its upper
variant, inside the finite region (6). The original four-block sufficient
condition remains its exact matrix-language equation with <=3 J letters.
Neither condition has been proved to succeed for every large odd h.

The next proof must use an invariant stronger than residue surjectivity,
stronger than merely having an independent tail parameter, and different
from the refuted same-target-consecutive-U condition. The bounded-bank density
budget proves why fixing a large but finite set of the present gap inputs
and then allowing arbitrarily long O/Z tails does not supply such an
invariant. These are precise stopping points, not claims that all other
anchoring or multi-root constructions are impossible.
