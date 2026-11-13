# Late Session 317: terminal readouts and contact capacity

**Status: source-pinned bounded review of inherited claims.** This note reads
five incoming files, checks their typed connections to the present packet, and
retains their stated scope. It does not rerun their replays, independently audit
the entire terminal compiler, or promote a new theorem. Earlier audit files stay
frozen. Chan–Pak remains open.

The read source is commit `20c9c7e0a500077380bc172d18f332d6c446caca`, which
was `origin/main` when inspected. All paths have prefix
`05-knowledge/results/session317_exact_anchor_atlas_20261113/`.

| File | Git blob |
|---|---|
| `README.md` | `2366020c4287fde647cc81b9e0c91b6e68a6c3f3` |
| `COALESCENT_MENUS.md` | `bbff16f1faf20e26b5e91d09f97b4e3c764d8dd4` |
| `SLOT_INVARIANTS.md` | `338f0edc213b457c9076f21f96d0df515d263f02` |
| `CONTACT_CAPACITY.md` | `df827febf649eb236ee2abd245a68b2c85643517` |
| `TWO_END_APERTURE.md` | `416edab6a7e23a4ef43fe17cd2e13c1a53db5c74` |

These exact versions can be recovered with `git show COMMIT:PATH` without
rebasing or reading a stale worktree copy. The README's other terminal-atlas,
inversion, return-star and radix files were not included in this bounded review;
its descriptions of them are retrieval leads, not additional reviewed proofs.

## 1. Same exact H can have arbitrarily many next readouts

`COALESCENT_MENUS.md` gives the exact operation `J(T,U)=C3[T,U,1]` for nonempty
blocks. With `P=Q_T` and `R=Q_U`, its slot polynomial is

    Q_J(z)=z[P(z)R(z)+sum_(epsilon=±1)
                    (P(z+epsilon)−P(z))(R(z+epsilon)−R(z))].

Its H readout is `(Q_T(2)−H(T))(Q_U(2)−H(U))+2H(T)H(U)`.
For r>=1, let T_i be an ordinal sum of i directed triangles and 2(r−i)
singletons. The source retains the ratios

    H(T_i)=3^i,
    Q_T_i(2)/H(T_i)=4^r,
    Q_T_i(3)/H(T_i)=11^i 9^(r−i).

For L>=r, the strong outputs `Y_i=J(T_i,T_(r−i)→TT_L)`, i=0,...,r,
have one exact H value

    H_(r,L)=3^r[(4^r−1)(4^r 2^L−1)+2],

and pairwise distinct Q(2). The only varying term in the latter is

    −2*3^r*4^r[2^L 11^i 9^(r−i)+3^L 11^(r−i)9^i].

The bracket decreases strictly because the ratio controlling its successive
changes is `(2/3)^L(11/9)^(2i−r+1)<1` for i<r. The constructions also have the
same order `5r+L+1`, directly by counting their blocks. This is an all-scale
mechanism for unbounded continuation ambiguity of **H alone**.

**Connection and boundary.** The terminal operation sees more than its current
scalar output. This is the same observer-closure issue as the current incidence
packet, but uses an explicit changing-size operator and selected Q evaluations,
not the current fixed-order dart matrix. It neither contradicts nor replaces the
full-Q closure computation: the menu deliberately varies Q(2), hence full Q.
The anchors H_(r,L) are sparse, not arbitrary target values. Moreover
`3^r | Q_Y_i(t)` at every integer t, and further cyclic wraps or ordinal singleton
padding preserve that common divisor. Changing a readout has not removed this
arithmetic obstruction to a prime-target grammar.

## 2. Preserved interfaces impose universal, not grammar-specific, barriers

`SLOT_INVARIANTS.md` uses `Q(t)=t^epsilon P(t²)` from odd-cycle packings. For
same-parity tournaments, equality at distinct integer slot counts t_i>=2 forces
`product_i(t_i²−1) | Delta H`. Different parity gives the corresponding squared-H
condition. In particular `Q(2) ≡ (−1)^n H (mod 3)`.

A +2 step preserving Q(2) must change parity and go from H=2 modulo3 to H=1
modulo3. Consequently **two successive such steps are impossible**, regardless
of the attachment grammar. Same-parity +6 is attained by a literal pair, but
cyclic return imposes further constraints. Under the actual cyclic operator

    (L p)(t)=t²p(t)+t[p(t+1)−p(t−1)],

persistent equality at t=2 forces recursively prescribed higher evaluations.
The source derives necessary H-gap multiples, after one/two returns, of
(60,90)/(3780,7560) for odd/even parity; for equal actual order they strengthen
to (120,90)/(7560,15120). These are necessary filters, not sufficient realization
or closure conditions. The nonzero forced evaluation sequence grows faster than
any fixed polynomial, obstructing indefinite return by a fixed difference p.

A useful bypass changes the preserved functional: source responses
`(H,Q2,Q3)=(99,1332,9153)` and `(103,1328,8877)` both wrap to H=1431, because the
next count reads H+Q2. Their new Q2 values are 23436 and 22860. This is a literal
coalescence at one exact target, not an infinite family of returns or a violation
of the fixed-Q2 obstruction. For a terminal recurrence the state must retain the
functional actually used by the next operation; preserving each summand is an
additional and sometimes unnecessary constraint.

## 3. Nonuniform contact is a limited arithmetic resource

`CONTACT_CAPACITY.md` gives two distinct bounds that should not be merged.
With uniform interblock arcs and k apex vertices, a physical Hamilton path uses
at most 2k marked block vertices adjacent to apexes. If
`g_i(s)=min_(|S|<=s) v_p(|Aut(B_i)_(S)|)`, then

    v_p(H(T)) >= min_(sum s_i<=2k) sum_i g_i(s_i).

The source's proof partitions paths by itinerary and actual contact names; the
product of pointwise stabilizers acts freely on each part. It exhibits strong
families where one apex genuinely removes four factors of three from an iterated
triangle rotation group. Thus a contact can consume several symmetry factors,
but its budget must be evaluated from marked stabilizers rather than vertex count.

Separately, for an **ordinal** block order and k apexes, an untouched block is
visited in at most k+1 runs. Therefore its contribution retains the divisor

    d_i(k)=gcd(F_j(B_i):1<=j<=min(k+1,|B_i|)).

At most 2k blocks are touched, so a lower p-valuation bound is the sum of the
`v_p(d_i(k))` after dropping its largest min(2k,r) entries. This survives even
when every block is rigid. The source's rigid H=9, Q2=84 block has d_i(1)=3;
a sharp one-apex family of r blocks has valuation exactly r−2. For r>=3 this
particular grammar therefore retains a factor3, excluding its large prime and
7p targets. This is not an unrestricted prime exclusion.

**Connection and boundary.** Endpoint incidence supplies the missing information
that full unpointed Q cannot recover automatically. One should carry pointed
slot counts or stabilizers at the actual contacts, together with untouched-block
slot divisors. The symmetry and arithmetic bounds cannot simply be multiplied
on the same factor; the source explicitly forbids that double counting.

The same note gives a separate homogeneous-substitution obstruction. If a strong
quotient Q on at least three vertices replaces one vertex by a nonempty n-vertex
U, then

    H(Q[v←U]) >= H(Q)H(U)+F_2(U) >= (n+2)H(U).

A two-run quotient itinerary supplies F_2(U), in addition to contiguous one-run
paths. Hence fixed-factor multiplication, such as by3 or9, cannot indefinitely
retain the whole growing old tournament as one homogeneous block. The bound is
uniform even when the quotient changes. Any proposed recurrence evading it must
change that physical interface; the statement does not prohibit nonuniform
contacts, interacting blocks, or other genuinely different constructions.

## 4. A new exact +2 aperture is not ready for repeated entry

`TWO_END_APERTURE.md` fixes a five-vertex terminal R with H(R)=9 and two ports
having the same optional-slot end response E=11 but Hamilton start counts1 and3.
For arbitrary nonempty M, the resulting strong counts are

    20 Q_M(2)+19 H(M),     20 Q_M(2)+17 H(M).

For transitive M=TT_L this gives adjacent odd counts
`20*2^L+19` and `20*2^L+17`. But switching the port changes the next two-slot
readout by

    Delta Q_out(2)=−128 H(M)−88 F_2(M)−8 F_3(M),

which is strictly negative for nonempty M. Thus the exact pair does not preserve
the required interface for a repeated fixed-Q2 +2 rule. The dependence on F_3 is
explained by two output slots and one apex permitting at most three runs through
M. The source determines the three coefficients from a triangular TT_1,TT_2,TT_3
basis and reports a TT_4 replay check; this review did not rerun that finite check.

## What changes in the present packet

The incoming work supplies concrete terminal operations and obstructions for
choosing a continuation-closed observer. It reinforces the need for jointly
retained incidence/slot data while showing that **which** combination is retained
may change from one operation to the next. Three quantities remain distinct:
the exact current anchor H, the next slot readout, and the physical contact data.

The cheapest decisive probes are: compare Q2 within an exact-H fibre; test the
mod3 fixed-Q2 readiness condition; compute the contact stabilizer/slot gcd; and
count how many runs an operation permits through the old block. None gives a
recurrence reconstructing unrestricted tournament incidence, a full target cover,
or an implication from the current n<=7 closure census to Chan–Pak.
