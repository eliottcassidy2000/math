# Strong arithmetic resets and a smaller critical-target grammar

**PROVED, scoped matrix identities and word reductions; FINITE-EXACT short-word
classification and literal witnesses.** This is a bounded follow-on to
[PROOF.md](PROOF.md), using exactly the same inherited thirteen generators.
It gives no coverage theorem and no complete classification of long products.

## Three different notions of a reset

The individual ranks, in inherited generator order, are

    1,3,2,2,2,2,3,3,3,3,3,3,3.

Only R=M1 has rank one. The other singular generators M3,M4,M5,M6 have rank two.
They do not all create permanent ordinal-sum cuts. In particular, M3,M4,M5
followed by C=M2 produce strongly connected tournaments with counts 5,9,17.
Their temporary sink vertices can be absorbed into later cycles because the
root remains on the other side of the cut.

Put N=M6. Its square is

    N=[1 4 4; 0 2 2; 1 2 2],
    N^2=(5,2,3)^t (1,4,4).                              (1)

The row ell=(1,4,4) is *N's own first row*. Therefore ell*s is an actual
tournament count H(Ns), not merely an arbitrary positive scalar divisor.
The column r=(5,2,3)^t is N*s0, hence also an actual tournament state.

For any prefix state s and any suffix matrix P, a word containing NN gives

    H(P N^2 s)=H(Ns) H(P r).                            (2)

Both factors are actual tournament counts at least five: H(Ns)=h+4E+4O>=5
because h>=1 and O>=1; H(Pr)>=5 by monotonicity. Thus such a word never yields
an odd prime or seven times an odd prime. A nontrivial factorization of 7p
into actual tournament counts would require the forbidden factor 7.

Unlike R, this arithmetic reset is compatible with strong connectivity. The
word NN produces the strongly connected order-seven tournament

    rows=(124,97,114,118,2,16,48), root=4,
    state=(25,10,15), coordinate gcd=5.                 (3)

There is no graph-level ordinal-sum factorization into two H=5 components.
An explicit isomorphism to the absorption lane's strong H=25 endpoint is

    old vertex -> endpoint vertex: (0,4,5,6,1,2,3).

It sends the root to endpoint vertex 1 and the rows to
(110,16,2,6,13,30,62). This identity was independently checked by the fibre
agent. Hence the count-preserving two-arc absorption switch actually reaches
this strong rank-one arithmetic reset, a concrete connection between the two
lanes. It does not by itself make the switch reusable in arbitrary contexts.

## A finite forbidden-pattern grammar for the critical lanes

For a in {1,2,4}, write

    X_a=[1 2 0; a a a; a a a],

so these are M3,M4,M5. Put

    ell_a=(1+2a,2+2a,2a)=e1^t C X_a.

Direct multiplication gives, for a,b in {1,2,4},

    X_b X_a=(1,b,b)^t ell_a,
    N C X_a=(5,2,3)^t ell_a.                            (4)

All factors have literal witnesses. The return ell_a*s is H(C X_a s)>=5.
The rays (1,b,b) are respectively C*s0, M3*s0, M4*s0 and have H=1.
The ray (5,2,3) is N*s0 and has H=5.

**Theorem.** Every critical count p or 7p represented by the thirteen-letter
language has a representation using the twelve letters other than R and
avoiding all the construction-order subwords

    X_a X_b,     NN,     X_a C N.                       (5)

**Proof.** First remove R by the prime/seven-prime extraction in PROOF.md.
Equation (2) excludes NN. The second identity in (4) excludes X_a C N,
because both factors of any continuation are actual counts at least five.

For X_a X_b, the first identity in (4) gives an attained return >=5 times
the suffix count starting from (1,b,b). If the suffix is nonempty, its first
letter increases H to at least 1+2b>=3; subsequent letters cannot decrease it.
This would again be a product of two nonunit attainable counts, impossible
for a critical target. Thus such a pair can only be terminal. Replace that
terminal X_b by C: the first row is exactly ell_a in both cases, so the
scalar H is unchanged and the forbidden pair disappears. QED.

This last replacement preserves H and gives an actual smaller tournament;
it need not preserve the full rooted state or future readiness. The theorem
is a representation reduction for critical targets, not an interchangeable
context identity. The resulting grammar is a sufficient search language if
the inherited thirteen-generator language can supply the critical lanes.
Its coverage remains OPEN.

## Coordinate gcd is a useful filter, not an attained factor

Integer matrices preserve divisibility of the three coordinates, so their gcd
cannot lose a divisor as a word grows. This is a valid necessary target filter.
But a coordinate gcd need not itself be an attained tournament count. The
literal R-free word (M7,M10,M13) gives the strongly connected order-twelve
witness

    rows=(4092,4081,4082,4022,3968,3856,3384,3936,
          2048,320,768,1536), root=10,
    state=(1757,994,1372), gcd=7, 1757=7*251.            (6)

Thus one of the required seven-prime targets has coordinate gcd seven,
although seven itself is forbidden as a tournament count. The failed
implication is “a common scalar divisor supplies a tournament factor.”
The repaired criterion in (1) and (4) retains both a realizable return row
and a realizable output ray. Neither rank one nor coordinate divisibility
alone supplies that information.

## Finite audit

Run

    python -B 05-knowledge/results/strategy_fibres_absorption_20261113/orbits/singular_check.py

The audit verifies all identities in (1) and (4). Among R-free words through
length three, the thirteen minimal rank-one subwords are exactly the nine
X_a X_b pairs, NN, and the three X_a C N words. “Minimal” means containing no
shorter rank-one contiguous subword. This finite classification is not asserted
for arbitrary word lengths.

It also checks all 22621 R-free words through length four. Among them, 6121
have prime or seven-prime counts; 471 require the terminal XX replacement.
The exclusions and replacements all pass. Independent path-subset enumeration
checks the full states of (3) and (6), and direct reachability checks establish
strong connectivity. The receipt retains their literal adjacency rows.

The contrast is now exact: R is a permanent graph cut; NN is a strong
arithmetic product; a gcd-seven state can lie on the essential seven-prime
lane without an attained seven-factor. Any larger reset classification should
retain these distinctions.
