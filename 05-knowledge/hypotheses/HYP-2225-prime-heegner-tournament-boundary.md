# HYP-2225: Prime-Heegner Windows Are Spine/Interior Boundary Carriers

**Status:** OPEN transfer hypothesis with exact arithmetic scout.

## Claim

For the Euler-Rabinowitsch quadratics

```text
f_p(x) = x^2 + x + p,
```

the lucky primes

```text
p in {2,3,5,11,17,41}
```

have an exact endpoint boundary:

```text
f_p(x) is prime for x=0..p-2,
f_p(p-1) = p^2.
```

The Heegner side channel is

```text
d = 4p - 1,
f_p(x) = ((2x+1)^2 + d)/4.
```

Thus the Euler lucky primes are the projection of the Heegner class-number-one
list through `d=4p-1`:

```text
d in {7,11,19,43,67,163}
<-> p in {2,3,5,11,17,41}.
```

The small Heegner numbers `{1,2,3}` are the exceptional part outside this
Euler-family shape: `d=1,2` are not `4p-1`, and `d=3` gives the degenerate
`p=1`.

The tournament-facing carrier is the indexing split:

```text
zero-based prime run length = p-1 = Hamiltonian spine length on p vertices;
interior inputs x=1..p-2 have count p-2 = Moon c3 floor for strong p-tournaments;
endpoint x=p-1 is a forced square sink p^2;
after fixing the spine, the free arc budget is C(p-1,2).
```

So the user's `n-2` signal is real if the endpoint/base terms are retained:
`x=0` is the source/base prime `p`, `x=1..p-2` is the interior run, and
`x=p-1` is the square sink.

## Evidence

S649 searches all primes `p<=500` for the exact boundary
`first_composite_x=p-1` and `f_p(p-1)=p^2`.  It recovers precisely:

```text
[2, 3, 5, 11, 17, 41]
```

The Heegner projection table is:

```text
d=1,2   -> outside d=4p-1 Euler-family shape
d=3     -> p=1 degenerate
d=7     -> p=2
d=11    -> p=3
d=19    -> p=5
d=43    -> p=11
d=67    -> p=17
d=163   -> p=41
```

The boundary/tournament ledger is:

```text
p=2:  run=1,  interior=0,  spine=1,  Moon floor=0,  off-path fiber=0
p=3:  run=2,  interior=1,  spine=2,  Moon floor=1,  off-path fiber=1
p=5:  run=4,  interior=3,  spine=4,  Moon floor=3,  off-path fiber=6
p=11: run=10, interior=9,  spine=10, Moon floor=9,  off-path fiber=45
p=17: run=16, interior=15, spine=16, Moon floor=15, off-path fiber=120
p=41: run=40, interior=39, spine=40, Moon floor=39, off-path fiber=780
```

The exact arithmetic identity is classical Heegner/Rabinowitsch structure; the
new repo-facing object is the source/interior/sink carrier that makes the
`p-2` count visible without forgetting the `p-1` spine.

## Tournament Analysis

The session explicitly challenged the vertex set.  Candidate vertices included
primes, Heegner discriminants, input positions, residues, Hamiltonian-path
arcs, directed 3-cycles, off-path arcs, and proof obligations.

The chosen Tournament Analysis uses carrier lenses, because they preserve the
boundary relation and expose what scalarization destroys:

```text
endpoint_square_failure_p2
heegner_class_number_one_side_channel
norm_line_d_equals_4p_minus_1
interior_length_p_minus_2_moon_floor
spine_length_p_minus_1
off_path_fiber_choose_p_minus_1_2
input_positions_as_vertices
raw_lucky_prime_list
```

The majority tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
hamiltonian_paths=1
```

The ranking puts the square endpoint and Heegner side channel above raw lucky
prime numerology.  That is the guardrail: the `n-2` echo is useful only while
the norm/class-number side channel and the spine/fiber tournament side channel
remain attached.

## Proof Use

The proposed transfer move is an endpoint-square carrier:

```text
long clean prefix + forced square/failure endpoint
-> source/interior/sink decomposition
-> compare interior count to strong-tournament floor invariants
-> compare full prefix length to fixed Hamiltonian spines
-> retain off-path deformation fibers as the hidden side channel
```

This can guide future work in three places:

1. LRC: treat `n-2` appearances as interior proof obligations, while `2n-1`
   remains the discriminant/root observer.
2. Unit distance: separate unit-spine length from bulk/failure endpoints in
   construction carriers, as in S648's Moser fixed quantum.
3. H-gap work: keep Moon's `c3>=n-2` as a lower-bound floor, not as a raw
   scalar equality.

## Limitations

This hypothesis does not prove the Heegner theorem or unique factorization.
It records the exact finite arithmetic bridge and adds a tournament carrier
dictionary that may be useful for transfer.  Any future proof must still retain
the class-number-one/norm side channel rather than replacing it with a
tournament count alone.

**See:** `04-computation/prime_heegner_tournament_boundary_s649.py`;
`05-knowledge/results/prime_heegner_tournament_boundary_s649.out`;
`07-reflections/prime-heegner-tournament-boundary-s649.md`; HYP-2224,
HYP-2223, HYP-2222, HYP-2217, HYP-2215, HYP-2200, THM-115.
