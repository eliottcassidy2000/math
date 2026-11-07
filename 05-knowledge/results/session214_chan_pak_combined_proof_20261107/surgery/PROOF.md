# Moving the comb exit: an all-order obstruction at 7,477

**Status: PROVED reduction plus FINITE-EXACT certificate, pending independent audit.** The exclusion is only from the explicitly defined class. The unrestricted Chan–Pak spectrum remains OPEN; an actual tournament with 7,477 paths is supplied below. No literature-priority claim.

## 1. The architecture and its relationship to Session213

Let H(T) denote the number of directed Hamiltonian vertex words of a tournament T.
Let F be the class of tournaments admitting a presentation with a transitive chain
x_0,...,x_(L-1) and at most two additional vertices a,b, where a has at most one
out-neighbour in the chain. The incidences of b to the chain and the a-b arc are
arbitrary. Zero or one additional vertex is permitted, and the empty chain is
permitted. The restriction concerns a's **chain** out-degree, not its total
out-degree; when a->b, a can have total out-degree two.

This contains both orientations of Session213's binary cycle comb. Reversing
its unique old exit a->x_1 to x_1->a and then reversing x_r->a to a->x_r relocates
the exit to any r. The present class is larger still: b need not beat x_0,
there is no distinguished x_1 attachment bit, and the exit can be absent.
It therefore tests a genuine route outside induced retention, not merely
another extension menu.

**Theorem.** The least positive odd integer other than 7 and 21 absent from F
is **7,477**. More precisely:

1. Every positive odd integer below 7,477 except 7 and 21 is H(T) for T in F.
2. No tournament in F, of any order, has H(T)=7,477.
3. A two-apex tournament in which a has exactly two chain out-neighbours has
   H(T)=7,477. Thus raising this chain out-degree allowance from one to two
   is necessary and sufficient for this particular target within the
   two-apex presentation hierarchy.

The exclusions of 7 and 21 in all tournaments are inherited from THM-1370;
they are not needed to prove the new class exclusion of 7,477.

## 2. A finite-to-all-order principle for a prime target

Two elementary facts drive the reduction.

**Hereditary closure.** Every induced subtournament of a tournament in F is
again in F. Restrict the transitive chain and retain whichever of a,b remain.
If only b remains, it is the unrestricted additional vertex. Deletion cannot
increase a's number of chain out-neighbours.

**Strong-component product.** If the strong components of T are ordered
T_1=>...=>T_k, every directed Hamiltonian word traverses each component in
that order and in one uninterrupted block. Hence

    H(T) = product_i H(T_i).

The product law is inherited standard tournament structure (THM-1370 and
THM-4094); the displayed block argument gives the needed direction directly.
A nontrivial strong tournament has more than one Hamiltonian path, for
example by taking the distinct rotations of a Hamiltonian cycle (Camion,
also inherited). Alternatively the cited bound below implies this directly.

Now 7,477 is prime: trial division through 86 excludes every possible divisor.
If H(T)=7,477, exactly one nontrivial strong component S has H(S)=7,477;
all remaining components have H=1 and must be singletons. Heredity gives S in F.

Use the **CITED** theorem of Busch, in the exact form recorded in
`01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`:

    H(S) >= f(n),
    f(3q)=3*5^(q-1),  f(3q+1)=5^q,  f(3q+2)=9*5^(q-1),

for strong S of order n>=3. This function is increasing, and f(18)=9,375>7,477.
Thus S has at most 17 vertices. The external theorem is inherited, not
re-proved or inferred from the finite tournament census.

## 3. Why the finite census includes every possible core

The actual census uses only the especially simple presentations

    chain length 1<=L<=17;
    a has exactly one chain out-neighbour x_r;
    both a-b orientations;
    all 2^L cuts of b.

Every S in F of order at most 17 can be embedded, without changing H, in a
presentation in that census, except the transitive cases which plainly have H=1.
The padding is always by global sink vertices, whose deletion leaves a
bijection of Hamiltonian words.

- If both a,b are present and a has exactly one chain exit, L<=15 already.
- If both are present and a has no chain exit, add a global sink z. It becomes
  a chain vertex after the old chain, with a->z and b->z. Now a has exactly
  one chain exit and the chain length is at most 16.
- If only one unrestricted apex b is present, add two successive global sinks
  a,z. Use z as the last chain vertex and a as the restricted apex. It beats
  exactly z, and b->a. The chain length is at most 17.
- If the single existing apex is a restricted a, it can be renamed as the
  unrestricted b and the preceding case applies. If no apex remains, S is
  transitive.

When a core already uses both apexes, strongness rules out an empty chain.
The small empty-chain or singleton boundaries are transitive and irrelevant
to the positive prime target. This padding argument is only a way to bound
an exact decision procedure; the padded tournaments need not be strong.

Consequently absence in the stated finite census proves absence from F
at every order. The proof does not assume a maximum order for arbitrary
nonstrong presentations, which may have unlimited source and sink padding.

## 4. Exact counting and audit surface

The census contains

    sum_(L=1)^17 2 L 2^L = 8,388,612

labelled parameter presentations, with no isomorphism or symmetry pruning.
The accepted odd counts through 7,477 are exactly all such counts other than
7, 21 and 7,477. There are 3,736 accepted values. These finite assertions
are reproduced by `verify.py` and its deterministic `RESULT.json`.

The counting transfer is **inherited**, in a short independent implementation,
from Session191's
`05-knowledge/results/session191_tournament_growth_20261104/two_apex/TWO_APEX.md`,
section 2. A Hamiltonian word has exactly one of the forms

    I,a,J,b,K;     I,b,J,a,K,

where I,J,K are increasing chain runs, possibly empty. Read chain vertices
in increasing order and assign each to one run. For each apex order retain
18 states: the empty/good-last/bad-last condition for I and J, and whether
K has started. The first vertex of J must be beaten by the preceding apex;
the first vertex of K must be beaten by the other apex. At the end, nonempty
I and J must point to their next apex; an empty J requires the corresponding
apex arc. This is a bijection, with no omitted or multiply counted words.

For fixed L,r and the apex arc, H is a multilinear polynomial of degree at
most two in the b-cut bits. In any Hamiltonian word, b has at most two
neighbours, so its validity depends on at most two such bits. Accordingly
values at the empty cut, every singleton and every pair determine the entire
polynomial. This attachment-quadraticity mechanism is inherited from
Session203. `coeff` obtains these coefficients using the run transfer;
`scan` traverses every Boolean cut in Gray order with exact integers.

Verification controls do not stand in for the all-order reasoning:

- Direct subset DP verifies every cut and both apex orientations for chain
  lengths 1 through 6 and zero or one restricted exit.
- Literal permutations give a third count through chain length 4.
- Direct transfer replays one physical witness for every accepted value,
  independently of its Gray traversal and quadratic interpolation.
- The two-exit 7,477 witness below is checked by transfer and subset DP.
- Checks use explicit exceptions and remain enabled with `python3 -O`.

The finite census is not a classification of tournament isomorphism types;
it is complete for these parameter presentations and this bounded scalar
range, with the padding proof supplying the needed all-order scope.

## 5. Positive control just beyond the barrier

Take a chain of length nine, set a->x_0,x_3 and x_i->a for the other chain
vertices, and set a->b. Let b beat precisely x_0,x_1,x_2,x_6,x_7. In chain
index order followed by a,b, the outgoing bit rows are

    [510,1020,1016,1520,2016,1984,896,768,1536,1033,199].

The subset Hamilton-path recursion and run transfer both give H=7,477.
This is an ordinary physical tournament, not a virtual subtraction of counts.
Since b beats x_0, it is obtained from an actual Session213 forward-apex comb
with the same b-incidences by three old-arc reversals: remove a->x_1 and add
a->x_0 and a->x_3. No claim is made that three is the minimum number of
arbitrary old-arc reversals from any comb; the sharp statement is the
chain-exit threshold within the stated presentation hierarchy.

## 6. Interpretation and remaining frontier

The closest proved mechanism is the binary cycle compiler. Its canonical
hostile is mistaking marked endpoint counts for total paths. The corrected
near miss here is the belief that merely relocating the bottleneck exit
could eliminate all total-path omissions. The least-used sidecar is prime
factorization of the output together with strong-component heredity.

The connection contract is explicit: comb -> old-arc surgery -> at-most-one
chain-exit class; the transitive spine and arbitrary second-apex cut survive,
the exact cycle code is discarded, and strong-component order is the sidecar
that converts a bounded failed search into an all-order certificate.

This result identifies a concrete target requiring more expressive apex
incidences than sole-exit relocation. It does not rule out two chain exits,
changing chain arcs, target-dependent general surgery, or unrestricted
Chan–Pak completeness. It neither proves a new ambient forbidden value nor
extends Session213's density-zero theorem to this larger class.

## 7. Reproduction

    python3 verify.py > RESULT.replay.json
    python3 -O verify.py > RESULT.optimized.json
    cmp RESULT.json RESULT.replay.json
    cmp RESULT.json RESULT.optimized.json

Run from this packet directory. Python 3 standard library only.
