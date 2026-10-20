# Marked tournament trijection, two views of ten, and reset composition

20 October 2026. **PROVED scoped elementary maps and composition laws;
FINITE-EXACT controls.** The inherited constructions are cited at their
point of use. No literature novelty, LRC(14), global Collatz, or unrestricted
tournament-metagraph connectivity claim is made. The closed blue `(n-4,3)`
endpoint layer belongs to the concurrent parent session and is not pursued
here.

## The exact trijection and its first lost coordinate

For a tournament with a marked Hamiltonian path `P=(v1,...,vn)`, number
vertices along P. For each nonadjacent position pair `i<j`, put an edge
`{i,j-1}` in a simple graph G on the fixed ordered set `[n-1]` exactly when
`v_j -> v_i`. These same bits orient a tournament U on `[n-1]` backwards
against its reference order. Thus

    isomorphism classes of marked pairs (T_n,P)
      <-> simple graphs G on ordered [n-1]
      <-> tournaments U on ordered [n-1].

Each inverse restores every arc, with the adjacent forward path forced.
The graph's order is essential: at n=4 the isomorphic one-edge graphs
with edges12,13,23 lift respectively to `C3 -> sink`, the strong
four-tournament, and `source -> C3`. Forgetting the original path instead
has exactly `H(T)/|Aut(T)|` fibres, inherited from LEM-003.

The translated predicates are explicit. T is transitive iff G is empty;
T is strong iff the intervals `[a,b]` of G's edges `{a,b}` cover every
path cut1,...,n-1; P closes iff G contains `{1,n-1}` (n>=3). Converse
together with path reversal becomes order reflection of G. A marked path
is blue exactly when G is fixed by that particular reflection.
See [the full proof, gauges and eight-row table](graph_audit.md).

The three simple graphs on four vertices with three edges are indeed
`P4`, `K1,3`, and `K3 + K1`, with labelled counts12,4,4. Graph complement
fixes the path type and exchanges the latter two. For a marked four-vertex
Hamiltonian path, its free-arc carrier `K4 minus P4` is itself P4,
specifically edges13,14,24. This explains exactly where the line appears.
The three graph types, the three arithmetic prime shapes, and the three
converse-merged tournament classes are distinct quotients; equality of
their counts does not provide an incidence-preserving identification.

## The two decompositions of ten occupy the same plane

The recovered `F=U+S` classification is `p,p^3,p^2qr`, where p,q,r in
the last term are distinct primes. F counts proper nontrivial divisors,
S the squarefree ones and U the prime ones. Since primes already belong
to S, this is a count equality. A tagged copy of primes maps to the
nonsquarefree divisors by `ell -> p*lcm(p,ell)`. Literal `p^q*r` with
only two distinct prime bases is a different family and does not qualify.

On the ten proper divisors of p^2qr use the faithful coordinates

    d=p^a q^b r^c  ->  (f,g)=(a-b-c,b-c).

Here is the entire object, with no arbitrary matching of equal-sized sets:

| g \ f | -2 | -1 | 0 | 1 | 2 | Row size |
|---|---|---|---|---|---|---:|
| +1 | | q | pq | p^2q | | 3 |
| 0 | qr | pqr | | p | p^2 | 4 |
| -1 | | r | pr | p^2r | | 3 |

Rows give **3-4-3**. The negative, zero and positive f regions give
**4-2-4**. Divisor complement flips both coordinates; exchange of q,r
flips g; their composite flips f. The central two points are exactly
the fixed points of the last operation. This is a precise reconstruction
of the requested4-2-4, since no earlier definition of that phrase was
recovered.

For p^kqr, the faithful chart is `(2a-k(b+c),b-c)`: only the excluded
bottom and top divisors collide. The sign split becomes `2k,2,2k` for
even k and `2k+1,0,2k+1` for odd k. A tetrahedral ten-point model cannot
make the full divisor symmetry geometric: complement has zero fixed
divisors, while every tetrahedral involution fixes at least two pair
letters. [Proof, inverse map and controls](arithmetic.md).

## The square arithmetic counts different deletions

For a 25-element set X,

    |X^2 minus diagonal|=600,
    |(X minus {x0})^2|=576.

The latter retains24 diagonal cells and552 off-diagonal cells. These
are literal Hamiltonian-path counts in an explicit ordinal product of
the strong four-tournament with five paths. The correct sum
`625+125+25+1=776` is realized by words in those five paths of lengths
4,3,2,0; excluding length1 is a supplied choice, not a consequence of
diagonal removal. [Explicit maps and independent proof audit](square_counts.md).

The inherited balanced two-return family has a separate actual deletion
map with fibres4,2,1: forgetting one Boolean bit in each coordinate gives
two lifts for finite states and one for the exceptional state. Thus
`q'=2q-1`, `H'=4H-4q+1`. This boundary state explains the correction;
the balanced-family q is `2^(k-1)+1`, so25 cannot be silently substituted
as one of its q values.

## Reset debt is stable after the gauge is fixed

For a lawful Collatz word write the shifted endpoint as

    z' = a z + P - D,       a>0, P>=0, D>=0,

where D is the discarded reset contribution weighted by later slopes.
For u followed by v,

    a_uv=a_v a_u,       D_uv=a_v D_u+D_v.

Therefore `D<=a-1` is composition-stable, but permits only expansion.
In reciprocal source coordinates `lambda=1/a, delta=D/a`, the law is

    lambda_uv=lambda_u lambda_v,
    delta_uv=delta_u+lambda_u delta_v,

and `delta<=lambda-1` is also stable; it now means `D<=1-a` and permits
only contraction. The two conditions cannot be mixed: native word13
from19 to11 concatenates one factor from each cone and lies in neither.

Legal domains compose by `C_uv=C_u intersection F_u^(-1)(C_v)`;
the original source and all first-exit inequalities must be retained.
Actual positive-charge payment is `D<n-F_w(n)`, a stricter demand than
descent alone. Paid authenticated descending blocks compose, but a first
exit cannot be split into earlier descending blocks. Universal entry
remains open. [Exact cylinders, proofs and short hostiles](reset_composition.md).

Concurrent mainline commit70cb917742 supplies a useful continuation of
the same forward cone: if `D_u<=c(a_u-1)`, `0<=c<=1`, and an actual final
odd step makes the endpoint smaller than the original source, then total
debt is less than2 and positive charge is paid. The generalized terminal
packet condition is `c(1-a_v)+D_v<=2`. Our note integrates and independently
checks that argument, including the c=0 boundary requiring a separate
dyadic-debt argument; incoming8dcc679ae3 independently made the same repair.
This certifies
authenticated descents; selection of such prefixes and terminal steps for
every source remains open.

## Portfolio, connection audit and verification

The inheritance pass used the tiling/free-action mechanism, the two
different divisor3-4-3 partitions, the balanced square/deletion map and
the positive-charge reset decomposition. Hostiles were missing path
order, complement exchanging a corner and midpoint, variable deletion
fibres, and the authenticated later descent55->53 with unpaid charge.
The less-used data were ordered cut intervals, the full divisor symmetry
action, the exceptional pointed state, and exact source pullbacks.

The portfolio was **Anchor:** marked trijection and graph invariants;
**Niche:** guarded reset composition; **Wildcard:** the two planar views
of ten and the square arithmetic. The five-concept board was revisited
after the recovered sources and exact checks:

| Concept | What the new object changes | What must remain attached |
|---|---|---|
| Marked path coordinates | Turns strength into interval coverage | Path order; arbitrary graph isomorphism loses it |
| Divisor exponent box | Two reflections give both requested splits | Prime labels and the full planar coordinates |
| Square path product | Separates diagonal deletion from point deletion | Whole paths; equal endpoints do not force blue |
| Exceptional deletion state | Explains fibres4,2,1 and growth corrections | Infinity/finite type in both coordinates |
| Affine reset composition | Gives exact stable debt inequalities | Multiplier convention, word order, source domain |

The shared research move is to transport an operation with its retained
coordinates. There is no asserted direct divisor-to-Collatz or
graph-class-to-prime-shape map. Existing META-PATTERNS cards used were
“Search the statement before the method”, “Type every analogy and every
implication”, and “Retain local incidence through counting”; no new
method card is needed for this repackaging.

[VALIDATION.md](VALIDATION.md) records reproduction, proof reviews,
bounded universes, corrected draft error and saved receipts. The all-order
proofs are structural; no large enumeration was used.
