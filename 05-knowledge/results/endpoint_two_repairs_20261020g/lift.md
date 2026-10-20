# Port-aware lifts and the last-port block repair

**PROVED:** every closed standard-blue path with endpoint scores `(n-3,2)`,
last-port exception `c->a`, and chord `c->x` has an explicit three-black
detour when n>=9. The proof retains a whole path block rather than assuming
that a singleton ledger commutes with substitution. A complementary result
repairs every aligned blue line of `C3[TT2,V,TT2]`: three literal black moves
when |V|>=2, and two black moves modulo ordinary isomorphism when |V|=1.

**PROVED:** an exact two-port parity criterion governs when a K ledger lifts
through path-block substitution. A seven-vertex hostile shows why the
unguarded lift fails even in the desired endpoint layer. **FINITE-EXACT:**
the small controls below. General black-core connectivity is not inferred.

## 1. Inheritance and retained coordinates

The closest proved mechanisms are the [two-block endpoint repair](../strong_blue_games_20261020f/detours.md)
and the [cyclic singleton-extension ledgers](../kernels_games_detours_20261020f/detours.md).
The first retains actual endpoint-preserving paths inside reversed blocks;
the second retains the middle and its two named exterior vertices. A middle
module alone is not the same data as the ports used by a sequence of moves.

Write R_P(T)=T^P for reversal of the path arcs and
`K_P(T)=(T^P)^op` for nonpath complementation. A sequence of R words
`Q,U,V` is the K sequence `Q,reverse(U),V`. A word is blue precisely when
its positional reversal is anti. A blackness argument below either uses
this full predicate or a sufficient nonzero endpoint-score excess.

The concept board is: exceptional endpoint incidences; two-port parity;
block path repair; positional anti-involutions; and small-order boundaries.
Anchor: the `(n-3,2)` endpoint layer. Niche: a correct substitution lemma.
Wildcard: a mirrored blow-up that initially defeats the obvious lift.

## 2. A whole-block repair of the last-port case

Let T have the blue directed Hamiltonian path

    P=(a,b,x,D,c,z),                 |D|=n-5,

where n>=7. Assume a has outdegree n-3, its two inneighbours are c,z, and
`c->x`. Blue reversal then makes the outneighbours of z exactly a,b.
In particular P closes by z->a. No hypothesis on the other interior arcs
is imposed. Write d,e for the first two vertices of D.

**Theorem A (PROVED).** At n=7 and at every n>=9 the following are three black K moves
ending exactly at K_P(T):

    Q1=(a,x,D,c,z,b),
    Q2=(D,a,x,c,z,b),
    Q3=(D,a,b,x,c,z).                                  (1)

**Direction and target.** Q1 uses the original path except for a->x and
z->b, both supplied by the endpoint incidences. Let S=K_Q1(T).
The new Q2 joins last(D)->a and x->c are complements of the nonpath
arcs a->last(D) and c->x. All its other joins are retained. After Q2,
the two original arcs a->b and b->x have each been complemented twice,
and the other Q3 joins are retained. Thus all words are directed.

Each internal D edge occurs three times. On the boundary the pairs ax,
last(D)a, xc and zb occur twice; the odd remainder is precisely the
boundary of P. Consequently `E(Q1) xor E(Q2) xor E(Q3)=E(P)`. Three K
moves also complement every pair an odd number of times, proving the
literal labelled target.

**Q1 is black for n>=7.** Its positional reversal is the same as that of
the legal cyclic cut `(b,x,D,c,z,a)`. Let rho reverse P and pi be its
cyclic successor. The cut is blue iff `pi rho pi^-1 rho=pi^2` is an
automorphism. For odd n this acts transitively, so all scores are
`(n-1)/2`. Comparing d(a)=n-3 forces n=5. For even n=2m the two parity
orbits each have size m. Within a's orbit its score is `(m-1)/2` by
transitivity. It beats at most m-1 vertices of the other orbit because
z belongs there and beats a. Hence

    2m-3 = d(a) <= (m-1)/2+(m-1),

which forces m<=3. Neither possibility occurs at n>=7.

**Q2 is black except possibly at orders 7 and 9.** In S,

    N^-_S(z)={a,c},       e->a.

The first equality follows by complementing z's original incidences,
retaining the Q1 joins c->z and z->b. The second follows because a->e
was a nonpath arc of Q1. The reversal of Q2 exchanges e,z, and sends a
to the vertex in zero-based position 4. If it were anti, e->a would
force `r_Q2(a)->z`, so that position would have to contain a or c.
Since D has n-5 vertices, this happens only for n=7 or n=9.

**Q3 is black except possibly at orders 8 and 9.** Check at its other
endpoint B=K_P(T), where a K word has the same colour. Here

    N^-_B(z)={a,b,c},       d->a.

Q3 pairs d with z and again sends a to position 4. Antisymmetry would
put that position in {a,b,c}, possible only at n=8 or n=9.

**The n=9 boundary is also black.** In both Q2 and Q3, a is at the fixed
middle position4. But d_S(a)=d_B(a)=3, whereas an anti-involution's fixed
vertex in an order-nine tournament must have score4. Thus neither word
can be blue. This fixed-point refinement was supplied and independently
audited by the companion theory lane. It completes the all-order assertion
n>=9. The remaining small possibilities are allowances of these witnesses,
not claims that the words are blue there.

**The n=7 boundary is black by a local arc.** Here D=(d,e), and Q2 fixes
x while exchanging b,d. In S the arc b->x is false (Q1 complemented
that omitted P arc), whereas x->d is true (Q1 retained it). These violate
the required anti equation. Q3 is already black by its position4 witness.
This additional refinement also comes from the companion theory audit.
Only n=8 remains an exception to this particular universal word. QED.

**The block move behind the proof.** More generally, write
`P=(a,b,x,D,U,z)`, with nonempty D,U and last(U)=c. If first(U)->x,
the words

    (a,x,D,U,z,b), (D,a,x,U,z,b), (D,a,b,x,U,z)

are directed and have the same exact XOR target. This extends the
singleton-D ledger developed in the companion lane. Choosing U=(c)
makes its guard automatic in Theorem A. The correction is to retain the
whole block and recalculate its positional reversal; a singleton reversal
identity is not silently reused for a longer block.

## 3. Mirrored two-vertex expansions have an explicit repair

Let V have r>=1 vertices and a blue path `M=(b,...,c)`. Set

    A=(a0,a1), Z=(z0,z1),
    a0->a1, z0->z1, A->V->Z->A,
    P=(a0,a1,M,z0,z1).                                  (2)

The source is `C3[TT2,V,TT2]`, and P is blue. Its endpoint scores are
`d(a0)=r+1=n-3` and `d(z1)=2`. Both endpoints are strong: P is a cycle
in the source, while K_P(T) has the cycle `(a0,a1,z0,z1,M,a0)`.

For r>=2 define

    H_R=(b,z0,reverse(M without b),z1),
    Q=(a0,M,z0,z1,a1),
    U=(a0,a1,H_R),
    Vword=(a1,reverse(H_R),a0).                          (3)

Here R denotes the displayed block `(M,z0,z1)`. After reversal of its
path, H_R is an endpoint-preserving Hamiltonian path: b->z0 is an
unchanged nonpath cross-arc, z0->c is reversed, the remainder of M is
traversed backwards, and its last vertex goes to z1 by an unchanged
cross-arc. **No openness assumption on M is needed.** The insertion of
z0 replaces the potentially absent b->c chord.

The inherited two-block identity now makes (3) an exact three-R-move
route to T^P. The three endpoint-score sums are

    2r+1, r+2, 2r,                 n-1=r+3.             (4)

For r>=3 the first two moves are black. The third is always black,
including its balanced r=3 boundary: just before it, z1->a0 is unchanged
but a1->b has become b->a1. Its reversal exchanges z1,b and a0,a1, so
these two arcs contradict antisymmetry. Thus every r>=3 has the literal
three-black detour (3), equivalently the K words Q,reverse(U),Vword.

At r=2, in the numbering P=(0,1,2,3,4,5), a separate literal three-R
certificate is

    (0,2,5,1,3,4), (4,0,1,5,2,3), (3,1,2,0,4,5).

At r=1, use the black R word `(1,2,3,4,0)` from T and the black R word
`(0,2,4,3,1)` from T^P. Their neighbours are isomorphic under the map
`(4,0,2,3,1)`, a two-black quotient detour. The latter two boundaries are
FINITE-EXACT fixed certificates. The general theorem is the block proof,
not an extrapolation from those checks.

The root lane's broader last-port `x->c` theorem subsumes the r>=3
geometric pattern (2). Its proof was independently audited here: the
same two-block path works and the common final-word arc witness avoids
the balanced score boundary. This section retains the simpler family
because it exposes the following failed lift particularly clearly.

## 4. Exact two-port criterion for lifting a K ledger

Replace each vertex v of a tournament S by a nonempty tournament B_v with
a fixed directed Hamiltonian word H_v. Initially all cross-arcs follow S;
within each block the first and last vertices of H_v are its entry and
exit ports. Lift a word by replacing every occurrence of v with H_v.

Let Q1,...,Q_l be a directed K ledger on S, l odd, ending exactly at
K_P(S). For each ordered pair u,v let C_uv(t) count the words among
Q1,...,Qt in which u immediately precedes v. Only pairs with at least
one nonsingleton block require the additional test below.

**Theorem L (PROVED).**

1. The lifted word Qi is directed in the formal lifted prefix state iff
   `C_vu(i-1)` is even for every adjacent ordered pair u,v of Qi touching
   a nonsingleton block.
2. The formal lifted ledger has the literal target K_lift(P)(lift(S)) iff
   `C_uv(l) mod2 = [u immediately precedes v in P]` for every ordered
   pair touching a nonsingleton block.

The first condition at every step makes the formal ledger a genuine
sequence of directed Hamiltonian moves. The theorem does not by itself
certify that those moves are black.

**Proof.** For a fixed block pair the two opposite macro directions use
two distinct physical pairs: exit(u)--entry(v), or entry(u)--exit(v).
They coincide exactly when both blocks are singletons. Before step i,
the macro orientation has been complemented i-1 times and toggled
`C_uv+C_vu` times. The physical pair selected by u followed by v has
been complemented i-1 times and toggled only C_uv times. Its difference
from the known legal macro orientation is precisely C_vu modulo two.
All internal H_v arcs remain directed because every K word retains them.
This proves (1). After an odd number of moves, the nonselected physical
cross-pairs and internal block pairs have the same global complement
phase as the requested target. Matching the two selected physical pairs
separately gives (2). Singleton pairs already agree by the given macro
endpoint identity. QED.

Thus the retained sidecar is **oriented port incidence modulo two at
every prefix**, not just the macro edge XOR at the final endpoint.

## 5. A seven-vertex hostile to the naive mirrored lift

Start with the cyclic singleton extension of TT3, on standard path
01234. Its inherited three-black K ledger is

    12340, 42310, 04231.

Replace its first and last vertices by TT2, with the other blocks single.
This is exactly (2) with r=3. Naive block expansion gives

    (2,3,4,5,6,0,1),
    (5,6,3,4,2,0,1),
    (0,1,5,6,3,4,2).

The first two words are directed, but the third needs 1->5 while the
actual current arc is 5->1. The macro arc from the last endpoint to the
first was used in the first word; its opposite occurrence in the third
word now selects different ports. Its final XOR also fails to match the
requested labelled target. This is not a disconnection: (3) repairs this
same source by three genuine black moves.

## 6. Reproduction and finite boundaries

`lift_check.py` imports no producer modules. Normal and optimized Python
check the following:

- Every blue middle presentation of orders 3..6 (86 cases), seven larger
  deterministic mirrored-chord middles of orders up to31, and both small
  family certificates. Every actual R/K word and labelled target is checked.
- Thirty-six module expansions check both directions of Theorem L's
  prefix and endpoint predicates, including the seven-vertex hostile.
- Complete last-port `c->x` universes at n=7,8,9 contain 8,32,256
  presentations. The fixed-tail ledger works in 8,30,256 respectively.
  Trying every legal D/U cut adds one n=8 presentation. Its sole remaining
  word failure has rows `(62,28,184,240,224,194,135,3)` and is passed to
  the independent finite-repair lane. A failed word is not a bridge.
- Sixteen deterministic chord samples at each n=10,11,16,25 check the
  proved large-order A construction; they are not a tournament census.

    python 05-knowledge/results/endpoint_two_repairs_20261020g/lift_check.py
    python -O 05-knowledge/results/endpoint_two_repairs_20261020g/lift_check.py

The exact results and source hash are in `lift_check.json`. This packet
does not claim that the other endpoint-layer branches or arbitrary strong
mixed blue edges are all repaired. Their status belongs to the separate
case-partition synthesis, after all branch proofs have been audited.
