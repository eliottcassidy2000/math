# What a hidden bit can preserve: exact ports and a projective path chart

2026-11-03. **PROVED** elementary propositions below, using the inherited
rooted-ear identity; **FINITE-EXACT** receipts in `verification.json`,
`audit.json`, `projective.json`. The independent computational path was
written and run by the same agent; this is not an external referee audit.
No external novelty claim. Universal odd-spectrum completeness remains **OPEN**.

## P1. A general criterion, relative to specified future operations

Let X carry a scalar observation pi:X→N and named partial operations.
Declare x~y when every finite operation word has the same legality and,
when legal, the same final scalar observation from x and y. Include the
empty word. Then an exact predictive representation `(pi, b extra bits)`
exists precisely when each pi-fibre has at most 2^b equivalence classes.

Proof: a faithful representation cannot merge distinguishable futures.
Conversely label the classes separately inside each scalar fibre. Future
equivalence is a congruence: appending any further word shows that each
legal operation has a well-defined output class. These labels give the
required state transition. This is a set-theoretic statement; it does not
supply computable class labels or an efficient decoder. Nor does it require
recovering differences that no permitted future operation can observe.

This is the standard observation-equivalence principle, not a new claim
about automata theory. The new application below supplies an unbounded,
strong-tournament witness using just one actual attachment.

## P2. Every port is distinguishable in a small strong family

Let E_n, n>=3, have vertices 0,...,n−1, all arcs increasing except n−1→0.
Its directed spanning cycle proves strongness. The inherited ordered path
decoder in [the endpoint/Pythagorean note](../endpoint_pythagorean_fibers_20261019.md)
gives h=H(E_n)=2^(n−2)+1.

Let K(X,q) count ordered spanning pairs (P,Q) of disjoint directed paths,
P nonempty ending at q, Q possibly empty. The inherited
[universal guard lift](../session164_emergent_20261103/complementary_guard/UNIVERSAL_GUARD_LIFT.md)
adds a,b with q→a→b→q and a,b→X−q. It preserves all old vertices and arcs and
has exact raw count h+2K(X,q)=3h+2J(X,q), where J=K−h>=0.

For E_n the complete rooted profile is

    K_0 = 2^(n−2)+1,
    K_q = 2^q + 3^(q−1) 2^(n−1−q) − [q=n−2],  1<=q<=n−2,
    K_(n−1) = 2^(n−1)−1.                              (1)

Proof. For an interior q, first count pairs whose two paths both increase.
The first can use any subset of the q earlier vertices, hence 2^q choices.
Exactly one choice fails when q=n−2: the second path would be (0,n−1),
whose required arc is absent. Next the backward edge belongs either to P
or to Q. In either case, each interior vertex below q has three allocations
(one path, or either increasing run of the path using the backward edge),
and each interior vertex above q has two. Each case gives
3^(q−1)2^(n−2−q). A pair cannot use the unique backward edge twice, so
these disjoint counts prove the middle formula. At q=0 there is the pair
beginning with the singleton0 and 2^(n−2) pairs using the backward edge in
P. At q=n−1 the first path must increase: all 2^(n−1) subsets are valid
except the first path (0,n−1). These prove both endpoints.

For every n>=4, all n numbers in (1) are distinct. Indeed K_1=K_0+1;
interior values increase, since their successive difference is

    2^q + 3^(q−1)2^(n−2−q) − [q+1=n−2] > 0.

At n=4 the list is (5,6,9,7); at n=5 it is (9,10,16,25,15).
For n>=6, K_2 < K_(n−1) < K_3: the two positive differences are
2^(n−3)−5 and 2^(n−4)+9. Thus the last port cannot collide either.
At n=3 all three values equal3, the sharp cyclic-symmetry exception.

Consequences, with explicit scope:

1. A representation of all marked E_n that retains even the exact base
   tournament and order, but stores its mark using a tag, needs at least
   ceil(log2 n) bits to predict this one ear operation. A uniform one- or
   two-bit repair therefore fails on this family. An integer port q is an
   efficient repair for this one step; closure under repeated attachments
   is a further problem.
2. E_n is asymmetric for n>=4: automorphisms preserve K and hence fix every
   vertex. This is an elementary alternative proof, not a novelty claim
   about this familiar tournament.
3. Substituting a canonical port without retaining the original port is
   observable immediately in the raw path count. H/Aut presentation counts
   play no role in any statement here.

## P3. Actual scalar roots: a better65 witness and a memory lower bound

On the six printed witnesses of
[the ternary forest](../triple_reentry_interfaces_20261103/RETURN_FOREST.md),
the ear responses distinguish respectively 1,6,6,8,7,7 marks. Thus each
nontrivial printed root already requires at least three bits for this
operation. This is a direct consequence of the inherited J profiles;
`audit.py` recomputes them by permutations and cuts, independently of subset DP.

There is a useful source choice at65. The old seven-vertex witness has

    J=(98,113,137,118,41,5,4),
    ear H=(391,421,469,431,277,205,203).

The standard eight-vertex E8, also with H=65, has

    J=(0,1,35,87,167,291,484,62),
    ear H=(195,197,265,369,529,777,1163,319).

These fifteen responses are distinct. If the retained scalar is only H=65,
four bits are necessary to distinguish this particular bank of fifteen
marked sources. If order and the appropriate source chart are also retained,
three bits suffice to identify the mark within either listed source.
This is a sharp bound for this one-step bank, not for all H=65 sources or
all future operations.

The second witness supplies actual first-ear outputs 3*65 and 3*65+2.
It repairs the earlier suggestion that a better65 witness still needed to
be found. It is NOT an induced repair of the original seven-vertex source:
that source contains R5 and has the inherited birth-only readiness
obstruction. Scalar-root selection permits choosing E8 from the outset;
preserving a supplied source does not permit that replacement.
There is no J=2 port in E8, so the same construction does not give199.
Neither195 nor197 comes with a proved repeated ternary return in this note.

There is nevertheless an exact small positive control: restrict E_n to the
two canonical ports0,1. Inserting a new interior vertex just before the
last vertex grows E_n to E_(n+1), preserving every old arc and both ports.
The scalar law is h→2h−1 and their slacks remain0,1. One bit identifies the
chosen port in this restricted, insertion-closed family. This does not
contradict P2, which permits every mark.

The first ear leaves this family for all r=n−2>=3: at port0 its output is
3*2^r+3, at port1 it is3*2^r+5. Subtracting1 gives respectively
2(3*2^(r−1)+1) and4(3*2^(r−2)+1); the parenthesized factors are odd and>1.
Neither is a power of2, so neither output equals H(E_m) for any m.
This is an explicit return obstruction for the attractive one-bit chart,
not a failure of its native insertion rule or of the larger proof target.

## P4. An exact 63/65 bridge, and the boundary of its symmetry

Put r=n−2, q=2^r, choose an F2-basis of F_q, and identify the subset of
interior labels before the backward edge with its binary field coordinate.
The exceptional increasing Hamilton path is infinity. The inherited decoder
is therefore an explicit bijection

    Hamilton paths of E_(r+2) ↔ P^1(F_q)=F_q union {infinity}.       (2)

Transporting PGL_2(F_q) through this decoder permutes actual Hamilton paths.
Subset complement is translation by the nonzero all-ones basis vector,
fixing infinity. Multiplication by a generator of F_q^* fixes0,infinity
and cycles the other q−1 addresses.

For the other cycle, view F_(q²) as a two-dimensional F_q-vector space.
Its norm-one subgroup has q+1 elements: the norm is onto because its
restriction to F_q^* is squaring, a bijection in characteristic two.
Multiplication acts on projective lines. A norm-one element stabilizing
a line must lie in F_q, where its norm is its square; hence it is1.
The action is free, and thus regular on the q+1 lines. Since a finite
field's multiplicative group is cyclic, this gives a (q+1)-cycle.
These are the usual split and nonsplit projective actions; the contribution
here is their concrete transport to the named Hamilton-path decoder.

At r=6, equation (2) gives65 paths, a63-cycle on nonzero finite addresses,
and a65-cycle from norm-one units of F4096. `projective.py` constructs these
using F64=F2[X]/(X^6+X+1), checks the quadratic extension, and decodes every
member of the65-cycle into a literal directed path in E8.

This is a genuine map with chosen coordinates, but it is a path-set action,
not a tournament vertex action. P2 proves E8 has no nontrivial vertex
automorphism. In particular this large artificial symmetry cannot be used
to move its eight rooted ear contexts while preserving their counts.
The original six scalar roots came from deleted nodes of a ternary tree;
this projective construction explains the pair63/65 in a different typed
object, not the entire root list or the missing reentry theorem.

## Reproduction

From the repository root:

```text
python3 05-knowledge/results/six_roots_hidden_bits_20261103/verify.py
python3 05-knowledge/results/six_roots_hidden_bits_20261103/audit.py
python3 05-knowledge/results/six_roots_hidden_bits_20261103/projective.py
```

Controls: E_n profile through14; literal ear outputs through n=9; independent
permutation/cut enumeration through8 and on all six roots; separately
memoized path counts on each seed's literal ear. Normal and optimized
Python runs retain explicit failure checks. The all-size proof is (1),
not extrapolation of these computations.
