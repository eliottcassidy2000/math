# A completely enumerated infinite source family and its internal moves

**PROVED in this note, using the reviewed all-k filtration theorem.** The
family consists exactly of tournaments whose strong components are K1 or
C3. The classification below covers their ordinary classes, automorphisms,
all within-family one-vertex births/deletions and single-arc flips, and
their arbitrary-extension common divisors. It does not classify every
outside extension count, native path-complement edge, or tournament.

## 1. Complete source invariant and all-height external arithmetic

An ordinary class is uniquely an ordered word in the components1 and3.
Let n be its vertex order, r its number of3s, and nu the sum of floor(L/2)
over maximal consecutive3-runs of lengths L. Then

```
H=3^r,          |Aut(T)|=3^r,
C_k(T)=3^max(r-k-min(k,nu),0),       L(T)=r-nu.          (1)
```

Here C_k is the gcd of all k-vertex induced extension counts, and L is
the least k with C_k=1. The first two equalities follow from ordered SCC
multiplication and independent triangle rotations. The last two are the
new [all-k half-contact theorem](../arithmetic/TRIANGLE_ALL_K_PROOF.md),
with independent [audit](../berggren/ALL_K_AUDIT.md). For r0 the source is
transitive, C_k=1 and L=0, including the empty-word convention.

Thus r and nu determine this entire external arithmetic profile; the
full SCC word remains necessary for actual source composition. A gcd1
profile does not mean that every desired extension count is attained.

## 2. Joint ordinary and labelled enumeration

Let z mark vertices, u triangles and v the matching size nu. The exact
ordinary generating function is

```
F(z,u,v)=(1+u*z^3)/(1-z-u*z^4-v*u^2*z^6).              (2)
```

Proof. With x=u z^3, a nonempty run of triangles has generating series
`R=(x+v*x^2)/(1-v*x^2)`: its matching number is floor(runlength/2).
Such runs are separated by nonempty singleton words; arbitrary singleton
padding is allowed at the two ends. Including the all-singleton class,

```
F=1/(1-z)+R/((1-z)(1-z-zR))=(1+R)/(1-z-zR),
```

which simplifies to (2). Setting v1 gives the familiar ordered-component
enumerator `1/(1-z-u*z^3)`. No different ordinary classes are identified
under converse.

Since the automorphism group has order3^r, a class contributes n!/3^r
labelled tournaments. Thus the corresponding labelled exponential
generating function is F(z,u/3,v). This is an exact labelled mass formula
for this subfamily of A000568, not a new total formula for all tournaments.

## 3. Exact lifetime layers and their size growth

Substitute u=w,v=w^(-1) in (2). There are no negative resulting powers,
because r-nu>=0. The full lifetime enumerator is

```
G(z,w)=(1+w*z^3)/(1-z-w*z^4-w*z^6).                    (3)
```

For lifetime ell>=1 the coefficient of w^ell is exactly

```
z^(4ell-1) (1+z^2)^(ell-1) (1+z^3)/(1-z)^(ell+1).     (4)
```

Lifetime0 has generating function1/(1-z). Formula (4) gives an exact
integer coefficient formula, eventually a polynomial in n of degree ell
with leading coefficient `2^ell/ell!`. It applies polynomially for all
n>=6ell; the coefficient form remains valid below that boundary too.
Its least order is4ell-1, attained uniquely by ell triangles separated
by singletons. Therefore the maximal lifetime within this family at
order n is exactly floor((n+1)/4). This is not asserted to be the maximum
over all tournaments.

To refine any fixed C_k layer, sum coefficients of (2) over pairs(r,nu)
with the desired exponent `max(r-k-min(k,nu),0)`. Thus the entire family
has explicit size enumerators for its arithmetic transition obstructions.

## 4. Every one-vertex birth that stays in the family

There are precisely two rules on the source SCC word:

1. Insert1 at any component boundary. This keeps H unchanged.
2. Replace two consecutive1s by3, using the newborn to close their
   directed triangle. This multiplies H by3.

Each component boundary gives one actual old-new cut. Each adjacent11
pair gives exactly one additional cyclic cut. Different positions can
produce the same ordinary target, and their multiplicities are retained.
If the source has c components and a adjacent11 pairs, exactly c+1+a
of its2^n one-vertex cuts stay in this family.

Completeness follows from SCC interval absorption: the newborn's strong
component either is a singleton inserted between old components, or has
size3 and absorbs exactly two adjacent old singleton components. It
cannot absorb an old triangle and still have size3. No interval containing
other components can be skipped. Deletion reverses these rules: remove
one singleton, or replace3 by11 (three labelled vertex choices).

A neutral singleton birth splitting a triangle run into lengths a,b
lowers nu by1 exactly when both a,b are odd. It then raises lifetime by1;
otherwise the complete C_k profile is unchanged. In the odd/odd case,
the exponent of C_k rises by1 precisely for
`old_nu <= k <= r-old_nu`. These are exact all-size transition rules,
not an extrapolation from the first H9 example.

## 5. Every single-arc flip that stays in the family

The rules, with labelled arc multiplicities, are exactly

```
3 -> 111:             three internal triangle arcs;
111 -> 3:             the arc between the first and third singletons;
adjacent11 -> 11:     one arc, reversing the order of those two names.
```

For the first rule every reversed triangle arc makes TT3. For the second,
the single reversed outer arc makes C3. The last rule is an ordinary-class
loop, but changes labels. Any other cross-component reversed arc merges
an interval with at least four vertices into a strong component and exits
the family. This supplies a complete internal single-arc metagraph, with
birth/deletion rules above and explicit numerical invariants at every node.

The [gap-transport macros](GAP_TRANSPORT.md) supply additional paths that
temporarily leave the family through S4, or stay in it through the
3->111->3 route. They preserve the final H and every old name while
changing the matching invariant. Their intermediate count costs are
explicit, so composition does not silently treat a macro as one edge.

## Reproduction and scope

`triangle_enumeration.py` independently compares (2) with every component
word through vertex order24, verifies (4), and tabulates ordinary and
labelled counts and all lifetime layers. It also constructs every source
word through order9 and checks ALL14,791 old-new cut masks and1,392 old
arc flips against the complete rules, including multiplicities. SCCs are
computed from literal reachability, independently of the word grammar.
Both normal and optimized runs pass; output is `triangle_enumeration.json`.

This is a fully controlled infinite type of metagraph source with exact
composition rules for the named operations. Its numerical source spectrum
remains {3^r}; it supplies structural and enumeration tools, not a proof
that all odd integers other than7 and21 occur as Hamilton-path counts.
