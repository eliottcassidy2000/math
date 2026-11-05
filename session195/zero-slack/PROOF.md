# Zero merge slack, root outdegree, and an unavoidable unbounded coordinate

Session 195, 2026-11-05. **PROVED all-order structural classification and grammar-specific sparsity obstruction.** The enlarged 67-gadget grammar's all-order coverage, and the full `{7,21}` Hamiltonian-path spectrum, remain **OPEN**. No finite coverage statement is used to prove either theorem here.

## 1. Conventions and inheritance

For a nonempty rooted tournament `(T,r)`, let `h=H(T)`, `d=H(T-r)`, and let e count Hamiltonian paths ending at r. Put `H(empty)=1`. Let R count ordered spanning two-path covers `(P,Q)` in which P may be empty, Q is nonempty, and Q starts at r. Define the nonnegative merge slack

    f = R+e-h-d.

The four-coordinate grammar and its exact transitions are proved in the companion combined-growth packet. This note identifies f=0 by an elementary graph property, then separates a universal consequence from a finite-bank grammar obstruction. The closest prior mechanisms are the Session 194 endpoint floor and module threshold, and Session 195's complete two-path boundary matrix. No literature novelty is claimed.

For a tournament S and distinct x,y, let `B_yx(S)` count ordered nonempty two-path covers `(P,Q)` with P ending at y and Q starting at x. These paths together cover S. This B is a boundary matrix entry, not the unrooted two-bin scalar used in a gadget signature below.

## 2. A positive two-endpoint extension lemma

**Lemma 1.** For any tournament S and any distinct vertices x,y,

    B_xy(S)+B_yx(S) >= H(S-{x,y}) >= 1.               (1)

**Proof.** Fix a Hamiltonian word `w_1,...,w_m` of the other vertices. If m=0, use `P=(x),Q=(y)`.

Otherwise, find the first index i such that at least one of x,y beats w_i. If such i exists, choose a in {x,y} beating w_i by a fixed tie rule, and let b be the other endpoint. Every preceding w_j beats both x and y, by minimality of i. Thus

    P=(w_1,...,w_(i-1),b),   Q=(a,w_i,...,w_m)

are legal nonempty paths, ending and starting at the two distinct named endpoints. If no such i exists, every w_j beats both endpoints; use `P=(w_1,...,w_m,x), Q=(y)`.

Deleting x,y from the two paths and concatenating the remainders recovers the original Hamiltonian word. The construction is therefore injective into the disjoint union counted on the left of (1). The usual insertion argument gives a Hamiltonian word on every induced tournament, including the unique empty word. QED.

This proof needs neither strongness nor a Hamiltonian path containing the arc between x and y. That stronger arc-containing assertion would be false in a transitive tournament with a nonconsecutive arc.

## 3. Universal exact classification of zero slack

**Theorem 2.** Let `S=T-r` and `X=N_T^+(r)`. Then

    f = sum_(x,y in X, x!=y) B_yx(S).                (2)

Consequently

    f=0  if and only if  outdegree_T(r)<=1.          (3)

**Proof.** Partition the R-covers `(P,Q)` according to whether they can be concatenated. The mergeable ones are exactly those with P empty or with its final vertex beating r. Concatenation is a bijection to Hamiltonian paths of T, inverted by splitting immediately before r; they contribute h.

In every remaining cover, P is nonempty and its last vertex y is beaten by r. Among these, the covers with `Q=(r)` contribute `d-e`: their first path is a Hamiltonian word of S whose final vertex does not beat r.

In each remaining cover Q has at least two vertices. Delete its first vertex r; the resulting second path begins at a vertex x beaten by r, and the first path ends at a different out-neighbor y. This is exactly an ordered nonempty two-path cover counted by `B_yx(S)`. Reinsert r before x to invert the correspondence. Thus `R=h+(d-e)+sum B_yx`, proving (2).

If |X|<=1 the sum is empty. If |X|>=2, choose two of its vertices; their two boundary terms have positive sum by Lemma 1. This proves (3). QED.

This theorem applies to **every** physical rooted tournament, not just a specified growth grammar. In particular, the f=0 image of the singleton root-switch A is explained by its actual new root having precisely one out-neighbor. A universal sink has `(h,e,R,d)=(d,d,d,d)` and also has zero slack.

## 4. The complete unit-deletion, zero-slack spectrum

**Theorem 3.** Among all rooted tournaments satisfying `d=1` and `f=0`, the possible Hamiltonian-path counts are exactly

    {1} union {2^m+1 : m>=1}.                         (4)

Every nontrivial displayed value has a strong representative.

**Proof.** A tournament has H=1 exactly when it is transitive. A transitive tournament has its unique increasing Hamiltonian word; a nontransitive tournament contains a directed triangle, whose three paths inject into the whole tournament by ordinary vertex insertion.

Thus `T-r` is a transitive chain, of order n say. By Theorem 2 the root r has no out-neighbor or exactly one. If it has none, it is a universal sink and H(T)=1.

Otherwise let its sole out-neighbor be the i-th vertex a in the chain. The exact singleton-birth path split gives

    H(T)=H(T-r)+R_a(T-r)-e_a(T-r)
        =1+2^(n-i)-[i=n].                            (5)

Here a path starting at a consists of a and any increasing subset of the later chain vertices, and the complementary path is unique. The sole chain Hamiltonian word ends at a only when i=n. Formula (5) is 1 in that case and `1+2^m` for m=n-i>=1 otherwise.

Conversely, take a chain of order m+1, let r beat only its first vertex, and let every other chain vertex beat r. The directed cycle consisting of r followed by the whole chain proves strongness, while (5) gives `2^m+1`. A singleton gives H=1. QED.

**Concrete normalization failure.** H15 is attained, for example, by the strong cyclic extension of C3 by two singleton blocks, whose outgoing rows are

    [10,12,9,16,7].

But 15 is neither 1 nor one more than a power of two. No tournament representing H15 can simultaneously have root-deletion count one and zero merge slack. This remains false even if arbitrary gadgets or arbitrary reconstructions are allowed. It is an obstruction to that normalization, not an absent Hamilton count.

## 5. A finite-bank bounded-deletion sparsity theorem

The remaining results concern a **specified finite grammar**, unlike Theorems 2-3. Start from the singleton state and permit K,D and a fixed finite bank of rooted root-switch gadgets, including A (the singleton gadget). In state `(h,e,R,d)`,

    K: (h+R-e, d, 2R-d, d),
    D: (h, 0, 2R, d).

A gadget with fixed physical signature `(g,s,L,E,Q,B,delta,d_G)` acts by

    h'=gh-se+LR,
    e'=Eh-delta e,
    R'=Qh+BR,
    d'=d_G h.                                      (6)

The bank may be the enlarged 67-gadget bank, but no enumeration of that bank is required for the theorem.

**Theorem 4.** Fix an integer D>=1. The number of distinct counts at most X attainable in this grammar with a marked representation satisfying

    d<=D, f=0

is `O((1+log X)^2)`, with the constant depending on D and the finite bank.

More precisely, these counts lie in a finite union of fixed-value sets, geometric families `d+C*2^j`, and families

    d+2^j(A+B*2^k),   A>=0, B>0, j>=0, k>=1,       (7)

with finitely many fixed integer parameter triples (d,A,B).

### Proof: bounded-count prefixes have only a dyadic unbounded coordinate

Take

    C_0=max(2, Q+B over all gadgets).

Every non-D operation satisfies `h'>=max(h,R)` and `R'<=C_0 h'`. The lower bounds follow either from the physical induced extension and the root-switch formula, or directly from `L>=g>=s`, `R>=e`, and `h>=e`. The upper bound follows from (6), and from `2R-d<=2h'` for K.

Hence any reachable input with h=t<=D is either a bounded state immediately following a non-D operation (or the seed), or that state followed by k neutral sinks. It therefore has

    R=r_0*2^k,  1<=r_0<=C_0 D,
    e=0 if k>=1,

while t,d and the e value at k=0 range over finite sets. This is an upper envelope; not every tuple in it is asserted to be reachable.

### Proof: the final K/D suffix when slack is zero

Write the nonnegative coordinates

    a=e, b=d-e, c=h-d, f=R+e-h-d.

K preserves f, and D sends it to `2b+c+2f`. Thus if the final state has f=0, every state after the last root switch already has f=0. Any D in that suffix is possible only from `b=c=f=0`, namely the state `(d,d,d,d)`.

K fixes this state. After one D the state is `(d,0,2d,d)`; a second D would create positive slack. After j>=1 subsequent K steps its count is `d+d*2^j`, and a later D would again create positive slack. Therefore a zero-slack suffix containing D gives only the bounded value d or that geometric family.

A suffix with no D is `K^j`. Starting from a zero-slack state z=(h_0,e_0,R_0,d), it gives

    h_j=d+2^j(h_0-e_0),    j>=1,                  (8)

while j=0 gives h_0. This follows from `R_0=h_0+d-e_0`; after the first K, the mark has e=d and R=h.

### Proof: combine the last root switch and the suffix

If no root switch occurred, the preceding K/D suffix analysis starts from the singleton and already gives a geometric family.

Otherwise, let G be the last root switch. Its output deletion count is `d=d_G t<=D`, so its input has t<=D. By the prefix envelope, its input R is `r_0*2^k`. If k=0, all input coordinates are bounded independently of X, so only finitely many output states occur; (8) gives finitely many geometric families.

If k>=1, the input endpoint count is zero. Formula (6) gives

    h_0=gt+Lr_0*2^k,
    e_0=Et,
    d=d_G t.

With no terminal K, its scalar value is

    d+(g-d_G)t+Lr_0*2^k.

With j>=1 terminal K steps it is

    d+2^j((g-E)t+Lr_0*2^k).

These have form (7): `g>=d_G`, `g>=E`, and `Lr_0>0`. There are finitely many gadgets, t<=D and r_0<=C_0D, so only finitely many coefficient triples. The allowed D-containing suffixes were already handled separately.

Finally, a value in (7) at most X has `B*2^(j+k)<=X`, hence `j+k=O(log X)`. There are `O((1+log X)^2)` possible pairs. The geometric and finite parts are smaller. This proves the theorem at all orders and all derivation depths. QED.

The finite-bank hypothesis is essential to this counting argument. The universal unit-deletion classification in Theorem 3 has no such hypothesis.

## 6. Infinitely many already-attained counts cannot be normalized this way

**Corollary 5.** For every fixed D and every fixed finite gadget bank containing A, infinitely many counts already attainable by A,D alone admit no representation in the full grammar with both d<=D and f=0.

**Proof.** Start from any state with `e=d` and `R=h`. For k>=1 the actual graph macro `D^k A` gives

    ((1+2^k)h, h, (1+2^k)h, h).                    (9)

It returns to e=d,R=h. The singleton starts in this class. Using k=1,2,4 therefore realizes every product

    3^a 5^b 17^c,   a,b,c>=0.

The factors 3,5,17 are distinct primes, so these counts are distinct for distinct exponent triples. Taking each exponent at most `floor(log X/(3log 17))` gives `Omega((1+log X)^3)` attained counts up to X. Theorem 4 allows only `O((1+log X)^2)` counts with the proposed bounded-deletion zero-slack normalization. Their difference is unbounded, proving the claim.

The nontrivial product witnesses can be chosen strong: from a strong input, the D^k A construction has a new root leading to the old root, while the added chain supplies a return route. The first such construction from the singleton is also strong. This does not assert that every possible representation of these counts has the same graph structure. QED.

Thus, within a fixed finite gadget bank, a cofinite or all-target proof cannot keep the root-deletion coordinate bounded while repeatedly returning to zero slack. It must handle unbounded deletion data, a different slack class, multiple nonreturning root moves, or reconstruction outside that normalization. This is a concrete obstruction to a natural induction architecture, not a disproof of the enlarged grammar.

## 7. Verification

`verify.py` independently computes subset Hamiltonian counts and ordered rooted-cover counts on every labelled tournament through order six; it checks (2)-(3) and the unit-deletion classification throughout that universe. Capacity bounds and literal two-endpoint injections are checked through order five, and literal nonmergeable-cover slack through order four. It also checks explicit H15 and multiplicative product witnesses, and a finite-instance implementation of the grammar envelope. Normal and optimized Python must agree. The all-order density statements follow from the proof, not the finite universe.
