# Growth recurrences and the exact deletion boundary

**Status: PROVED deductions, with explicitly identified FINITE-EXACT seed
inputs. Same-agent proof review; no claim of literature priority.**

The preceding packet's regular-parent theorem R1 and clone-corrected
count R2 are inherited from
[regular_parent_symmetry/PROOFS.md](../regular_parent_symmetry_20261029/PROOFS.md).
Here a prime tournament means one of order at least three with no
nontrivial proper module; its vertex count need not be a prime integer.
All attachments preserve the old literal orientation.

## T1. Every odd carousel supplies a rigid family of size b^k+1

Let R_b be the cyclic tournament on Z/b, b>=3 odd, with arcs
x->x+j for 1<=j<=(b-1)/2. It is regular and locally transitive.

Its full automorphism group is C_b: in the transitive out-neighborhood
of x the unique source is x+1, so an automorphism preserves the
successor cycle. Rotations supply all b such maps.

It is prime. If M is a proper module, an outside vertex puts M inside
one in- or out-neighborhood, hence M is transitive. Regularity of the
whole tournament makes every module regular internally (its vertices
have the same external incidence). A transitive regular tournament is
a singleton.

Let T_k=R_b^[k]. Color a word x in (Z/b)^k by

    a(x)=sum_i 1[x_i=0] mod 2.

The number of ones is

    A_k=(b^k-(b-2)^k)/2.                                  (1)

This is the inherited parity construction, now applied uniformly to
all odd carousel bases. Its two complementary colorings have distinct
weights. At each node a color-preserving automorphism must fix the
distinguished child; the quotient C_b action then becomes trivial.
Induction eliminates every child stabilizer. Thus adding one vertex
v that beats the color-one vertices gives a rigid unmarked tournament
U_k, by the regular-parent theorem.

The maximum transitive size of R_b is t=(b+1)/2: use the source of a
transitive set for the upper bound and a consecutive interval for
equality. The inherited weighted substitution formula gives

    |U_k|=b^k+1,   Aut(U_k)=1,   tau(U_k)=t^k+1.           (2)

For the last equality, the product of the transitive sets
{1,...,t}, avoiding the distinguished symbol zero, is all color zero
and beats v. It attains the general upper bound tau(T_k)+1.

The parent group has order b^((b^k-1)/(b-1)). For prime b=p it is
a p-group; this does not assign that group to U_k, whose group is trivial.

## T2. All-order orbit recurrences for these children

Let F_k count binary coloring orbits of T_k, and R_k those with
trivial stabilizer. For prime p=b, the usual necklaces give

    F_0=R_0=2,
    F_(k+1)=(F_k^p+(p-1)F_k)/p,
    R_(k+1)=(R_k^p-R_k)/p.                               (3)

A regular coloring requires every child coloring regular and a
nonconstant quotient word. A word of prime length has a nontrivial
rotation precisely when it is constant. This proves (3).
For general odd b the corresponding classical formulas are

    F_(k+1)=b^-1 sum_(d|b) phi(d) F_k^(b/d),
    R_(k+1)=b^-1 sum_(d|b) mu(d) R_k^(b/d).               (4)

Thus the Euler-totient necklace count and the Möbius primitive-necklace
count are two different observables of the same growth operation.

The parent is vertex-transitive, so its unmarked child count is F_k-1.
For prime p, the number of rigid children is R_1-1 at k=1 and R_k at
k>=2. The point stabilizer being subtracted has order
p^((p^k-1)/(p-1)-k), trivial at k=1 and nontrivial thereafter.

| p | parent p^2 | all unmarked children | rigid children |
|---|---:|---:|---:|
| 3 | 9 | 23 | 2 |
| 5 | 25 | 6559 | 1554 |
| 7 | 49 | 182857159 | 87460002 |

The binary regular palette two is a fixed point of (3) exactly when
p=3: 2^(p-1)=p+1. For p>=5, 2^(p-1)>p+1 and the iterates grow strictly.
The exceptional ternary constancy was already proved in the previous
packet; this explains its boundary among prime carousel bases.

These are complete counts for the stated regular-parent extension
families, not counts of every tournament of that size.

## T3. Mixed bases retain their order even when scalar products agree

For prime regular bases B_i with distinguishing cuts of size s_i, let b_i be their
orders and suppose each cut complement retains a maximum transitive
set of size t_i. The same parity construction gives

    n=product_i b_i,
    imbalance=product_i (b_i-2s_i),
    cut weight=(n-imbalance)/2,
    tau(child)=product_i t_i+1.                          (5)

Full group orders instead obey the ordered recurrence

    g_outer=|Aut(B_outer)| * g_inner^(b_outer).

For R3[R5] and R5[R3] the old size is15, the cut weight is6, and
the child transitive maximum is7 in both cases. The core groups have
orders375 and1215. The two rigid16-vertex children are nonisomorphic,
verified natively. Also, their new vertex has outdegree6, whereas clone
attachments to a regular15-parent have degrees7 or8. Their regular-card
root is consequently unique, so the nonisomorphic cores already prove
the children differ. Products in (5) do not reconstruct the ordered LCA
tree.

## D1. A sharp vertex-deletion law for doubly regular tournaments

Let H be any doubly regular tournament of order q=4t+3>=7, and put
h=(q-1)/2. Then:

1. Every induced subtournament with more than (q+1)/2 vertices is prime.
2. The smallest number of vertices whose deletion makes H decomposable
   is h.
3. Exactly q(q+3)/2 deletion sets of size h attain this boundary.
4. Each boundary subtournament has exactly one nontrivial proper module:
   either a pair, or a regular h-vertex neighborhood beneath a source
   or above a sink.

This is about **vertex deletion**. It is not the classical arc-reversal
simplicity index. The literature distinction is recorded in RECOVERY.md.

**Proof of the inequality and equality module sizes.** Let S be the
skew adjacency matrix, with S_ij=1 when i->j. Double regularity gives

    S1=0,   S^T S=qI-J.                                 (6)

For a set M of size m>=2, let U be the outside vertices that see all
of M in one direction, and D the remaining outside vertices. If M is
to be a proper module after deletions, U is nonempty and all of D
must be deleted. An outside uniform witness puts M inside one
neighborhood, so m<=h.

Project the coordinate space on M orthogonally to its sum-zero
subspace. Equation (6) gives trace q(m-1) for the restricted squared
operator. Each row in U contributes zero. If d_u is the signed row
sum inside M and e_w the signed row sum from w in D to M, then exactly

    q(m-1)=m(m-1)-sum_(u in M)d_u^2/m
                      + |D|m-sum_(w in D)e_w^2/m.        (7)

For odd m, every e_w is odd. For even m, every d_u is odd. Therefore

    |D| >= m(q-m)/(m+1),                      m odd;
    |D| >= [(m-1)(q-m)+1]/m,                  m even.     (8)

Subtracting h from these bounds gives, respectively,

    (m-1)(q-1-2m)/(2(m+1)),
    (m-2)(q-1-2m)/(2m).

They are nonnegative for 2<=m<=h, and equality is possible only at
m=h (odd) or m=2 (even). This proves the lower deletion bound and the
classification of possible module sizes at equality.

For a pair {u,v}, its outside separating set D(u,v) has size h by
the doubly regular intersection numbers. Deleting it makes the pair
a module. For a vertex u, retaining {u} together with N^+(u) or N^-(u)
also gives a module of size h after h deletions. These neighborhoods
are regular, by the common-neighbor count t. Hence the lower bound
is attained.

If a boundary module has size h, an outside witness forces it to be
exactly a neighborhood, and there is only one outside surviving vertex.
If it has size two, every deleted vertex belongs to its separating set.
These account for all boundary deletion sets.

**Distinctness and the count.** The 2q in/out-neighborhood deletion sets
are distinct. Equal in/out sets at two different vertices would give
zero common out-neighbors, contradicting t>0; equal same-direction
sets are already incompatible with their mutual edge.

A cone over a regular h-tournament cannot have a twin pair: two twins
inside a regular tournament have different degrees, and a pair including
the cone point would force a source or sink in that regular tournament.

Two pair modules in the same boundary subtournament cannot overlap.
Their union would be a three-vertex module. For h>3 its size is forbidden
by (8); for h=3 equality forces it regular, whereas two overlapping
twin pairs make a transitive triple.

Nor can the pairs be disjoint. For f=e_u-e_v and g=e_a-e_b, the two
vectors Sf and Sg overlap only on D and have entries +/-2 there.
The set D has odd size h, so their dot product is a nonzero multiple
of4. But (6) makes this dot product zero. Contradiction.

There are consequently binom(q,2) different pair cuts, disjoint from
2q cone cuts, giving binom(q,2)+2q=q(q+3)/2. No second nontrivial module
is possible in a boundary subtournament. QED.

**Local consequence.** Every in- and out-neighborhood of H is itself a
regular prime tournament. A proper module inside such a neighborhood
would persist inside its cone. The boundary classification excludes
every size except two, and regularity excludes twins.

**Exact controls.** Complete candidate-module universes are checked at
q=7,11,15,19,23,27. Completeness follows because each candidate is a subset
of an in/out-neighborhood of an outside witness. The boundary cut
counts are35,77,135,209,299,405. An independent induced-subset/module
enumeration checks q=7,11. No random sample is called an exhaustive test.

The order15 control is deliberately not a finite-field Paley construction.
Adjoin a source to P7, form its skew matrix K, and set H=I+K.
The classical skew-Hadamard doubling

    H' = [[H,H],[-H^T,H^T]]

satisfies H'+H'^T=2I and H'H'^T=16I. Switch K'=H'-I so its first
row is all positive off the diagonal, then delete that source. The
remaining15-tournament is doubly regular. Its full native automorphism
group has order21 and its transitive maximum is5. This tests the
theorem outside the prime-power field family. Preservation of the group
under this doubling is already the known Hanaki theorem, recovered as
[THM-4557](../../../01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md);
the order21 computation is a control, not a new symmetry result.

## D2. A uniform multiport guard and two exact attachment radii

Let B be a fixed prime tournament of size b and let H be doubly regular
of size q. Every extension of B^[k] by r arbitrary vertices remains
H-free whenever

    q-r>b  and  r<(q-1)/2.                              (9)

An induced copy of H would meet the old core in at least q-r vertices.
By D1 that intersection is prime and has more than b vertices. A prime
tournament larger than b cannot embed in B^[k]: at a substitution node
it must be contained in one block, or meet each block at most once.
The second case is too small, and the first descends to a singleton.

The bound can be sharpened by retaining the actual prime isotypes.
Let M_B(H) be the largest size of an induced subtournament of H that
embeds in some B^[k]. The stable minimum number of new vertices needed
to introduce H is exactly |H|-M_B(H). The upper bound is realized by
adding the missing vertices of H around an embedding of its old part.
The lower bound follows by intersecting any alleged copy with the core.

Two complete finite obstruction banks yield:

| old family | forbidden H | M_B(H) | arbitrary attachments always safe | first possible introduction |
|---|---|---:|---:|---:|
| R5^[k] | P7 | 4 | 2 | 3 |
| P7^[k] | P11 | 5 | 5 | 6 |

These thresholds hold for every k>=1, since the lower witnesses already
embed directly in B. For every5-subset of P7, and every6-subset of P11,
the certificate supplies an induced prime obstruction absent from B.
Primeness makes it absent from the whole substitution tower. Direct
4- and5-vertex embeddings supply the matching lower witnesses.
No assumption about random attachments or their automorphism groups
enters these exact radii.

## D3. New local bases at 9,11,13,29

For Paley H=P_q, label the out-neighborhood of zero by powers of a
generator of the square subgroup. If its order is b=(q-1)/2 and g is
that generator, the induced cyclic connection set is

    C={d in Z/b minus zero: chi(g^d-1)=1}.               (10)

Indeed chi(g^i)=1 factors out of chi(g^j-g^i). If q=p^e,
Frobenius acts on this cyclic index as d->pd. This is a literal
orientation-preserving map with its coordinate choice retained.

D1 proves these local bases prime and regular at every q. The following
additional data are FINITE-EXACT full-group and transitive-set calculations:

| ambient q | base b | full group order | cut size s | tau(base) | b-2s |
|---|---:|---:|---:|---:|---:|
| 7 | 3 | 3 | 1 | 2 | 1 |
| 11 | 5 | 5 | 1 | 3 | 3 |
| 19 | 9 | 9 | 1 | 4 | 7 |
| 23 | 11 | 11 | 1 | 4 | 9 |
| 27 | 13 | 39 | 2 | 4 | 9 |
| 59 | 29 | 29 | 1 | 8 | 27 |

Each displayed cut has trivial stabilizer and its complement retains
a maximum transitive set. Thus it supplies all-depth rigid families
with size b^k+1, cut weight [b^k-(b-2s)^k]/2, and transitive maximum
tau(base)^k+1, by the proved parity construction.

In particular the11- and13-vertex bases give the common imbalance9^k
and transitive maximum4^k+1 at different sizes11^k+1 and13^k+1.
The base13 group is C13 semidirect C3: the action (10) supplies39
automorphisms, and native enumeration confirms there are no others.
This recovers the existing <3>-cyclotomic13 tournament from the
cycle101 work by a new concrete local realization, not a newly claimed
isomorphism class.

For a prime-sized base with full group C_p semidirect C_ell, ell an odd
prime dividing p-1, set c=1+(p-1)/ell. On Q named colors the exact
ordinary orbit and regular-orbit polynomials are

    F(Q)=[Q^p+p(ell-1)Q^c+(p-1)Q]/(p ell),
    R(Q)=[Q^p-p Q^c+(p-1)Q]/(p ell).                    (11)

There are p order-ell point stabilizers. Nonconstant colorings fixed
by two distinct ones would be fixed by the whole group, which contains
all translations, an impossibility. This proves the second formula;
the first is Burnside. Iterating R counts regular-color orbits in
uniform towers.

For P7, P11 and local13 the binary regular-orbit counts are4,36,200;
their unmarked one-vertex child counts are11,43,231. The minimum
permutation indices are4,8,8. The preceding packet's higher-node
coupling proof therefore extends: with n=p^k, k and Q tending to infinity,
n/Q^d->lambda gives residual group C_ell^Poisson(lambda). At the n/p
bottom nodes the probability of an order-ell stabilizer is
p(Q^-d-Q^-(p-1)); full-group bottom stabilizers vanish in this regime.
The inherited higher-node error tends to zero by the same cycle-index
bound. Ports are fixed individually, their incidence bits are independent
and fair, and only bounded observables are covered. For Q=2^r this is
a subsequence limit with the stated ratio, not a continuous-offset law.
