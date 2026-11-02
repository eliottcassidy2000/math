# All zero completions of ordinal sums of directed triangles

**Status: PROVED, all orders.** This classifies the nonstrong hostile used in
Session137. It is not a solution of zero rigidity for strong quotients.

Let `T=C3_0 => C3_1 => ... => C3_(k-1)`, with all edges directed from earlier
to later blocks. Delete one vertex from block r, leaving the homogeneous
ordered pair u->w. For a new vertex x let a,b be its beaten indicators on
u,w and t_i in{0,1,2,3} its number of beaten vertices in intact block i.

## Theorem

The completions with H=H(T)=3^k are exactly:

1. **Canonical:** `(a,b)=(1,0)`, t_i=0 before r and t_i=3 after r.
2. **Previous-block defect, when r>0:** `(a,b)=(0,1)`, t_(r-1)=1,
   with every other intact-block indicator canonical.
3. **Next-block defect, when r<k-1:** `(a,b)=(0,1)`, t_(r+1)=2,
   with every other intact-block indicator canonical.

Each noncanonical block profile represents exactly three labelled cuts.
Thus the exact number of labelled restoring cuts is1 when k=1; for k>=2 it
is4 at an endpoint block and7 at an interior block. The noncanonical
completions leave the C3-substitution family by the inherited canonical-card
theorem.

## 1. Explicit ordinal-sum capacities

For k>=2 put

    h=H(T-v)=3^(k-1),       L=2*3^(k-2).

The broken-pair capacity is3L, and every internal capacity of an intact
triangle is L. For two intact blocks i<j, with d=j-i,

    c_ij = h*4^(d-1) + [d=1]*3^(k-3).                   (1)

There are two intact blocks only when k>=3, so the second term is integral.
For an intact block i after r,

    c_ui = 2h*4^(i-r-1),
    c_wi =  h*4^(i-r-1)+[i=r+1]*3^(k-2).              (2)

For an intact block i before r,

    c_ui =  h*4^(r-i-1)+[i=r-1]*3^(k-2),
    c_wi = 2h*4^(r-i-1).                               (3)

These capacities are per individual pair of vertices, not block totals.
In particular the cheap incident capacity in (2),(3) is2L for an adjacent
block and at least6L for a nonadjacent block. The other adjacent incident
capacity is3L. Every intact/intact capacity is at least(5/3)L.

### Derivation

Every valid Hamiltonian word of an ordinal sum traverses its blocks in order.
For a fixed intact triangle vertex, the sum over splitting that triangle
between a prefix and suffix with that vertex fixed as the suffix start is3;
the analogous fixed-prefix-end sum is also3. The unrestricted split sum is

    sum_(A subset C3) H(C3[A])H(C3-A)=12,

whereas its ordinary H is3. For TT2 the unrestricted split sum is4 and its
ordinary H is1. Thus every intervening block, whether C3 or the broken TT2,
contributes the same split/ordinary ratio4.

For intact i<j, a reverse exposed gap from block j to block i has count
`h*4^(d-1)`: the two endpoint factors3 exactly cancel their two ordinary
block factors3. A forward exposed gap is possible only for adjacent blocks,
when its count is `3^(k-3)`. Their sum proves (1).

For the broken pair, the suffix-start split factors of u,w are2,1, while the
prefix-end split factors are1,2. The direct adjacent valid gap additionally
uses only w as a path end or u as a path start. This gives (2),(3).
Internal entries also follow from the identities in README Sections1–2.

## 2. Budget proof of the classification

The exact positive cut formula says the new H minus h is

    3L*a(1-b)
    + L*#{intact blocks with t_i=1 or2}
    + sum_(intact i<j)c_ij*t_i*(3-t_j)
    + sum_(i<r)t_i[c_ui*(1-a)+c_wi*(1-b)]
    + sum_(i>r)(3-t_i)[c_ui*a+c_wi*b].                 (4)

The target increment is exactly2h=3L. All terms in (4) are nonnegative.

First, no two intact blocks i<j can have t_i>0 and t_j<3. If both are mixed,
their internal cost is2L and their crossing cost is at least(5/3)L, totaling
more than3L. If either is not mixed, the crossing multiplicity is at least3,
so that term alone is at least5L. Hence intact profiles consist of zeros,
possibly one mixed block, then threes; there are no intact-block cut terms.

* If `(a,b)=(1,0)`, the broken pair alone spends3L. Positivity of the incident
  capacities forces the canonical orientation on every intact block.
* If `(a,b)=(0,0)`, any selected vertex before r would cost at least5L, since
  `c_ui+c_wi>=5L`. Thus none occurs. All incident terms then vanish, and the
  remaining internal cost is0 or L, not3L. The case `(1,1)` is dual.
* If `(a,b)=(0,1)`, every wrong incidence before r costs c_ui and every wrong
  incidence after r costs c_wi. Each costs at least2L; a nonadjacent one
  costs at least6L. Thus the budget allows exactly one wrong incidence in
  exactly one adjacent block. Its incident cost2L plus the mixed-block costL
  is3L, giving precisely cases2 and3. With no wrong incidence the cost is0.

This exhausts the four pair states. When k=1, C3 has three paths and its
three transitive completions have one, giving case1 directly.

## 3. Why these zeros do not lift by cyclic strongification

For a tournament U let f_U(j) count **ordered** covers of its vertices by j
nonempty directed paths. The local canonical six-vertex module
`U=C3=>C3` has profile

    (9,126,684,1656,1800,720).

Either of its two noncanonical zero types, up to reversal and triangle
rotation, has profile

    (9,150,756,1704,1800,720).

Thus the exact difference is

    (0,24,72,48,0,0).                                  (5)

These tiny coefficients can be checked by720 permutations, choosing cuts
at every invalid adjacency and optionally at valid adjacencies. The literal
`module_profile_probe.py` does this; `pathcover_probe.py` supplies an
independent subset-DP route.

For completeness, ordinal-join path-cover profiles have the positive
composition law

    f_(A=>B)(ell)
      = sum_(p,q,j>=0; ell=p+q-j)
          f_A(p) f_B(q)
          (p+q-j)!/[j!(p-j)!(q-j)!],                   (6)

with `0<=j<=min(p,q)` and positive p,q for nonempty blocks. To see this,
forget the orders on the p and q local path components; select a matching
of j A-components to j B-components to concatenate, then order the
remaining p+q-j components. Every global path crosses from A to B at most
once, so the construction is bijective.

Formula (6) preserves entrywise profile inequalities. Adjoining a C3
multiplies a difference supported first at component count2 by12 at that
coordinate (the p=2,q=1,j=1 and p=2,q=2,j=2 terms). Consequently every
noncanonical zero completion V of the theorem satisfies

    f_V(j)>=f_T(j) for all j,
    f_V(1)=f_T(1),
    f_V(2)-f_T(2)=24*12^(k-2)>0.                       (7)

Indeed V differs only in one adjacent six-vertex module, and its remaining
ordinal triangle factors apply the positive law (6).

Every fixed substitution context has H equal to a nonnegative linear
functional of the substituted module's ordered path-cover profile: contract
its maximal runs as in README Section3. Therefore replacing T by V cannot
decrease H in any such context. It **strictly increases** H whenever that
context admits a valid run-word using the module twice. In particular, for
any two nonempty tournaments A,B,

    H(C3[V,A,B]) > H(C3[T,A,B]).                        (8)

The word `module,A,B,module` already supplies a positive coefficient of
f_module(2). Thus the known nonstrong zero controls cannot be turned into
strong-quotient counterexamples simply by embedding their ordinal-triangle
module into this cyclic context. This is a specific strongification
obstruction, not a proof for arbitrary cuts or arbitrary strong quotients.

## 4. Audits

`ordinal_zero_audit.py` compares every explicit capacity above with the
independent quotient derivative compiler for k<=6 and every root. It then
checks all 669,924 block profiles for transitive quotients k<=8 (36 deletion
rows) against the **already proved** classification. This is a bounded audit
of a closed nonstrong theorem, not a larger strong-quotient search.

`module_profile_probe.py` directly enumerates720 local permutations for each
six-vertex control. `pathcover_probe.py` independently computes full cover
profiles by a subset DP, where adding a valid adjacency multiplies by1+t
and an invalid adjacency by t. Its k=2,3 controls agree with (5)--(7).
