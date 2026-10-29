# The exact quotient-walk endpoint response, and its closure boundary

**PROVED elementary, scoped; same-agent audit.** Local results E1-E5 belong
to this file. These concern ordinary unlabelled source/partner observations
of supplied-run assembly. They do not claim that the resulting response
quotient closes every later operation on a retained Hamilton presentation.

## E1. State the interface before minimizing it

Let A be an n-vertex tournament, and F a directed spanning path forest in A,
with k nonempty path components. Reversing all arcs except the forest arcs
gives B. The aligned pair (A,B) recovers F from its agreement arcs.

An input orders the supplied components as R_1,...,R_k. Put
s_i=first(R_i), t_i=last(R_i). Assembly takes other supplied children and a
legal quotient walk, traversing each supplied run intact, and forms a full
Hamilton path P in the substituted tournament S. The observed partner
reverses every arc except those of P. Observe its ordinary isomorphism class;
the ordinary source class may also be observed.

Throughout a comparison over a fixed aligned pair, the source class does not
depend on the ordering of runs. Crucially, this observation does not retain
the Hamilton path, a cursor inside a run, or a common source/target transporter.

There are two versions:
- **Class-only response:** the observations just stated.
- **Length-aware response:** additionally observe the ordered run lengths.

The latter is often the appropriate supplied-run interface. A distinction
between them is essential for the smallest hostile example.

## E2. Exact all-size response theorem

For k>=2, the complete class-only quotient-walk response of the child is
exactly the rooted tournament

    [B; (s_1,t_1),...,(s_k,t_k)].

For the length-aware response, include |R_i| with each ordered pair.
Over a fixed aligned pair, two orderings have the same response under EVERY
such context iff these rooted records are isomorphic (and lengths match
in the second version).

**Sufficiency.** Within the child, the partner is B, independent of component
order. All additional dependence on order lies in the seam arcs connecting
a run's endpoints to other children. An isomorphism of the displayed endpoint
records extends by the identity on the other children, preserving every
cross relation and seam. It therefore identifies the unlabelled partners
for any same quotient walk and any auxiliary data. If adjacent runs of the
same child are allowed to merge, their new seam is also determined by the
ordered endpoint record; the argument still applies.

**Necessity.** The following single explicit C3 context recovers B and every
ordered endpoint pair. Its two auxiliary children are strongly connected.
Its order is

    N = n + k(n+4) + 2.

No minimal-context-size claim is made.

### A recognizable degree construction

Set L=n+3, b=1+kL, c=k+1. Use quotient A -> U -> V -> A.
The auxiliary child U on b vertices is the transitive tournament in increasing
order with only the edge between its first and last vertices reversed.
Use the same construction on c vertices for V. Both are strong because their
increasing Hamilton paths close into directed Hamilton cycles.

Partition U into an initial singleton u_0 and k successive length-L runs
U_1,...,U_k. Give V its c singleton runs v_0,...,v_k. Assemble

    u_0, v_0, R_1, U_1, v_1, R_2, U_2, v_2, ...,
    R_k, U_k, v_k.

This is a legal Hamilton path in the C3 substitution. Number U's vertices
internally as u_r, 0<=r<b. In the final partner their outdegrees are

    deg(u_0)=n+2,
    deg(u_r)=n+r                 (1<=r<=b-2),
    deg(u_(b-1))=n+b-2.

Every old A vertex has degree in [c,c+n-1]. Every V vertex has degree at
least b. These follow from

    deg_partner(x)=N-1-deg_source(x)
                   +1_(x=path start)-1_(x=path end).

In U_1 use anchor u_c. In U_j for j>=2 use the second vertex of that run.
Their indices r satisfy c-1<r<b-n, so their degrees are globally unique.
They are intrinsic vertices of the unlabelled output.

Recover each U_j in both directions from its anchor. Before its first vertex,
do not continue. At an interior u_r, its preceding run vertex is the unique
vertex of degree n+r-1 pointing to u_r. Before its last vertex, its next run
vertex is the unique vertex of degree n+r+1 beaten by u_r, except at
u_(b-1), where the expected degree is n+b-2.

Why no false candidate survives: backward recovery stays below the V degree
band; every A vertex loses to an interior U vertex. Forward recovery stays
above the A degree band; an interior U vertex loses to every V vertex.
The exceptional u_0 loses to all these interior vertices. At the final
degree repetition, the preceding vertex itself has no loop and is excluded.
These observations also handle coincidences in the low-degree band.

Remove the k recovered runs. The remaining vertices of degree at least b
are exactly V. Among the remaining low-degree vertices, u_0 is uniquely
characterized by having only one outgoing arc to V: an A vertex has c or
c-1 outgoing arcs there, at least two. What remains is the old child, with
its full relation B.

Now v_0 is the unique V vertex beaten by u_0. Each v_i, i>=1, is the unique
V vertex beaten by the last vertex of U_i. The unique A vertex beaten by
v_(i-1) is s_i. The unique A vertex beating the first vertex of U_i is t_i.
Thus every ordered endpoint pair is recovered intrinsically. QED.

The exception k=1 has no component-order question over a fixed pair. The
theorem's construction is only asserted for k>=2.

## E3. The smallest response fiber can be smaller than the faithful order fiber

Let G=Aut(A,B). Its action on F's k components is faithful: a directed path
has no nonidentity automorphism fixing its component setwise.

Let H be the subgroup of Aut(B) that permutes the set of ordered endpoint
pairs {(s_R,t_R):R in F}. Let Hbar be its permutation image on those pairs.
For length-aware observation, restrict H to preserve their length colours.

Then the exact numbers over this fixed pair are

    faithful ordered-cover classes:       k!/|G|,
    quotient-walk response classes:        k!/|Hbar|.

Indeed complete orderings have free permutation-group actions; E2 identifies
exactly the Hbar orbits. G embeds in Hbar, so each response cell contains
|Hbar|/|G| faithful classes. Restoring the lost information requires a relative
coset, not just the order of a stabilizer. Unless the inclusion is normal,
that coset space need not itself be a group.

This resolves the ordinary partner-response question. It does NOT prove that
Hbar-orbit data is closed when the next operation may inspect internal run
positions, split a run, retain the full assembled path, or demand a common
source/target transport. Those operations can test information discarded
here. Nor does it cover a separate native blue-flag observation.

## E4. Sharp small boundaries and infinite phase losses

**No loss when all runs have at most two vertices.** Every forest edge is
then exactly one non-loop endpoint pair, so any H automorphism preserves F.
Thus Hbar=G. This includes all-singleton covers, where the81/27 control keeps
its full symmetry distinction.

**The first class-only loss is at n=7.** For n<=6, if there are at most two
interior vertices, every odd-order H action fixes them individually. Moving
a component containing an interior vertex would require at least three
non-singleton components, hence at least six endpoints plus an interior.
If there are at least three interiors, at most three endpoints remain; a
non-singleton component uses two, leaving no possible odd nontrivial action
on three components. Therefore every nontrivial endpoint action comes from G.

The bound is attained by a family. For odd m>=3, take the cyclic tournament
C_m with steps1,...,(m-1)/2. Replace each vertex i by an edge x_i -> y_i,
and add z with x_i -> z -> y_i for every i. Call the result B, of order2m+1.
Take forest

    (x_0,z,y_0), (x_1,y_1), ..., (x_(m-1),y_(m-1)),

and obtain A by forest complement.

The degree classes of B are X,Y,{z}. Aut(C_m) is cyclic of order m:
fixing a vertex fixes its transitive out- and in-neighbourhoods pointwise.
Relations to X identify every Y vertex, so Aut(B)=C_m, rotating both lists
together. Its endpoint action has order m. Preserving F must fix its unique
three-vertex run, killing the rotation. Hence G=1. Moreover A is regular
of degree m and therefore strong. B is strong too.

Thus m! faithful orders collapse to (m-1)! responses. At m=3, n=7:

    A=(112,49,11,67,44,76,22),
    B=(78,12,120,48,99,3,42),
    F=((0,6,1),(2,3),(4,5)),

and six orders collapse to two. The cyclically shifted order changes the
length list, so this particular loss disappears if lengths are observed.

**The first length-aware loss is at n=8.** For n<=7 and at most two
interiors, an odd action cannot move a length>=3 component among three
same-length components: that would require at least three interiors.
It fixes the interior vertices and the long components, hence their paths.
With at least three interiors there are at most four endpoints, too few
for a nontrivial odd component action while supporting a nonsingleton path.
Thus Hbar=G in the length-aware interface through7.

For sharpness, start with C_m with each vertex replaced by x_i -> y_i,
adjoin a universal source s and sink t, and take

    (s,x_0,x_1,...,x_(m-1),t), (y_0),...,(y_(m-1)).

Here n=2m+2, G=1 and Hbar=C_m. Rotating the singleton components preserves
the length list and all endpoint responses, while it moves the hidden X
positions inside the long run. A is strong: its long-run vertices and the
arc t->s form a directed cycle, and each Y vertex receives from t and beats s.

At m=3 the literal data are

    A=(2,104,99,37,15,145,57,95),
    B=(254,156,152,240,224,198,134,0),
    F=((0,1,3,5,7),(2),(4),(6)).

There are24 faithful order classes and8 response classes. Keeping the long
run first already gives6 versus2, with the same length list on both sides.

## E5. What Frobenius/LCA transport must retain

The lost coordinate here is an actual rotation. A target-side endpoint
transporter exists, but it is not a simultaneous automorphism of A and B:
it moves a hidden part of F. This is the same precise compatibility issue
as the inherited Frobenius/LCA port work, with an explicit new realization.

LCA multiplication is sound for actual maps on a shared free frame.
Replacing those maps by marginal endpoint responses, their group orders,
or independently selected phases is not sound. The appropriate restoration
over the response fiber is the relative source-compatible coset described
in E3. A marked interior vertex kills the cyclic phase in both displayed
families, but it is extra information.

The new endpoint minimum is useful precisely because it also states its
limit. The smallest carrier CLOSED under all intended future run operations
still requires that operation menu to be fixed; E2 settles the one-stage
ordinary response and supplies both a decoder and sharp counterexamples to
faithful component-order recovery from that response alone.
