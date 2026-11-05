# Exact physical replacement branches and a complete rebasing failure at H=43

Session 198. **PROVED:** two response identities give target-dependent physical
binary children, with explicit infinite primitive source families. **FINITE-EXACT
with proved all-order cutoffs:** all C2/C3 source-pair presentations of H=43 fail
to produce either 85 or 87 by source extensions and family switches, although both
children have explicit fresh physical constructions. Exact general coverage of
the remaining odd integers is still **OPEN**.

## 1. Actual replacement, not universal word transport

A source q=(D,E,U,d) is an actual nontrivial rooted tournament with one Hamilton
path starting at its root. Counts and cover conventions are those of the frozen
BILINEAR_C2 and C3_TWO_SOURCE dependencies. Set t=(1,0,2,1), the rooted transitive
edge, and let b=(3,1,3,1) be a rooted directed triangle. Use these physical moves:

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

O appends a sink. Z appends a vertex beating only the old root and retains that
root. J makes the same birth and retains the newborn root. Every move preserves
root-start count one. These are actual graph operations, not merely matrices.

Start with the physical C2 count

    n=N(t,Oq)=3D+12U-2d.                             (1)

Replace the source pair and gate, constructing M(b,q). Its count is

    M(b,q)=9D-4E+22U-4d
          =2n+(3D-4E-2U).                          (2)

Thus the physical response condition

    3D-4E-2U=epsilon, epsilon in {-1,+1},           (3)

produces exactly the binary child 2n+epsilon, with C=2. The move discards the
sink in Oq, replaces the left transitive edge by a triangle, and changes the gate
from N to M. It is not an extension of each original source, and no
scalar-preserving intermediate rebase is claimed.

There is a second actual replacement, using the shared-root three-cycle family:

    J0(t,Jq)=15D-3E+15U
            =2n+(9D-3E-9U+4d).                     (4)

Whenever the final parenthesis equals one, this constructs the upper binary
child with C=3. Equations (1)-(4) follow by substituting the fixed source into the
reviewed physical response rows. All sources are nontrivial, so no singleton
exception is being used.

The response conditions are essential. The unrestricted cone is not asserted to
contain a physical q for every desired n. The benefit is a genuine exact branch
on a specified physical fiber, rather than another congruence or a purported
identity on every pair of sources.

## 2. Infinite primitive fibers supplying both C2 children

Words below are applied from left to right, starting at t.

For k>=1, put r=2^(k-1) and take

    q^-_k = t . O Z^k J O Z.

Its actual profile is

    (D,E,U,d)=(36r-1,6r-1,42r+1,6r-1).             (5)

To check this at every k, the prefix O Z^k has profile
(6r-1,1,6r+1,1). Substitution in J,O,Z gives (5).
It satisfies 3D-4E-2U=-1. It is primitive: any common divisor of its four
coordinates divides D-6E=5 and U-7E=8, hence is one. Equations (1)-(2) give

    n^-_k=600r+11,       M(b,q^-_k)=1200r+21=2n^-_k-1.  (6)

For the positive fiber take

    q^+_k = t . O J O Z^k J O Z.

Its actual profile is

    (D,E,U,d)=(108r-7,18r-3,126r-5,18r-3).          (7)

The prefix O J O Z^k has profile (18r-3,1,18r+1,1), after which J,O,Z gives
(7). Now 3D-4E-2U=1. Any common divisor divides D-6E=11 and U-7E=16, so this
source is also primitive. Therefore

    n^+_k=1800r-75,      M(b,q^+_k)=3600r-149=2n^+_k+1. (8)

These constructions work for every k. They are infinite exact families, not
claims that the fibers cover all parent counts. The final output in each case
has one more vertex than the displayed N-parent construction, but Section 5
shows that one extra vertex does not imply an induced extension of that parent.

A particularly simple infinite family for (4) comes from

    q_k=t . O^(k-1) J=(x,1,x,1), x=2^k+1, k>=1.

These are strong primitive one-apex sources, and (4) gives

    n_k=15*2^k+13,      J0(t,Jq_k)=30*2^k+27=2n_k+1. (9)

The first case is the upper-child replacement 43->87.

## 3. A complete exact inverse at the parent value 43

The complete C2 theorem gives {15}, N(v,w), and M(v,w). Its physical bounds
force U_v,U_w<=floor((43-1)/6)=7. The complete finite source grammar gives
exactly these seven nontrivial profiles under that cap:

    (1,0,2,1), (3,1,3,1), (1,0,4,1), (5,1,5,1),
    (3,0,6,1), (5,2,6,3), (5,1,7,1).

Testing the unordered pairs allowed by the proved product bounds gives nine N
candidates, with counts 13,23,25,41,43,41,53,41,45, and two M candidates, with
counts 29,49. Thus, up to swapping its two sources, the only C2 source-profile
presentation of 43 is

    43=N(t,w),        w=(3,0,6,1)=O(b).             (10)

The complete C3 inverse has the conservative source bound U<=floor(43/3)=14.
There are exactly twenty source profiles under this cap, including the
singleton. Exhausting all 2,964 admissible ordered pairs across its eight forms,
including every permitted singleton boundary, gives scalar support through 43
exactly

    {15,27,29,31,33}.

Consequently no C3 source-pair presentation of 43 supplies an alternative rebase.
These inverse checks are all-order conclusions: the source bank has its proved
monotone cutoff, not a graph-size or construction-depth cutoff. The separate
regular-five C2 exception is 15, so it causes no exception at 43.

## 4. All source extensions and all family switches still miss both children

Use the positive source coordinates

    y=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E),

and the ten positive matrices P_F of the reviewed finite-grammar packet. Its
exact matrix inequalities are

    P_F >= P_N entrywise for every one of the ten forms,
    N(Xv,w) >= (4/3)N(v,w), X in {O,Z,J},

with the same bound on the right. Thus, after k total source letters, every
final C2/C3 form has count at least (4/3)^k*43. At k>=3 this exceeds 87.
Only the 34 pairs of source words having combined length at most two need
checking. Across ten final forms these give 340 tests, or 298 distinct readout
matrices. Their complete output spectrum through 87 is

    {43,69,77,79,81,83}.                            (11)

In particular neither 85 nor 87 occurs. The reversed ordering in (10) gives the
same scalar set because the source-word set and the ten forms are closed under
transpose. Sources never become singleton under these operations, so no omitted
singleton final boundary is reachable from (10).

Equation (11) applies even after choosing any source-pair presentation with the
same parent count 43: Section 3 exhausted all of them. This is stronger than a
failure at one arbitrarily chosen marked representation. It still does not
forbid replacing a source by a smaller or unrelated source, as (2)-(4) do.
It is not an exclusion of either scalar child from the tournament spectrum.

## 5. Fresh physical children and the stronger one-vertex boundary

For q=b, (1) gives 43, (2) gives 85, and (4) gives 87. Thus

    M(b,b)=85, C=2;
    J0(t,Jb)=87, C=3.                              (12)

There is also a C2 construction of the upper child:

    N(t,t . O Z J)=87,
    profile(t . O Z J)=(11,4,12,5).                (13)

All these children have eight vertices. One literal parent in (10) has rows

    [114,124,113,37,8,16,56].

The verifier constructs the reverse ordered-source presentation and explicitly
identifies it with the converse of this graph. There are no hidden rooted graph
choices in the two source profiles: H=1 forces a transitive source, and U=2
forces its order to be two. A weak start-one H=3 source has a triangle initial
component and only singleton tail components. A strong component with H=3 has
at most three vertices, since cutting any Hamilton cycle yields as many distinct
paths as vertices. Its rooted cover count is three; U=6 therefore forces exactly
one tail vertex. The triangle's root choices are isomorphic. The complete C2
normal form consequently places every C2 graph with H=43 among the two displayed
ordered constructions, which are converses. Section 3 excludes C3.

Enumerating all 128 orientations of a new vertex against the displayed parent,
without any restriction on the child's cycle count, produces no H=85 or H=87.
Converse gives the same conclusion for the other parent presentation. This is a
finite-exact exclusion of every one-vertex induced extension of these parent
classes. It is not an arbitrary-depth induced-extension theorem, and says
nothing about parent graphs with more than three cycles.

Thus the eight-vertex constructions (12)-(13) actually replace the parent
structure. They cannot be obtained merely by adding a vertex to either of its
C2/C3 source-pair graphs. This is the concrete geometric distinction between a
fresh target construction and a returning source extension.

## 6. Reproduction and remaining arithmetic task

Run `python3 verify.py` and `python3 -O verify.py` in a copied packet. The two
receipts agree. The verifier separately checks raw and positive form formulas,
all 11 C2 inverse candidates, every admissible C3 pair under the proved cap,
all short extension readouts, physical parent/child graphs by literal
permutations and subset DP, the converse map, every arbitrary one-vertex birth,
and the two replacement identities. It checks the ray profiles for k=1..100
and literal ray graphs for k=1..3. Those finite controls supplement the displayed
all-order word calculations.

positive_model.py and the dependencies/*.py files are unchanged copies from the
reviewed finite-grammar packet. Their hashes, and those of the named proof
snapshots, are recorded. No dirty checkout, external graph bank, or numerical
saturation assumption is used.

The constructive task now exposed is specific: for a desired child, find an
actual source on a response fiber such as (3) or (4), with the required exact
parent readout (1), or find a different physical replacement identity. The
explicit rays show these fibers are nonempty and unbounded even on primitive
sources. They do not prove that the required fibers meet every remaining target.
