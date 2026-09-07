# The equality sextic: a marked infinity-plumbing obstruction

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** The conclusion concerns one explicit whole
nonproperness curve, not the planar Jacobian conjecture in general. No
projective-to-affine centrality assertion is used.

## 1. Actual object, inheritance, and conclusion

Let C be the affine image of

\[
 U(t)=t^4+t^3+t^2,\qquad
 V(t)=2t^6+3t^5+2t^3+2t^2.
\]

The independently audited [resolution budget](planar_jc48_sep06_resolution_budget.md)
and its [audit](planar_jc48_sep06_resolution_budget_audit.md) give a birational
normalization, one finite ordinary (2,5) cusp, four ordinary nodes, and one
point at infinity. Its projective sextic has singularities A4+A8+4A1;
the infinity branch has local pair (2,9) and contact six with the line at
infinity. The resolved Nori inequality is exactly at equality, D²=2N=8.

**Marked obstruction.** Every homomorphism from the actual infinity-link
complement group to A6 that sends a curve meridian to a single three-cycle
fixes a label. Consequently no surjective homomorphism

\[
 \pi_1(\mathbb C^2\setminus C)\longrightarrow A_6
\]

can send a curve meridian to a single three-cycle.

If C were the whole nonproperness curve of a nonautomorphic planar Keller
map, the [odd-cusp passport](planar_jc48_sep06_odd_cusp.md) and the
[alternating supplier](planar_jc48_sep06_alternating.md) would force exactly
this surjective six-sheet representation. Therefore **this explicit C is
excluded as a whole Keller nonproperness curve**. The
[independent topology, algebra and source audit](planar_jc48_sep06_boundary_plumbing_audit.md)
and root analytic read both pass. The statement does not assert that the affine complement group
is abelian.

The live concept board and its proof coordinates are:

| Inherited object | Retained information | Decisive test here |
|---|---|---|
| Strict resolved Nori mechanism | Actual affine complement | Its strict inequality fails at equality for C. |
| Semilocal cusp/node A6 passport | Marked finite meridians, retained sheets | It lacks the infinity group relations. |
| Tangent infinity divisor | Weights, adjacency, curve meridian | Seven-fibre plumbing presentation below. |
| Projective sextic classification | Quotient after killing the infinity meridian | It cannot alone identify the affine group. |
| Two-pair splice determinant method | Full boundary and a one-preimage hypothesis | Its additional hypothesis remains unnecessary for the new obstruction. |

The hostile to an overly strong conclusion is an exact marked plumbing
representation with image A5 fixing just one label. Thus the relation system
permits a nonabelian image and the one-fixed-label conclusion is sharp.

## 2. The actual infinity divisor and its marked presentation

Only the six infinity blowups are used here. Start with the projective
line L, of square +1. In order, the centres lie on

\[
 L;\quad L\cap E_1;\quad L\cap E_2;\quad E_3;
 \quad E_4;\quad E_4\cap E_5.
\]

The fourth and fifth centres are free points of the indicated component.
These are the centres in the explicit charts of the resolution-budget
proof. The total infinity divisor has edges

\[
 L-E_3-E_2-E_1,\qquad E_3-E_4-E_6-E_5,
\]

and the strict curve meets E6 transversely at a point away from all these
edges. The self-intersections, in the order L,E1,E2,E3,E4,E5,E6, are

\[
 (-2,-2,-2,-2,-3,-2,-1).
\tag{1}
\]

All seven compact components are rational. A small plumbed neighbourhood
of the total infinity divisor is the blowup of a neighbourhood of the
original line, and its outer boundary is the same large sphere. Removing
the transversal curve disc removes the corresponding knot tube. Hence
its marked boundary complement is the actual link-at-infinity exterior.
Blowups take place inside the neighbourhood; they do not replace the
affine complement or change the meridian of the strict curve.

Here is the presentation with its sign and access convention. Use positive
complex normal-circle fibres on the compact components. On a rational
component of square -w, remove small base discs at the incident plumbing
points; a section of its normal-circle bundle gives a product of the
incident boundary loops equal to the w-th power of the normal fibre.
At a plumbing crossing the two normal fibres commute, and the boundary
loop in one base is the neighbour's fibre. At the curve arrow the boundary
loop is its positive meridian μ. Choose the base access paths along the
tree, rooting it at E3. There are no graph cycles and therefore no extra
stable-letter generators. At E3 choose the order L,E2,E4; at E6 choose
E4,E5,arrow. This fixes the section/access gauge. Changing access conjugates
the marked generators and does not change their cycle types or whether the
whole image fixes a label.

Write the fibres on L,E1,E2,E3,E4,E5,E6 as l,a,b,c,d,e,f respectively.
The vertex relations and crossing commutators are

\[
\begin{aligned}
 l^2&=c,&a^2&=b,&b^2&=ac,\\
 c^2&=lbd,&d^3&=cf,&e^2&=f,&f&=de\mu;
\end{aligned}
\tag{2}
\]

\[
 [l,c]=[a,b]=[b,c]=[c,d]=[d,f]=[e,f]=[\mu,f]=1.
\tag{3}
\]

The sign in (2) is minus the self-intersection, not the self-intersection.
Reversing the orientation used to view the outer three-manifold does not
reverse individual relations. As a check, setting μ=1 in (2)–(3) fills
the knot tube and gives the three-sphere group; the exact A6 check below
then has only the identity assignment. The minus-intersection matrix of
(1) has determinant -1, as it must for this total infinity tree.
An additional marked sign check in the abelianization is
\[
 (l,a,b,c,d,e,f)=(-6,-4,-8,-12,-10,-9,-18)\mu.
\]
In particular the original infinity-line meridian is -6 times the curve
meridian, consistent with the degree-six projective relation.

Eliminating b gives c=a³. Thus every represented solution satisfies

\[
 c=l^2=a^3,\quad b=a^2,\quad
 d=(lb)^{-1}c^2,\quad f=c^{-1}d^3=e^2,\quad
 \mu=(de)^{-1}f.
\tag{4}
\]

The proof below uses only necessary relations (2)–(3), with this actual
marked meridian; it does not assume local transitivity.

## 3. Every marked A6 image fixes a sheet

The possible element orders of A6 are 1,2,3,4,5. An element that is both a
square and a cube has type identity, (2,2,1,1), or (5,1). This follows
immediately by taking powers of the six possible even cycle types on six
letters. In particular, no element of order four is a square in A6.

First consider c=l²=a³.

* If c has order two, a=c and l has order four. Consequently b=1,
  d=l⁻¹ and f=c⁻¹d³=l⁻¹. This is impossible because f=e² would have
  order four.
* If c has order five, l=c³ and a=c². Then b=c⁴, d=1, f=c⁻¹,
  e=c² and μ=c². This contradicts the prescribed order three of μ.

It follows that c=1. Now f=d³=e² is again both a cube and a square.
If f has order two, d=f and e has order four, whence μ=e⁻¹ has order
four. If f has order five, d=f² and e=f³, whence μ=f has order five.
Therefore

\[
 c=f=1.
\tag{5}
\]

The remaining relations give

\[
 l^2=a^3=d^3=e^2=\mu^3=1,\qquad
 d=(la^2)^{-1},\qquad de\mu=1.
\tag{6}
\]

Both H1=⟨l,a⟩ and H2=⟨d,e⟩ are quotients of the triangle group
⟨x,y : x²=y³=(xy)³=1⟩=A4. For completeness, conjugate x by the
three powers of y, obtaining involutions x0,x1,x2. The relation
(xy)³=1 gives x0x1x2=1; thus x0 and x1 commute and generate a
normal Klein group or the trivial group. The resulting group is A4,
C3, or trivial. The usual permutations (12)(34),(123) realize the
order-twelve possibility, proving the asserted triangle presentation.

We need its actual permutation action, rather than its abstract name.
An A4 subgroup of A6 containing a **single** three-cycle acts on at
most four labels. Indeed a transitive six-label A4 action has point
stabilizer of order two; no order-three element can fix a point in that
action. In a nontransitive faithful action, an orbit of size four leaves
two fixed points, since A4 has no quotient C2. If every orbit has size
at most three, the normal Klein subgroup acts trivially on each orbit,
contradicting faithfulness. A C3 image generated by a single three-cycle
acts on its three labels.

H2 contains the prescribed single three-cycle μ, so its support has size
at most four. Also d≠1, since d=1 would force μ=e⁻¹ of order dividing
two. Every order-three element in the resulting natural A4 action or
C3 action is a single three-cycle. Thus d is a single three-cycle.
H1 contains this same d and also has support of size at most four.
Their supports share the three labels moved by d. Consequently

\[
 |\operatorname{supp}\langle H_1,H_2\rangle|
 \le 4+4-3=5.
\tag{7}
\]

All generators in (2) lie in ⟨H1,H2⟩, so the entire marked plumbing
image fixes at least one of the six labels. This proves the obstruction
without an enumeration or a source compactification hypothesis.

## 4. The actual affine consumer: no generic-infinity assumption

The needed inclusion

\[
 \pi_1(S^3_R\setminus(C\cap S^3_R))
 \twoheadrightarrow \pi_1(\mathbb C^2\setminus C)
\tag{8}
\]

is explicitly proved for an arbitrary affine plane curve in the proof of
Leidy–Maxim, *Higher-order Alexander invariants of plane algebraic curves*,
Theorem 4.7, author-PDF page 11. Their argument uses a generic projective
line inside a small neighbourhood of the chosen infinity line and the
Lefschetz surjection for the complement of the union of both divisors.
It does not assume the chosen infinity line is generic or transverse.
That additional hypothesis first appears in the following Corollary 4.8.
[Primary author PDF](https://people.math.wisc.edu/~maxim/hoA.pdf),
[publication DOI](https://doi.org/10.1155/IMRN/2006/12976).

The map (8) takes the arrow μ to a meridian of the same irreducible
affine curve, up to conjugacy. Compose it with the necessary actual A6
monodromy from Section 1. Surjectivity would make the infinity image A6,
whereas (7) says that image fixes a label. This is the contradiction.

No abstract passport is promoted to an actual source. Only the necessary
representation of a hypothetical source is used. In particular neither
smoothness of a finite normal envelope nor the public one-dicritical
smoothness claim is a dependency.

## 5. Complete exact universe and reproducibility

The [standalone source](../../04-computation/planar_jc48_sep06_boundary_plumbing.py)
uses only Python's standard library. It constructs the infinity graph by
the six blowups, checks its exact rational determinants, and enumerates
**all A6 assignments** to (2)–(3) with arrow type (3,1,1,1), with no
transitivity filter or quotient by conjugacy. Equation (4) is the complete
elimination: l and a range over every compatible square/cube root of c;
e ranges over every square root of f; all relations are finally rechecked.

There are exactly 4,000 labelled assignments, 100 for each of the forty
possible arrows. Every assignment has c=f=1. The numbers fixing exactly
one, two, or three labels are respectively 2,160, 1,800, and 40. The
independent triangle-action control checks all 400 relevant labelled
triangle assignments. No topological claim is delegated to these counts.

A positive control, written as zero-based permutation tuples, has

\[
 l=(0,1,3,2,5,4),\quad a=(0,1,2,4,5,3),\quad
 e=(0,2,1,4,3,5).
\]

Its other fibres follow (4), its arrow is (0,2,4,3,1,5), and its image
has order 60, fixing exactly label zero. This refutes the stronger
claim that these relations force a cyclic or abelian image. It is a
representation of the infinity group, not an actual Keller cover.

Reproduce from the worktree root:

```bash
python3 -B 04-computation/planar_jc48_sep06_boundary_plumbing.py
python3 -B -O 04-computation/planar_jc48_sep06_boundary_plumbing.py
```

The frozen [output](planar_jc48_sep06_boundary_plumbing.out) records 20,419
always-active gates. Normal and optimized runs are byte-identical (1,936
bytes). The complete labelled-assignment semantic SHA256 is
`6a7b588be824405ada82187b35012650dc96d58da6c793cfaa2ec32e9d84fe9d`.

* Source SHA256: `5d0cf9812c611bbc1d467800e8902fee22bb1019c6ced34f1d1fbb03ce838f8a`.
* Output SHA256: `054c3d5613625f4e306883e989104a87951b7bea11f18fd6f160c81da496db6b`.

## 6. Literature recovery and the paths that do not prove the conclusion

**Projective group, CITED.** Akyol–Degtyarev's *Geography of irreducible
plane sextics*, Corollary 2.8, pp.8–9, gives
π1(P²\Cbar)=C6 here. The hypotheses are checked as follows: μ=16;
the set A8+A4+4A1 has weight three, excluding every torus class; it is
absent from the D10/D14 special lists in Theorems 2.13–2.16 and from
the nonabelian exceptions in Corollary 2.8. Section 5.2 establishes the
statement for every curve in these strata.
[Primary institutional PDF](https://archive.mpim-bonn.mpg.de/id/eprint/2160/1/preprint_2013_63.pdf).

Filling the infinity line kills the normal closure of its meridian λ:
π_aff/⟨⟨λ⟩⟩=C6. Thus any surjective A6 representation would have
nontrivial ρ(λ), normally generating A6. This follows because its quotient
by that normal closure is an abelian quotient of the perfect group A6.
The fact that the projective group is cyclic does not make λ central.
Section 4 instead uses the actual infinity link surjection.

**Two characteristic pairs at infinity, VERIFIED.** The graph in Section 2
has two trivalent vertices, E3 and E6, counting the curve arrow. Its
splice-arm determinants are (2,3,1) at E3 and (9,2) at E6. This is the
two-node graph in Orevkov's convention, even though the plane singularity
at the infinity point has the single local pair (2,9). Independently,
formal expansion gives

\[
 V-2U^{3/2}+\frac{15}{4}U=\frac{35}{8}t^3+O(t^2).
\]

Writing U=r⁴, the essential descending exponent gcd sequence is
4,6,3, hence 4→2→1. Local singularity data and the marked affine
boundary are different objects.

**Conditional determinant route, not a dependency.** Orevkov's
*Counter-examples to the “Jacobian Conjecture at Infinity”*, §2.4,
pp.11–13, also assumes that some irreducible component L_j of the
**target infinity divisor** has exactly one noncontracted irreducible
component in its full preimage. Contracted preimage curves are allowed.
His equations give N=m²xd1d2²R1S1, with positive integral factors,
R1≥2 and x=a+n≥n+1. The degree-six passport m=1,n=3 would then
require N≥8. [Primary author PDF](https://www.math.univ-toulouse.fr/~orevkov/jci-e.pdf).

The public [two-pair elimination note](https://github.com/alok/jacobian-two/blob/main/docs/two-pair-infinity-elimination.md)
already records this specialization and keeps its extra hypothesis
explicit; it is an antecedent, not an established independent source for
a stronger theorem. Our unique deleted divisor lies above the **affine**
curve C and does not by itself imply the target-infinity condition.
The new marked-group obstruction avoids that unpaid implication.

## 7. Scope and next question

The source-to-target map in the successful connection is the genuine
infinity-link inclusion (8). It preserves the marked curve meridian and
the whole hypothetical monodromy image. The boundary graph retains
information lost by the projective quotient and the semilocal cusp/node
passport. Its two tetrahedral pieces share a three-cycle and therefore
cannot collectively move six sheets.

The same conclusion holds for any actual whole-support curve with this
marked infinity tree and the inherited (2,5)-cusp-plus-at-least-two-nodes
passport. The tree and arrow must be supplied, not inferred from the local
infinity cusp label alone.

This argument leaves the full affine complement group unidentified and
does not state a general exclusion of every Nori-equality curve. A useful
next question is which other actual resolved infinity trees force a
proper support for three-cycle representations; that question must retain
the marked arrow and an actual curve supplying the tree.
