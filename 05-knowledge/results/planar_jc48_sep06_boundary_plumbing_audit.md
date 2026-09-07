# Independent audit: the equality sextic infinity-plumbing obstruction

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS; FINITE-EXACT replay PASS.**
Auditor: `three_ray_geometry`, 2026-09-06. This audits the main theorem in
[the producer note](planar_jc48_sep06_boundary_plumbing.md), its actual
infinity topology, and its exact finite controls. It does not rely on the
optional literature-recovery appendix or assert an affine group computation.

## Scope and actual curve

The curve is the image of
\(U=t^4+t^3+t^2\), \(V=2t^6+3t^5+2t^3+2t^2\).
The [independently audited resolution supplier](planar_jc48_sep06_resolution_budget.md)
provides its birational normalization, finite (2,5) cusp and four ordinary
nodes, and infinity (2,9) branch. The strict Nori inequality is unavailable:
the resolved value is exactly \(D^2=8=2N\).

I checked that the present proof uses the actual total infinity divisor and
its marked curve meridian. It proves that every homomorphism from that
infinity-link group into \(A_6\) carrying the curve meridian to a single
three-cycle fixes a label. The infinity epimorphism then forbids the
transitive monodromy required by the inherited sole-support Keller
passport. It does not conclude that the affine complement is abelian.

## Independent infinity and sign reconstruction

The six centres are successively on
\(L\), \(L\cap E_1\), \(L\cap E_2\), a free point of \(E_3\), a free
point of \(E_4\), and \(E_4\cap E_5\). They give precisely
\(L-E_3-E_2-E_1\) and \(E_3-E_4-E_6-E_5\), with the curve arrow on
\(E_6\). The weights in the order \((L,E_1,E_2,E_3,E_4,E_5,E_6)\) are
\((-2,-2,-2,-2,-3,-2,-1)\).

I reconstructed these centres from the supplier's infinity charts. Writing
\(X=U(1/z)/V(1/z)\), \(Z=1/V(1/z)\), the branch has initial orders
\((2,6)\); successive orders are \((2,4),(2,2),(2,2),(2,1),(1,1),(1,1)\).
The strict line leaves the branch at the third blowup: in
\(\rho=Z/X^3-4\), it meets at \(\rho=-4\), while the branch meets at
\(\rho=0\). With \(\tau=\rho/X+30\), the final centre is the satellite
\(E_4\cap E_5\); the last transverse curve point has coordinate
\(X/\tau^2=1/2450\), away from the divisor intersections. The finite
cusp blowups do not belong in this infinity graph.

Blowing up inside a tubular neighbourhood of the original projective line
does not change its boundary. Removing the final transverse curve disk
marks the actual knot exterior in the large sphere. Thus the graph is not
being substituted for an unproved affine or projective presentation.

For positive complex normal fibres, a rational vertex of square \(-w\)
has incident base-loop product equal to its fibre to power \(w\). The
crossing identifications make adjacent fibres commute. Tree access paths
introduce no stable letters. In the declared cyclic orders this gives
exactly the producer's seven vertex relations and seven commutators,
including \(f=de\mu\), with positive curve meridian \(\mu\).

An independent marked abelianization check is
\[
(l,a,b,c,d,e,f)=(-6,-4,-8,-12,-10,-9,-18)\mu.
\]
It agrees both with every displayed relation and with the blowup valuation
recurrences: start \(L=-6\mu\), then add the incident divisor valuations
and the strict-curve multiplicity \(2,2,2,2,1,1\) at the six centres.
In particular the original line meridian is \(-6\mu\), the correct
projective degree-six relation. This check detects the load-bearing arrow
and self-intersection signs. The determinant of the negative intersection
matrix is \(-1\), also consistent with the boundary of the blown-up line.

## Independent group argument

I checked the reductions \(c=l^2=a^3\), \(b=a^2\),
\(d=(lb)^{-1}c^2\), \(f=c^{-1}d^3=e^2\),
\(\mu=(de)^{-1}f\). The common square/cube cycle types in \(A_6\)
are identity, double transposition and five-cycle. There are six even
cycle types on six letters; the producer repaired an initial editorial
miscount before freezing.

If \(c\) is a double transposition, \(a=c\), \(b=1\), \(l\) has
order four, and \(d=f=l^{-1}\); this is impossible since no order-four
element is a square in \(A_6\). If \(c\) is a five-cycle, its cyclic
centralizer forces \(l=c^3,a=c^2,b=c^4,d=1,f=c^{-1},e=c^2\), so
\(\mu=c^2\) has order five. Hence \(c=1\). Only after this step is
\(f=d^3=e^2\) both a cube and a square. Its nonidentity possibilities
force respectively an order-four or order-five \(\mu\), so \(f=1\).

The remaining two groups \(\langle l,a\rangle\) and
\(\langle d,e\rangle\) are quotients of the (2,3,3) triangle group.
The elementary three-conjugate-involutions proof identifies this group
as \(A_4\), with possible images \(A_4,C_3,1\). A faithful six-letter
\(A_4\) action containing a single three-cycle cannot be transitive:
a transitive action has order-two stabilizers and no order-three element
fixes a point. A nontransitive faithful action has its natural four-letter
orbit and two fixed points. The \(C_3\) case has support three.
The second triangle group contains the prescribed \(\mu\), forcing
its support to have size at most four; its nontrivial \(d\) is a single
three-cycle. The first triangle group contains this same \(d\), so its
support also has size at most four. Their supports overlap in the three
letters moved by \(d\), and their union therefore has size at most five.
No transitivity assumption was used in proving this obstruction.

## Primary infinity epimorphism

I independently read Leidy--Maxim, *Higher-order Alexander invariants of
plane algebraic curves*, Theorem 4.7 and its proof on PDF page 11 of the
[author PDF](https://people.math.wisc.edu/~maxim/hoA.pdf).
The statement begins with any affine plane curve. The proof explicitly
uses the epimorphism from the link-at-infinity group to the affine
complement group: a generic projective line sufficiently close to the line
at infinity lies in its tubular neighbourhood, and the Lefschetz
surjection factors through that neighbourhood complement. The extra
transversality hypothesis appears only in the following Corollary 4.8.
Thus the tangent infinity line in the present sextic does not invalidate
this input. The marked curve meridian maps to the affine meridian up to
conjugacy. No claim of centrality for the killed projective meridian is
needed or inferred.

## Exact source and frozen replay

I read the complete standalone source. Its universe is all 360 even
permutations on six raw labels. It enumerates every common square/cube
value \(c\), every square root \(l\), every cube root \(a\), then
reconstructs \(b,d,f\), enumerates every square root \(e\), and
reconstructs \(\mu\). All original vertex relations and commutators
are checked before retaining the single-three-cycle arrow. This is an
exhaustive elimination, with no local-transitivity or orbit filter.

The 4,000 retained assignments have fixed-label counts
\(40,1800,2160\) for respectively three, two and one fixed labels.
Every one satisfies \(c=f=1\), with the two triangle supports and their
union checked directly. The 400 separate triangle controls, graph
reconstruction, determinant, marked abelianization, and an actual
nonabelian \(A_5\) image of order 60 fixing exactly one label provide
independent positive and hostile controls. The \(A_5\) witness prevents
an unsupported abelianity conclusion.

Independent normal and optimized runs were both byte-for-byte equal to
the saved output: **20,419 always-active gates, 1,936 output bytes**.

```text
python3 -B 04-computation/planar_jc48_sep06_boundary_plumbing.py
python3 -B -O 04-computation/planar_jc48_sep06_boundary_plumbing.py
source SHA256 5d0cf9812c611bbc1d467800e8902fee22bb1019c6ced34f1d1fbb03ce838f8a
output SHA256 054c3d5613625f4e306883e989104a87951b7bea11f18fd6f160c81da496db6b
```

The proof's main topology, algebra, exact source and affine consumer pass.
The optional appendix's separate classification and Orevkov recovery are
not dependencies of this conclusion and are outside this audit's
independent literature-verification scope. No unresolved correction
remains in the audited theorem.
