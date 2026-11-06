# A finite test for one-state affine source grammars

Status: PROVED, independently audited with exact verification. This concerns append-only O/Z/J source words, not arbitrary tournament rewrites, multistate encodings, nonlinear source profiles, or the unrestricted spectrum.

Use the inherited positive source cone coordinates

    (a,f,e,p)=(D-E-d,U-D+E-d,E,d-E)>=0.

The physical count after closing is F=2a+f+3e+3p. The generators are

    O(a,f,e,p)=(a+e,a+2f+e+2p,0,e+p),
    Z(a,f,e,p)=(2a+f+e+p,f,e+p,0),
    J(a,f,e,p)=(f+p,0,a+e+p,e).

These are inherited from Session200/synchronized-joint-inverse/PROOF.md. The diagnostic below is a direct deduction and makes no literature-priority claim.

## Growth bounds the complete word search

Every nonnegative v satisfies F(Mv)>=(3/2)F(v), for each generator M. In fact

    2 F(Ov)-3F(v)=f+3e+p,
    2 F(Jv)-3F(v)=f+3e+p,
    2 F(Zv)-3F(v)=2a+3f+e+p.

Therefore any nonzero nonnegative eigenray M_W v=qv obeys

    (3/2)^length(W)<=q.

This is an exact finite bound on every word with that eigenray, including nonminimal words. It is not a bound on source order for arbitrary scalar representations.

## Exclusion of four one-state digit grammars

Fix a source family s(u)=u v+w for all sufficiently large integers u, with v nonzero and with every s(u) in the nonnegative cone. Then v>=0. Suppose for base q there are physical source words W_d satisfying

    M_(W_d) s(u)=s(q u+d)

for every sufficiently large u and for a digit set containing one representative of every residue modulo q. The family s is the same for every branch. Comparing coefficients forces M_(W_d)v=qv. For a fixed actual matrix M, comparing constant terms gives Mw=w+dv, so that matrix supports at most one digit d.

For q=3,5,7,9 the preceding bound limits word lengths to2,3,4,5 respectively. Exact rational enumeration gives at most1,1,0,2 distinct actual matrices, respectively, on any common nonnegative eigenray. Each is less than q. Hence no such one-state affine source digit grammar exists at these four bases.

The exact ray lists and matrix word classes are in SOURCE_RAYS_VERIFIED.json. `verify_source_rays.py` checks every word in the proved length universe using integer matrices and rational nullspaces, rejects mixed-sign eigenvectors, and confirms that all nonzero kernels encountered are one-dimensional. It performs no floating-point interpolation. Normal and optimized Python output agree.

This is not a barrier to several source types with different rays, several parameter coordinates, non-affine profiles, or edits/deletions of old tournament arcs and vertices. It does not rule out an unconditional eventual spectrum theorem by another mechanism.
