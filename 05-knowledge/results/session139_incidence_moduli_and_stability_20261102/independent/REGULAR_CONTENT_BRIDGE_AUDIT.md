# R11/Paley11 coarse-spectrum versus context-content bridge

**Independent audit: PASS.** Only these new audit files were created; frozen
source files and earlier receipts were preserved. No Git mutation.

Direct adjacency checks use the standard cyclic connection set
`{1,2,3,4,5}` for R11 and the quadratic-residue set `{1,3,4,5,9}` for P11.
Both graphs have order11, every outdegree5, and exactly55 directed cyclic
triangles (counted once per unordered triple). The latter also follows from
`binom(11,3)-11*binom(5,2)=55`.

Both have principal adjacency eigenvalue5 and every nonprincipal adjacency
eigenvalue has real part -1/2. Indeed regularity gives A1=5*1 and
A^T1=5*1. On the invariant subspace perpendicular to1,
`A=(-I+K)/2`, where K=A-A^T is real skew-symmetric. Its eigenvalues there
are therefore `-1/2+i*mu/2` for real mu. This is an exact structural proof,
not a floating-point eigenvalue estimate; the script checks the matrix
identities, K1=0, and normality directly for both graphs.

The previously independently audited literal R11 profile has context
content33. Its exact receipt is pinned by SHA256
`5ee1af479cff6ecef4a34874a94768e3d86432d24da7bed4e0c8d0dc28c18605`.
The55 distinct affine maps x->a*x+b with a in `{1,3,4,5,9}` are checked
again for closure and arc preservation on P11. The accepted free-action
context theorem therefore gives55|g(P11). Thus R11 is not contextually zero
modulo5 (or55), while P11 is contextually zero modulo55.

**Conclusion:** order, regular degree, cyclic triangle count, and the
nonprincipal adjacency critical-line property, even jointly, do not
determine arithmetic context content or zero-mod5 contextual behavior.
This does not assert that full spectra fail to distinguish this binary
pair: direct fourth traces are880 for R11 and660 for P11, so they are
explicitly **not cospectral**. No exact value g(P11)=55 is claimed.

`audit_regular_content_bridge.py` imports no research implementation. Its
normal and `-O` outputs agree byte-for-byte; every check raises an explicit
exception. This is a two-graph bounded bridge, not a new circulant census or
a general classification of arithmetic context states.
