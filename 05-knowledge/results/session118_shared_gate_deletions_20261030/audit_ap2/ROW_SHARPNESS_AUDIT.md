# Independent final audit of AP1 and GP sharpness

**Verdict: both proposed constructions pass for every M>=1.** This is an
independent algebraic and DAG-incidence audit, not another computational run
or certification of the full provisional elimination theorem.

## AP1: three deletions are sharp

Use f=pw+rz+IP_M on n=2M+4 coordinates, and the marked pair (p,r) excluding
(0,0). The circuit is G=pw, H=rz, K=G+H, output K+IP_M, with the standard
2M-1-gate inner-product circuit. It has 2M+3=n-1 gates. Its quadratic polar
form is nondegenerate, so every constant affine flat has dimension at most
M+2; equality is attained by fixing p=r and all IP x-coordinates to zero.
Thus d=M+2 is valid and sharp on the full cube.

Both paired coordinates have distinct unique children G,H; they occur in no
affine equation. There are no dead, degenerate, unary, or preparable affine
gates. G may be first in a topological order. Its other parent w is unpaired.
The safe fixing r=1 leaves p free, and w=0 kills G. H reroutes to z and K
reroutes to z, so exactly the three old gates G,H,K disappear. The remainder
is z+IP_M and has 2M gates on 2M+1 essential inputs, hence is minimum.
The surviving p is unused; no other variable usage disappears. Therefore
Delta I=-1, Delta q=-1, Delta a=-1.

The initial safe-core dimension is 2M+3>=d+2, and the final one is
2M+2>d, including M=1. All restrictions and the terminal nonconstancy guard
are therefore legitimate. The initial core is minimum even on the marked
product domain: p has a legal sensitivity witness with r=w=1; r has one
with p=z=1; w with p=1 and r=1; z with r=1. Each IP coordinate has its usual
partner-one witness with p=r=1 and w=z=0. Hence every input is essential
on legal same-coordinate fibers and the n-1 dependency bound applies.

## GP: two deletions are sharp

Use f=xy+z+IP_M on n=2M+3 coordinates, no initial pairs, with G=xy,
K=G+z and output K+IP_M. This is a live reduced circuit with n-1=2M+2
gates and no preparable affine gate. G has exactly one child, K.
GP excludes the corner (x,y)=(1,1), makes G zero, and reroutes K to z.
Exactly G,K disappear. The residual z+IP_M is again a minimum 2M-gate
circuit. The paired coordinates x,y become unused, all other input usages
survive, so Delta I=-2, Delta q=1, Delta a=-1.

The disperser-dimension parameter is d=M+1. For precision, the quadratic
polar form has rank 2M+2 and radical spanned by the z-coordinate. A direction
of a constant affine flat is isotropic. Its dimension is at most M+2; if
it had dimension M+2, it would necessarily contain the radical, because its
image in the nondegenerate quotient has dimension at most M+1. But f changes
by one along z, forbidding constancy. Thus constant flats have dimension at
most M+1. Equality is attained by setting x=z=0 and all IP x-coordinates to
zero while leaving y and the M IP y-coordinates free.

The initial safe-core dimension 2M+3 is at least d+2, and the final dimension
2M+2 exceeds d. Output nonconstancy and every legal-domain condition hold.
All n initial coordinates are essential, so the initial n-1 circuit is
minimum, not a padding-based sharpness example.

These prove sharp local deletion counts for the stated chosen restrictions.
They do not show that every possible elimination strategy encounters these
rows at a given density, or that the examples are asymptotically hard.

## Minimum-circuit sharpness of the antichain factor

The antichain sharpness family AND_i((a_i b_i) XOR x_i) has 3r essential
inputs and 3r-1 gates: r AND seeds, r XOR boundary gates and an r-1-gate
AND tree. To witness any coordinate in one factor, set every other factor
to one and use the elementary sensitivity witness inside that factor.
Consequently the initial circuit is minimum. After all a_i=0, the residual
conjunction of the r independent x_i has minimum size r-1 (zero gates for
r=1). Exactly 2r gates disappear. Thus the uniform antichain factor two is
sharp even when both endpoint circuits are minimum.
