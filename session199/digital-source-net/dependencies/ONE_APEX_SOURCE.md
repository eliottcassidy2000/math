# Two fixed-left source-pair families have every admissible modular residue

Status: **PROVED**, conditional only on the named reviewed physical graft laws and
one-apex gadget compiler. Session 198. This is an all-order congruence theorem,
not an exact-target coverage theorem.

For each positive integer m and each residue r compatible with an odd integer,
there are infinitely many distinct Hamilton-path counts congruent to r modulo m
among strong tournaments having exactly two Hamilton cycles. The same assertion
holds separately for exactly three Hamilton cycles. The witnesses in each class
use one fixed two-vertex left source and an actual one-apex right source.
Thus the assertion includes mixed even moduli, not merely separate odd-modulus
and power-of-two results.

## 1. Two explicit physical families

A source is a rooted tournament with exactly one Hamilton path starting at its
root. Its response is (D,E,U,d), comprising the Hamilton-path count, root-ending
count, reserved-root-start two-cover count, and root-deletion count. All are
measured on that one physical graph.

Fix the transitive two-vertex source v=(1,0,2,1). For every L>=2 let G(L) be the
reviewed ending-one gadget: its base is transitive, its incidence word is the
complemented binary digits of L-1 in low-to-high order, and its root is 0. Its
converse W(L) is a nontrivial source with response

    w(L)=(g(L), L mod 2, L, 1).

This is an explicit graph family, not an arbitrary tuple assignment. The source
has deletion count one because its root-deletion is transitive.

Apply the reviewed C2 M graft to v,w(L), or apply the reviewed C3 Type-B graft to
the same pair. The respective left response rows (h,E_A,R_A,K_A) are

    C2: (5,3,13,2),       C3: (5,2,14,2).

The arbitrary-gadget scalar formula therefore gives

    F2(L)=5g(L)+13L-3(L mod 2)-2,
    F3(L)=5g(L)+14L-2(L mod 2)-2.                 (1)

Each is the path count of a literal strong graph with cycle count 2 or 3,
respectively. Both graphs have bit_length(L-1)+5 vertices. The fixed left
adjacency rows can be taken to be

    C2: [10,12,1,4], selected entries {0,3};
    C3: [2,12,1,5], selected entries {0,3}.

All cross arcs go from the fixed left graph to G(L), except that its root beats
both selected entries. Its ending count is one. The selected starting-count sums
are 2 and 3, so the reviewed cycle law gives the asserted cycles and strongness.

Only even L will be needed. In that subfamily,

    F2(L)=5g(L)+13L-2,       F3(L)=5g(L)+14L-2.    (2)

## 2. The inherited joint modular compiler, with its mechanism

For f(n)=g(n+1), the physical gadget compiler gives

    f(0)=1, f(2x+1)=f(x), f(2x)=f(x)+x+(x mod 2). (3)

Consequently, for every even modulus M and every pair of residues u,z with z
odd, there are arbitrarily large L satisfying L=u and g(L)=z modulo M.
This is the exact theorem in JOINT_MODULAR_GADGET.md. The following recalls its
proof so that the role of physical realizability is explicit.

Write M=2^k q, k>=1 and q odd. Choose a>=2 with 2^a=1 modulo q, taking a=2 if q=1.
Appending the binary word W=1^a 0^k to a current n=0 modulo M returns n to zero
modulo M and adds

    Delta=(2^a-1)(2^k-1)+1

to f modulo M. The first zero contributes the additional 1 in this expression.
One has Delta=1 modulo q and gcd(Delta,M)=2. For k>=2, Delta=2 modulo 4; for k=1,
Delta=2^a and its gcd with the odd q is one.

Let s=u-1 modulo M, and append the fixed-width binary word for s after j copies
of W. Starting at n=0 modulo M, that final word changes f by f(s)-1. Thus the
result has L=u and g(L)=f(s)+j Delta modulo M. Choose j modulo M/2 to obtain z.
Increasing j by positive multiples of M/2 makes the literal binary word, hence
L, arbitrarily large without changing either residue. Every resulting integer
L is fed to the actual gadget compiler.

For M=2^k the same reviewed theorem gives a word of length at most 4k^2+2k+1 for
any joint residue. Its neutral initial word makes L nontrivial. No uniform
bound on the numerical size of an exact-target witness follows.

## 3. Explicit residue selection for both cycle classes

Fix m>=1 and an admissible residue r. Choose an odd integer R congruent to r
modulo m; for odd m an even representative can be replaced by r+m. Put

    t=(R-3)/2,        M=2m.

For C2 choose

    u=-6t modulo M,       z=1+16t modulo M.

Then u is even and z is odd, and

    5z+13u-2 = 3+80t-78t = R modulo M.            (4)

For C3 instead choose

    u=-2t modulo M,       z=1+6t modulo M,

so

    5z+14u-2 = 3+30t-28t = R modulo M.            (5)

These are integral Bezout identities: 5*8-13*3=1 and 5*3-14*1=1. They require
no invertibility of 3,5,13,14 or m. Apply the joint compiler with (u,z). Since M
is even, its resulting actual L is even. Equations (2), (4), and (5) produce the
requested strong physical graphs and prove the theorem.

The compiler supplies arbitrarily large L. Since F2(L)>=13L+3 and
F3(L)>=14L+3 for even L, their counts are unbounded and yield infinitely many
distinct integers in the requested residue class.

For a target modulus 2^k with k>=1, the selection equations may be performed
modulo 2^k directly. The reviewed short joint compiler therefore gives each
class a witness of at most 4k^2+2k+6 vertices. This is a modular witness-order
bound only. For general m the recalled compiler yields O(m^2) vertices.

## 4. Why this does not close the exact-value problem

The result removes every obstruction depending only on one fixed congruence
modulus, even within these two highly restricted source-pair families. It does
not construct a graph with a prescribed exact count, an interval, or a cofinite
spectrum. It does not imply that the pair (g(L),L) can be selected independently
over the integers. They are independent only as residues.

The small bounds make this distinction concrete: F2(L)>=29 and F3(L)>=31 in our
even-L families. Both families nevertheless represent 19 modulo every positive
modulus. The numerical integer 19 itself is absent from these families.

A nearby family does not automatically have the same congruence scope. The C2 N
graft with the same fixed v and even L has 3g(L)+6L-2, always 1 modulo 3. The C3
Type-A graft has 5g(L)+15L-6, always 4 modulo 5. The two selected rows in (1), and
their coprime coefficients, are load-bearing.

## 5. Reproduction and provenance

Run `python3 verify.py` and `python3 -O verify.py`. The standard-library verifier
checks the source responses and both actual graft graphs by independent subset
path/cycle counting for a stated finite L universe; verifies the joint compiler
by direct bit transitions for every admissible residue in its stated finite
modulus universe; and checks the two nearby modular obstructions. Those tests
supplement the proof and are not a finite-to-infinite inference.

The dependencies are frozen reviewed snapshots in dependencies/: BILINEAR_C2.md
for the M scalar/cycle law; C3_TWO_SOURCE.md for Type B; and
JOINT_MODULAR_GADGET.md for the joint gadget theorem and short dyadic compiler.
Their hashes are recorded by the verifier. The existing all-modulus C1 result
already has broader scope than separate odd/dyadic C1 arguments. The new
observation here is its direct fixed-source-pair transfer to C2 and C3, not a new
C1 theorem or a claim of literature priority.
