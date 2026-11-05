# Independent source counts modulo every modulus with a fixed unit interface

Session 198. **PROVED, physical congruence theorem.** This is an interface
addendum to the mixed-modulus source-pair result. It gives no exact-value
coverage or inverse closure of a formal source cone.

## 1. A literal converse source

Use the reviewed one-apex gadget G(L), L>=2. Put n=L-1, s=bit_length(n),
and B=2^s. Its transitive base is x_0,...,x_(s-1), in increasing order.
Its root b is beaten by x_i exactly when bit i of n is one. Write

    H(G(L))=g(L), e_b(G(L))=1, L_b(G(L))=L,
    s_b(G(L))=L mod 2.

Here L_b counts ordered spanning covers (X,Y) with X nonempty and ending
at b and Y possibly empty. The last base vertex points to b because the
leading bit of n is one.

Add a vertex a with all arcs a->G(L) except b->a, giving T(L). Its root
deletion is the transitive order a,x_0,...,x_(s-1). Every path starting at b
must next visit a, because a has no other incoming arc, and then traverses
that transitive order. Every path ending at b traverses the entire deletion
before b. Thus s_b(T)=e_b(T)=H(T-b)=1. The word
b,a,x_0,...,x_(s-1),b is its unique Hamilton cycle, so T is strong C1.

Split Hamilton paths according to their use of b->a. Those avoiding this
arc start at a and give g(L)-(L mod 2) paths. Those using it are in bijection
with the L ending-b covers of G(L), by inserting a immediately after b.
Consequently

    H(T(L))=g(L)+L-(L mod 2).                       (1)

There is an especially simple exact cover response. In an ending-b cover
of T, choose the base vertices preceding b. The empty choice contributes
one. If the largest chosen base vertex is x_i, its incidence to b must
have bit i of n equal to one, and the i+1 earlier vertices
(a,x_0,...,x_(i-1)) are independently chosen or omitted. The complementary
base vertices form Y in their unique increasing order. Therefore

    L_b(T)=1+sum_(i:bit_i(n)=1) 2^(i+1)=2L-1.       (2)

Take the literal converse V(L)=T(L)^op and retain b. Converse exchanges
starting and ending counts and covers, while preserving H and the deletion
count. It follows that the actual source V(L) has

    (D,E,U,d)=(g(L)+L-(L mod 2), 1, 2L-1, 1),
    s_b(V(L))=1, C(V(L))=1.                        (3)

The graph has s+2 vertices and a one-vertex feedback set. In particular,
(3) is an actual source profile with a fixed unit interface, not a formal
choice of four integers.

## 2. Independent residues in the fixed interface

The reviewed JOINT_MODULAR_GADGET theorem says that, for every even M>=2,
any residue u modulo M and any odd residue z modulo M occur simultaneously
as L=u and g(L)=z, for arbitrarily large actual L. Its proof uses physical
binary gadget words and is recalled in the accompanying main packet.

**Theorem.** For every even M>=2 and any odd residues D_0,U_0 modulo M,
there are infinitely many distinct strong C1 sources with

    D=D_0, U=U_0 (mod M),     E=s=d=1.

Proof. Choose odd representatives D_0,U_0 in {1,3,...,M-1}, and put

    u=(U_0+1)/2,
    z=D_0-u+(u mod 2) (mod M).                     (4)

This is ordinary integer division of the chosen odd representative plus
one, not inversion of 2 modulo an even modulus. The residue z is odd.
Choose L=u and g(L)=z modulo M using the joint compiler. Because M is even,
L and u have the same parity. Formula (3) then gives exactly the two desired
congruences. Indeed 2L-1=U_0 modulo 2M, a stronger cover congruence than
needed. Arbitrarily large compiler outputs have D>=L, so the sources have
unbounded Hamilton counts and give infinitely many distinct profiles. QED.

For an arbitrary positive modulus m and any two residue classes each
compatible with an odd integer, choose odd lifts modulo 2m and apply the
theorem with M=2m. This includes simultaneous mixed even moduli. There is
no assumption that independent constructions modulo coprime factors can
be combined without keeping the same graph.

For M=2^k, the inherited short joint compiler uses at most
4k^2+2k+1 bits in L-1. Hence these fixed-interface sources require at most

    4k^2+2k+3 vertices.                            (5)

The general construction gives O(M^2) vertices. These are modular witness
bounds, not exact-target graph-order bounds.

## 3. Interpretation and boundary

The earlier DYADIC_PROFILE_UNIVERSALITY theorem prescribed D,U modulo
powers of two by controlling the opposite cover Q(L) through leading bits.
The converse in (2) removes that requirement: U is exactly 2L-1, so the
older joint (g,L) compiler already gives the stronger arbitrary-modulus
interface theorem. No new claim about independent g,L,Q or B residues is
made. In particular the physical H and U cannot be prescribed independently
over the integers: (3) only asserts their congruence independence.

The hypothesis L>=2 matters. At L=1 the added graph is transitive, the
root-ending count of T is zero (and the root-starting count of its converse
is zero), and there is no Hamilton cycle. The
modular compiler supplies arbitrarily large L, so this boundary is avoided.

The source retains s=e=d=1 exactly. The bilinear graft's output is a
separate graph and is not asserted to retain a start-one interface.

## 4. Dependencies and verification

The proof uses only the reviewed gadget compiler in ROOTED.md and its
JOINT_MODULAR_GADGET theorem. Equations (1)-(3) have direct counting proofs
above. `check_interface.py` independently builds the literal graphs,
counts both endpoint and cover responses by subset recursions, counts
cycles, and checks the modular selection against two binary evaluators.
It also checks the excluded singleton boundary. The finite controls are
not the premise of the all-order proof.
