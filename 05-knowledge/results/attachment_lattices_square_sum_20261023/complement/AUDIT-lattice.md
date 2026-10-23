# Independent audit of the cyclic attachment lattice packet

23 October 2026. Reviewer: complement/trijection lane, independent of the
lattice lane's derivation. **PASS** for the stated cyclic carrier. Reviewed
[lattice/README.md](../lattice/README.md) and its checker.

The inverse-image subgroup lattice has basis (h,0),(s,v), with h|m, v|n and
h|s(n/v). The two proposed congruences are necessary and sufficient. The
one-sided kernels have the stated sizes m/h and n*gcd(h,s)/(hv). The quotient
by the single shared incidence map, its Smith/minor index test, and the global
kernel count follow from the finite abelian first isomorphism theorem. A
global offset is tested against the integer column lattice, not its rational
span. Unequal orders and isolated variables are correctly retained.

The stronger any-degree elimination theorem is valid. Each fixed neighbor
value supplies a coset of a fixed subgroup of the cyclic middle group.
Generalized CRT makes pairwise overlap equivalent to total overlap. Unary
projections and **all** fill-in pairs therefore preserve the boundary relation,
and gcd of the kernel sizes is the constant number of hidden middle witnesses.
The text correctly absorbs a middle unary restriction first, retains the old
relations for owner back-substitution, and restricts the counting claim to
uniform weights with a common scalar multiplier. Arbitrary weighted functions
and noncyclic abelian domains do not inherit this shortcut.

The actual C3[R_r,R_r,R_r] cyclic action and equality-mod-three blocks check
directly; the triangle on orders 3,9,15 has image size 9 and kernel 405/9=45.
The order-21 matching example correctly prevents silently numbering nonabelian
charts as an additive cyclic group. Abstract HNF cosets do not by themselves
prove realizability by every fixed actual cell action; the packet now states
that scope explicitly.

An initial saved-receipt comparison mixed Python tuples with JSON lists. This
was a prepublication serialization error, not a mathematical mismatch. The
author normalized the record, and the repaired default checker passed 57,909
checks independently in this lane.

The additional [audit_lattice.py](audit_lattice.py) uses direct subgroup
generation by pairs, rather than HNF enumeration, to recover all 683 subgroup
normal forms for cyclic endpoint orders through 9. It directly intersects all
163,309 triples of cyclic cosets for moduli through 18, testing Helly-2 and the
gcd intersection cardinality. The three affine lines x=0, y=0, x+y=1 in F3^2
are a positive pairwise / negative joint hostile. Normal and optimized Python
both passed **173,149** checks. This independent computation supplements, rather
than replaces, the all-size CRT proof.

```text
python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/complement/audit_lattice.py
python -B -O 05-knowledge/results/attachment_lattices_square_sum_20261023/complement/audit_lattice.py
```
