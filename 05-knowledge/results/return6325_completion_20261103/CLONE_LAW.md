# The deletion coordinate controls marked cloning

**PROVED for all finite A below; FINITE-EXACT applications to the supplied
6325 sources.** The ordered-run substitution mechanism is inherited from
[COMMON_GRAMMAR G4](../six_output_return_20261103/COMMON_GRAMMAR.md) and
[BC1/BC7](../boundary_carry_returns_20261103/PROOFS.md). The marked deletion
formula below refines that mechanism.

Let A be a nonempty tournament with named vertex p. Add r beaten only by p.
Replace r by a transitive module of m>=1 vertices r1,...,rm, in that order,
retaining r1 as the live root. Write

    a=H(A), s=start_p(A), d=H(A-p), F=E_p,A(2), H(empty)=1.

E counts two ordered bins, the first nonempty and ending at p, the second
possibly empty. For the resulting tournament T_m,

    H(T_m)=a+F-s+(2^m-2)(F-d).                         (CL1)

**Proof.** A path has at most two visits to the clone module: every later
entry uses p, which can occur only once. One-visit paths have the inherited
newborn count a+F-s. A two-visit path starts in the clone, then traverses a
nonempty old path ending at p, then the other clone run, then the remaining
old path, possibly empty. The middle old path cannot be the singleton p,
because the clone loses to p. All other such covers are legal. Their count
is F-d. Assign a nonempty proper subset of the m clones to the first run;
each run has its unique increasing order, giving 2^m-2 choices. The
restriction/concatenation maps are inverse. QED.

The marked defect obeys the particularly simple rule

    E_r1,T_m(2)-H(T_m) = s,  m=1;
                         d,  m>=2.                   (CL2)

For an indegree-one root the inherited defect identity is q-start_r1,
which equals start_p(T_m-r1): Hamilton paths not starting at p in the
deleted graph are exactly those after a legal initial r1. For m>=2 a path
starting at p must immediately traverse all the remaining clones in their
increasing order, since no later old vertex can enter that module. It then
traverses any Hamilton path of A-p. This gives d, including A={p}.

## Actual source-return boundary

The three supplied graphs all have H=6325 at m=1. Their laws are:

|Literal graph|H(T_m)|Defect at m=1|Defect for every m>=2|
|---|---|---:|---:|
|Inherited source|8130*2^(m-1)-1805|0|153|
|Returned +4 source|10814*2^(m-1)-4489|2|27|
|Returned +2 source|8460*2^(m-1)-2135|1|3|

Each row follows from its exact parent (a,s,d,F), retained in
`clone_probe.json`; independent endpoint algorithms check m=1,2,3,4.
The +4 repair changes H along the cloned family by
2684*(2^(m-1)-1); the +2 repair changes it by 330*(2^(m-1)-1).
Thus both are count-preserving at m=1 only.

Every pair of native commands in the supplied return has the same entry
and exit at the retained root. Expanding that root by a selected child
Hamilton path in both guards preserves the literal repair relation:
inside the child the two identical complement involutions cancel; the
external exceptions agree. The permanent birth in the +2 return is uniform
on the child. Hence the repair programs lift as relation-level programs,
but CL1 proves that their scalar H-neutrality does not lift. This is a
specific failure boundary, not a failure of the helper compiler.

Intrinsic halving of an even clone block retains its pair fibres and
reduces m to m/2. Recovering a marked reset requires d as well as the
original one-clone readiness s. A controller storing only H and its one
visible readiness bit discards this distinction. The open task is an
actual operation returning these larger defects to a common input grammar;
neither scalar interpolation nor field phase supplies that operation.
