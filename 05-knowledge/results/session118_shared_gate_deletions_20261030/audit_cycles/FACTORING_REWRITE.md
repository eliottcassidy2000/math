# A fanout-safe factoring rewrite closes the tight-prefix example

**PROVED classical distributive rewrite with explicit DAG boundary conditions.**
This is a semantic circuit optimization, not a new automatic deletion row and
not a universal surplus theorem for GP/P2 sequences. It supplies an actual
improvement for the counterexample family rather than merely identifying its
nonminimality.

## 1. Exclusive-cone factoring theorem

Let a binary XOR tree combine r>=2 selected AND outputs

    L_i = b AND u_i,    i=1,...,r,

and t>=0 other leaf values v_1,...,v_t. Assume:

1. All internal XOR-tree gates have no users outside this tree, except that
   the tree's root may have arbitrary downstream users and may be the output.
2. Each selected AND output L_i has exactly its displayed occurrence as a
   tree leaf and has no other user or independent output reference.
3. All b,u_i,v_j wires are retained boundary values. They originate outside
   the gates to be removed. Their existing outside fanout is left unchanged.
4. The selected AND gates and the internal XOR gates are distinct. The
   subgraph is acyclic and is the indicated tree, not merely a collection of
   expressions whose shared internal values have been counted repeatedly.

Then the cone can be replaced without changing its Boolean function or any
outside signal, saving exactly r-1 gates in the displayed implementations.

**Proof.** Associativity and commutativity of XOR allow the selected products
and the other leaves to be regrouped without increasing the number of XOR
operations. Distributivity over F_2 gives

    (b u_1)+...+(b u_r)+v_1+...+v_t
       = b(u_1+...+u_r)+v_1+...+v_t.

The old cone has r selected AND gates and r+t-1 XOR gates, for 2r+t-1 gates.
The replacement has r-1 XOR gates to compute the parity of the u_i, one AND
with b, and t XOR gates to combine the other leaves. Its total is r+t gates,
a saving of r-1. Any further simplification only improves that upper bound.
This includes t=0, where the last combination stage is empty.

All deleted non-root outputs were exclusive to the old cone. Therefore no
outside user loses a required value. Redirect every old root use to the new
root. Shared input/boundary wires, including b and the u_i, remain intact;
they may have arbitrarily many users elsewhere. The new gates depend only
on these retained predecessor values, so they can be placed before the old
root's downstream users without creating a cycle. This is a DAG-preserving
rewrite, not an invalid direct-sum lower-bound argument.

**Why the exclusivity conditions matter.** If a product b u_i is also used
elsewhere, its gate cannot simply be removed. One may retain it, but the
advertised saving then requires a fresh count and can disappear. Likewise a
shared internal XOR value may need to survive for another consumer. Equality
of the root's Boolean expression alone is not enough to justify deleting all
of its syntactic constituents.

Repeated or coincident retained boundary values cause no correctness problem;
the identity still holds. The stated count is for the displayed construction,
not a claim that either circuit is minimum.

## 2. Application to the actual GP/P2 residual blocks

After one prescribed tight cycle, the six-gate residual block is

    B=yb, T=ac, H=ab,
    V=z+B, W=V+T, O=W+H.

The selected products are B and H, sharing b. The XOR tree combines leaves
z,B,T,H. Its intermediate V,W have no external users. B,H have no users
outside the tree. The root O can feed the unchanged global XOR connector.
The retained a wire is also used by T=ac, which is allowed: shared *parents*
are retained, while the product outputs being removed are exclusive.

Reassociate and factor:

    O = z+yb+ac+ab = b(y+a)+ac+z.

An explicit replacement is

    U=y+a, P=bU, T=ac, V=P+T, O=V+z.

It has five gates. The three old XOR gates are replaced by three XOR gates:
one forming y+a and two combining P,T,z. The AND count falls from three to
two. Thus the saving is exactly one in the exhibited circuits, with no hidden
cost for reassociation and no assumption of free input-coordinate changes.

For k disjoint blocks, each exclusive cone is separate and the k-1 global XOR
connector gates are retained. The original residual has 6k+(k-1)=7k-1 gates;
the rewritten circuit has 5k+(k-1)=6k-1 gates. Exactly k distinct gates are
saved. Starting from the original 11k-1-gate presentation, the prescribed k
GP/P2 cycles remove 4k gates and the additional semantic rewrite removes k
more. Relative to the strengthened potential W=s+4I+2q, this final rewrite
preserves m,q,I,a and contributes an additional potential loss of k.

This demonstrates recoverable linear surplus for this specific hostile
presentation. It does not prove that every long tight prefix exposes these
exclusive common-factor cones, nor that the resulting circuit is minimum.
That structural coverage statement remains a separate open obligation.

## 3. A further legal M/D4 step, with coordinates charged honestly

There is an additional improvement in the *affine-state representation*, which
must not be confused with the original-coordinate gate saving above. After
factoring one block, move its affine gate U=y+a into the fair system (operation
M). Its equation is

    U+y+a=0.

Now y is an unpaired parameter whose unique computational child is this
equation. Its only prior use was the gate U=y+a. Apply D4: exchange parameter y
with auxiliary U, naming the new parameter w=U. The old y becomes an auxiliary
satisfying

    y=w+a.

The fair auxiliary matrix was [1] before the exchange and remains [1] after
it: the column of old parameter y is also [1]. The now-auxiliary y appears only
in this equation and is not the output. Delete that auxiliary and equation.
The remaining empty auxiliary matrix is invertible (determinant one). All
surviving gates now compute

    P=bw, T=ac, V=P+T, O=V+z.

These four gates represent the original residual function through the explicitly
retained affine parameter map

    (w,a,b,c,z) -> (y=w+a,a,b,c,z).

This map is a bijection, its inverse being w=y+a. Relative to the six original
block inputs, the complete embedding also retains x=0:

    (w,a,b,c,z) -> (x=0,y=w+a,a,b,c,z).

It has rank five. It changes neither the safe-core dimension nor the domain:
there are no marked pairs after the GP/P2 cycle, and all five new parameters
range freely. The new w, a, b, c, z each have a computational child, so I stays
five. M preserves s; D4 reduces s by one. The additional change is exactly

    (Delta s, Delta I, Delta q, Delta a)=(-1,0,0,0).

Across k blocks the auxiliary matrices are identity matrices. Each exchange
replaces a column by the identical unit column; deleting the corresponding
row/column leaves an identity matrix of one smaller size. The affine embeddings
are block diagonal and have full column rank 5k. Thus the entire continuation
is covered by the source's legal affine preparation and role-exchange rules.

The final state has 4k block gates and k-1 unchanged XOR connector gates:
**5k-1 charged gates**, compared with 7k-1 immediately after the tight prefix.
Every one of its 5k parameters is essential: for each quadratic factor, set
its partner to one, and for each z bit simply toggle it. The essential-input
bound therefore shows this new-coordinate circuit is minimum. The combined
factoring and M/D4 continuation contributes exactly **2k** further charged-object
loss, with I,q,a unchanged, hence the same 2k loss in W=s+4I+2q.

**Coordinate-cost boundary.** In the original fixed surviving input coordinates
(y,a,b,c,z), one must compute each w=y+a, costing k XOR gates in the displayed
implementation. That gives the earlier 6k-1-gate original-coordinate circuit.
The 5k-1 count is an affine-state representation with a recorded input map, not
a free-GL wiring claim or an upper bound of 5k-1 gates in the original coordinates.
The source's elimination machinery permits that representation change; ordinary
circuit complexity does not silently give it for free.

## 4. Literal checks

The independent `check_independent.py` builds the old and replacement residual
DAGs separately. It checks their exact gate-count difference for k=1,...,25
and equality on all 32 and 1,024 surviving-input assignments for k=1,2. It also
checks the general factoring identity on all inputs for r=2,...,8, totaling
1,016 assignments. The M/D4 continuation is also checked on all 32 and 1,024
new-coordinate assignments for k=1,2, including the explicit inverse map. For
k=1,...,25 it checks embedding rank 5k, every auxiliary identity-matrix stage,
all essential-input witnesses, and the exact final inventory. No producer
implementation is imported.

Normal and optimized Python executions agree byte for byte. These tests are
finite regression certificates; distributivity, boundary preservation and the
explicit gate counts above prove the all-r and all-k statements.
