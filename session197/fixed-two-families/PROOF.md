# Three symmetric two-integer families in the low-cycle spectrum

Session197. **Explicit all-order constructions and exact bounded inverse tests for these families.** Every output below is a strong tournament with exactly two Hamilton cycles and feedback vertex number at most two. These families are not claimed to exhaust that graph class, the full C<=2 scalar class, or the Hamilton-path spectrum.

## 1. Physical words and the arithmetic function

Begin with the edge a->b and retain these two physical ports forever. A mask S in {0,1,2,3} appends a vertex that beats precisely the selected ports and loses to every other old vertex; bit0 selects a and bit1 selects b. Every newborn is forgotten as a port. Thus deleting a,b leaves a transitive tournament.

For u>=1 let g(u) be the inherited one-apex gadget path count:

    g(1)=1;
    g(2l)=g(l);
    g(2l-1)=g(l)+l-(l mod2), l>=2.                 (1)

Write p=u mod2 and t=v mod2. When u>=2, let w be the binary digits of u-1, most significant first. When v>=2, let z be the binary digits of v-1, least significant first. The first w digit and last z digit are1. Empty words are used at parameter1.

Define three actual words:

- N(u,v), u,v>=2: map w digits1,0 to masks1,3; append mask2; append z as masks0,1.
- M(u,v), u,v>=2: the same prefix; append masks0,3; append z as masks0,1.
- P(u,v), u,v>=1: map w digits1,0 to masks2,3; append masks0,2,3; append mask 2 z_i once for each z digit, giving masks0,2.

The claimed path counts are

    N(u,v)=g(u)g(v)+(u-p)g(v)+(v-t)g(u)
            +3uv-u-v-p-t+pt;                       (2)

    M(u,v)=g(u)g(v)+2u g(v)+2v g(u)
            +6uv-u-v-u t-v p-p-t;                  (3)

    P(u,v)=g(u)g(v)+4v g(u)+4u g(v)
            +6uv-u t-v p.                          (4)

These formulas will be derived from physical response counts. In particular the parameters are not arbitrary formal profiles.

## 2. One-apex prefix responses

At a root r, use h=H, ending count e, ordered two-cover count R with the second path nonempty and starting at r, and d=H(T-r). Empty first paths are allowed and H(empty)=1. The audited retained-root births are

    O(h,e,R,d)=(h,0,2R,d),
    Z(h,e,R,d)=(h-e+R,d,2R-d,d),                    (5)

where O appends a universal sink and Z appends a vertex beating only r. The singleton one-seam J marks the newborn instead and gives

    J(h,e,R)=(h-e+R,h-e,h+R).                       (6)

For V(u), take one apex b and a transitive chain whose incidence digits [b->x_i] are w, or just b when u=1. Its root-b profile is

    (h,e,R,d)=(g(u),p,u,1),  s_b=1.                 (7)

For a nonempty chain the root begins its unique path through the entire increasing chain, because the first digit is1. Appending a digit1 is O and a digit0 is Z. Starting from the singleton and applying (5) proves (7): R is one plus the binary code, e is one minus the last digit, and the h recurrence agrees with (1). The leading digit is1, so no exceptional root-start loss occurs at a singleton Z birth.

For N and M, the prefix graph U(u) is V(u) with a final vertex a pointing only to b. In the port order a,b its word is exactly the prescribed1/3 prefix. With D=g(u), the exact responses are

    h=D+u-p;
    e_a=D-p, R_a=D+u;
    e_b=1,   R_b=2u-1;
    s_a=s_b=1.                                     (8)

These follow from (5)-(6). U is strong C1: deleting b leaves the transitive chain x_1,...,x_m,a, with b beating its first vertex and losing to its last.

## 3. The first gate

For N, mask2 appends a clone v after a: both a and v beat only b among the old source vertices, and a->v. Keeping a as root gives

    h0=D+3u-p-2, e0=u-1,
    R0=2D+4u-p-1, d0=D+u-p, s_a=2.                (9)

This is the audited adjacent-pair first-clone identity, applied to (8). Deleting a leaves a copy of U, which explains d0. The graph remains C1. Its later suffix will produce C2.

## 4. The second gate and the other actual root

Start instead from U, append the same clone v, then append z beating only the old a. Keep v as the root. Relabel the graph into the order

    ports (v,b), then x_1,...,x_m,z,a.

The fixed-port word is precisely the prefix followed by0,3. Thus this is M's actual graph, not a scalar replacement. Its profile at v, now port a, is

    h0=3D+6u-2-2p,
    e0=D+u-1,
    R0=3D+8u-1-p,
    d0=D+2u,  s_a=2, C=2.                          (10)

Here are direct cover checks for the less immediate rows. Write W for the graph after deleting v. It is the singleton one-seam extension of U at a, hence H(W)=D+2u. In W a path ending at b must finish z,a,b and has a unique preceding transitive chain; thus e_b(W)=1. A path ending at z contributes h(U)-e_a(U)=u. The old root v has exactly the two outgoing neighbors b,z, so e_v=H(W)-1-u=D+u-1.

For R_v, split the reserved path by its beginning. If it is the singleton v, the contribution is d0. If it begins v,z, it must next take a and then b unless it stops: stopping at z contributes h(U), stopping at a contributes D, and continuing through b contributes R_b(V)=u. If it begins v,b, distribute z,a among the two path slots of a cover of V reserved at b. There are four placements per cover; precisely one is forbidden, namely appending a to the reserved singleton b. The contribution is4u-1. Their sum is d0+h(U)+D+u+(4u-1)=3D+8u-1-p.

The count h0 follows from the singleton transfer at the first clone using (9), or from the audited nonadjacent-pair two-birth identity. The same marked isomorphism identifies v as the desired LAST clone, with start count2. The source snapshot contains explicit hostile examples showing why exchanging the two clone roots would be wrong.

## 5. The third gate is a previously audited transitive cell

Before P's gate, the graph is the ordered join V(u)->b. After masks0,2 it is V(u)->C3(b,z,y). The final mask3 appends r beating a,b. Delete b from the result. The remaining graph W is obtained from V by the transitive three-vertex gadget z->y->r, with the sole reverse seam r->a, and marked exit z.

The audited distinct-terminal transfer for this entry-last/exit-first cell is

    Dcell = [[1,0,4],[0,0,1],[4,-1,6]].             (11)

Consequently the profile of W at z is

    (D+4u, u, 4D-p+6u).

The full graph is W with b appended as a new root whose only outgoing neighbor is z. Applying J gives

    h0=5D+9u-p,
    e0=D+3u,
    R0=5D+10u-p,
    d0=D+4u.                                       (12)

This establishes the e and R rows at all orders without fitting them from data. The two z-start paths in W are z,y,r followed by the unique a-start V path, and z,r followed by that V path then y. Thus s_z=2. The one-seam cycle formula makes the final graph strong C2 with s_b=2. This remains valid at u=1, when V is a singleton.

## 6. Exact scalar normal form for the binary suffix

For any rooted input (h,e,R,d) with at least two vertices, append O/Z according to z, ending in Z. Encode it by v=1+sum z_i2^(i-1), so v>=2. Let Q(v),B(v) be the inherited gadget constants accompanying g(v). Iterating (5) gives

    h' = h+(g(v)-1)d+(v-1)(R-d)+(1-t)(d-e),
    e' = d,
    R' = Q(v)d+B(v)(R-d).                          (13)

This can be proved by stripping the first letter. O corresponds to v=2l-1 and Z to v=2l, and (1) together with the corresponding Q,B recurrences gives the displayed identities. The one-letter Z case is (5). At v=1 there is no suffix; the scalar formula still equals h, but e must remain e rather than being overwritten by d.

Substitution of (9), (10), and (12) into (13) gives (2), (3), and (4), respectively. In P, v=1 is allowed and uses the original gate graph. In all nonempty suffixes the root's start count remains2: a root-start path must traverse the old graph before the new tail vertex. Cutting a new Hamilton cycle at the last Z vertex gives exactly that root-start count. Hence every asserted output has C=2 and is strong, even though O steps may make intermediate graphs weak.

## 7. Symmetry is physical converse duality

Converse the graph, swap the two ports, and reverse the order of all nonport vertices. A mask transforms by

    0<->3, 1->1, 2->2.

In each family this exchanges the prefix and suffix binary words and maps its central gate to itself:2,03, and023, respectively. Thus the resulting graph is exactly the same family's construction with parameters v,u. The N gate is fixed, M's two gate vertices are exchanged, and P's gate is reversed. This proves the scalar symmetries by an actual graph isomorphism and proves that every diagonal u=v output is self-converse. No assertion that the distinguished response root is preserved is needed or made.

## 8. Exact inverse regions and limitations

For x>=2, g(x)>=1+2(x mod2). Substituting this into (2)-(3) gives

    N(u,v)>=3uv+1,
    M(u,v)>=6uv+u+v+1,                             (14)

with equality when u and v are powers of two. For P, g>=1 and p,t<=1 give

    P(u,v)>=6uv+3u+3v+1.                            (15)

Therefore membership in each family at a specified integer h has an exact finite test. By symmetry take u<=v and enumerate uv<=(h-1)/3 for N, uv<=(h-5)/6 for M, or uv<=(h-7)/6 for P, checking the exact formula. Every solution directly compiles into its physical word. These inverse regions do not rely on a graph-order cutoff, a guessed seed bank or profile saturation.

Examples are N(2,2)=13, M(2,2)=29, M(4,2)=55, P(1,24)=231 and P(8,3)=247. Each is a physical C2 output. They illustrate why different gates matter; no claim that their union covers all desired integers is made.

## Dependencies and reproduction

The basic pair/clone identities and J/Z/O rows use the independently reviewed pair-entry-bridge snapshot, PROOF.md SHA fcff7ca192c4055243aeab2344582686214aa7aa25504d1233f3b6b33465a856. The distinct-terminal cell is the audited ribbon-terminals/ribbon-sign-flip entry-last/exit-first transitive cell; its general transfer implementation is SHA9ada1cc1cfaf33030a89122dcf4ed65103aa23f5a2ad01e5dd5c288cde198a3e. Equation(11) also has direct literal controls in this packet. The g,Q,B recurrences come from the audited ROOTED.md snapshot ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e.

The included standard-library verifier checks actual graph counts, all prefix/gate response rows, actual converse adjacency maps, parameter boundaries, inverse searches, and the asserted examples. Finite controls corroborate the preceding all-order proofs and do not promote the union to a complete scalar spectrum.
