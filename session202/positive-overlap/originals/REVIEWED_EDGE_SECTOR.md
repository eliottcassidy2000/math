# A freely prescribed edge-order sector and exact two-overlap arithmetic

Status: PROOF CANDIDATE, pending independent audit. The construction and formulas below are all-order claims; no eventual odd-spectrum threshold is proved.

The one-apex source function and its physical binary construction are inherited from Session201/central-arithmetic/PROOF.md. The contribution here is a direct marked-edge interpretation, the exact two-overlap scalar/sector laws, and an explicit separation of free sector control from the still-missing scalar extraction. No priority claim is made for the underlying one-apex family.

## 1. An explicit graph with an arbitrary positive forward sector

For u>=2 let m=bit_length(u-1). Start with an increasing chain P=(p_0,...,p_(m-1)). Add vertices r,c. Put r->c, P->r, and orient c->p_j exactly when bit m-1-j of u-1 is one. Call the resulting tournament E(u), with marked arc r->c. Set E(1) to be the edge r->c, with no private vertices.

Every u>=2 gives a strong tournament with exactly one Hamilton cycle. The first bit of u-1 is one, so r,c,p_0,...,p_(m-1),r is a spanning directed cycle. Every Hamilton cycle must use r->c because c is r's sole out-neighbor; the remaining private vertices must occur in their unique increasing order. Its order is exactly 2+bit_length(u-1).

Let f(u) count Hamilton paths with r before c, and let a(u) count paths with c before r. Then

    f(u)=u,
    a(u)=g(u)-(u mod 2),
    H(E(u))=u+g(u)-(u mod 2),                         (1)

where

    g(1)=1,
    g(2v)=g(v),
    g(2v-1)=g(v)+v-(v mod 2), v>=2.                 (2)

Thus one actual boundary-order sector is an entirely free ordinary integer. The other sector is correlated with it by (2); it is not another free parameter.

**Forward-sector count.** A path with r before c must use r,c consecutively. The private vertices split into increasing paths before r and after c. There is one choice with the latter path empty. If its first vertex is p_j, c must beat p_j; all earlier private vertices lie before r, and each later vertex chooses either side independently. Therefore

    f(u)=1+sum_(c->p_j) 2^(m-1-j)=1+(u-1)=u.

**The other sector and the recurrence.** In a path placing c before r, r must be last, again because its only out-neighbor is c. Let D(u)=H(E(u)-r). Of these deletion paths, the only ones that cannot be followed by r end at c. There is precisely one such path when the last private vertex beats c, and none otherwise. For u>=2 this indicator is u mod2, because the last binary digit is that of u-1. At u=1 the singleton c also contributes one. Hence a(u)=D(u)-(u mod2).

Appending a final private vertex with c beating it gives a universal sink in E(u)-r, so D(2u)=D(u). Appending one beating c gives two disjoint kinds of paths: old deletion paths followed by the new vertex, except the (u mod2) old paths ending at c; and paths using the new vertex immediately before c, of which the forward-sector argument gives u. Thus D(2u-1)=D(u)+u-(u mod2) for u>=2. Together with D(1)=1, this proves D=g and (1).

## 2. Two literal shared-edge operations

Keep the marked arc r->c and all old edges. Let p=u mod2 and write a=a(u).

**Move A.** Add z with z->r, c->z and P->z. The new three-vertex piece on r,c,z is a directed triangle, and it overlaps E(u) in exactly the edge r,c. Its incidences to every private vertex are explicit. The enlarged graph is E(2u), with the same marked arc, and

    (f,a) -> (2f,a+p),
    H -> 2f+a+p = 2u+g(u).                          (3)

**Move B.** Add z with z->r,c and P->z. The new three-vertex piece on r,c,z is transitive. For u>=2 the enlarged graph is E(2u-1), and

    (f,a) -> (2f-1,a+f-1),
    H -> 3f+a-2 = 3u+g(u)-(u mod2)-2.               (4)

Both graph identities follow by appending, respectively, a one or zero to the binary word u-1. The sector identities then follow from (2). Starting from E(2), these two operations produce every E(u), u>=2: recursively halve an even target and replace an odd target u>2 by (u+1)/2.

These are exact two-overlap scalar laws, not merely the overlap injection's lower bound. They are not a product of the constituent counts. For example the triangle move A has second input count3, but its output is 2u+g(u), rather than 3H(E(u)).

In fact, for u=2^k+1 with k>=1, recurrence (2) gives g(u)=u. The input and output counts under A are respectively 2u-1 and 3u, so

    H(output)/(H(input)*3)=u/(2u-1) -> 1/2.

This is a positive nonproduct composition with arbitrarily large E(u) inputs; the other input is the fixed triangle. It does not establish a universal half-product lower bound for arbitrary two-overlap gluings.

## 3. Why the free sector is not yet a free Hamilton count

Moves A and B stay inside the same one-apex family. They control f exactly but retain the arithmetic correlation a=g(f)-(f mod2). No independent choice of a is licensed.

For a concrete hostile, E(u) never has H=13. Since H(E(u))>=u, a possible witness must have u<=13. At u=1,...,13 the exact counts are

    1,3,5,5,9,9,11,9,17,15,19,15,23.

There is an unrestricted five-vertex H13 tournament, so this is a limitation of this explicit family, not an omitted value of the unrestricted spectrum. In particular sector surjectivity cannot replace a scalar-extraction argument.

## 4. A precise sufficient extraction target

A concrete constructive goal is to turn the input E(u), including its binary adjacency recipe, into actual tournaments with counts

    6u-1 and 6u+1, for every u>=2.                   (5)

This is OPEN; no operation producing (5) has been constructed here. If such two graph constructors are found, positive existence for every positive odd count other than7 and21 follows from a short exact arithmetic argument:

- (5) covers all positive odd integers coprime to3 other than1,5,7. Counts1 and5 have the edge and E(3) as witnesses.
- An ordered join multiplies Hamilton-path counts. Thus multiplying a covered odd count by any power of3 is realized by adjoining the corresponding number of triangle factors.
- The only remaining odd numbers are7*3^a. The cases a=0,1 are7,21. For a>=2, E(39) has count63, and an ordered join with a-2 triangles gives7*3^a.

The global exclusions of7 and21 are separate inherited results and are not proved here. This is an implication from the explicitly stated extraction problem, not a proof of that problem or an assertion that a bounded graft can solve it. The extraction would need to remove or control the correlated g contribution without silently replacing the input by an unrealizable formal sector vector.

## 5. Verification and scope

The accompanying standalone verifier constructs literal adjacency rows, counts both sectors and total paths independently by subset DP, checks unique Hamilton cycles, and compares the two exact graph extensions with their larger-index graphs. It also verifies the H13 hostile and the E(39) repair. The all-order identities above follow from the path decompositions, not from interpolation.

The local scalar-extraction probes are exploratory only and are not part of this proof. They found no successful constructor in their stated small tested graft classes. A failed finite probe is not a theorem excluding other two-overlap constructions.
