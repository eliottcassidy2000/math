# Binary growth cannot branch indefinitely under singleton append/re-mark

**Status: PROVED for the precisely specified two-operation alphabet; FINITE-EXACT controls.** The unrestricted odd Hamiltonian-count spectrum and the existence of a richer actual binary constructor remain **OPEN**. This is not a homogeneous-module theorem and not a fixed-affine-ray theorem.

## Inherited actual interface

For a nonempty tournament A and an external marked vertex x with cut S=N+(x) intersect A, use the inherited source coordinates

    h=H(A), F=H(A+x),
    q=# ordered two-path covers (P,Q) of A with P possibly empty, Q nonempty, and first(Q) in S,
    g=# Hamiltonian paths of A ending in S;
    C=(F-h)/2, L=(q-F+h-g)/2, t=h-g.

The four positive coordinates (C,L,g,t) are nonnegative integers and h=g+t>=1. This positivity and the graph operations below are inherited from Session134 constructor/PROOF.md and Session137 constructors/THEORY.md, pinned by path, commit and SHA-256 in SOURCE_PINS.json.

Append a new vertex z with A->z->x. There are two authorized choices of the marked vertex of the resulting *same full tournament*:

    S(C,L,g,t) = (2C+L+g, L, 0, g+t)      [keep x marked],
    R(C,L,g,t) = (C+L+g, 0, t, 2C+g)      [mark fresh z].

These are S0 and R1 in the inherited notation. In particular every finite word in S,R has an actual tournament reconstruction. Old vertices are not arbitrarily re-marked; there are no other graph edits or neutral steps in this alphabet.

The two laws also have a short self-contained counting proof. In the extended graph, a Hamilton path with z last contributes F-h+g, since the forbidden old paths ending x number h-g. If z is not last, its next vertex must be x. A nonempty A-suffix after x contributes q ordered covers; an empty suffix contributes h. Thus F_new=F+q+g. Keeping x gives parent height h, endpoint count g'=0, and q'=2q by assigning z to either path and deleting its forced last occurrence. Re-marking z gives parent height F, g'=h-g (old paths ending x), and q'=q+h by splitting a cover whose second path starts x according to whether that path is the singleton (x). These native laws give the displayed positive-coordinate formulas. In particular P may be empty throughout; omitting that case would invalidate the interface.

Let F(v)=2C+g+t and I(v)=2L+g-t. Direct multiplication gives

    F(Sv)=F(Rv)=2F(v)+I(v),
    I(Sv)=I(v)-2g,
    I(Rv)=2L-2C-I(v).                         (1)

A **binary history** is a word in S,R for which every step satisfies F_next=2F+1 or F_next=2F-1. Thus its current source is eligible exactly when I=+1 or -1. The operation chosen determines the next source state; at a fixed current state both operations have the SAME next count, not two independent digits.

## The finite algebraic trap is an all-length theorem

**Lemma.** For any nonnegative integer v=(C,L,g,t) with g+t>=1:

1. If v, Rv, SRv are all eligible, then F(v)=1. The complete list is
   (0,0,0,1), (0,0,1,0), (0,1,0,1).
2. If v, Rv, RRv are all eligible, then F(v) is 1 or 3. The complete list is
   (0,0,0,1), (0,0,1,0), (1,0,0,1).

**Proof.** Put d=I(v), e=I(Rv), with d,e in {-1,1}.
For part 1 let f=I(SRv) in {-1,1}. Equation (1) yields

    t=(e-f)/2,        2C+g=t-e.

Nonnegativity leaves only (e,f,t)=(-1,-1,0) or (1,-1,1). The former forces C=0,g=1, and d=2L+1 forces L=0. The latter forces C=g=0 and d=2L-1 allows L=0 or1. This is exactly the displayed list.

For part 2, direct substitution gives

    I(RRv)=-2L-g-t.

Eligibility and g+t>=1 force L=0,g+t=1. If (g,t)=(0,1), e=1-2C permits C=0 or1. If (g,t)=(1,0), e=-1-2C forces C=0. This gives exactly the second list. QED.

**Theorem (all-depth binary thinning).** Fix any nonnegative integral seed v with g+t>=1 and F(v)>=5. In every binary history of length k>=4, all of its first k-1 operations are S. The last operation can be S or R, but it has the same scalar output in either case. Therefore for each k>=4 there is at most one output, namely F(S^k v), and if this S-history is not binary there is no such length-k binary history.

**Proof.** Binary growth is strictly increasing when F>=5. By the lemma an R at any of the first k-2 positions is impossible: there would be three binary steps beginning with that R. Thus an R could only occur in the penultimate or last position. Suppose it is penultimate. Its preceding prefix consists of j=k-2>=2 copies of S. After two copies of S,

    C_2-L = 4C+2L+2g,

which is nonnegative and even. It cannot be zero for an eligible seed of F>=5: zero would force C=L=g=0, and then I=-t=+-1 forces F=t=1. Hence C_2-L>=2. Subsequent S steps preserve g=0 and increase C-L, so the state preceding the putative penultimate R has C-L>=2. But by (1), for I in {-1,1}, its R image has I'=2L-2C-I<=-3 and is not eligible for the final binary step. Contradiction. Only the final R is possible, and both choices there have equal F by (1). QED.

**Corollary (finite actual seed bank).** Starting from any finite bank of actual sources with F>=5, the set of outputs of all binary histories has O(log X) distinct elements at most X. In particular it is not cofinite among the odds and cannot supply both binary descendants indefinitely, even with changing profiles, arbitrary operation depth, or a finite-state choice rule.

There are finitely many histories of lengths 0..3. At each later length there is at most one output per seed. Since every binary step has F_next>=2F-1, a history ending below X has length at most log_2((X-1)/(F_seed-1)). This proves the count without enumerating a depth cutoff. Equivalently, the long tails are explicit dyadic rays. If the seed is eligible with I=d, its S-tail is indefinitely eligible exactly when either g=0, or g=1 and d=1; the respective outputs are

    F_k=2^k(F_0+d)-d,                 [g=0],
    F_k=2^k F_0+1, k>=1.             [g=1,d=1].

### Specified neutral moves cannot repair the language

The conclusion survives arbitrary interleavings of the two inherited H-neutral operations

    N(C,L,g,t)=(C,C+2L+g,g+t,0)      [S1],
    Z(C,L,g,t)=(0,0,0,F)            [append a universal sink and mark it].

Now allow S or R either as a binary-growth step or as an H-neutral step. For F>=5 define the binary-dead set

    D={v : I(v)>=F(v)} union {(0,0,0,t):t>=5}.

Both N and Z send every nonnegative source of F>=5 into D: I(Nv)=F+4L+2g, and Z is on the second ray. No source in D is binary-eligible. A neutral S or R is possible exactly when F+I=2(C+L+g)=0, so C=L=g=0. From that empty-cut ray, S is the identity and R gives (0,0,t,0), whose I=F=t. Thus every allowed neutral transition preserves D, and no future binary growth can escape it. Consequently N and Z can only occur after the final binary growth; they add no scalar outputs. This includes an unbounded number of these neutral moves, not a depth-limited control.

This corollary names its neutral operations exactly. It does not cover arbitrary other neutral gadgets, two or more fresh vertices per **growth** step, nonuniform fresh-to-parent edges, old-edge edits, additional marked vertices, or arbitrary legal re-marking of an older vertex. Any of these changes the operation language and requires a new proof.

## Positive actual rays, and the missing operation

The obstruction does not deny the elementary one-digit constructions. Let A be the transitive tournament on vertices 0<1<2 (arcs from smaller to larger).

* Cut S={0} gives native (h,F,q,g)=(1,5,4,0), positive state (2,0,0,1), and I=-1. Repeating S gives F_k=2^(k+2)+1.
* Cut S={0,1} gives native (1,5,6,0), positive state (2,1,0,1), and I=1. Repeating S gives F_k=6*2^k-1.

These are literal four-vertex H=5 tournaments, not formal source vectors. The verifier reconstructs the graphs after S and R and counts paths and ordered covers independently. The first ray is the inherited reversed-extreme-arc family (THM-4749), not a new construction. The second is also inherited in the h=1,L=1 exact fibre of Session137 BOUNDED_INTERIOR.md. What is new here is the all-depth two-operation obstruction, not either ray.

For a binary proof strategy, a concrete missing resource is now visible: a legal operation outside {S,R} must restore or change I=+-1 while keeping an actual sufficiently rich family. The usual one-step re-marking does not do this indefinitely. A Jacobian or inverse calculation on an affine state space does not by itself furnish this missing graph operation.

## Verification boundary

`verify_binary.py` uses only the Python standard library. It checks the matrix identities, the complete symbolic finite cases in the lemma by solving the sign equations, a finite nonnegative-state history census, exact actual witnesses and their continuations, and hostile controls. Its normal, -O, and relocated outputs must agree. The proof above, not the finite census, supports the all-length and O(log X) statements. No unrestricted negative claim or new all-integer positive realization is made.
