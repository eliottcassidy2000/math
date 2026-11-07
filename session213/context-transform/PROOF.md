# Why the complete binary cycle code has a sparse total-path spectrum

Status: PROVED elementary statements below; independent proof audit PASS; exact frozen-file review recorded in the Session213 audit packet. Finite controls are separated from the all-order proof. No literature-priority claim. The unrestricted Hamiltonian-path spectrum remains OPEN.

## 1. Scope, inheritance, and the decisive change of observable

Write H(T) for directed Hamiltonian vertex words and C(T) for directed Hamiltonian cycles modulo rotation. The input is exactly Session212's family B+(m,q,tau), B−(m,q,tau), with m>=1, 0<=q<2^m and tau in {0,1}. Its order is m+4. On chain vertices x0,...,x(m+1), all arcs point from lower to higher index. The two further vertices are a,b, with:

- a->x1, and xj->a for j!=1;
- b->x0;
- x1->b iff tau=1;
- x(i+2)->b iff bit i of q is one, for 0<=i<m;
- a->b in B+, and b->a in B−.

The inherited exact law is C(B+)=q+tau+1 and C(B−)=q+tau. Thus these cycles realize every nonnegative integer when orders vary. The new statement is that the *total path counts of this very same family* form a zero-density set, even with every leading-zero padding allowed. The proof identifies the three-way allocation that defeats the binary cycle register.

Inheritance at source snapshot `25fde962779e`:

- `session212/cycle-bridge/PROOF.md`: complete cycle code and marked binary growth. Its cycle bijection is used, not claimed anew.
- `session212/cycle-bridge/ENDPOINT_BRIDGE.md`: the marked endpoint core U=B−−a and s_x1(U)=q+tau.
- `01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`: H is nondecreasing under induced extension, with the exact insertion/orphan loss ledger.
- `05-knowledge/results/spine_grammars_cutwidth_apexes_hankel_seams_20261104.md`, section on two-apex words: the special all-one count 2^(m+2)+2*3^m−1 is already present. We recover it as a control, not a new family.
- `01-canon/theorems/THM-4890-one-apex-collision-balance.md`, inherited Objects/(A): the one-apex scalar g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2). The underlying word recurrence is older than THM-4890.

Anchor: can the complete cycle code force an odd path tail? Niche: exact one-apex deletion and path-cover response. Wildcard: count the possible outputs below X rather than seek an affine invariant. The concept board is cycle cardinality / total paths / chain allocations / path covers / induced extensions / padding. The canonical hostile is confusing a marked endpoint count with total H. The least-used sidecar is the location of a chain vertex in one of three path runs.

## 2. The ternary contribution inside the binary family

Put epsilon_i=bit_i(q) and

    V(m,q)=sum_(i=0)^(m-1) epsilon_i 3^i 2^(m-1-i).

**Theorem 1.** Both apex orientations satisfy

    H(B±(m,q,tau)) >= (1+tau)2^m + 2V(m,q).             (1)

Proof. Construct two disjoint classes of spanning paths. Chain vertices within each displayed run occur in increasing order.

Class I has form I,b,J,a,x1,K, where x0 is in J and the m bit vertices are partitioned among I,J,K. If I is empty, each bit vertex can independently lie in J or K, giving 2^m paths. If I is nonempty and its last vertex is x(i+2), its arc to b requires epsilon_i=1. Each earlier bit vertex has three choices I,J,K, and each later bit vertex has only J,K. This gives 3^i2^(m-1-i). All other boundary arcs are automatic: b enters J at x0; J is nonempty and ends at a vertex other than x1, hence beats a; a->x1; and x1 precedes K. In particular this class never uses the apex arc.

Class II has form I,a,x1,J,b,K, where x0 is in K and the bit vertices are partitioned among I,J,K. If J is empty, x1->b is needed, so this gives tau*2^m paths. If J ends at x(i+2), epsilon_i=1 is needed and the same three-choice/two-choice count gives 3^i2^(m-1-i). The other boundaries are automatic: the optional I ends at a vertex beating a, b enters K at x0, and x1 precedes J. This class also never uses the apex arc.

Their apex orders differ, so the two classes are disjoint. Every parameter choice is recovered uniquely from its path. Adding the counts proves (1). QED.

If q>0 has bit length l, write m=l+r with r>=0. The highest one bit alone gives

    H(B±) >= 2*3^(l-1)*2^r.                            (2)

This lower bound is about total paths, while the cycles are only q+tau or q+tau+1. The missing coordinate is the third run: a bit vertex below the last selected bit has three choices in these path classes, but only two in the cycle bijection.

## 3. An all-order zero-density theorem, including all padding

Let F be the set of the four-variant comb presentations above, at all orders, and let

    A(X)=#{H(T): T in F, H(T)<=X},
    alpha=log(2)/log(3).

**Theorem 2.** For every real X>=2,

    A(X) <= [8*2^(-alpha)/(1-2^(-alpha))] X^alpha
              +4 floor(log2 X)
           <16 X^alpha +4 floor(log2 X).               (3)

In particular A(X)=O(X^alpha)=o(X). The total-path spectrum of this complete cycle family omits infinitely many odd integers and cannot contain an odd tail.

Proof. Counting presentations instead of distinct outputs gives an upper bound. For q>0, put j=l-1. For each j,r there are 2^j choices of q and four choices of (tau, apex orientation). By (2), a presentation with H<=X must satisfy 2^(r+1)3^j<=X. For a fixed r with Y=X/2^(r+1)>=1, its presentation count is at most

    4 sum_(j=0)^floor(log3 Y) 2^j < 8 Y^alpha.

Summing over all r>=0 only enlarges this bound and gives the geometric series in (3). For q=0, (1) gives H>=2^m, so there are at most four presentations for each 1<=m<=floor(log2 X). This includes even the zero-cycle B− boundary, so no positive-cycle restriction has hidden a case. The integer inequalities 2^5>3^3 and 2^8>3^5 give alpha>3/5>log_2(3/2). Thus one has 2^(-alpha)<2/3, and therefore 8*2^(-alpha)/(1-2^(-alpha))<16. Finally log X=O(X^alpha) and alpha<1. QED.

The exponent is an upper bound for the number of distinct attained values, not an asserted sharp asymptotic. We do not claim every remaining odd integer is forbidden in all tournaments. Indeed the omitted values may occur in completely different tournament families.

**Corollary 3 (finite-menu extension barrier).** Fix K deterministic recipes Phi_1,...,Phi_K, each defined on every comb presentation, and suppose Phi_i(T) is a finite tournament containing T as an induced subtournament. The union of their total-path spectra has at most K times the right side of (3) values at most X, and has density zero.

Proof. Every Hamiltonian path of a tournament has a legal insertion site for any added vertex: use its first/last vertex or a 0-to-1 transition in the adjacency word. Different old paths have disjoint sets of inserted paths. Iteration gives H(Phi_i(T))>=H(T). Therefore an output count <=X can arise only from an input counted in the proof of (3), with at most K outputs per input. QED.

The recipes may depend on the entire marked presentation and may add any finite number of vertices. The essential restrictions are at most K outputs per input and retaining the input as an induced subtournament. This does NOT rule out deleting vertices or changing old arcs, unboundedly many outputs for one input, combining several independently varying codes, or an unrelated constructive grammar. It is not a theorem against cycle-to-path transformations in general.

## 4. Exact total-path readout for the sole-exit orientation

The lower bound is enough for Theorem 2, but the full B− readout shows the extra data explicitly. Let U=B−−a, let W be the induced tournament on b and the m bit vertices, and let P_W(2) count ordered pairs of possibly empty directed paths partitioning V(W). Put

    Q=P_W(2)/2,
    R=1+sum_i (1-epsilon_i)2^(m-1-i),
    e=epsilon_(m-1).

Here Q is an integer: the nonempty vertex b lies in exactly one of the two bins, and swapping bins is a fixed-point-free involution.

**Theorem 4.** For every m>=1,

    H(B−)=H(U)−(2−tau)e +(1+tau)2^m+2V+2Q−(2−tau)R.    (4)

Proof. Every path either ends at a or uses the forced successor a,x1. Paths ending at a are H(U)−e_x1(U), where

    e_x1(U)=(2−tau)e.

Indeed a path ending x1 must put every bit vertex before b, in increasing order. It then ends b,x0,x1, or, only if tau=0, has x0 before b and ends b,x1. In both cases the largest bit vertex must beat b.

For the other paths, let E be an optional increasing bit-vertex run whose last vertex beats b; S an optional increasing run whose first vertex is beaten by b; and f an arbitrary optional increasing bit run. An empty run is permitted in E,S,f. The number of disjoint allocations covering every bit vertex to (E,S,f) is Q: insert b between E,S in the bin containing b of an ordered two-path cover. The number of allocations to (S,f) is R: S is empty in one case, or its first index i has epsilon_i=0 and the later vertices have two choices.

There are four disjoint cases for paths with forced segment a,x1:

- Apex order b,a and x0 after b: 2^m+V, as Class I.
- Apex order b,a and x0 before b: Q−R. Here the pre-b bit run E must be nonempty because x0 itself loses to b. The post-b run is S and the run after x1 is f.
- Apex order a,b and x0 after b: tau*2^m+V, as Class II.
- Apex order a,b and x0 before a: Q−(1−tau)R. Here the run between x1 and b is E; its empty choice is allowed exactly when tau=1. The post-b run is S and the pre-a run is f.

The location of x0 exhausts the possibilities, since it cannot occur in a chain run after x1 but before b. Summing proves (4). QED.

As controls, when tau=1 and q=2^m−1 one obtains the inherited all-one ray

    C=2^m,  H=2*3^m+2^(m+2)−1.

When tau=1,q=0 one obtains

    C=1,  H=3*2^m+3.

Thus even the fixed cycle count one has unbounded H along the padded family. The all-one ray also rules out any uniformly affine-in-C induced extension, but Theorem 2 is stronger: it excludes an odd tail without assuming any scalar readout formula.

## 5. What deleting the sole-exit vertex actually returns

In the tau=1 slice, let k=q+1=C(B−). The recovered one-apex word of U has code 4q+2 in the convention that bit j is one when xj->b. Padding adds source/sink vertices and does not change this scalar. The inherited one-apex recurrence gives

    H(U)=g(4q+2)=g(q)+2q+2=2k+g(k−1).                  (5)

The elementary inherited bounds 1<=g(v)<=v+1 also follow directly by induction from the recurrence. Consequently

    2k+1 <= H(U) <= 3k.                                (6)

The lower equality holds at k=2^s, since g(2^s−1)=1. The upper equality holds at k=2^s+1 for s>=1, since g(2^s)=2^s+1 in that range. These are sharp infinitely often; k=1 gives H=3. This is a constructive nearby-value statement, not exact coverage: for instance 11 is absent from (5), since H<=11 forces k<=5 and those outputs are 3,5,9,9,15. Deletion escapes Corollary 3 but loses the clean cycle register and retains the nonconstant g term.

## 6. Reproduction and independent validation boundary

Run `python3 verify.py` and `python3 -O verify.py`. The self-contained checker uses explicit exceptions, no imported research implementation, and checks:

- all four variants for every q at m=1,...,7 by independent subset path DP;
- the two literal path classes at m<=4, checking cardinality, disjointness, and every arc;
- all Hamiltonian permutations at m=1,2 as a second exact counting route;
- every B− identity (4) through m=7 using separately computed two-bin cover counts;
- the deletion formula (5), the cycle formula and both special rays;
- exact presentation-count bounds over a finite grid, independently enumerating the integer budget region in (2).

The proof of zero density is the all-order allocation/counting argument above, not extrapolation from this finite range. An independent audit must review both the allocation bijections and the quantifiers of Corollary 3.
