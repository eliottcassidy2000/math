# Unit determinants: a physical sign switch and a prime-amplitude barrier

Session 196, 2026-11-05. **All-order exact constructions and restrictions.** The full Hamiltonian-path spectrum remains open. No literature-priority claim is made.

## 1. The synthesis interface

The named companion source [cloning synthesis proof](../spectral-moves/PROOF.md), Theorem 4A, constructs a tournament F_m(X,Y) from two nonempty rooted pieces identified at c and a positive transitive m-chain Z. Write A=X-c and B=Y-c. Its cross arcs are

    A -> Z -> B, A -> B, Z -> c,

and all arcs within X and Y are retained. Put

    h=H(X), e=ends_c(X), d=H(X-c), R=R_c(X),
    U=R-d, sigma=R-h-d+e,
    g=H(Y), t=H(Y-c), L=L_c(Y).

Here R counts two-path covers whose second path starts at c, and L counts two-path covers whose first path ends at c; the other path may be empty. Then, for every m>=1,

    H(F_m(X,Y))=U L 2^m + d g-sigma t.               (1)

Every positive binary transitive-cloning amplifier arises this way, after taking a converse when needed. Its amplitude and offset are

    alpha=UL>0, beta=dg-sigma t.                     (2)

The formula is only for positive clone multiplicities. It does not interpret m=0 as deletion.

We use the physical inequalities `e<=d<=h<=R`, `L>=g>=1`, and the reviewed S195 zero-slack theorem: sigma>=0, and sigma=0 exactly when c has at most one out-neighbor in X. Also sigma is even. Indeed the singleton root-switch gives a physical tournament of count h+R-e; tournament Hamilton counts are odd, so R-e is even. Both h,d are odd, and consequently R+e-h-d is even.

## 2. A fixed physical compiler switches either unit sign

Suppose (T,r) is any positive binary amplifier with offset b in {+1,-1}. Thus for positive m,

    H(T[r:m])=a 2^m+b, a>0,

where T[r:m] replaces r by a transitive m-chain. Let H=H(T). Make Y=T[r:2], and mark either of its two clones as c. Deleting the marked clone returns a tournament isomorphic to T. Hence

    t=H, g=2H-b.                                    (3)

Take the fixed transitive three-vertex X with rows `[6,4,0]`, rooted at its source c=0. Its profile is

    (h,e,R,d)=(1,0,4,1), U=3, sigma=2.

Substituting (3) into (1) gives

    H(F_m(X,Y))=3L_c(Y) 2^m-b.                       (4)

**Theorem 1 (unit-sign switch).** Every actual marked unit-offset binary amplifier admits an explicit induced extension to an amplifier with the opposite offset. At m=1 it adds four vertices, changes the count to

    N=6L_c(Y)-b,

and gives a new marked vertex z whose repeated cloning performs `N -> 2N+b`. All old vertices and arcs of T are retained after choosing one of the two c-clones. The operation may be repeated indefinitely, switching the offset sign each time.

**Proof.** Equation (4) is an identity on actual graphs. Its positive amplitude 3L means the new z is a positive binary amplifier, with offset -b. The construction adds one clone to T, two vertices of X-c, and one z, so the total order increases by four. Choosing one clone in Y recovers the original induced T. Reapply the construction to the new marked amplifier to iterate it. QED.

If T is strong, every output is strong. Cloning a vertex of a strong tournament preserves strongness. In F_1, the new vertices a_1,a_2,z and c have the directed cycle

    c -> a_1 -> a_2 -> z -> c.

This strong subgraph intersects the strong Y at c and includes every newly added vertex. The whole tournament is therefore strong. Larger positive clone multiplicities preserve strongness too.

This is **not** a literal opposite binary step at the old count H. Since L>=g=2H-b,

    N>=12H-7b.                                      (5)

The opposite operation is available at the new count N, after the template has changed. It is also not two fixed commuting clone buttons: the marked root is moved to z and new non-clone cross arcs are introduced. The companion commutation and finite-template barriers therefore do not apply to this compiler.

Every switched amplitude is divisible by three and is composite: a positive amplifier contains a directed triangle, so H>=3, g>=5 and L>=5. This simple compiler supplies sign control, not arbitrary amplitudes or all-target coverage.

## 3. Complete slack arithmetic when d=1

Assume H(X-c)=1. Then X-c is transitive: a nontransitive tournament contains a directed triangle, whose three Hamilton paths inject into those of the full tournament by deterministic insertion. Order the transitive chain A source-first, and let its incidence word have bit 1 exactly where c beats the corresponding chain vertex.

Let q be the nonnegative integer represented by that binary word, allowing leading zeros, and let f(q)=H(X). The empty word has q=0 and f(0)=1. Direct rooted-cover counting gives

    R=1+q, U=q, e=1-(q mod 2),
    sigma(q)=q+1-(q mod 2)-f(q).                    (6)

For R, the singleton second path c contributes one. A nonempty continuation first using the i-th chain vertex exists exactly when its bit is 1; the later chain vertices choose either path independently, contributing its binary place value. The endpoint count e is one precisely when the last chain vertex beats c. Empty A obeys the displayed convention as well.

Appending bit 1 adds a universal sink and preserves f. Appending bit 0 adds a vertex beating only c, so the elementary K law gives the increment q+(q mod2). Thus

    f(0)=1,
    f(2q+1)=f(q),
    f(2q)=f(q)+q+(q mod2).                          (7)

Leading zeros do not change these values. Equations (6)-(7) prove the complete slack recursion

    sigma(0)=0,
    sigma(2q)=sigma(q),
    sigma(2q+1)=q+sigma(q)+(q mod2).                (8)

This uses **all** rooted tournaments with deletion count one, not just transitive X. X itself is an arbitrary one-apex tournament.

**Theorem 2 (fixed-slack classification).** For every fixed nonnegative even s, all positive U attainable with d=1 and sigma=s are exactly

    U=v 2^j, j>=0,
    v positive odd, v<=2s+1, sigma(v)=s.             (9)

The test sigma(v)=s uses (8). In particular,

    sigma=0: U in {1,2,4,8,...},
    sigma=2: U in {3,5} * {1,2,4,8,...}.             (10)

The additional U=0 occurs only at sigma=0 and gives no positive amplifier.

**Proof.** Strip powers of two from U using the first recursion in (8), leaving an odd v=2q+1. If sigma(v)=s, then `q+sigma(q)+(q mod2)=s`, so q<=s and v<=2s+1. Conversely every displayed v and every j gives that slack by (8), and its binary word is an actual incidence assignment on a transitive chain. For s=0 the odd core is v=1. For s=2, q can only be 1 or2, giving v=3 or5; q=0 gives zero and q>=3 gives at least four. QED.

Thus fixed slack and deletion one produce only finitely many dyadic amplitude rays, even though the actual piece order is unbounded.

There are genuine missing slack values. For example, **sigma=10 is impossible when d=1 at every order**. The odd-core bound reduces it to the eleven odd v from 1 through 21; their slack values from (8) are respectively

    0,2,2,6,4,8,8,14,8,14,12,

none equal to 10. The smaller nonnegative even values 0,2,4,6,8 all occur. This finite test is an all-order exclusion because Theorem 2 supplies the exhaustive odd-core bound.

For a concrete determinant obstruction, take the actual partner

    Y=[14,12,24,16,3], root c=4,
    (g,t,L)=(11,1,13).

With d=1, offset +1 would require `11-sigma=1`, hence the impossible sigma=10. No rooted X of deletion count one can supply that unit-offset synthesis, however large its order. The numerical equation has its obvious coprime solution; the missing ingredient is physical realization of the slack.

**Corollary 2.1 (finite partner banks still give only dyadic rays).** Fix finitely many rooted partner signatures `(g,t,L)` for Y. Allow arbitrary rooted X with d=1 and U>0, arbitrary m>=1, and retain only unit offsets beta in {+1,-1}. The resulting positive binary amplifiers have Hamilton counts forming a finite union of rays

    v L 2^k + beta, k>=1.                           (11a)

In particular there are O(log M) distinct resulting counts at most M.

Indeed each partner and sign force `sigma=(g-beta)/t`. If this is not a nonnegative even integer, there is no solution. Otherwise Theorem 2 leaves only finitely many odd cores v, and U=v2^j. Substituting in (1) gives (11a) with k=j+m. Every valid odd core gives its whole displayed ray. This is not merely a fixed-template statement: X is allowed arbitrary order and arbitrary incidence words, and it is the unit determinant that reduces those possibilities to finitely many arithmetic rays. It does not apply when the partner signatures themselves grow without bound.

## 4. Odd prime amplitudes force the nontrivial determinant case

The two factors of alpha in (2) have physical restrictions beyond integer factorization.

**Lemma 3.** For any physical rooted X, U=1 forces

    (h,e,R,d)=(1,0,2,1), sigma=0.                   (11)

**Proof.** Since `U=h-e+sigma`, sigma is nonnegative even, and h-e>=0, U=1 forces sigma=0 and h-e=1. By zero-slack geometry, c has at most one out-neighbor. It cannot have none, because a universal sink has h=e. Let a be its unique out-neighbor in S=X-c. The singleton-birth and terminal-removal identities give

    h=d+R_a(S)-ends_a(S), e=d-ends_a(S).

Thus `R_a(S)=h-e=1`. Since R_a(S)>=H(S)=d, d=1. Now R=d+U=2, and the odd count h<=R forces h=1 and e=0. QED.

**Theorem 4 (prime-amplitude reduction).** Suppose alpha=p is an odd prime and beta in {+1,-1} in the complete physical star synthesis. Then

    U=p, L=1, g=t=1, sigma=d-beta.                  (12)

In particular Y is transitive with c as its source. If d=1, the only possibilities are

    p=3 or5, beta=-1.                              (13)

Both occur physically. There are no positive-unit-offset odd prime amplitudes when d=1, even allowing arbitrary Y.

**Proof.** Since UL=p, either U=1 or L=1. In the first case Lemma 3 gives d=1 and sigma=0, so beta=g. A unit offset forces beta=1 and g=1. Then Y is transitive. If c has j predecessors in its source-first order, L=2^j; this cannot be the odd prime p. Thus U=1 is impossible.

We must have L=1 and U=p. Since L>=g>=1, g=1, so Y is transitive. Its formula L=2^j forces j=0; deleting c leaves a transitive tournament, including the empty case, so t=1. Equation (2) gives beta=d-sigma, proving (12).

For d=1, beta=1 requires sigma=0; Theorem 2 then makes U a power of two, impossible for odd prime p. Beta=-1 requires sigma=2, so U is 3 or5 times a power of two. Primality forces p=3 or5. Conversely choose q=3 or5 in Section 3 and Y a singleton. Equation (1) gives respectively `3*2^m-1` and `5*2^m-1`, proving physical attainment. QED.

The restriction concerns the amplitude p, not a forbidden Hamiltonian count. A target may admit another sign, a smaller clone exponent with composite amplitude, or a completely different operation. In particular this is not a global obstruction for any prime Hamilton count.

For an odd prime amplitude p>=7, every unit-offset synthesis therefore requires d>=3, with nontransitive X-c. Mere coprimality of d and sigma does not provide the needed rooted pieces. The sign-switch compiler in Section 2 also cannot supply such prime amplitudes, because its amplitudes are divisible by three.

The restriction can genuinely be escaped by allowing larger deletion data. Take

    X=[16,1,3,7,14], root c=2; Y=[0].

Then `(h,e,R,d)=(9,2,16,5)`, so U=11, sigma=4 and beta=5-4=1. The first strong output has

    rows=[6,48,34,39,45,1], marked clone vertex=5,
    H(F_m)=11*2^m+1,
    counts 23,45,89,177,353,... .

This is a positive contrast to the d=1 barrier, not a proof that every larger prime amplitude has such a presentation.

## 5. Verification and scope

`verify.py` is self-contained and uses explicit exceptions, not removable assertions. It independently computes Hamilton counts and endpoint-cover responses by subset dynamic programming, checks small profiles by literal permutations and cuts, constructs the sign-switch graphs at adjacency level, and verifies unit coefficients by actually cloning the newly marked vertex. It checks the complete d=1 incidence-word recurrence over a finite word universe and the fixed-slack odd-core rule; those controls are not substituted for the proof of (8)-(10).

Run `python3 verify.py > receipt.json` and `python3 -O verify.py > receipt.optimized.json`, then compare the outputs.

Named sources are the current Session 196 transitive-cloning packet, and S195 `combined-growth/PROOF.md`, `zero-slack/PROOF.md` at commit `ee07a5aac9cdfd8d3834fa834df4fb357fcaf769`. The d=1 Hamilton recurrence is inherited one-apex arithmetic (also explicitly present in S193 `modular_sync/PROOF.md`); this packet does not claim it as new. The deductions developed here are the fixed-slack classification and finite-partner sparsity corollary, the physical unit-sign compiler, and the prime-amplitude reduction. No dirty checkout was read or changed.
