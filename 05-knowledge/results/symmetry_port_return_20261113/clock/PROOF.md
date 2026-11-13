# A wrapped neutral refill that preserves the old graph

**PROVED:** an arbitrary-module two-arc transfer, indefinite changing-strong-module refill, protected marks, rigidity, and an exact inverse for the resulting word language. **FINITE-EXACT:** the small independent controls in `receipt.json`. The arbitrary-word family remains sparse; no all-target, prime-distribution, or Collatz convergence claim is made.

## 1. Inheritance and the new interface

Write Q_T(t) for ordered t-slot path covers, with empty slots allowed, and H(T)=Q_T(1). The inherited Π frame has ordered sites

    a,v,c,u,d0,d1,x0,f0,f1

and all arcs forward except d0,d1→a and f0,f1→c. Its two-module law is

    H(Π[M,N])=60 Q_M(3) Q_N(2)+H(M)H(N),                 (1)

with M at u and N at v. See [the clock boundary](../../target_collatz_transitions_20261113/collatz/CLOCK_BOUNDARY.md). [Session309](../../session309_changing_module_returns_20261113/MIGRATING_CLOCK.md) already moves whole final ordinal factors and refills an empty v-site. Its [fixed-bank theorem](../../session309_changing_module_returns_20261113/FIXED_BANKS.md) rules out a universal supply of changed-intercept prime refills from a finite SCC repertoire.

The [Session310 neutral pair](../../session310_positive_clocks_20261113/NEUTRAL_TAILS.md) is

    A=(2,28,25,81,97,79,7),
    B=(2,80,3,7,173,143,189,15),
    Q_A(t)=a(t)=t^7+20t^5+60t^3+14t,
    Q_B(t)=t a(t).

Both A and B are strong and rigid; the rigidity certificate is in [Session311](../../session311_recursive_symmetry_20261113/WREATH_ATLAS.md). Thus B can be replaced neutrally by A▷{w} in every uniform context. This is inherited, including its loss of selective inside-port information.

The new step wraps the reservoir. Put

    W(T)=C3[T,x,y], with T→x→y→T,
    G(M,T)=Π[M,W(T)].

The outer ports a,c and wrapper ports x,y remain named. A released sink of T has one extra incidence preventing the old one-arc transfer. Reversing that incidence as well gives the reusable move below. The v-module stays strong throughout the completed cycles; its site is never vacated.

The concept board is full-profile equivalence / actual transfer sites / protected ports / strong-component growth / arithmetic entry residues. The source-to-target map retains all old vertices and old arcs, not just their count, but does not keep the count constant throughout a cycle.

## 2. The complete wrapper observation

For every tournament T, including the empty boundary with Q_empty(t)=1,

    Q_W(T)(t)=t Q_T(t+1)+t² Q_T(t)−t Q_T(t−1).          (2)

This is a special case of the inherited C3 substitution kernel [THM-3134, full-profile transform](../../../../01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md). A short direct proof clarifies the empty correction. If x,y have different colors there are t(t−1) choices, contributing Q_T(t). If they share a color, there are t choices; that color's T-subset contributes Q_subset(2)+H(subset), except that the empty subset must contribute one instead of two. Summing these three terms gives Q_T(t+1)+Q_T(t)−Q_T(t−1). The displayed identity follows. Equality for positive integer t proves the polynomial identity.

For nonempty T, write (h,p,w)=(H(T),Q_T(2),Q_T(3)). Then

    H(W(T))=p+h,
    Q_W(T)(2)=2w+4p−2h.                               (3)

For arbitrary nonempty M let m=H(M), z=Q_M(3). Combining (1) and (3),

    H(G(M,T))=120zw+(240z+m)p+(m−120z)h.                (4)

The negative coefficient is an exact inclusion correction, not a negative path count. For M=TT_j, the complete nonempty-run signature at T is

    (720·3^j+3, 600·3^j+1, 120·3^j),                  (5)

and all higher sectors vanish. Indeed the only exterior in-neighbors of T are a and y, so every path has at most three T-runs. The checker counts those finite exterior words directly for j=1,2 and independently checks (4) by literal subset path DP on six small graphs. Thus (4) does not rely on fitting a few growing modules.

## 3. The two-arc transfer and the recurring neutral step

Suppose the inner module is T▷{w}. The vertex w loses to all T vertices and to y, and beats x. As a vertex of the outer v-slot, it also beats c. Reverse exactly

    w→x and w→c.

Every other incidence already agrees with the target. The result is precisely

    G(M,T▷{w}) → G({w}▷M,T).                          (6)

This is an equality of named graphs after the two specified reversals. All vertices and all other arcs remain. Both x and c remain in place. In particular, W(T) is still strong for every nonempty T, and the transfer does not consume the wrapper.

Using Q_(T▷1)(k)=kQ_T(k), the exact count change is

    H_after−H_before=(240z−m)p−240zh.                  (7)

The dependence on p and h is essential: this is not an ordinary or minus Collatz step on H alone.

Now start with any nonempty M,T. Choose any integer d≥0. Carry out:

1. Insert d new cyclic wrappers inside T, obtaining U=W^d(T).
2. Add a fresh copy of B as a final ordinal factor inside the wrapper: U▷B.
3. Replace only that fresh B by A▷{w}. Full-profile equality makes this step count-neutral in the entire graph.
4. Perform (6), moving w into the outer left slot.

The endpoint is

    (M,T) → ({w}▷M, W^d(T)▷A).                       (8)

Each cycle adds exactly 8+2d vertices. Every old vertex and every old-old arc survives. The neutral replacement changes only arcs among the eight fresh B vertices; the transfer changes only two arcs incident to its fresh w. The named a,c,x,y and the complete old induced graph are therefore available for the next cycle. Every endpoint's right module is strong, since it is W of a nonempty module.

Any chosen mark in the initial T survives forever. During the neutral step, all rooted or multi-marked path-cover responses with marks outside the fresh B agree: each exterior itinerary sees the same ordered nonempty cover profile in the fresh module. This includes a protected internal mark in old T. It does **not** say that births leave its full rooted response unchanged, or that unmarked Q alone determines arbitrary future selective operations at that mark. For the next cycle, the literal partition and four named ports are the sufficient interface.

This differs from a fixed finite reserve B^r: a new strong module W(W^d(T)▷A), of growing order, is produced at every endpoint. It also differs from Session309's vacant-site refill: the right site is never empty, and the old graph embeds induced throughout.

## 4. All endpoints and neutral presentations are rigid

Collapse the two transitive pairs d0,d1 and f0,f1 in Π. The seven-site quotient, in the corresponding inherited order, has rows

    (110,124,56,112,97,64,4).                          (9)

It is prime and rigid. The finite certificate tests all 119 nonsingleton proper subsets for the module condition and all 5040 permutations for automorphisms; no proper module and only the identity occur.

For a prime strong quotient, every proper module in a substituted graph lies within one quotient block. To see this, the blocks met by a crossing module form a module in the quotient, hence must be all blocks. Any excluded vertex then sees both orientations to those intersections because the quotient is strong, a contradiction. Consequently these blocks are canonical maximal proper modules. Automorphisms act through quotient automorphisms and the internal automorphisms of the blocks.

Take M=TT_j, T initially A, and iterate (8). An ordinal join of rigid tournaments is rigid because its SCC order is canonical. Also W(T) is rigid whenever T is rigid and |T|>1: its three canonical cyclic blocks have sizes |T|,1,1, so no nontrivial rotation survives. Induction proves rigidity of every inner T and of W(T). The prime rigid quotient (9) then proves rigidity of the whole G. The prepared U▷B and neutral U▷A▷1 presentations are rigid for the same reason. No claim is needed about the graph between the first and second transfer flips.

Thus the protected marks are not being used to hide an unremoved automorphism. The endpoint tournaments themselves are asymmetric. A separate simultaneous neutral symmetry erasure can take place inside an unmarked uniform factor; its compatibility with a primitive inside mark is addressed in the sibling ports lane.

## 5. The noncommuting word language

Let f=Q_T and define

    (Wf)(t)=t f(t+1)+t² f(t)−t f(t−1),
    (Af)(t)=a(t)f(t).

Cycle d acts by f→A W^d f and j→j+1. Arbitrary finite words (d1,...,dr) in nonnegative integers are physical, with the full polynomial carried forward. They do not reduce to the fixed three-moment state once wrappers are allowed: evaluating the next wrapper at t requires the old value at t+1.

The operators genuinely differ in order:

    (WA−AW)f(t)
     =t[a(t+1)−a(t)]f(t+1)−t[a(t−1)−a(t)]f(t−1).     (10)

For every integer t≥1 and nonempty T this is strictly positive, since a has positive coefficients and f(k)>0 for k≥1. At t=1 the difference is exactly1181 Q_T(2). Thus a rearrangement that preserves only the multiset of wraps and A-factors loses an observable coordinate.

Addresses are recoverable from the actual graphs. Remove the final A SCC in T. If the remaining prefix is decomposable, the last d was zero. If it is strong and larger than the seed A, peel its uniquely identifiable C3 wrappers until reaching a decomposable predecessor or the seven-vertex seed. This gives d and the previous T. The canonical outer quotient and wrapper recover j and T. Different words therefore give nonisomorphic rooted and unrooted endpoint graphs. No scalar-H injectivity is asserted.

The exact generating function by added vertices is

    Σ_words z^(added vertices)=(1−z²)/(1−z²−z⁸),      (11)

because a cycle costs 8+2d. This is ordinary composition counting, the same denominator pattern as the [earlier chain-language analysis](../../h17_exact_repairs_20261112/CHAIN_LANGUAGE.md) after a change of length unit; the new content is the physical address interpretation.

### Exact finite target inverse and its limitation

Fix the initial j0≥1 and T=A. For a word of r cycles with total D=Σd_i,

    j=j0+r,  |T|=7+7r+2D.

Every tournament X has Q_X(2)≥2^|X|, by summing the positive path counts on a subset and its complement. Therefore

    H(G(TT_j,T))≥240·3^j·2^|T|
                =K0·384^r·4^D,
    K0=240·3^j0·2^7.                                 (12)

Given Y, enumerate all r,D with K0·384^r·4^D≤Y, all nonnegative compositions of D into r entries, apply the exact polynomial operators, and test (4). The empty word is tested separately. This terminates and supplies every actual construction address in this language having target Y. It assumes no unrestricted module-realizability oracle. The `select_target` routine implements precisely this bound, with no heuristic pruning; its positive controls include the empty address and one actual refill, and its nearby hostile lies outside the language.

The six compensated outputs of the current determinant-growth lane have also been tested, allowing every j0≥1 as well as every word. The bound K0≤Y leaves only14 or15 initial j0 values per target, and380 word states in total. None enters this fixed A-seed language, including384,666,081,601. This exact gate is saved in `determinant_entry.json`; reproduce with `python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/entry_probe.py`. It does not exclude a different actual starting T, which (8) explicitly allows, or any other clock construction. The two positive lanes therefore cannot yet be concatenated at these six targets without a new entry interface.

Nevertheless the language is not a complete target compiler. For s=1/4,

    Σ_words (384^r4^D)^(-s)
       =1/[1−384^(-s)/(1−4^(-s))]<∞.

Consequently at most C(Y/K0)^(1/4) words, and hence at most that many distinct counts, lie below Y, where C is the displayed finite sum. This bound does not depend on scalar injectivity. Arbitrary changing strong modules and a growing observation depth have escaped the earlier fixed five-mode ray, but this particular fixed rule still has zero-density arithmetic image.

## 6. Concrete arithmetic gates, retaining the sidecar

If d=0 always, independent parameters j,n≥1 give T=A^n and

    H_(j,n)=3^j C_n+D_n,
    C_n=120(8709^n+2·1276^n−95^n),
    D_n=1276^n+95^n.                                 (13)

Along one repeated cycle j,n both increase by one. Its five diagonal modes are26127,3828,285,1276,95; this formula is a restricted sublanguage, not the arbitrary-word state.

Every value in (13) is1 modulo5, so the fixed divisor5 of Session310's plain Π reservoir is absent. It is divisible by3 exactly when n is odd. Also7 divides it precisely when

    n≡0 (mod3), j≡3 (mod6), or
    n≡1 (mod3), j≡0 (mod6).                           (14)

These are divisibility gates, not prime or prime-cofactor assertions. They follow directly from the periods of2,3,4 modulo7 in (13).

There is no fixed prime divisor of the whole language: the seed at j0=1 has count4,021,131, and its first zero-wrap refill has count85,423,165,841; their gcd is1. This exact local fact coexists with the zero-density bound in §5. It does not imply prime outputs.

For arbitrary words there is a useful exact ternary entry automaton. On (h,p) modulo3,

    W:(h,p)→(h+p,h+p),    A:(h,p)→(−h,p),             (15)

because Q_T(3) is divisible by3 for nonempty T. The outer H is h+p modulo3. Wrapping when h+p=0 sends the pair to the permanent zero state, so every later target is divisible by3. A successful nonzero wrap followed by A produces opposite coordinates and a current target divisible by3; a subsequent zero-wrap cycle restores equal nonzero coordinates. Starting at A=(−1,1), a first positive d enters the permanent trap. This is a precise physical entry guard, not a claim about ordinary Collatz.

For any modulus q, the entire residue observer is finite. The positive odd-cycle polynomial has integral monomial coefficients, so store f(0),...,f(q−1) modulo q, with cyclic indices. Formula (2) and multiplication by a act exactly on this vector. Together with3^j modulo q, (4) gives the target residue. Taking q=7 and49 distinguishes a seven-divisible target from exact valuation one. The state remains only a residue observer: it does not replace the exact polynomial for prescribed integer targets.

### A certified infinite seven-adic phase inside the word language

There is a positive all-parameter application of that observer. Take j0=1 and the words

    (0,d,0),   d=4+336t,   t≥0.                       (16)

The endpoint is G(TT4,T_d), where

    T_d=W^d(A▷A)▷A▷A,
    Q_Td=a² W^d(a²).

Its order is49+672t. Every graph is strong and rigid by §4. From t to t+1, insert336 cyclic wrappers around the existing strong prefix W^d(A²), leaving the two terminal A blocks, left module, exterior ports and all old arcs intact. Thus this subfamily itself has an induced all-parameter growth rule.

For q=3,5,49, let v_d be the q-entry evaluation vector of W^d(a²) modulo q. The explicit finite certificate in `phase_receipt.json` verifies the exact linear identities

    v_338=v_2   for q=3,5,49.                         (17)

The transition v→Wv is deterministic and linear. Hence (17) implies v_(d+336)=v_d for every d≥2, not just the sampled indices. Multiplying the phase-four vector by a² and applying (4) at j=4 gives, respectively,

    H≡2 (mod3),   H≡3 (mod5),   H≡35 (mod49).         (18)

Consequently every graph in (16) has exact seven-adic valuation one, and its cofactor is coprime to105. The first literal construction has49 vertices and count

    1,482,747,909,335,627,020,487,025,233
      =7·211,821,129,905,089,574,355,289,319.

No primality of this cofactor, or of infinitely many cofactors, is asserted. The period336 is a proved return period of these three particular finite-state orbits, not a claimed universal period for arbitrary polynomials or all primes. The certificate retains both equal state vectors and computes the first exact count independently from the integer polynomial. Reproduce with `python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/phase_certificate.py`.

## 7. Controls and stopping boundary

Run:

```bash
python3 -B 05-knowledge/results/symmetry_port_return_20261113/clock/refill.py
```

The checker uses four small inner modules for direct slot-transform checks, six literal outer Hamilton counts, complete finite run certificates for j=1,2, the prime-quotient module/automorphism certificate, and four finite words including positive wrappers and the ternary hostile. It compares named adjacencies for every cycle, the old induced graph, the exact neutral step, transfer drift, and modular observers at3,7,49. The inherited A/B full-profile equality is a named dependency, not re-inferred from these controls.

The achieved return is structural and indefinite. Its count-neutral part is exact, but the whole cycle changes the target. The remaining task is a target-directed choice of new rule or module that overcomes (12), while retaining the protected interface. This packet does not turn any of its graph transitions into an edge of ordinary positive Collatz.
