# A minimal four-coordinate grammar for growing rooted tournaments

Session 195, 2026-11-05. **PROVED exact all-order interfaces and target-bounded decision procedure; FINITE-EXACT support results.** The unrestricted Hamilton-path spectrum with exceptions `{7,21}` remains **OPEN**.

## 1. The result and the boundary

A rooted tournament can grow by any sequence of the following operations:

1. attach a fresh rooted gadget across one reversed cross arc and move the root into the gadget;
2. add a vertex whose only old out-neighbor is the root, keeping the old root;
3. add a universal sink, keeping the old root.

The exact state `(h,e,R,d)` closes under all these operations, where `h=H(T)`, `e` counts Hamilton paths ending at the root, `R` counts ordered two-path covers whose second path is nonempty and starts at the root, and `d=H(T-r)`. The empty tournament has H=1. This grammar genuinely returns to its input interface after every operation, with unbounded numbers of pieces and root changes. No endpoint-one hypothesis or fixed graph-order bound is needed.

The fourth coordinate is necessary for this operation set: the old `(h,e,R)` state does not determine a two-step continuation. More strongly, the continuation rank is exactly four, among homogeneous linear representations.

With every rooted gadget of order at most five and the singleton seed, the grammar's **complete arbitrary-depth** support through 20,000 misses precisely `7,21,67,173`. This is not a global spectrum theorem, nor a proof that the same four omissions persist above 20,000. It is an exact classification within the defined grammar at the stated target bound.

## 2. Physical operations and their exact states

All rows in the certificates are outgoing-neighbor bitmasks on zero-based named vertices.

### A: new root with one outgoing old arc

Adjoin x with `x -> r` and `u -> x` for every old u other than r; mark x. Then

    A(h,e,R,d) = (h+R-e, h-e, h+R, h).

**Proof.** Paths not using x->r are old Hamilton paths followed by x, except those ending at r, giving h-e. Paths using x->r correspond exactly to the ordered covers counted by R: put x between the two paths, the second starting at r. A Hamilton path ending at x has count h-e. For the new rooted R, the second path is either the singleton x, contributing h, or begins x,r, in bijection with an old R-cover. Deleting x leaves T, giving d'=h. QED.

### K: the same vertex birth, retain the old root

Use the same graph as A but continue marking r. Then

    K(h,e,R,d) = (h+R-e, d, 2R-d, d).

**Proof.** A Hamilton path ending at r must end x,r, preceded by any Hamilton path of T-r, giving d. This includes a singleton old T, with empty prefix. In an R-cover rooted at r, x can lie in the first path, where it is the last vertex: removing x gives any old R-cover. Or it can lie in the second path, again last, except that an old second path consisting only of r cannot be extended because x->r. The excluded old covers number d. Thus R'=R+(R-d). Deleting r leaves T-r with a universal sink x, whose H is d. QED.

### D: neutral sink, retain the old root

Adjoin a universal sink x. Then

    D(h,e,R,d) = (h, 0, 2R, d).

Every old Hamilton path extends uniquely by x. No Hamilton path ends at the old root. In each R-cover, x is appended to either path, including an empty first path. Root deletion also adjoins a sink, so d is unchanged. All three operations preserve every old arc.

### General root-switch gadget

Take any nonempty rooted gadget (G,b). Set

    g=H(G), s=s_b(G), E=e_b(G), Q=R_b(G),
    L=sum_{b in S} e_b(G[S]) H(G-S),
    B=sum_{S subset G-b} H(G[S]) H(G-b-S),
    delta=[|G|=1], d_G=H(G-b).

Form `T -> G`, reverse only `r -> b` into `b -> r`, and mark b. The Session 195 root-switch theorem gives

    h' = g h - s e + L R,
    e' = E h - delta e,
    R' = Q h + B R.

For self-containment: paths not using the reversed seam contribute gh-se; those using it split into the rooted covers R and L. An endpoint-b Hamilton path cannot use the seam; it has the ordinary Eh count except the forbidden concatenation when G is the singleton. An R'-cover not using the seam has all of T in its first path and contributes Qh. A cover using the seam has a rooted R-cover of T and an arbitrary ordered two-bin cover of G-b, giving BR. This separation includes the cases of empty bins. There cannot be two G-to-T transitions because only b->r goes backwards.

Deleting the newly marked b removes the only reverse seam. The remainder is the ordered join `T -> (G-b)`, so

    d' = d_G h.                                      (1)

Thus the whole gadget bank combines with A,K,D in **the same four coordinates**. A itself is the singleton gadget.

## 3. A nonnegative arithmetic state, not a formal signed triple

Every physical rooted tournament satisfies

    0 <= e <= d <= h <= R,   R+e >= h+d.

The first bounds follow by deleting the final root from an endpoint path, by ordinary Hamilton-path insertion, and by splitting every Hamilton word immediately before its occurrence of r. For the last bound, partition the R-covers by whether the final vertex of the first path beats r. The mergeable covers, including an empty first path, count h. Among the other covers, those with singleton second path r count d-e. Consequently `R-h >= d-e`.

Use nonnegative coordinates

    a=e, b=d-e, c=h-d, f=R+e-h-d.

They recover the original state by

    e=a, d=a+b, h=a+b+c, R=a+2b+c+f.

The three elementary operations become positive integer matrices:

    A(a,b,c,f) = (b+c, a, 2b+c+f, 0),
    K(a,b,c,f) = (a+b, 0, 2b+2c+f, f),
    D(a,b,c,f) = (0, a+b, c, 2b+c+2f).             (2)

The general gadget also has nonnegative coefficients in these coordinates. In its final f-coordinate put

    U=Q+E-g-d_G >=0, V=B-L >=0, W=s-delta >=0.

Here U is the gadget's own f-coordinate. Deleting b from an L-cover injects into the ordered two-bin covers of G-b, so V>=0; W>=0 is immediate. Then

    a'=(E-delta)a+Eb+Ec,
    b'=(d_G-E+delta)a+(d_G-E)b+(d_G-E)c,
    c'=(g-s+L-d_G)a+(g+2L-d_G)b+(g+L-d_G)c+Lf,
    f'=(U+V+W)a+(U+2V)b+(U+V)c+Vf.               (3)

All c' coefficients are nonnegative because g>=d_G and L>=g>=s. All other coefficients were justified above. These are actual graph transformations valid on every physical input, not an assumption that every point in the nonnegative cone is realizable.

## 4. Why four coordinates really are needed

Take

    C=[2,12,9,17,7].

Roots 2 and 3 both have `(h,e,R)=(13,2,18)`, but their d-values are respectively 5 and 3. Applying K gives

    root 2: (29,5,31,5),
    root 3: (29,3,33,3).

Following with A yields H=55 and H=59. Thus no deterministic rule on the old triple alone can implement all K/A continuations.

There is also a linear minimality certificate. Four reachable inputs, with operation words read left to right, are

    empty: (1,1,1,1), A: (1,0,2,1),
    AA: (3,1,3,1), AAD: (3,0,6,1).

Observe H after suffixes empty, A, AA, KA. Their continuation-by-input matrix (rows are suffixes, columns are the listed inputs) is

    [1,1,3,3]
    [1,3,5,9]
    [3,5,9,15]
    [1,5,9,19].

Its determinant is -16. Any homogeneous linear continuation representation of dimension k factors this matrix through a k-dimensional space; therefore k>=4. Equations (1)-(3) supply dimension four. This proves exact linear dimension, not any assertion about nonlinear encodings or arbitrary tournament operations.

## 5. Complete arbitrary-depth support below a target bound

This grammar admits a finite exact target test without any graph-order cutoff.

For A, K, and every root-switch gadget,

    h' >= R.

For a gadget this follows from `g>=s`, `h>=e`, and `L>=1`. For K and A it is `h'=R+(h-e)`. Also R never decreases:

    R_A'=h+R>R,
    R_K'=2R-d>=R,
    R_D'=2R>R,
    R_gadget'=Qh+BR>R,

using d<=R, Q>=1 and B>=1. A K-step with no state change is harmless because visited states are deduplicated.

**Target-bounded closure algorithm.** Start at (1,1,1,1). Retain every reachable state with R<=X. For each retained state, apply each generator; record its h' if h'<=X, even when R'>X, but continue exploring it only when R'<=X. States have `0<=e<=d<=h<=R<=X`, so finitely many can be retained.

**Completeness proof.** Remove trailing D operations from a scalar witness, since they do not change h. If a derivation ever first reaches R>X, any subsequent non-D operation has h'>X, so it cannot produce a target <=X. All earlier prefixes therefore have R<=X. The possible last non-D operation is precisely one of the recorded terminal outputs. Replacing a prefix by a previous occurrence of the same state preserves every continuation count because the transitions depend only on that state. Hence saturation is complete for arbitrarily long derivations, including redundant neutral vertices. No Busch bound, graph census, or guessed maximum depth is used.

The supplied gadget inventory represents all rooted tournaments of orders 1 through 5: all 1,099 labelled graphs and 5,405 roots yield exactly 65 distinct **eight-entry** signatures `(g,s,L,E,Q,B,delta,d_G)`. Although the first seven already distinguish these 65 entries in this small bank, the eighth is explicitly included when checking completeness. The general theorem works with any specified finite physical gadget bank.

Exact computational results:

- A,K,D alone, singleton seed: through 32,768, precisely 29 odd values are absent:
  `7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647`.
- Root-switch gadgets of order at most five alone, singleton seed: the companion packet certifies through 20,000 the seven omissions `7,21,67,163,173,247,511`.
- Combined bank plus K,D, singleton seed: through 20,000, 2,559,611 retained states and precisely four omissions `7,21,67,173`.

The combined result is **not merely a union**: 163 and 247 are in neither of the two simpler grammars but arise by interleaving their operations. Every claimed omission is within the defined grammar, not among unrestricted tournaments. The absence of 67 and 173 from this grammar is already decided by the finite target calculation at 173, regardless of how long one grows a graph.

## 6. Explicit complementary constructions and infinite families

The following rows and states are independently recomputed by subset DP, with small profiles additionally checked by literal permutation/cut counting.

An H29 seed within the gadget grammar is

    T=[54,48,2,7,12,28], r=3,
    (h,e,R,d)=(29,3,49,5).

Then K and A give

    (29,3,49,5) -> (75,5,93,5) -> (163,70,168,75).

The final H163 graph is

    [246,240,194,7,204,220,136,8], root 7.

An H13 seed within the gadget grammar is

    U=[54,36,0,7,14,28], r=3,
    (h,e,R,d)=(13,0,36,3).

Then K,K,A give

    (13,0,36,3) -> (49,3,69,3)
                 -> (115,3,135,3) -> (247,112,250,115).

The final graph is

    [502,484,448,7,462,476,392,264,8], root 8.

The packet gives full singleton-to-seed gadget descriptions as well, so neither seed is an unexplained formal state.

**All-order K-spine identity.** For any rooted physical input `(h,e,R,d)` and k>=1,

    K^k: h_k=h+(2^k-1)(R-d)+d-e,
         e_k=d, R_k=d+2^k(R-d), d_k=d.

Induction from K proves the formula. After one further A the scalar identity, also valid at k=0, is

    H(A after K^k) = h+(2^(k+1)-1)(R-d)+d-e.       (4)

Thus the H29 seed gives `88*2^k-13`, and the marked H13 seed gives `66*2^k-17`, for every k>=0. The latter includes 247 (k=2), 511=7*73 (k=3), and prime 1039 (k=4). The former includes primes 163 and 691 (k=1,3). Small primality assertions are checked by division through the integer square root. No infinitude of prime or 7p values is asserted.

These one-parameter families illustrate the mechanism; the theorem permits arbitrary interleaving of different gadgets, K, and D, so its structural complexity is not bounded by a fixed number of clocks or runs.

## 7. Hostiles: return of an interface is not ternary closure

Every generator is integer-linear on the four-coordinate state. Consequently the gcd of its coordinates divides every descendant coordinate, regardless of how gadget choices depend on the current target.

For example, two consecutive TT2-sink root-switches produce

    [28,17,26,2,8], root 3,
    (h,e,R,d)=(9,3,9,3).

Every combined-grammar descendant of this specific marked graph has H divisible by 3. It therefore cannot have H29 or H31, despite both being globally attainable and already present elsewhere in this grammar. This rules out literal per-input ternary rearming even with unlimited future gadget choices.

All operations here preserve the old graph induced. They do not erase the S194 endpoint floors, module obstruction, or induced-ancestry exclusions. Different target derivations may branch from different earlier states. The contribution is an exact, composable, minimal linear state and a decidable structural grammar with useful complementary constructions, **not** a proof that it covers every required prime and 7p atom. The residual all-order representation/coverage problem remains open.

## 8. Reproduction

Run from this folder:

    python3 verify.py --bound 1000 > receipt-1000.json
    python3 -O verify.py --bound 1000 > receipt-1000.optimized.json
    cmp receipt-1000.json receipt-1000.optimized.json
    g++ -O3 -std=c++17 census.cpp -o /tmp/s195_four_state
    /tmp/s195_four_state inventory.txt 20000 combined > receipt-20000.json
    /tmp/s195_four_state inventory.txt 32768 near > receipt-near-32768.json

The Python verifier checks every supplied physical gadget, independently regenerates the full order-five signature set, tests all elementary operations on every rooted labelled tournament through order four, checks general gadget gluing against actual graphs, verifies positive-coordinate matrices, explicit complementary witnesses, the rank minor and gcd hostile, and independently saturates the combined grammar through the selected bound. The C++ program uses the same proved exact state recurrence but an independent implementation of saturation. At the common bound 1000 both state counts and support must agree. No checks use removable Python assertions.

Named inherited sources: S194 guide at `188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec`, especially its ternary-fan, endpoint-floor, weak-parent, and module-barrier proofs; S193 prime-induction and one-arc-join packets as linked from its guide at `169792fb5eeb483e351112b5a09706e6d7bd89f0`. The Session 195 companion rooted-interface packet supplies the independently audited general root-switch theorem and 65-gadget inventory. All new files were written in a session-specific workspace, without modifying the dirty repository.
