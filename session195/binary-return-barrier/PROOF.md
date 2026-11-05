# An arbitrary-gadget obstruction to binary-ready reset induction

Session 195, 2026-11-05. **PROVED all-order structural obstruction for a specified induction architecture, with finite controls.** This is not a counterexample to the enlarged tournament grammar or to the `{7,21}` spectrum conjecture. There is no new finite target census here.

## 1. A natural exact doubling interface

Use the physical four-coordinate state `(h,e,R,d)` from the companion combined-growth packet, with `d=H(T-r)`. Its operations are

    K(h,e,R,d)=(h+R-e,d,2R-d,d),
    D(h,e,R,d)=(h,0,2R,d).

Call a physical rooted state **binary-ready with sign σ** if

    e=d, R=h+d+σ, σ in {-1,+1}.                    (1)

A K birth then gives

    h'=2h+σ, e'=d, R'=h'+d+σ, d'=d.              (2)

Thus K implements an actual repeated binary digit step and returns the same marked interface. A plausible induction would use a rooted gadget to refresh or change this interface and then read another run of binary digits. We prove that no rooted gadget of any order can do this after the deletion count reaches three, even with arbitrarily long K/D padding following the gadget.

## 2. Structural bounds for endpoint-complete roots

**Lemma 1.** If a physical rooted tournament has e=d, then either

    (h,e,R,d)=(d,d,d,d),

with the root a universal sink, or `h>=3d`.

**Proof.** Since e=d>=1, the root is in the last strong component C of T. Hamilton-path multiplicativity across the ordered SCCs gives `h=a H(C)`, `e=a e_r(C)`, and `d=a H(C-r)`, where a is the product of the earlier component counts.

If C is a singleton, r is a universal sink. Its rooted second path can only be the singleton r, so R=h=d.

Otherwise C is nontrivial and strong. If C-r were strong, each of its vertices could end a Hamilton path. Choose an out-neighbor of r inside C; a Hamilton path of C-r ending there cannot extend by r. Thus `e_r(C)<H(C-r)`, contradicting e=d. Consequently C-r is weak. The proved Session 194 weak-deletion ratio gives `H(C)>=3H(C-r)`, and multiplication by a gives h>=3d. QED.

The only external structural input here is [the sharp weak-deletion ratio](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/weak-deletion/STRONGER_RATIO.md), sections 1-3. The positivity of all prescribed endpoint counts in a strong tournament follows by the elementary endpoint-preserving insertion lemma in that session. No strong-order cutoff is used.

The physical cone also gives `R+e>=h+d`. Under e=d this says `R>=h`, hence

    τ:=R-h-d >= -d.                                (3)

A binary-ready root with d>=3 cannot be the sink case of Lemma 1, because a sink has τ=-d, not ±1. Thus every input considered in the main theorem satisfies h>=3d.

## 3. Arbitrary padding reduces to two short return tests

For any endpoint-complete state e=d, let τ=R-h-d. K preserves τ and sends h to `2h+τ`. For m>=1, the padding block `D^m K` gives

    h'=(1+2^m)h+2^m d+2^m τ,
    τ'=(2^m-1)h+(2^m-2)d+2^m τ.                  (4)

When h>=3d and τ>=-d, equation (4) implies

    τ'>=(3*2^m-5)d>=d.                            (5)

**Lemma 2 (no delayed padding repair).** Fix d>=3. Starting from any physical endpoint-complete state e=d with τ not equal to ±1, no K/D word can end binary-ready.

**Proof.** If the starting root is a universal sink, its entire state is `(d,d,d,d)`. Every K/D transition is integer-linear, so every descendant coordinate is divisible by d, and therefore τ cannot become ±1. Otherwise Lemma 1 gives h>=3d. K preserves τ. At the first D-block followed by K, inequality (5) gives τ'>=d. The resulting endpoint-complete state is not the sink case, so the same argument applies to every later D-block. A word ending in D has e=0 and is not ready. QED.

Now apply a root-switch gadget to a binary-ready input. Write the resulting physical state as z=(h_z,e_z,R_z,d_z). To return to e=d after a nonempty padding word, its first K must follow `D^m` for some m>=0. At that first K,

    τ_K=R_z-h_z+e_z-2d_z                 if m=0,
    τ_m=2^m R_z-h_z-2d_z                if m>=1.   (6)

For m>=2, physical inequalities `R_z>=h_z>=d_z` imply

    τ_m>=h_z>=d_z.                               (7)

If d_z>=3, neither (7) nor any later padding can return to ±1, by Lemma 2. Therefore all possible returns after arbitrary padding are decided by just three tests:

- the gadget output itself is ready;
- the output after K is ready;
- the output after DK is ready.

This is an all-length reduction, not a bounded suffix search.

## 4. Every physical gadget excludes all three tests

Let (G,b) be an arbitrary nonempty rooted physical gadget, of any order. Use its usual profile

    (g,s,L,E,Q,B,delta,d_G),

where g=H(G), s and E count Hamilton paths starting and ending at b, Q=R_b(G), L counts two-path covers whose first path is nonempty and ends at b, B is the ordered two-bin path-cover count of G-b, delta indicates a singleton gadget, and d_G=H(G-b). Its actual root-switch operator is

    h_z=gh-se+LR,
    e_z=Eh-delta e,
    R_z=Qh+BR,
    d_z=d_G h.                                    (8)

Put

    U=Q+E-g-d_G, V=B-L.

**Lemma 3 (the missing first-bin identity).** One has

    U>=0,
    V=sum_{u in N_G^+(b)} L_u(G-b).                 (9)

Thus V>=d_G if b has any out-neighbor. If b is a sink, then

    U=V=0, g=E=Q=d_G, s=delta.

**Proof.** U>=0 is the universal physical cone inequality `R+e>=h+d`, applied to G. For V, delete b from the first path of an L-cover. This injects into the ordered two-bin covers of G-b. A cover is in the image precisely when its first bin is empty or its first path ends at an in-neighbor of b. The missing covers have nonempty first path ending at a unique out-neighbor u of b, and their number is the sum in (9). Each `L_u(G-b)>=H(G-b)=d_G`, by splitting every Hamilton word immediately after u. This proves the positive bound on V. If b is a sink, every cover is in the image and the other displayed identities follow directly from the definitions. For the singleton, s=delta=1; otherwise both are zero. QED.

Substitute a binary-ready input e=d, R=h+d+σ, with σ=±1, into (8). The offset after K is

    τ_K=(U+V-d_G)h+(V+s-delta)d+Vσ.               (10)

If b has an out-neighbor, then delta=0 and Lemma 3 gives

    τ_K>=V(d-1)>=d_G(d-1)>=2                      (11)

whenever d>=3. If b is a sink, equation (10) instead becomes

    τ_K=-d_G h<=-3.                              (12)

So a gadget followed by K never returns. A direct binary-ready gadget output would remain binary-ready under K, so this also rules out every direct return without any separate profile restriction.

For the gadget followed by DK, direct substitution gives

    τ_DK=A h+(B+V+s)d+(B+V)σ,
    A=(g-E)+(L-E)+2U+2V.                         (13)

Here A>=0 because `L>=g>=E`, U>=0, and V>=0. Also B>=1. Therefore

    τ_DK>=(B+V)(d-1)+sd>=2                       (14)

for d>=3. All three short tests from section 3 are universally excluded.

**Theorem 4 (no one-gadget binary reset).** Let T be any physical binary-ready input with d>=3. Choose **any** nonempty rooted tournament G, apply its one-reversed-seam root-switch operation, and then apply any finite K/D word. The result is not binary-ready with either sign.

**Proof.** Lemma 1 gives h>=3d. The output deletion count is `d_z=d_G h>=h>=9`. Direct, K, and DK returns are excluded by (10)-(14). Longer initial D-blocks are excluded by (7), and every delayed repair after the first K is excluded by Lemma 2. A padding word with no K is either empty, already excluded, or ends at e=0. QED.

The gadget may depend on the entire input, on the desired target, or on any arithmetic side information. There is no gadget-size cutoff, finite-bank assumption, or target cutoff in this theorem. The single-root-switch restriction is essential to its stated scope; several root switches while outside the ready class are not included.

### Checks against the enlarged bank

The packet retains the 65 rooted gadgets of order at most five plus the two order-six exception-closing gadgets in `inventory67.json`. All their profiles are independently recomputed. `coefficient_certificate.json` records the 156 direct/K/DK affine forms for this bank and independently verifies the universal cone bounds; their positive lower bounds are at least 5 and negative upper bounds at most -9 on h>=3d,d>=3. These checks are controls for the structural theorem, not its all-order justification.

## 5. Why the small-deletion exception matters

The hypothesis d>=3 cannot be dropped. A rooted directed triangle has

    (h,e,R,d)=(3,1,3,1), σ=-1.

The literal sequence A,D,K gives

    (3,1,3,1) -> (5,2,6,3) -> (5,0,12,3)
              -> (17,3,21,3).

The last state is binary-ready with σ=+1 and d=3. Thus a first reset can work, but Theorem 4 blocks every later reset of the specified form. This is an independently checked physical example, not a formal state.

K alone does remain reusable. From the H17 state it supplies the entire infinite sequence `18*2^k-1`, k>=0. What fails is the proposed mechanism for changing the ready state again.

## 6. An all-order sparsity consequence for the proposed induction

Fix any finite bank of binary-ready physical seeds, each with h>=3. Consider only derivations that consist of repeated K digit steps and **resets**, where each reset uses a gadget from any fixed finite bank of nonempty rooted tournaments, followed by an arbitrary K/D word, and where the state must be binary-ready both before and after every reset.

**Corollary 5.** Such a derivation has at most one reset. The set of its distinct attainable counts up to X has size `O((1+log X)^2)`. It consequently misses infinitely many ordinary-prime atoms and infinitely many 7p atoms.

**Proof.** All deletion counts d are odd. If a ready input has d>=3, Theorem 4 forbids a reset. If d=1, any reset changes it to `d_G h>=3`, and no subsequent reset is possible. Thus there is at most one.

For no-reset paths, K gives `h_k=2^k(h_0+σ)-σ`, and h_0>=3 bounds k by O(log X). For a reset after k steps, its output deletion count is at least three. Sections 3-4 imply that a successful post-gadget padding word reaches readiness either immediately, after K, or after DK, and thereafter contains only K steps: a later D-block would destroy readiness and could not be repaired. Thus the non-doubling part belongs to a fixed finite collection of at most three short forms per gadget. After it there are l K-steps, with l=O(log X), since Hamilton counts never decrease under these induced additions and every ready count is at least three. There are O(log X) choices for each of k and l, and finitely many seeds, signs, and short forms. This proves the bound.

A set with this polylogarithmic counting bound has convergent reciprocal sum, by dyadic grouping. The sum of reciprocals of odd primes diverges; removing the finitely many excluded small primes does not change that. The same holds for numbers 7p with p an odd prime at least five and p!=7. Hence neither atom lane can be covered up to only finitely many exceptions. QED.

This corollary concerns the specified ready-state strategy, not the full 67-gadget grammar or an arbitrary growing-gadget construction. In particular it does not constrain multiple root switches while outside the ready class, differently marked interfaces, old-arc reconstruction, or a target-dependent proof that uses witnesses unrelated to one another. Those remain genuine routes for the all-order problem.

## 7. Reproduction and remaining obligation

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The verifier independently recomputes all 67 physical gadget profiles by subset DP and literal permutation/cut counting, checks the missing-first-bin identity (9) and the nonnegative coefficient decomposition (13), checks the three coefficient tests without any h,d target cutoff, and checks the d=1 physical positive control and subsequent binary K family. `coefficient_certificate.json` lists every tested affine form and its universal bound. No tournament census is extended here.

The enlarged grammar still has verified finite coverage, but this attempted binary-ready induction cannot prove its all-order coverage. Simply enlarging the rooted gadget bank cannot repair this strategy. A successful argument must use a richer returning state, more than one root switch before returning, or an essentially different descent. The full `{7,21}` spectrum remains open.
