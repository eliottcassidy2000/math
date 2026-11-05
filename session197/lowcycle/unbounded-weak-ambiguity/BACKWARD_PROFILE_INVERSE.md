# A backward arithmetic inverse for the complete rooted start-one class

Session197. **Exact all-order profile membership and witness reconstruction.** The complete rooted grammar and its physical transfers are inherited from `ROOTED.md`; the necessary response cone is inherited from `ARITHMETIC_INVERSE.md`. This packet adds a numerical last-step test and a terminating backward membership algorithm. It does not prove cofinite scalar coverage.

Write a rooted profile as (h,e,R), where h counts Hamilton paths, e counts paths ending at the root, and R counts ordered path covers whose second, nonempty path starts at the root. Let M_L be the actual ordinary one-seam gadget indexed by L>=1. Its constants g,Q,B satisfy

    g(1)=Q(1)=B(1)=1;
    (g,Q,B)(2l)=(g(l),Q(l),2B(l));
    (g,Q,B)(2l−1)=(g(l)+l−(l mod2),Q(l)+B(l),2B(l))

in the last line for l>=2. Its transfer is

    h'=g(L)h−(L mod2)e+L R,
    e'=h−[L=1]e,
    R'=Q(L)h+B(L)R.                                 (1)

Starting at (1,1,1), these gadgets give precisely every strong rooted tournament with start count one, together with the singleton and the rooted two-vertex transient. Each accepted recipe below therefore constructs a literal graph, not merely a formal integer state.

## 1. The inherited cone and the exact last-step test

Every nontrivial strong start-one profile satisfies

    h<=R<=2(h−e),  e>=1.                              (2)

In particular e<=h/2. For clarity, the cone is not asserted for arbitrary roots of a C1 graph. The ordinary transfer sends every nonempty physical rooted input into this cone whenever the output is nontrivial strong; the singleton seed itself is outside its upper inequality.

The gadget inequalities used here are Q>=g, B>=L, and, for L>=2, B<=2L−2 and 2g−Q>=1, with 2g−Q>=3 for odd L>=3. They follow directly from the binary recursion. The lower cone bound follows by subtracting h' from R'. For L=1, the upper-cone slack is R−h>=0. For L>=2, using R>=h and e<=h in its upper-cone slack proves nonnegativity. The elementary R>=h inequality follows by cutting each Hamilton word immediately before its root into an ordered cover.

**Last-step theorem.** Suppose (h,e,R) is a reachable strong start-one profile with h>3. Its last ordinary gadget is determined as a singleton or nonsingleton by:

- e>h/3: it is L=1.
- e<h/3: it has L>=2.
- e=h/3 and R>h: it is L=1.
- e=h/3 and R=h: it is L=2.

No reachable profile has e>h/3 and R=h.

**Proof.** If the last gadget is a singleton, its predecessor is nontrivial strong: the only smaller transient inputs give output h<=3. Write that predecessor (a,b,c). It satisfies c<=2(a−b), so

    h=(a−b)+c<=3(a−b)=3e.

Equality forces c=2(a−b), and then R−h=b>0. Thus a singleton cannot produce the balanced equality boundary.

For L>=2 the predecessor has path count e. Since every physical predecessor has R_old>=e and e_old<=e, even L gives h>=e+L e>=3e. Odd L>=3 has g>=3 and gives h>=3e−e+3e=5e. Equality h=3e can therefore occur only at L=2 with R_old=e. Its output has R=e+2R_old=h. These cases prove the test. QED.

At h=3 the triangle profile(3,1,3) has both a singleton transient predecessor(1,0,2) and a nonsingleton seed predecessor. This is why the theorem states h>3.

The scope restriction is essential. The strong C1 graph[8,1,3,6], rooted at vertex2, has (h,e,R)=(5,1,9) and start count2. Its singleton output is(13,4,14), for which e<h/3. Thus neither arbitrary strong inputs nor arbitrary roots of C1 graphs satisfy the last-step theorem.

## 2. The forced singleton inverse

Solving (1) at L=1 gives the unique possible predecessor

    (h_old,e_old,R_old)=(R−h+e, R−h, h−e).            (3)

This inverse is used only in the singleton branch of the last-step theorem, followed by a recursive reachability test. It is not an assertion that all integer points in the cone are reachable or that the cone is inverse-closed.

For any nontrivial strong predecessor, h_new>=3h_old/2: the cone gives e_old<=h_old/2 and R_old>=h_old. Hence a necessary singleton-branch check is

    3h_old<=2h.                                       (4)

### The deletion response needs no extra coordinate in this class

For every reachable strong start-one profile, the deletion count d=H(T−root) is already determined by its three responses:

    d=R−h+e  if 3e>h, or if 3e=h and R>h;
    d=e      otherwise.                              (5a)

At either count1 base state, set d=1. The triangle also satisfies (5a), with d=1.

Indeed deleting the newest gadget root leaves its predecessor followed by the transitive gadget tail, so d equals the predecessor's h. Apply the last-step theorem and (3), or e=h_old in the nonsingleton branch. The ambiguous triangle predecessors have the same count1. Consequently different start-one graphs with the same (h,e,R) cannot have different deletion counts. This is a nonlinear, piecewise-linear readout on the stated reachable class, not a linear extension of the universal rooted interface and not a claim for arbitrary marked roots.

## 3. The bounded nonsingleton inverse

For L>=2 the predecessor path count is exactly h_old=e. From the output R equation,

    R_old=(R−Q(L)e)/B(L).                             (5)

This must be an integer at least e. Because Q>=1 and B>=L,

    2<=L<=floor(R/e)−1.                              (6)

Thus there are finitely many candidate gadgets. Also h>=3e, so every recursive predecessor has path count at most h/3.

For odd L, the predecessor endpoint response is forced:

    e_old=g(L)e+L R_old−h.                            (7)

For even L, first require h=g(L)e+L R_old. The scalar and R equations then do not retain e_old; test its possible integers0,...,e recursively. This explicit existential choice is necessary. One cannot invert an even gadget by inventing a missing endpoint response.

## 4. Complete recursive decision and graph reconstruction

Define Accept(h,e,R) on integer triples, with a recipe output, as follows.

1. At h=1 accept only (1,1,1), with the empty recipe, and (1,0,2), with recipe[1]. At h=3 accept only (3,1,3), with recipe[2].
2. For h>3 reject out-of-range states, even h, or failure of e>=1 and (2).
3. If 3e>h, or if 3e=h and R>h, take only the forced singleton predecessor (3). Reject unless (4) holds and that predecessor is accepted. Append1 to its recipe.
4. Otherwise take the nonsingleton candidates (5)-(7), including the even-L endpoint alternatives. Accept exactly when one predecessor is accepted, and append its actual L.

No graph census or strong-order estimate appears in this procedure. Memoization is optional for correctness.

**Theorem.** Accept recognizes exactly the singleton, the rooted two-vertex transient, and all strong start-one profiles. Every returned recipe constructs the requested triple by physical gadgets.

**Proof.** Every recursive call at h>3 has strictly smaller path count: at most2h/3 in the singleton branch and at most h/3 otherwise. Therefore recursion terminates. The base profiles are exhaustive. Indeed an ordinary output of count1 can only be the first singleton step; an output of count3 must come from one of the two count1 profiles, since every strong step grows h by at least3/2. Directly substituting the count1 inputs in (1) gives only(3,1,3) at count3.

For soundness, every accepted predecessor is physical by induction and (1),(3),(5)-(7) give exactly the requested output. Appending its literal gadget preserves the rooted grammar. For completeness, take an actual grammar presentation. The last-step theorem selects its correct branch; (3) is its exact singleton predecessor, or its nonsingleton index lies within (6) and its old responses appear in (5)-(7). Induction accepts that smaller physical predecessor. QED.

All nontrivial reachable h are odd, also directly by the gadget recursion and the invariant R=e modulo2. This justifies the parity rejection without an additional source assumption.

The recursion depth is logarithmic in h, apart from the two explicit small base states. Memoized target membership is finite and polynomial in the numerical target: there are O(h^3) triples within the bounded cone, and at most O(h) candidate predecessor tests per nonsingleton triple, up to logarithmic-size integer arithmetic. This is a numerical bound, not polynomial time in the input bit length.

By the separately proved C1 classification, scalar H belongs to the strong C1 spectrum exactly when at least one triple

    1<=e<=floor(H/2), H<=R<=2(H−e)

is accepted. This is an exact scalar inverse criterion with witness reconstruction. It does not supply a uniform choice of e,R for every large odd H, which remains the obstacle to a cofinite-coverage proof.

## 5. Checks and boundaries

`verify_backward_profiles.py` independently generates the bounded forward grammar through H61, then checks every81,344 integer triple with odd h<=61,0<=e<=h,1<=R<=2h against the backward algorithm. The two sets agree exactly on159 reachable profiles. Every accepted recipe is also reconstructed as a literal tournament and counted by subset path/cover DP, including its actual root-deletion count and formula (5a). Both the transient triangle boundary and the wrong-root hostile are checked explicitly. Four previously reviewed larger output triples, with H=5099741,5477461,8231473,9551177, are also inverted and their complete returned recipes replayed. These are full-profile tests, not a new scalar coverage claim.

The previously known physical profile(5,1,7) is accepted through L3. Its formal singleton inverse(3,2,4) is rejected; the last-step test correctly prevents the false physical descent in the first place. Finite tests supplement the induction and never replace it.

Run `python3 verify_backward_profiles.py` and its `-O` form. The verifier writes `backward_profile_certificate.json`. The proof dependencies are the frozen `ROOTED.md`, `ARITHMETIC_INVERSE.md`, and `C1_COMPLETE.md`; their hashes are bound in the receipt. No stronger claim about C2/C3, arbitrary marked roots, all integer cone points, or all odd counts is made.
