# From native Zeckendorf guards to exact seed-predecessor languages

**Status: PROVED transfer and inverse-guard identities; FINITE-EXACT new-interface controls.** Historical results are inherited at their stated scopes, not re-audited. This is an addition to the global [nine-state seed compiler](../seed_decision/PROOF.md), not a second claim of its discovery. Universal slope coverage remains OPEN.

## 1. The recovered repository advances

The relevant source is substantially richer than ordinary uniqueness of Fibonacci digits.

1. `05-knowledge/results/inverse65_zeckendorf_20261103/ZECKENDORF.md`, Z1–Z5, recovers the actual SCC prefix code `1↔a,21↔b,22↔c` with a terminal-2 flag, and proves **two compiler-model restrictions**. A prefix-correct same-source append-only compiler cannot realize both legal Fibonacci children. A final-only finite-state transducer emitting the fixed 65-symbol bank must have eventual zero-cycle length divisible by three and at least twelve. The latter comes from actual product matrices with characteristic polynomial `(X−1)(X²−L_r X+(−1)^r)`, not from a generic digit analogy.
2. The companion `INVERSE.md`, I1–I3, supplies a genuine exact inverse: preserve `(h,x,y)`, invert each nonsingular bank matrix rationally, then test integrality, the physical cone, and strict decrease of the positive norm `6h+9x+5y`. The scalar fibre alone is insufficient. Its I5 distinguishes six H=55 presentations with no 111 attachment from a different H=55 source that has such an attachment.
3. The companion `CHARGES.md`, C1–C4, classifies actual conserved signed charge lines and gives the norm-13 recurrence with **a retained consecutive pair and phase**. The backward step `(U,V)→((40U−V)/13,U)` needs integrality and arrival at the specified seed. A quadratic norm equation alone is not an orbit-membership certificate. Modulo 13^e, the actual pair has transient exactly e and period exactly 13^e; the coprime-modulus action is invertible.
4. `05-knowledge/results/zeckendorf_guard_automaton_20261003.md`, “Exact two-register digit arithmetic,” gives the native reader

       (X,Y,p) -> (Y+d, X+Y+2d, d),  with pd=0,

   and the exact charge map `Q(X)=(2Y−3X,2X−Y)`. Reduction modulo m gives at most 2m² arithmetic guard states. The explicit 2/8 example proves that value residue plus the color does not determine the next shifted residue. Its guarded-route construction preserves the actual prior route, not merely convergence somewhere.
5. `05-knowledge/results/duck_zeckendorf_20260925.md`, Z4 and Z7–Z12, supplies the full two-coordinate charge and the addition-tree carry cocycle

       b(n)=floor((n+1)/phi²), Q(n)=(n−2b(n),b(n)),
       delta(a,b)=b(a)+b(b)−b(a+b) in {−1,0,1}.

   Relaxed distinct Fibonacci representations preserve Q; overlapping supports need the correction `(-2delta,delta)`. Neutral charge at 6 is real. These XOR charges are not the user's separate red/black/blue stripe or the typed extra-unit markers.

The actual transferable object is therefore **a legal transition together with the additional register and the predecessor/phase guard**. We now give its explicit seed-machine counterpart.

## 2. New dual compiler: propagate target requirements backward, not guessed scalar predecessors

Fix exact multiplier q and cap K=q+2. In the global seed compiler let a block L have its nine positive local counts L_ab and let the old state be z_ij. Write T_L for its 9×9 nonnegative matrix, in row-vector convention. The exact destination is

    d(i,j,a,b)=(a if a!=0 else i, j if j!=0 else b).

Thus `z'=z T_L`. This matrix need not be nonsingular. Applying the 65-bank's rational inverse mechanically would therefore be unjustified. Replace it by the following **contravariant observation pullback**.

For any nonnegative column f define

    (T_L f)_ij = sum_(a,b) L_ab f_(d(i,j,a,b)).          (D1)

All operations may be capped at K. Define terminal observation columns

    f_H=e_01+e_02,  f_W=e_22,  f_U=e_02+e_12+e_22.

For a literal suffix word w=L1...Lt, retain the triple of columns

    Phi(w)=(T_L1...T_Lt f_H,
            T_L1...T_Lt f_W,
            T_L1...T_Lt f_U),                         (D2)

coordinatewise capped at K. This has 27 nonnegative coordinates. Its precise predecessor update is

    Phi(Lw)=(T_L Phi_H(w),T_L Phi_W(w),T_L Phi_U(w)).   (D3)

Given a capped prefix z, set `(h,w,u)=(z Phi_H,z Phi_W,z Phi_U)` in the capped semiring. Then this particular physical suffix completes z to a seed **iff**

    u=q, h=w+1, h<=q+1.                               (D4)

Proof: associativity gives `(z T_w) f=z(T_w f)` for each observable. The cap is a semiring homomorphism. Condition D4 keeps all observed values below K, so the parent proof's no-overflow acceptance lemma applies. This is an exact requirement on the same physical prefix and suffix, not on unrelated scalar witnesses.

**Finite inverse-language theorem.** There are at most `(q+3)^27` distinct suffix profiles. Exhaustive backward expansion from the terminal triple therefore constructs every possible target-predecessor guard in finite time. Equal profiles have equal seed acceptance on every capped prefix. Intersecting the forward nine-state graph with these backward guards gives a sound and complete bidirectional constructor. Store an actual suffix word with each profile; the emitted graph is its ordinal concatenation with the actual prefix.

This is not claimed faster in the worst case than forward search. It offers a different mathematical certificate: a set of suffix guards closed under every block pullback, with no guard accepting the initial state, certifies absence. A successful meeting supplies a literal word. Minimum forward-plus-backward depth gives least SCC-block count, and positive vertex costs give the corresponding weighted search.

The direct relational version is also exact: for a target set R of capped states,

    Pre_L(R)={z : cap(z T_L) in R}.

D1–D4 compress suffix observations rather than enumerating arbitrary R. This is the precise replacement of inverse65's “integral predecessor in a physical cone” when block maps are singular. Well-foundedness comes from finite graph distance, not an invented strictly shrinking H.

## 3. Native Fibonacci charge repairs the actual three-channel seed sum

The unit multiplier observable is not one abstract number chosen independently of endpoints. It is

    q=x+y+z,   (x,y,z)=(z_02,z_12,z_22).               (F1)

These count disjoint endpoint classes. Disjoint combinatorial classes do **not** imply disjoint Fibonacci supports of their cardinalities. The duck-note correction is therefore directly load-bearing.

Define

    kappa=delta(x,y)+delta(x+y,z)
         =b(x)+b(y)+b(z)−b(q).                        (F2)

Then

    Q(q)=Q(x)+Q(y)+Q(z)+(2kappa,−kappa).               (F3)

The inherited cocycle shows that this corrected charge is independent of how the three channels are grouped. It is not an analogy: F3 is the exact arithmetic interface between the three physical endpoint counts and the native two-register target reader.

There are only four possible three-channel carries:

    kappa in {−2,−1,0,1}.                             (F4)

Indeed, writing alpha=1/phi² and the three fractional parts as r_x,r_y,r_z,

    kappa=−floor(r_x+r_y+r_z−2alpha),

and `−1<−2alpha<r_x+r_y+r_z−2alpha<3−2alpha<3`. All four occur in the new literal small-block controls. Modulo two, only the second charge coordinate changes, by kappa mod2. The H=W+1 predicate has its own carry-corrected charge test

    Q(H)=Q(W)+Q(1)+(2delta(W,1),−delta(W,1)).           (F5)

A charge-only test without F2/F5 can wrongly reject a genuine seed. Conversely charge equality is merely a necessary condition, not sufficient for D4.

This creates an exact guard intersection: the canonical target word is read by `(X,Y,p)`; its terminal arithmetic state supplies Q(q) and q modulo any selected modulus. The candidate physical state supplies the three U channels; F2–F3 compare them to the same charge without losing carries. At accepted or still viable states U<=q these three channels are exact and below the cap, so the carry calculation is meaningful. Overflow markers are not integers with a defined Fibonacci charge and must not be fed into this test.

No new finite-state Fibonacci-addition transducer is assumed here. The native reader is finite after fixed-modulus reduction, and the fixed-q physical graph is finite; exact integer carry checks on its bounded relevant coordinates define a finite filter. For unbounded targets the second register and carry remain explicit, and no constant-size arithmetic filter is claimed.

## 4. Why this extension escapes the old no-go theorem without weakening it

The inverse65 obstruction concerns one retained physical source, a fixed 65-symbol emission bank, and prefix-correct or final-only append-only digit compilation. Our global seed compiler instead:

- reads the entire target before selecting the finite SCC alphabet A_q;
- retains the nine physical endpoint channels, not only scalar H;
- may search and backtrack, or reconstruct a wholly new source;
- uses a target-dependent finite graph whose size and alphabet grow with q;
- certifies literal final adjacency and selected vertices after reconstruction.

It therefore does not satisfy the hypotheses of the twelve-state spectral obstruction. This is a precise change of compiler model, not a dismissal of that obstruction. Nor does the quadratic norm-13 clock prove coverage: D3 retains the actual predecessor word and D4 tests its real target, just as the historical quadratic inverse required a seed and divisibility in addition to a norm.

An interesting remaining uniform obligation is now visible: find a finite description of a family of backward profiles Phi indexed by native Fibonacci registers, closed under suitable legal block pullbacks, that always meets the initial state at some multiplier c2^j. Such a **profile invariant** would contain the missing endpoint data absent from a plain digit induction. This is a precisely specified OPEN goal, not a claim that the current finite search already supplies a uniform invariant.

## 5. Controls and deliverables

`history_bridge_check.py` imports the sibling seed compiler only for its physical block counts and composition map. It independently constructs D1, the native Fibonacci reader, exact Beatty charge, and carry checks. It uses 22 transfer vectors from marked singleton, triangle, and one explicit four-vertex tournament; 507 literal prefixes; and cap8 (q6).

Normal and optimized Python runs agree:

- 33,462 exact forward-versus-pullback tests;
- 11,154 three-channel charge and native-reader interface tests;
- 8,000 addition-tree cocycle controls;
- the retained 2/8 missing-register and neutral-charge hostile examples.

These are small controls for the new bridge, not a re-audit of the inherited Zeckendorf results, a full inverse-profile census, or new universal slope coverage. Reproduce from this directory with

    python3 history_bridge_check.py
    python3 -O history_bridge_check.py

The deterministic output is `history_bridge_receipt.json`. No external packages are needed.
