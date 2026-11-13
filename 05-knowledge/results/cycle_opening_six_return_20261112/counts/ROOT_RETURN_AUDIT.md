# Independent review of the physical module return and L/T macros

**VERIFIED / PROOF AUDITED:** the claims explicitly listed below. This
review is by the full-count lane, separately from the compiler and return
search implementations. It does not promote an indefinitely repeatable
ready grammar or Chan–Pak from OPEN.

## 1. Root and incoming-port preserving module return

Reviewed [MODULE_RETURN.md](../compiler/MODULE_RETURN.md), including its
all-k-to-arbitrary-module step, complete path count, actual named map,
strongness, and stated next-iteration boundary. The proof passes.

For the tail `1^k0001`, at m=2 the inherited L coefficient is

    1 + (2^k-1) + 2^(k+3) = 9·2^k.

At m=1, D is zero because the first tail bit is one. B is nine: one
all-before assignment, seven assignments with a nonempty selection of
the three zero vertices after the root and the final one before it, and
one with all three zeros and the final one after it. The triangle has
H=3, E(2)=S(2)=3, and P(2)=12. Hence the funnel parent has

    H = 27·2^k+27,       E(2) = 108·2^k-3,

and has exactly one Hamilton path starting at its new root. The indegree-
one child therefore has full H=135·2^k+23 and defect one.

The arbitrary-module lifting is valid: each fixed uniform-module context
has a finite expansion in ordered nonempty run-cover counts K_j(M).
The functions K_j(TT_k)=j!S(k,j) are linearly independent, including K1.
The all-k identity therefore fixes the entire coefficient vector. It
preserves H(M)=u and P_M(2)=v as two separate coordinates, giving full child
H=23u+135v, defect u, and root-deletion count 27u+27v. No property of an
arbitrary M is inferred from its size alone.

`cross_audit.py` independently constructs the source, the two births, and
the funnel target. On TT1, TT2, C3, C4, and the five-vertex H9 core, it
checks all complete counts and defects and every individual comparison.
The pulled-back map changes exactly six fixed arcs and three complete
module boundaries, hence 6+3|M| arcs. The newest physical root and its
sole predecessor are fixed, as is the entire newest-root star. Internal
module arcs are also fixed. This validates realizability, not only scalar
count agreement. All five modules also pass for the sibling's distinct
[sparse module-migration target](../returns/MODULE_MIGRATION.md).

The return is a simultaneous graph replacement; neither the proof nor
this review claims that each intermediate arc flip preserves H or is a
cycle-opening operation.

## 2. Next state and arithmetic barriers

The returned state `(x,b,d)=(H-d,delta,d)` is

    (108v-4u, u, 27v+27u).

Two ordinary indegree-one births have full count `3H+2b` and new defect
`x+b`. Consequently the next target and defect are exactly

    405v+71u,       108v-3u.

For u=1, re-entry into the initial family would require
`135v+23 = 45v'+7`, impossible modulo45. Re-entry of the following target
into the returned family would require `405v+71 = 135v'+23`, impossible
modulo135. These are correctly scoped family barriers, even allowing
arbitrary integral v'. They do not exclude another ready template.

## 3. Independent L/T derivation

For any indegree-one root r put `d=H(T-r)`, `b=E_r(2)-H(T)=d-s_r`,
and `x=H(T)-d`. Let L append one vertex beaten only by r and mark it.
Let T append a and then b, both beaten only by the old r when each is
inserted, so the second vertex beats the first; mark the second vertex.
The exact closed state maps are

    L(x,b,d) = (x+2b, d-b, x+d),
    T(x,b,d) = (2x+2b, d, 2x+2b+d).

The L law follows from the full-count birth bijection. After the first
birth, a path starting at old r must begin r,a, followed by an arbitrary
Hamilton path on T-r. Thus that old-root starting count is d. Deleting
old r leaves a global source over T-r, also with count d. Its defect is
therefore zero, and the second birth gives the stated T law directly.
This checks the orientation of the two new vertices, a potential source
of a silent indegree-two error.

Writing q=2x+2b+d, the first two T steps give q(2,1,3). The eight-vertex
macros A=TTTLL and B=TTLTL then give respectively

    q(24,9,27),       q(22,12,29).

Both have full count51q. The scalar functional for the next q evaluates
to93q and97q. This proves the complete **L/T-state** diamonds ABA=BAA and
ABB=BAB. The observation scope matters: these are equalities of the
closed three-coordinate state for future L/T operations, not a claim of
equality of the entire rooted path-cover polynomial or every star query.

On four literal sources, including the two different H9 sources and the
H97 source, `cross_audit.py` verifies L, T, and both macros directly by
subset path counts. In each A/B pair exactly four arcs differ, and the
latest named root and its entire star agree. The proof is universal;
these finite checks independently validate its implementation.

For the special one-shot retarget discussed by the return-search lane,
T and L² have equal H exactly when `x+2b=2d`. A second equality would force
`x+b=0`, equivalently H(T)=s_r(T), which is impossible: if a is r's sole
predecessor, prepend a,r to any Hamilton path on T-{a,r}. This constructs
a full path not starting at r. The obstruction is correctly specific to
that retarget, not to the new larger module replacements.

## Reproduction

```
python 05-knowledge/results/cycle_opening_six_return_20261112/counts/cross_audit.py
python -O 05-knowledge/results/cycle_opening_six_return_20261112/counts/cross_audit.py
```

The script uses explicit failure checks, so optimization cannot disable
validation. `cross_audit.json` records the modules, literal states, changed
arc counts, and script hash. The primary count compiler is separately
checked by permutation enumeration and by all225 chart controls in
`audit_counts.py`; this review does not replace those checks.
