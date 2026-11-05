# An exact five-times window for equal-response descendants

Session 196, 2026-11-05. **PROVED all-order, architecture-specific theorem.** The unrestricted Hamiltonian-path spectrum and the full 67-gadget grammar's eventual coverage remain open.

## 1. Statement

For a rooted tournament `(T,r)`, write

    h=H(T), e=ends_r(T), R=starting_r_two_covers(T), d=H(T-r),

with an empty first path allowed in R and H(empty)=1. Suppose R=h. The physical inequalities force e=d.

Start from this *particular marked tournament*. Allow arbitrarily many stationary-root operations K and D, and root-switch one-seam attachments with **arbitrary nonempty finite rooted gadgets**, not merely the 67-piece bank. K adds a vertex beating only the old root and retains that root. D adds a universal sink and retains the root. A denotes the singleton root-switch operation; the new vertex beats only the old root and becomes the new root. Words are applied left to right.

**Theorem 1 (exact five-times window).** The complete set of descendant Hamiltonian-path counts at most 5h is

    {h,3h,5h},                          if h=d;
    {h,2h-d,3h,4h-3d,5h},              if h>d.       (1)

In the second case physical structure forces h>=3d. Thus the displayed values are in nondecreasing order, with a repeated value only when h=3d. The theorem applies at every tournament order and operation depth; finite checks are controls, not the source of its exhaustiveness.

**Corollary 1.1 (exact ternary escape).** Suppose h>=3 and h>d. The complete descendant spectrum through 3h+4 is

    {h,2h-d,3h}
      union ({4h-3d} if h-3d<=4, otherwise empty).    (2)

For delta in {2,4},

    3h+delta is a descendant count
      if and only if h-3d=delta.                    (3)

When it exists, the extra value is realized by KK (or KA). For h=d>=3 both literal children are excluded. Therefore **no physical equal-response state with h>=3 can reach both literal ternary children** under this full operation architecture.

Stationary operations really can escape the inherited moving-root barrier, but only at the precise deletion gaps in (3). This does not exclude either scalar target from other starting graphs or other operations.

## 2. Physical preliminaries

The Session 195 combined interface proves

    A(H,e,R,d)=(H+R-e,H-e,H+R,H),
    K(H,e,R,d)=(H+R-e,d,2R-d,d),
    D(H,e,R,d)=(H,0,2R,d),

and for a rooted gadget with signature `(g,s,L,E,Q,B,delta,dG)` proves

    H'=gH-se+LR, e'=EH-delta*e,
    R'=QH+BR, d'=dG*H.                              (4)

Every physical state satisfies `0<=e<=d<=H<=R`, `R+e>=H+d`. For every gadget, `L>=g>=s`, `Q>=g>=1`, `B>=1`. Consequently H never decreases, R never decreases, and every non-D operation has `H'>=R`. Terminal D operations can therefore always be discarded when deciding scalar reachability; if R exceeds a target ceiling, no subsequent non-D operation can return below it. Neutral loops do not create any unlisted scalar.

### Lemma 1: every nontrivial physical equal-response state has h>=3d

From `R=h`, the inequality `R+e>=h+d` and `e<=d` give e=d. The zero-slack identity from Session 195 states that

    R+e-h-d=0 iff outdegree_T(r)<=1.

Here its left side is zero. If r has no out-neighbor, it is a universal sink and h=d, contrary to h>d. Thus r has exactly one out-neighbor a. Put S=T-r. Removing the last root from a Hamilton path gives

    e=d-ends_a(S).

Hence ends_a(S)=0. In the linearly ordered strong components `C_1 -> ... -> C_k` of S, a lies in a nonterminal component C_j, so j<k. Indeed every vertex of a strong tournament is a Hamilton-path endpoint.

Let P(U;2) count ordered spanning two-path covers, with both slots allowed empty. Directly splitting the components of a rooted cover gives

    R_a(S)= [product_(i<j) H(C_i)] R_a(C_j)
            [product_(i>j) P(C_i;2)].               (5)

Earlier components must all lie in the first path; every later component independently splits between the two ordered paths. Since `R_a(C_j)>=H(C_j)` and `P(C_i;2)>=2H(C_i)`, equation (5) gives `R_a(S)>=2d`.

The singleton-birth identity, applied to S with new root r beating only a, gives

    h=d+R_a(S)-ends_a(S)>=3d.                       (6)

This argument covers weak T as well as strong T; no strongness assumption or external vertex-order bound is used.

For completeness, the zero-slack fact used above follows by partitioning rooted covers into mergeable covers (h of them), nonmergeable covers with second path `(r)` (d-e), and the remainder. The latter are precisely two nonempty paths of S ending and starting at distinct out-neighbors of r. For any distinct x,y, at least one of those two endpoint orders has a spanning two-cover: insert x,y around a Hamilton word of the remaining vertices at its first position beaten by either x or y. Therefore the remainder is empty exactly when r has at most one out-neighbor.

## 3. Excluding large and nontransitive gadgets below 5h

We first classify counts **strictly below 5h**, assuming h>d. Throughout this proof h,d denote the initial numbers, and put `h_1=2h-d>=5h/3`. Since H never decreases, a word reaching H>=5h cannot contribute a later count strictly below 5h. At the end we add the attainable boundary value 5h.

A nontransitive gadget has g>=3: it contains a directed triangle, whose Hamilton words inject into Hamilton words of the whole tournament by deterministic insertion. If such a gadget is the first non-D operation and no D precedes it, (4) and d<=h/3 give

    H' >= g(h-d)+Lh >= 5h.

If at least one D precedes it, R>=2h and e=0, so its count is still larger. After any first non-D operation remaining below 5h, its R is at least h_1: this holds for A,K directly and for every gadget since `R'>=H+R>=2h`. Every later nontransitive gadget has `H'>=LR>=3R>=5h`. Therefore no nontransitive gadget occurs in a word whose final count is strictly below 5h.

For a transitive m-vertex gadget, mark its vertex j in source-first order, 0<=j<m. Its exact constants are

    g=1, s=[j=0], L=2^j,
    E=[j=m-1], Q=2^(m-1-j), B=2^(m-1).

These follow by assigning the vertices before or after the root independently to the two cover slots.

If j=0, its scalar output is the same as A. For m>=2 its output has e'=0 and `R'=2^(m-1)(H+R)>=4h`. It also has H'>=h_1: if it is the first non-D step this is direct (preceding D only increases the value), while otherwise H was already at least h_1. Every subsequent non-D operation from this output has scalar at least `H'+R'`: for A,K this is exact because e'=0, and for a general gadget it follows from `g,L>=1`. Inserting D preserves e=0 and only increases R. Consequently every subsequent non-D output is at least

    h_1+4h=6h-d>5h.

A source-marked transitive gadget of size at least two is therefore terminal among non-D operations in the strict window. Replacing that terminal gadget by A preserves its scalar, so it adds no new count.

If j>=1, its scalar output is `H'=H+2^j R`. On the initial state, staying strictly below 5h forces j=1 and yields 3h. After any preceding D this expression is at least 5h. After a preceding non-D step, both H and R are at least h_1, so it is at least `3h_1>=5h`. Thus a nonsource-marked transitive gadget can only be the first non-D step, with no preceding D, and contributes exactly 3h. Its output has H=3h, R>=3h and e<=h; any next A,K or source-marked transitive operation has scalar at least `H+R-e>=5h`, and a nonsource-marked one is larger. Nontransitive gadgets have already been excluded. Inserting D cannot lower any of these bounds.

It remains to exhaust A,K,D.

## 4. Exact elementary exhaustion and the boundary

Both A and K first give h_1=2h-d, with states

    A: (h_1,h-d,2h,h),
    K: (h_1,d,h_1,d).

The possible two-letter non-D words give

    AA: (3h,h,4h-d,h_1),
    AK: (3h,h,3h,h),
    KA: (4h-3d,2h-2d,4h-2d,h_1),
    KK: (4h-3d,d,4h-3d,d).                         (7)

From AA, AK, KA and KK, the next A or K has scalar respectively

    6h-d, 5h, 6h-3d, 8h-7d,

all at least 5h, using d<=h/3. Thus three non-D elementary steps are impossible strictly below 5h. Section 3 already rules out other continuations except replacing a terminal source-marked gadget by A.

If D comes first, one D followed by A or K gives exactly 3h. Two or more initial D moves followed by a non-D operation have e=0, H=h and R>=4h, giving scalar at least 5h. After DA the next A/K count is 5h; after DK it is `7h-2d>=5h`. A D inserted after an initial A or K makes e=0 and doubles R. A subsequent A/K then has count `H+2R>=3h_1>=5h`; a larger transitive operation is no cheaper. More D moves only increase this bound. A D inserted after two non-D steps cannot improve on the third-step bounds above. Trailing D moves leave the scalar unchanged and can be discarded.

This exhausts arbitrary-depth words strictly below 5h. The empty word, A, AA and KK realize respectively

    h, 2h-d, 3h, 4h-3d.

All four are strictly below 5h. Finally DDA realizes 5h. No larger intermediate H can return below the ceiling, so this proves the second line of (1), including its boundary.

If h=d, the input is h times the singleton state `(1,1,1,1)`. Every operation is integer-linear, so each descendant is h times the state obtained by applying the same word to a singleton. The latter is physical and has an odd Hamilton count. Thus the only possible scalar values through 5h are h,3h,5h. The empty word, DA and DDA realize them. This proves the first line of (1), including h=1.

For Corollary 1.1, h>d implies h>=3d>=3. Hence `3h+4<5h`. Restrict (1) to that smaller ceiling. The value `4h-3d=3h+(h-3d)` survives exactly when h-3d<=4. The other two smaller values cannot equal either literal child, giving (2)-(3). In the diagonal case h=d>=3, neither child is an odd multiple of h (including h=3, where 11 and 13 are not multiples). QED.

## 5. Explicit controls and interpretation

The primitive strong state `(h,d,h,d)=(5,1,5,1)` has outgoing rows

    [12,9,2,4], root=2.

Its exact descendant counts through 19 are `{5,9,15,17}`. KK realizes 17, so the inherited moving-root exclusion of 17 does not survive adding stationary K. The same input still cannot reach 19 in the combined architecture, even with arbitrary-sized gadgets.

The primitive strong H49 state inherited from the Session 195 moving-root barrier is

    [122,112,67,102,108,4,32], root=5,
    (h,e,R,d)=(49,13,49,13).

Here h-3d=10. Its exact combined-architecture descendants through 151 have counts `{49,85,147}`. Neither 149 nor 151 is attainable from this particular marked input, despite its state content being one. This strengthens the inherited primitive hostile from moving-root-only operations to arbitrary root-switch gadgets interleaved with K and D.

The verifier constructs positive examples of both escape gaps from small physical tournaments. It independently counts paths and covers, verifies the displayed states, constructs KK literally, and exhausts the frozen 67-gadget bank through 5h and through 3h+4 for each enumerated nontrivial input state. Those bank checks are finite controls; the arbitrary-gadget conclusion follows from Sections 2-4.

### An infinite primitive hostile family, and a presentation contrast

For every m>=3, take a transitive chain on m+1 vertices and a root r beating only its first vertex. This is strong, since r followed by the whole chain and back to r is a directed spanning cycle. The root-deleted count is d=1, the root endpoint count is e=1, and the singleton-birth formula gives

    h=2^m+1, R=h.

Its coordinate gcd is one. Since `h-3d=2^m-2>=6`, neither `3h+2` nor `3h+4` is a descendant count, allowing any finite rooted gadgets and any number of interleaved K,D operations.

For example, H17 has such a marked representation with d=1 and cannot reach either 53 or 55 from that input. In contrast, the verified six-vertex presentation

    rows=[56,17,51,54,32,2], root=5,
    state=(17,5,17,5)

does reach 53 by KK. Likewise

    rows=[56,49,51,22,32,8], root=5,
    state=(19,5,19,5)

reaches 61 by KK. The exact adjacency rows of both resulting graphs are in `receipt.json`. The H17 comparison is a concrete warning against replacing a marked-state induction hypothesis by scalar attainability alone.

## 6. Source boundary and reproduction

Named sources at commit `ee07a5aac9cdfd8d3834fa834df4fb357fcaf769`:

- `session195/combined-growth/PROOF.md`: the four-coordinate interface.
- `session195/rooted-interfaces/EQUAL_RESPONSE_BARRIER.md`: the moving-root-only predecessor theorem and primitive H49 graph.
- `session195/zero-slack/PROOF.md`: the zero-slack/outdegree classification, reproduced above as needed.
- `session195/gadget-extension/bank67.txt`: the frozen finite control bank, copied unchanged here.

The preliminary fixed-modulus idea was discarded as new work: Session 193's `modular_sync/PROOF.md` already proves a stronger theorem and attributes basic modular completeness to Session 25. No modular-completeness novelty is claimed here.

Run `python3 verify.py` and `python3 -O verify.py` in this folder and compare their JSON outputs. The verifier uses explicit exception checks, not removable Python assertions. No dirty checkout was read or changed.
