# Every all-even-child quotient converse has modulo-sixteen ambient transfer

**Status: structural proof, independent review pending.** Quotient order and orientation are unrestricted. Every child must have even order at least four; corresponding target children must preserve complete ordered-cover profiles. This theorem supplies transfer, not local nonvanishing or arbitrary-profile reconstruction.

## 1. Main theorem and physical conventions

Let Q be any tournament on m>=1 fixed child labels. Let A_i and B_i be tournaments on the same physical child sets, each of even order at least four, with q_k(A_i)=q_k(B_i) for every k. Matching is childwise between A_i and B_i: profiles and even orders may differ freely between different children i and j. No equality across distinct child positions is assumed. Here q_h counts ordered spanning covers by h nonempty directed paths. Compare T=Q[A_i] and U=Q^op[B_i].

Use actual ports x->y, with incidences constant on every child. For any one fixed physical block query z, let Delta q_h(z) be the difference of the complete ordered-cover counts of the two augmented cores. Then, modulo sixteen,

    Delta q_2(z)=2 Delta q_1(z),
    Delta q_h(z)=0 for every h>=3.                       (1)

Consequently EVERY nonzero local Hamilton contrast modulo sixteen survives every fixed ambient root substitution, when both ports inherit the original core root's exterior incidences. The same core query works in every ambient. The ambient contrast is

    [c_1+2c_2] Delta q_1(z) modulo sixteen,               (2)

where c_1+2c_2 is odd. It need not equal the original residue, but cannot annihilate it.

All child names, actual physical vertices and port masks stay aligned. The theorem does not assert that such a nonzero local contrast exists for every Q or every child assignment. It is not the false generic modulo-sixteen transfer principle for arbitrary pairs of modules.

## 2. A universal path-slot parity lemma

Let R be any tournament on a fixed label set containing a distinguished label i and at least one other label. For k>=1 define K_h(R;i^k) to count ordered h-tuples of nonempty directed words that collectively contain i exactly k times and every other label exactly once. An uncut word cannot have equal adjacent labels; cuts separate the h words. Thus this is the exact repeated-child run kernel, before multiplying by child profile coordinates.

For k in{2,3,5}, the complete higher-cover parity table is independent of every arc of R:

    k=2: K_h=1 modulo2 exactly at h=3, for h>=2;
    k=3: K_h=1 modulo2 exactly at h=3, for h>=2;
    k=5: K_h=1 modulo2 exactly at h=5, for h>=2.          (3)

All unlisted higher h have parity zero. At h=1, K_1(R;i^3) and K_1(R;i^5) are even; K_1(R;i^2) is not asserted even.

**Proof.** Permute the h ordered path slots. Two identical paths can contain no label other than i, because every other label occurs once globally. A path using only i must be the singleton word i, since equal adjacent occurrences would be one run rather than two. Therefore, if d slots are singleton-i paths, the stabilizer is exactly S_d and the orbit size is h!/d!.

There is at least one nonsingleton-i slot because another label occurs. If h-d>=2, the product(d+1)...h contains an even factor, so that orbit contributes zero modulo two. The sole surviving possibility for h>=2 is h=d+1 odd: there are h-1 identical singleton-i paths and one remaining path containing every other label, together with r=k-h+1 copies of i.

If r<0 there is no such path. If r=0 or1, its count is an ordinary tournament Hamilton count and is odd. If r>=3, its r-1 intermediate return segments have disjoint nonempty supports in the other labels. Permuting those segments is free and preserves all joins, so the count is divisible by(r-1)! and hence even. The only unresolved multiplicity is r=2. For k in{2,3,5}, the equation r=2 forces h=k-1, which is either h=1 or even and therefore cannot be a surviving higher-cover case. This gives exactly(3). The same return-segment argument directly proves the h=1 claims. QED.

For the two-port kernels used below, x and y are unique other labels, so the lemma's nonempty-complement condition always holds. No path-slot action is applied after forgetting cuts, and no fixed-point-free assertion is made when singleton-i slots coincide.

## 3. The complete even-child domain and its affine contrast

The inherited normalized profile domain for every actual even-order child of order at least four is

    q(A_i)/H(A_i)=(0,1,14+8e_i,4+8f_i,8,8g_i,0,...),
    e_i,f_i,g_i in{0,1}, modulo sixteen.                (4)

Every H(A_i) is odd. The upper domain contains every actual profile; it need not consist entirely of realizable vectors. Its source is the independently reviewed S188 `even_cyclic_mod16/PROOF.md`, at immutable commit `8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7` under `05-knowledge/results/session188_structural_bridges_20261104/`.

Let g_k=(-2)^(k-1). Modulo sixteen its coordinates are(1,14,4,8,0,0,...). This is the formal geometric profile from the independently reviewed Session189 positive-Seidel-energy theorem. That theorem proves, for every quotient Q and every child-uniform two-port query, the complete geometric contrast identity

    Delta_Q^g(t)=(1+2t) delta_Q^g,                      (5)

where t^(h-1) records q_h. It holds for all quotient orders; no finite prime-seven or prime-eight classification is used here.

Write Lambda=product_i H(A_i), a common odd unit. Expand each profile(4) around g in the exact cut-aware run functional. Products of two perturbations have factor64 and vanish modulo sixteen. In a first-order term, every other child must occur once, because any repeated other child supplies an additional even profile factor. Therefore, for every h,

    Delta q_h / Lambda
      = Delta_Q^g[h]
        +8 sum_i [e_i Delta K_h(i^2)
                    +f_i Delta K_h(i^3)+g_i Delta K_h(i^5)]
          modulo sixteen.                              (6)

Here each kernel is computed on the actual bare quotient-plus-port tournament at z, in the two orientations. Formula(6) concerns contrasts: individual higher-cover derivatives need not vanish. For example, the parity-one entries in(3) occur on each side and cancel only after subtraction.

The parity lemma makes every derivative contrast in(6) zero for h>=2. At h=1, the k=3,5 terms also vanish. Thus for an integer eta determined by the q2-bit corrections,

    Delta q_1=Lambda(delta_Q^g+8 eta),
    Delta q_2=2 Lambda delta_Q^g,
    Delta q_h=0 for h>=3,                               (7)

modulo sixteen. Multiplying the first line by two eliminates its correction, proving(1). The common odd Lambda has not been discarded.

The formal geometric series causes no uncontrolled tail: modulo sixteen, only bounded total run excess contributes. The affine expansion handles all actual normalized coordinates. Independently,16 divides h! for h>=6, so the actual higher-cover tail vanishes as well.

## 4. Ambient consequence and boundaries

Give both ports the original root incidence at each fixed exterior vertex. The port-augmented core remains a common module. Exact run substitution into the ambient writes its Hamilton count as sum_h c_h q_h, where c_h are nonnegative integers determined by the exterior and c_1=H(R) is odd. Equation(1) gives(2).

In fact c_1+2c_2=H(R[TT2]), since a transitive pair has profile(1,2,0,...), so this multiplier is itself a tournament Hamilton count and is odd. Odd multiplication cannot kill a nonzero residue modulo sixteen.

The result is restricted to quotient converse with common even-child profiles and child-uniform queries. It does not apply to arbitrary same-order module pairs, singleton or pair children, odd children, nonuniform marked ports inside a child, unequal source/target child profiles, or independently relabeled responses. The known generic modulo-sixteen transfer failure is therefore untouched.

If a local nonvanishing theorem is available in this profile domain, (1) supplies its arbitrary fixed ambient extension without a new cut-aware quotient enumeration. It supplies no local nonvanishing by itself. Root-containing original-module families may then be treated by the separate negative-dual monotonicity argument, with that argument's original-module and fixed-exterior hypotheses retained.

## 5. Reproducible corroboration and data scope

Run `python3 verify.py` or `python3 -O verify.py` with Python 3 and a C++17 compiler. No third-party package is used. The executable retains 16,212 complete exact repeated-child cut-kernel vectors, each q1 through q9, on every labeled tournament of order two through five and every distinguished label, for k=2,3,5. Every possible nonzero coordinate is included; it checks the universal table including h beyond five.

The 944 literal physical controls retain the actual source and target augmented rows, physical port masks, child sets, port and exterior names, and exact q1 through q5. Covers q_h for h>=6 are zero modulo sixteen by the free physical path-slot action, so these five exact integers determine the complete modulo-sixteen profiles. Controls include nontransitive even children, unequal even child orders four and six, and profile-preserving interior converses. There are 400 local queries, eight actual cyclic-root ambient queries with multiplier five, and all 64 coarse queries for an outside-domain hostile example. Maximum augmented order is eighteen; the integer bound 18! binomial(17,4)<2^64 protects every stored exact coordinate.

The outside-domain example is the directed triangle inflated by transitive children of orders 2,2,3. At block query (a,c)=(0,1), with the physically aligned interior converses specified in the data, its exact first five contrasts are (116,2520,15960,47040,73920). Modulo sixteen these are (4,8,8,0,0); the nonzero third coordinate violates (1), illustrating why the domain restriction cannot simply be dropped.

`kernels.txt` and `physical_records.json` retain complete finite data, including the hostile bank rather than only its first witness. The theorem is the slot-permutation and affine proof above; finite checks corroborate it rather than replace it.

No new claim of historical priority is made for the free permutation of return segments. Its use is inherited from the S188 return-segment work. The new argument explicitly analyzes the stabilizers of global ordered path slots, which may contain identical singleton-i paths.
