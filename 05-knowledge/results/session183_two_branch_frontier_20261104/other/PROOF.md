# Exact classification of clock-uniform digits in the transitive suffix chart

2026-11-04. **PROVED, chart-scoped.** A nonempty transitive suffix with physical ports admits no clock-uniform +2 transition to another transitive suffix. Its only clock-uniform +4 transition has binary port masks 7·2^a and 21·2^b, with a compensating prefix-clock shift a−b. This classifies all suffix lengths and all port sets, not merely a bounded edit menu. It does not classify arbitrary tournaments or certify native guarded operations.

## Setup and antecedents

Use the actual attachment M(TT_t;TT_m,W) from S182's [ordered-component proof](https://github.com/eliottcassidy2000/math/blob/3f2fa002801868c280fafed0c05a1d9060df6c86/05-knowledge/results/session182_physical_returns_20261104/induction/PROOF.md): root r beats the prefix, the prefix beats the suffix, W beats r and r beats the other suffix vertices. Number the suffix vertices 0,...,m−1 in their transitive order, and encode the port set by f=Σ_{i∈W}2^i. Strongness is equivalent to the terminal suffix vertex being a port. Thus m is precisely the bit length of f. Leading zero positions at the low end mean actual initial nonport vertices; they are not fictitious ports.

Write j(f) for the two-path joint response of this suffix and port set. The full count and physical root defect are

    H_f(t)=1−j(f)+2^t f,       E_r(2)−H_f(t)=j(f).

S182 constructs the +4 example 7→21 and its low-zero padding. The present result proves these are all such +4 maps and that +2 has no representative anywhere in this transitive chart. The immediately preceding [transitive-chart result](https://github.com/eliottcassidy2000/math/blob/3f2fa002801868c280fafed0c05a1d9060df6c86/05-knowledge/results/multiport_defect_grammar_20261103/TRANSITIVE_CHART.md) proves one digit per fixed source and scans odd masks below 2^16, explicitly leaving this all-mask classification open. Its pair formula and clock-shift reduction are antecedents reused here; the new content is the exact three-carry recurrence, its global bounds, and the resulting all-length classification. Its companion two-digit entry has a nontransitive TT2→C3 suffix on the +2 branch and is not contradicted. This is not an unrestricted spectrum obstruction.

## 1. A physical binary recursion for the joint response

Splitting off the lowest bit gives

    j(2n)=j(n),
    j(2n+1)=j(n)+n+(n mod 2),                  (1)
    j(0)=0.

For completeness, each pair of ports at positions i<k contributes 2^(k−i−1) to j, plus one extra contribution when k=i+1. The first contribution counts covers in which the first slot ends at k and the second starts at i: every intermediate vertex independently belongs to either slot. The extra contribution counts the opposite boundary order, possible only for adjacent positions. Removing a nonport at position0 leaves this pair count unchanged. Removing a port at position0 contributes Σ_{k≥1,b_k=1}2^(k−1)+b_1=n+(n mod2), proving (1). Hence low-zero padding does not change j.

## 2. An exact carry system and bounds

For c=0,1,2 define D_c(n)=3j(n)−j(3n+c), and set p=n mod2. All three D_c(0) vanish. Equation (1), applied to the physical port masks, gives

    D_0(2n)   = D_0(n),
    D_0(2n+1) = D_1(n)+4p−2,
    D_1(2n)   = D_0(n)−3n−p,
    D_1(2n+1) = D_2(n)+3n+3p,
    D_2(2n)   = D_1(n),
    D_2(2n+1) = D_2(n)+2p−2.                 (2)

These are exact identities, not numerical fits. Simultaneous induction proves

    −2n ≤ D_0(n) ≤ 2n,
    −3n ≤ D_1(n) ≤ 4n,
    −2n ≤ D_2(n) ≤ 2n.                     (3)

Here is the full induction check. At an even argument2n, the D0 bound follows directly; D1 lies between −5n−p and −n−p, which lies inside [−6n,8n] (n=0 is immediate, and p≤n otherwise); and D2 lies in [−3n,4n]⊆[−4n,4n]. At argument2n+1, D0 lies in [−3n+4p−2,4n+4p−2]⊆[−4n−2,4n+2]; D1 lies in [n+3p,5n+3p]⊆[−6n−3,8n+4]; and D2 lies in [−2n+2p−2,2n+2p−2]⊆[−4n−2,4n+2]. All recursive arguments are smaller except the immediate zero base identities.

## 3. Solving the two digit equations at every length

Since D0(2n)=D0(n), remove powers of2 and consider odd masks. For l≥0, equation (2) yields

    D_0(4l+1)=D_0(l)−3l−(l mod2)−2
             ≤−l−(l mod2)−2 < 0,

    D_0(4l+3)=D_2(l)+3l+3(l mod2)+2
             ≥l+3(l mod2)+2.                         (4)

For l≥5 the second bound is at least7. The five remaining values, evaluated exactly from (1), are

    l:                 0  1  2   3  4
    D_0(4l+3):         2  6  8  12  8.

It follows for every positive integer f that

    D_0(f)=4 has no solutions,
    D_0(f)=6 if and only if f=7·2^a, a≥0.             (5)

## 4. Classification with arbitrary fixed clock shifts

Suppose actual strong attachments satisfy

    H_g(t+k)=3H_f(t)+d

for all sufficiently large integer t, with fixed integer k and d∈{2,4}. Comparing the independent constant and 2^t coefficients gives

    2^k g=3f,          3j(f)−j(g)=d+2.              (6)

Write f=2^a u with u odd. The first identity forces g=2^b(3u), b≥0, and k=a−b. Padding invariance reduces the second identity to D0(u)=d+2. By (5), d=2 is impossible. For d=4 it forces u=7, and conversely f=7·2^a, g=21·2^b, k=a−b satisfy both identities because j(7)=6 and j(21)=12. Both suffixes have terminal ports and hence strong full attachments. Choose t≥max(1,1−k) so both prefixes are nonempty.

Thus the displayed family is exhaustive, including suffix sizes that were not searched and arbitrary fixed positive or negative clock shifts. In particular its +4 output, whose odd port mask is21, cannot take another +2 or +4 step while staying within this chart. No two-step digit word, let alone a reusable two-branch class, exists entirely here.

This boundary requires transitive suffixes and transitive prefixes. Nontransitive components, changes to the attachment structure, or nonuniform-in-clock transformations remain outside it. In particular the actual S182 lineage has a nontransitive fourteen-vertex suffix and is not ruled out. No edge-reversal or native-guard program is asserted by this classification.

## Verification

`verify.py` checks the pair-count interpretation against direct enumeration of actual two-path covers; checks the six recurrences and simultaneous inequalities; checks the exact digit classification on a large finite range as a control (the proof above handles the infinite tail); and independently counts actual full tournament Hamilton paths for every nonempty port mask through suffix order7 and prefix clocks1,2. The displayed +4 family is checked with actual matrices. Explicit exceptions survive Python optimization. See `verification.json`.
