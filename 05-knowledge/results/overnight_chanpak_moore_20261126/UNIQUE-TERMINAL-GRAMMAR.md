# The Fibonacci grammar of uniquely terminal tournament roots

2026-11-26. Continuation of [the mixed-spine result](TRANSFERRED-BOUNDARIES.md). This gives an all-order classification of the admissible right inputs to the address-9 unit lane, independent of the search for its second exit.

## 1. One terminal path determines a recursive normal form

Let $e(T,r)$ count directed spanning paths ending at the marked root $r$. Put $W=T-r$, let $C$ be its final strongly connected component, and let $P$ be the product of the spanning-path counts of all earlier components. The condensation of a tournament is a total order. Therefore

    e(T,r) = P ∑_{u in C with u→r} e(C,u).

Every vertex in $C$ contributes a positive terminal count when it is an in-neighbor of $r$: this is trivial for a singleton and follows by cutting a directed Hamilton cycle at that vertex for a nontrivial strong component. Every nontrivial earlier strong component has at least three spanning paths for the same reason.

Consequently $e(T,r)=1$ holds exactly for the following recursive grammar, besides the singleton:

1. Start with a strongly connected rooted tournament $(C,u)$ with exactly one spanning path ending at $u$.
2. Put any ordered chain $Q$ before $C$, with every arc from $Q$ to $C$.
3. Adjoin the new root $r$. Put $u\to r$ and $r\to C\setminus\{u\}$; orient contacts from $r$ to $Q$ freely.

Indeed $P=1$ makes every earlier component a singleton, and the sum must have exactly one term, equal to one. Conversely every constructed object has that unique terminal path. The grammar is unambiguous on rooted isomorphism classes: deleting $r$ recovers the ordered prefix and final component, and the unique in-neighbor in that component recovers $u$. The order of the chain vertices fixes the positions of the free contact bits.

Call the strong and unrestricted rooted class counts $A_n$ and $B_n$. For a nonempty $Q$, the result is strong precisely when $r$ beats its first vertex. If it does not, that vertex is a source; if it does, $r$ reaches every vertex and every vertex reaches $r$ via $u$. For an empty $Q$, the result is strong precisely when $C$ is nontrivial. Thus for $n\ge2$,

    B_n = Σ_{c<n} A_c 2^(n−1−c),

with $(A_2,B_2)=(0,1)$ and, for $n\ge3$,

    A_n = A_(n−1)+B_(n−1),
    B_n = A_(n−1)+2B_(n−1).

Writing $F_0=0,F_1=1$ for Fibonacci numbers gives **$A_n=F_{2n-4}$ and $B_n=F_{2n-3}$** for $n\ge2$. The total class generating function, including the singleton, is $(z-2z^2)/(1-3z+z^2)$. The two-state matrix is the square of the Fibonacci matrix in the coordinate order $(A,B)$; its dominant eigenvalue is $\varphi^2$. This golden-ratio growth has a structural source, not a numerical coincidence.

Every automorphism fixing $r$ must fix the unique terminal path pointwise. In particular the grammar's rooted classes account for exactly $(n-1)!B_n$ labeled orientations with the last label as root, matching the independent labeled census through order seven. The unrooted automorphism group acts freely on the set of terminal-unique vertices.

The same factorization makes a prime-sensitive extension possible: if $e(T,r)$ is a power of two all earlier components are still singletons, since tournament spanning-path counts are odd. For an odd positive target, earlier nontrivial components may only carry factors of that target. The branching occurs only in the terminal component and is the actual partition of the terminal count among its in-neighbors.

## 2. Exact second-exit census from the grammar

At left address 9 with its fixed mixed-spine edit, obtaining the opposite unit exit by toggling a right arc requires

    e=1,  e' ≡ 8 (mod 9),  9[(h'+R')−(h+R)]−e'+1 = 2,

where $R$ and $R'$ are the two rooted cover counts before and after the toggle. Enumerate the admissible starting right cores by the grammar, not all tournaments. For every possible toggle, the terminal count after editing is the sum of endpoint counts of the induced nonroot tournament over the new in-neighbor set. This makes the residue filter exact and cheap relative to full cover evaluation.

A root-incident flip can contribute only by making a different vertex of the final component an in-neighbor of the root (an earlier chain vertex does not affect $e$, while removing the unique terminal predecessor makes $e'=0$). That new predecessor must have terminal count 7 mod 9 inside $C$. A flip wholly within $Q$ multiplies the old terminal count by the spanning-path count of that newly edited ordered chain. Such a flip of an arc of gap $d\ge2$ gives $e'=1+2^{d-1}$, by the two-run subset choice at the reversed edge, so $d\equiv5\pmod 6$ is necessary in that case. The implementation still checks every arc.

**No opposite exit exists in this lane with any right core through order 13.** At order 13 the grammar generated 28,657 rooted isomorphism classes. Every one of their 78 possible single-arc edits was tested for the terminal residue, and the remaining 118,828 edits were fully profiled, without a hit. This is still a finite obstruction to one specific construction, not an all-order bound or a proof of Chan–Pak.

Replay: `g++ -O3 -std=c++17 unique_terminal.cpp -o unique_terminal && ./unique_terminal 13`. All counts in the search are exact integers.
