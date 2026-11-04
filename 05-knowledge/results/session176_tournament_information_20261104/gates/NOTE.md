# Shared gates, tournament dependence, and the cost of a boundary state

Session 176. **Proved elementary consequences of the inherited one-feedback constructor; no general Boolean circuit lower bound or novelty priority is claimed.** The packet makes a precise distinction between three quantities: pairwise dependence, arithmetic interface rank, and exact streaming memory.

## 1. Source and model

The tournament T(b), for b in {0,1}^n, has base vertices 1,...,n with i→j iff i<j, and apex z with z→i iff b_i=1. Session 175 `construction/README.md`, §1, already proves its path-count formula and three-register recurrence. We use those results and recheck them independently by literal subset DP. All subsequent identities are over the integers or rationals unless explicitly stated otherwise.

Let H_n(b) count directed Hamilton vertex words of T(b). The inherited formula expands to

    H_n(b) = 1 + sum_{i=1}^{n-1} 2^(n-i) b_i
               - sum_{1≤i<j≤n} 2^(j-i-1) b_i b_j
               - sum_{i=1}^{n-1} b_i b_(i+1).                 (1)

To derive this without relying on a numeric fit, delete z from a Hamilton word. Its left and right parts are increasing. With both parts nonempty, let i be the left maximum and j the right minimum. When i<j there is exactly one choice, j=i+1, contributing (1-b_i)b_(i+1). When j<i, the i-j-1 intervening vertices may be freely assigned left or right, contributing 2^(i-j-1)b_j(1-b_i). The two endpoint contributions are b_1 and 1-b_n. Expansion gives (1). For n=0 the convention is H_0=1.

## 2. Complete pairwise interaction, linear arithmetic circuit

For i<j, the discrete mixed derivative (difference after setting b_i and b_j to 0 or 1, with all other bits fixed) is

    Δ_i Δ_j H_n = -2^(j-i-1) - 1_{j=i+1}.                 (2)

It never vanishes. The pairwise interaction graph is therefore the complete graph K_n, with n(n-1)/2 edges, for every n≥2. No higher derivative survives, since (1) is quadratic.

Nevertheless, start (P,L,h)=(0,0,1) and update on each input bit c by

    h' = h + (1-c)(P+L),
    P' = 2P+c,
    L' = c.                                               (3)

Then h=H_n(b). This is the inherited Session 175 recurrence, now written as a branch-free straight-line circuit. With binary arithmetic operations +,-,×, and only constants 0,1,2, each step uses at most six gates: P+L, 1-c, their product, addition to h, 2P, and addition of c. Thus at most 6n gates suffice. The outputs P',L',h' use the old register values simultaneously.

**Barrier.** A proposed lower-bound principle that charges a distinct arithmetic gate to every nonzero pairwise derivative is false, already for an explicit tournament counting invariant. The shared prefix gate P serves many later interactions. This does not refute a properly proved gate-elimination theorem: it identifies exactly the multiplicity a charging argument must control. Also, an arithmetic gate with growing integer operands is not a constant-size Boolean gate. A straightforward binary implementation has polynomial size, for example O(n²) with linear-size ripple arithmetic on O(n)-bit registers; no B2 lower bound follows from (2).

## 3. Every ordered cut has only two bilinear interaction channels

Cut the input between k and k+1. The cross part of (1) is exactly

    -A_k(b_1,...,b_k) B_k(b_(k+1),...,b_n) - b_k b_(k+1),

where

    A_k = sum_{i=1}^k 2^(k-i)b_i,
    B_k = sum_{j=k+1}^n 2^(j-k-1)b_j.                     (4)

The k×(n-k) cross-coefficient matrix therefore has rank at most two over Q. If k≥2 and n-k≥2, it has rank exactly two: the minor on rows k-1,k and columns k+1,k+2 is

    [ -2  -4 ]
    [ -2  -2 ],

with determinant -4. If one side has size one, the nonzero matrix has rank one. This proves the exact number of channels needed **in a representation as a sum of products of linear forms in the two sides**. It is not minimality over arbitrary circuits, arbitrary nonlinear encodings, or arbitrary bit protocols.

Thus a dense derivative graph does not imply a large ordered-cut rank. The long-range interactions are geometric and collapse into one weighted channel; the second channel records the adjacent-bit correction. This makes the shared gate in (3) algebraically visible.

## 4. Constant arithmetic width still needs linear exact memory

For a prefix word w of length k, let P(w) be its binary integer value. Applying (3) twice gives

    H_(k+2)(w00) = h(w)+3P(w)+L(w),
    H_(k+2)(w01) = h(w)+ P(w)+L(w).

Consequently

    H_(k+2)(w00) - H_(k+2)(w01) = 2P(w).                (5)

**Distinguishability theorem.** All 2^k prefixes of length k have different exact residual output functions on suffixes of length two. Indeed, if two prefixes had the same outputs for every such suffix, subtracting the outputs for 00 and 01 in (5) would give the same binary value P, hence the same prefix.

Any deterministic one-pass exact algorithm computing H_(k+2), with no access to the prefix after the cut, must therefore have at least 2^k distinct states at that cut and hence at least k bits of working information. This assumes the final output is produced from the state and suffix; an earlier output tape readable by the suffix processor is part of the state. Equivalently, deterministic one-way communication from the prefix holder to the two-bit suffix holder, who must output H exactly, requires k bits in the worst case. Sending the prefix achieves k, so this one-way complexity is exactly k.

This is a state-count argument, not an entropy heuristic and not a Boolean circuit-size lower bound. It coexists with the rank-two bilinear interface and the linear-size arithmetic circuit. The register P in (3) is literally a growing k-bit encoding of the entire prefix: constant register count is not constant information.

The same theorem excludes any exact finite-state automaton with a fixed finite state set for the full family of output functions H_n under this streaming model. It does not exclude a fixed-dimensional linear or arithmetic representation with unbounded entries; (3) provides one.

## 5. What this transfers to the larger program

The useful common object is an **extension response**, not the isolated value of a tournament statistic. Two partial histories may be merged only if their responses to every allowed continuation agree (or have the precisely declared compensated equivalence). Here two probe continuations already distinguish every prefix. A small number of algebraic coordinates can encode that large distinction space.

For circuit elimination, the lesson is equally precise: first identify which gates are shared by multiple proposed charges, then prove that the selected deletions are distinct. Dense dependency alone cannot supply that proof. The older Session 117/118 gate-elimination work should remain separate: its B2 gate basis and legal-restriction coverage hypotheses are not established or transferred by this arithmetic example.

Useful scope table:

- All pairwise interactions present: proved by (2).
- Ordered bilinear cut rank: exactly two away from end cuts, by (4).
- Exact streaming state bits: at least k after k inputs with two inputs remaining, by (5).
- Arithmetic circuit size: at most 6n in the specified basis, by (3).
- General B2/NP circuit lower bound: not supplied.

## 6. Reproduction and limits

`python3 check.py > certificate.json` and `python3 -O check.py > certificate.optimized.json` require only the Python standard library. The checker independently compares literal tournament subset DP, the quadratic polynomial and the streaming circuit for every word through length eight; checks mixed derivatives on all base assignments through length seven; verifies cut factorizations and rank witnesses through n=100; and checks prefix response injectivity through k=12. Its checks raise explicit exceptions and remain active under -O. The all-length claims rely on the proofs above, not extrapolation from these finite tests.
