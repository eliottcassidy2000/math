# Crossing-matching continuants: a reduced O^m J^2 gate, with lossless pair decoding

Status: PROVED elementary identities and consequences below; FINITE-EXACT literal checks. Independent audit [PASS](../audit/GATES_AUDIT.md). This is an exact reduced sublanguage of the inherited O/Z/J construction, **not a new primitive and not Chan–Pak completeness**. No literature-priority claim.

## 1. Inheritance, scope, and what survived

The selected input is main `38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817` (materialized read-only in `session221-stage`; newer root-change packet inspected from this exact commit using `git show` in `session219-stage`).

Closest proved mechanism: ordered root amalgamation and the sole-outneighbor insertion identity. Session216 `paired-reset/FINITE_SEARCHES.md` gives the seven-coordinate amalgamation; Session220 `session220_binary_envelopes_and_selectors_20261107/exploration/REPORT.md` recovers a one-seam four-coordinate transfer. THM-4904, `01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md`, contains the newer O/Z/J chamber and bounded-count work. The present gate is exactly its `O^m J^2`, after swapping the two last birth labels.

Hostile / corrected near miss: commuting Q/Z counters do not encode independent digits; Session220's constant ribbon family has a doubled mode coupled to a much larger physical eigenmode. The least-used sidecar pursued here is the root-starting two-path-cover count after moving the root to the penultimate birth. The output is a universal **two-coordinate** macro with independently chosen dyadic parameter, exact inverse on pairs, and a sharp warning against scalar continuation.

Other candidates were tested and retired as recovered machinery: two disjoint reverse arcs (crossing and nested formulas verified on 64 and 27 parameter triples, respectively), and a general bounded-active-port transfer. THM-4239 `all-n-strong-two-reversal-presentation-sphere-and-tail-fiber-product-boundary`, THM-4748 `two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity`, and THM-4878 `graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses` already cover that territory. The six-state exploratory specialization is not promoted as a new automaton theorem and is omitted from this portable release.

THM-4744 `01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md` concerns reverse arcs sharing an endpoint and gives a product law. Here all reversed-edge endpoints are distinct, and consecutive intervals cross. Their counts therefore need not be Moon products. This distinguishes the presented literal family, but not the inherited macro calculus.

Concept board: (i) multiple-seam physical adjacency; (ii) moving-root continuation; (iii) finite-state reduction; (iv) programmable dyadic arithmetic; (v) inverse sidecars; (vi) support density. The positive signal was the reset relation in Section 2. It replaced the initial generic automaton target by a small exact arithmetic carrier. The density and collision tests delimit, rather than erase, that survivor.

## 2. Universal physical macro

For a nonempty rooted tournament `(T,a)`, let

- `h=H(T)`;
- `e=end_a(T)`;
- `R=R(T,a)` be the number of ordered spanning two-directed-path covers `(P,Q)` with `P` starting at `a`. The path `Q` may be empty, and `P` is nonempty. Set `H(empty)=1`.

For an integer `m>=0`, form `G_m(T,a)` as follows.

1. Append `m` global sinks in transitive order.
2. Append one further global sink `b`.
3. Append `c`, with `c->a` and every other earlier vertex pointing toward `c`.
4. Mark `b` as the new root.

This is complete physical adjacency, adding exactly `m+2` vertices. The only backward arc in the appended-order description is `c->a`. Equivalently, append the `m` sinks, perform two inherited sole-exit/new-root `J` births, then swap their birth labels. Thus this is `O^m J^2`, not an independently new primitive.

**Universal gate theorem.** With `q=2^(m+1)`,

    h' = h + q R,
    R' = 2h + q R - 1_[m=0] e.

In particular, for **every** rooted input and every `m>=1`,

    (h,R) -> (h+qR, 2h+qR),   q in {4,8,16,...}.             (1)

No extra endpoint or deletion coordinate is needed for iterating these positive-m macros. This is an exact all-input formula, not a fit to a seed orbit.

### Proof, including the boundary correction

If a new vertex `c` has exactly one outneighbor `a` in an old tournament `U`, paths ending at `c` are old Hamilton paths not ending at `a`, contributing `H(U)-end_a(U)`. Every other Hamilton path has the unique form `Q,c,P`, with `P` starting at `a`; these are counted by `R(U,a)`. Therefore its count is `H(U)+R(U,a)-end_a(U)`.

Appending `r` global sinks preserves `h` and multiplies `R` by `2^r`: independently assign each new sink to one of the two labeled paths, placing each assigned subsequence at that path's end. If `r>=1`, no Hamilton path ends at the old root.

Apply these facts to `U=T+TT_(m+1)` before adding `c`. This gives `h'=h+2^(m+1)R` for every `m>=0`.

For `R'`, the path starting at `b` has exactly three possibilities:

- It is just `b`. The other path is Hamiltonian in `T+TT_m+c`, contributing `h+2^m R-1_[m=0]e`.
- It is exactly `b,c`. The other path lies in `T+TT_m`, contributing `h`.
- It starts `b,c,a` and continues. Removing `b,c` leaves a root-starting two-cover of `T+TT_m`, contributing `2^m R`.

The alternatives are disjoint and exhaustive because `b` has sole outneighbor `c` and `c` has sole outneighbor `a`. Their sum proves the theorem.

The `m>=1` condition is necessary for this two-coordinate reset. A singleton input has `(h,R,e)=(1,1,1)`. At `m=0` the output is a directed triangle with `(h',R')=(3,3)`, not `(3,4)`.

## 3. Literal crossing-matching compiler and variable continuants

Start from the singleton with `(h,R)=(1,1)`. Choose any finite program `(m_1,...,m_k)` of positive integers and apply the gates in that order. Write `q_i=2^(m_i+1)`.

The resulting tournament has

    n = 1 + sum_i(m_i+2) = 1 + sum_i(log_2(q_i)+1)           (2)

vertices. Its explicit order is

    a_0, F_1, a_1, c_1, F_2, a_2, c_2, ..., F_k, a_k, c_k,

where `F_i` is a transitive string of `m_i` vertices. All arcs point forward except

    c_i -> a_(i-1),  i=1,...,k.                              (3)

The `k` backward arcs form a matching. Consecutive intervals cross; nonconsecutive intervals are disjoint. At most two backward arcs cross any cut, even when their spans are unbounded. The tournament is strongly connected: initially it is a point; given a strong old tournament, every old vertex reaches `a`, then the new vertices, while the new vertices reach `c->a`. This also proves strongness of the first nontrivial stage.

Let `h_j` be the count after `j` stages. Equation (1) gives `R_j=h_j+h_(j-1)`. With the formal initial values `h_(-1)=0, h_0=1`, the exact programmable recurrence is

    h_j = (q_j+1) h_(j-1) + q_j h_(j-2).                    (4)

The coefficient `q_j` is independently chosen at every stage. This is a positive variable generalized continuant, not merely a constant-block spectral ray. Equivalently the pair `(h,R)` is multiplied by

    M_q = [[1,q],[2,q]].                                    (5)

The matrices are noncommuting. For programs `(1,2)` and `(2,1)`, the literal counts are respectively `53` and `49`. A constant `m_i=1` gives

    1,5,29,165,941,5365,30589,174405,994381,...,
    h_j=5h_(j-1)+4h_(j-2).

Every output has an explicit short witness. Since `h_j>q_j h_(j-1)` and `log_2 q_j>=2`, (2) yields

    n <= 1 + (3/2) log_2 h_k.                               (6)

This bound is only for outputs of this compiler, not all odd targets.

## 4. Lossless pair decoder; the scalar is insufficient

Every reachable pair `(h,R)` records the previous count exactly:

    p = R-h = h_(k-1).

If the pair is the singleton `(1,1)`, stop. Otherwise `p>=1`. Since `0<=h_(k-2)<p`, (4) gives

    q <= (h-p)/p < 2q.

There is exactly one power of two `q>=4` in this dyadic bracket. Compute it by integer arithmetic, and recover

    previous_previous = [h-(q+1)p]/q.

Require this number to be an integer in `[0,p)`, and repeat on `(previous_previous,p)`. Accept exactly when the terminal pair is `(0,1)`. This proves that each reachable pair has **one program**, and provides an exact terminating recognition/decoding algorithm. The successive counts decrease by a factor greater than four. No scalar inversion theorem follows.

An independent two-adic inverse is also available. Beyond the first gate, `R_old≡2 (mod 4)` because both adjacent counts are `1 (mod 4)`. Therefore from a noninitial pair, with `p=R-h>1` and `d=2h-R`,

    q = 2^(v_2(d)-1),  h_old=p,  R_old=d/q.

The first gate is handled separately by `p=1`, `q=h-1`, `R_old=1`. The verifier checks both inverse algorithms on the full finite census below.

**Hostile scalar collision.** The two programs `(2,2)` and `(3,1)` have the same count `89`, on nine vertices, but their states are `(h,R)=(89,98)` and `(89,106)`. Appending the same `m=1` gate gives respectively `481` and `513`. Thus `h` alone is not a lawful continuation coordinate even within this constrained family. The exhaustive monotone search through count `10^6` finds no smaller collision; in particular the claim of minimal collision at `89` is FINITE-EXACT and complete below `89`, not an isomorphism claim.

## 5. Sharp lane boundaries: congruence and zero density

Every output is `1 (mod 4)`. More precisely,

    h_k≡5 (mod 8) if m_1=1,
    h_k≡1 (mod 8) if m_1>=2.

Indeed the first count is `q_1+1`. Thereafter `R≡2 (mod 4)`, so every increment `qR` is divisible by eight. This permanently excludes every odd target that is `3 (mod 4)`.

The entire infinite-parameter family also has zero density. Put `L=floor(log_2 X)`. An output at most `X` requires `sum_i log_2 q_i<=L`. The number of compositions of an integer `s` into parts at least two is `F_(s-1)`; including the empty program, the number of programs with this sum at most `L` is `F_(L+1)`. Consequently

    #{distinct outputs <=X} <= F_(floor(log_2 X)+1)
                           = O(X^(log_2(phi))),
    log_2(phi)=0.69424... <1.                                (7)

Collisions can only lower the support size. This is a theorem about this exact reduced language, not all O/Z/J words, all cutwidth-two tournaments, or all matching-reversal tournaments. It explains why a genuinely programmable, invertible pair carrier still fails to supply target coverage.

## 6. Reproduction and finite universe

Run from this directory:

    python3 verify_gate.py > VERIFY.json
    python3 -O verify_gate.py > VERIFY_OPT.json
    cmp VERIFY.json VERIFY_OPT.json

Only Python's standard library is required. All seventeen verification and input-precondition tests use explicit `require`/`raise` checks, so Python optimization cannot remove them. Normal and optimized runs produce byte-identical receipts, preserving the exact counts below. The verifier uses independently implemented subset dynamic programs on literal adjacency, including all-subset root-start counts to compute `R` by convolution. It does not use the six-profile exploratory transfer.

Passed:

- Every labeled tournament of orders `1..4`, every root: `285` rooted presentations.
- For each, all `m=0,1,2,3`: `1,140` literal gate checks, including `285` corrected zero-padding controls.
- In every gate case, actual adjacency is compared to `O^m J^2` after swapping the last two labels.
- Every word of length at most four over `m in {1,2,3,4}` whose physical order is at most sixteen: `110` literal program checks. Counts, both decoders, strongness, matching endpoints, and cutwidth at most two are checked.
- Exact monotone recursive census of **all** positive-m programs with count at most `1,000,000`, with no separate length or digit cap: `4,047` words including empty, `2,737` distinct scalar outputs. Both inverse algorithms agree on every word. The loop stops at each node only when the next increasing dyadic coefficient forces the new count above the bound; every later stage increases the count.
- The count-89 collision and its two different literal continuations are independently checked by subset DP.

No repository files were edited, committed, or published by this worker. The deliverable is `REPORT.md`, `verify_gate.py`, `VERIFY.json`, `VERIFY_OPT.json`, and `SHA256SUMS`. Retired exploratory code is not needed for the stated results or included in this portable release.
