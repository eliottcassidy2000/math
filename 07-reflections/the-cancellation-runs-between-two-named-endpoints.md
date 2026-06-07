# The cancellation runs between two named endpoints

*monad-explorer-2026-06-07 (deep-research, 12th session). On THM-438's cycle-rank
triangle `t(k,m)` and the meaning of "the Catalan is a cancellation, not a count."*

---

## The slogan, made literal

Across nine addenda, THM-438 kept circling one thesis: the appearance of the Catalan
numbers `C_k` in the Paley cluster integrals is a **cancellation**, not a count. The
even-series patterns over-count wildly (`A088368 ~ e·k!`); the Möbius signs collapse the
over-count to `C_k`. ADDENDUM-3 sharpened it ("the count is so unstructured it is
uncatalogued"); MISTAKE-062 made it sharper still (the unsigned support count `1,3,13,67,
383,2351` is in no OEIS sequence).

This session closes the loop with a clean, almost visual statement. Attach to the cycle-rank
triangle every natural one-variable sequence and ask OEIS which are catalogued. The answer is
stark — **exactly two**, and they are the two *ends* of the cancellation:

```
   WILD END                                            TAME END
   diagonal  t(k,k) = A088368(k)  ~ e·k!     ────►     signed sum  Σ_m(−1)^m t(k,m) = (−1)^k C_k
   "partitions into noncrossing lists"                 Catalan A000108
   the all-pairings / moment over-count                the free-cumulant answer
```

Everything strictly between them is structureless: the top residues `P_m(1)=1,3,20,181`,
the sub-diagonal `9,72,580,4845`, column 3 `13,72,230,560`, the *unsigned* row sum
`1,4,23,160,1262,10944`, the #-lines counts — all return *No results*. The cancellation is a
bridge whose two abutments are famous landmarks and whose every intermediate pier is
anonymous. That is what "a cancellation, not a count" means, stated as a fact about OEIS.

## The wild end has a name now

The diagonal — long described in this thread only by its mechanism ("the bigon-tree /
all-pairings over-count") — is **A088368: partitions of `[n]` into sets of noncrossing
lists** (Callan, 2007). It is not uncatalogued at all. It carries two functional equations,
`A(x)=Σ n!x^n A(x)^n` and `A(x/F(x))=F(x)` with `F=Σ n!x^n`, and a rigorous asymptotic
`a(n) ~ e·n!` (Kotesovec). The thread's own model recovers Callan's exactly: max-rank
patterns are **doubled plane trees**, the Euler tour visits each vertex `deg(v)` times so the
weight is `∏(deg(v)−1)!` over noncrossing tours = sets of noncrossing lists.

This matters beyond bookkeeping. It means the *over-count itself* is a well-understood object
with a clean asymptotic; the mystery is entirely in the cancellation, which sends a sequence
growing like `e·k!` to one growing like `4^k/k^{3/2}`. The moment method's two faces —
Gaussian (all pairings, `e·k!`) and semicircular (non-crossing, `C_k`) — are *both* named
here; the Möbius inversion is the only unnamed thing, and it is unnamed because its
intermediate data is genuinely without arithmetic structure.

## A correction, and its mirror

ADDENDUM-9 had flagged `A088368 ~ e·m!` as FALSE, replacing it with an empirical `~m!(m+2)/2`,
because `a(m)/m!` for `m≤7` rises *past* `e`. It does — and then it turns around. The ratio
**peaks at m=8 (≈4.36) and descends** back to `e` (3.14 at m=20, still falling). ADD-9 saw
only the rising side of a slow, overshooting limit. The original slogan stands; the
"correction" was the error (MISTAKE-063).

This is the exact mirror of MISTAKE-062 from three sessions earlier. There, five common terms
were over-trusted as a sequence *identity* (the A215257 coincidence). Here, six small terms
were over-trusted as an asymptotic *refutation*. Both failures are the same epistemic shape:
**a factorial-scale quantity still visibly moving at `n≤8` has converged to nothing.** Slow
asymptotics overshoot and turn; a monotone prefix is not a limit; and OEIS already knew.

## `e` at the wild end

A small bonus. Writing the column as `t(k,m)=(k)_m h_m(k)` (ADD-9), the bridge polynomial
`h_m` has two asymptotic ends:

```
   h_m(m)  = A088368(m)/m!  → e        (wild end, the noncrossing-lists overcount)
   h_m(−1) = (2^m−1)/((−1)^m m!) → 0   (tame end, Mersenne, super-exponentially)
```

So `e` is the wild-end limit of the bridge. The project's "everything is the triangle" motif
lists four constants (`√2, π, e, γ`); `e` already appeared there through Stirling/Gamma growth
in Burnside. It now appears a second, independent time — as the limiting value of the bridge
polynomial at its wild root-cluster, via Callan's `e·k!`. When a constant shows up twice by
unrelated routes inside one object, that is the object telling you the two routes are not
unrelated: both are the factorial (the symmetric group's order) seen once through orbit
counting and once through the all-pairings over-count. The cancellation that takes `e·k!` to
`C_k` is, in free-probability terms, exactly the passage from the classical (Gaussian,
`e`-flavored) moments to the free (semicircular, Catalan) cumulants — the same passage the
loop equation `xF²+F=1` encodes (ADD-4). The two named endpoints are the two probability
theories; the anonymous middle is the functor between them.

## What this buys the program

The diagonal (wild end) is now closed-form and named. The two open handoffs both live at the
**tame end**: handoff #1 `(k)_m | t(k,m)` and handoff #2 `g_m(−1)=(−1)^m(2^m−1)`. The cleanest
new sub-target is finite and number-theory-free: prove the bijection *doubled plane trees ↔
Callan noncrossing lists*, upgrading `t(k,k)=A088368(k)` from VERIFIED+OEIS-match to PROVED,
with `h_m(m)→e` falling out as a corollary of Kotesovec. That is a self-contained foothold on
one abutment of the bridge — and once one end is proved, the cancellation has somewhere solid
to start from.
