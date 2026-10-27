# ROOT transport: a positive affine carry and a minimal clocked norm

## Native coordinate, with its constant retained

For a positive odd primitive seed `s≡1,3,7 (mod 8)`, put `B(x)=4x+1`.
Take `a=2` for `s≡1 (mod8)` and `a=1` otherwise, and write

    (3s+1)/2^a = B^j(P(s)),       t=a+2j.

With `X=3s+1`, `Y=3P(s)+1`, the exact identity is

    2^t Y = 3X + 2^a.

Indeed `B^j(u)=4^j u+(4^j−1)/3`; substitute and multiply by three.
This positive carry is part of the dynamics, not a removable error term.
For example the seeds `33,41,49,57` have successors `25,31,9,43` and
clocks `2,2,4,2`. The root construction is multiplicative, but its branch
conditions and this carry are additive data.

The original finite source charges are Eisenstein integers

    A_N(s)=Σ_(j:B^j(s)<4N) ω^(B^j(s))
             ·[(1+2ω) if B^j(s)<2N, else ω],
    μ_N(s)=|A_N(s)|².

Zero charges must not silently erase the identities of original sources.
The following construction works for any nonnegative source weights; it
therefore does not depend on a claimed all-scale estimate for these charges.

## Why the original weighted quotient can fail

A pullback on a weighted function space is defined modulo zero norm only
if a positive-weight source never pulls from a zero-weight target.
There is a concrete exterior obstruction for the earlier forward-closed
627-state frontier. For `N≥3360`, choose seeds in `[3N,4N)` congruent to
`7,11,15 (mod16)`. There are at least `3 floor(N/16)>627` candidates.
Each has `a=1,j=0`, successor `(3s+1)/2>4N`, source weight one and target
weight zero. The successors are distinct. At least one lies outside the
frontier; forward closure then forces its source outside as well. Thus the
exterior quotient does not repair the problem. At `N=4096`, the earlier
explicit witness is `12295→18443`, clock one. The frontier and that
membership check are cited finite inputs from the
[preceding ROOT checkpoint](https://github.com/eliottcassidy2000/math/blob/318d10ff81b7e9f147e22162d32762059ef15931/05-knowledge/results/session75_overnight_black_root_spectral_20261027/README.checkpoint3.md),
not data included in this small package.

## The clock lift and its least occupation weight

Fix an integer budget `H`. Lift an orbit to `(x,τ)`, starting at clock zero;
perform a full native edge only when its arrival clock is at most `H`.
Otherwise retain the current state as terminal. Boundary hits are terminal
as well. No fractional edge or uncharged holding step is introduced.
Because every clock increment is positive, this is a finite directed
acyclic graph with at most `H+1` states per original source. Preserve source
labels and terminal ports whenever subsequent observables can name them.
A merger is valid only for observables constant on its fibres.

Let `i` inject `μ(s)` at `(s,0)`, and let `K` push weights forward along
nonterminal edges, adding incoming weights at shared states. Then
`K^(H+1)=0` and

    W = Σ_(j=0)^H K^j i,          W=i+KW.

This is the unique solution of that equation and the least nonnegative
weight satisfying `w≥i+Kw`: iterate the inequality and use nilpotence.
For the terminal-zero pullback `Cf(q)=f(Tq)` at a nonterminal state and
`Cf(q)=0` at a terminal one,

    ||f||_W²−||Cf||_W² = Σ_s μ(s)|f(s,0)|².

This follows simply by grouping the pullback norm by its target. Positive
predecessor weight implies positive target weight, so the weighted quotient
is now well-defined. The operator is contractive; it generally has norm
one, not a strict contraction. A coordinate at a positive-clock target
with positive incoming weight realizes equality.

Because injection occurs only at clock zero and there are no incoming
edges there, `w≥i+Kw` is equivalently `Kw≤w` together with the initial
lower bounds. That equivalence is false for arbitrary injection: on the
edge `0→1`, injection `(1,1)` has occupation `(1,2)`, whereas `(1,1)`
dominates the injection and satisfies `Kw≤w` without satisfying the sum
inequality. The included checker retains this hostile example.

If source `s` traverses `L_s` edges, conservation gives

    Σ_q W(q)=Σ_s μ(s)(L_s+1)≤(H+1)Σ_s μ(s),
    terminal mass = Σ_s μ(s).

These are nonnegative occupation identities. Coherent addition of signed
Eisenstein currents is a different quadratic operation and may increase
energy at a merged port. None of the identities proves the open all-scale
ROOT primitive-source debt bound or a universal finite arrival budget.
The small clock checker illustrates the theorem on a shared-target graph;
it does not replay the omitted native frontier audit.
