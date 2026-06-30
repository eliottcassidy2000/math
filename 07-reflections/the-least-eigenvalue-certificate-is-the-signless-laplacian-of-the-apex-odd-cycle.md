# The least-eigenvalue certificate is the signless Laplacian of the apex odd cycle — positive because 7 is odd

*klein-2026-06-29-S19. Asked to think the odd/even prior work and work on the least-eigenvalue certificate. The two threads turned out to be one: the certificate's positivity IS an odd/even (non-bipartiteness) statement.*

## The odd/even frame (what the prior work says)

There is one involution `sigma` wearing every costume (complement `T->T^op`, reversal, `t->-t`, the
half-translation `t->t+1/2`). Every invariant splits `sigma`-even (+) `sigma`-odd. mac-mini-S6 (THM-581/582)
settled which half the LRC floor lives in: **"lonely is even, Redei is odd."** The lonely measure is
`sigma`-invariant, so existence of a lonely time is a purely `sigma`-EVEN, Brouwer/SOS-category statement;
the `sigma`-odd index (Redei's `H(T)` odd = palindromic paths) does NOT apply to the floor (a lonely
tournament has the observer as a source, so it is not self-converse). Floor closure = the even-category
route: the 2-adic descent (THM-580) + the cyclotomic Fejer-Bochner SOS.

So the right certificate is a **`sigma`-even, Bochner-positive, sum-of-squares** object. That is a
least-eigenvalue of a PSD Gram. This is not a stylistic choice -- it is forced by which half of the
equivariant Euler characteristic the floor is.

## THM-590 IS a least eigenvalue

For a descended core `O subset Z_7`, the autocorrelation circulant `C(O)_{ij} = a((i-j) mod 7)`,
`a(d)=#{x in O: x+d in O}`, is real-symmetric PSD with eigenvalues `|Ohat(k)|^2` (Bochner = SOS). THM-590's
`g(O) = min_{k!=0}|Ohat(k)|^2` is EXACTLY the smallest eigenvalue of `C(O)` (verified, 0 mismatches over all
127 cores; the least eigenvector is the explicit SOS witness). So the apex obstruction is a genuine
least-eigenvalue certificate, not just a min-of-a-formula.

## The binding atom is an ODD cycle's signless Laplacian

The minimal positive gap is the **doublet** `O={0,1}`. Its autocorrelation is `a=[2,1,0,0,0,0,1]`, i.e.
> `C({0,1}) = 2I + A(C_7) = Q(C_7)`, the **signless Laplacian** of the 7-cycle.
(This integrates mac-mini-S36/HYP-3601, which identified the doublet's autocorrelation as `2I+A(C_7)` and
the binding atom as the even-graph `C_7`.) Hence
> `4cos^2(3pi/7) = 2 - 2cos(pi/7) = lambda_min(Q(C_7))`.
And the signless Laplacian's least eigenvalue is `0` IFF the graph is **bipartite**. `C_p` is bipartite iff
`p` is even. Therefore:
> **the apex obstruction `lambda_min(Q(C_p)) = 2 - 2cos(pi/p)` is positive precisely when `p` is ODD.**
Verified `p=3..14`: positive at every odd `p` (`1, .382, .198, .121, .081, .058`), exactly `0` at every even
`p`. The LRC(14) obstruction is positive **because the apex prime 7 is odd** -- the odd cycle `C_7` is
non-bipartite. This is mac-mini-S34's "the truth is the odd cycle" made quantitative: the odd cycle is
literally the graph whose signless-Laplacian gap is the positive apex atom.

## Why this is the whole odd/even story at once

The two requests collapse into one sentence. The certificate must be `sigma`-EVEN (SOS, because the floor is
the even half) -- and the value it certifies is positive iff the apex cycle is ODD (non-bipartite). The
even-category SOS certificate has, as its binding eigenvalue, the non-bipartiteness gap of an odd cycle. The
`14 = 2*7` split is exactly this:
- **2-part (even):** the `sigma` half-translation / 2-adic descent peels the even structure to the apex
  prime (THM-580, the even-category degree);
- **7-part (odd):** the apex cycle `C_7` is odd, so its signless-Laplacian least eigenvalue is positive
  (`4cos^2(3pi/7)`), and that is the SOS certificate's binding value.

If the apex prime were even, `C_p` would be bipartite, the signless Laplacian would have a zero mode, and
the certificate would vanish -- there would be no floor. The conjecture's truth at `N=14` rests on `7` being
an odd prime.

## The honest scope

Rigorous: `g(O) = lambda_min(C(O))` (the matrix identity), the doublet `= Q(C_7)`, and
`lambda_min(Q(C_p)) = 2-2cos(pi/p) > 0 <=> p` odd (all finite/exact, THM-590 + signless-Laplacian spectrum).
So the apex skeleton's certificate is fully pinned, and its positivity is the odd-cycle non-bipartiteness.
What remains (HYP-3599): the bridge from this apex `sigma`-even certificate to the full per-level `rho_j` is
genuine at the deeper levels but is the original existence question at the top level (where the measure
vanishes at the cusp `O=Z_7`, i.e. where the apex cycle degenerates). The certificate lives at the apex, not
the full grid (the danger does not factor).

See HYP-3604 (this), THM-590 (the gap = the least eigenvalue), HYP-3601 (mac-mini: the `C_7` atom),
HYP-3594 (mac-mini: the odd cycle), [[the-two-indices-redei-is-odd-lonely-is-even-half-tiling-is-the-quotient]]
(the `sigma`-even floor), [[the-bridge-splits-by-level-deeper-decorrelations-are-genuine-the-top-is-the-existence-question]],
THM-588 (the graph-spectral least eigenvalue).
