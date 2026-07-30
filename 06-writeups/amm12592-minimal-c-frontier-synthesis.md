# The AMM 12592 minimal-C frontier: what is now known (death-star, 2026-07-30)

> **Living synthesis — refreshed 2026-07-30 during the coinC session.** Canon
> truth lives in THM-2160, THM-2225, THM-2966, THM-2967 and HYP-9061; this
> document is the connective narrative. Lane C/D/E results are integrated as
> they land.

## 0. The problem

A coin shows `0` with unknown probability `p in (0,1)`. The critical value
`n` of a nonconstant flip stream is the length of its maximal constant
initial run. AMM Problem 12592 (Elliot Glazer, 2026) asks for exactly fair
deterministic extractors whose pathwise stop time is at most `2n` (part a)
and `max(2, 2n-1)` (part b); repo canon THM-2225 proves both, THM-2160
sharpens to `T(1)=2, T(2)=3, T(n) <= 2n-2` for `n >= 3`. The frontier
question (HYP-9061):

> **(Q)** What is `C* = inf{C : exists D and a fair extractor with
> T(n) <= Cn + D}`?

All classical schemes give `C = 2`. Before this session nothing below 2 was
known impossible beyond the trivial `T(n) >= n+1`.

## 1. The spine normal form (THM-2966, PROVED)

A `T`-deadline fair extractor exists **iff** integer vectors
`0 <= w_{m,k} <= binom(d_m,k)`, `d_m = T(m)-m-1`, satisfy

```text
sum_{m>=1} p^m q W_m(p) + sum_{m>=1} q^m p V_m(p) = 1/2   on (0,1),   (S)
```

with `W_m` the depth-`d_m` decided-tree polynomial with coefficients
`w_{m,k}`. Hence `C* = 1 + gamma*` where `gamma*` is the minimal linear
degree-growth rate for which (S) is solvable. Everything below is a study
of (S).

## 2. Sublinear excess is impossible (THM-2967, PROVED, audit finalizing)

**No fair extractor satisfies `T(n) = n + o(n)`.** Mechanism: with
`d_m = o(m)`, the window series `F(p) = sum p^m (1-p) W_m(p)` has integer
coefficients and radius 1; the identity `F(p) + G(1-p) = 1/2` continues `F`
across `p=1`; Polya-Carlson forces rationality; Fatou + Kronecker force all
poles to be roots of unity; the `z -> 1-z` symmetry then pins them to the
primitive sixth roots (`|z| = |1-z| = 1`); and integrality of the resulting
quasi-polynomial coefficient tails collapses the identity to
`A(p) + B(1-p) = 1/2` with `A, B` integer polynomials — an integer equals
1/2 at `p = 0`. The method is rate-free: it says nothing about
`d_m ~ gamma m`, so it kills `C = 1 + o(1)` envelopes without bounding
`C*` away from 1.

## 3. The carry mechanics of (S)

Writing `W_m = 1/2 + Delta_m` turns (S) into a zero-sum deficit flow.
Facts (THM-2966 corollaries + lane A verification on THM-2225):

1. Cells with odd budget `binom(d_m,k)` carry forced half-integer deficits
   (Lucas). Dyadic-tail shells confine the parity failure to one **corner
   word** `0^n 1^l` per shell: a single `+-1/2` at monomial `(n, l)`.
2. **Ratio 2 is self-annihilating:** at `l = n` the 0-side and 1-side
   corners share the monomial `(n,n)` and cancel — THM-2160's middle pair,
   explaining why every classical scheme has `C = 2`.
3. **Dyadic split-jump:** over `F_2`, `(1+u)^N = (1+u)^{N-2^a}(1+u^{2^a})`:
   a half-deficit at `(z,o)` pushed `2^a` levels leaves halves at
   `(z+2^a, o)` and `(z, o+2^a)`, everything interior being integer —
   absorbable only where integer capacity exists.
4. **The capacity desert:** at rate `gamma < 1`, capacity lives in the two
   cones `o <= gamma z + O(D)` and `o >= z/gamma - O(D)`; between them no
   cell exists. Near-cone-edge capacity is polynomial (edge corridor width
   `O(D)`), deep-cone capacity exponential.
5. **The C=2 fixed point (exact lemma):** the naive single-jump routing of
   consecutive corners costs additive slack exactly `D = l` at tail scale
   `l` (frozen referee), i.e. envelope `T(n) = (1+gamma)n + (1-gamma)n =
   2n` at every `gamma` — shell gain and routing cost cancel identically.
6. **The two-bias certificate:** the externally supplied gate
   `(2457/6592) log(8847357/2974400) - log(1285/896) > 1/25` (decoded and
   exactly refereed this session) matches the shape of the desert's
   entropy max-flow condition; its biases are the Legendre duals of lattice
   rays of slopes `q_A/p_A ~ 0.697` and `q_B/p_B ~ 0.336`, one in the
   desert, one just inside the 0-cone edge (`2457/6592 ~ 0.3727`). Whether
   it certifies a lower bound (plausibly `C* >= 9049/6592 ~ 1.3727`) or
   gates a construction remains open: both readings are recorded in
   HYP-9061 with the lower-bound reading currently favored.

## 4. Open tension

Hand analysis of smarter routings (along-the-cone-edge hopping with
same-side opposite-sign corner pairing, debris merging, escape cascades
along the capacity-1 edge corridor) neither closes nor excludes rates
arbitrarily close to `gamma = 0`: the accounting requires the exact LP that
the certificate presumably dualizes. The decisive finite experiments are
lane C (exact truncated feasibility of `T(4)=5`-type envelopes with
rational Farkas certificates) and lane D (exact ledger simulation of
corner routing at `gamma in {1/2, 0.6, 0.75, 0.9}`).

## 5. Lane results

*(integrated as they land)*

- **Lane A (DONE, VERIFIED-EXACT):** checksum spine ledger closes
  intra-block; the only cross-spine flux is the boundary pair; budget binds
  exactly at block openers `m = 2^k`; `C = 2` is pinned by openers.
  Artifacts: `04-computation/amm12592_checksum_spine_ledger_laneA_deathstar.py`.
- **Lane B (DONE, referee 12/12; draft secured):** THM-2967 above.
- **Lane C (pending):** small-`n` envelope frontier.
- **Lane D (pending):** carry-ledger sustainability at `gamma < 1`.
- **Lane E (pending):** LRC(14) mixed-sector audit (anchor duty).

## 6. Classical context

von Neumann (1951) extraction, Hoeffding–Simons, Stout–Warren tree
procedures, and Mitzenmacher's survey optimize expected flips; the
critical-run pathwise deadline is new to AMM 12592. The parity mechanics
here share the char-2 Pascal/Hasse algebra of the repo's LRC thirteen-sheet
fibre (THM-2160 §7, THM-2201) — the same `F_2[C_h] = F_2[eps]/(eps^h)`
triangularization, used here for shells, there for owner reconstruction.
