# The anatomy of a tight runner set — eight conditions and the animal they describe

*kind-pasteur-2026-06-22-S38. A free-ranging dissection of what AP {1,…,13} and Goddyn–Wong
{1,…,11,13,24} secretly have in common — and what any tight set of their kind must be.*

The owner asked: get abstract and convoluted; find the shared DNA. Here is the dissection. Each
layer is a *necessary condition* that AP and GW both pass and the loose near-misses fail. Stacked,
they nearly pin the animal down.

## The eight layers

**1. The divisibility threshold is exactly 14.** Define `q(S) = min{d : S has no multiple of d}`.
Then `M(S) ≥ 1/q(S)` — this is the q-witness lemma re-read as a threshold. A tight set has
`q(S) = 14`: it is **13-covering** (a multiple of every `d ≤ 13`) and **14-avoiding** (no multiple
of 14). The loose look-alikes give themselves away here: `{1..11,13,26}` has the *same residues mod
14 as the AP* but no multiple of 12, so `q = 12` and `M = 1/12`. Residues lie; divisibility tells the
truth. (The danger zone `M < 1/14` is `q > 14` — exactly the covering sets. The threshold *is* the
covering reduction.)

**2. It is a perfect one-hole tiling of ℤ/14.** At the optimum `t* = a/14`, the AP residues are
`{1,…,13} = ℤ/14 ∖ {0}`: the lonely runner sits in the unique empty residue. GW is the same with one
residue doubled and one (12) vacated. So a tight set is a *near-perfect tiling of ℤ/14 with exactly
one structural hole* — the observer's seat. (Fuglede whispers here: the unsafe arcs tile [0,1).)

**3. It covers the ±units.** For every unit `a` mod 14, some runner has `s·a ≡ ±1` — else `t = a/14`
would leave the observer lonelier than 1/14. (mac-mini's necessary condition.)

**4. Its phase gaps obey Steinhaus.** At the optimum the residues form a ≤3-gap configuration on
ℤ/14; `g(14) ≤ 3` is the open rigidity. (mac-mini.)

**5. It is reachable from the AP by valid doublings.** Among single-swaps, tight ⟺ the swap replaces
`v` by `2v` (the **doubling map** `x ↦ 2x mod 14`, never `3v`) and the swap is *valid*. GW is AP with
12 ↦ 24. Double-doublings are never tight. The AP is the root; tight sets are its doubling-children.

**6. The doubling is gated by the Jacobsthal function.** A site `v` admits a valid doubling iff the
window `[14−v, 27−2v]` lies entirely inside a **coprime-gap of `v`** (an interval between consecutive
integers coprime to `v`). For `v = 12` the window `[2,3]` sits in the gap `(1,5)` between coprimes 1
and 5 of 12. *Only `v = 12` passes* — so the census is `{AP, GW}`, and the Jacobsthal function is the
gatekeeper. The acceleration site is "the highly-composite number just below `n−1`," divisible by the
small primes the window contains.

**7. It is rigid against its Farey neighbor.** The denominator 41 has haunted this project — the
bounded-denominator worst case `{1..11,13,84}` (witness `D = 41`), and the nearest near-miss
`M = 3/41`. The reason is one determinant: `det[[1,3],[14,41]] = −1`, so **3/41 is the Farey neighbor
of 1/14** — the first Stern-Brocot child where a competing optimum can hide. Tightness demands that no
runner-gap opens at 3/41. 41 is the *second resonance* of n = 14, Farey-adjacent to the apex.

**8. It factors as 2-adic ⊗ 7-adic.** Everything above splits along `14 = 2·7`. The **2-adic** factor
is the doubling map `x ↦ 2x` (layers 5–6, the acceleration moves). The **7-adic** factor is the apex:
`7` is the *unique residue whose double is 0 mod 14`, sitting at the edge of the doubling map — the
pre-image of the observer's hole. The ±units, the binding pair, the sectors all live in the 7-adic
half. The Cayley–Dickson split of the project's foundation governs the tight locus too.

## The animal

Read top to bottom, the eight layers describe a single object:

> **A tight runner set is the canonical perfect one-hole tiling of ℤ/14 (the AP, Dirichlet-extremal),
> together with its orbit under the Goddyn–Wong doubling operad `x ↦ 2x`, where each doubling move is
> admitted only when the Jacobsthal function leaves a coprime-gap wide enough — and stays rigid against
> its Farey neighbors. For `14 = 2·7`, the operad's gate opens at exactly one site (`v = 12`), so the
> orbit is two points: `{AP, GW}`.**

This is why the census is so small and so stubborn. The AP is forced (it is the only perfect tiling).
The single exception GW exists because 12 — uniquely divisible by both primes in its Jacobsthal window
— lets the 2-adic doubling fire once without tearing the tiling. Every other site is gated shut. The
"convoluted" conditions are not independent decorations; they are the same fact (the AP-tiling is
2-adically rigid except at the Jacobsthal-gated site) seen through divisibility, tiling, doubling,
coprimality, and Farey geometry in turn.

## Why this does not (yet) finish LRC(14)

The eight layers characterize the *bounded* census — single doublings of the AP. The open core is
whether some *unbounded or exotically-structured* set sneaks past all eight: a set that is
13-covering/14-avoiding, ±unit-covering, ≤3-gap, Farey-rigid, yet *not* a doubling-child of the AP.
mac-mini's `g(14) ≤ 3` is the same wall from the Steinhaus side. The animal is well-dissected; proving
no *other* animal of its kind exists is exactly LRC for 13 runners.

## The pointer beyond

If the tight locus is "the AP-tiling under a Jacobsthal-gated doubling operad," then the right
generalization is functorial: for each `n`, the census is the AP plus the doubling-orbit admitted by
the Jacobsthal gate of `n`. The conjecture to chase — **the tight-set census of LRC(n) is computed by
the Jacobsthal function of the acceleration sites of `n`**, with the prime-factorization of `n`
choosing which adic doublings are available (2-adic for even `n`, generally `p`-adic acceleration
`x ↦ px` for each `p | n`). The lonely-runner extremal problem would then be a *Jacobsthal-gated
operad census* — and the apex prime `7` of `14`, the pre-image of the hole, is why the gate is so
narrow. The runners have been tiling a clock built from `2 × 7` all along; the tight ones are the
tilings the doubling map cannot improve.

— Related: [[lrc14-thread]], HYP-2917 (q(S)=14), HYP-2918 (doubling operad + gcd-window),
HYP-2919 (Jacobsthal / Farey-41 / tiling), HYP-2914 (Dirichlet-extremal), HYP-2893 (GW = Jacobsthal
acceleration), mac-mini S53 (g(n)≤3); `the-lonely-runner-is-a-random-round-tournament.md`,
`everything-is-the-triangle.md` (the 2·7 Cayley–Dickson foundation).
