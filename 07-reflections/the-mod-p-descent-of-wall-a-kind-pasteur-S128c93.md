# The mod-p descent of Wall A: descend, prove in F_p, cage-lift back

**kind-pasteur-2026-07-19-S128c93** (HYP-8020) · owner: *"find an even more creative new
proof synthesis attack angle."*

## 0. The angle in one paragraph

Wall A — the inverse/rigidity complex (near-floor ⟹ AP/GW structure) — has resisted
every ℤ-side attack because heights are unbounded and the extremizers are arithmetic
(the ~45% mistake genus; the certificate-rung mechanism; the escape tails). **But Wall A
has a mod-p shadow: "every improper 13-tuple mod p is near an AP/GW dilate mod p" — and
three things changed this week that make the shadow the right place to attack.** (1) The
c92 census MEASURED the shadow to be essentially TRUE at p ~ 300–460 (minimal covers
collapse to ~10¹–10², all size exactly 13, within Hamming ≤ 3 of the ansatz at unsmooth
primes, with exact AP/GW dilates reappearing). (2) The mod-p statement lives in the ONE
arena where inverse-theorem technology actually exists — finite-field additive
combinatorics: sum-product, Bourgain–Glibichuk–Konyagin exponential-sum bounds,
expansion in groups — machinery built precisely for "small sets in F_p cannot have
approximate multiplicative-additive structure unless structured." The ℤ-side literature
gap (klein HYP-6100: "the literature covers bounded-max, NOT the spread gap") does not
exist mod p. (3) The HYP-7940 cage already knows how to LIFT per-prime structure back to
ℤ-height rigidity — that is literally what it does. So the pipeline is:

> **DESCEND** (Wall A → the mod-p shell theorem) → **PROVE IN F_p** (expansion /
> counting starves non-shell covers) → **CAGE-LIFT** (per-prime shell membership + the
> J-separator elimination ⟹ ℤ rigidity + micro-gap emptiness to explicit heights).

No route in the repo's 50-day history has attacked Wall A from inside a finite field.
The census makes it empirically safe; the cage makes it consequential.

## 1. The precise target (new, verified unposed)

**The mod-p shell theorem (target).** There exist k₀ and a class P* of primes (the
unsmooth/safe class) such that for p ∈ P*: every 13-tuple of residues mod p with no
witness time a/p (improper; equivalently its 13 folded danger-APs A(w) = {fold(jw) :
j ≤ ⌊p/14⌋} cover [1, (p−1)/2]) has folded speed-set within Hamming distance k₀ of a
dilate of {1..13} or {1..11,13,24} mod p.

Fragments already in canon (credited): the mod-13/19/23/25 spread lemmas
(klein/boxeph, kernel-pure) are NECESSARY CONDITIONS of exactly this statement — each
says an improper-regime tuple must antipodally cover — but none is an inverse (they
constrain, they do not classify). opus-S107's "sum-product coincidence" (the extremizer
= additive minimality ∩ residue completeness) is the ℤ-side slogan of the same
intersection; boxeph-S90's function-field port was the adjacent finite-model instinct
(F_q[t], different model, abandoned); the descent to literal F_p residues with a lift
mechanism is new.

## 2. Why F_p technology plausibly closes it (the mechanism, from the S6 face)

Covering [1, h] by 13 danger-APs A(w) = w·±[1, dk] (dilates of one interval) requires
overlap efficiency: waste = Σ|A| − |∪A| must land exactly at 13·dk − h (the S1
identity). Overlaps are RESONANCES: r(w, w′) = #{(j, j′) ∈ [1,dk]² : jw ≡ ±j′w′} —
the number of small-ratio representations of w′/w. Structured covers manufacture
resonance patterns; generic w-sets waste. At safe primes p = 2q+1 the multiplicative
group is C₂ × C_q — no small subgroups, so ratio-sets cannot cluster: the ratio graph
is expansion-constrained, resonances equidistribute, and waste-efficient 13-covers are
starved down to the shell. That is the S6 face made mechanistic — and "small
multiplicative structure forces additive spreading" is THE sum-product phenomenon.
The diagnostic below measures its signature (resonance-tail concentration, safe vs
smooth). The proof-shaped question for the descent: **bound the number of
waste-efficient 13-covers via resonance-moment counting (2nd/4th moments of r over the
chosen pairs), with BGK-type input at the point where a cover's ratio-multiset must
concentrate.** This is a counting argument in a finite group — bounded, checkable
per p, and of exactly the shape finite-field combinatorics proves.

## 3. The lift (what the cage already does, and the one named upgrade)

HYP-7940's machinery: a primitive family with M below the rung at p is improper mod p;
shell membership mod each caging prime + the J-separator two-stage elimination + the
Newton tower pin V to the terminal families up to explicit heights (281,577 at the
{1,2,4,8}-grid; klein's sandwich already consumes this as its bottom layer). With
EXACT-dilate shells the lift is built. With Hamming-k₀ shells the elimination needs the
**bounded-defect Newton upgrade** (power sums matching up to k₀ substitutions): the
defect positions vary per prime, so either (i) a k₀-defect version of the R_m
elimination (structured: defects contribute O(H^{2m}) correction terms with bounded
support — the same shape opus's 11-even branch handled), or (ii) prime-majority
voting (defects at ≥ (1 − ε) of caging primes force ℤ-defect structure — pigeonhole
over defect patterns, C(13, k₀)-bounded). Named obligation, not hand-waved: this is
the pipeline's step-2-to-step-3 glue, one focused session to draft.

## 4. The computational twin: the safe-prime sieve

Independent of the analytic descent, the S6 face redesigns the S–T finite check:
**run the sieve preferentially at unsmooth primes, where the I(13,p,1) grind is
10–100× thinner** (measured: draw-success 3/800 at 457 vs 111/800 at 307; the c92
abundance table). The S–T architecture never chose primes by sea-thinness — prime
selection was availability-ordered. Diagnostics this session: (a) whether the caging
budget (Σ ln p > 657.6) is satisfiable inside the safe/unsmooth classes alone; (b) sea
samples at safe {719, 839, 983} vs smooth controls {631, 991, 1009}; (c) the
resonance-concentration statistic. **Diagnostic results (frozen out):** (a) BUDGET FEASIBLE IN-CLASS: 85 safe primes
(227..7523, sum ln = 665.5 > 657.6) or 89 unsmooth primes (227..4517) each carry the
full caging product alone. Cost caveat stated honestly: the S-T p^((k+1)/2) lifting
heuristic explodes at the class's large tail (7523^7), BUT that heuristic models
lifting of SURVIVORS and does not apply to near-empty seas -- the safe-prime sieve's
true cost model must be rebuilt around sea size, not p-size (named follow-up). (b) THE
REFINEMENT THAT SHARPENS THE TARGET: at p ~ 700-1000 the greedy-reachable sea is PURE
SHELL at every prime -- safe primes {719, 839, 983} sample ONLY distance-0 ansatz
dilates (1-3 distinct covers per 400 draws, all exact AP/GW dilates); smooth controls
{631, 991, 1009} nearly so (991: distances {0:2, 1:2, 2:1}). The collapse completes
universally by ~700; smoothness modulates the RATE in the 200-600 transition band (the
c92 contrast). Consequence: **pose the mod-p shell theorem at p >= ~700 with k0 <= 2,
and put the caging set there** -- the empirical k0 is far smaller than feared, which
keeps the bounded-defect lift upgrade light. (c) HONEST NULL: the random-pair
resonance statistic does NOT distinguish safe from smooth (tail-ratio ~3.5x mean
everywhere; variances mixed) -- the S6 mechanism's signature lives in the STRUCTURED
pairs that covers actually use, not in random pairs; the right instrument is the
internal r-distribution of found covers (one more session). The S6 correlation itself
(greedy difficulty vs unsmoothness) stands as measured; its pairwise-tail MECHANISM is
not confirmed by this instrument.

## 5. Convergences logged while claiming (the fleet's same-hour work)

opus-S407 (HYP-8005) landed the **tournament regularity bridge**: the half-turn
observer tournament's 3-cycle count c₃ separates AP (112) from GW (108) at the floor —
a TOURNAMENT-side separator, arriving hours after the moment-side J-separator
(HYP-7940) and the profile-side identity window (MISTAKE-093). Three separators, three
instruments, one pair of families: the sliver's owners now form their own face-family
(J ≡ c₃ ≡ profile-kink as AP/GW discriminants — worth one canonical statement). And
their ghost result ("extremality = forbidding your own best rational approximation" —
the excluded element is the penultimate convergent) is exactly the kind of per-prime
structural fact the descent's F_p proof would want as an anchor.

## 6. Honest obstacles

(1) The census measured the GREEDY-REACHABLE shell; the mod-p theorem needs ALL covers
— the F_p proof must handle what sampling cannot see (that is its job; the counting
frame of §2 is reachability-free). (2) k₀ might grow with p (census says ≤ 3–6 at
tested p; if unbounded, the lift's defect machinery fails — the FIRST thing the F_p
analysis must pin). (3) BGK-type bounds have absolute constants that may be too weak
at 13 sets of density 1/14 — the counting argument may need the specific interval
structure (Fourier on ±[1,dk] is explicit — sin kernels — which helps). (4) The cage
lift inherits its conditional character until the mod-p theorem replaces the un-run
sieve's terminal-class hypothesis — that replacement is the point: **the descent
would make the cage UNCONDITIONAL at its caging primes.**

## 7. Handoffs

- **anyone with finite-field appetite (opus's pinning genre):** the §2 counting
  question is one literature session (BGK/sum-product statements at interval dilates
  in F_p) + one drafting session; the census + diagnostics supply the constants to
  aim for.
- **boxeph:** the resonance-moment counting is your cover-debt frame in F_p; the
  spread lemmas you proved are its first-order necessary conditions — the inverse is
  the same machinery run backwards.
- **death-star:** the safe-prime sieve twin needs your ghost evaluator at scale; the
  budget table (part a) says whether prime selection alone re-prices the finite check.
- **opus:** your c₃ separator + my J-separator + the profile window = the separator
  face-family; and your ghost/convergent result is a candidate F_p anchor lemma.
- **mac-mini:** your small-witness law (S125, recommended to me) and this descent are
  siblings — both are "denominator-bounded verification with an lcm-channel residue";
  the descent's F_p counting may supply your (a)-part's 78 BAD_q analyses wholesale.

**Files:** diagnostics script + frozen out (`lrc14_modp_descent_diagnostics_kps_S128c93`),
HYP-8020, SESSION-LOG entry, backlog lead.
