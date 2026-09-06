# Separated missing intervals remove three more connected-complement clocks

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** Exactly the three authorized clocks14904,14898,14892 are excluded. Removing them from the previously audited necessary set leaves **8198 clocks, maximum14886**. No other physical clock is newly scanned. LRC(14) remains OPEN.

The mechanism is a uniform near-resonance lemma for the danger intersection of a pair(1,q). It replaces independent interval minima by the relative positions of the intervals in which a grid point is missing. A strict endpoint condition is essential, even for a pair in the actual inert atlas.

The resulting physical corollary retains the general connected-complement scope. Let n be a primitive row of thirteen distinct positive integer speeds. Choose any six labels A and their seven-label complement B, and put t=gcd(A), g=gcd(B), V=A/t, U=B/g. Suppose the strict atlas graph on U is connected. Then n is weakly safe whenever t is outside the new necessary set; in particular,

`gcd(A)>14886 implies a phase of clearance at least1/14`.

Neither decoder-span equality, a height box, nor connectivity of A is required. The lower-dimensional six-component phase supplier remains **CITED**, exactly as in the inherited grid theorem.

[Independent proof and exact audit](continuing6_20260906_lrc_near_resonance_audit.md) passes.

## 1. Inheritance, hostile, and connection contract

The inputs, with their proof and citation scopes unchanged, are:

- [third_20260906_grid.md](third_20260906_grid.md) and its audit: marginal danger ceilings, exact overlap credit, and physical grid lifting. It inherits THM-2060, [THM-2060-crt-tail-coset-saturation](../../01-canon/theorems/THM-2060-crt-tail-coset-saturation.md), and THM-2064, [THM-2064-multitail-sheet-capacity-and-dyadic-seam](../../01-canon/theorems/THM-2064-multitail-sheet-capacity-and-dyadic-seam.md).
- [third_20260906_decoder_profile_graph.md](third_20260906_decoder_profile_graph.md): complete projected words and connectedness force some actual atlas edge with sheet gcd e<=6. It does not prescribe that edge.
- [third_20260906_grid_full_words.md](third_20260906_grid_full_words.md), [third_20260906_grid_refined.md](third_20260906_grid_refined.md), and [continuing5_20260906_lrc_clock16704_complete.md](continuing5_20260906_lrc_clock16704_complete.md): complete word costs and the current necessary set of8201 clocks, maximum14904.
- THM-3818, [THM-3818-scaled-inert-cubeclass-support-two-pair-packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md): the5855 coprime positive ratios p<q, with p+q<=356 and every prime factor of p+q congruent2 modulo3 with exponent at most2.
- The complete projected profile JSON of [overnight12_20260906_lrc_gcd_semigroup.md](overnight12_20260906_lrc_gcd_semigroup.md), pinned below.

Common-translate overlap counting is inherited, including the parallel [long_frontier2_sep06_lrc.md](long_frontier2_sep06_lrc.md) adaptive overlap work. This note does not claim a new general gluing algorithm. Targeted searches for the three clocks, near-resonance counts, missing gaps and the exact q355 consumer found no earlier version of the present separated-missing-interval lemma or three-clock exclusion. No external priority claim is made.

The board is: a compatible full word; a supplied actual edge; one common translate; open danger endpoints; missing-interval positions; and the physical lift. The corrected near misses are separately minimizing each tooth, or taking a pair's maximum cost from an incompatible word. The frozen16704 controls show those repairs are independent. Here every scalar survivor has the same ratio(1,355), so a geometric lemma becomes more useful than another location table.

The source-to-target map takes all placed teeth of the same pair to their closed absence intervals on the phase circle. It preserves the exact overlap count, including endpoint values. Replacing them by their lengths would destroy their simultaneous incidence. The needed sidecar is their center progression modulo1. The cheapest decisive test is the equality case where the first and last closed gaps touch. The resulting credit then maps to the physical grid with its multiplicity e; retaining the complete word supplies the compatible ceiling budget.

## 2. A general near-resonance lemma

Write D_v={x mod1: ||vx||<1/14}. Let q>=2 be an integer not divisible by14, let R=floor(q/14), and let s be a positive integer such that n=7q-s>0. Define the translated grid

`X_alpha={(alpha+j)/n mod1:0<=j<n}`.

If

`s(14R+1)<7q`,                                                   (1)

then, over all real alpha,

`min |X_alpha intersect D_1 intersect D_q|=2R`,

`max |X_alpha intersect D_1 intersect D_q|=2R+1`.                 (2)

This is an all-parameter elementary lemma, not an extrapolation from the three clocks.

**Proof.** Put q=14R+r,1<=r<=13. Inside the unwrapped D_1 interval, the intersection is exactly the2R+1 full teeth

`I_k=(k/q-1/(14q), k/q+1/(14q)), -R<=k<=R`.                     (3)

Their extreme endpoints satisfy14R+1<=q, while the next tooth does not have a nonempty intersection. Equality at a D_1 endpoint is harmless because both danger sets are open. The condition14 not dividing q is retained: in the divisible case the extreme teeth are clipped halves.

Each tooth has length1/(7q), so its n-fold image is an open circle arc of length1-s/(7q)<1. Consequently it contains either zero or one grid point. The phases at which it contains zero form a **closed** absence interval of length

`ell=s/(7q)`

centered at

`c_k=1/2+n*k/q=1/2-s*k/q mod1`.                                (4)

For R>=1, a difference d=|k-l| lies between1 and2R. Condition(1) implies0<sd/q<1. The circular center distance is therefore at least

`min(s/q,1-2Rs/q)>s/(7q)=ell`.

Thus the closed absence intervals are pairwise disjoint, with a positive gap between them. At any phase at most one tooth is absent, giving the lower bound2R. A phase at an absence center attains it. There are phases outside all the separated absence intervals, attaining2R+1. For R=0 there is one absence interval of length less than1, and the same conclusions follow directly. This proves(2).

The proof also supplies a reusable sufficient certificate from the actual positions, rather than a general equivalence between separation and a physical safe phase. No converse beyond the exact extrema under(1) is claimed.

## 3. Equality really loses one point

Condition(1) cannot be weakened to <=. For every R>=1, take

`q=14R+1, s=7, n=98R`.

The absence intervals have length1/q. The centers belonging to k=+R and k=-R are1/(2q) and1-1/(2q); their closed intervals touch at alpha=0. All other pairwise center distances exceed1/q. Hence exactly two teeth are absent at that phase, and no phase misses more than two:

`minimum=2R-1, maximum=2R+1`.                                  (5)

The maximum follows also from the absence-union measure(2R+1)/q<1. Thus the first failed implication is strict separation -> at most one absence when strictness has been discarded. An open-cell-only test misses the touching phase.

The smallest generic example in this family is q15,n98, with minimum1 and maximum3. That ratio(1,15) is not in the strict atlas. The stronger actual-ratio control is

`q43,n294, pair(1,43), minimum5 atalpha0, maximum7`.

Here p+q=44=2^2*11 satisfies the exact atlas. This is an actual allowed pair, not a claimed thirteen-speed failure or decoder realization. The nearby strict control q43,s6 has minimum6 and maximum7. The verifier reconstructs all three geometries and checks every critical wall and cell by literal modular danger predicates.

## 4. The complete word and candidate universes

For each authorized t, let F_t consist of seven-state multisets d with entries among the42 inherited allowed divisors of t, gcd(d)=1, and every proper nonempty positional subset I satisfying the full profile condition

`(c,sort(gcd(c,d_j):j outside I)) in profile level7-|I|`,

where c=gcd(d_i:i in I). Repeated states retain their positional multiplicity. There are126 proper positional subsets. These are necessary projected words in a hypothetical failure, not converse physical realizations.

The source enumerates every sorted multiset in each divisor domain without branch pruning. For each valid word it evaluates

`E(d)=sum_i d_i ceil(t/(7d_i))-t`.

One enumeration also updates the exact maximum M_(a,b) over every unordered pair of states occurring with the required multiplicity. The complete valid-word lists, costs, all margin types including infeasible ones, and every attaining owner are retained in the JSON. Completeness of the enumeration proves the upper bounds; the owners alone do not.

| t | Allowed divisor count | All multisets | Valid full words | Realized pair types | M(t) | M_(6,6) |
|---:|---:|---:|---:|---:|---:|---:|
|14904|16|170544|9593|127|113|74|
|14898|4|120|61|10|12|12|
|14892|9|6435|681|39|38|32|

The exact divisor domains are respectively

```
14904: 1,2,3,4,6,8,9,12,18,23,24,27,36,46,54,72
14898: 1,2,3,6
14892: 1,2,3,4,6,12,17,34,51
```

Every5855 atlas pair and every e<=6 dividing t is retained. Put n=t/e. The forced margins are e*gcd(n,p) and e*gcd(n,q). The complete clipped component-length numerators over14pq are2p once and min(2p,p+q-14k) twice for1<=k<(p+q)/14. Every list is independently reconstructed by literal interval intersection in the producer. The separate credit is

`C_sep=e sum_components(ceil(n*length)-1)`.

There are29275,23420,29275 base candidates, respectively. In each clock, **every candidate except(e,p,q)=(6,1,355) already satisfies C_sep>M(t)**. Thus81967 of81970 candidates are excluded at that stage.

The full pair-conditioned tables are retained but are not load-bearing for those scalar exclusions. The alternate conditional-first accounting is3976/25298/1,3720/19699/1,2679/26595/1, where the columns are infeasible, separate, located. Both orderings retain the complete candidate universe and agree that only one location calculation is needed at each clock.

## 5. One structural location certificate closes all three

For the sole remaining ratio, q355, R25, and

`n=t/6=2484,2483,2482=7*355-s`, with s=1,2,3.

Each satisfies s*351<2485. The lemma gives exactly50 to51 pair-danger points on every normalized grid, so the actual overlap credit is

`C_loc=6*50=300`.

The separate credit is zero: all51 separate teeth have n*length<1. Their zero minima cannot occur simultaneously.

The exact compatible(6,6) owners and surpluses are

| t | Owner | Cost | Located credit minus cost |
|---:|---|---:|---:|
|14904|(1,3,6,6,8,46,54)|74|226|
|14898|(2,2,2,3,3,6,6)|12|288|
|14892|(1,1,3,6,6,12,17)|32|268|

The weaker global-cost comparisons also suffice:300 exceeds113,12,38. For each of the three, the complete profile has102 projected endpoint walls; the verifier checks every wall and every open cell by both the interval count and literal modular predicates on all n grid points. These are exact controls for the analytic lemma, not a finite approximation to a minimum over real phases. Additional literal t-grid controls check the multiplicity6 map and nontrivial coprime multipliers.

## 6. Physical transfer and actual necessary-set deletion

Take a primitive seven-shape U with sheet word in F_t and connected strict atlas graph. The inherited graph certificate supplies some actual pair(Dp,Dq) with e=gcd(t,D)<=6. Its reduced ratio is in the complete atlas above. Writing D=e*h gives gcd(h,t/e)=1, which proves both the forced margins and the map of a translated t-grid to a translated(t/e)-grid with multiplicity e. The unit h changes the translate and permutes its points; it does not change the uniform minimum.

Each individual danger count is at most d_i ceil(t/(7d_i)). At a grid point with m dangerous labels, the selected pair indicator is at most(m-1)_+. Therefore a valid pair credit C implies at least C-E(d) weak-safe grid points. Every supplied candidate above has a strictly positive surplus. This proves the profile-conditioned grid supplier at the three clocks.

Now use the thirteen-speed row and chosen A,B from the opening statement. If no weak-safe phase existed, the inherited global joint-shadow profiles would hold; gcd(t,g)=gcd(V)=gcd(U)=1, and the sheet word of U would lie in F_t. The **CITED** proper six-component LRC input supplies a V-safe phase eta. Its t physical lifts(eta+j)/t preserve safety of A; multiplication by g permutes the translated t-grid in the U clock. The grid supplier gives a lift weak-safe for B too, a contradiction. Open danger endpoints yield weak safety, and no strict-clearance conclusion is inferred.

The current8201-element necessary set is read in full from the audited16704 certificate, with semantic pin

`ddbe0e091d36d54c8f6a7c8ea631bbf363d799b9adaa2ec4fe4cf56250d11a76`.

Only14904,14898,14892 are deleted. The source checks actual membership and actual set difference, obtaining8198 elements and maximum14886. The complete resulting array has semantic SHA256

`c146d73c149c448e744e138aa3bb6f8c286748be81353b8a6605adc47a129117`.

No adjacent-clock monotonicity is assumed. Membership in this necessary set does not assert an unsafe physical realization. The new maximum gives the corollary in the opening statement, with all inherited obligations unchanged for clocks that remain.

## 7. Reproduction, pins, and stopping boundary

```
python ../../04-computation/continuing6_20260906_lrc_near_resonance.py
python -O ../../04-computation/continuing6_20260906_lrc_near_resonance.py
```

Outside the repository, run the same filename from this directory. Filed paths prefer adjacent profile data and sibling results; the outside-worktree fallback is retained. No producer implementation is imported.

Both runs pass99305 always-active gates, with byte-identical actual LF output and JSON. Capture rejects carriage returns rather than normalizing them. The complete certificate includes10335 valid words and costs, every conditional type and owner, all81970 candidate decisions, exact wall/cell profiles, endpoint hostiles, native multiplicity controls, and the new necessary array.

- Source SHA256:`6027344aba5f9bda6c2637d7a992cb66367ec05f03c9cbde26761c65f438f93f`.
- Output SHA256:`9a650dd5227063c06f62c88bdd6036632761167e3e29475340e9a35bf6719f61`.
- JSON SHA256:`7692b34aa305701be0df279435d281ab498c680a7c1195c455f381b0a9f35661`.
- Inherited profile JSON SHA256:`935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.

After the independent proof and exact audit passed, the certificate status field was changed from a pending proof candidate to `FINITE-EXACT complete certificate`; this changed no arithmetic or producer logic. Both execution modes were rerun, with the same99305 gates and all mathematical data unchanged. The pins above identify that final data-status version. Parent owns the report's audited-status promotion and link.

The transferable new operation is to count simultaneous absences after the near-resonance projection. It succeeds because their relative centers separate; it fails sharply when the closed gaps touch. The complete three-clock consumer is the bounded consequence established here. No further clock is tested, no universal exclusion over the remaining8198 is asserted, and no new theorem ID is claimed. The independent audit above supplies status promotion; the manifest pins the filed artifacts.
