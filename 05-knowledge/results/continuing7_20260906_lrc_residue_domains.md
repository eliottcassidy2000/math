# Exact divisor-domain costs remove three further clocks

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The three newly tested clocks14886,14880,14874 are all excluded. Their actual deletion leaves **8195 necessary clocks, maximum14868**. No other whole clock is newly tested. LRC(14) remains OPEN.

There is also a uniform selected-edge consequence: the supplied edge(e,p,q)=(6,1,355) can never be the remaining obstruction at any clock t>=14868. This does not exclude every other edge at an untested clock, and in particular does not establish whole-clock closure at14868.

The physical conclusion retains the general scope of the inherited theorem. In a primitive row of thirteen distinct positive speeds, choose any six-label subset A with complement B. Put t=gcd(A),g=gcd(B),V=A/t,U=B/g. If the strict atlas graph on U is connected, then t outside the new necessary set implies a phase of clearance at least1/14. In particular,

`gcd(A)>14868 implies weak full-row safety`.

No decoder-span equality, physical height bound, or connectivity of A is added. Proper six-component LRC remains the same **CITED** input as in the inherited grid theorem.

[Independent proof and exact audit](continuing7_20260906_lrc_residue_domains_audit.md) passes.

## 1. Inheritance and the two arithmetic coordinates

The proof inherits:

- `third_20260906_grid.md` and its independent audit: marginal ceilings, pair credit, multiplicity excess and physical grid lifting; its underlying exact marginal input includes THM-2060, `01-canon/theorems/THM-2060-crt-tail-coset-saturation.md`.
- `third_20260906_decoder_profile_graph.md`: full projected words and connectedness supply some actual edge with sheet gcd e<=6, without prescribing the edge.
- THM-3818, `01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`: all5855 coprime ratios p<q with p+q<=356 and every prime factor of p+q congruent2 modulo3 with exponent at most2.
- `third_20260906_grid_full_words.md`: complete projected words, the residue-cost identity, and the sharp universal full-word caps210,270,192,239,197,224 in residues1 through6. Its proof of the upper bounds, not just its displayed owners, is used below.
- `continuing5_20260906_lrc_clock16704_complete.md`: pair-conditioned costs attached to one compatible full word, with multiplicity retained.
- `continuing6_20260906_lrc_near_resonance.md`: the strict missing-interval lemma, its equality hostile, and the full8198-element necessary set, maximum14886.

The board is: the exact clock divisor domain; its residue modulo7; the same full word and forced margins; one actual atlas edge; relative missing-interval positions; and the physical lift. The canonical endpoint hostile remains q43,n294: two closed absence gaps touch and its minimum drops from6 to5. The corrected cost near miss is keeping the residue while forgetting the divisor domain. The least-used sidecar in this cycle is the full family of costs on one fixed domain, rather than one maximizing word at one clock.

Targeted searches and the full-word source revealed that the residue mechanism and modulus5388292800 are already explicit in the inherited theorem. They are credited as inputs. The new finite tables below couple all residues to the three exact domains, and the uniform edge argument combines the inherited global270 cap with near-resonance positions and the forced-margin filter. No new general gluing algorithm is claimed.

## 2. Exact cost transport on a fixed divisor domain

Let A_0 be the42 inherited allowed sheet states, with

`Lambda=lcm(A_0)=5388292800`.

Every state is coprime to7. For a clock t define

`D(t)={d in A_0:d divides t}={d in A_0:d divides gcd(t,Lambda)}`.

The full-word set F_D consists of sorted seven-state multisets with entries in D, gcd1, and all126 proper nonempty positional subsets I satisfying

`(c,sort(gcd(c,d_j):j outside I)) in inherited profile level7-|I|`,

where c=gcd(d_i:i in I). Repetitions are allowed and positional multiplicities are retained. These are necessary projected words, not converse physical realizations.

For a fixed D, this set is independent of the clock. If tau=t mod7 and every state divides t, its exact ceiling excess is

`E_tau(d)=(1/7)sum_i d_i*((-tau*d_i^(-1)) mod7)`.

The numerator is an integer multiple of7 because there are seven summands, each congruent-tau modulo7. Thus the entire global and pair-conditioned cost functions depend only on D and tau. Domain inclusion gives inclusion of feasible full words and cannot decrease the maximum; keeping tau alone does not determine the budget.

One complete enumeration of each domain updates every pair-conditioned maximum at all seven residues. A repeated margin(a,a) is inserted only when two copies occur. The all-residue table contains every feasible type and an attaining owner; absence from a table means no full word has those margins.

The three domains are exactly those with gcd(t,Lambda)=18,480,6:

```
D18={1,2,3,6,9,18},
D480={1,2,3,4,5,6,8,10,12,15,16,20,24,30,32,40,48,60},
D6={1,2,3,6}.
```

Their exact maxima are:

| Domain | tau0 | tau1 | tau2 | tau3 | tau4 | tau5 | tau6 |
|---|---:|---:|---:|---:|---:|---:|---:|
|D18|0|33|30|12|44|27|24|
|D480|0|94|88|68|116|95|92|
|D6|0|6|12|12|14|15|20|

These tables apply at every height whenever the clock has the declared exact domain and residue. The theorem is about the cost quotient; it does not itself provide a safe phase at all such clocks.

The domain sidecar is necessary. At tau5 the small domain cap is15, but the actual14880-domain full word(1,3,4,6,32,48,60) costs95 and passes every inherited full profile. The failed implication is transporting a budget across a change of domain merely because the residue is unchanged. This full-word witness is not an unsafe physical row.

## 3. Complete enumeration at the three permitted clocks

The clocks were selected from the actual pinned necessary array, in descending order. The source exhausts every sorted seven-state multiset in each divisor domain, without branch pruning. Every accepted word is tested against all126 profile conditions. The complete accepted lists and costs are retained.

| Clock | All multisets | Valid full words | Realized unordered pair types | Global maximum |
|---:|---:|---:|---:|---:|
|14886|792|221|21|44|
|14880|346104|17826|158|95|
|14874|120|61|10|20|

There are347016 multisets and18108 valid words in total. For each valid word the native ceiling formula is checked against the residue formula. Every conditional owner is checked with its multiplicity; exhaustive enumeration, rather than those owner checks alone, supplies the upper bound. All seven residue tables are generated during these same enumerations; no extra physical clocks are searched.

For every5855 atlas pair and each e<=6 dividing the clock, put n=t/e. The forced margins are

`d_p=e*gcd(n,p), d_q=e*gcd(n,q)`.

There are23420,35130,23420 candidates, respectively, for81970 in total. In particular e5 is retained at14880. The exact clipped component lengths are reconstructed from literal pair-danger intersections. Their separate credit is

`C_sep=e sum_components(ceil(n*length)-1)`.

Every candidate except(6,1,355) has C_sep greater than the clock's global full-word maximum. The least noncritical credits are150,148,150, respectively. They even exceed the maximum over all seven residues on their fixed domains:44,116,20. These statements follow from the complete atlas computations, not from examining the preferred edge.

An alternative conditional-first accounting gives:

| Clock | Infeasible margins | Separate credit closes | Needs location |
|---:|---:|---:|---:|
|14886|0|23419|1|
|14880|5557|29572|1|
|14874|1848|21571|1|

Every candidate and decision is recorded. There is no survivor after the following located test.

## 4. The remaining pair, with its changed forced margins

For p1,q355,e6 the normalized grid lengths are n2481,2480,2479. These equal7q-s with s4,5,6. The inherited lemma applies because R=floor(q/14)=25 and s*351<2485. The51 full pair-danger teeth have exact normalized overlap minimum50 and maximum51. Consequently

`C_loc=6*50=300`.

Every individual tooth has separate credit zero. Their individual worst translates cannot be imposed simultaneously. Each complete102-wall profile is checked at every endpoint and every open cell, including a second path using the literal modular predicates on all grid points.

The actual forced margins and compatible costs are:

| Clock | Forced margins | Conditional maximum | Attaining full word | Located surplus |
|---:|---|---:|---|---:|
|14886|(6,6)|32|(1,2,2,6,6,9,18)|268|
|14880|(6,30)|77|(1,3,6,20,30,32,48)|223|
|14874|(6,6)|20|(1,2,2,3,6,6,6)|280|

At14880, n2480 is divisible by5, so the second margin changes to30. Treating it as6 would attach the wrong owner. The weaker global comparisons300>44,95,20 also suffice for these three clocks, so the conditioned tables are retained without claiming they are necessary for this particular exclusion.

## 5. A uniform selected-edge exclusion, including exact resonance

The inherited full-word theorem gives E<=270 for every clock; in residue0 it gives E=0. Combine this global cost with the same actual edge(e,p,q)=(6,1,355). Suppose t>=14868 and e divides t, so n=t/6>=2478.

- If2478<=n<=2484, then n=2485-s with1<=s<=7. The strict inequality s*351<2485 holds even at s7. Hence C_loc=300>270.
- If n=2485, the normalized pair really can have zero intersection count: all absence points coincide at alpha1/2. However, its forced q margin is6*gcd(2485,355)=2130, outside the entire allowed state set. Thus this formal pair cannot be the actual supplied edge of a full-profile row. A false positive count is not used to bridge the resonance.
- If n>=2486, each of the51 full teeth has length1/2485 and separate count at least1. Therefore C_sep>=6*51=306>270.

This proves the selected-edge exclusion at every t>=14868. The only arithmetic boundary that does not give overlap is removed by its exact forced margin. The local resonance check is not a complete test of another clock. No other edge at an untested clock is excluded by this argument.

The map here retains three different coordinates: one common translate near resonance, native divisibility at exact resonance, and monotone separate counts above resonance. None may substitute for another without its condition. This is the uniform consequence extracted from the finite positive signal.

## 6. Physical transfer and the new necessary set

For a primitive seven-shape U with word in F_t and connected strict graph, the inherited graph theorem supplies an actual pair(Dp,Dq) with e=gcd(t,D)<=6. Writing D=e*h gives gcd(h,t/e)=1. This proves the forced-margin formulas and sends any translated t-grid to a translated(t/e)-grid with multiplicity e. The unit h changes the phase and permutes the reduced grid. The verifier includes literal multiplicity controls.

The sum of individual danger counts is at most t+E(d). At a grid point with m dangerous labels, the selected pair indicator is at most(m-1)_+. Therefore a uniform pair credit C leaves at least C-E(d) weak-safe grid points. Every candidate at the three tested clocks has positive surplus.

For the thirteen-speed row in the opening statement, hypothetical failure supplies the global inherited profiles. Primitivity gives gcd(t,g)=gcd(V)=gcd(U)=1. The **CITED** proper six-component theorem supplies a V-safe phase eta. Its lifts(eta+j)/t preserve A safety, while multiplication by g permutes the t-grid in the U clock. The profile-conditioned supplier gives a lift weak-safe for B, contradicting failure. This argument uses neither decoder equality nor a height box, and it does not infer strict clearance from weak endpoints.

The old8198-element array is pinned by

`c146d73c149c448e744e138aa3bb6f8c286748be81353b8a6605adc47a129117`.

Removing exactly14886,14880,14874 leaves8195 elements, with maximum14868. The complete new array is stored and has semantic SHA256

`c67f5e98eff3fe406a4416057c6063095290330a50f039e731ced0fc2ca4657a`.

The source verifies the actual difference and maximum, rather than assuming adjacency of clocks. Membership remains necessary, not evidence of a realizable failure. In particular the uniform selected-edge result above must not be mistaken for closure of the new maximum clock.

## 7. Reproduction and frozen evidence

```
python ../../04-computation/continuing7_20260906_lrc_residue_domains.py
python -O ../../04-computation/continuing7_20260906_lrc_residue_domains.py
```

The outside version runs by filename from this directory. Filed paths use adjacent profile data and sibling results, with the outside-worktree fallback retained. The inherited universal-cap output is read and its exact six maxima checked; its independent optimization proof remains the upper-bound supplier. No producer implementation is imported.

Both execution modes pass253341 always-active gates with identical actual LF output and JSON. The source was adapted from our preceding producer; the separate referee supplies the independent path. The JSON retains every valid word and cost, all seven conditional tables per domain, all81970 edge/sheet decisions, literal located controls, inherited cap values, the uniform selected-edge parameters, and the full necessary array.

- Source SHA256:`1afc8c3438c9db2ad6bad53b771d84a9b6e46d4f6464722aa4ce240d73084e1b`.
- Output SHA256:`f7554bdb450c2f6a042c0caef9270eba4832934682d72b1c8cfaaff5db223858`.
- JSON SHA256:`4764d3def5d24702465f702180bad715b14402e9c45ddcd72ae1968e676a8086`.
- Full inherited profile JSON:`935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.
- Six inherited global-cap values, semantic SHA256:`6d57ff60ed6a7ed80af2fb8c1efd854eda75ec8131bbab82f41b8701ac89059f`.

The residue operation is inherited; its exact restricted-domain tables and the all-height selected-edge consumer are the reusable objects recorded here. No fourth whole-clock test, repository edit, Git mutation, or scarce theorem-ID claim was made. The independent audit above supports the filed status.
