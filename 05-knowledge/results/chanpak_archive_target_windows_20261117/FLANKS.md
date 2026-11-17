# Parity-balanced separators and an exact two-axis unit repair

Sampled eight old reflections uniformly without replacement from the sorted local reflection list, seed 861774. The algebraic JC2 finite-place obstruction, archived LRC mobile-centre and harmonic-support calculation, relaxed SOS Toeplitz test, and the explicitly refuted cusp analogy all press the same methodological point: keep the measured sidecar, not the suggestive scalar. Here the relevant harmonic is just the signed parity of a selected path slot, and its vanishing can be proved exactly. See the reproducible sample list below.

## General marked interval law

Let the tournament obtained on deleting a vertex v be the ordinal join of nonempty blocks D_i (they need not be strong). Set h_i=H(D_i), P_i=P₂(D_i) with optional labelled slots. In D_i, a chosen directed path runs in slot one and its complement in optional slot two. Let L_i^e,L_i^o be its weighted counts by length parity when slot one is nonempty and starts at a vertex beaten by v; let R_i^e,R_i^o be the corresponding nonempty end counts when that end beats v. Write L_i=L_i^e+L_i^o, R_i=R_i^e+R_i^o, and let `L_i^-=L_i^e−L_i^o`, `R_i^-=R_i^e−R_i^o`. Write d_i=H(D_i+v)−h_i for the increment on the isolated block with its actual tags.

Expand Rédei–Berge at one colour: each pack of disjoint directed odd cycles contributes 2 per nontrivial cycle. The packs missing v contribute ∏h_i. Otherwise the unique cycle containing v is v→old path→v, with even nonzero old path length. Its old vertices occupy a first and last ordinal block, i≤j. Any untouched block contributes h. If i=j the contribution is d_i∏_(k≠i)h_k. If i<j, the first and last nonempty runs contribute L_i and R_j with parity tracked; intermediate runs are optional. In any nonempty block, the weighted counts of an optional first slot having even versus odd length are exactly equal (each odd-cycle pack atom, including a singleton, is assigned to either slot; assigning any fixed atom flips parity without changing weight). Thus an intermediate nonempty block kills the parity character. For i<j the contribution is

    outer ∙ L_i R_j ∏_(i<k<j)P_k,                   if j≥i+2;
    outer ∙ (L_i R_j + L_i^- R_j^-),               if j=i+1,

where outer means ∏_(k<i or k>j)h_k. The factor 2 from the v-cycle cancels the 1/2 from selecting even parity. For adjacent blocks the sign sidecar survives. This formula includes all paths counted by the odd-pack identity, and signed terms are only an exact rewriting of nonnegative counts.

For a transitive left flank with binary v-beats word u, the weighted start count is the binary integer u itself: if a marked vertex starts the run, each subsequent vertex is independently assigned either to it or to the complementary slot. The signed count cancels except for the final vertex, so it vanishes if the final bit is 0. Dually for a transitive right flank its weighted end count is the reverse-complement integer w, and signed count vanishes if the initial bit is 1.

## A real ±2 landing gate at two independent arbitrary binary arguments

Take D=(56,33,35,6,14,24) in outgoing-mask notation. It is strong and has h=45,P=432. The apex v may use either outgoing mask 6 or 9 on D. Both settings have L=129,R=258, but H(D+v)=137 and 139 respectively. The profile is checked by independent subset DP. Every nonnegative even u,w is allowed. For u>0 use its standard binary word as a left TT flank; for w>0 use the reverse complement of its standard binary word as a right TT flank; omit a flank at 0. The required signed counts are zero by the preceding lemma; the left apex contacts its first vertex, and the right flank's last vertex reaches the apex, ensuring strength. If either flank is absent, D is strong and the v masks are genuinely mixed, filling the corresponding role.

Let F(z) be the Hamilton count of a transitive spine plus apex with binary v→vertex word z. Leading zero padding does not change F. Set F(0)=1. If a word has an inversion i<j, w_i=1,w_j=0, the contribution is 2^(j−i−1)+[j=i+1], by choosing an even subset of intermediate spine vertices in the v-cycle. Hence

    F(2z+1)=F(z),      F(2z)=F(z)+z+[z odd].

The right flank has H=F(w) by reversing the tournament and relabeling its spine backwards. The exact paired reads are therefore

    H=κ +45(F(u)+F(w)−2)+258u+129w+432uw,    κ∈{137,139}.

The same physical flanks produce a literal +2 repair at every u,w. The outputs of the first and second contact choices are respectively 2 and 1 modulo 3; thus a prescribed target prime other than 3 fixes κ rather than leaving it a free phase. There is no claim that every target has legal u,w.

An explicit two-height subfamily uses m copies of the transitive word 10 at left and n at right. Then u=2(4^m−1)/3, w=2(4^n−1)/3. The recurrence for F on successive 10 words verifies the closed forms

    H=192·4^(m+n)−86·4^n+30(m+n)+31       and the same +2,

for all m,n≥0. These are not just frozen-order switches.

## Strong fixed-alphabet obstruction on both prime axes

For another application take an arbitrary ordered string of nontrivial strong deletion blocks drawn from the two graph types C3 and the strong four-vertex tournament with H=5, and allow *arbitrary* contacts to a single apex v. A strong output must have v beating some vertex of the first block and losing to some vertex of the last. If there are at least three blocks, the first-to-last term gives H≥∏ P_i over all the strictly intermediate blocks. The per-block contact alphabet has at most four masks up to rotation in C3, and at most sixteen for the four-vertex block; their optional two-slot counts are 12 and 32. For s sufficiently close to one from below,

    4/12^s + 16/32^s < 1

since the sum is 5/6 at s=1. The count of intermediate words of P-product at most X is at most X^s times the sum of their inverse s-weights. The latter sums over all lengths as a convergent geometric series; the endpoints cost only a fixed factor. Thus the grammar realizes O(X^s) strong values ≤X for some s<1, and misses relative density one of primes both as direct values p and after normalizing by seven. This does not treat arbitrary variable strong blocks or arbitrarily many non-apex vertices joining the existing SCCs. A disconnected output in this grammar cannot generate more prime or seven-times-prime targets except at bounded length: any such count has only one nontrivial SCC, and old SCCs cannot split on adding a vertex. With two or more old blocks, all must merge through the apex.

## Random archive provenance

Use `find math/07-reflections -name '*.md' | sort` and sample eight with Python `random.Random(861774).sample(...,8)`. The sequence in this checkout:
1. jc2-degree18-bc-algebraic-finite-place-genus-floor-opus-20260725
2. merged-line-parity-even-odd-s675b
3. the-equi-axis-equicontinuity-is-equidistribution-is-dodge-not-patch-kps-S27
4. lrc-mobile-common-centre-rank-fixed-zero-gap-and-harmonic-support-codex-20260815
5. sos-toeplitz-certificate-breaks-the-riesz-product-stall-macmini-S169
6. the-lrc14-obstruction-is-the-first-cusp-form-scaled-cores-clocks-are-cusps-boxeph-S220 (refuted transfer; not used as a theorem)
7. u21-is-57-and-the-3N-crossover-closes-to-a-tie-at-27-s710
8. the-mediant-attainer-trichotomy-and-the-parity-of-N-macmini-S28

## Prescribed inverse in O(sqrt(N) log N) candidate evaluations

The two-flank family admits an exact bounded decoder, implemented in `inverse.py`. Select κ by N mod3, reject nonpositive range and parity, and try both orientations for which flank argument x is smaller. Since all other summands of the count after κ are nonnegative, x≤sqrt((N−κ)/432), including x=0. For fixed x the other argument w must solve A w+45F(w)=T, where A=129+432x or 258+432x depending on orientation, and T is the exact remainder after moving the smaller argument terms left.

For each w=2^t q+r, 0≤r<2^t, the inversion formula across the word boundary gives

    F(w)=F(q)+RC_t(r) q + [r has leading zero in width t](q mod2)+F(r)−1.

Here RC is reverse complement with the exact suffix width, and the edge indicator is defined to vanish for width zero. Choose the smallest t≥0 with A2^t≥floor(sqrt(N)). For each suffix r the resulting equation is C q+45 F(q)+45g(q mod2)=D with C=A2^t+45RC_t(r), g∈{0,1}. Induction on the binary recurrence gives 1≤F(q)≤q+1 (if q is odd then F(q)≤q). Consequently only the integer interval from max(0,ceil((D−90)/(C+45))) to floor((D−45)/C) needs checking. Its length is O(1) since C is at least floor(sqrt(N)) and D≤N+O(1). Check exact equality, parity and w≥x. The number of suffixes for fixed x is O(1+sqrt(N)/(x+1)), and summing over smaller x yields O(sqrt(N) log N) candidate evaluations, each with polynomial bit cost in log N. Stable +2 physical repair therefore comes with a genuine absolute *family membership* decoder, not a guarantee of acceptance.

The replay recovers all 1,250 input-profile examples and rejects both 65537 and 458759=7·65537 from this family; it accepts 1000000007 via (u,w,κ)=(3533710,0,137). These are exact family-membership statements, not claims about other tournament constructions or the whole spectrum.
