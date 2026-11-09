> Independent audit PASS: ../bare_coverage_audit/AUDIT.md. The complete bounded union is verified by a different transfer-tree implementation and a byte-for-byte bitset match. This is a fixed-core language, not the unrestricted spectrum.

# Bare-context placement: endpoint defects and exact occupancy collapse

Status: PROVED identities and exhaustive target bound; FINITE-EXACT census in RESULT.json. No unrestricted spectrum or all-target claim. Parent retains commit ownership.

## Inheritance and concept board

Closest mechanism: Session229 dynamic/PROOF.md seven-state composition and pure-T repair. Source: Session231 rerooted-grafts/PROOF.md §4 bare family. Hostile: dropping endpoint sidecars falsely predicts common-context unit spacing. Corrected near miss: the old pure-T family had constant S,E, whereas these smaller bare bodies do not. Least-used sidecar: endpoint S=25 versus45 and E=0 versus14. Live concepts: endpoint defects; strong normalization; last-zero bounds; interval placement versus occupancy; absorption into the singleton language. This artifact is the assigned niche census, not a new unrestricted cutoff.

## Exact legal contexts

Let B_j=X^j⇒C⇒Y^(k−j), with X=101000111000, Y=111011100110, C outgoing rows [4,1,18,7,75,31,175,63], and cut mask124. Write its seven-state (45,q,A,G,T_j,S_j,E_j). The inherited formulas give

q=1440·4096^k, A=878·4096^k+171630J, G=577+2383200J,
J=(4096^k−1)/4095,
Z=(13617152·4096^k+1602510k+717259)/24843,
S_j=25 if j=0 and45 otherwise, E_j=14 if j=k and0 otherwise,
T_j=A+45−E_j−Z+2j.

The C state is (45,1440,878,577,332,25,14), independently verified by literal subset endpoint DP.

Let P be an arbitrary prefix integer, represented without leading zeroes (P=0 means empty prefix). Let R be a length-b suffix integer, including leading zeroes; b=0 means empty suffix. Write rev_b(R) for bit reversal, first_b(R) for its first bit, and last_b(R) for its last bit, all zero when empty. Let f(word) be the sum over selected positions i<j of 2^(j−i−1)+1_(j=i+1). Equivalently f is the seven-state T of the singleton word. Exact apex Hamilton count is

F=(qP+A)(2^b−rev_b(R))+45(R+1−f(R)−last_b(R))
  −PG−45f(P)−T_j−(P mod2)S_j−E_j c_b(R),

where c_b(R)=first_b(R) for b>0 and c_0=1. This follows directly by composing prefix, body and suffix, then using F=A_total−T_total+45−E_total. In particular the empty-suffix formula is

F=Z−2j+(q−G)P−45f(P)−(P mod2)S_j.

Leading unselected prefix vertices are global sources and trailing selected suffix vertices global sinks in the apex tournament. Removing these neutral singleton SCCs preserves F. Thus normalize nonempty prefix to start1 and nonempty suffix to end0. Every normalized graph is strong: the selected first prefix vertex, or selected vertices of the initial strong C, are reachable from apex and reach everything later; the final unselected suffix vertex, or unselected vertices of terminal C, are reachable from everything earlier and reach apex. X starts1 and Y ends0 handle the remaining bare boundary cases.

For a fixed common context the j-step is −2, plus −20(P mod2) at j:0→1, plus14(1−c_b(R)) at j:k−1→k. At k=1 these two defects occur on the same step. Since their possible sum is0,−20,14,−6, full unit spacing for k≥1 holds exactly when P is even and c_b(R)=1. For normalized nonempty suffixes this means they begin1 and end0. For example prefix1 with empty suffix gives 5759665,5759643, not twins. Prefix10 and suffix0 gives24838551,24838563, not twins. Empty context gives2245227,2245225; empty prefix and suffix10 gives9781153,9781151.

## Target-complete bounds, including empty suffix

Fix target cap N. Insertion of a new tournament vertex injects Hamilton paths: choose a deterministic valid insertion position, and deleting the added vertex recovers the input. Hence every context has F≥Z_k−2j≥Z_k−2k. This bounds k.

If normalized suffix is nonempty, its final vertex is unselected. Just before that vertex, write state A_*,T_*. Appending0 gives F=2A_*−T_*+45≥A_*+45, because the covers counted by T_* form a subset of those counted by A_*. Moreover A_*≥2^(b−1)(qP+A). Thus

2^(b−1)A≤N−45,
P≤floor(((N−45)//2^(b−1)−A)/q).

For empty suffix, f(P)≤P because again T_prefix≤A_prefix=P. Therefore

F≥Z−2j+(q−G−45)P−S_j,
P≤floor((N−Z+2j+S_j)/(q−G−45)).

The denominator is positive for every k. These bounds enumerate a finite superset of all normalized contexts below N, without an arbitrary vertex cutoff. Applying the exact F filter makes the census complete for this construction universe.

## All-scale occupancy collapse

There is a more decisive exact reduction. For every k,j and context (P,R), absorb X^j into the prefix and Y^(k−j) into the suffix. The resulting k=0 context is the IDENTICAL cut-labelled tournament, not merely an equal-count representative. Conversely k=0 is part of the family. Therefore the union of all bare interval contexts equals the union of arbitrary singleton contexts around C alone, at every target scale.

This map preserves the entire graph and count. It forgets the designated interval decomposition and hence loses the information that nearby graphs form a consecutive window. The extra sidecar needed to recover that organization is (k,j) and the location of the absorbed X/Y words. The window theorem supplies correlated witnesses, not additional occupied targets. Larger k cannot increase occupancy of the arbitrary-context language. This is stronger than a finite observation and does not rule out new cores, non-singleton surrounding blocks, or an independent theorem that the k=0 language covers some tail.

## Reproduction and finite scope

PREREGISTRATION.md records universe, hypothesis, explicit complete bounds, and controls before the census. A failed proposed suffix0 control was corrected before the coverage run, with exact hostile values retained. Run `python3 -I -B census.py controls`, then `python3 -I -B census.py 100000000`. RESULT.json includes full k≥1 witness records [k,j,P,R,b,order], representation counts, distinct target counts, missing fractions among all positive odds, and the ten largest INTERNAL gaps. A gap tuple (w,a,b) means endpoints a,b attained and the open interval has integer width w=b−a−2, or w/2 missing odd targets. Boundary omissions below the first or above the last attained value are excluded from the internal-gap ranking but included in missing counts.

Positive checks compare the closed formula with seven-state transfer for k≤2 and short contexts, and with independent literal full-graph Hamilton DP for k=0 and short contexts. The finite census includes all k with Z_k−2k≤cap; at cap10^8 these are k=0,1, since Z_2−4=9196067473. There is no claim about arbitrary tournaments or the unrestricted spectrum.

### Census results

At cap10^7: 75,076 attained positive odd targets out of5,000,000, leaving4,924,924 missing (98.49848%). Ten targets have a k≥1 witness; every one already has a k=0 witness.

At cap10^8: 907,984 attained out of50,000,000, leaving49,092,016 missing (98.184032%). There are916,809 legal representations at k=0 and160 at k=1; the latter160 distinct targets all also occur at k=0. The longest internal gap has854 missing odds between attained987,933 and989,643. The second longest at this cap has754 missing odds between77,196,595 and77,198,105. These percentages are finite facts, not asymptotic density claims.

Explicit witness examples, in [k,j,P,R,b,order] format:
- 577: [0,0,0,0,0,9].
- 2,245,227: [1,0,0,0,0,21], also [0,0,0,3814,12,21].
- 2,245,225: [1,1,0,0,0,21], also [0,0,2616,0,0,21].
- 62,532,633: [1,1,0,30,5,26].

The source k=1 window uses12 added vertices to organize two nearby counts that the general k=0 context language already contains. Thus there is no occupancy benefit from increasing the interval parameter inside this precise placement grammar, even though the bare interval theorem remains useful as an interval certificate.
