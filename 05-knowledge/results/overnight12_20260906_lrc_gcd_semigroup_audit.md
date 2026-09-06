# Independent referee: exact signed boxes and actual LRC entry

**Status: INDEPENDENT ANALYTIC PASS; FINITE-EXACT independent controls.**
The all-height claims in Sections 2–5 of the
[producer report](overnight12_20260906_lrc_gcd_semigroup.md)
are accepted in their stated scopes. No mathematical repair remains pending.
LRC(14) is not proved. The coefficient lemmas are elementary; the actual-entry
closure additionally inherits THM-3818 and cited eleven-speed LRC.

The audited final producer source is
`15a3fa8f511640c2fee3404cb5cc5b377af4164f8db3dcec7aa1702c2f15e4c3`,
and its output is
`b3f47f0f979c9b126dd6a67b2a83c0dc8dc30fe46b9834d389d3425f82f66ab8`.
The referee did not import or reuse the producer computation. The proof was
derived independently before reading its final prose. A harmless illustration
outside the lemma's hypotheses was replaced by Q=3, (a,b)=(2,3). The producer
also updated the incoming hereditary caps to 90/30 and explicitly excluded
the opposite crossing orientation from its collection of 110 tests.

## 1. Exact central radius and the lost support information

For coprime 1≤a<b≤Q, put L=Q(a+b) and R=L−(a−1)(b−1). Given 0≤x≤R,
the complement M=L−x satisfies M≥ab−a−b+1. Choose u∈[0,b−1] with
au≡M mod b and v=(M−au)/b. A negative v would give M≤ab−a−b;
thus v≥0. Also v≤L/b<2Q and u≤Q. Hence the signed coefficients
Q−u,Q−v are in the required box. At x=R+1, the complement ab−a−b
has no nonnegative representation, so the point is absent. For a=1 this
reduces to the full support endpoint L+1. Symmetry handles negative x.

The strict inequality 2R>L is valid throughout the stated range: using Q≥b,
2R−L≥b(b−a)+2a+2b−2>0. This is the important extra fact beyond a generic
semigroup conductor. The central interval does not describe the whole signed
support: with Q=3 and (a,b)=(2,3), 14 is missing while 15 is present.

The producer's integer interval test is exact for positive, zero, and negative
x. All solutions of ar+bs=x are r=r₀+kb and s=s₀−ka. Intersecting the four
coefficient inequalities gives precisely its displayed ceiling/floor bounds.
The independent verifier checks the constructed witness as well as membership.

## 2. The minimal outside coefficient is sufficient and necessary

Write A=da, B=db, δ=gcd(d,Y), c=d/δ, x=Y/δ. A positive outside coefficient
C in CY=rA+sB must be C=jc. If c>Q there is no eligible C. If x belongs to
the signed box, j=1 works. If x is absent, x>R>L/2, and every j≥2 exceeds
the entire positive support L. Thus a larger external coefficient cannot
rescue a missing minimal point. This proves the claimed iff.

Both restrictions have independent hostile controls. With Q=2, a=1,b=6,
x=3, the point x is absent while 2x=6 is present; b≤Q cannot be omitted.
With physical core 3,6 and outside 1 at Q=2, the normalized point is present
but c=3>Q. These are generic coefficient examples, not claimed decoder entries.

For actual entry, any such crossing has a nonzero outside partial sum. It lies
in W but outside the componentwise decoder span, contradicting W=V_dec. This
uses the actual decoder hypothesis, not a displayed split or a rank heuristic.
All 110 tests concern two core labels and one outside label. The eleven supports
with one core label and both outside labels remain outside this test collection.

## 3. Discrete scale and phase conclusion

From x=gs/δ>R, set z=gcd(δ,s), d=δ/z, e=s/z. Since gcd(d,e)=1 and δ|gs,
g=dh and x=eh. The first scale compatible with these two retained conditions is
d(floor(R/e)+1), including the strict inequality correctly. Its stated optimality
is relative to these two conditions only. It need not realize the defining gcd
δ, avoid every exact box gap, or satisfy actual entry.

The gluing proof preserves the physical scales and gcd(t,g)=1. An eleven-speed
1/12-safe phase gives a closed 1/14-safe arc of total length 1/(42K), where K
is the maximum of the entire primitive core. The g pair lifts become a complete
shifted g-grid in the core clock. Therefore g≥42K suffices, including equality.
The complementary branch 11t≥21q uses the inherited pair arc and the t-grid.

I checked the external supplier in the primary
[Sungkawichai–Trakulthongchai paper, v2, Theorem 1.3](https://arxiv.org/html/2604.23906v2).
Its convention is k nonzero speeds with clearance 1/(k+1), and its theorem covers
k≤12. The use at k=11 is correctly typed. This remains a CITED dependency;
the referee does not reproduce that paper's computer-assisted proof.

The endpoint cutoff follows without the new gcd hierarchy. Failure of the first
branch implies t≤677; D≤H=floor(Q/(42·177)) gives c≤677H<Q. Since
R≥QK/D and p≤177, actual entry forces g>QK/(Dp)≥42K. The integer values
Q=567,869,252,041 and H=76,388,115 agree. No optimality of H as an LRC
failure boundary is asserted. An interior pair can improve R, but it cannot
replace the whole-core maximum in the Lipschitz arc. The (1,…,10,85) hostile
verifies that distinction exactly.

## 4. Actual controls and current lineage

The independent verifier reconstructs all 5,855 atlas pairs and uses union-find
to recover exactly the 11+2 graph components in both producer controls. It checks
the physical Q² box and g>2QK. This last inequality excludes every support-three
crossing in both orientations: a nonzero outside partial sum is an integer
multiple of g, while at most two core terms have absolute sum ≤2QK. A zero
outside partial sum cannot balance a nonzero single core term. Connected internal
decoder rows already span dimension eleven, establishing W=V_dec. This supplies
the necessary global entry check missing from a selected-support test alone.

The referee independently finds tiny-denominator core phases, then selects a
nearest lift of pair phase 1/4. Literal thirteen-coordinate checks verify the
resulting phases, without enumerating an enormous cyclic grid. The unitless
control has best interior radius 29Q−182 at (14,15), versus 22Q−84 at (14,30)
among maximum pairs. It already satisfies the endpoint criterion; no additional
class closed only by interior pairs is claimed.

The incoming `lrc14_joint_shadow_empty_core_next_sep06.md` is the completed,
audited dependency for current subset caps 1,2,4,9,30,90 at sizes12,…,7. Its
finite classifier is inherited, not rerun here; no RESERVED theorem is used.
The primorial example violates these restrictions and is retained solely as
a normalization hostile. The core gcd bound t≤2 applies to a hypothetical
strict failure, not to every safe decoder entry.

## 5. Independent reproduction

[Referee source](../../04-computation/overnight12_20260906_lrc_gcd_semigroup_audit.py)
and [frozen output](overnight12_20260906_lrc_gcd_semigroup_audit.out):

```
python -B 04-computation/overnight12_20260906_lrc_gcd_semigroup_audit.py
python -B -O 04-computation/overnight12_20260906_lrc_gcd_semigroup_audit.py
```

Both runs pass **149,274 always-active gates**, with byte-identical LF outputs.
The independent universes are all 268 primitive pairs at Q=2,…,13 and all 63,962
integers in their support ranges plus two margins; 18,630 literal full outside-
coefficient cases at Q=2,…,7, common scale1,…,9 and Y=1,…,45; 8,325 first-scale
cases R=1,…,37 and δ,s=1,…,15; plus the exact atlas, entry, phase, score, and
hostile controls above. The universal proof is separate from this finite bank.

Referee source SHA256:
`3016fa28cc5c6ca0e9109582014861153308d66f5e5f4966cc6a4776b8dec8bb`.
Output and optimized-output SHA256:
`8b2950db233f2a5806da82ba166e094ac9676aea19b37eb1314e9e695ae5c673`.
No Git, shared navigation, or producer file was modified by this referee.
The input lookup supports adjacent outside-folder artifacts and the eventual
repository layout with sources in `04-computation` and outputs in
`05-knowledge/results`; producer raw-byte hash checks are unchanged.

**Filing:** root integrated these independently audited artifacts in the twelfth
checkpoint. Reproduction commands are relative to the repository root.
