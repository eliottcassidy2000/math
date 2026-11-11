# Carmichael twins, Catalan carry, and the 99·2^k+23 ready ray

Status: PROVED elementary all-k statements below, conditional only on the explicitly inherited physical source formulas where indicated. FINITE-EXACT targeted controls in check.py and receipt.json. OPEN: every-odd-slope unit-offset seeds, and Chan–Pak. This packet does not re-audit historical research or claim a literature novelty search.

## Inheritance and concept board

Anchor: actual odd Hamilton-path count coverage and endpoint readiness. Niche: Carmichael's *totient multiplicity* mechanism transferred to physical source fibers. Wildcard: Catalan's power-gap projection as an explanation of why scalar arithmetic loses the endpoint carry.

The closest proved mechanism is Session200's physical O/J twin, not a gcd identity. Its hostile is two equal F readouts whose next Z outputs disagree. The corrected near miss is “many words/classes imply every count”: word, numerical-profile and tournament-class fibers are different. The least-used sidecar is Session200's joint (m,h) descent set, not just membership of m in the balanced image.

Live board: (1) physical endpoint profile, (2) twin fibers and forbidden chamber, (3) ordered carry versus raw gap, (4) first25-prime phase automaton, (5) the dyadic ready ray with triadic continuation mass. Each new calculation below is typed against these five objects.

## 1. Recovered historical source map

All paths below are repository-relative, recovered from HEAD 13268e7459bb. Some are sparse in the current worktree; `git show HEAD:path` reads their committed contents. This is a source map, not a request to re-verify its literature.

- `session200/carmichael-fibers/PROOF.md`, §§2–6: physical twins, distinct numerical-profile fibers, two-run shifted Mersenne products, and canonical joint descent. This is the direct source of the present arithmetic obstruction.
- `session199/physical-adjoint/PROOF.md`: chronological word reversal with O/Z exchange, preserving the scalar F but requiring transported boundaries for additional observations.
- `01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md` and `05-knowledge/results/root_change_carmichael_chamber_20261107/CARMICHAEL_MAP.md`: fixed root-change-budget singleton chamber and the carefully separated classical Carmichael literature. The main theorem counts words, not profiles/classes. Its O_K(sqrt(X)(log X)^floor(K/2)) singleton bound does not show absence of empty fibers.
- `05-knowledge/results/definability_lenses_spectrum_diophantine_carmichael_20261105.md`, §§0,2,4: strong-class multiplicities; explicit correction of finite missing multiplicities being read as global; physical versus merely integral profile distinction. Hostile formal profile (1,0,6,1) is not a physical source.
- `05-knowledge/results/session218_smith_repunit_totient_20261107/provenance/PROVENANCE.md`: disambiguates totient conjecture, Carmichael numbers, lambda and Camion pseudoprimes. These are not one invariant.
- `05-knowledge/results/catalan_elliptic_20260921_catalan.md`, §§3–5: exact ordered carry B, raw power gap Delta, realization denominator q=|Delta|/gcd(B,|Delta|), and primitive-word versus repeated-word distinction. Catalan classifies a unit-gap subfamily, not the cancellation gate.
- `01-canon/theorems/THM-4516-perfect-power-clocks-are-fermat-catalan-identities.md`: clock-to-power-identity route. This is a restriction of arithmetic clocks, not a target-coverage compiler.
- `05-knowledge/results/berggren_target_return_20261111/reset/TARGET_RETURN.md`, §§1,4,5: actual ready representatives at 99·2^k+23; complete marked response; endpoint-ready failures at 667 and 1997; explicit growing surplus after prescribed digits.
- `05-knowledge/results/berggren_target_return_20261111/arithmetic/PROOF.md`, §§2,5,6: first25-prime wheel has prior context; actual pure-H source pairs and right-length phase retain source identity and baseline. Its 30000 translation is not the present exponential sequence.

Classical literature is inherited only as credited in these files. No use of Carmichael's conjecture or general Catalan's theorem is necessary for the new proofs below.

## 2. Carmichael gives partner construction, not missing-value construction

Use source state s=(D,E,U,d), F(s)=D-E+U and G(s)=D+2U, with chronological operations

    O(s)=(D,0,2U,d),
    Z(s)=(D-E+U,d,2U-d,d),
    J(s)=(D-E+U,D-E,D+U,D).

Starting at (1,0,2,1), these are inherited literal graph births. Session200 proves

    F(O(s))=F(J(s))=G(s),
    F(Z^j O(s))-F(Z^j J(s))=(2^(j+1)-2)(U-d).

The two sources are genuinely distinct because their E values are 0 and D-E>0. This is stronger than a numerical coincidence and weaker than unrestricted source substitution. The invariant needed for subsequent Z is U-d; for arbitrary future words retain the full physical profile and recipe.

Exact source-to-target map: tagged parent recipes (s,O),(s,J) map to two source profiles in one F fiber. Preserved predicate: actual F count. Lost data: endpoint response. Restoration sidecar: U-d, or full state. Cheapest hostile: at the rooted edge, twins have F=5, then one Z yields 11 versus9. Thus a Carmichael pairing increases multiplicity only at already reached targets. It does not fill an empty fiber, which is what odd-count coverage needs.

Session200's deeper descent is already a useful task-specific obligation. For a represented nonbalanced target n, a canonical terminal O Z^j suffix supplies

    n=m+2(2^j-1)h,  m<n/2,
    (m,h)=(D+2U,2U-d) for an actual common source.

The joint source condition is indispensable. For given n,j, the permitted interval is

    n/2^(j+1)<h<n/(2^(j+1)-1),  m=n-2(2^j-1)h.

The historical hostile443 has profiles (219,1,225,1) and (191,1,253,1), canonical (j,h,m)=(5,7,9),(2,63,65), but no j=1 representation in the inherited complete fiber. This suggests a *joint endpoint descent atlas*, indexed by physical (m,h), not a proof by scalar balanced parents.

## 3. New exact obstruction: the entire ready ray avoids the two-run chamber

Let a_k=99·2^k+23, k>=1. The inherited two-run sources O^a Z^b, a,b>=1, have

    F=2+(2^(a+1)-1)(2^(b+1)-1).

**Theorem.** For all k>=1, a_k is not a two-run value above. In particular the known physical ready representatives of a_k cannot be obtained from this particular Carmichael chamber.

**Proof.** Write r=a+1,s=b+1, so r,s>=2. Equality would give

    (2^r-1)(2^s-1)=99·2^k+21.

For k>=3 the right side is5 modulo8. A Mersenne factor with exponent2 is3 modulo8; with exponent>=3 it is7 modulo8. Thus exactly one exponent must be2. The other, say s, would satisfy

    2^s=33·2^k+8.

At k=3 this is272=16·17, not a power of2. For k>=4 its exact2-adic valuation is3, and its quotient by8 is the odd integer33·2^(k-3)+1>1, again impossible. For k=1,2 the products are219=3·73 and417=3·139. Their only possible Mersenne factor>=3 is3 (73+1 and139+1 are not powers of2), and the partner is not Mersenne. QED.

This is an endpoint-grammar exclusion, not a nonrealizability theorem for the count. The incoming R_k construction realizes every one of these targets and is1-ready. Their coexistence shows exactly why a new grammar/endpoint state is needed beyond the simplest twin chamber. Longer O/Z words and root-change words are not excluded.

## 4. What the first25 primes really say about a_k

The primes below100 are exactly25, ending at97. Their count does not impose a relation between99 and23; the precise relation is a finite congruence sieve on the *index k*.

For p odd, p not dividing99·23,

    p divides a_k iff 2^k=-23/99 mod p.

If this target lies in the cyclic subgroup generated by2, there is one forbidden residue k=r_p mod ord_p(2); otherwise none. For p=3,11 the count is23 mod p, nonzero. For p=23 it is99·2^k mod23, nonzero. For k>=1 it is odd, so p=2 never divides it.

The entire nonempty first25-prime phase table is

    p : forbidden k residue / period
    5:3/4; 13:1/12; 17:1/8; 19:9/18;
    29:27/28; 37:23/36; 53:25/52; 59:35/58;
    61:0/60; 67:56/66; 83:15/82.

All other first25 primes have empty forbidden sets. In particular7 never divides the family: 99·2^k+23=2^k+2 mod7 and the powers of2 are1,2,4. The indicator of first25-wheel survival has valid period428468040, the lcm of the displayed periods. We do not claim this is its minimal period. This exact automaton can be used to route factor searches, but its accepted states are not prime certificates or graph endpoint states.

The first wheel survivors are k=2,4,5,6,8,10,12,14. All eight give primes, respectively419,1607,3191,6359,25367,101399,405527,1622039. The next survivor is a hostile:

    k=16: a_k=6488087=1301·4987.

Both factors exceed97. This is the first wheel false positive for1<=k<=16; receipt.json contains complete factorization of every term in that range. It is a plausible explanation of an initially striking prime pattern, not an extension to all k. Trial division by the first25 primes alone certifies a surviving a_k only while a_k<101^2, which includes k<=6. Later primality claims above are individually certified by integer trial division.

The modular observations also prove, for every k>=1,

    a_k=2 mod3,   a_k=5 mod9.

Therefore a_k is never a square or cube, and hence never any perfect power whose exponent is divisible by2 or3. This does not rule out prime exponent>=5. Additionally v_3(a_k+1)=1, so a_k+1 is never any nontrivial perfect power: any such power divisible by3 would have3-adic valuation at least2. These elementary power obstructions do not require Catalan's theorem.

## 5. Catalan's useful import: retain the ordered carry and the observed mode

The recovered Catalan note takes an ordered halving word to (K,L,B), with Delta=2^K-3^L. Its smallest same-clock hostile is

    word(1,3): (Delta,B)=(7,5),
    word(2,2): (Delta,B)=(7,7).

Only the latter has integral fixed point B/Delta at parameter1. Repetition multiplies both B and Delta by the same geometric factor, leaving q=|Delta|/gcd(B,|Delta|) unchanged. This is why a raw power-gap test loses the actual realization predicate.

On the ready ray, the inherited actual R_k state in the marked gauge is

    (H,E,S,U,V,d,K)=
    (99·2^k+23,3,27·2^k+2,
     450·3^k+204·2^k-117,99·2^k+24,27·2^k+3,
     450·3^k+240·2^k-126).

It is physically1-ready since V-H=1, but its full response contains3^k. The analogy is precise: projection to H preserves the target integer and discards a continuation mode, just as (K,L,B)->(K,L) discards the translation. The carrier is different, so the Collatz/Catalan integrality theorem is not imported as a tournament theorem.

**New elementary consequence.** There is no real polynomial P with U(R_k)=P(H(R_k)) for every k>=1. Indeed U/H tends to infinity, whereas U/H^2 tends to0, because the dominant rates are3^k,2^k,4^k. A nonzero polynomial has an integer degree, which would have to be strictly between1 and2. Even an exact scalar dyadic law therefore cannot polynomially recover this sidecar uniformly. On this ray k can of course be recovered by a logarithm, so no information-theoretic impossibility is asserted.

Important nuance: the specific next two births use (H,V,S), not U alone. Thus the nonpolynomial U result obstructs *full response compression*, not by itself every next-digit operation. The actual inherited failure is sharper: after the known ready returns, the k=1 target1997 has no strong1- or2-ready representative; the next repeated two-birth continuation has growing defect. The exact surplus over the ideal ternary digit target is an ordered weighted sum2 sum_i3^(j-1-i)(b_i-1). That is the genuine target-side analogue of the carry: deleting it changes the target.

## 6. Total-count exceptions are not forbidden seed slopes

The parent task distinguishes two statements: Hamilton totals N=7,21 are forbidden; slope parameters c in unit-offset ray seeds N=1+c·2^s need not avoid7 or21. In particular slope c=7 gives N=15 at s=1 and slope c=21 gives N=43. The known all-c<=4097 bank is inherited task context, not recomputed here. HYP9658(A)'s seed predicate H-W=1 with U=c·2^j is a different interface from1-readiness V-H=1.

**Exact ray-coordinate result.** For k>=2,

    v_2(a_k-1)=1,  oddpart(a_k-1)=99·2^(k-1)+11.

For k=1, a_1-1=220=4·55. Consequently the infinitely many ready totals a_k do not provide a single fixed-c unit-offset ray as k varies: their canonical odd slopes vary and their2-adic heights stay1 for k>=2. Nor does V-H=1 alone verify H-W=1. These are distinct coordinates and seed predicates. Proof: a_k-1=2(99·2^(k-1)+11); the bracket is odd when k>=2. No endpoint existence theorem is used in this arithmetic observation.

## 7. Concrete next seed problem, not a gcd surrogate

A useful count-coverage certificate for each odd total N outside7,21 must contain: a literal tournament recipe, H=N, a chosen endpoint/root, the continuation observations used by the digit rule, and an applicable actual return/compensation map. The inheritance C_(n-1)=1 is compatible with every obstruction above because gcd of extension counts says nothing about any particular total N being present.

The new theorem isolates one whole already-represented target ray outside the smallest Carmichael chamber; the phase automaton selects which arithmetic shortcuts are available without erasing the actual source. The historical joint (m,h) descent gives a more promising finite-per-target specification than a scalar prime sieve. The missing universal theorem is that every remaining odd target has either a realizable joint descent parent or another explicit endpoint-ready grammar certificate. Neither the twin multiplicity theorem nor Catalan's unit-gap theorem supplies that existential step.

## Reproduction

    python3 check.py --write
    python3 check.py
    python3 -O check.py

The script uses explicit exceptions, no assertions. Controls cover every first25 prime over four index periods; complete factorizations for k1..16; the two-run obstruction for k1..64 against every possible Mersenne divisor exponent;343 physical twin continuation cases; the443 witnesses; the ordered-carry hostile; and the inherited ready-state formula at k1..12. They are narrow controls for the new mechanisms, not blanket historical verification. The all-k conclusions rest on the proofs in this note, not these finite boxes.
