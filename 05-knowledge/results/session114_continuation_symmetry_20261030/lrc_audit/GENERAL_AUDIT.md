# Independent audit: all-q Farey core

**Verdict: PASS, scoped.** `../lrc/GENERAL_CORE.md` proves the complete maximum-clearance spectrum for C_q={1,...,q−1,q+1} with any positive integer tail, q>=3. It also proves the two-maximizer description and terminal slope on the strict branch2m>s(q), and the stated primorial characterization of tight multiples. The fixed q12 audit remains independently valid.

## Universal proof review

- The determinant-one neighbor basis proves the exact AP Farey tent, not merely an upper envelope. A represented fraction with denominator<=q−1 cannot have both neighbor coefficients positive; the resulting vertical-distance inequality is valid also when one coefficient is zero.
- Neighbor denominator sums>=q. Above1/(q+2), only sumq and sum(q+1) tents matter. The added speedq+1 vanishes at the latter centers, forcing clearance<=1/(q+1+e), where e>=1 is the actual side slope. Thus they disappear from the open top chamber.
- At primitive a/q, the AP slopes are the inverse residue p=a^(-1) modq and q−p. For a other than±1 both slopes>=2; the extra phase moves less than1/q and remains strictly safe. At a=±1 the added wall replaces owner1 withq+1; the other wall cannot bind, including q3. Therefore the minimum remaining owner is the least integer>=2 coprime toq.
- The zero-phase tail argument is applied only after proving the complete windows. It yields the exact strict-branch radius, precisely two maximizing endpoints, and their reciprocal-owner slope.
- On2m<=s(q), the top-chamber inequality excludes higher radii. Primitive(q+2)-grid witnesses attain the floor because(q+2) cannot divideqm under m<=(q−1)/2. The equality boundary s=2,m1 is included without asserting the strict branch's maximizer count.
- The least absent unit s is prime. For m>=2,2m is composite, so2m<=s exactly when every prime<=2m dividesq. This proves the tight-family criterion without asserting a classification of all tight sources.

## Independent finite controls

`independent_general.py` imports no producer functions. It constructs literal safe-interval intersections for q3..80 and compares them with inverse-residue owner windows at four rational radii per core, including the top endpoint. It separately verifies the least-owner identity, floor witnesses, and strict-branch witnesses for m1..s(q)+3. The output is `general_receipt.json`: 1,308 explicit exception-based gates pass identically under normal Python and `python -O`. Reproduce by running `python independent_general.py` and `python -O independent_general.py` in this directory.

These tests audit the formulas. The all-q conclusion uses the proof above; it is not inferred from a finite range. No unrestricted Lonely Runner theorem or literature-priority conclusion follows.

## Completeness clarification and historical gate reconciliation

The final proof explicitly excludes a wrapped extra component in the a=1 AP window. Its original interval is[r,(1−r)/(q−1)]. Below raw phase1, the added-speed distance is already less thanr because(q+1)r>1−r. Raw phase stays below2. The remaining safe branch requires(q+1)t>=1+r, exactly the stated clipped endpoint. Reflection treats a=q−1. This clarification is necessary for a fully explicit containment argument; it changes no formula or test result.

The tightness criterion is inherited from the Goddyn–Wong acceleration criterion, as recorded in `THM-709-doubling-tight-locus-singleton.md`; its elementary derivation here is not a novelty claim. THM-709 has an ID collision, so retain its exact filename when citing.

The historical THM-539 statement that this F-family stops at level4 is already superseded by later repository results. THM-1286 records6/1271 and proposes the primorial cascade. THM-1271 certifies7/16183 and leaves uniform gate completeness open. THM-1276 records further exact instances throughD10. Those examples and the conjectured gate mechanism are inherited. The all-q proof establishes the previously missing uniform statement:

For q=N−1>=3 and D>=3, the canonical family has exact slack-one value D/(D(N+1)−1) if and only if s(q)=2D−1. Equivalently2D−1 is prime, every smaller prime dividesq, and2D−1 does not divideq. This follows by equating the proved strict-branch denominators; the target is strictly above the floor, so the max-floor branch cannot create an exception.

For each primep>=5, let q be the product of primes belowp and m=(p+1)/2. Then s(q)=p=2m−1 and the exact clearance is m/(m(q+2)−1). Its excess above1/(q+2) is exactly1/((q+2)(m(q+2)−1)). With k=q+1 speeds, the normalized excess is less than1/m and tends to zero. Thus explicit attained positive excesses have no uniform positive lower bound of the formc/k². This does not assert the existence of a global second minimum or classify the whole spectrum, and is not a new discovery of the previously certified finite rungs.
