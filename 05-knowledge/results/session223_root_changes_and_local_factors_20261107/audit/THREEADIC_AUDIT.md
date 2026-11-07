# Independent audit of the fixed-length mod-3 theorem

## Verdict: PASS for the supplied theorem; useful additional dependence theorem

Audited `../threeadic/PROOF.md`, sections 1–8, as supplied on 2026-11-07. All fixed-length exactly-two-J probabilities, the one-J comparison, the O(1/L) error, and the every-precision counting identity are correct. No claim about height-window asymptotics follows. The report already states this limitation correctly.

I independently built all 81 states and their O/Z directed graph directly from the physical formulas, found communicating classes by reachability, checked closed-class uniformity by incoming edge counts, and computed the projection from one-step entry into those classes. This does not import the author's projection formulas or code. The result is six closed classes of sizes 1,1,1,6,18,18 and 36 transient states. Exact rational contraction gives

- zero J: (1/3,1/3,1/3);
- one J: (41/108,8/27,35/108);
- two J: (34/81,47/162,47/162).

An independent integer DP on (state,J-count) counts every length through 200 with at most two J. At selected lengths, exactly-two-J totals equal binom(L,2)2^(L−2). The observed L-scaled errors remain bounded in the displayed sample, consistent with but not a proof of the proved error estimate. Both scripts use explicit checks active under `-O`.

Files: `audit_threeadic.py`, `THREEADIC_CHECKS.json`, identical `THREEADIC_CHECKS_OPT.json`.

## Proof review

1. The eighteen-state irreducibility argument is valid. On E=0, O gives S(D,U)=(D,−U) and ZO gives T(D,U)=(D+U,U+d). Explicitly T^(-1)(D,U)=(D−U+d,U−d), whereas STS(D,U)=(D−U,U−d); this supplies a nonzero D-translation. Together with T it supplies U-translations. A self-loop exists, proving aperiodicity.
2. The six-state d=0,U!=0 class has period exactly two; P^2 has a self-loop and a nonzero D-translation on each U slice. The three zero-U classes are absorbing. Arbitrary E enters a closed class in one step.
3. The coefficient 1/12 on delta_0 and 1/36 on delta_(-a) after pi_a J Pi is correct. Enumerating the two E choices at D=U=0 is important; overlooking that branch would incorrectly give uniform residues.
4. The length-triple weighting is exactly uniform: every triple a+b+c=L−2 has 2^(L−2) words. The leading z=1 pole has order three with coefficient pi_1 R Pi R Pi. All other unit-circle poles have order at most two because the initial d=1 row has no z=−1 pole. Finite stochastic matrices have semisimple peripheral spectrum because their powers are bounded. This justifies the O(L)/Theta(L^2) error without suppressing period two. The one-J argument is the same with two resolvents.
5. The modulo-3^k transfer is an exact word-count identity, not a higher-adic mixing theorem. The distinction is correctly stated.

## Additional theorem: the modular law is correlated with a real-state partition

This directly sharpens the stated real/3-adic gap. It is NOT safe to assume independence from arbitrary real-state sidecars.

Write the sampled word A J B J C, with B nonempty. Let sigma be the last letter of B, and let s=(D,E,U,d) be the actual integer state just before the second J. Physical inequalities give

    sigma=O  iff  D/U <= 1/2,
    sigma=Z  iff  D/U > 1/2.

Indeed an O gives D'/U'=D/(2U)<=1/2 since D<=U. A Z gives

    D'/U'=(D−E+U)/(2U−d)>1/2,

since 2(D−E)+d>0, using D−E>=d>=1. These two alternatives exhaust B's last letter. Thus this is a real ratio partition, not merely an artificial modular condition.

**Claim (PROVED).** Conditional on B nonempty and sigma, under the same fixed-total-length ensemble as the audited theorem, the final F-residue law tends to

    sigma=O: (71/162, 91/324, 91/324),
    sigma=Z: (65/162, 97/324, 97/324),

with error O(1/L). Their average is the audited unconditional limit. In particular the limiting zero-residue probabilities differ by 1/27.

### Exact derivation

Let O and Z now denote their deterministic modular transition matrices. Put

    v=pi_1 R Pi=(pi_1+pi_2)/3+(2/9)eta+(1/12)delta_0+(1/36)delta_2.

The conditional leading row is v O R Pi or v Z R Pi. Under pi_a O, D,U are uniform and E=0. Under pi_a Z, D,U are uniform and E=a. Therefore

    pi_a O R Pi=(pi_1+pi_2)/3+(2/9)eta+(1/9)delta_0,
    pi_a Z R Pi=(pi_1+pi_2)/3+(2/9)eta+(1/18)delta_0+(1/18)delta_(-a).

Also eta O=eta Z=eta and delta_c O=delta_c Z=delta_c. For the O case the absorbing mass is 17/108 at zero and none elsewhere; the remaining 91/108 has uniform F. For Z the absorbing masses are 13/108 at zero and 1/54 at each nonzero residue, again leaving uniform mass 91/108. This yields the displayed laws. The independent graph-based code verifies these contractions exactly.

For the fixed-length justification, write B=B' sigma. The remaining O/Z lengths satisfy |A|+|B'|+|C|=L−3, uniformly. The normalized generating expression is the original three-resolvent expression with the fixed sigma matrix inserted immediately before the second R. The initial resolvent still has no pole at −1; hence the same leading-pole argument proves O(1/L). The event B empty has probability 2/L and does not affect the unconditional limit. Each specified sigma event has probability (L−2)/(2L).

### What it does and does not obstruct

This disproves a blanket assertion that final residue is asymptotically independent of the real history of the rooted state. It does NOT prove dependence on the final scalar height F/2^L alone, and does NOT refute any correctly scoped scalar-height equidistribution statement. A joint renewal kernel can retain this two-piece real partition rather than discarding it. Actual height-window factors still require a theorem for that kernel.

## Other concrete limitations and constructive next step

For exactly two J at length L, the inherited height inequalities imply

    2^(L−2)+2 <= F <= 2^(L+2).

The lower bound is W>=2^(#O+#Z); the upper bound follows from U at most doubling per letter, U(t)=2, and F<=D+U<=2U. Thus F/2^L lies in a fixed compact interval (1/4,4]. A height window only involves a bounded band of nearby length scales, but it still selects a nonuniform part of each real distribution. Compactness alone does not remove this selection.

No fixed modular quotient can determine every height cutoff: t and O^(2*3^(k−1))t coincide modulo 3^k, while their F-values differ. This follows from 2^(2*3^(k−1))=1 mod3^k and the explicit states (1,0,2^(h+1),1). It rules out a residue-only threshold test, not a richer real/modular operator.

The disciplined next target is a joint transfer/renewal theorem on normalized real state and finite modular state, explicitly retaining the pre-J ratio partition above. The fixed-length theorem supplies its exact modular marginal and its conditional control cases. One must then integrate the correct height-window indicator and show the effect of lower-J terms in M_2, rather than substitute fixed-length probabilities directly into HYP-9494(E).
