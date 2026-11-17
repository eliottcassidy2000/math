# A nonuniform twelve-vertex directional source and exact growth obstructions

The original one-arc reader observes source profile (H,P₂) with weights (a,b)=(12h−3,15h−4), where 5a−4b=1; its canonical direction is (10,−8). Source discovery below keeps the original switch edge {5,6} and changes its surroundings. All graph masks are zero-based and outgoing.

## Literal linked growth through order twelve

Start A=(174,124,120,240,65,80,1,118) and its mate formed by reversing {5,6}. In turn append vertices of old-vertex outgoing masks 43, 75, 612, 392. The masks refer to all current vertices, and neither second vertex in a two-step block beats its predecessor. The exact lower final rows are

    (3246,3196,2936,1264,3905,2640,2305,1910,1579,2123,2660,392).

The exact mate replaces rows 5,6 by 2576,2337. Their ordinary and optional-two-slot counts are (271787,6325688) and (271797,6325680). Hence the original reader obtains the exact +2 pair

    98146764 h−26118113, 98146764 h−26118111

at every permitted h. At h=1 these are 72028651 and 72028653. At all h the latter is divisible by 3, so this new phase is not a twin-prime family. To test an exact target in either track, solve its affine equation for h by integer division and then apply the original exact ABC height decoder.

This stage was selected by an exhaustive specialized inverse for two appended vertices. Fix the first new contact mask, then inspect a quadratic cut polynomial in the second mask. If C is the old graph with the first added vertex, an additional vertex s contributes

    P_q(C+s) = q P_q(C) + 2q Σ_(i beaten by s, j beating s) W_q(i,j),

where W_q(i,j) counts positive even-length directed paths i→…→j in C, each weighted by P_q of its complement. Indeed every new odd pack either leaves s singleton or puts it in one directed odd cycle; deleting s from that cycle gives the indicated even path. This gives exact q=1 and q=2 responses from subset path DP and complementary path counts. For each first mask use a subset recurrence on the second cut responses; no untested probabilistic assumptions enter.

Of the two established order-ten descendants (12,31) and (43,75), the above is the unique order-twelve continuation preserving the oriented (10,−8), over all 2·1024² pairs of masks. Its full set of 4096² two-vertex extensions has no order-fourteen continuation preserving (10,−8). Hence the exact two-vertices-at-a-time preservation scheme starting at this initial source stops at twelve. This does NOT prohibit a larger batch which passes through a noncanonical intermediate source, a different initial source, or a different arc.

Replays: `python verify.py` checks both literal profiles by independent full subset DP. `extend.cpp` and `extend14.cpp` recheck the exhaustive first- and second-extension scans using the cut polynomial. Orienting the first new vertex toward the second costs nothing in an exhaustive test because relabeling can always choose that orientation.

## Exact unequal-width obstruction for the two ten-vertex templates

Keep the original eight-vertex source, and replace the two successfully added vertices of masks (12,31) or (43,75) by uniform transitive blocks of widths t,u, respectively; they point from the first block to the second. There is no canonical pair in either template except t=u=1, even if the heights are entirely independent.

Here is the finite algebraic certificate. Contract maximal runs in either module to its module label. A projected path consists of eight fixed old vertices (some possibly omitted in an optional slot) and two labels with bounded visits, because in each stretch between fixed vertices the module labels can only go from first to second. The skeleton DP in `tandem.py` counts all walks in this quotient by fixed-vertex mask, numbers of runs of each label, and last label. For q=2 split the fixed-vertex set between two ordered paths, and multiply their skeleton counts at corresponding run totals. Inserting a physical TT_t into its r runs uses exactly Q_r(t)=r!{t choose r}_Stirling assignments; each assigned subset per run has a unique increasing path. Thus the response is exact finite bivariate exponential polynomial in i^t j^u (convert Q by inclusion–exclusion).

For (12,31), the leading corner term of Delta H−10 is +18·4^t 4^u; every other base is bounded by that monomial times (3/4)^min(t,u). Summing absolute coefficients gives dominance when min(t,u)≥16. For (43,75), the maximal terms are 48(4^t 5^u−5^t 4^u). If t≠u their difference is at least a fifth of the larger product, and all remaining terms relative to that product are bounded by coefficient sums times (4/5)^min(t,u), giving exclusion for min≥23. On the diagonal the maximal terms cancel but the remaining leading base is −84·16^t, so an exact univariate dominance bound excludes the large diagonal too.

For each fixed small t or u, the coefficient table gives a univariate exponential polynomial, which is nonzero in at least one of the two coordinates. Cap the remaining height at the first point where the highest nonzero term exceeds all lower terms in absolute value; division by the highest base proves the inequality persists forever. The maximum such strip thresholds are 16 and 23, respectively. `width_certificate.py` independently evaluates the finite strips and returns exactly {(1,1)} in each case. The obstruction is for these two uniform two-module templates, not arbitrary antennae around an arbitrary base.

## A forced third-grade signal in every matching canonical pair

This argument needs no specific template or one-arc assumption. Same-order tournaments with directional difference (10,−8) necessarily have even order by P₂≡(−1)^nH mod3. In the odd-pack expansion, the slot polynomial uses only powers q^d with d≡n mod2, because every atom has odd size. Thus its difference is an even-coefficient polynomial D(z) in z=q², with zero constant. The two readouts impose D(1)=10 and D(4)=−8. In Z[z] the monic factors of the difference give

    D(z) = 14z−4z² + 2z(z−1)(z−4) J(z),   J∈Z[z].

Indeed (D−14z+4z²)/2 vanishes at 0,1,4, so each monic linear factor divides it successively over Z. In particular

    Delta P₃ ≡ −198 (mod 720),      |Delta P₃|≥198.

For any integer q≥3, the representative 14q²−4q⁴ lies strictly inside the half-modulus, giving |Delta P_q|≥4q⁴−14q². The exact two-coordinate match necessarily carries a larger higher-slot signal. A context with one additional visible profile coordinate cannot simply assume that the source is still quiet there. This is a scoped obstruction to transport, not a proof that a specifically compensated higher reader cannot work.
