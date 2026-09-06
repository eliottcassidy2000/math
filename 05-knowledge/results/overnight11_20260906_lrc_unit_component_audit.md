# Independent audit: actual-entry primitive-unit eleven-component closure

**Status: PASS.** The all-parameter proof in
[the primary report](overnight11_20260906_lrc_unit_component.md)
is correct relative to its explicitly inherited actual-entry hypotheses and
the cited eleven-speed LRC supplier. No repair is requested. The numerical
controls below are **FINITE-EXACT**; they are not the basis for the universal
conclusion and do not enlarge its decoder-entry scope.

## 1. Sources and exact dependency boundary

I read the complete primary proof and standalone source, and independently
read **THM-3818**, Sections 6.4--6.5, especially (15o)--(15q) and (15aa)--(15ab):
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`.
I also checked the pigeonhole supplier **THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`.
The physical box, full bounded-relation span, and actual decoder equality
are essential hypotheses; the result does not follow from an abstract
rank-eleven relation subspace or a displayed 11+2 decomposition alone.

For the external safe-arc supplier I read the current CORE-PAPERS entry and
the primary source: Touch Sungkawichai and Tanupat Trakulthongchai,
*Eleven, twelve, and thirteen lonely runners*, arXiv:2604.23906v2,
Theorem 1.3, printed page 2. Its convention is k nonzero speeds, and it states
LRC(k) for k<=12. Thus eleven nonzero speeds have a phase with clearance
1/12. This is a **CITED preprint theorem**; I did not independently audit its
computer-assisted proof. [Primary PDF](https://arxiv.org/pdf/2604.23906v2).

## 2. Independent proof audit

Write the actual primitive physical row as `tV union g(p,q)`, with eleven
distinct positive entries in V, a literal unit in V, `gcd(t,g)=1`, and the
remaining hypotheses stated in the primary report. Put `Q=91^6` and
`K=max V`. Eleven distinct positive coordinates imply K>=11, so the unit
and maximum are distinct physical columns.

The inherited internal-height bound is correctly applied to the primitive
pair (1,K), giving K<=Q. For clarity, its actual-entry mechanism is that
THM-2052 supplies a bounded relation on the two internal coordinates and one
outside coordinate. Decoder equality forbids a crossing relation, so the
outside coefficient vanishes. The remaining internal relation is a nonzero
multiple of the primitive pair relation; its coefficient height forces
K<=Q. Merely knowing gcd(V)=1 would not give this conclusion.

The clamped quotient `a=min(floor(x/K),Q)`, `r=x-aK` correctly represents
every integer `1<=x<=Q(K+1)` using `0<=a,r<=Q`. In the clamped branch,
`r<=Q` follows from the endpoint assumption, not the ordinary Euclidean
remainder bound. At the inclusive endpoint the coefficients are both Q.

If `11t>=21q`, the pair has a closed safe arc of length at least
`11/(21q)>=1/t`. Its elementary one-third-clearance center is obtained by
solving `pj=floor((p+q)/2) mod(p+q)`. A complete shifted t-grid in the pair
clock, obtained from physical lifts of a core-safe phase, meets this closed
arc. The equality boundary is valid for the weak target 1/14.

Otherwise `t<=floor((21q-1)/11)<=677<Q`. With `delta=gcd(t,p)`, the integer
`c=t/delta` is at most Q. If `x=gp/delta<=Q(K+1)`, bounded division gives

    c(gp) - a(tK) - r t = 0.

This is a literal integer relation on at most three distinct physical
coordinates, with all coefficient magnitudes at most Q. Zero a or r is
allowed; a and r cannot both vanish because x>0. Its restriction to the
pair component has the nonzero physical sum cgp. Every vector in the
decoder span has zero physical sum separately on each component. Therefore
the displayed row belongs to W but not V_dec, an actual contradiction to
the stated equality. There is no inference from the rank of selected rows.

The surviving inequality is consequently

    g > delta Q(K+1)/p > 42K,

using delta>=1, p<=177, and `Q>42*177`. Eleven-speed LRC and the
K-Lipschitz clearance function supply a closed 1/14-safe core arc of length
`2(1/12-1/14)/K=1/(42K)`. For a fixed safe pair phase z, the g physical
lifts `(z+j)/g` map in the core clock to a complete shifted g-grid because
gcd(t,g)=1. Its spacing is strictly smaller than that arc length. This proves
the asserted common physical phase. A coprime-scale sidecar cannot be dropped.

The d=3 consumer is a correct weaker consequence of the scale inequality:
with `g=3h`, `K=max(3M,E)`, and p<=177, the bound implies `14h>29K`, hence
both attachment inequalities. The theorem itself needs neither that chart
nor the attachment theorem. The virtual-family discussion correctly tests
entry and does not assert that all displayed virtual rows actually enter.

## 3. Nonvacuity and hostile controls

The primary nonvacuity control

    V=(1,4,6,8,10,12,14,15,16,18,22),
    t=1, g=2^45, (p,q)=(1,3)

has an actual decoder graph with exactly these 11+2 components. My referee
reconstructs it using multiplicatively generated admissible pair sums and a
separate breadth-first component search. It checks the finite physical box,
primitivity, and distinctness. The inequality `g>2Q max V` excludes every
crossing support-at-most-three relation: with one pair term, the two core
terms are too small; with two pair terms, their nonzero sum is an integer
multiple of g and cannot balance one bounded core term. If that pair sum
were zero, a nonzero core term would still be impossible. Connected internal
decoder graphs then have the full dimension eleven, so W=V_dec. A rational
physical safe phase is independently constructed and checked.

The inherited primitive nonunit hostile is also reproduced exactly: it has
gcd one and maximum `237907127334685115>Q`, while every primitive internal
pair height is at most 127. This locates the first failed implication in an
unlicensed extension: gcd one cannot replace the literal normalized unit
when deriving K<=Q. The audit also includes two small independent hostiles:
the signed height-two combination of coordinates 1 and 6 cannot represent
3, and the multiplier-two image of a four-grid misses the closed arc
[1/8,3/8] despite its length being 1/4. These test the coefficient and clock
sidecars separately.

## 4. Independent exact engine and reproduction

[The referee source](../../04-computation/overnight11_20260906_lrc_unit_component_audit.py)
uses only the Python standard library, no producer imports, no floating
arithmetic, and no height census. Its **131,958 always-active gates** cover:

* 4,017 full bounded-division interval points, independently compared with
  Minkowski sums of coefficient boxes, and sharp exterior points;
* 9,085 scale-cleared crossing identities with nonzero component sums;
* all 19,314 coprime pairs p<q with p+q<=356, a universe strictly larger
  than the inherited atlas, checking exact pair centers and cutoff bounds;
* an independent multiplicative reconstruction of the 5,855 atlas entries;
* 47 literal rational physical clock witnesses in both branches, including
  a `11t=21q` control, plus 714 complete closed-grid transports;
* the actual-entry nonvacuity control and the three sidecar hostiles above.

The generic physical clock bank checks only the gluing mechanism. It is not
claimed to certify decoder entry for each member. Its exact finite universe
does not substitute for the all-parameter argument in Section 2.

Run from this directory:

    python -B D-computation/overnight11_20260906_lrc_unit_component_audit.py
    python -B D-computation/overnight11_20260906_lrc_unit_component_audit.py

Normal and optimized transcripts are byte-identical LF. Source SHA256:
`37d2d41531770f047cb3979c214ea0a57b5a8ef8574d1c2b3f4570734b2ba51d`.
Output and optimized-output SHA256:
`22e866526bc150cff6aac2ce41b5610946f2ea799fdaf98f600c64855544dac7`.

Primary source reviewed:
`d2fe7245e8d304557b036b72a03a76a46c06e871f60057e5a137892d1ce27bea`.
Primary output reviewed:
`33ee43062b8f3035d2fa150972d9e817bfa72ec210f71b92db5f4fb282ea2608`.

**Verdict:** the precise primitive-unit actual-entry subclass closes at
clearance at least 1/14. Arbitrary primitive eleven-component shapes and
LRC(14) remain outside this theorem. A separate proposed signed endpoint
extension is not part of this frozen primary or its audit.

**Filing:** root integrated these audited artifacts in the eleventh checkpoint;
reproduction commands above are relative to the repository root. Outside-worktree
locations preserve author provenance, not the present reproduction location.
