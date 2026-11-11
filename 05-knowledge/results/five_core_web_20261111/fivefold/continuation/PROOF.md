# Pursuing shell transport to actual target selection

**2026-11-10. PROVED IN THIS NOTE:** internal-bit derivative and a centered
doubling inequality, with exact equality boundary. **FINITE-EXACT:** actual
inverse witnesses and joint failure tests on the inherited322 small-parent
holes through `2^24-1`; no new coverage range. **REFUTED:** an axis-only
H9 repair and a constant signed chart-multiplicity identity, in the precise
forms below. **OPEN:** a uniform five-core target selector and Chan–Pak.

This continuation pursues the positive shell identity in
[DYADIC_SHELLS.md](../DYADIC_SHELLS.md), rather than stopping at the fact
that every marginal jump occurs. Two proposed completion mechanisms were
tested: short-flank repairs to the smaller-core holes, and a signed
partition-of-unity formula for the five chart languages. A third probe
keeps Pollock's joint-defect relation on actual compiler witnesses.

## 1. Exact inverse witnesses for the inherited hard targets

The target universe is exactly the322 extra holes of the first three cores
K1,C3,S4, inherited from
`chanpak_five_core_bridges_20261110/cover/support_24.json`.
The largest is16500351. The ten targets requiring one particular H9 core
are retained separately. We did not extend the bound or recalculate the
complete five-core support.

The terminating positive inverse exhausts all H9 witnesses to these
targets, in the inherited44-row quotient by equal signatures and flank
swap. It returns3560 witnesses. Every one is retained in `witnesses.json`
as `[target,row,core_index,A,B]`; the row table fixes literal core masks,
apex cut and the whole response signature. `analysis.json` expands the
ten indispensable-core targets and both worst short-flank targets into
named coordinates. These are actual tournament instructions, not bare
solutions to a desired recurrence.

### Idea1: repair the old holes on an axis or a uniformly short flank

The positive signal is strong but has a precise boundary. Among the322
targets, the minimum of `min(A,B)` over **all** H9 witnesses has histogram

|Minimum short flank|Number of targets|
|---:|---:|
|0|303|
|1|9|
|2|5|
|3|1|
|4|1|
|6|1|
|13|2|

The first axis-only failure is53707. Its complete inverse has just two
witnesses, modulo the declared row quotient:

    S5a, cut5,  A=141, B=4;
    S5b, cut14, A=34,  B=18.

Thus an axis repair fails, while these actual two-flank repairs survive.
The first example requiring short-flank value13 is1759233. Its complete
inverse has just three witnesses:

    S5a, cut16, A=70,   B=297;
    S5b, cut1,  A=357,  B=58;
    S5b, cut12, A=1534, B=13.

Target10826727 also has minimum13, attained by

    S5a, cut24, A=13, B=9464.

Its other seven witnesses are in the receipt. These two targets refute
every proposed cap below13 for this inherited hard-target universe.
They do not prove that a constant cap13 works at all heights, or that
any constant cap is impossible. The stopping gate is now explicit:
**either prove a bounded short-flank theorem, or describe how the needed
short flank grows.** Marginal jump fullness supplies neither.

## 2. A useful positive inequality survives target selection

Write `D_k(n)=2^k-rev_k(n)-topbit_k(n)` as in the shell theorem.

**Proposition1.** For `k>=1,0<=n<2^k`,

    D_k(n)>=g(n)-1,
    g(2^k+n)>=2g(n)-1.                                (1)

Equality holds exactly when the padded k-bit word n starts with1 and
contains at most one0.

Proof: unfold the recurrence for g. Every zero at digit j contributes
`q_j+(q_j mod2)`, where `q_j=floor(n/2^(j+1))`. At a non-top zero,

    q_j+(q_j mod2)<=2^(k-j-1),

the corresponding term of D. Equality requires the entire higher prefix
to consist of ones. A zero in the top position contributes0 to g−1 but2
to D. Thus equality requires a leading1 and no second zero; these
conditions also suffice. QED.

For a fixed actual core/cut and B, put

    alpha=LB+i,
    c_B=h+d(g(B)-1)+oB+e(B mod2),
    Psi_B(A)=Phi(A,B)-c_B-s(A mod2)
             =d(g(A)-1)+alpha A.

Then (1) gives the centered-height inequality

    Psi_B(2^k+n)>=2Psi_B(n).                          (2)

When `alpha>0`, this is a proper height on the fixed-B address axis.
If `alpha=0`, the degenerate axis can have arbitrarily long all-one words
at height zero; the inherited trailing-one normalization remains necessary.
The operation preserves the core, cut, B and apex; it is an actual
shortening rule inside one response chart. It is not an induction on unmarked scalar H:
the centering and available continuation depend on B and the marked chart.

## 3. The joint defect probe distinguishes one step from a whole trajectory

Let D denote the inherited322 small-parent holes, and start with each of
their3560 H9 witnesses. Remove the highest1-bit of one nonzero flank,
provided its bit position is at least1. This is the reverse of the exact
physical high-bit extension, with padding retained.

**FINITE-EXACT:** all5612 such predecessor comparisons land outside D.
Every resulting lower count has a smaller-core realization in the inherited
range. That is a genuine joint statement about the source and increment,
not just their marginal histograms.

The tempting strengthening fails sharply. Allow removal of any internal
1-bit, retaining a strictly smaller H. Among31092 bit deletions there is
one defect-pair receipt:

    core S5b=(2,4,25,3,11), cut12,
    Phi(2,204)=42223,
    Phi(2,196)=40685.

Both totals are in D. The right addresses differ by clearing bit3.
Allowing several high bits to be removed also returns to D:58270 suffix
truncation comparisons give134 incidences. For example

    S5a, cut2: Phi(0,69870)=965301,
               Phi(0,46)=663,

where46 is the six-bit suffix of69870. Thus the observed one-step property
cannot be silently iterated into a defect-free whole trajectory.

### The missing coordinate has an exact derivative

**Proposition2.** For `p>=0,k>=1,0<=s<2^k`, let

    A0=2p*2^k+s,       A1=(2p+1)*2^k+s.

Then

    g(A1)-g(A0)=D_k(s)-p-(p mod2).                    (3)

Proof: before appending the k suffix digits, changing `2p` to `2p+1`
changes g by `-p-(p mod2)`. Each appended zero contributes the appropriate
power of two to the difference of the running prefixes. The first
appended digit contributes one additional unit when it is zero, because
the two initial prefixes have opposite parity. Their total is precisely
`2^k-rev_k(s)-topbit_k(s)=D_k(s)`. QED.

Hence the full compiler's literal single-apex-arc flip has derivative

    Phi(A1,B)-Phi(A0,B)
      =(LB+i)*2^k+d*(D_k(s)-p-(p mod2)).              (4)

The high-bit shell law is p=0. Internal moves have a prefix penalty.
Equation(4), including the symmetric right-flank version, gives a direct
way to investigate target-preserving switches or exact small jumps while
retaining their required source. The42223/40685 pair shows why the prefix
must not be dropped.

**Stopping gate.** Even if the observed5612-case one-step exclusion were
proved uniformly, it starts from an already existing H9 witness. It does
not supply a witness for a hypothetical global hole. Moreover a different
small-parent realization of the same predecessor total need not preserve
its root response. The needed addition is a target-driven entrance into
one of these correlated shell pairs, or an actual response-preserving
rewrite. This is the precise remaining burden of the Pollock transfer.

## 4. Idea2: can signed chart multiplicities prove coverage directly?

There is another way to use a finite atlas: its multiplicity functions
might admit a signed partition of unity. For each of the44 response rows
let `R_y(N)` be the number of normalized address pairs realizing N.
Uniform degenerate axes use the inherited even-address normalization;
all other address pairs are counted literally. These counts are finite
and independent of the chosen larger target bound.

The tested proposal was to choose fixed rational coefficients w_y with

    sum_y w_y R_y(N)=1 for odd N other than7,21,
                    0 for N=7,21.                    (5)

No positivity of the weights was imposed. A second proposal gives a
separate weight to each row and pair of input parities,176 coefficients.

Exact rational linear algebra on odd N through1023 gives:

|Weights allowed|Unknowns|Coefficient rank|Augmented rank|First inconsistent prefix ends at|
|---|---:|---:|---:|---:|
|Row only|44|43|44|47|
|Row and both input parities|176|168|169|309|

Thus both forms of (5) are REFUTED. This is a finite contradiction to a
proposed all-target identity; it is not an asymptotic inference.
`chart_counts.json` records the exact observation rows needed for each
contradiction, and the script uses primitive-integer elimination, computing
rational rank without a potentially misleading reduction modulo a prime.

The strongest survivor is a more informative weighting or a genuinely
nonlinear overlap correction: coefficients may depend on shell, reversed
address, carry, or intersections of chart images. The test does not reject
such formulas. The point of the probe was to see whether the geometry of
five charts could bypass target-by-target source selection; constant
chart/parity weights cannot do so.

## 5. Reproduction

From the repository root:

```bash
c++ -O3 -std=c++17 05-knowledge/results/five_core_web_20261111/fivefold/continuation/inverse.cpp -o 05-knowledge/results/five_core_web_20261111/fivefold/continuation/inverse.bin
./05-knowledge/results/five_core_web_20261111/fivefold/continuation/inverse.bin 05-knowledge/results/chanpak_five_core_bridges_20261110/cover/signatures.txt 05-knowledge/results/five_core_web_20261111/fivefold/continuation/targets.txt 05-knowledge/results/five_core_web_20261111/fivefold/continuation/witnesses.json
python3 -B 05-knowledge/results/five_core_web_20261111/fivefold/continuation/analyse.py
python3 -B 05-knowledge/results/five_core_web_20261111/fivefold/continuation/chart_counts.py
```

The C++ inverse evaluated100103627 target-bounded H9 address pairs and
retained the3560 matching witnesses. Python independently reevaluates
each retained witness through the inherited compiler, checks the exact
scope of the high-bit, internal-bit and suffix probes, and exercises
2046 doubling/equality cases and32736 internal-bit derivatives.
Source hashes are in `analysis.json`. The compiled executable is transient
and is not part of the research artifact. The shell and derivative proofs
are symbolic; no inherited coverage claim is being re-audited here.
