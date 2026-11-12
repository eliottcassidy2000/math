# Positive-sector granularity and certified real error

Status: **PROVED all-order bounds; FINITE-EXACT inherited-cache comparisons.**
This is an error-correction architecture, not a complete structural classification
of tournaments, and not a proof of Chan–Pak. No outside priority is asserted.

## Inherited mechanism and changed question

Session271 `forecast/PROOF.md`, §§1–3, gives the odd-cycle Burnside formula,
positive support/defect expansion, and relative support-deficit envelope. It
credits THM-4845 and Session98 for the defect hierarchy. Session270 gives the
exact engine and the independently checked cache through 250. We retain those
results rather than auditing or renaming them. The hostile example is that a
tiny relative error multiplied by an exponentially enormous count is not a small
absolute error. Here we quantify a stronger obstruction: **at large order, every
individual positive symmetry sector is already too large to omit for rounding**.
The constructive response is a real interval plus an exact low-residue sidecar,
implemented in the sibling `engine/` lane; the method does not require declaring
an asymptotic forecast exact.

## 1. Granularity theorem: no proper positive-sector truncation rounds exactly

Let λ be any partition of n into odd parts, let r be its number of parts,
m_j their multiplicities, and

    e(λ) = Σ_i (λ_i−1)/2 + Σ_(i<j) gcd(λ_i,λ_j),
    z(λ) = ∏_j j^m_j m_j!,       t(λ)=2^e(λ)/z(λ).

Burnside says a(n)=Σ_λ t(λ). All these are positive rational numbers. Since
gcd≥1,

    e(λ) ≥ (n−r)/2 + r(r−1)/2 = (n+r²−2r)/2.

Also ∏m_j!≤r! and AM–GM gives ∏λ_i≤(n/r)^r. Consequently
z(λ)≤(n/r)^r r!≤n^r, whence completing a square yields

    log₂ t(λ) ≥ n/2+r²/2−r(1+log₂n)
               ≥ [n−(1+log₂n)²]/2.                       (G)

For n≥64, the last expression is increasing and is at least 15/2. For example,
log₂n+1≤n/8 on this range, so the derivative of n−(1+log₂n)² is at least
1−1/(4 log 2)>0. Therefore **every sector t(λ)≥2^(15/2)>181** for n≥64.

It follows that at n≥64, any sum over a proper subset of odd cycle types
underestimates a(n) by more than181. No such sum can recover a(n) merely by
rounding to the nearest integer, however its retained terms are chosen and
however tiny its relative error. This includes support, defect, orbit-loss,
maximum-cycle-length and adaptively ranked positive-sector truncations.

Scope is important. This does NOT obstruct a signed resummation, an exact
modular completion, an interval narrowed by independent information, or a new
identity evaluating the aggregate omitted contribution. It rules out treating
a proper *positive-sector partial sum itself* as an accurate integer predictor.
It also does not claim that each sector is an integer. The bound is on its real
magnitude.

### Fixed-depth corrections remove only linearly many bits

If a fixed support cutoff s leaves an odd cycle of fixed length j>s, its
omitted one-cycle sector has exact contribution

    2^[C(n−j,2)+(n−j)+(j−1)/2] / [j (n−j)!].             (W)

Relative to b(n)=2^C(n,2)/n!, it is

    (n)_j/j · 2^[−(j−1)n+(j²−1)/2].

Thus its absolute logarithm is log₂b(n)−(j−1)n+O(log n). Any fixed number of
support corrections leaves Θ(n²) absolute-error bits. A forecast can gain
linearly many correct leading bits without coming remotely close to correct
low bits. This explicit survivor, not merely a big-O upper bound, explains
why the signed errors in `ERRORS.csv` grow while relative errors shrink.

## 2. A sharper exact-rational support interval

Let μ consist of the nonunit cycles, with support m. For odd-cycle permutations,
every moved unordered pair lies in an orbit of odd length at least3. If there
are r nonunit cycles, each fixed vertex loses m−r≥2m/3 pair orbits. Therefore
its deficit from N=C(n,2) satisfies the inherited inequality

    D ≥ ceil[m(2n−m−1)/3] = d(n,m).

Instead of the earlier crude n^m count of supported permutations, count all
admissible supported permutations exactly. Define

    c_m = [x^m] exp(Σ_(odd j≥3) x^j/j),
    c_0=1,     m c_m = Σ_(odd 3≤j≤m) c_(m−j).

Then m! c_m is the number of fixed-point-free permutations on m vertices
whose cycles are odd. There are n! c_m/(n−m)! permutations of this support.
If L_s(n) is the exact Burnside contribution from all support≤s, positivity
and the deficit inequality give

    L_s(n) ≤ a(n) ≤ L_s(n)+B_s(n),
    B_s(n)=Σ_(m>s, m≥3)^n c_m 2^[N−d(n,m)]/(n−m)!.     (I)

Every quantity in (I) is rational with integer powers of2; there are no
floating-point endpoints or fractional-power rounding problems. Taking
ceil(L_s) and floor(L_s+B_s) gives a certified integer interval. The upper
bound is valid even when loose. Support restrictions c_1=c_2=c_4=0 are handled
automatically. Full support s≥n gives B_s=0 and exact Burnside enumeration.
`certified_interval.py` implements this theorem directly.

This is a tightening and exact-arithmetic implementation of Session271's
support argument, not a newly discovered Burnside formula. It computes the
number of admissible symmetry types/permutations, not classes with a prescribed
Hamiltonian-path count.

## 3. Exact recovery needs a second coordinate

Given any certified integer interval [A,B] containing a(n), and an exact
residue a(n)≡r (mod M), if B−A<M, then at most one integer in [A,B] has this
residue. Existence follows from the certified inputs. Its value is

    A + ((r−A) mod M).

Checking that it is ≤B closes the machine certificate. An interval alone need
not be narrow enough to round; a residue alone need not determine magnitude.
Together they can certify all digits. The sibling exact-engine lane realizes
the residue as a dyadic sector calculation and prunes sectors vanishing modulo
2^b. Its arithmetic identities, pruning proof, experiments and costs belong
there, not to the present real-error theorem.

The conceptual connection is a **two-coordinate completion principle**:
retain coarse magnitude and fine congruence, rather than refining one lossy
observable indefinitely. In Chan–Pak work the analogous missing coordinate is
the physical endpoint/root interface beyond a scalar count. This is a design
analogy, not a transfer theorem from enumerative accuracy to spectrum coverage.

## 4. Reproducible error portrait

`ERRORS.csv` contains every n=0,…,250 and support cutoff0,3,6,9: 1,004 rows.
The signed error is always `a(n)−L_s(n)`, so direction is exact or
underestimate. Both absolute and relative log₁₀ magnitudes are displayed to
32-digit decimal working precision. These logarithmic strings are NOT proof
certificates. Exact numerator/denominator are included through24; beyond that,
the exact rational is specified by the inherited cache and the executable
support formula, with a SHA-256 of its numerator/denominator serialization to
make regenerated values unambiguous. The denominator is recorded at every size.
No exact current values beyond the inherited cache's range250 are asserted.

`CONTROLS.json` checks (I) by exact rational comparisons at14 selected orders,
including14,64,100,150,200,250 and the empty/small boundaries, at all four cuts.
`certified_interval.py` separately computes the full odd-partition Burnside sum
through18, testing each of these cuts. This validates new formulas and boundary
handling, not an expansion of the old numerical evidence for Chan–Pak.

Reproduce from this directory:

    python3 certified_interval.py
    python3 error_table.py

The inherited cache path is Session270 `enumeration/prefix.json`; its hash is
in `CONTROLS.json`. All-digit accuracy for arbitrary giant n still needs an
executed exact recovery, not an enormous-precision prediction relabelled exact.
Writing all digits alone costs Ω(n²) output bits, since the identity lower
bound and a(n)≤2^N imply log₂ a(n)=Θ(n²). Neither these results nor the exact
engine establishes instantaneous enumeration or a finite prime alphabet for
all recursive tournament structures.
