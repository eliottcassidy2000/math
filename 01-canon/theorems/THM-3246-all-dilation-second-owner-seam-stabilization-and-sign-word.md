---
id: THM-3246
title: "All-dilation second-owner seam stabilization and sign word"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the THM-3224 lane (P,Q;e,f)=(3,5;1,2), every one of the 168 cleared
  cell masses is one explicit quadratic for every integer dilation g>=1.
  The 156 interior cells are certified by exact affine-ray stability; the
  twelve formerly open seam cells have a direct one-alignment formula.
  Consequently the second owner corrector has exactly 156 positive and 12
  negative entries at every dilation, with nonzero Hodge sum 1/24696.  This
  upgrades the old finite scout and its Singer-line counterfeit, but it is
  one second-corrector lane and gives no LRC row exclusion.
source: root/2026-08-03
audit: >
  The assertion-independent companion pins and replays promoted THM-3224;
  checks 9,360 exact interior term-branch certificates and 1,008 direct
  interval controls; verifies the seam component and unaligned-pair bounds;
  reconstructs all 168 quadratics, limits, first corrections, and second
  correctors; and exhausts the 672 projective Singer-gauge classes modulo
  radial F_13^* scaling.  Normal, optimized, and stored transcripts and
  LF-normalized hashes agree.  Independent hostile audits rederived the
  cleared centre-difference bounds, unique aligned seam family and endpoint
  contributions, reflected seam, strict unaligned bounds, ray/head
  quantifiers, complete quadratic and corrector tables, Hodge sum, and
  projective Singer census; a separate implementation matched 5,040 direct
  controls through g=30, both pinned digests, and the THM-3224 complete-orbit
  sum.
depends_on:
  - THM-3224-complete-lrc-orbit-bernoulli-gcd-carry-and-owner-hodge-splitting
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
related:
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
  - THM-3211-uniform-lrc-channel-limit-bernoulli-cubic-and-sharp-floor
script: 04-computation/lrc_second_owner_all_dilation_seam_thm3246.py
output: 05-knowledge/results/lrc_second_owner_all_dilation_seam_thm3246.out
script_sha256: e23b098b38aa2199a348f48f8ab4ac0ce5913c870ead972bd31296494fc25a4b
output_sha256: d7f7dd96b01c597113e78f903cad36246cb47b10e9a1758cb831aa0e83e8cebc
hash_basis: LF-normalized bytes
---

# THM-3246 -- all-dilation second-owner seam stabilization and sign word

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3224 proves the owner-dependent second-corrector/Hodge formalism but
leaves one attractive `168`-owner word as a finite exact scout.  The interior
branch compiler was not the problem: only twelve cells meeting the cut of the
unit interval lacked a stabilization certificate.  Treating that cut as a
seam, rather than asking the interior engine to cross it, closes the word for
every dilation.

## 1. Fixed lane and cleared masses

Fix

```text
(P,Q;e,f)=(3,5;1,2),                L0=168.              (1)
```

For cell owner `j in {0,...,167}` and integer `g>=1`, use THM-3224's exact
cell mass `I_(g,j)` and put

```text
D_g=(504g-1)(840g-2),
N_(g,j)=D_g I_(g,j).                                      (2)
```

The residue period is

```text
M=|Qe-Pf|=1.                                               (3)
```

Thus there is one eventual branch per cell.  The result below proves that it
starts at `g=1`; no hidden pre-stabilization head remains.

## 2. Complete all-dilation quadratic table

For every `g>=1`, the cleared numerator in `(2)` is exactly

| cells `j` | `N_(g,j)` |
|---|---|
| `0<=j<=5` or `162<=j<=167` | `12096g^2-1032g+2` |
| `6<=j<=23` or `144<=j<=161` | `12096g^2-24g` |
| `24<=j<=71` | `(16044-168j)g^2+48g` |
| `72<=j<=95` | `4032g^2+96g` |
| `96<=j<=143` | `(168j-12012)g^2+48g` |

Call this table `(4)`.

In particular

```text
N_(g,167-j)=N_(g,j).                                      (5)
```

### Interior certificate

For `6<=j<=161`, THM-3224's exact floor-moment compiler has period one.
Each cell has

```text
3 residue indices * 10 shift indices * 2 triangle kernels = 60           (6)
```

affine endpoint terms.  At the first stable ray step, the compiler verifies the
affine endpoint slopes and all max/min branch inequalities; those inequalities
then persist for every `g>=2`.  The separately checked `g=1` head agrees with
the same quadratic.  The resulting

```text
156*60=9360                                                   (7)
```

term certificates give the four interior pieces of `(4)`.  This is the same
infinite-ray lemma already used in the proved parts of THM-3224, now exhausted
over the complete interior.

### Direct seam certificate

It remains to prove `(4)` for `0<=j<=5`.  Set

```text
z=504g-1,                 w=840g-2.                       (8)
```

An `e`-tooth with index `k` and an `f`-tooth with index `l` have cleared
centre difference

```text
[w(j+168k)-z(2j+168l)]/168
 =g(168d-j)-(k+d)/3,          d=5k-3l.                   (9)
```

The legal ranges are `0<=k<=3g`, `0<=l<=5g`.  If `d>=1`, the right side of
`(9)` is at least

```text
162g-1/3 > 96g-3/14;
```

if `d<=-1`, its absolute value is at least

```text
168g-1/3 > 96g-3/14.                                    (10)
```

The quantity on the right of `(10)` is exactly the cleared overlap radius
`12(z+w)/168`.  Hence no unaligned pair intersects.  When `d=0`, coprimality
forces

```text
(k,l)=(3t,5t),                0<=t<=g.                   (11)
```

For these pairs the `f`-interval is nested inside the `e`-interval.  After
clearing by `zw`, the three contributions are

```text
t=0:          z(2j+12),
sum_(t=1)^(g-1): (g-1)24z,
t=g:          z max(10-2j,0).                            (12)
```

Their sum is independent of `j`:

```text
z[(2j+12)+24(g-1)+max(10-2j,0)]
 =(504g-1)(24g-2)
 =12096g^2-1032g+2.                                    (13)
```

At `j=5`, the last term is zero and `(13)` is unchanged.  The map `x -> 1-x`
identifies cell `j` with cell `167-j`, proving the opposite seam.  Open/closed
endpoint conventions affect only a null set and not the masses.

## 3. Exact second-owner sign word

Use THM-3224's cell limit `L_j`, first correction `c_j`, and definition

```text
q_j=[kappa_j-2L_j+1848c_j]/423360,                       (14)
```

where `kappa_j` is the constant coefficient in `(4)`.  Since `(3)` has one
residue, `q_j` is independent of `g`.  Substitution of the exact Bernoulli
endpoint values into `(14)` gives

```text
#{j:q_j>0}=156,       #{j:q_j<0}=12,       #{j:q_j=0}=0, (15)

q_j=-17/3087000
  iff j in {0,1,2,3,4,5,162,163,164,165,166,167},        (16)

1/6174000 <= q_j <= 751/666792000       when q_j>0.      (17)
```

Reflection gives `q_(167-j)=q_j`, and the exact Hodge sum is

```text
sum_(j=0)^167 q_j = 1/24696.                             (18)
```

Thus the finite sign scout in THM-3224 was correct, but finite interpolation
was not its proof.  The proof is the interior ray atlas plus the separate
one-alignment seam chart.

## 4. All-dilation Singer-line counterfeit

Let

```text
N={0,1,2,3,4,5,162,163,164,165,166,167}.                (19)
```

These are the negative owners for every `g>=1`.  Their residues modulo `14`
are twelve distinct classes.  Use THM-3234's identification
`F_169~=F_13^2`, fix a primitive `alpha in F_169^*`, and let
`c in F_169^*` and `a in (Z/168Z)^*`.  Under the Singer-equivariant gauge

```text
j |-> c alpha^(aj),              gcd(a,168)=1,           (20)
```

the nonzero points of a vector line through zero are one congruence class
modulo `14`.  Multiplication by the unit `a` and the phase shift `c` preserve
the distinctness of the twelve residues in `(19)`.  Therefore

```text
max_(c in F_169^*, gcd(a,168)=1, vector lines L through 0)
|{c alpha^(aj):j in N} intersect (L\{0})| = 1.         (21)
```

Multiplying `c` by `F_13^*` only rescales a vector line.  Thus the complete
projective census has

```text
phi(168) * 14 = 48 * 14 = 672
```

classes, exactly the number exhausted by the companion.

This promotes THM-3234's finite sign-line hostile to an all-dilation theorem:
even after adjoining the completion point, the negative owner class cannot
be an affine/vector line.  The obstruction is cyclic placement, not its
cardinality `12` (or `13` after adjoining a head).

## 5. Scope

The word `negative` in `(15)--(21)` refers to the second asymptotic corrector
`q_j`, not to the positive cell mass `I_(g,j)`.  This theorem handles one
ordered lane only.  It supplies no safety implication, canonical endpoint
current, Heisenberg clutch, owner-to-root map, row exclusion, or decrement of
the `LRC(14)` ledger.

Its reusable gain is methodological: a circular seam should be given its own
chart.  Forcing a nonwrapping affine-ray compiler across the cut created the
old stabilization gap even though the exact boundary geometry was simpler
than the interior.

## 6. Exact companion

Run

```text
python3 04-computation/lrc_second_owner_all_dilation_seam_thm3246.py
python3 -O 04-computation/lrc_second_owner_all_dilation_seam_thm3246.py
```

and compare LF-normalized bytes with the declared output.  The companion
pins and replays THM-3224; checks the symbolic seam decomposition, every
interior ray certificate, and `1008` direct interval controls; reconstructs
the complete `q` word and pins the row-table and `q`-word digests; and exhausts
all projective unit/phase Singer-gauge classes.  It uses exact rational
arithmetic and no optimization-sensitive assertions, floating point,
randomness, or discovery cache.

QED.
