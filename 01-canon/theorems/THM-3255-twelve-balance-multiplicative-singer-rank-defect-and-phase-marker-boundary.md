---
id: THM-3255
title: "Twelve-balance multiplicative Singer rank defect and phase-marker boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The positive all-dilation owner-mass word of THM-3253 and the signed Hodge
  word of THM-3252 both have multiplicative C_168 translation rank exactly
  157.  Their common eleven-dimensional defect is the augmentation
  representation of C_12: precisely the nonzero frequencies 14,28,...,154
  vanish.  Every primitive Singer decimation preserves this defect, so all
  8,064 relabellings remain in the same 157-space despite their additive
  13 by 13 nonsingularity. Each of the twelve individual norm-phase slices
  instead has multiplicative rank 168, although its additive support has rank
  at most 7. A centered norm-fibre marker is the exact rank-11 complement;
  every nonzero nonnegative filling marker has rank at least 12, sharply
  attained by one Boolean norm-fibre indicator.
source: root/2026-08-03
audit: >
  The assertion-independent companion pins THM-3246, THM-3252 and THM-3253;
  reconstructs the positive mass coefficient triples and the signed Hodge
  word; verifies twelve equal residue sums; builds all cyclotomic polynomials
  over Z without a CAS; proves exact divisibility and all complementary
  nonvanishing by two integral remainders and eight finite-field quadratic
  certificates; and checks the rank, unit-decimation invariance, phase-atom
  repair and constant hostile. A secondary SymPy companion identifies the
  C_12 quotient with the exact F_169/F_13 norm phase, proves all 24 individual
  Hodge/mass phase-slice unit statements, checks the additive support bound,
  and constructs the sharp signed and nonnegative markers. Independent
  finite-field DFT and field-norm reconstructions reproduced the ranks,
  fibres, gauges and positivity boundary. Normal, optimized and stored
  transcript replay and all declared LF hashes pass.
depends_on:
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
  - THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word
  - THM-3252-singer-compactified-owner-hodge-word-universal-charged-cyclicity
  - THM-3253-positive-owner-mass-newton-cyclicity-and-maximal-common-heisenberg-module
related:
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
  - THM-2543-augmentation-norm-relative-phase-local-system-dichotomy
script: 04-computation/lrc_twelve_balance_singer_rank_defect_thm3255.py
output: 05-knowledge/results/lrc_twelve_balance_singer_rank_defect_thm3255.out
script_sha256: e1b42874879ae1057418bb9aa0f95bf8d5af2140af415e49ea8a0c7f72cfd35f
output_sha256: 7a229b92d63577a3d79eb78b34418a39b42d083ab0a2f065a6bc1106978d6e45
secondary_script: 04-computation/lrc_multiplicative_singer_twelve_balance_discovery_20260803.py
secondary_output: 05-knowledge/results/lrc_multiplicative_singer_twelve_balance_discovery_20260803.out
secondary_script_sha256: 1066ed22d8f6e2ac8a424f62a2d8ea2649133642ebf9c133853893d7b580e143
secondary_output_sha256: cb4a89fe66d7ae38dd11266b66d5a2c6509999cabcef9a7e02e252df71462d4f
hash_basis: LF-normalized bytes
---

# THM-3255 -- twelve-balance multiplicative Singer rank defect and phase-marker boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3252 and THM-3253 prove a striking additive statement: after placing
the 168 owner coefficients on the punctured Singer plane, every primitive
cyclic relabelling is a nonsingular `13 by 13` matrix.  This theorem asks the
orthogonal multiplicative question.  Do phase translation and primitive
decimation of the owner word recover all functions on the 168 cyclic owners?

The answer is no, for a rigid reason.  Both the positive mass word and its
signed second-corrector companion are exactly balanced over the twelve residue
classes modulo 12.  Their common missing quotient is eleven-dimensional.

## 1. The two owner words in `Q[C_168]`

Let

```text
M_g(x)=sum_(j=0)^167 N_(g,j)x^j in Z[x]/(x^168-1),     (1)
```

where `N_(g,j)` is THM-3246's cleared positive mass numerator in the ordered
`(3,5;1,2)` lane.  Thus `M_g` is positive coefficientwise for every integer
`g>=1`.

Let `q_j` be THM-3252's signed Hodge word and clear its common denominator:

```text
scale=32006016000,
Q(x)=sum_(j=0)^167 (scale q_j)x^j in Z[x]/(x^168-1).   (2)
```

The compactification value at the additive-plane origin is deliberately not
part of `(2)`: this section studies the multiplicative phase group on the 168
punctured owners.

## 2. Exact twelve-balance

Writing each mass numerator as a quadratic in `g`, direct summation gives,
for every residue `r mod 12`,

```text
sum_(j == r mod 12) N_(g,j)=120960g^2-528g+2.          (3)
```

The independently reconstructed Hodge word satisfies

```text
sum_(j == r mod 12) scale q_j=108000                  (4)
```

for every `r`.  Consequently both words are divisible by

```text
C_12(x)=1+x+...+x^11
       =(x^12-1)/(x-1)
       =Phi_2 Phi_3 Phi_4 Phi_6 Phi_12.                (5)
```

Equivalently, their Fourier coefficients vanish at exactly the candidate
frequencies

```text
k=14,28,42,56,70,84,98,112,126,140,154.               (6)
```

These are the eleven nontrivial characters of the quotient `C_12`.

## 3. No further zero for the positive word

It remains to show that `(5)` is the full gcd with `x^168-1`, uniformly for
all integer `g>=1`.  The trivial character is

```text
M_g(1)=24(60480g^2-264g+1)>0.                          (7)
```

For primitive roots of orders 8 and 24, exact reduction gives

```text
M_g(x) = -2(504g-1)(x-1)(x+1)^2             mod Phi_8,
M_g(x) =  4(504g-1)(x+1)(x^6-x^2-1)         mod Phi_24. (8)
```

The displayed residual factors have degree below the corresponding
cyclotomic polynomial and are nonzero, while `504g-1!=0` for integral `g`.

For each remaining order, the following table chooses a prime `p`, a
primitive root `z` of that order in `F_p`, and records

```text
M_g(z)=Ag^2+Bg+C,      Delta=B^2-4AC.                   (9)
```

| order | `p` | `z` | `(A,B,C)` | `Delta` |
|---:|---:|---:|---:|---:|
| 7 | 29 | 16 | `(14,0,16)` | 3 |
| 14 | 43 | 27 | `(33,34,18)` | 27 |
| 21 | 211 | 180 | `(173,17,163)` | 167 |
| 28 | 113 | 81 | `(8,91,46)` | 29 |
| 42 | 127 | 27 | `(47,2,76)` | 67 |
| 56 | 113 | 9 | `(98,90,71)` | 43 |
| 84 | 337 | 227 | `(162,270,146)` | 197 |
| 168 | 673 | 625 | `(147,130,250)` | 462 |

In every row `A!=0` and `Delta` is a quadratic nonresidue modulo `p`, so
`(9)` has no solution for any residue `g mod p`.  If an integral specialization
of `M_g` vanished at a complex primitive root of that order, the monic
cyclotomic polynomial would divide `M_g` over `Z`; reducing modulo `p` would
contradict the selected row.  Equations `(7)--(9)` cover every divisor of 168
outside `(5)`.  Hence, for every integer `g>=1`,

```text
gcd(M_g(x),x^168-1)=C_12(x).                           (10)
```

The two exceptional integral remainders in `(8)` are load-bearing.  A
finite-field root-free certificate cannot exist there because `504g-1` has
a root modulo every prime not dividing 504, although it has no integral root.

## 4. The Hodge word has the identical defect

Exact integer division of the reconstructed word `(2)` by all sixteen
cyclotomic factors of `x^168-1` gives zero remainder precisely for

```text
d=2,3,4,6,12.                                          (11)
```

The full remainder-bank digest is

```text
de93c88360c6e17d5cb249ff220eff9ea5a34e499dc54bf9d42d63c5cdc7ed96. (12)
```

Thus

```text
gcd(Q(x),x^168-1)=C_12(x).                             (13)
```

The signed Hodge correction repairs additive `13 by 13` determinants, but it
cannot repair the missing multiplicative characters.

## 5. Exact rank and all 8,064 relabellings

Over characteristic zero, the translation span of a word on `C_168` has one
dimension for each nonzero Fourier coefficient.  Equations `(10)` and `(13)`
therefore give

```text
dim span_C168(M_g)=dim span_C168(Q)=168-deg C_12=157.   (14)
```

A primitive Singer relabelling sends a word `W(x)` to

```text
x^b W(x^a),       a in (Z/168Z)^*, b in Z/168Z.        (15)
```

Translation multiplies Fourier coefficients by units.  Decimation permutes
characters without changing their orders, and hence permutes the complete
set `(6)`.  All `48*168=8064` relabellings of either word therefore lie in,
and individually generate, the same 157-dimensional ideal.  Even allowing
all dilations `g`, all positive mass coefficient words, and the signed Hodge
word cannot leave this ideal because every one is divisible by `C_12`.

This gives the promised contrast:

```text
additive 13x13 matrix rank = 13 in every gauge,
multiplicative C168 phase-orbit rank = 157 < 168.       (16)
```

Additive nonsingularity does not preserve the full multiplicative phase label.

## 6. Norm phase and representation-rank reversal

In THM-3234's deterministic field F_169/F_13, the chosen generator alpha
satisfies

~~~text
Norm(alpha)=6,          ord_(F_13^*)(6)=12.             (17)
~~~

Hence

~~~text
Norm(alpha^j)=6^j,
Norm(alpha^j)=Norm(alpha^k) iff j=k mod 12.             (18)
~~~

The twelve residue classes are exactly the twelve multiplicative norm phases,
each with fourteen points. Equation (3) or (4) says precisely that the
weighted norm pushforward is uniform. Its nontrivial eleven characters are
the vanished modes in (6).

The negative seam owners

~~~text
0,1,2,3,4,5,162,163,164,165,166,167
~~~

form one complete norm-phase transversal. Each Hodge phase therefore contains
thirteen positive owners and one negative owner.

For a residue r, restrict the two words to one fibre:

~~~text
Q_r(x)=sum_(j=r mod 12) scale*q_j*x^j,
M_(g,r)(x)=sum_(j=r mod 12) N_(g,j)*x^j.               (19)
~~~

The secondary companion proves, for every r and every integer g>=1,

~~~text
gcd(Q_r,x^168-1)=gcd(M_(g,r),x^168-1)=1.               (20)
~~~

Thus every individual signed Hodge slice and every fourteen-point positive
mass slice has multiplicative translation rank 168. Cancellation appears
only after the twelve individually cyclic slices are summed, when the rank
drops to 157.

The additive Singer-plane carrier behaves oppositely. The support of one norm
phase meets either eight rows and seven columns or seven rows and eight
columns, alternating with r. Every matrix supported on one phase therefore
has additive rank at most seven, while the sum of all phases has additive
rank thirteen by THM-3253. The same decomposition consequently gives

~~~text
one norm phase: multiplicative rank 168, additive rank at most 7;
sum of phases:   multiplicative rank 157, additive rank 13.
~~~

This is a representation change, not an intertwiner: phase disintegration can
raise one rank while lowering the other.

## 7. Necessary and sufficient sidecar

Let `V in Q[C_168]` be any additional coefficient-plane word and allow all
its phase translations.  Since `M_g` is already nonzero on the other 157
characters,

```text
span_C168(M_g,V)=Q[C_168]
iff Vhat(k)!=0 for every k in (6).                     (21)
```

Let e_r be the Boolean indicator of one norm fibre and put

~~~text
C_r=12*e_r-1.                                          (22)
~~~

The centered marker C_r has Fourier support exactly on the eleven missing
characters. Its translation orbit has rank 11 and is the exact direct
complement to the rank-157 owner ideal. It has zero total and is necessarily
signed: fourteen entries are positive and 154 are negative.

The Boolean marker e_r is nonnegative and has orbit rank 12. It contains the
same eleven missing characters plus the trivial character, so adjoining it
fills the regular module with exactly one trivial-line overlap. This is the
sharp positivity tax. Any nonzero nonnegative rational word has positive
total mass and therefore a nonzero trivial Fourier coefficient; if it also
fills all eleven missing modes, its orbit rank is at least twelve. The fibre
indicator attains that lower bound.

A single labeled phase atom `V=delta_j` has a root of unity, hence a nonzero
coefficient, at every character and repairs all eleven missing modes.  This
is sharp as one additional orbit-generator.  By contrast, a constant word
has only the trivial Fourier coefficient and repairs none of `(6)`.

The theorem therefore identifies the cheapest missing datum: not another
unlabeled mass or asymptotic corrector, but an absolute phase marker (or any
sidecar meeting the eleven tests in `(21)`).

## 8. Scope and holotopy boundary

This is a coefficient-plane theorem.  It neither constructs a canonical
owner-to-Singer identification nor produces a physical LRC endpoint marker.
The one-atom repair is an algebraic sufficiency statement, not an existing
positive observable.  No row/state-to-owner-phase map, oriented edge clutch,
Markov transport, owner semantic alignment, NC2 theorem, Gaussian Moment
Conjecture, or `LRC(14)` decrement follows.

What is ruled out is precise: phase averaging, primitive decimation, changing
the dilation, or adjoining THM-3252's signed Hodge word cannot recover the
lost `C_12` augmentation directions.  Any proposed transport through the
THM-3253 frame must supply a sidecar that breaks twelve-balance.

## 9. Exact companions

Run

```text
python3 04-computation/lrc_twelve_balance_singer_rank_defect_thm3255.py
python3 -O 04-computation/lrc_twelve_balance_singer_rank_defect_thm3255.py
python3 04-computation/lrc_multiplicative_singer_twelve_balance_discovery_20260803.py
python3 -O 04-computation/lrc_multiplicative_singer_twelve_balance_discovery_20260803.py
```

and compare LF-normalized bytes with their declared outputs. The primary
companion uses a CAS-independent exact cyclotomic certificate. The secondary
companion uses exact SymPy polynomial arithmetic for the phase-slice and
marker strengthening. Neither has floating point, randomness, a discovery
cache, or an optimization-sensitive assertion.

QED.
