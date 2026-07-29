---
id: THM-2917
title: "Diameter-four three-slot resultant atlas and sharpness corollary"
status: >
  PROOF-COMPLETE COROLLARY CANDIDATE + VERIFIED-EXACT; AWAITING
  INDEPENDENT HOSTILE AUDIT.  This is an explicit finite-pattern
  resultant atlas inside the stronger arbitrary-support THM-2824, not a
  new SFC(3) range.  The six translated diameter-at-most-four patterns
  have coefficient-positive one-variable resultants.  Combined with
  THM-2173, the bound is sharp with full support and every positive-depth
  two-charge lift has exact uniform Gaussian detection depth six.
source: root/three-slot-diameter-four-boundary-2026-07-29
depends_on:
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2812-consecutive-three-slot-factorial-moment-six-detection
  - THM-2908-consecutive-four-slot-projective-resultant-closure
script: 04-computation/gmc_all_three_slot_diameter_four_detection_thm2917.py
output: 05-knowledge/results/gmc_all_three_slot_diameter_four_detection_thm2917.out
script_sha256: 4e7880ee2700803ec7205cfcfaeb79198b4343c8ec55a7a26dc1c713b709cd1b
output_sha256: a414aa9760b4e56766510900bf0238a4efa41b08e2d23c254174a63e04e82497
hash_basis: LF-normalized bytes
---

# THM-2917 -- diameter-four three-slot resultant atlas

**PROOF-COMPLETE COROLLARY CANDIDATE + VERIFIED-EXACT; AWAITING
INDEPENDENT HOSTILE AUDIT.**

THM-2824 already proves first-window SFC(3) on **every** arbitrary
three-slot support.  The purpose of this corollary is narrower and
computationally different: it gives a compact, explicit elimination atlas
for all coordinate faces of a five-consecutive-slot window.  Those
resultants are useful as independently replayable boundary certificates
for four-slot Macaulay proofs.  No enlargement of THM-2824's support range
is claimed.

Put

```text
L(s^k)=k!.                                             (1)
```

For every integer `n>=0`, every three-element set

```text
E subset {n,n+1,n+2,n+3,n+4},                         (2)
```

and every nonzero

```text
H=sum_(e in E) z_e s^e,                               (3)
```

at least one of

```text
L(H),                  L(H^2),                  L(H^3) (4)
```

is nonzero.  Equivalently, first-window SFC(3) holds for every
three-slot support of diameter at most four.

The bound is exact on every support: there is a full-support `H` for
which the first two quantities in `(4)` vanish and the third does not.

## 1. Translation-normalized support atlas

Subtract the smallest exponent in `(2)` and absorb it into `n`.  Up to
translation, the complete support list is

```text
(0,1,2), (0,1,3), (0,2,3),
(0,1,4), (0,2,4), (0,3,4).                            (5)
```

There is no reflection quotient here: factorial weights distinguish
`(0,1,4)` from `(0,3,4)`, so both remain in the atlas.

A common zero of `(4)` cannot lie on a coordinate boundary.  Indeed, if

```text
H=A s^u+B s^v,                   0<=u<v,              (6)
```

and `L(H)=0`, then

```text
B/A=-u!/v! in R_<0.                                  (7)
```

Thus `H=A R` for a nonzero real polynomial `R`, and

```text
L(H^2)=A^2 integral_0^infinity R(s)^2e^(-s)ds!=0.    (8)
```

It remains to exclude exact three-term support.

## 2. Uniform one-variable elimination

Fix offsets `(0,p,q)` from `(5)` and write

```text
H=a s^n+b s^(n+p)+c s^(n+q).                          (9)
```

For `m=1,2,3`, normalize by `(mn)!`.  Direct multinomial expansion gives

```text
F_m:=L(H^m)/(mn)!

 =sum_(i+j+k=m) m!/(i!j!k!)
    a^i b^j c^k (mn+1)_(jp+kq),                      (10)
```

where `(x)_r=x(x+1)...(x+r-1)`.  In particular,

```text
F_1=a+(n+1)_p b+(n+1)_q c.                           (11)
```

On exact support set `c=1`, `b=t`, and eliminate `a` through `(11)`.
Let `Q_(p,q)(n,t)` and `C_(p,q)(n,t)` be the resulting second and third
normalized moments.  Their exact profiles, followed by the profile of
their `t`-resultant, are

| offsets | `Q:(deg_t,deg_n,terms)` | `C:(deg_t,deg_n,terms)` | `Res:(deg_n,terms)` |
|---|---:|---:|---:|
| `(0,1,2)` | `(2,4,12)` | `(3,6,22)` | `(18,19)` |
| `(0,1,3)` | `(2,6,15)` | `(3,9,28)` | `(24,25)` |
| `(0,2,3)` | `(2,6,18)` | `(3,9,34)` | `(30,31)` |
| `(0,1,4)` | `(2,8,18)` | `(3,12,34)` | `(30,31)` |
| `(0,2,4)` | `(2,8,21)` | `(3,12,40)` | `(36,37)` |
| `(0,3,4)` | `(2,8,24)` | `(3,12,46)` | `(42,43)` |

Every resultant has the form

```text
u prod_(j=1)^4 (n+j)^e_j P_d(n),                     (12)
```

with `u>0`, `e_j>=0`, and the following nonlinear degree and constant:

| offsets | `u` | nonzero linear powers | `d` | `P_d(0)` |
|---|---:|---|---:|---:|
| `(0,1,2)` | `4` | `(n+1)^7(n+2)^6` | `5` | `29` |
| `(0,1,3)` | `16` | `(n+1)^7(n+2)(n+3)^6` | `10` | `1453192` |
| `(0,2,3)` | `4` | `(n+1)^5(n+2)^7(n+3)^6` | `12` | `39189888` |
| `(0,1,4)` | `4` | `(n+1)^7(n+2)^2(n+3)(n+4)^6` | `14` | `765149683104` |
| `(0,2,4)` | `16` | `(n+1)^5(n+2)^7(n+3)(n+4)^6` | `17` | `16074139599576` |
| `(0,3,4)` | `4` | `(n+1)^5(n+2)^5(n+3)^7(n+4)^6` | `19` | `744201921369600` |

All coefficients of every `P_d` are strictly positive.  The degree
`5`, `10`, and `12` factors agree with the consecutive and
diameter-three controls.  For completeness, the three new
diameter-four factors are

```text
P_14
 =3906250000n^14+153465625000n^13+2141463070625n^12
  +16110905735898n^11+76343348840566n^10
  +245796457278026n^9+559841937053220n^8
  +920753590518838n^7+1101018143554598n^6
  +952814256369302n^5+586833288165003n^4
  +248960417472216n^3+68604780704868n^2
  +10962908594064n+765149683104,                      (13)

P_17
 =186624000000n^17+7959824640000n^16
  +126793629979200n^15+1117460812693568n^14
  +6344299442392480n^13+25067401842485800n^12
  +72047735338528457n^11+154671856554509376n^10
  +251871626451044850n^9+313490106547563094n^8
  +298596370241889104n^7+216611859219731894n^6
  +118304708288976394n^5+47660117398643490n^4
  +13686966390095559n^3+2642193415750482n^2
  +306348272267076n+16074139599576,                  (14)

P_19
 =6034392250000n^19+274851324405000n^18
  +4819671803445425n^17+47598028826544462n^16
  +307673310883599975n^15+1407171549205784070n^14
  +4767066189073755466n^13+12309036752506420908n^12
  +24666885804400963070n^11+38783147945655616356n^10
  +48101045175484840293n^9+47094970755564618678n^8
  +36265629151836952555n^7+21777280851261165678n^6
  +10050382807734748496n^5+3483927192296471952n^4
  +874836500125302000n^3+149796914758435296n^2
  +15604425653224320n+744201921369600.                (15)
```

Therefore `(12)` is strictly positive for every integer `n>=0`.  The
polynomials `Q_(p,q)` and `C_(p,q)` have no common complex root, so no
exact three-term polynomial can make `(4)` vanish.  Together with the
boundary argument, this proves the theorem.

## 3. Full-support sharpness and Gaussian lift

THM-2173 (or the same Krull-height argument) gives, in every prescribed
three-dimensional support envelope, a nonzero `H` satisfying

```text
L(H)=L(H^2)=0.                                        (16)
```

Section 1 forces every such witness to have all three coefficients
nonzero, and the theorem just proved forces

```text
L(H^3)!=0.                                            (17)
```

Thus factorial detection depth three is exact on every support in `(2)`.

If the smallest exponent is positive, write `H=s h(s)`, put `s=ZW` for
a standard complex Gaussian `Z` with `W=conj(Z)`, and define

```text
P=alpha W+Z h(s),                    alpha!=0.         (18)
```

Charge balance gives

```text
E[P^(2j+1)]=0,
E[P^(2j)]=binom(2j,j)alpha^j L(H^j).                 (19)
```

Hence every such genuinely two-charge lift is detected by moment at
most six.  Applying `(19)` to `(16)--(17)` gives a four-monomial,
full-support lift whose moments one through five vanish and whose sixth
is nonzero.  Its exact uniform Gaussian detection depth is therefore

```text
6=(4-1)*2.                                            (20)
```

## 4. Scope and exact verification

This explicit atlas is confined to three slots, diameter at most four, and
the first moment window.  Its detection statement is subsumed by
THM-2824; its new payload is the six-cell coefficient-positive resultant
table and the corresponding explicit sharpness certificate.  It does not
prove a shifted window, SFC(4), an arbitrary-charge effective Gaussian
bound, or a new nullcone classification.

The exact companion reconstructs `(10)` independently for all six
supports, performs every resultant and factorization over `Z[n]`, checks
the complete profiles and coefficient positivity, and hash-pins each
nonlinear coefficient vector.  It has explicit `require` gates and no
floating point.  Run

```text
python 04-computation/gmc_all_three_slot_diameter_four_detection_thm2917.py
python -O 04-computation/gmc_all_three_slot_diameter_four_detection_thm2917.py
```

Both executions byte-match the stored output.  Immutable artifact hashes
and the final status promotion await an independent hostile audit.

**QED (candidate pending independent audit).**
