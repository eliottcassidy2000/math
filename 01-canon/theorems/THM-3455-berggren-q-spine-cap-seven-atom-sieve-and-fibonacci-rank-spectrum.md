---
id: THM-3455
title: "Berggren q-spine cap-seven atom sieve and Fibonacci rank spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Intersecting
  THM-3453's fifteen divisor atoms with the parabolic Berggren labels
  q_t=4t(t+1)+3 leaves exactly the atoms 9,11,51 and gives a complete exact
  rank-4/6/7/>7 residue law of minimal period 1683.  Fibonacci sampling has
  minimal rank period 360 and exact labelled-index densities
  1/6,1/4,11/180,47/90.  The associated natural and harmonic coefficients are
  proved, with recurrence index, spine index/rooted depth, and nonlinear-label
  subsets kept distinct.  This remains a literal transverse half-twist theorem and does
  not prove LRC(14).
source: codex-2026-08-15-berggren-q-spine-cap-seven-sieve
audit: >
  independent atom, quadratic-residue, exact-rank, overlap-priority, period,
  Pisano, density, harmonic, indexing, and repo-novelty audits; normal,
  optimized, and stored exact replay; dependency/hash and scope gates
depends_on:
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3415-zero-mode-cochain-global-rank-five-support
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3453-global-literal-half-twist-cap-seven-support-classification
  - THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order
script: 04-computation/berggren_q_spine_cap7_fibonacci_spectrum_thm3455.py
output: 05-knowledge/results/berggren_q_spine_cap7_fibonacci_spectrum_thm3455.out
script_sha256: 1b2c01462bf31844b90deb87f4552ae260a6abeade0f01eb8baad612fbc77306
output_sha256: bb092c7185d7f369ec77154c361c645c44388f027e1b62d657cd1880b7a81344
semantic_sha256: 330499ea2bcbf3d2a0da6d3512870ebfae83c9e8268c78e7291ed00f0d95e652
hash_basis: LF-normalized bytes
---

# THM-3455 -- Berggren q-spine cap-seven atom sieve and Fibonacci rank spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is a repository-new synthesis and corollary of the proved all-modulus
classification.  No literature-priority claim is made.

## 1. Typed setup and inheritance

THM-3334's depth-`d` parabolic Berggren point has scalar

```text
Q_d=4d^2+12d+11.                                        (1)
```

Following THM-3454's parameter indexing, write

```text
q_t=Q_(t-1)=4t(t+1)+3=(2t+1)^2+2,          t>=1.        (2)
```

The corresponding point is `P_t=U^(t-1)(3,4,5)`.  Let `rho_H(M)` be
THM-3453's least number of distinct transverse strict literal half-twist masks
covering every sheet modulo `M`.

THM-3453 proves

```text
rho_H(M)<=7
iff
some d in D divides M,                                  (3)

D={8,9,10,11,12,13,14,15,23,25,29,38,51,68,148},
```

with atom ranks

```text
rank 4: 8,9;
rank 5: 10,12;
rank 6: 11,15,23,25;
rank 7: 13,14,29,38,51,68,148.                          (4)
```

THM-3415/3416 give the complete literal support through rank six, so the
minimum rank of any modulus in `(3)` is the minimum rank of its dividing
atoms.  This minimum convention is essential when several atoms divide one
label.

The closest mechanism is THM-3453's exact divisor descent.  Its canonical
hostile here is

```text
t=20,             q_t=1683=9*11*17.                     (5)
```

The label is divisible by `9`, `11`, and `51`, but its exact rank is four,
not six or seven.  A bare membership list without rank priority is false.

## 2. Only the atoms `9,11,51` meet the spine

Put

```text
x=2t+1,                 q_t=x^2+2.                       (6)
```

First, every `q_t` is odd, so no even atom in `D` divides it.  Modulo five,
`-2` is not a square, so the atoms `15` and `25` are impossible.  For an odd
prime `p`,

```text
(-2|p)=1  iff  p=1 or 3 mod 8.                          (7)
```

The primes `13,23,29` are respectively `5,7,5 mod 8`, so their atoms are
impossible as well.  The only survivors are therefore

```text
9,11,51.                                                 (8)
```

Solving `(6)` in the complete residue rings gives

```text
9 | q_t   iff t=2,6 mod 9;
11 | q_t  iff t=1,9 mod 11;
51 | q_t  iff t=3,20,30,47 mod 51.                      (9)
```

For example, modulo `51` the four roots of `x^2=-2` are
`x=7,10,41,44`, which become the four `t` classes in `(9)` after division by
two.  The exact companion checks every residue for all fifteen atoms.

Combining `(3)--(4)` with `(9)` proves the complete spine law

```text
rho_H(q_t)=4
  iff t=2,6 mod 9;

rho_H(q_t)=6
  iff t is not rank 4 and t=1,9 mod 11;

rho_H(q_t)=7
  iff t is neither rank 4 nor rank 6
      and t=3,20,30,47 mod 51;

rho_H(q_t)>7
  otherwise.                                             (10)
```

In particular, exact rank five never occurs on this odd branch.

The first four geometric labels realize all four outcomes:

```text
(t,q_t,rho_H(q_t))
=(1,11,6),(2,27,4),(3,51,7),(4,83,>7).                  (11)
```

Oddness alone therefore does not force a negative result.

## 3. Exact period and harmonic density on the full spine

The joint period in `(9)` is

```text
lcm(9,11,51)=1683.                                       (12)
```

Within one period, the exact counts are

| exact rank | residue classes | density |
|---:|---:|---:|
| `4` | `374` | `2/9` |
| `6` | `238` | `14/99` |
| `7` | `72` | `8/187` |
| `>7` | `999` | `111/187` |

Here the rank-six count is

```text
(2/11)(1-2/9)=14/99.                                    (13)
```

For rank seven, first remove the rank-four classes.  The four roots modulo
`51` meet the two roots modulo `9` in four classes modulo `153`, because the
required congruences agree modulo three.  Thus the surviving density is

```text
4/51-4/153=8/153.                                        (14)
```

The rank-six predicate modulo `11` is independent of modulus `153`, so
removing it leaves

```text
(8/153)(1-2/11)=8/187.                                   (15)
```

The complete rank word has minimal period `1683`: the exact companion checks
all eleven proper divisors of `(12)`.  The cap-seven support density is

```text
2/9+14/99+8/187=76/187.                                  (16)
```

For each rank class `r`, periodicity gives both natural and harmonic density:

```text
#{t<=T:rho_H(q_t)=r}=delta_r T+O(1),
sum_(t<=T,rho_H(q_t)=r) 1/t=delta_r log T+O(1),          (17)
```

where `delta_r` is the corresponding table entry.  In particular the union
of ranks at most seven has harmonic coefficient `76/187` in the **spine-index**
harmonic series.

Actual rooted depth is `h=t-1`, not `t`.  After omitting the root `h=0`, the
shifted periodic set also satisfies

```text
sum_(1<=h<=H,rho_H(q_(h+1))=r) 1/h=delta_r log H+O(1).  (17a)
```

Thus index and rooted-depth carriers have the same asymptotic coefficient,
but they are not literally the same subset or the same reciprocal series.

The nonlinear labels have different asymptotics.  Since `q_t=4t^2+O(t)`,

```text
#{q_t<=X:rho_H(q_t)<=7}=(38/187)sqrt(X)+O(1),            (18)
```

and

```text
sum_(rho_H(q_t)<=7) 1/q_t                               (19)
```

converges by comparison with `sum 1/(4t^2)`.  Neither the spine-index nor the
shifted rooted-depth harmonic set is interchangeable with the branch-label
harmonic set.

## 4. Fibonacci sampling and the period `360`

Now put

```text
r_n=rho_H(q_(F_n)),                 n>=1.                (20)
```

The relevant Pisano periods are

```text
pi(9)=24,              pi(11)=10,              pi(51)=72.
```

Consequently the labelled rank word has period dividing

```text
lcm(24,10,72)=360.                                      (21)
```

The exact return calculation modulo `1683` gives

```text
(F_360,F_361)=(0,1) mod 1683,                            (22)
```

and no proper divisor of `360` gives this return.  The rank word itself also
has minimal period `360`.

The residue classification can be stated without listing all `360` states:

```text
r_n=4 iff F_n mod 9 is in {2,6};
r_n=6 iff not rank 4 and F_n mod 11 is in {1,9};
r_n=7 iff neither lower rank and
             F_n mod 51 is in {3,20,30,47}.             (23)
```

The individual indicator periods and residue representatives are:

```text
rank 4, period 24:
  n=3,16,20,21 mod 24;

rank 6, period 120:
  n=1,2,9,11,12,19,22,29,31,32,39,41,42,49,52,
    59,61,62,71,72,79,81,82,89,91,101,102,109,111,119
    mod 120;

rank 7, period 360:
  n=4,14,28,43,76,86,100,115,134,148,158,
    173,187,206,220,230,244,245,278,316,317,350
    mod 360.                                              (24)
```

Thus one period splits exactly as

| exact rank | indices | density/harmonic coefficient |
|---:|---:|---:|
| `4` | `60` | `1/6` |
| `6` | `90` | `1/4` |
| `7` | `22` | `11/180` |
| `>7` | `188` | `47/90` |

The cap-seven union has labelled-index density and harmonic coefficient

```text
1/6+1/4+11/180=43/90.                                   (25)
```

This is an index-subset statement: `n` remains a label.  The values
`F_1=F_2=1` and hence `q_(F_1)=q_(F_2)=11` repeat.  If the carrier is instead
an unlabelled subset of Fibonacci values or branch labels, start at `n=2` as
required by MISTAKE-209.  The finite duplicate does not change `(25)`, but it
does change exact set-versus-multiset claims.

The harmonic behaviors bifurcate:

```text
sum_(n:r_n<=7) 1/n              diverges with coefficient 43/90;
sum_(n>=2:r_n<=7) 1/F_n         converges;
sum_(n>=3:r_n<=7) 1/(F_n-1)     converges;
sum_(n>=2:r_n<=7) 1/q_(F_n)     converges.               (26)
```

The latter three follow already from `F_(n+2)>=2F_n`.  The rooted-depth sum
starts at `n=3` because `F_1-1=F_2-1=0`.  A subset may therefore
have positive density in its recurrence **indices** while being exponentially
sparse in both its selected spine-index/rooted-depth values and nonlinear
labels.

## 5. Exact companion and controls

The exact companion uses only integral and rational arithmetic.  It checks:

- every residue modulo every one of THM-3453's fifteen atoms;
- every `t` in the complete period `1683`, with exact rank counts and all
  proper-period hostiles;
- the four consecutive positive controls `(11,27,51,83)` and the triple-atom
  priority hostile `q_20=1683`;
- the exact Pisano return modulo `1683` and all proper divisors of `360`;
- all `360` Fibonacci rank states, the individual periods `24,120,360`, and
  every residue in `(24)`; and
- all rational densities and support sums in `(16)` and `(25)`.

Normal and optimized runs reproduce the stored transcript byte for byte after
LF normalization:

```text
python3 -B 04-computation/berggren_q_spine_cap7_fibonacci_spectrum_thm3455.py
python3 -B -O 04-computation/berggren_q_spine_cap7_fibonacci_spectrum_thm3455.py
```

The script freezes LF-normalized dependency hashes and uses explicit
exceptions rather than optimization-sensitive assertions.

## 6. Scope and non-consequences

| field | boundary |
|---|---|
| source | parabolic Berggren `q_t` labels and their Fibonacci spine-index sample (at rooted depths `F_n-1`) |
| predicate | THM-3453's transverse strict literal half-twist cover rank |
| preserved | atom divisibility, rank priority, periodic index, harmonic coefficient |
| retained existence | every finite-rank grade has a labelled witness at `c=1/(2q_t)` with zero complete mode cochain |
| lost | the compressed word forgets owner masks and widths; other centres, nonzero current, phase transport, and the LRC exit remain absent |
| positive controls | `q_1=11`, `q_2=27`, `q_3=51` |
| hostile | `q_20=1683` forces lower-rank priority |

This theorem classifies only literal transverse half-twist masks on the named
spine labels.  Every finite-rank positive does realize the fixed common time
`c=1/(2q_t)` and a zero complete mode cochain after one labelled witness is
reattached; the periodic rank word alone forgets that witness.  It does not
classify other centres or construct a nonzero LRC current, decrement
certificate, spectral closure, or Jacobian map.  Neither the density `76/187`
nor `43/90` is an LRC success probability.  LRC(14) remains open.
