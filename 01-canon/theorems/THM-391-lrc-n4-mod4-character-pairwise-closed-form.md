---
id: THM-391
name: lrc-n4-mod4-character-pairwise-closed-form
status: PROVED
date: 2026-06-01
session: monad-researcher-2026-06-01-S1
depends_on:
  - HYP-2004
  - THM-369
---

# THM-391: The n=4 pairwise safe measure is a mod-4 odd-character closed form

## Statement

Work at LRC threshold `1/n` with `n=4` (3 runners). For a single positive
integer speed `s`, let the **safe set** be

```text
S_s = { t in [0,1) : ||s t|| >= 1/4 },     |S_s| = 1/2.
```

Let `chi4` be the unique non-principal (odd) Dirichlet character mod 4:
`chi4(1)=+1`, `chi4(3)=-1`, `chi4(even)=0`. Then for distinct positive integers
`a, b` with `d = gcd(a,b)`, `a' = a/d`, `b' = b/d`,

```text
|S_a ∩ S_b|  =  1/4  +  chi4(a') chi4(b') / (4 a' b').
```

Equivalently the pairwise resonance correction is
`R2^(a,b) = chi4(a') chi4(b') / (4 a' b')`.

### Consequences
- The correction is **nonzero iff both cofactors `a', b'` are odd**; it has
  magnitude `1/(4 a' b') <= 1/12` (max at `{a',b'}={1,3}`).
- `|S_a ∩ S_b| >= 1/4 - 1/12 = 1/6`, with equality iff `{a', b'} = {1,3}`
  (e.g. `(a,b)=(1,3)` or any `(d, 3d)`).

This is the exact **mod-4 odd-harmonic analogue** of the n=3 Legendre-symbol
closed form `|B_a ∩ B_b| = 4/9 + (2/9) chi3(a) chi3(b)/(ab)` (HYP-2004), and it
answers handoff item (A) of HYP-2004.

## Proof

The safe indicator `g(x) = 1[||x|| >= 1/4]` is the indicator of the half-circle
`[1/4, 3/4]`. Its Fourier coefficients on `[0,1)` are

```text
g_0 = 1 - 2/4 = 1/2,
g_k = ∫_{1/4}^{3/4} e^{-2πi k x} dx = -sin(πk/2)/(πk)   (k ≠ 0).
```

The arithmetic key is the elementary identity, valid for every integer `k`,

```text
sin(πk/2) = chi4(k):     k even -> 0,   k≡1 (4) -> +1,   k≡3 (4) -> -1.
```

Hence `g_k = -chi4(k)/(πk)`; in particular **all even harmonics vanish**.

By orthogonality (Parseval for the product of two characters of a single
variable `t`),

```text
|S_a ∩ S_b| = ∫_0^1 g(at) g(bt) dt = Σ_{ k a + l b = 0 } g_k g_l .
```

The term `k=l=0` gives `g_0^2 = 1/4`. For nonzero solutions, `gcd(a,b)=d` gives
`k a + l b = 0 ⇔ k a' + l b' = 0 ⇔ (k,l) = (b' M, -a' M)`, `M ∈ Z\{0}`. Thus

```text
correction = Σ_{M≠0} g_{b'M} g_{-a'M}
           = Σ_{M≠0} [ chi4(b'M) chi4(-a'M) ] / [ π^2 (b'M)(-a'M) ].
```

Using `chi4` multiplicative and `chi4(-1) = -1`:

```text
chi4(b'M) chi4(-a'M) = chi4(-1) chi4(a') chi4(b') chi4(M)^2
                     = - chi4(a') chi4(b') chi4(M)^2,
```

and `chi4(M)^2 = 1` for `M` odd, `0` for `M` even. The denominator is
`-π^2 a' b' M^2`. The two minus signs cancel, leaving

```text
correction = ( chi4(a') chi4(b') / (π^2 a' b') ) Σ_{M odd} 1/M^2 .
```

Finally `Σ_{M≠0, odd} 1/M^2 = 2 · π^2/8 = π^2/4`, so

```text
correction = chi4(a') chi4(b') / (4 a' b').    ∎
```

(The `chi4(a')chi4(b')` factor automatically encodes the parity condition: if
either cofactor is even the character is 0 and there is no correction, matching
the vanishing of even harmonics.)

## Verification

`04-computation/lrc_n4_mod4_character_monad.py`: the closed form matches the
exact rational measure `|S_a ∩ S_b|` for **all** pairs `a<b<=14` with **0
mismatches**. The character identity `sin(πk/2)=chi4(k)` is checked for
`|k|<=20`.

## Context / use

This is the pairwise (support-2) term of the full n=4 covering decomposition
`|SAFE| = 1/8 + R2 + R3` (HYP-2004 / HYP-2040). The genuine 3-term resonance
`R3` does **not** collapse to a single character (it is a triple character sum
over a rank-2 resonance lattice) — this is the same wall as larger `n`. See
HYP-2040 for how the n=4 conjecture is nonetheless reduced to a single
measure-gap statement plus one boundary witness.
