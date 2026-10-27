# Additive supplement: most escaped prime targets really are touchable

2026-10-27. **PROVED corollaries of cited classical results; FINITE-EXACT controls.**
This does not alter the frozen Session73 arithmetic packet. Its A1 is used at
the report hash recorded in `SOURCE_PINS.json`. No all-prime classification,
effective all-height cutoff, or new best Goldbach exponent is claimed.

## S1. Almost every prime target is touchable, unconditionally

Let `P(X)` count primes at most X. Then

`#{prime m<=X: F(m) is nonempty} / P(X) -> 1`.

**Proof.** The classical Goldbach exceptional-set estimate supplies some
`theta<1` such that at most `O(X^theta)` even integers at most X have no
representation as a sum of two odd primes. For a concrete published input we
use `theta=0.914` from the first-page theorem of
[Hongze Li, *The exceptional set of Goldbach numbers (II)*, Acta Arithmetica 92 (2000), 71–88](https://matwbn.icm.edu.pl/ksiazki/aa/aa92/aa9217.pdf).
This older estimate is sufficient; it is not stated as the current record.

For prime `m>=11`, an ordinary partition `m-1=p+q` gives an aliquot source
`pq` unless the two primes coincide. Thus any prime target without a distinct
partition is either associated to an ordinary Goldbach exception, or has the
form `m=2p+1` with p prime. The first set has size `o(P(X))`. The second also
has relative prime density zero: apply the fixed-affine residue argument of
base-report A1b to `m=2p+1` (S=2,T=1). For every auxiliary odd prime ell,
the primitive residue `m=1 mod ell` forces `p=ell`; finitely many such
residue exclusions and fixed-modulus PNT give a product tending to zero.
Every remaining prime has a distinct-prime product source. QED.

This is a derived classical corollary, not a new Goldbach theorem. In particular
we did **not** silently equate ordinary and distinct Goldbach exceptions.
The estimate leaves room for infinitely many exceptional prime targets; it
does not prove that 2 and 5 are the only ones.

## S2. Genuine unbounded-core witnesses, not merely missing templates

Fix a finite core bank A as in the base report. Combining S1 with A1 gives

```
#{prime m<=X: F(m) nonempty and m not in K_A} / P(X) -> 1.
```

Hence for almost every prime target, there is a genuine aliquot witness, but
**every** witness lies outside the finite bank and is nonsquarefull. Every
witness therefore has a prime factor of exponent one; all its admissible core
values lie outside A. If A is `2..C`, all such cores exceed C.
This strengthens the base report's necessary growth coordinate without using
the full distinct Goldbach conjecture. The fixed-bank quantifier remains essential.

## S3. A constructive progression that defeats any fixed core bank

Take a nonempty finite A (the empty-bank case is immediate). For each `a in A`,
put `S_a=s(a)>0`, `T_a=sigma(a)>0`. Choose mutually distinct
odd primes `ell_a` not dividing `S_a T_a`. By CRT choose

`r = T_a mod ell_a` for all a, with `Q=product_a ell_a`.

Then `gcd(r,Q)=1`. If a target `m=r mod Q` belongs to the a-branch
`m=S_a p+T_a` with p prime, divisibility forces `p=ell_a`. Consequently
**every** prime target in this progression beyond
`B=max_a(S_a ell_a+T_a)` evades the entire core bank.

Fixed-modulus PNT, S1, and the base squarefull bound now yield

```
#{prime m<=X: m=r mod Q, F(m) nonempty, m not in K_A}
    ~ P(X)/phi(Q).
```

Thus one explicit positive-relative-density progression supplies infinitely
many **touchable** prime targets which cannot be witnessed by this bank or by
any squarefull source. The residue assignment is a certificate of template
failure, not of untouchability. It is also not an algorithm giving the witness
for every prime in the progression.

Example for `A={2,3,4,5}`: use auxiliary primes `(5,3,11,7)`. Then
`Q=1155`, `r=1063`, `B=40`. The first chosen prime target 1063 has the
explicit source `11*1051=11561`, since `1063=1+11+1051`; both factors exceed
the bank. The checker computes its **complete** fiber by the bound
`n<=(1063-1)^2=1127844` and verifies no squarefull source exists.
An independent guarded-core/squarefull partition reproduces all 37 sources.

For `A=2..10`, the deterministic auxiliary-prime rule gives
`Q=3234846615`, `r=2900072083`. The explicit prime target `9369765313`
has source `13*9369765299`, beyond the bank. This large example certifies the
prime witnesses and bank obstruction only; no complete fiber or absence of
squarefull sources is claimed for it.

## S4. Reproduction and limits

```
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json
python3 check.py --negative-control
python3 -O check.py --negative-control
```

The first two outputs match exactly. The last two must fail: they deliberately
replace `s(3^2)=4` by the false equal-prime value `1+3+3=7`.
The checker also verifies CRT certificates, primes by exact trial division,
and every prime in the `1063 mod 1155` progression through 200000, preserving
the actual source pair. These finite observations are controls, not the proof
of S1–S3. Original-source theorem and hypotheses were checked; the published
analytic proof was not independently re-proved here.
