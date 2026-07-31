# Independent audit of the all-`k` support-transfer ladder

**Audit verdict: PASS, with scope wording required.**  No mathematical or
computational defect was found in the promotion candidate
`04-computation/lrc14_aligned_drift_support_transfer_ladder_thm2928.py` at
source SHA-256
`8db781fb3e7dc8fdc4df2bf3c6d83869a9ffe52f41c7d70c25bbd0a9b0122bea`.
This differs from the audited scratch source SHA
`207dc56b6ddd551d2e9dc75cded60a8fc486b297979b989a8d53a1b864c1cfc6`
only by changing the two user-facing words `scout` to `referee`; all output
after the title line is byte-identical.
The finite filters are necessary conditions, not a cover census, and the
reported denominator objects are unordered denominator multisets, not
physical speed/phase packets.

## 1. Uniform support-transfer lemma

Let `F` be a six-element subset of `{1,...,14}`, put

```text
L=14 lcm(F),
```

and let `J` be the set of full `1/L`-cells safe from the six body combs.
For `D|L`, put `S_D=J mod D`.  Suppose `k` of the seven tail speeds are
aligned, so their normalized common safe set is `R_A`, and let `p=7-k`
tails be genuine drifts.  Write each drift in lowest terms as

```text
z_i/L=c_i/d_i,       gcd(c_i,d_i)=1,       d_i>1,
```

and put

```text
D=lcm(d_1,...,d_p),          a_i=c_i D/d_i.
```

Then `D|L`.  Moreover the quotient speeds `a_i` are distinct: indeed
`a_i=D z_i/L`, so equality of two `a_i` would imply equality of the
corresponding physical tail speeds.

The body-address projection of the aligned safe carrier is

```text
Y_D(A)=union_(r in S_D)(r+R_A)/D.
```

The pieces occupy distinct `1/D`-cells up to null endpoints, hence

```text
mu(Y_D(A))=(|S_D|/D)mu(R_A).
```

If the tails cover, then

```text
Y_D(A) subset union_(i=1)^p D_(a_i).
```

Let `u_j` be any uniform safe-mass floor for `j` distinct radius-`1/14`
comb speeds.  Since the `a_i` are distinct,

```text
mu(R_A)>=u_k,
mu(union_i D_(a_i))<=1-u_p.
```

Therefore every cover necessarily satisfies

```text
|S_D|/D <= (1-u_p)/u_k.                              (ST)
```

The inequality in `(ST)` is non-strict.  No equality row occurs in the
present finite body/divisor census, but that computational fact is not a
license to change the theorem to a strict inequality.

For `k=0` use `R_empty=T` and `u_0=1`.  For `p=0` the right-hand drift union
is empty; this all-aligned case is handled separately and should not be
fed to the positive-arity denominator formula below.

Using THM-2928 `(25a)` and its cited providers,

```text
u_0,...,u_7 =
1, 6/7, 66/91, 55/91, 558/1183, 478/1365, 61/273, 15/154,
```

the exact cutoffs `(1-u_(7-k))/u_k` are

```text
k             0        1        2        3       4       5       6    7
cutoff     139/154  106/117  887/990  125/143  26/31  375/478  39/61  0.
```

These filters must not be described as nested or monotone: the `k=1`
cutoff `106/117` is slightly larger than the `k=0` cutoff `139/154`.

## 2. Exact-lcm denominator-multiset formula

For `p>=1`, let `a_p(D)` be the number of nondecreasing `p`-tuples
(equivalently, multisets of cardinality `p`) of divisors `d_i>1` whose lcm
is exactly `D`.  Let

```text
B_p(e)=C(tau(e)+p-2,p).
```

There are `tau(e)-1` allowed divisors greater than one, so stars and bars
shows that `B_p(e)` counts all such multisets whose entries divide `e`.
Partitioning them by their exact lcm gives the divisor-poset zeta identity

```text
B_p(e)=sum_(d|e) a_p(d).
```

Möbius inversion therefore gives

```text
a_p(D)=sum_(e|D) mu(D/e) C(tau(e)+p-2,p).
```

This permits repeated denominators, as it must: distinct physical or
quotient speeds can have the same reduced denominator.  It correctly gives
`a_p(1)=0`, `a_p(q)=1` for prime `q`, and `a_1(D)=1` for every `D>1`.

The number called `denominator_shapes` is

```text
sum_(surviving D) a_p(D),
```

with each modulus counted once.  The number called `raw_occurrences` is

```text
sum_(surviving body,D rows) a_p(D).
```

Neither number counts numerator choices, physical slopes, phases, or
realized covers.

## 3. Independent computation

The promoted audit script

```text
04-computation/
lrc14_aligned_drift_support_transfer_ladder_independent_audit_thm2928.py
```

has SHA-256
`417830cff7a767227d93bbcee42ad57b75adf2b335dc5fa8fe50e85a972bb792`.
It reconstructs every safe-cell word by a fresh integer endpoint sweep,
projects support by merged cyclic arcs rather than the target bitset, and
counts exact-lcm multisets by downward divisor-poset subtraction rather
than the target Möbius sum.  It checks:

```text
3,003 bodies,
251,536 body/divisor rows,
360 actually occurring moduli,
2,520 recurrence-versus-Möbius formula cases,
420 direct brute-force multiset cases,
220 hostile safe-cell endpoint checks.
```

Its output SHA-256 is
`43a73d9daa2beafb69541db6dc9bf205d9f4d9e4ac0ebf83a268c868551923f4`.
Ordinary and optimized outputs are byte-identical.  The target ordinary and
optimized outputs are also byte-identical, with output SHA-256
`808ec922a881e1d6d9541539ee51ff520e44c4fa7c98208b315c82f91df59e81`.

The independently reproduced ledgers are:

```text
k       0       1       2       3       4       5      6   7
rows 27210   27240   27163   26970   13778   10976   6237   0
bodies
      3003    3003    3003    3003    3003    3003   3003   0
divisors
       219     219     219     217     206     193    171   0
```

```text
k  denominator shapes             raw occurrences
0  161,535,777,082,757             1,504,842,061,942,849
1    3,095,010,121,875                38,954,725,590,760
2       50,874,159,718                   951,545,890,235
3          694,921,995                    21,357,714,101
4            7,483,350                       298,255,882
5               56,419                         3,066,274
6                  171                             6,237
7                    0                                 0
```

All eight semantic digests agree with the target.

Reproduction:

```bash
python3 04-computation/lrc14_aligned_drift_support_transfer_ladder_thm2928.py
python3 -O 04-computation/lrc14_aligned_drift_support_transfer_ladder_thm2928.py
python3 04-computation/lrc14_aligned_drift_support_transfer_ladder_independent_audit_thm2928.py
python3 -O 04-computation/lrc14_aligned_drift_support_transfer_ladder_independent_audit_thm2928.py
```

## 4. Canonical wording

Recommended theorem status:

> **FINITE-EXACT necessary support census.**  For each `k=0,...,7`, every
> literal six-body cover with `k` aligned and `p=7-k` genuine drift tails
> obeys `(ST)`.  Exact enumeration of all `3,003` bodies and all `D|L_F`
> gives the displayed body/divisor ledgers.  For `p>=1`, Möbius inversion
> gives the displayed counts of unordered denominator multisets and their
> body/divisor occurrences.  These are necessary quotient-address counts,
> not counts of numerator/phase packets or realized covers.  The `p=0`
> branch is handled separately and is empty.
