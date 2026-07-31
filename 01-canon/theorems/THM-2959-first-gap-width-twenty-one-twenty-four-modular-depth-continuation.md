---
id: THM-2959
title: "First-gap width twenty-one through twenty-four modular depth continuation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For each
  first-gap support (0,1,2,M),
  21<=M<=24, THM-2949's fixed rank-35 cofactor is nonzero at
  every nonnegative integer depth.  Width 21 has one explicit
  surplus negative-root factor n+17 beyond THM-2957's floor-law
  factor; after removing it, all four primitive cores have frozen
  finite-field root-free gates.  This is not a complete-width,
  arbitrary-width, or positive-real-ray theorem.
source: codex-gmc-first-gap-modular-continuation-2026-07-29
audit: >
  An independent hostile audit checked FLINT coefficient order,
  integral primitive normalization, all four degree invoices, the
  width-21 surplus, primality and every residue in the four modular
  gates, all twelve outside-grid determinants, the THM-2947 rank-gap
  handoff, and the integer-only scope.  Fresh normal and optimized
  replays reproduced the stored transcript, hashes, and record digest
  exactly.
depends_on:
  - THM-2957-first-gap-width-fifteen-twenty-modular-depth-ladder
  - THM-2949-fixed-rank-thirty-five-cofactor-newton-atlas
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
related:
  - THM-2955-width-twenty-fixed-fifth-compound-mod-97-gate
  - THM-2956-koszul-gale-fixed-fifth-compound-exchange
script: 04-computation/gmc_first_gap_width_twenty_one_twenty_four_modular_depth_continuation_thm2959.py
output: 05-knowledge/results/gmc_first_gap_width_twenty_one_twenty_four_modular_depth_continuation_thm2959.out
script_sha256: 1311076e49d164581bc9ef32ac7cba890c9ec8d2dc72581254e2d5cd257f8605
output_sha256: 82634f764eca493f4354be254c482ddcb760a0404db4075b7bb60b37828c6b83
thm2949_dependency_sha256: 9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54
hash_basis: LF-normalized bytes
---

# THM-2959 -- first-gap widths twenty-one through twenty-four modular continuation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

For \(21\leq M\leq24\), let \(P_M(n)\) be the fixed
\(20Q+10C+5F\), rank-thirty-five cofactor constructed in THM-2949 for
the normalized support

```text
(0,1,2,M).
```

Then

```text
P_M(n) != 0 for every integer n>=0 and every 21<=M<=24.       (1)
```

Consequently the PROVED conjugate-pair rank gap of THM-2947 gives full
rank of the degree-seven Macaulay map at every physical depth for each
of these four supports.  In particular, first-window SFC(4) holds for
every translate of

```text
(0,1,2,M),  21<=M<=24.
```

The certificate is discrete.  Each sampled cofactor has two real-sign
runs before its final positive integer tail.  The proof of nonvanishing
is a finite-field gate on a primitive high-degree core, not positivity
on the full real ray.

## 1. The inherited floor-law factor and the width-21 seam

Write

```text
a=floor(M/3), b=floor(M/2), c=floor(2M/3).
```

Define the THM-2957 floor-law factor

```text
B_M(n)
 =(2n+1)^5 (n+1)^6 (n+2)^21
  prod_(3<=r<=a)     (n+r)^26
  prod_(a<r<=b)      (n+r)^24
  prod_(b<r<=c)      (n+r)^20
  prod_(c<r<=M)      (n+r)^19.                         (2)
```

Let

```text
q200=[x0^2]Q,  c300=[x0^3]C.
```

For widths \(22,23,24\), exact polynomial division gives

```text
P_M = c300*q200^5*B_M*core_M.                          (3)
```

Width \(21\) is the first seam where `(2)` undercounts one harmless
negative-root factor:

```text
P_21 = c300*q200^5*B_21*(n+17)*core_21.                (4)
```

The companion verifies the factor `n+17` by exact division and by the
exact root `-17` of the quotient before division.  This explains why
the uncorrected quotient can have no root-free prime witness: the
residue `-17 mod p` is a root for every prime \(p\).  The root is
negative and does not threaten a physical depth.

For all four widths the companion also verifies

```text
deg q200=M-1,             every coefficient of q200 is positive,
deg c300=2M-1,            every coefficient of -c300 is positive.
```

Thus `q200`, `c300`, every factor in `B_M`, and the width-21 surplus
are nonzero for every real \(n\geq0\).

## 2. Primitive cores and four exact finite-field gates

Primitive-normalize each remaining `core_M` with positive leading
coefficient and call it \(R_M\).  Before integer conversion, the
companion explicitly verifies that every quotient coefficient is
integral.  The exact bank is:

| \(M\) | \(\deg P_M\) | \(\deg B_M\) | surplus | \(\deg R_M\) | root-free prime | SHA-256 of primitive coefficient vector |
|---:|---:|---:|:---|---:|---:|:---|
| 21 | 1120 | 447 | `(n+17)` | 531 | 127 | `3826e18fe51b3f8f3f7359257355823fc107cdb414b54447728ec813dac0e1fe` |
| 22 | 1175 | 470 | none | 557 | 107 | `5fa6bb6af27d164d4d69e40f3da1be5a7f7885968068701f25a1e29c8f7d8bd7` |
| 23 | 1230 | 490 | none | 585 | 113 | `4aa61de9c1793b40354540dbdc06bb771aa49795ddd2e5d694a5d90743bc0c8b` |
| 24 | 1285 | 516 | none | 607 | 137 | `0379f6fdc92f0380698f7e48fc3b8b4a9a6aeba4cbdd52f3765ca23a710992a2` |

For each row, exact Horner evaluation verifies

```text
R_M(r) != 0 mod p_M for every r in F_(p_M).             (5)
```

The SHA-256 digests of the complete residue-value tables are:

| \(M\) | complete gate-table digest |
|---:|:---|
| 21 | `365eeeafc09b137555bae4711574baa3b14c7d2b36c73bab9bf027c41535b0e7` |
| 22 | `9432915e09a30033831510b7832068ec47a2635147962af06dd5a5425f36b488` |
| 23 | `deb2e21c89809291ca19d728dc6a8682cbfd04b689d44936c1bd88006cd4c31c` |
| 24 | `315d7cea1008f6546286f580889b1d65efe4d144140bda83c5629482e8982129` |

If \(R_M(n)=0\) for an integer \(n\), reduction modulo \(p_M\) would
contradict `(5)` at \(n\bmod p_M\).  Combining `(2)--(5)` proves
`(1)`.

The witnesses

```text
127, 107, 113, 137
```

are frozen certificates, not a prime-selection law.

## 3. Exact sign geometry beyond the `4M` window

Direct exact evaluation through \(M^2\), together with a one-sign
Gregory--Newton vector at the displayed final base, gives:

| \(M\) | exact sign runs | GN positive-tail base | relation to \(4M\) |
|---:|:---|---:|:---|
| 21 | `0..70:+,71..133:-,134..441:+` | 134 | fail |
| 22 | `0..79:+,80..150:-,151..484:+` | 151 | fail |
| 23 | `0..89:+,90..169:-,170..529:+` | 170 | fail |
| 24 | `0..99:+,100..189:-,190..576:+` | 190 | fail |

The complete Newton-vector digests are:

```text
M21 d8d073a344a7649afee6d4373c65d6bf70d52662585aa6816850ae5ed41b8fca
M22 13494fc018c276ea8fedb0fc83e9aca50e66d533696f84ee982818d7a848fee7
M23 36c93dd5e27fb71436ceae3de415fe2b7e01d895bcfe743bc30498af008e8e7c
M24 081beb5f76430d613afabe15c1a98bf17b37163c937bab9ecfb37df4f5320945
```

Thus every new tail begins after \(4M\).  The modular gate, rather than
a bounded Gregory--Newton cutoff, is the uniform proof engine.

## 4. Exact replay and audit boundary

Run

```text
python -B 04-computation/gmc_first_gap_width_twenty_one_twenty_four_modular_depth_continuation_thm2959.py
python -O -B 04-computation/gmc_first_gap_width_twenty_one_twenty_four_modular_depth_continuation_thm2959.py
```

The two runs byte-match the stored transcript.  For every width the
companion verifies:

1. all interpolation values and three direct determinants outside the
   interpolation grid, at depths
   `deg(P_M)+4M-2`, `deg(P_M)+4M-1`, and `deg(P_M)+4M`;
2. positivity of `q200` and the fixed orientation of `c300`;
3. exact division by `c300*q200^5*B_M`, plus the exact width-21
   surplus `(n+17)`;
4. integral quotient coefficients, primitive content, core degree,
   and frozen coefficient digest;
5. every residue of the frozen root-free gate and its digest;
6. every integer sign through \(M^2\);
7. the complete one-sign Newton vector at the displayed tail base.

There are twelve direct outside-grid determinant controls in total.
The four-record digest is

```text
c45265c9537ba68f8ee8c73ff4a23dde3a634fd6e590dcb9ba96b72687655e66.
```

This theorem is independent of the retracted all-width
`nextprime(4(n+M))` finite-field rank shortcut.  That shortcut fails at
`M=10`, support `(0,5,6,10)`, `n=0`, `p=41`, where the modular rank is
only `34`.  Here each frozen prime tests one explicit characteristic-zero
cofactor after exact factor removal; no claim is made about generic
modular rank at a next prime.

The result covers exactly four first-gap supports and their translates.
It does **not** prove every support of widths \(21\)--\(24\), arbitrary
width, irreducibility or real-root-freeness of the cores, a uniform
prime rule, NC(2), or GMC(2).

The separate stable two-full-chart computation at width \(21\) is a
stronger characteristic-zero signal: its full determinant gcd has
nonnegative coefficients and positive constant term.  It is not used
here and remains outside this theorem pending its own audit.

**QED.**
