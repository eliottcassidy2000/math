---
id: THM-2957
title: "First-gap width-fifteen-to-twenty modular depth ladder"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On each of
  the six first-gap supports
  (0,1,2,M), 15<=M<=20, THM-2949's fixed rank-35 cofactor is
  nonzero at every integer depth.  Exact division by a stable
  floor-law negative-root factor leaves a primitive core with a
  frozen root-free finite-field gate.  This is not a complete-width
  or arbitrary-width theorem.
source: codex-gmc-first-gap-modular-ladder-2026-07-29
audit: >
  Independent hostile audit accepted the floor-law factor, exact
  divisions, integral primitive cores, coefficient order, six complete
  finite-field gates, direct outside-grid determinants, sign/tail
  boundaries, normal/optimized replay, LF hashes, and first-gap-only
  scope.
depends_on:
  - THM-2949-fixed-rank-thirty-five-cofactor-newton-atlas
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
related:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2955-width-twenty-fixed-fifth-compound-mod-97-gate
  - THM-2956-koszul-gale-fixed-fifth-compound-exchange
script: 04-computation/gmc_first_gap_width_fifteen_twenty_modular_depth_ladder_thm2957.py
output: 05-knowledge/results/gmc_first_gap_width_fifteen_twenty_modular_depth_ladder_thm2957.out
script_sha256: d8157df94c9c14f6607905819dae00ea47ae260344d331a974e90bef152655ed
output_sha256: 1472caa06836ed6e4784315b16b2377f006589a21db4bb422cc3d83a17253584
thm2949_dependency_sha256: 9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54
hash_basis: LF-normalized bytes
---

# THM-2957 -- first-gap width-fifteen-to-twenty modular depth ladder

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

For \(15\leq M\leq20\), let \(P_M(n)\) be the fixed
\(20Q+10C+5F\), rank-thirty-five cofactor constructed in THM-2949 for
the normalized support

```text
(0,1,2,M).
```

Then

```text
P_M(n) != 0 for every integer n>=0 and every 15<=M<=20.       (1)
```

Consequently the PROVED conjugate-pair rank gap of THM-2947 gives full
rank of the degree-seven Macaulay map at every physical depth for each
of these six supports.  In particular, first-window SFC(4) holds for
every translate of

```text
(0,1,2,M),  15<=M<=20.
```

The certificate is deliberately discrete.  Each sampled integer
sequence has two sign changes.  What excludes integer zeros is a
finite-field gate on the remaining high-degree core, not real
positivity.

## 1. The stable negative-root factor

Write

```text
a=floor(M/3), b=floor(M/2), c=floor(2M/3).
```

Define

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

The exact companion constructs \(P_M\) from the live THM-2949
constructor and verifies exact polynomial division

```text
P_M = c300*q200^5*B_M*core_M.                          (3)
```

It also verifies

```text
deg q200=M-1,             every coefficient of q200 is positive,
deg c300=2M-1,            every coefficient of -c300 is positive.
```

Thus `q200`, `c300`, and `B_M` are nonzero for every real \(n\geq0\).

Primitive-normalize `core_M` with positive leading coefficient and call
the result \(R_M\).  The exact degree and digest census is:

| \(M\) | \(\deg P_M\) | \(\deg B_M\) | \(\deg R_M\) | SHA-256 of primitive coefficient vector |
|---:|---:|---:|---:|:---|
| 15 | 790 | 313 | 378 | `b6fea710768bfa6be0c7c4c529c1698f1a0ddff0dbf2e95e568272a020823423` |
| 16 | 845 | 336 | 403 | `505e5c935a11ec2beff10d22d16768a0bfe1878b37c50ac19f8ada852b32950b` |
| 17 | 900 | 356 | 431 | `447b5b527c3b230aebb088a9733ff0480e1c8822f3ffe32deecbcaec463669bc` |
| 18 | 955 | 382 | 453 | `75a1081c87b3322fa3f81cd6bd35520cf963223e141a1dc84ca020395f801148` |
| 19 | 1010 | 401 | 482 | `86c1a256e878afa7b9890818430b77d87134bcdb68654699cb0bb385fc5a6acf` |
| 20 | 1065 | 425 | 506 | `a0176d5c76c67c931de4217da581d1599e29d573514917f416638c1d2497355a` |

The cofactor degree is uniformly \(55M-35\).  Formula (2), rather than
an expensive factorization of \(P_M\), isolates every forced
negative-root factor needed by the argument.

For these widths its degree and the remaining core degree have the
closed floor forms

```text
deg B_M = 19M+2 floor(M/3)+4 floor(M/2)+floor(2M/3)-20,
deg R_M = 29M-9-2 floor(M/3)-4 floor(M/2)-floor(2M/3). (4)
```

These are degree identities only.  They do not promote the observed
factor pattern to widths outside the stated six-member bank.

## 2. Six finite-field root gates

For each row below, exact Horner evaluation verifies

```text
R_M(r) != 0 mod p_M for every r in F_(p_M).             (5)
```

| \(M\) | \(p_M\) | SHA-256 of the complete residue-value table |
|---:|---:|:---|
| 15 | 83 | `5c3f2d2af0bf23808b9a7b0fc023a38eb3d18d57c49e86e58f6c5bdd074dc26e` |
| 16 | 71 | `a1f03010ddecd4296d9ef80ccc6cb9503c3d4a7dfc33fb153dc23cd5ce9a5208` |
| 17 | 107 | `bce5b1c385b75d3dc865c89672f64a120ca926d3a5bb4885f2f49d37b3f75242` |
| 18 | 73 | `7ed43f14e61981c1c0cec8fd100ffe811639fc09c6e2247d65840d61429007b2` |
| 19 | 79 | `0738bb31e61f4d2b2d77bc91a2c14599beb3f64ecfc1dfdb28cdc03bf515f59c` |
| 20 | 97 | `81743f0f0a65a8c3d3b92e313b8ee73b3376cd656b9c9943d93be34f07181978` |

If \(R_M(n)=0\) for an integer \(n\), reducing modulo \(p_M\) would
contradict (5) at the residue \(n\bmod p_M\).  Hence \(R_M\) has no
integral root.  Combining this with (2)--(3) proves (1).

The primes are witnesses, not a claimed law: the sequence is

```text
83, 71, 107, 73, 79, 97.
```

In particular width 17 already refutes a naive monotone or
approximately-\(5M\) prime rule.

## 3. Exact sign geometry and the failure of a uniform `4M` cutoff

Direct exact evaluation through \(M^2\), together with a one-sign
Gregory--Newton tail vector at the displayed final base, gives:

| \(M\) | exact sign runs | GN positive-tail base | relation to \(4M\) |
|---:|:---|---:|:---|
| 15 | `0..28:+,29..51:-,52..225:+` | 52 | pass |
| 16 | `0..33:+,34..62:-,63..256:+` | 63 | pass |
| 17 | `0..39:+,40..74:-,75..289:+` | 75 | fail |
| 18 | `0..46:+,47..87:-,88..324:+` | 88 | fail |
| 19 | `0..53:+,54..101:-,102..361:+` | 102 | fail |
| 20 | `0..61:+,62..116:-,117..400:+` | 117 | fail |

Thus width 17 is the first failure of the `4M` envelope in this
six-width ladder.  Width 20 is a sharp hostile to any attempt to keep
that cutoff: every base \(r\leq80\) lies before the eventual positive
tail and the cofactor takes both signs after \(r\), so no one-sign
Newton vector can start there.

This identifies the structural division of labor:

```text
explicit B_M  = stable forced negative-root geometry;
R_M           = integer two-sign-change high core;
mod-p gate    = integer-depth nonvanishing;
GN tail       = optional sign geometry, not the uniform proof engine.
```

## 4. Replay and audit boundary

Run

```text
python -B 04-computation/gmc_first_gap_width_fifteen_twenty_modular_depth_ladder_thm2957.py
python -O -B 04-computation/gmc_first_gap_width_fifteen_twenty_modular_depth_ladder_thm2957.py
```

The two runs byte-match

```text
05-knowledge/results/gmc_first_gap_width_fifteen_twenty_modular_depth_ladder_thm2957.out
```

and the six-record digest is

```text
cc826193b445e6ff27ee990d06287df417d6fe8242fd7116d86f8e478d842017.
```

The companion imports

```text
04-computation/gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py
```

and no scratch constructor.  For every width it verifies:

1. all interpolation values and three direct determinants outside the
   interpolation grid, at depths
   `deg(P_M)+4M-2`, `deg(P_M)+4M-1`, and `deg(P_M)+4M`;
2. exact division by `c300*q200^5*B_M`;
3. primitive content, core degree, and frozen coefficient digest;
4. every residue in the root-free prime table and its frozen digest;
5. every integer sign through \(M^2\);
6. the complete one-sign Newton vector at the displayed tail base.

This is a finite exact certificate for the six stated supports.  It
does **not** prove that the high cores are irreducible, real-root-free,
or members of a standard orthogonal-polynomial family, and it does not
claim a prime-selection formula beyond the frozen witnesses above.

An independent reconstruction replayed the canonical companion against
the stored transcript, rebuilt all six quotients from the promoted
THM-2949 constructor, and checked separately that every core is an
`fmpz_poly`.  Thus primitive normalization does not truncate rational
coefficients.  For each width it also compared direct evaluation of the
primitive quotient with Horner evaluation at residues `0`, `1`,
`floor(p_M/2)`, and `p_M-1`; all `24` comparisons agree.  The audit
rechecked the ascending coefficient convention, all exact degrees and
digests, and the distinction between integer nonvanishing and real-ray
positivity.

**QED.**
