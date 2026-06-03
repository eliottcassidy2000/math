# LRC n=14 Carry-Fiber Conservativity, S611

The latest n=14 improvements have a common shape: they are all quotient maps,
and each one is useful only to the extent that it keeps the data its predicate
needs.

HYP-2166 packages the current n=14 proof surface as a quotient tower.  Its
carry seam is the piece this S611 diagnostic isolates.

HYP-2162/THM-407 folds the raw `C=27` shells by `G=<2,-1>` to the three
prime-3 strata `{gcd 1, gcd 3, gcd 9}`.  HYP-2164 classifies the
least-positive `Res_27` section and leaves only AP, `V*`, and nonprimitive
`2*AP` at the floor.  HYP-2165 says the 64 fixed self-converse classes carry
Redei parity, while the `C=27` owner layer carries cheap-pair and
positive-measure certificates.

The missing structure is the carry fiber.

If a lifted speed is written

```text
v = r + 27k,
```

then `27 == -1 (mod 14)`, so

```text
v == r - k (mod 14).
```

This is the whole lift/CRT problem in one line.  The `Res_27` shadow `r` is
not what the `n=14` clock sees.  The clock sees the residue together with the
carry.  So the least-positive representative is a section of the coimage, not
the conservative object itself.

S611 turns that into a small exact audit.  All `36` unit scalar lifts of AP
and `V*` remain at `M=1/14`, as they must by scaling invariance.  But after
projecting each one to the least-positive `Res_27` section, only three shadows
remain at the floor:

```text
AP,
V*,
2*AP.
```

The other `33` section shadows are strict.  That both validates HYP-2164 and
shows exactly why HYP-2164 cannot be the full lifted theorem: the section can
turn a floor lift into a strict representative by forgetting carry data.

The local carry check is also encouraging.  Over AP and `V*`, adding one or
two isolated `+27` carries always makes the row strict:

```text
AP  weight 1: min M=1/13
AP  weight 2: min M=1/12
V*  weight 1: min M=2/25
V*  weight 2: min M=1/12
```

So there is no obvious new local floor family hiding in the carry fiber.  A
lifted floor row must use a globally coherent carry pattern, plausibly one
cohomologous to a scalar AP or `V*` carry, or else it should route to the
HYP-2165 owner certificate or the Cprime/multiple CRT contradiction.

From the anti-Poisson viewpoint, this is a rigidity statement about
all-orders cancellation.  The base quotient says almost every section shadow
has positive safe measure above the floor.  The local carry audit says small
fiber defects also reopen the ground cell.  Therefore the residual cancellation
is not a generic shell-support feature; it is a rigid section-plus-descent
phenomenon.

The proof route I would try next is a carry cocycle theorem.  Treat

```text
K_i = (v_i-r_i)/27
```

as a `Z`-valued cocycle over the `Res_27` shell/owner quotient.  The base
residue gives pair-sum shell data, the carry gives the `n`-clock correction,
and endpoint owners impose local compatibility inequalities.  Then prove that
any cocycle preserving floor cancellation is either scalar AP/`V*`, a known
nonprimitive normalization, or it creates a bounded witness/certificate.

Artifacts:

- `04-computation/lrc_n14_carry_conservativity_s611.py`
- `05-knowledge/results/lrc_n14_carry_conservativity_s611.out`
- `05-knowledge/hypotheses/HYP-2167-lrc-n14-carry-fiber-conservativity.md`
