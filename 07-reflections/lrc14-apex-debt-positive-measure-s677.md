# LRC14 Apex-Debt Positive Measure (S677)

Prompt: spend a long session trying to improve the next target

```text
prove no new normalized p_0=0 wall exists in the Res_27 carry/owner fiber.
```

The target improved by splitting the wall problem into the piece that matters.

## The Split

AP and Vstar are `p_0=0` walls, but they are no-multiple rows.  So they are
not dangerous for LRC14: `t=1/14` and the other unit `n`-clock times already
witness them.

`2AP` is also a `p_0=0` wall, but it is nonprimitive.  Its witnesses are just
the doubled preimages of the AP witnesses, and it disappears after gcd
normalization.

So the residual is not "all p_0=0 rows."  It is:

```text
primitive rows with a multiple of 14 and p_0=0.
```

In carry coordinates `v=r+27k`, this is the apex-debt congruence:

```text
14 | v  <=>  k == r (mod 14).
```

That is much better.  The proof target becomes:

```text
primitive apex debt => p_0 > 0.
```

This is the old C' branch in sharper language: the multiple of `14` kills the
`n`-clock endpoint witnesses, but it should open positive Lebesgue measure
elsewhere.

## Exact Audit

S677 adds `04-computation/lrc14_apex_debt_lebesgue_s677.py` and stores
`05-knowledge/results/lrc14_apex_debt_lebesgue_s677.out`.

It checks exact rational `p_0` on four coherent carry families over AP and
Vstar:

- one coordinate carried to its first multiple of `14`;
- that apex debt plus one side carry;
- contiguous carry blocks of height `1` or `2`;
- small affine carry laws `k=a*r+b (mod 14)`.

Across `720` probes, the primitive multiple branch has:

- `414` primitive multiple rows;
- `0` `p_0=0` walls;
- `414` positive-measure rows.

The minimum primitive-multiple safe measure is

```text
181/28028
```

at the AP row

```text
(1,2,3,4,5,7,8,9,10,11,12,13,168).
```

This is a genuinely useful witness for the shape of the residual: it has a
multiple of `14`, gcd `1`, no surviving `n`-clock endpoints, no gcd-scaled
endpoint rescue, and still has four positive safe intervals.

## Why This Improves HYP-2252

HYP-2252 said "no new normalized wall."  S677 says most walls are not relevant
to the proof risk:

- no-multiple wall: handled by the old clock;
- nonprimitive wall: handled by gcd reduction;
- primitive multiple wall: must be ruled out.

So the next theorem should be narrower and more plausible:

> Primitive apex debt pays positive measure.

The exact audit does not prove the theorem, but it removes a lot of fog.  It
also gives the right side channel: the proof should track where the carry
congruence `k_i=r_i mod 14` lands in the owner ledger.

## Tournament Analysis

Vertices are proof filters, not runners.  I considered runners, clocks,
gcd-scaled clocks, apex congruence sites, carry vectors, safe intervals, owner
obligations, and proof obligations.

The filter tournament is transitive:

```text
apex_congruence_debt
> primitive_multiple_positive_measure
> owner_private_derivative
> gcd_scaled_endpoint_wall
> no_multiple_n_clock
> raw_res27_shadow
> raw_first_moment
```

This order is a good guide.  First name the apex congruence.  Then prove
positive measure on the primitive multiple branch.  Use owner-private deletion
as the side channel.  Only then drop down to the old endpoint-wall cases.

## Next Move

Attach HYP-2241's owner-private bit directly to apex-debt sites:

1. For every `k_i=r_i mod 14` site, record which D/U/N obligations it owns
   before and after deletion.
2. Test whether primitive apex-debt rows with the same visible `Res_27` shadow
   but different owner-private ledgers ever have small `p_0`.
3. Try to prove a deletion lemma: if the apex site is private for any
   obligation, its deletion opens a safe interval; if it is never private, the
   remaining owners force a cheap-pair or block-channel positive-measure
   certificate.

That would turn the computational split into a proof mechanism.
