---
id: HYP-2674
title: LRC14 uniform Delta tail via same-sign packet alignment
status: OPEN; exact scout after HYP-2653d correction
source: codex-2026-06-20-S46
depends_on:
  - HYP-2673
  - HYP-2671
  - HYP-2653d
  - HYP-2653
related:
  - HYP-2672
  - HYP-2670
  - HYP-2661
  - HYP-2644
  - OPEN-Q-108
---

# HYP-2674 - Uniform Delta Packet Alignment

## Claim

The corrected HYP-2653d far-tail target

```text
sup_{max(E')>B} Delta_w(E',w) <= cap_k - Q(k-1)
```

should be attacked by classifying same-sign one-missed-sector packet alignments,
not by a generic discrepancy estimate first.

For `E=E' union {w}`, the KPS/codex bridge writes

```text
Delta_w = p0(E) - Phi(E')
        = (1/w) sum_{1-missed cells c}
          [G0(w*b_c-s_c/7)-G0(w*a_c-s_c/7)].
```

HYP-2674 records the packet sign word over missed sectors `s=1..6`.  The known
dangerous rows are not hidden cancellation failures: every packet is positive.
Thus the finite obstruction is a `++++++` alignment pocket.  The tail proof
should show that after this pocket, either sign changes appear or the total
positive packet mass is too small to threaten the uniform `Delta_w` margin.

## Exact Scout

Script:

```text
04-computation/lrc14_uniform_delta_packet_alignment_codex_s46.py
```

Output:

```text
05-knowledge/results/lrc14_uniform_delta_packet_alignment_codex_s46.out
```

The k=9 far-tail margin is

```text
cap_9 - Q(8) = 129643/980980 ~= 0.132157.
```

Named rows:

```text
finite B13 packet:
  E'=(0,1,2,4,6,7,8,10), w=12
  Delta_w = 997/5880 ~= 0.169558
  sign word = ++++++
  note: exceeds the k=9 far margin, but it is a finite max(E')=10 pocket row.

HYP-2671 dyadic s=4 alignment:
  E'=(0,1,2,4,8,12,16,20), w=24
  Delta_w = 457/3920 ~= 0.116582
  margin gap = 12223/784784
  sign word = ++++++

non-shell warning row:
  E'=(0,2,3,5,6,15), w=18
  Delta_w = 11/315 ~= 0.034921
  sign word = ++++++
```

Dyadic family scan for

```text
E_s={0,1,2,4,8,3s,4s,5s},  s=3..120
```

gives:

```text
w=6s:
  global best        s=4,  Delta=457/3920
  best after s>20    s=22, Delta=2539/64680 ~= 0.039255
  best after s>40    s=60, Delta=127/3675 ~= 0.034558

w=10s:
  global best        s=16, Delta=673/47040 ~= 0.014307
  best after s>20    s=28, Delta=841/82320 ~= 0.010216
  best after s>40    s=44, Delta=19/2310 ~= 0.008225
```

The important point is not just that these are below margin.  The `w=6s`
family has its only threatening spike at `s=4`, while the post-20 tail has
more than `0.0929` absolute k=9 margin.  This supports HYP-2653d's empirical
`B=20` safety as a finite-alignment cutoff, not a mysterious numerical cliff.

## Incoming KPS Addendum

After HYP-2674 was first pushed, incoming KPS work added the finite-half and
per-sector-Koksma side of the same picture.

Finite half:

```text
05-knowledge/results/lrc14_finite_half_span14_kps.out
```

certifies `p0(E)<=cap_k` for all primitive `E` with `max(E)<=14`,
`k=8..12`, with zero violations.  For example:

```text
k=8:  sets=3431, maxp0=0.32721, cap=0.38146
k=9:  sets=3003, maxp0=0.41616, cap=0.49426
k=10: sets=2002, maxp0=0.50448, cap=0.60440
```

Per-sector Koksma:

```text
04-computation/lrc14_ck_per-sector-koksma_kps-Sx-wf.py
05-knowledge/results/lrc14_ck_per-sector-koksma_kps-Sx-wf.out
```

proves and verifies the exact per-sector telescope

```text
w*Delta_w = sum_s sum_{[a,b] in exact-{s} runs}
            [G0(w*b-s/7)-G0(w*a-s/7)].
```

It also verifies that the tempting full "sector-s-missed" telescope is wrong:
cells with `|miss|>=2` are double-counted.  The rigorous Koksma/BV bound is

```text
w|Delta_w| <= (6/49) sum_s |R_s| <= (6/7) sigma(E').
```

Most importantly, KPS proves that any standalone `w|Delta_w|` constant is
false: for

```text
E'_M={0,1,2,3} union {M,M+1,M+2,M+3}, w=22M,
```

the diagnostic `w|Delta_w|` grows like `0.08M`.  This does not kill the LRC14
route because the wide base has small plateau/p0.  The real analytic closer is
a joint statement:

```text
p0(E)=Plat(E')+Delta_w,
wide E' => small Plat/p0,
bounded E' + large w => sigma(E')/w decay.
```

So HYP-2674 should be read as the bounded-near-plateau packet-alignment
subproblem.  It classifies why the finite pocket and the HYP-2671 dyadic row
are dangerous.  The very wide resonant branch may have huge `w*Delta_w`, but it
is routed through the KPS joint plateau/sigma argument rather than through a
standalone packet-alignment cap.

## Proof Route

1. Define the one-missed-sector packet sign word
   `sgn_s(E',w)=sgn(sum_{c:miss(c)=s} G0(w*b_c-s/7)-G0(w*a_c-s/7))`.
2. Prove a finite alignment lemma: rows with large positive `Delta_w` must have
   the `++++++` packet word plus a dyadic/fold address matching the B13 or
   HYP-2671 pockets.
3. Prove a tail lemma after `max(E')>B`:
   either some packet is nonpositive, giving signed cancellation, or the
   same-sign packet mass has the dyadic-tail bound seen above.
4. Splice this with HYP-2673 and the KPS per-sector result: finite pocket below
   `B`, shell-full `p1` taxes where applicable, bounded-base `sigma(E')/w`
   decay, and wide-base small-plateau control.

## Tournament Analysis

Vertices:

```text
finite_B13_packet
dyadic_s4_alignment
nonshell_warning
dyadic_w6_tail_after20
dyadic_w10_tail_after20
```

Pairwise observable: larger positive `Delta_w` points to the riskier proof
obligation.

Hamiltonian path from the script:

```text
finite_B13_packet
> dyadic_s4_alignment
> dyadic_w6_tail_after20
> nonshell_warning
> dyadic_w10_tail_after20
```

Score histogram is `[0,1,2,3,4]` and there are no directed `3`-cycles in this
risk order.

Challenged assumption: the useful tournament vertices are runners or residues.
For the corrected far-tail lemma, the better vertices are proof obligations
indexed by packet sign words and finite alignment pockets.  This preserves the
LRC predicate `Delta_w <= cap_k-Q(k-1)` while deliberately forgetting raw
runner labels.

## Honest Status

LRC(14) is not proved.  HYP-2674 is a proof-route improvement: it turns the
near-plateau part of the corrected `Delta_w` problem into a finite same-sign
alignment problem plus a signed tail problem.  Incoming KPS per-sector work
adds the complementary wide-base route: huge `w*Delta_w` resonances are harmless
when the plateau itself is small.
