# LRC14 Fixed-Margin Labelled Packets

## Claim

The HYP-2961 family/sporadic classifier needs one more structural layer:

```text
family = fixed-margin labelled packet class
sporadic = bounded singleton after family parameters are removed
```

I added HYP-2962 and a small classifier script:

```text
04-computation/lrc14_fixed_margin_labelled_packet_classifier_codex_s150.py
05-knowledge/results/lrc14_fixed_margin_labelled_packet_classifier_codex_s150.out
```

## External Pattern

I used arXiv:2606.22636 only as a proof-pattern import.  That paper proves a
fixed-margin binary swap-chain spectral-gap result by comparing to a two-row
heat-bath chain, reducing to three rows, and splitting the three-row model into
a scalar count sector plus Johnson harmonic non-scalar sectors:

```text
https://arxiv.org/abs/2606.22636
https://arxiv.org/html/2606.22636
```

For LRC14 the analogous split is:

```text
fixed margins            -> labelled packet families
swap moves               -> packet-preserving mutations
scalar count sector      -> qdiv / q-cover data
Johnson-like sectors     -> boundary owners, source spectrum, C27/K33, gK8
```

No result from that paper is being used as an LRC theorem.

## Packet Signature

For each row the script records:

```text
qdiv bucket
q-cover vector for d=2..14
U14 apex zero profile
U14 boundary profile
C27 shell state counts
packet atom keys
exact M and argmax denominators
strict Haar-open mass
C27 transfer and S145 packet route
```

Those margins define a labelled packet class.

## Run Result

Representative run:

```text
rows audited                         = 18
fixed-margin classes                 = 16
shared family signatures             = 2
singleton sporadic signatures        = 14
strict counterexamples found         = 0
qdiv>14 representatives              = 5
qdiv>14 zero-open representatives    = 0
```

Bucket counts:

```text
D0 direct q-witness family                 4
D2 positive Haar-open row                  2
D2 positive covering Haar-open family      5
D3 positive unit-petal/GW strip            3
D4 positive K33/state-lift labelled row    2
S0 equality sporadic AP/GW                 2
```

Shared packet classes:

```text
12->28 and 12->56:
  direct-q apex decoy family.

12->84 and 12->168:
  qdiv>14 covering lcm-tail family, discharged by positive open mass.
```

The second one is the signal: covering rows should be classified as fixed-margin
packet fibers, not as isolated rows.

## Live Buckets

No row in this bank landed in the dangerous zero-open covering core.  In
HYP-2961 terms:

```text
L1 apex-multiple residual:
  not seriously tested; needs a many-14-multiple bank.

L2 genuine-wide zero-moment:
  moment image not computed here.

L3 bounded covering core:
  represented by covering repairs and lcm tails; all positive-open.

L4 K33 zero-open state lift:
  no zero-open representative found.

L5 unnamed source kernel:
  no zero-open non-K33 representative found.
```

## Theorem Target

Build `P(S)` with rows/features:

```text
q-cover obligations,
U14 apex contacts,
C27 shell states,
Haar open/boundary front,
source-spectrum labels,
K33/state-lift flags,
gK8/L_y moment coordinates.
```

Then prove:

```text
Every strict LRC14 counterexample packet lies in L1, L2, L3, L4, or L5.
Every packet class outside those live buckets has a named discharge.
```

The next real search should range over packet-preserving swaps inside a fixed
margin class and add the missing gK8/L_y coordinate.

