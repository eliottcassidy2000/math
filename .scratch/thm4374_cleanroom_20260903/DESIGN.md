# THM-4374 clean-room design freeze

**Date:** 2026-09-03  
**Status:** independent derivation fixed before reading
`04-computation/lrc14_seventeen_step_metric_exit_observability_thm4374.py`.

## Scope and inheritance

The source universe is exactly the THM-4365 family with odd `P>=11019`.
Only its first missing-component metric coordinate is observed.  The closest
proved mechanisms are THM-4365's centered-residue law and THM-4367's
primitive-ray scale fibres.  The canonical hostile is THM-4367's equal-exit
four-tuple that splits under `P -> P+2`.  The corrected near miss is that a
quotient for one output need not be a congruence for the operation.  The
least-used sidecar is the first return of the centered phase to the strict
active interval.

Anchor / Niche / Wildcard:

- Anchor: audit the proposed horizon-17 observer and every quantifier.
- Niche: audit open-tooth boundary phases separately from strict activity.
- Wildcard: reconstruct `P` directly from two suitably placed active outputs,
  without using the proposed program's enumeration strategy.

Concept board:

```text
centered phase | active interval | primitive scale fibre | return time
metric word | direct rational decoder | shift congruence
```

The applicable method cards are "Controlled forgetting and unlabeled
quotients require a sidecar" and "Re-evaluate a certificate after every
fibre-changing operation".  The scale/return-time sidecar is exactly what the
static metric quotient loses.

## Independent derivation

Set

```text
A=3371, S=1303, M=14A=47194, N=M/2=23597.
```

For the centered remainder `SP=Mn+rho`, `-M/2<rho<=M/2`, put `c=A-rho`.
For odd `P`, `c` is even.  With `x=c/2 mod N`, strict activity is precisely
`I={1,...,3370}` and the shift is the rotation `x -> x-S mod N`.

Inside a current active THM-4367 fibre, write

```text
c=ag, P=bg, a+Sb=A*kappa, g*kappa=1 (mod 14).
```

For `c>=2608`, time one remains active and

```text
E(bg+2)=g*kappa / (14(bg+2)).
```

Cross multiplication makes this injective in `g`.  For `c<=2606`, times
1 through 15 are inactive.  The once-wrapped returns are

```text
c_16=c+5498, c_17=c+2892.
```

Thus time 16 is active iff `c<=1242`, while time 17 is active for all
`c<=2606`.  At a common return time `t` the output is

```text
(g*kappa+14)/(14(bg+2t)).
```

Equality at two distinct scales reduces to `t*kappa=7b`, impossible for
`t=16,17` because `kappa` is a unit modulo 7.  This gives the claimed
partitions without any bounded search.

The fibre-size bounds follow from one residue class of `g` modulo 14:

```text
ag<=6740                    -> at most 241;
ag<=2606                    -> at most 94;
1244<=ag<=2606              -> at most 49;
horizon at least 17         -> 1.
```

The proposed witness `(a,b,kappa)=(2,47595,18397)` has
`g=1,15,...,3361`, so it realizes all four numbers.

For the global claim, the preimage intervals `I+tS`, `0<=t<=16`, have
unwrapped union `[1,3370+16*1303]=[1,24218]`, which covers all `N=23597`
phases.  The first-hit census should be `3370, 15*1303, 682`.  If first entry
is delayed, its phase `y` lies in `[2068,3370]`, so the following phase
`y-S` is active too.  Two consecutive active outputs with shifted parameter
`Q` and excesses `d_0,d_1` above `S/M` recover

```text
Q = 2(S+M*d_1)/(M(d_0-d_1)).
```

If the current phase is active but its successor is not, the next active
time is 16 or 17 after one wrap.  From current excess `d_0` and return excess
`d_t`, recover

```text
P = (M-2tS-2tM*d_t)/(M(d_t-d_0)).
```

The possible zero denominator is exactly the already excluded equation
`t*kappa=7b`.  Hence `W_17` determines `P` directly.  The stated
`(253031,258645)` pair will be checked as a horizon-16 collision and a
horizon-17 split.

Finally, the kernel of any output-compatible forward-shift state relation is
an equivalence relation stable under iteration.  Iterating 17 times gives
equal `W_17`, hence equal parameters.  This is the congruence corollary; it is
logical, not a finite-enumeration inference.

## Independent verifier plan

The verifier will use only Python's standard library and no repository
imports.  It will:

1. enumerate all 23,597 rotation phases, including both inactive boundaries;
2. construct an actual coprime `(b,kappa)` representative for every finite
   structural type `(even a, unit kappa mod 14)`;
3. exhaust all 281,073 unordered scale pairs in the resulting stronger full
   structural fibres and compare literal rational output words;
4. compute every horizon partition and the sharp witness;
5. exercise the direct decoder on two complete phase periods;
6. check the horizon-16 hostile and exact return fractions;
7. run normal, `-O`, `-I`, and fixed-hashseed modes and require identical
   raw-LF byte streams.

No canon or navigation file is in scope for editing.  LRC(14) remains open.
