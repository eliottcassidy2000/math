---
id: HYP-1989
status: OPEN
source: codex-2026-06-01-S518
related:
  - HYP-1988
  - HYP-1987
  - HYP-1986
  - HYP-1981
  - THM-385
  - THM-384
  - THM-381
  - THM-380
  - THM-369
---

# HYP-1989: n=145 CRT aperture proof route

## Statement

For total LRC denominator

```text
n = 145 = 5 * 29,
```

the most promising whacky proof route is a CRT unit-wall aperture argument.

At any unit wall

```text
t = a/145,  gcd(a,145)=1,
```

every moving speed whose residue is nonzero modulo `145` is already safe from
the observer, because its circular distance is at least `1/145`.  The only
blockers at such a wall are speeds divisible by `145`.

Thus:

1. If no speed is divisible by `145`, THM-369 with `q=145` gives a lonely time.
2. If a `145`-divisible speed exists, the `144` moving-runner budget prevents
   the speed set from also occupying all `144` nonzero residues modulo `145`.
3. Therefore at least one antipodal nonzero residue pair `{r,-r}` is not fully
   occupied.
4. Choosing a unit wall with `a^{-1}=r` makes one of the two boundary sides
   unpinned, creating a one-sided aperture adjacent to the unit wall.

The open proof target is the zero-residue embryo:

```text
speeds divisible by 145, scaled down by 145.
```

The conjecture is that for every primitive n=145 speed system, some unit-wall
aperture lets this zero-residue embryo move out of the forbidden cap before the
pinned side collapses.  If it cannot, the obstruction should be expressible as
a labelled endpoint-pressure core contradicting the THM-380 trilemma.

## Evidence

S518 adds the atlas script:

```text
04-computation/lrc_n145_whacky_reframes_s518.py
05-knowledge/results/lrc_n145_whacky_reframes_s518.out
```

The script records the basic n=145 arithmetic:

```text
145 = 5 * 29
moving runners = 144
threshold = 1/145
phi(145) = 112
nonzero nonunits = 32
```

Model-family audits show the denominator-trap behavior:

```text
initial_1_to_144:
  q=145 blockers=0, so THM-369 solves at denominator 145.

unit_wall_with_one_zero:
  q=145 blockers=1, but one nonzero residue is missing; antipodal pairs
  full/half/empty = 71/1/0.

lcm_spike_units:
  sieve_complete=True with q=145 blockers=1; the same one-sided aperture
  profile appears.
```

The whacky route tournament ranks the proof languages:

```text
1. unit_wall_aperture
2. zero_residue_embryo
3. almost_source_tunnel
4. crt_5_29_two_moons
5. observer_score_descent
```

This ranking is not a proof, but it exposes the concrete local problem:
source reachability at n=145 may be reduced to escaping a unit-wall aperture,
not to enumerating any A000568-scale target.

## Whacky Reframes

### 1. Unit-Wall Aperture

The wall `a/145` is almost a source by construction.  Nonzero residues are
safe; zero residues are blockers.  A missing antipodal boundary side creates
the aperture.

### 2. Zero-Residue Embryo

The speeds divisible by `145` form a smaller speed system sitting at the
observer at every unit wall.  The local perturbation problem is to move this
embryo out of the forbidden cap while preserving the nonzero residue runners.

### 3. Almost-Source Tunnel

By THM-385, almost-source means exactly one blocker.  If the aperture can
reduce the zero-residue embryo to one blocker, the remaining state may be
attackable by side-defect and pressure repair.

### 4. CRT Two Moons

Because `145=5*29`, there are two large denominator moons: mod `5` and mod
`29`.  The `q=5`, `q=29`, and `q=145` rational traps are nested CRT gates.
Counterexample-shaped systems must block all three.

### 5. Paley Failure as Data

Both `5` and `29` are `1 mod 4`, so the Paley tournament shortcut does not
directly apply.  That failure may be useful: the CRT residue space is naturally
bipartite/signed rather than Paley-regular.

## Predictions

1. Any n=145 near-counterexample must be denominator-sieve complete and contain
   at least one `145`-divisible speed.
2. Such a system always has a unit-wall aperture in the residue sense.
3. The real obstruction is not residue occupancy; it is dynamic aperture width
   versus the scaled zero-residue embryo.
4. Observer-score descent should be the right meter: the local task is to move
   from a positive blocker layer to the source layer.
5. If the embryo cannot escape, its failed escape should produce labelled
   endpoint-pressure data, not just a scalar near miss.

## Next Tests

1. Implement an exact local aperture solver around `a/145` walls for moderate
   speed families.
2. Add n=145 model rows to the observer-score and endpoint-pressure scripts.
3. Separate the zero-residue embryo into quotient speeds `v/145` and run the
   THM-369 sieve recursively on that quotient.
4. Search for a dual certificate: aperture-pinned boundary runner plus
   zero-embryo endpoint core.
