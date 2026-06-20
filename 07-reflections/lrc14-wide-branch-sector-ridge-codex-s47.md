# LRC14 wide-branch sector ridge - codex s47

The useful correction from this session is that `span>14` is too coarse.
The exact sector scout separates boundary rows from true-wide rows.

For `k=9` in the exact `{0}+8` box through `20`, the leader with `span>14`
is not a wide-base row at all:

```text
E=(0,2,4,6,8,10,12,14,15)
p0=437/1176
cap_9-p0=20627/168168
second-largest=14
```

The true-wide leader in the same box is lower:

```text
E=(0,4,6,8,10,12,14,15,16)
p0=321/980
cap_9-p0=11681/70070
second-largest=15
```

That changes the proof shape.  The KPS comfortable-margin route should not be
phrased as a single spread lemma first.  It wants an ordered split:

```text
finite span<=14 check
boundary collar: second<=14, one-far AP-like rows
true-wide sector deficit: second>14, Freiman/GAP/state-word scaffold
post-25 packet tail: KPS packet-mass decay
```

The true-wide rows that remain high are still not random.  They have small
sumset excess and visible GAP structure: even progressions, shifted intervals,
or short cluster packets.  That matches the user's Freiman/additive-energy
prompt: high `p0` is high additive energy plus low state-word entropy, not a
generic large-spread phenomenon.

The dyadic Delta rows also stop looking mysterious after direct `p0`
evaluation.  The HYP-2671 full dyadic row has

```text
p0=29/112
cap_9-p0=3769/16016
```

so it is dangerous for a decoupled `Delta_w` proof but harmless for the direct
sector predicate.  This is exactly the Plat/Delta entanglement: the same packet
alignment that raises Delta lives on a small-plateau base.

I do not think this proves LRC(14), but it narrows the crux.  The next theorem
should either prove a compression monotonicity for the boundary collar, or prove
that true-wide rows with `second>14` and high `p0` must compress to one of the
finite GAP scaffolds already visible in the exact box.

Incoming THM-546 slots into this cleanly.  Its

```text
|Delta_w| <= kappa V(E')/(pi^2 w)
```

is the rigorous gapped one-far estimate.  The rows found here are the ungapped
or near-ungapped residues where `w` is not large compared with the base, so the
proof still needs a compression or dimension argument.  In other words:
THM-546 handles distance in scale; HYP-2675 identifies what is left when scale
distance is absent.
