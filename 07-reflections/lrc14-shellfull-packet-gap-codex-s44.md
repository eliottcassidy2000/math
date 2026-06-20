# LRC14 Shell-Full Packet Gap - Codex S44

The useful shift in this session was to stop asking whether the shell-full
`2p1/5` tax merely survives bigger bounded boxes.  HYP-2669 already made that
look likely through `B=24`.  The better question is why the frontier does not
move.

The B30 scout says the answer is not hidden in raw speed size.  Once shell 1 is
full, the single high row above `1/3` is still the old dyadic-even B13 packet.
Rows with new speeds beyond `14` drop below `1/3`; rows beyond `24` drop below
`1/4` in the exact scan.

That makes the shell-full proof feel less like a global scalar inequality and
more like a packet-exclusion theorem:

```text
finite dangerous packet
new-speed packet decay
far-tail packet decay
```

The packet data also says what should be retained before scalarizing.  The B13
leader has small-denominator phase mass and a thick fold-target profile.  The
new-speed leader still has small denominators, but its fold reciprocal mass is
much lower.  So the missing invariant is probably not "small denominator" by
itself; it is small denominator plus fold-target concentration plus Glaisher
tower word.

This integrates cleanly with the repo's older fold and relation work.  HYP-2643
said fold multiplicity is target transport, not count.  HYP-2670 says the same
thing in the p1-tax quotient: the leader is not just relation-rich, it keeps the
right targets alive while the extended dyadic tower spreads them.

Concrete next proof move:

1. Certify `max(E')<=14` shell-full rows as a finite packet ledger.
2. Prove a new-speed inequality `max(E')>14 => Delta^+ <= p1/3`.
3. Try to strengthen the far-tail statement to `max(E')>24 => Delta^+ <= p1/4`
   or to a packet-dependent variant that implies it.

The shell-damaged side remains HYP-2661/HYP-2666.  The shell-full side now has
a visible shape.
