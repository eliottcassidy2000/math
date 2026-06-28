# LRC14 Expanded Residue / Owner-Support Repair

HYP-3311's residue-word success on the curated `31`-row bank was real, but it
was not the next theorem.  The enlarged HYP-2963 bank shows why.

The first failure is the one recent work prepared us for:

```text
P10+GW
GW-shell alias 12->132
```

same coarse packet, same nonunit residue word, different theorem exits.  This
is a genuine height leak, and `v2` / exact nonunit height repairs it.

But the more interesting failure appears immediately after that:

```text
petal 13->26
single swap 1->26
single swap 3->26
single swap 5->26
```

These rows agree on the nonunit side completely.  Same residue word.  Same
`v2` word.  Same exact nonunit height word.  Yet one row is
`unit-petal-named` and the others are `positive-Haar-open`.

So the enlarged-bank obstruction is not "height is missing" in any simple
sense.  The real missing coordinate is which endpoint-owner channels are
active at the positive boundary.  Once that owner-support word is added, the
mixed fibers disappear on every scanned bank up through
`single_limit=72`, `two_swap_limit=20`.

The frontier extension adds useful pressure: at `(72,20)` there are now two
height-persistent mixed fibers.  The old `petal 13->26` family remains, and a
new `petal 10->20` fiber collides with positive-open two-drop/add-20 rows even
after exact nonunit height is retained.  The owner-support word still splits
both.

This feels like the right hidden connection to preserve:

```text
residue data tells us which covering layer we are on;
height data tells us how far the covering branch has moved;
owner-support tells us which boundary channels actually carry the theorem exit.
```

The recent story was:

```text
HYP-3311: residue is enough on the curated bank
HYP-3402: expect height or owner-current to be next
HYP-3403: same-residue height debt is the stress target
HYP-3404: residue-word breakpoint theorem is the top finite-lemma lead
HYP-3405: AP-collar finite lemma certifies the unit-height companion case
HYP-3406: yes, but the stronger stable leak is owner-side, not height-side
HYP-3407: turn this sidecar chain into recursive cut/signature tests
```

In other words: tropical height walls explain the first leak, but endpoint
owners explain the first leak that survives tropical repair.
