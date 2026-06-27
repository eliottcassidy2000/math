# LRC14, THM-079, and the STAR obligation

The THM-079 analogy is real, but it should not be allowed to hide the last
logical step.

The useful template is:

```text
THM-079:
  component/multiplicative reduction
  -> one strong atom
  -> forcing adds odd-cycle mass, forbidding H=21.

LRC14:
  q-witness / dilation / large-speed boundary-state reduction
  -> one bounded covering atom
  -> prove the atom is AP/GW boundary-only, or lifts to forbidden H=7/K3.
```

The apex-7 observation is strong but conditional.  A covering set contains a
multiple of `14`, so every denom-14 point is blocked.  This excludes AP/GW-style
apex equality rows from the covering core.  It does not by itself prove the
covering core has no off-apex witness or no sub-tight row.

S120 correction: THM-568 changes the status of the structural half, but only
to the local apex-shell statement.  At a tight opposite-binder point
`t=a/D`, it proves `14|D` and `D|(u+v)`.  The denom-14 location still needs
shell collapse `D=14h => h=1`; what remains is proving covering rows are
strict on every shell, or state-lifting shell height `h>1` to the forbidden
packet.

The exact audit makes the distinction visible:

```text
AP {1..13}:              M=1/14, non-covering, denom-14 witness.
GW {1..11,13,24}:        M=1/14, non-covering, denom-14 witness.
{1..11,13,84}:           covering, denom-14 blocked, M=7/89 > 1/14 at 37/89.
{1..11,13,168}:          covering, denom-14 blocked, M=14/173 > 1/14 at 72/173.
{1..11,13,504}:          covering, denom-14 blocked, M=42/509 > 1/14 at 212/509.
```

Before THM-568, the theorem target had to be stated as STAR+:

```text
For every primitive reduced 13-speed atom, M(S) >= 1/14,
with equality only on AP/Goddyn-Wong, and equality witnessed at denom 14.
```

Or as STAR0 plus a boundary-forcing lemma:

```text
M(S)=1/14 implies apex AP/GW,
and every hypothetical M(S)<1/14 covering family has a reduced boundary atom
with M=1/14.
```

Without that second line, equality-locus classification alone did not rule out
a non-apex sub-tight covering row.

After the corrected THM-568/HYP-2929 read and incoming THM-571, the preferred
residual is sharper:

```text
S = R union Q, Q subset 14*Z, |Q|<=6, R 14-free
  => M(S) > 1/14.
```

The 14-free core has a `1/13` margin by smaller LRC.  The multiples of 14 must
fail to cover the core's margin interval on every shell `D=14h`.  The
`|Q|>=7` apex-majority branch is now closed by THM-571 via `14 -> 7` gamma
descent.  The `|Q|<=6` branch is where the S31aa/S31v comb-teeth union bound
still needs a uniform finite-core / scale-separation input.  KPS-S31ab is useful
evidence for the AP/GW-to-multiple-of-14 perturbation family, but not yet a
general proof of this residual.

There are now three clean attacks.

1. THM-568 residual route: prove the `<=6` multiples-of-14 branch above.
2. Tight-locus route: prove no reduced non-AP/GW atom has `M<=1/14`, extending
   THM-560 and the Goddyn-Wong census into a full finite Node-2 theorem.
3. Forbidden-H7 route: prove the reduced covering atom trying to saturate or
   beat the apex bound yields a tournament-conflict-realizable activity-2
   packet graph with connected `I(.,2)=7`; HYP-2908 then turns it into the
   forbidden `K3` atom.

This is the sharper lesson from THM-079: the final move is not the numerical
identity `14=2*7`, but a category-correct atom-forcing theorem.  THM-568 now
supplies the apex-denominator forcing; the last open piece is the
multiples-of-14 covering strictness or an equivalent conflict-packet lift.

Artifacts:

- `05-knowledge/hypotheses/HYP-2909-lrc14-thm079-template-star-obligation.md`
- `04-computation/lrc14_thm079_star_audit_codex_s119.py`
- `05-knowledge/results/lrc14_thm079_star_audit_codex_s119.out`
