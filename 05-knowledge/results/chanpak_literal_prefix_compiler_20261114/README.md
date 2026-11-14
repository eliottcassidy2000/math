# Actual tournament prefixes at the prescribed coordinates

**Two previously missing kinds of physical source are now explicit.**
The construction still does not establish universal Chan–Pak coverage.

| prime cofactor \(p\) | real prefix state \((a,b)\) | prefix order | ordinary count after the seven seam |
|---:|---|---:|---:|
| 1,137,067 | \((1105,2634)\) | 15 | 7,959,469 |
| 1,512,281 | \((23,4736)\) | 12 | 10,585,967 |

In the first row the least arithmetic point \((469,3176)\) is
not a physical terrace state, but another point on the **same target
line** works. In the second row the chosen coordinates are
impossible for *any* terrace prefix. An arbitrary-root preparation
gets them from a different kind of six-vertex core. Both prefixes
were checked by literal path-subset dynamic programming; the exact
extension proof gives the full counts. The compiler outputs the complete
adjacency masks of both the prefix and the final strong tournament.

- [Arbitrary-root bridge](ARBITRARY_ROOTED_BRIDGE.md): the signed
  rooted subset-path profile becomes exactly ready after a positive
  bridge. Thereafter the same corridor matrices act literally, even
  when the old core was not a terrace. This unlocks the actual
  deletion pair \((H(C),2^r F_2(C))\) for arbitrary cores \(C\).
- [Physical boundary](PHYSICAL_BOUNDARY.md) and
  [target compiler](TARGET_COMPILER.md): at maximal outer depth
  test all at most 49 points on the target line, with exact slope and
  predecessor filters for the terrace. A separate strong-component
  factorization chooses a rooted nonterrace source from a real core
  bank and checks the actual pair, not a relaxed cone.
- [Guarded growth](GUARDED_GROWTH.md): an unbounded bank with distinct
  inner structures is selected back by the target's own maximal
  inverse; modular valuation peels its source directly.
- [Two-context switch](SOURCE_SWITCH.md): swap two small corridors
  rather than continue both states in the same way. A commutator
  gives exact cofactor adjustment \(\pm2a\). For \(a=1\) this is the
  desired unit adjustment, with explicit prime cofactors. The
  [combined switch](COMBINED_SWITCH.md) applies it to a literal new
  nonterrace source \((5,64)\), giving a prime cofactor 3,539,307,127
  next to a cofactor divisible by 3.
- [Changed equal-coordinate sources](CHANGED_EQUAL_U.md): four
  literal rank-13 pairs, including an inverse-derived pair at
  common value 26,235 and a reset-free pair that produces
  the ordinary prime 26,360,017 beside 26,360,019.

Run `python construct_target.py 1137067 1512281` for literal witnesses.
Run `python bridge.py`, `python physical_seed.py`,
`python physical_prefix.py`, `python replay_source_switch.py`,
`python replay_combined.py`, and `python replay_changed_equal_u.py`
for the independent checks.

**Remaining obligation.** No argument says that every prime cofactor's
line meets one of the source classes, or that every ordinary prime is
covered. The source search now handles whole target lines and arbitrary
rooted cores, but choosing a universally successful physical pair is
still the absolute-target problem.

**Follow-up: two contacts.** [The new endpoint correction](TWO_CONTACT_REPAIR.md)
gives actual two-width ordinary unit repairs, including adjacent
prime 2053 and 2051 = 7 times 293. It also constructs a third
exact target-line source: \((1349,2718)\) for prime cofactor 1,229,903,
a nonterrace state with small even defect rather than a forbidden
dyadic slope. Run `python two_contact.py` for the direct replay and
the target-line inverse for this expanded bank.
