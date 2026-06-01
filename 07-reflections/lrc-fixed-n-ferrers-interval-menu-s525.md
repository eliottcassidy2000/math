# Fixed n means a tiny tournament menu (S525)

The user's phrasing is exactly the clean version of the LRC tournament
program: for a particular `n`, the problem is not "search the whole torus."
It is "prove the arithmetic clock hits a particular finite menu of tournament
isomorphism classes."

The size shift is:

```text
total LRC denominator n  ->  observer-marked tournament on n vertices
lonely/source state      ->  observer is a source
delete the observer      ->  runner tournament on m=n-1 vertices
```

So the target base is `A000568(n-1)`, not `A000568(n)`.  But the target is much
smaller than that base.  In an open source state all `m` moving runners sit in
the safe arc `[1/n,1-1/n]` of length `1-2/n`.  Once they are ordered in that
arc, a pair reverses relative to the order exactly when its separation exceeds
`1/2`.  Those reversed pairs form a Ferrers filter in the triangular pair grid.

That is the fixed-`n` structure:

```text
A000568(n-1)
  contains the Catalan-many Ferrers half-turn patterns
    cut down by the safe-arc length 1-2/n
      quotient by isomorphism
        plus compactified boundary/tie classes
```

The S525 computation separates the open and boundary pieces.  The open Ferrers
menus have isomorphism-class counts

```text
n=4..9: 1, 2, 4, 6, 10, 16
```

This is the fixed-safe-arc refinement of HYP-1996's circular menu.  From
`n>=5` the open counts agree with the circular `2*Fib` pattern; at `n=4`, the
safe arc has length exactly `1/2`, so the cyclic class is pushed to the wall
rather than appearing as an open tournament class.

Paired with HYP-1997 and HYP-1998, this gives the moving pieces of the
fixed-`n` statement: the LRC trajectory is an arithmetic-constrained walk on a
pointed metagraph, the open movie lives in the round/A000016 body, and the
target is this Ferrers interval menu plus its compactified wall boundary.

inside

```text
A000568(n-1): 2, 4, 12, 56, 456, 6880.
```

Bounded speed clocks hit exactly those open classes for `n=4..8`.  The older
S512/S520 menus are sometimes larger because they include equality-wall
witnesses.  That is not a contradiction; it is the THM-383 compactification
showing up in class language.

This makes the proof target sharper.  At fixed `n`, LRC asks for a forced visit
to a tiny Ferrers interval menu, or to its wall boundary, by every primitive
arithmetic clock.  The raw A000568 universe is the ambient quotient.  The
Ferrers interval menu is the real target.
