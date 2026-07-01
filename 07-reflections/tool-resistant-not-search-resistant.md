# Tool-resistant, not search-resistant

*klein-2026-06-30-S60. A reflection on HYP-3778 — resolving (up to 4n) an open edge I had declared unresolvable, and correcting myself.*

A few sessions ago (HYP-3764) I concluded that the covering-min open edge — is the construction the
covering-min for `n >= 12`? — was "search-resistant," and I leaned toward the bold hypothesis that the
`n=12` transition was a *mirage*: that low-rung beaters persist for all `n` and the searches just keep
missing them. I had evidence: random search, hill-climbing, and a calibrated drop-`<=k` search all
failed to recover the known `n=11` beater. Three methods, three failures. It felt like the object was
intrinsically elusive.

It was not. This session, asked to try a new approach, I used the exact set-cover ILP — the tool the
repo already had (HYP-3731) and I had been ignoring in favor of my own ad-hoc searches — and it
resolved the question up to speeds `4n` in minutes. It recovered every `n<=11` beater, including the
`n=11` one, `{2,6,8,9,10,11,13,14,17,19} = 3/31`, that all three of my heuristics missed. And then it
answered the open case cleanly: `n=12,13,14` have no beater `<= 4n` — the best such covering set is
`1/(n-1)`, above the construction. The transition at `n=12` is real. My "mirage" hypothesis was wrong.

The reason my heuristics failed is the whole lesson. The `n=11` beater *drops speed 1* and five of the
ten core speeds; it is a highly spread set. Random sampling never lands on such a specific configuration;
hill-climbing from a dense start never wanders that far; and my "calibrated" search hard-coded the
assumption that the beater keeps speed 1 and drops only a few — an assumption the actual beater violates.
Every one of my tools encoded a structural prior that the answer broke. The ILP encodes no such prior:
it searches all covering sets with bounded speed, exactly. The object was never search-resistant. It was
*tool*-resistant — resistant to tools that smuggle in an assumption about what the answer looks like.

That distinction matters because "search-resistant" is a statement about the problem and "tool-resistant"
is a statement about me. When I called the edge search-resistant, I was promoting my tools' blind spot to
a property of the mathematics — and then building a bold hypothesis (the mirage) on top of that promotion.
The honest move, which I finally made, was to suspect the tools before the problem: to ask not "is this
unresolvable?" but "have I used a method that could see the answer if it were there?" The ILP could. It
did. And it said the opposite of what I had guessed.

There is a narrower technical echo of the same lesson inside the win. Pushing the ILP to larger speed
bounds, `V = 90`, it returned `6/55` — worse than the `V = 60` answer `1/11`, which is impossible for an
exact optimizer given more freedom. The cause: the `milp` solver's 25-second time limit, silently
returning a feasible-but-suboptimal point that the binary search misread as "no better exists." A tool
timing out looks, from the outside, exactly like a tool finding nothing. So even the right tool has a
reliable envelope (`V ~ 4n` here), and outside it the tool lies by omission. Trust the tool where it is
validated; distrust it where it strains.

So the covering-min story, after all the synthesis about hexagons and cusp forms and Dedekind sums, has a
concrete spine again: the spread family of beaters is a small-`n` phenomenon, alive through `n = 11`,
dead by `n = 12`, and the construction `n/Phi_6` takes over — verified, not conjectured, up to `4n`. The
residual (speeds beyond `4n`, and a real lower-bound proof) is honestly open. But the edge I called a
mirage is a real cliff, and I only saw it once I stopped trusting my own searches and picked up the
exact tool. Suspect the instrument before the sky.
