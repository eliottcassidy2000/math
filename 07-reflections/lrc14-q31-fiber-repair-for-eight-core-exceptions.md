# LRC14 Q31 fiber repair for the eight-core exceptions

The first below-nine-core layer did not behave in the most convenient way, and
that is probably good news.

HYP-2470 is the clean theorem for the LRC proof route:

```text
retain at least eight core speeds in the carry window
=> Q27 witness or plain q<=41 witness.
```

This note records the companion fiber view.  The naive target after HYP-2465
was:

```text
delete four core speeds + add five window speeds
=> still cannot block Q27.
```

That statement is false.  Two deletion shapes really can block Q27:

```text
(28,42,56,84)
(42,56,70,84)
```

HYP-2470 repairs the false statement by adding the missing plain shells through
`41`.  HYP-2471 repairs the same false statement in a different language: keep
the divisor/fiber ladder `{1,2,7,14}*m`, but extend it from `m<=27` to `m<=31`.
Both exceptional deletion shapes become set-cover infeasible on Q31.

So the proof grammar from HYP-2469 got more concrete:

```text
scalar Q27 quotient passes in two finite shapes
+ retained plain shell side-channel
=> shell-41 theorem

scalar Q27 quotient passes in two finite shapes
+ retained ladder side-channel
=> Q31 fiber repair
```

The bounded statement is now two-layered:

```text
inside 1..1092, retaining >=8 core speeds forces Q27 or plain q<=41.
inside 1..1092, retaining >=8 core speeds also forces Q31.
```

The first line is the one to use when proving LRC14 directly.  The second line
is the one to use when trying to understand why the obstruction cannot survive
inside the carry/fiber grammar.

This changes the next proof target.  We should not say "below-nine-core" as if
the whole region is still live.  The bounded live region is now below-eight-core:
rows retaining at most seven speeds of `7*{1,...,12}`, plus anything that leaves
the carry window.

The two Q27 exception shapes are also telling us which typed budget matters.
They delete high even core speeds and cover heavily with `7`-multiple /
`13`-clock / `q=91`-fiber addresses.  That suggests the next `e=5` computation
should not be a blind larger MILP.  It should split the added speeds by role:

```text
primitive key,
7-multiple support mass,
13-clock debt,
q=91/q=161 divisor fiber,
shell-27 antipodal class,
Bprime owner.
```

The wild but plausible theorem shape is now:

> A bounded no-Q27/no-plain-41 row with eight or more retained core speeds is
> impossible.
> A bounded no-Q27/no-plain-41 row with seven or fewer retained core speeds must
> either have concentrated support load, a low-clock opening, or a deletion/owner
> descent that lowers the number of non-core packets.

Q31 is not a new magic number.  It is the first fibered place where the hidden
support that Q27 forgot becomes visible.
