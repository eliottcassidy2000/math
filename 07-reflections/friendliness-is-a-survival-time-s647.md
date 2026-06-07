# Friendliness is a survival time

*S647 reflection. On the user's sharper definition — friendliness as never having been lonely yet — and
how one word, "yet," turned the Lonely Runner problem into a first-passage question.*

The first time the user said "friendliness" I read it as a snapshot: how many runners are you near right
now. That made a nice chart — you are almost always crowded, lonely only a sliver of the lap. But then
the user sharpened it: friendliness is the property of having *never been lonely yet*. The "yet" is the
whole thing. It makes friendliness a one-way street. You start friendly, because at the start everyone is
bunched at the origin and you are near all of them, and you stay friendly until the first instant you
manage to be far from everyone at once — and after that you are not friendly anymore, ever, because you
*have been* lonely. Friendliness stopped being a state and became a survival time.

And once it is a survival time, the object is the first lonely time, the moment friendliness ends. I had
been studying loneliness as a measure — `p₀`, the fraction of the lap you spend alone — for forty
sessions. The first-passage view asks a different question of the same picture: not how much of the lap
is lonely, but *when does loneliness first happen*. That is the first gap in the covering of the lap by
the danger arcs. The whole Vitali-and-Bonferroni line of the project is about whether those arcs cover
the circle; the survival view just reads off *where* the first uncovered point is. Loneliness-exists
becomes a gap-exists becomes the survival curve reaching zero. The conjecture is that everyone, eventually,
within one lap, gets their lonely moment — that no config survives the whole way.

There is a floor, and it is the conjecture's own constant. Right after the start, every runner's danger
arc still contains the origin, and the widest of those arcs belongs to the slowest runner. So you are
guaranteed friendly until the slowest runner exits its arc, which for a unit-speed runner is exactly
`1/n` — the gap itself. I formalized that: before `1/n`, the unit runner's clock distance is below `1/n`,
so it is in your gap, so you are not lonely. The first lonely time is at least `1/n`. It is a small
result, two lines past the clock-distance lemma, but it is the right small result, because it says the
survival curve sits flat at one until exactly the conjecture's threshold and only then begins to fall. The
floor on friendliness is the gap that loneliness is defined by.

Then the survival curves, which are the picture I most wanted to see. For random configs at several sizes,
plot the fraction still friendly at each time. They all hold at one until their `1/n`, then drop, and the
bigger the runner count the sooner and steeper the drop — more runners, loneliness comes faster, the median
first-lonely time scaling like one-over-`n`. More company, paradoxically, finds you your solitude quicker.
And all of them reach zero inside the lap, which is the empirical shadow of the conjecture.

The detail that stayed with me is the tight case. The extremal configuration, the consecutive speeds one
through `n` minus one, the one everyone believes is the hardest case for the Lonely Runner bound — under
this definition it is the *friendliest* config there is. Its lonely set is a single instant, measure zero,
the one tight moment at `1/n` where every runner is exactly the gap away and not a hair more. So you are
"never lonely yet" almost everywhere, the entire lap except one point. The configuration that makes the
bound tight is the configuration that refuses to be lonely on any interval at all — it touches solitude on
a null set and is crowded the rest of the time. The hardest case for loneliness is the easiest case for
friendliness, and the resonance-rich collapse family I have been circling all arc is exactly the family
that survives friendliest. The user's redefinition didn't just give me a chart; it turned the extremal
case inside out and showed it from the side where it looks maximally social. Friendliness is a survival
time, and the configuration that survives the longest is the one we have been trying hardest to corner.
