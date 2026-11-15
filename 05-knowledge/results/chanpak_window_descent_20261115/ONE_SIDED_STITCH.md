# A five-count window by facing two one-sided palettes toward each other

A full two-sided switch at spine count 3 does not exist (see
[the CRT window note](CRT_WINDOWS.md)). That does not forbid using 3
as the *first* block: it need only keep its outward observables.
If there are just two varying blocks, the right block likewise need
only keep its observables on the inward side. This lets two larger
physical palettes face one another across an arbitrary common
transitive tape.

## Literal palettes

A word is a string of ordinal singletons marked 0 or 1; the one apex
beats precisely the marked vertices. Rows are outgoing bitmasks; a
core is specified by its rows and its mark mask. Each row of the table
represents word → core → word, with the apex added after the complete
spine is assembled. Empty words are written −.

| Core rows | Prefix | Mask | Suffix | Local apex count |
|---|---:|---:|---:|---:|
| (2,4,1) | 10 | 6 | 11 | 17 |
| (2,4,1) | − | 6 | 1011 | 23 |
| (2,4,1) | 1010 | 7 | − | 27 |
| (2,4,1) | 1 | 4 | 011 | 29 |
| (2,4,1) | 101 | 4 | 1 | 43 |
| (6,8,10,1) | 100 | 9 | 0 | 233 |
| (6,8,10,1) | − | 0 | 111000 | 255 |
| (6,8,10,1) | − | 8 | 11000 | 259 |
| (6,8,10,1) | − | 6 | 1000 | 261 |
| (6,8,10,1) | 10 | 2 | 00 | 267 |

In the first bank every right state is $(h,S,d)=(3,129,3)$.
The second bank as written has right state $(5,280,0)$; **reverse
all arcs and complement the mask** before using it on the right.
Its local apex count is unchanged, and its left state, the
unmarked-end total and signed bias, is $(h,t,f)=(5,280,0)$.
Literal subset DP verifies these identities in the replay.

## Stitch and invert

Between a chosen left and right member insert $k\ge1$ zero singletons.
The left state immediately before the right bank is
$(3,129\cdot 2^k,0)$. Its apex count has grown by
$129(2^k-1)+3$. The [cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md)
now gives exactly
\[
N=36765\cdot2^k-645+5a+3b,
\qquad
a\in\{17,23,27,29,43\},\quad b\in\{233,255,259,261,267\}.
\]
The offsets $5a+3b=910,912,914,916,918$ come respectively from
$(29,255),(27,259),(43,233),(23,267),(27,261)$. Thus every $k\ge1$
has the five-count window
\[
36765\cdot2^k+\{265,267,269,271,273\}.
\]

This example additionally gives **every odd residue modulo 30 at
every such fixed height**, not merely when the height may vary.
Indeed the five $a$ cover all three odd classes modulo 6 and the
five $b$ cover all five odd classes modulo 10. Reducing the formula
first modulo 3 and 5 chooses the two arms independently. After that
choice, test the integer quotient: that is a separate constraint.
For example the offset $5a+3b=920$ is in the interval between the
minimum and maximum and its class is present, but it is absent from
the bank itself.

The zero tape can in fact be replaced by **any** common transitive
marked tape, shifting all 25 outputs by the same amount. The first
block's right state stays matched through the tape; the last block
reads only that state and its own fixed incoming data. The completed
right state of the *whole* construction need not match, so this
cannot simply be nested as another independent internal digit.
That is the useful design lesson: broaden the palettes at the two
outer ends, and ask full matching only of the blocks genuinely in
the interior. To cover an arbitrary target one still has to arrange
the common translation and its carry, not just pass the two residue
tests.
