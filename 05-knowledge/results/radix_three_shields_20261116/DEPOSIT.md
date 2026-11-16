# A repeatable independent deposit behind an exact ×3 step

The [corrected fixed-root bridge law](BRIDGE_BOUNDS.md) supplies
a genuine two-choice deposit. It multiplies every old count
difference by 3 and deposits an independent new bit without
changing any coordinate used by the next bridge. Its deposit
is much too large to be a ternary odometer.

Use right modules V0=C3→TT2, rooted in the triangle, and
V1=TT2→C3, rooted at the second singleton. Both have ordinary
count 3, two-slot count 48, and weighted root-start profile
(E_s,O_s)=(6,6). Their weighted root-end profiles are
(1,2) and (3,3), respectively. These values follow by reading
the consecutive modules into the optional before/after slot,
or by literal weighted path enumeration. For example a root
in V0 reaches the right transitive pair with both parity sums
2; the root-ended paths stay in the triangle. In V1 the right
start reads the complete two-slot triangle with both parity
sums 6; the only possible predecessor run before its root is
the previous singleton, and the untouched triangle contributes
3 to each of the two root-end parities.

On the left state (H,E,O,g), either choice has the identical
continuation

    E' = O' = 24(E+O)−6g,     g'=3g.

But its output count is

    V0: H'=3H+4E+2O,
    V1: H'=3H+6E+6O.

This is a literal triangular transfer, not a product of
independent local counts. It is repeatable for any number
of choices, since retaining the original root uses exactly
the same state coordinates. All old count differences are
scaled by 3, while the new right choice contributes
2E+4O at that layer. We make no strongness claim for these
right modules, which have nontrivial ordered components.

For a concrete start, take a rooted C3 with
(H,E,O,g)=(3,1,2,1). The first two outputs are 17 and 27,
with shared state (E',O',g')=(66,66,3). At the next step
the four outputs are 447,477,843,873: both old settings
survive independently of the right choice. After the first
step write E=O=W. At later depths

    W'=48W−6g,   g'=3g,   H'=3H+(6+6b)W,   b=0,1.

Hence b contributes 6W at the moment of insertion, then
accumulates a factor 3 for each later layer. At a fixed
output depth the ratio between successive (newer over older)
weights is (48−6g/W)/3, greater than 15 in this example,
except for the first ratio 396/30>13. This gives a
superincreasing binary scale, so the gaps do not shrink.
It is a positive independent **deposit** mechanism, but
a decidedly non-ternary one. To realize width-four radix
digits, one needs more than an exact transport of differences:
the new bounded deposit must also be calibrated independently
of the growing endpoint response.

The [replay](check_bridge_bounds.py) checks the exact
profiles, the first 17/27 switch, and all four next outputs.

There is also a general parity check on trying to force the
endpoint growth ratio itself to be 3: for every nonempty V,
P_V(2) is even (swap the two ordered slots, a fixed-point-free
involution), while H(V) is odd. Thus P_V(2) cannot equal 3H(V).
A proportional scale copy alone cannot supply an odd ratio
three; any successful correction would have to use another
coordinate or a non-proportional return.

### All count-three right profiles: why this deposit cannot be tuned small

For fixed P_V(2), a right count-three module is a triangle
plus a fixed number of singleton SCCs. A root lies in the
triangle, before it, or after it. Let u,v be the singleton
counts before and after the root's component, excluding
the root itself if it is a singleton. The exact total start
and end readouts, and their signed even-minus-odd differences,
are as follows:

| Root location | Start total | Start sign difference | End total | End sign difference |
|--|--:|--:|--:|--:|
| triangle | 3·2^v | −1 if v=0, else 0 | 3·2^u | −1 if u=0, else 0 |
| before triangle | 12·2^v | 0 | 3·2^u | −3 if u=0, else 0 |
| after triangle | 3·2^v | −3 if v=0, else 0 | 12·2^u | 0 |

These follow from factoring along the ordered components;
a subsequent or preceding nonempty optional component cancels
a signed parity sum. At fixed h_V=3, two-slot count and
start profile, each root type has at most one end profile.
If the singleton-before and singleton-after types match,
their end profiles coincide (their u values differ by two,
and both end parity differences are zero). If a triangle
root also matches, its end total is exactly half of the
singleton-root end total; the latter has equal parities.
Thus there are at most **two** distinct end profiles for
any synchronized right count-three palette, and any
nonzero change increases both end coordinates in one direction.
The left deposit 2(E ΔO_V+O ΔE_V) cannot be a cancellation
against the fast mode. This rules out a four-state bounded
deposit in this fixed-root H=3 language, not a different
root or a shield read by different coordinates.
