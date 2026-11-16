# A bridge cannot conceal the continuation drift just by changing notation

Let (H,X,D,g) be the retained-root state: X=E+O, D=E−O,
g the Hamilton count after deletion of the root. Follow the
[exact reverse-bridge law](BRIDGE_RESET.md). In the following
three experiments all probes retain this root. Normalize each
read by subtracting its own scalar multiple of the initial H.

| Future probe | Normalized count |
|--|--:|
| one right singleton | X+D |
| one right directed triangle | 3X+D |
| two consecutive right singletons | 3X+D−2g |

For a singleton the update is
H↦H+X+D, X↦2X−g, D↦−g, g↦g. For the triangle
the update is H↦3H+3X+D, X↦12X−3g, D↦−g, g↦3g.
The three normalized reads recover X from the first two,
then D, and g from the second and third. The row determinant
has absolute value 4, so this separation works modulo
every odd modulus as well as over the integers. In
particular a continuation quotient that preserves the
normalized scalar completions of these very short probes
cannot forget g. The usual root-change drift is therefore
observable, not merely bookkeeping in a padded state.
This does not refute a language which excludes these
probes or has a different boundary contract.

Even if the language uses only repeated triangles, the
first three normalized scalar rows are

    (3,1,0),  (45,3,−10),  (567,9,−168),

with determinant 648. Thus the full continuation is again
recovered over the rationals, although not uniformly at
3-adic precision through those three tests. This is the
finite observer idea from the planar Jacobian archive:
the relevant quotient is defined by the outputs of
actual legal continuations, not chosen because it
forgets a troublesome coordinate.

## Mixed policies: exact synchronization pins two ratios at once

For one right choice R retaining the old root and another
Y moving to the new root, let h_R=H(R),
h_-Y=H(Y−y), a_R,a_Y their total weighted start counts.
Assume Y−y nonempty. Equal g' and X' require

    h_R g = h_-Y H,
    (P_R(2)−P_(Y−y)(2))X = a_R g + a_Y H.

All counts on the right are positive. Thus two fixed
companions fix H/g and X/H whenever they synchronize.
A proposed reusable mixed-policy palette cannot tune either
ratio independently by just increasing the baseline.
Any surviving bounded difference must exploit D with
the remaining start/end constraints, or change the
companions with the target. The [actual count-three
case](BRIDGE_RESET.md) has synchronized +4, but no +2.

One small exact warning about weakening the ready state:
start with two scalar-only matched states with a positive
H difference, and do a common move to the right root.
If that root has a nonzero start-even count, then at
the new root the E difference is strictly positive.
The next conventional reverse bridge necessarily adds
positive mixed weight to the difference, because its
right end-odd response includes the singleton root and
is positive. Thus the raw product of the local H scales
no longer describes even the next scalar difference.
A triangle at its own root always has this behavior.
