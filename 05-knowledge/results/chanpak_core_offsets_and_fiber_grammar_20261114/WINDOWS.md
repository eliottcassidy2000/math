# Exact windows by changing the pair-path core

**Proved with literal certificates:** for every integer t>=1, there are
strong tournaments for every one of the counts

    288 * 4^t + {-11,-9,-7,-5,-3,-1,1,3,5},
    13790 * 2^t + {885,887,889,891,893,895,897,899,901,903,905}.

That is nine and eleven consecutive odd totals, respectively, with a
physical core selected from the residual. Neither statement says that
an arbitrary target lies inside one of these windows.

## Factoring the height out of the offset

For notation and the direct path-count argument see the inherited
[bilateral formula](../chanpak_movable_exact_windows_20261114/MOVABLE_TRIPLES.md).
A marked core A has an apex contact set U. Its counts are h=H(A);
f=ordered complementary path pairs (empty paths allowed); P and Q count
such pairs with a specified nonempty path starting in U and ending outside U,
respectively (for a single restriction either choice of the indexed path
is equivalent by exchange). S counts whole-core paths ending in U;
T counts those starting outside U. E counts pairs with first path ending
and second starting in U, both nonempty. The core plus apex has K=h+P-S-E.

Attach ordered transitive wings before and after A, then add the apex.
Its left contact word is the binary representation of x; its right
contact word is the reverse complement of the binary representation of y.
The core contact is U. The inherited count is

    fxy + Qx + Py + h(g(x)+g(y)-2) + K + T p(x) + S p(y),

p(x)=x mod 2, g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+p(v).
All these outputs are strong: the apex beats the first vertex and is
beaten by the last, every earlier vertex reaches the last, and the first
reaches all other old vertices.

Take u=2^t with t>=1. Iterating the trailing bits gives

    g(au)=g(a)+(u-1)a+p(a),       g(bu-1)=g(b-1).

For x=au and y=bu-1, the exact shape is D u²+J u+delta with

    D=fab; J=bP+a(Q-f+h);
    delta=h-E+h(g(a)+g(b-1)-a+p(a)-2).

The quadratic bank allows h, f and wing factors to change while holding
(D,J)=(288,0). This is distinct from a moving chart where h is frozen.

For x=au and y=b (positive fixed right address), the shape is D u+delta:

    D=a(fb+Q+h);
    delta=h+P(b+1)+S(p(b)-1)-E+h(g(a)+g(b)-a+p(a)-2).

The bank holds D=13790, but uses nine values of h and several different
actual f values; it is not a menu of scalar offsets chosen after the fact.
Literal rows, contacts and multipliers are in [window.py](window.py). For
example delta=895 uses rows (2,12,17,37,11,87,31), contact 27 and
(a,b)=(5,3). Its h is 49, not the final count nor the correction alone.

## Absolute decoding, and what is still missing

Given N, for each listed delta try (N-delta)/D and require exactly a
positive power of 2 (power of 4 for the quadratic family), with exponent
respecting t>=1. A success returns an actual tournament at exactly N,
with O(log N) vertices. The first windows are 1141..1157 and 28465..28485
(odd values only). On a specified height the eleven totals meet all odd
residues modulo 22, within distance 20 of the left endpoint; equality to
the target does not follow from residue matching. The union of the two
charts has only O(log X) distinct counts through X. The remaining selector
must let the baseline itself depend on N, with a residual menu that
overlaps the next baseline.

Replay with python3 window.py. It computes local profiles from subset
path DP and assembles enlarged tournaments for t=1,2, independently
recounting H and checking strong connectivity. Run python3 window.py 28475
to select an exact target. Outgoing rows use bits from vertex 0 upwards.
The all-height result depends on the certificates and identities, not
on completeness of the bounded small-core search that found them.

## Why filling between the dyadic heights needs another bank

The exact bit invariant rho(x)=g(x)-x+p(x) obeys rho(2x)=rho(x).
For odd v, g(v)=g((v-1)/2) and g(w)<=w+1 by induction, so
rho(v)<=(3-v)/2. For the even-height attempt x=a u with arbitrary
even u (right address still b), a row's actual value is

    13790 u + delta + h (rho(a u)-rho(a)).

The rows delta=885 and delta=889 have h=67 and 17 respectively, and
their a=10 and a=5 have the same odd part 5, with rho(5)=-1. If v is
the odd part of u, the first row minus the second is exactly

    -4 + 50 (rho(5v)+1).

For v>=3 this is at most -254, greater in absolute value than the
diameter 20 of the eleven-point window. Thus using the same eleven
rows with an arbitrary even height does not retain the window,
even if its absolute center is allowed to shift. This only rejects
that direct interpolation; varying the core bank or the right
address is still available.
