# Two square roots replace the inverse search inside an exact depth

The [square-sum ternary family](SQUARE_TREE.md) puts both addresses into
a fixed-sum interval at each depth. This supplement gives a direct
absolute decoder on that interval, without binary search.

Fix S=x+y and center j (determined by an odd target N mod 3), and restrict
1≤t=min(x,y)≤floor(S/3). The two orientations have quadratic baselines

    b_-(t) = −432t²+(432S+129)t+129S+135+2j  (x=S−t, y=t),
    b_+(t) = −432t²+(432S−129)t+258S+135+2j  (x=t, y=S−t).

The omitted term is between 0 and 45S+45. For steps t→t+1
within the interval, both baselines increase by at least
144S−129: use S−2t−1≥S/3 when t+1≤floor(S/3). For S≥3 this
exceeds 45S+45. Thus for each orientation the uncertainty bands
are disjoint; at **most one candidate t** need be checked.
The reversed-order comparison from the square tree note shows that
for its depth bands the two orientations together cannot collide
either. No monotonicity of the checksum is needed.

For b(t)=−432t²+B t+C, take h=floor(S/3). If N<b(1), reject;
if N≥b(h), test h; otherwise compute the exact candidate

    t = floor((B−ceil(sqrt(B²−1728(N−C))))/864).

This is the greatest t with b(t)≤N on the smaller-root branch.
The radicand is nonnegative because N<b(h) and h is on the
increasing branch. Replacing the radical by its integer ceiling
is exact after flooring, since the divisor and the other term are
integers. Test the literal checksum; then test the physical
square-sum root and branch addresses using their unique ceiling
inverse. This yields an exact, target-first decision on the
picked tree using at most two square roots per candidate depth
and no orientation scan. For square-sum seed (1,3), depth zero
is handled separately.

The method and the scaled-swap inverse are instances of the same
criterion: find a monotone leading height, prove that all remaining
graph data fit in disjoint height bands, then enforce the remaining
physical conditions exactly. A successful height decoder diagnoses
a missed target but does not manufacture a branch for it.
