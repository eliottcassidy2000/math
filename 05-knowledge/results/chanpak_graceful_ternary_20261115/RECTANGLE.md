# Sharp span of four distinct counts from two one-apex interval choices

**Theorem.** In the one-apex transitive-spine family, take two independent
two-setting controls on disjoint contiguous intervals, keeping all other
bits fixed. If the four Hamilton-path counts are distinct, their span is
at least **8**. In particular they cannot be four consecutive odd counts,
even if mixed interaction is allowed. The bound is attained.

This strengthens the additive obstruction in the graceful matching
decoder. The theorem does not require either control to occupy the whole
word, nor does it restrict the fixed prefix, middle block, or suffix.
It does not cover interleaved controls.

## Finite reduction, without a bound on the spine

Split the full word anywhere between the two changed intervals. Its
count is

    J(a,b) = h(a)+h(b)−1 + ab + (a mod 2)(b mod 2),

with ordinary binary left address a and reverse-complement binary right
address b. Here h(t) is the one-apex count of a binary presentation of
t, and leading zero padding does not matter. Relabel settings so the
left addresses are a,a+A and the right addresses b,b+B with A,B≥1;
write the four readouts p,q,r,s in rectangle order, with p=J(a,b)
and s=J(a+A,b+B).

If the four readouts are consecutive odds, their mixed difference
M=p+s−q−r lies in {0,±4,±8}. But J gives M=AB plus the product of
the two parity differences, hence M≥AB−1≥0. M=0 is precisely the
adjacent-address opposite-parity case. The graceful decoder has an
elementary proof ruling out four consecutive odds in that additive
case; the finite check below also covers it directly. In all cases
under consideration AB≤9, so A,B≤9.

The one-apex formula shows

    h(t+1)−h(t) ≥ −2 ceil(t/4)    for every t≥0.

At even t this is exact: flipping the final bit deletes the cycles
ending at that site. At odd t a carry creates a new 1 followed by a
run of zeros. The new 1 contributes positively; all existing 1s to
its left have nonnegative net gain because the new zero run's
weighted sum is at least the lost single zero's contribution.
Thus the odd increment is positive (explicit coefficients appear in
the graceful decoder). Summing the bound gives

    q−p ≥ A[b−(a+A+2)/2] − 1.

Since four consecutive odds have q−p≤6, it follows that
2b−a≤A+2+14/A≤17, the last inequality using 1≤A≤9.
Symmetrically 2a−b≤17. Together these force 0≤a,b≤17.

This leaves a finite exact universal check, *not* an enumeration
assuming that the two intervals or the entire spine are small.
The adjacent verifier checks all a,b≤17 and all A,B≤9 with AB≤9,
including M=0. In every case the four counts either repeat or have
span larger than 6. Since they are odd, the minimum possible
four-distinct span is then at least 8.

**Sharp witness.** Use a,a+A=0,1 and b,b+B=5,6. The four exact
counts are (3,9,5,11), with span 8. One literal choice uses left words
0,1 and right words 010,100. The verifier separately counts
Hamiltonian paths on these four five-vertex tournaments directly,
not just from J.
