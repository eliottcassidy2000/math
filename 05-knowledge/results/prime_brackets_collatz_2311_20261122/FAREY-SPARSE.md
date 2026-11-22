# Any-depth signed Collatz exclusion with at most two odd steps

**FINITE-EXACT.** At spine 4 and at spine 8, the two positive E61
Farey flanks cannot be related in either orientation by a shortened
Collatz script of positive length containing at most two odd steps.
At each odd step the sign in \((3x\pm1)/2\) is independent. This
extends the [bounded-depth fixed-sign search](../chanpak_farey_two_three_clocks_20261122/COLLATZ-MAHLER.md) in a
different direction. The unrestricted Farey search remains open.

For \(d\) total steps and \(m\) odd inputs at positions \(t_j\)
(starting at zero), let
\[
 C=\sum_{j=0}^{m-1}\epsilon_j2^{t_j}3^{m-1-j}.
\]
The exact relation is \(2^d y=3^m x+C\). The congruence at \(2^d\)
is necessary and sufficient for the chosen parity word: reducing
modulo 2 forces its first input parity, and stripping the first
step repeats the same argument. Thus \(x=B+2^d z,\ 
y=(3^mB+C)/2^d+3^m z\) with
\(B=-C(3^m)^{-1}\bmod2^d\). Require both addresses positive.
Along this affine ray the last solution to the quadratic inequality
\(P(x,y)\le N\) is the only address pair worth testing, since
\(H(x,y)<P(x+1,y+1)\), as in the previous search.

The height and positivity of the opposite flank imply
\(x\le\lfloor(N-k-c)/(q+t)\rfloor\) in physical forward orientation;
transpose \(t,c\) for the other orientation. An odd step increases
the state by a factor of at most 2 and every other step halves it.
So the positive endpoint implies \(d\le 2m+\lfloor\log_2 x_{\max}\rfloor\).
Both orientations then have total depth at most 103 or 97
respectively. This is a proved absolute cutoff, not a depth cap
set by runtime.

| Spine | Signed presentations of positive length, both orientations | Parity-corrected envelope | Surviving modulo \(2h\) | Exact hits |
|---:|---:|---:|---:|---:|
|4|1,477,634|3|0|0|
|8|1,235,194|1|0|0|

Run `python3 farey-sparse-signed.py` (standard library). The replay
recovers the mixed-sign planted chain \(29\to44\to22\to11\to16\to8\)
at its own chart height in both orientations and both spines.
