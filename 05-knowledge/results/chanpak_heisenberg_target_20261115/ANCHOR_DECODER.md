# An absolute-target decoder for changing binary-cut hosts

**Proved, for a specified physical class.** The class varies the host at all orders and has two independently chosen transitive cores of positive orders. It is not the full Hamiltonian spectrum and there is no coverage claim for this class. One complete search in this class takes only O(sqrt(N)) candidate words, each processed in O(log N) simple integer operations. A narrower subclass gives an exact uninterrupted high-digit window and identifies the missing resource as low-digit checksum repair.

## 1. Physical host and its census

Take the one-outlet two-core gluing with exterior R and anchor q: R→W→U, U→q, R−q→U. Here R−q is a transitively ordered spine, and w is its bit word with bit 1 iff q beats that vertex. Let u=val(w) (first bit most significant), A(w)=H(R), and B(w) count the two ordered path slots covering R with first path starting at q. K=H(R−q)=1. Empty spine is allowed (R is a singleton). If the spine is nonempty, the result is strong exactly when w begins with 1: then q reaches the first spine vertex and hence every vertex of R. A word beginning with 0 cannot reach that first vertex; earlier initial zero bits could instead be treated as neutral ordinal prefixes for objectives that do not require strength.

Each anchored first path selects no spine vertices (one possibility) or has a first selected spine vertex carrying bit 1, with all later selections free. Thus B(w)=1+u. The odd-cycle packing expansion gives
```
A(w) = 1 + sum of f(j−i−1) over pairs i<j with w_i=1,w_j=0,
f(0)=2, f(d)=2^d for d≥1.
```
Every directed cycle uses q. Its spine vertices are listed in increasing order, starting in a 1 and ending in a 0; an odd cycle uses an even number of spine vertices. For given endpoints, the number of allowed interior selections is 1 at distance d=0, or 2^(d−1) otherwise, and the odd-cycle expansion gives each cycle weight 2. For each i, replacing all later bits by 0 bounds its contribution by 2^(n−i), where n is the spine length; hence 1≤A(w)≤u+1. In particular A is not an independently chosen coordinate even though B ranges over a contiguous binary band.

For fast replay, append bit b to a word whose last bit is t, with t=0 at an empty word. Then
```
u' = 2u+b
A' = A + (1−b)(u+t).
```
The terms are exactly the inversions into the new final 0, with the adjacent term counted twice.

For transitive U=TT_m, W=TT_s, set x=2^m−1, C=2^s x. The exact physical path count by the two-core gluing is
```
N(w;m,s) = A(w)+1+x(2^s(u+1)−1).
```

## 2. A product with just one boundary correction

Split w into prefix p and nonempty suffix v of length k. Put a=val(p), b=val(v), r=rev_k(b) (reverse the k bits, retaining zero digits). Then
```
A(pv) = A(p)+A(v)−1 + a(2^k−1−r) + e,
e = [p ends in 1 and v starts in 0].
```
The cross inversions factor: a prefix start contributes its position weight to a, and a suffix end contributes its reversed-position weight to 2^k−1−r. Every cross inversion has this product weight except the adjacent one, whose true contribution is 2 rather than 1. This correction e is essential; the empty prefix has e=0.

Equivalently,
```
N(pv;m,s) = N(v;m,s) + aD + A(p)−1+e,
D = (C+1)2^k−1−r.
```
As the suffix varies over its 2^k words, D visits every integer in [C·2^k,(C+1)·2^k−1] exactly once. This is an independently selectable slope before the low correction is evaluated.

### Exact inverse for the whole class

For specified m,s and total spine length n≥2, split k=ceil(n/2). Enumerate each fixed-width suffix v. The residual 0≤A(p)−1+e≤a+1≤2^k is strictly less than D≥C·2^k (C≥2). Therefore the prescribed target determines the only possible prefix:
```
a = floor((N−N(v;m,s))/D).
```
Reject if a does not fit the prefix width with first digit 1. Otherwise rebuild the prefix and check the exact residual. This is a complete inverse, not a modulus hit. The n=0,1 cases use their unique valid strong anchor words. Since N≥C·2^(n−1) for a nonempty strong anchor, summing 2^ceil(n/2) over feasible lengths for each m,s is O(sqrt(N/C)). Then sum over m,s using C≥2^(m+s−1) to get O(sqrt(N)) candidate suffixes. The n=0,1 loops contribute O(log² N) candidates and do not change the bound. Checking (rather than listing) each witness takes O(log N) simple integer operations, with O(log N)-bit words up to constant factors. This search deliberately produces only this class.

### Gap-separated exact high-digit window

Use the fixed prefix p=1 0^g for g≥1. Here A(p)−1=2^g and e=0. Thus
```
N(pv;m,s) = 2^g[(C+1)2^k−r] + [Cb + A(v)+1+C−x].
```
For g≥k+ceil(log₂(C+2)), the positive bracket is less than 2^g: by A(v)≤b+1 it is at most (C+1)2^k+1−x≤(C+1)2^k. Therefore the quotient floor(N/2^g) takes every value in [C·2^k+1,(C+1)·2^k] exactly once. Given that absolute quotient, choose the physical suffix by reversing r=(C+1)2^k−floor(N/2^g); the remaining necessary and sufficient condition is equality with the explicit low checksum. This really selects every leading window position, but not every target within those windows. Even complete freedom of the quotient has left a single low readout per position.

## 3. Two variables that can be changed separately

For an arbitrary tournament in the place of a transitive W, the same anchor (K=1) and transitive U with x=2^m−1 give
```
N = [A(w)+1−x] h(W) + x B(w) Q(W).
```
This entire slice, conditional on the second core W, still has the same target-driven inverse. Write t=h(W), D_W=xQ(W). Then
```
N(pv)=N(v)+a [D_W·2^k+t(2^k−1−r)]+t[A(p)−1+e].
```
The residual is at most t·2^k for a half-word split, and the bracketed slope is at least 2t·2^k since Q(W)≥2h(W). Thus again take the floor to recover the only admissible prefix for each suffix. Also N≥D_W·2^(n−1), so for fixed m and fixed W the enumeration takes O(sqrt(N/D_W)) candidates plus the two endpoint cases. This conditional decoder does not enumerate arbitrary W at that cost.

The gap-separated form is particularly transparent:
```
N(1 0^g v)=2^g [(D_W+t)2^k−t r]+N(v),
0<N(v)≤(D_W+t)2^k+(1−x)t.
```
For sufficiently large g, its quotient is an arithmetic progression of step t=h(W). If that step divides the target's difference from the top endpoint and the index lies in range, the high clock selects the literal host; exact equality is then the checksum test. Changing W changes both the step and the baseline. The open problem is to arrange translated high-quotient bands and their physical low checksum sets to satisfy every prescribed N, not just to solve a congruence within a fixed band.

## Replay

Run `python replay_anchor_decoder.py`; it checks 3,000 seeded seam/quotient instances and independently enumerates actual Hamiltonian paths in 120 literal glued tournaments.
