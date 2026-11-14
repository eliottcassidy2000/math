# PRIMEGAME's hidden copy, return guards, and exact timing

**Proved:** the complete markerless-entry map, sharp first-epoch poison
test, and exact step counts. This is not a claim that the usual start
is flawed. Conway's [original paper](https://gwern.net/doc/cs/computable/1987-conway.pdf),
§§1 and 5, already proves the canonical prime-generation theorem and
explains the divisor loop. Here we extend the entrance and timing
rather than promote that inherited theorem as new. S323 already
compiled the Chan–Pak bounded search into a more general FRACTRAN
orbit; this is a precise look at PRIMEGAME itself.

## 1. Copy instruction and the concealed natural-number section

Use Conway's exact fraction order, with labels A through N:

    17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23,
    77/19, 1/17, 11/13, 13/11, 15/2, 1/7, 55.

A state C(n,d)=5^n 7^d·13 is a checkpoint. The candidate is n,
not the magnitude of the displayed integer, and d is the current
descending test divisor.

Start with a **markerless** 2^u 3^v 5^w 7^t, nonnegative exponents.
The L loop (15/2) copies each factor 2 into *both* the 3 and 5
registers. It preserves the difference between those two valuations.
Then M=1/7 drains t, N=55 introduces one factor 5 and a marker 11,
each EF pair moves one 3 to 7, and K swaps marker 11 for 13. The
first checkpoint is exactly

    n=u+w+1,     d=u+v.

Thus d−(n−1)=v−w, and the entrance is canonical exactly when
the old 3- and 5-exponents agree. Initial state 2^(n−1) enters
C(n,n−1), a clean copy of the candidate n (e.g. 2^9 starts at
candidate 10). Naked integer 10=2·5 is different: it enters (3,1),
still correctly gets the prime 3 but has skipped trying divisor 2.

A cofactor coprime to every participating prime rides along inertly;
claims about literal pure powers assume no such cofactor is present.
We do not assert a classification for arbitrary inputs carrying
spurious simultaneous state markers.

## 2. Exact return for all divisor states, and where poison matters

If d=0, J K alternate forever between marker 13 and marker 11, never
changing the 5-register. This is a literal two-cycle outside the
canonical orbit. Otherwise write n=qd+r with 0≤r<d. A full successful
subtraction of d copies from the 5-counter follows (AB)^d J (EF)^d K:
each AB moves 5 and 7 to 2 and 3, J moves 13 to 11, EF restores the
7s from the 3s, and K returns the marker 13. There are q such full
passes. They do not revisit the clean section because 2^qd is now
present (when q>0).

At the first incomplete subtraction the machine does (AB)^r A and
is at

    2^n 3^r 7^(d−r−1)·17.

If r>0, C consumes one 3 and changes 17 to 19; (DG)^n transports
all 2s back to 5; H adds a 7 and moves to marker 11; (EF)^(r−1)K
restores the 7-counter to d−1 and completes the failed trial.
The exact cost of the trial is

    (4d+2)q + 4r+2n+2 = 6n+2q+2.

If r=0, I removes the 17 marker to give 2^n 7^(d−1).
L^n M^(d−1) N (EF)^n K completes the restart at C(n+1,n).
The exact cost of this successful trial is

    (4d+2)q + 3n+d+3 = 7n+2q+d+3.

These formulas also cover oversize divisor windows d≥n: at d>n
one gets q=0 and declines until d=n, then succeeds there. For
n≥2 and d>0, the first candidate emits a pure 2^n precisely when
d<spf(n), the least prime factor of n. Above that threshold the
descending window includes a divisor other than one (including n
itself if d≥n). For n=1, d>0 always eventually emits the old
sentinel 2. In every positive-d case, after this first candidate the
next entrance is canonical and the normal candidate sequence persists.

Thus for composite n≥4, an entrance defect is falsely read as a
prime only if 1≤d<spf(n). The smallest such raw markerless seed
is 50=2·5², which emits 16. For n=4, d must be 1 and the minimal
2^u 3^v 5^w is (u,v,w)=(1,0,2); the only other choice is (0,1,3).
For n≥6 composite, d<spf(n)≤sqrt(n), so w=n−1−u≥n−sqrt(n)>3,
giving a larger factor of 5 alone. The ignored 7-factor cannot
decrease the minimal seed. This is scoped to the markerless lane.

The [paired-reset variant](https://research.vu.nl/ws/portalfiles/portal/2705740/243603.pdf) uses 15/14,15/2,55
in place of 15/2,1/7,55. In the markerless entrance its n is
unchanged, but d becomes v+max(u,t), because 7 is removed only
while consuming a factor 2. A canonical successful trial has
t=d−1≤n=u, so it agrees on canonical outputs and checkpoints,
saving d−1 steps each accepted trial. On other seeds the programs
differ: naked 7 reaches a d=0 two-cycle in the original, but d=1 in
the paired variant.

## 3. The hidden program clock has an arithmetic signature

At the canonical entrance C(n,n−1), put p=spf(n), g=n/p.
The exact time to the next canonical entrance is

    C(n)=Σ_(d=g+1,...,n−1) [6n+2 floor(n/d)+2]
           + 7n + 2n/g + g + 3.

This is Θ(n²) uniformly for n≥2: n−1−g failed probes, with
g≤n/2, each costing order n, plus the one acceptance. More
precisely,

    C(n)/n² = 6(1−1/p)+O(log n/n).

Thus the full accumulation set of normalised candidate costs is
{6(1−1/p): p prime} ∪ {6}. Each finite p is attained along p^k,
and increasing primes attain 6. Any other convergent subsequence
either stabilises some least prime factor on a subsequence or
has least prime factors tending to infinity. Canonical candidate
10 takes 338 steps and emits no 2^10.

There is also a cumulative constant. The integers with least
prime factor p form a periodic set of density
δ_p=(1/p) product_(l<p prime)(1−1/l), and their weighted second
moment up to N is δ_p N³/3+O_p(N²). The large-p tail in the
weighted sum for 1/spf(n) is bounded by 1/P after normalisation
for any threshold P, so finite periodic averaging and then P→∞
give

    lim_(N→∞) N^(−3) Σ_(n=2)^N C(n)
      = 2[1 − Σ_p (1/p²) product_(l<p prime)(1−1/l)].

This is timing from the hidden divisor layer, not information
visible in the sequence of prime readouts alone.
