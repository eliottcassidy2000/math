# Complete physical offset inventory at slope 908 and core count nine

**Finite-exact after an all-order reduction.** Consider the unrestricted
left-address chart with one common apex over a marked inner core, a
left transitive address u>=1 and a fixed right address b>=1. Allow
arbitrary fixed low bits of the original left address: these belong
to the inner core once they are absorbed, as in [LINEAR_WINDOWS](LINEAR_WINDOWS.md).
Keep h=9 and a linear coefficient D=908 in front of u, before
adding 9g(u). Then the list of absolute responses is precisely

    908u+9g(u)+{624,706,752,754,756,758,772,784}.

There is no parity sidecar in any admitted marked core and right
address. Thus the four-point window already found is maximal within
this precise chart, regardless of original core size or its fixed
left padding.

The [bilateral readout](../chanpak_movable_exact_windows_20261114/MOVABLE_TRIPLES.md)
gives D=f b+Q and parity sidecar T p(u); the constant is
P b+9(g(b)-2)+K+S p(b). Since 0<=Q and f is a sum over all
subsets of products of positive path counts, f>=2^n. Thus n<=9
for any physical core (before adding the two free/fixed wings).
At h=9, ordered SCC factorization leaves only one nontrivial H9
component or two H3 components; arbitrary singleton SCCs are free.
The [Busch bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
bounds either H9 component to at most five vertices and H3 to
at most three. Through order four no class has H9. The two strong
order-five H9 classes are the known
[Session 224 birth fiber](../session224_first_occurrence_20261107/README.md);
their literal rows are included in the replay. At H3 the only
strong component is the directed triangle.

For each such ordered word and each contact mask, the exact pair
DP determines f,Q and hence the **only** possible b=(908-Q)/f.
Only seventeen marked choices survive the positive integrality test;
their responses collapse to the eight offsets above. They are physical,
not polynomial extrapolations. The population counts marked choices
in this enumerator; it is not a count of unpointed tournament classes.
Masks related by a block automorphism may both be counted. The
all-order size argument and exact class inventory establish completeness,
unlike an unrestricted brute-force tournament-size trial.

The same local-band proof in [COMPLEMENT_GAPS](COMPLEMENT_GAPS.md)
applies to all eight offsets: their span is 160, and for A>=256 the
lower gap at fixed low-tail class is greater than 2A+160. Thus this
entire slope and core-count slice supplies less than 80A/908 of the
A-1 odds in the indicated band near A². This rules out extending
coverage in that band by finding other hidden masks at these same
slopes, without ruling out a new slope, different h, both free
addresses, or a non-affine selector.

Run python3 nine_slope.py to replay the complete chart inventory.

## All slope banks through 4095: the apparent alternatives are aliases

The same exact cutoff with variable D<=4095 allows inner orders only
through 11. This gives 86,523 marked choices after including the right
address parameter, 6,286 distinct pairs (D,T), and no five-window.
The full four-window list, recording (D,T,first offset), is

    (908,0,752),   (1816,0,1660), (1825,9,752),
    (3632,0,3476), (3641,9,1660), (3650,0,2586), (3659,9,752).

Each of the six larger-slope entries is exactly the original window
with one or two fixed low bits of u. To verify, a chart written as
D u+9g(u)+T p(u)+delta turns into:
on u=2v+1, D'=2D, T'=0, delta'=delta+D+T;
on u=2v, D'=2D+9, T'=9, delta'=delta.
These rules give precisely the other six entries and do not add
any absolute targets to the widest-window output. This identifies
the minimal slope 908 in the sweep and rejects an apparent family
of independently movable same-width charts. It does not classify
slopes above 4095 or narrower windows in this range.

Replay the inclusive sweep with python3 nine_slope.py --sweep.
