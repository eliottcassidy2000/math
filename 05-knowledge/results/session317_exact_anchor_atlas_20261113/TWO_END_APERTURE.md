# An exact two-ended aperture giving adjacent absolute counts

Let R have outgoing masks (10,28,25,16,1), with H(R)=9. The marked options b=0 and b=2 both have the optional-second-slot end response E=11, but have start counts s_b=1 and 3, respectively. Form the ordinal sum of a singleton a, an arbitrary nonempty tournament M, and R. Add apex x, beating a and every R vertex except b, losing to all of M and to b. This output is strong: R has the directed cycle 0,1,2,3,4,0, all vertices of the ordinal join reach b, and x enters at a.

Write h=H(M), q=Q_M(2), h_R=H(R). With x removed, a Hamilton path has two optional slots on the sides of x. Endpoint and same-block contributions sum to (h_R+E)h: the left can only start after x and provides h_R h; the right allows the ending port b with optional suffix, giving E h. The only feasible forward seam goes from M to R and contributes (h_R-s_b)h. A backward seam must start the second slot at a. If the first slot ends in M, this contributes (q-h)h_R; if it ends at b, this contributes E q. Hence

```
H = (E+h_R-s_b)h + (E+h_R)q.
```

Thus the two counts are 20Q_M(2)+19H(M) and 20Q_M(2)+17H(M), and transitive M of arbitrary length L gives the exact adjacent odd counts 20*2^L+19 and 20*2^L+17.

The next two-slot interface strictly separates the choices. Switching b=0 to b=2 changes Q_out(2) by

```
-128 H(M) - 88 F_2(M) - 8 F_3(M).
```

Indeed an ordered pair of output slots has at most three runs in M, because deleting the single apex yields an ordinal quotient, so at most one slot can visit M twice. A fixed run itinerary expands by F_j(M) if it has j runs. For M=TT_1,TT_2,TT_3 the triangular basis matrix consists of rows (1,0,0),(1,2,0),(1,6,6). Literal subset DP gives differences -128,-304,-704, hence the coefficients, with an independent check at TT_4 in the replay. The pair is absolute but does not preserve the two-slot state needed for reentry.
