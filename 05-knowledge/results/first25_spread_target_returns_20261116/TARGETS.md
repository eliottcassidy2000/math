# Exact inverse on the old +4 row, using three different cores

Write g for the count after deleting the named root, E and O for
the counts of root-start paths of even and odd length in vertices,
each weighted by the Hamilton count of its complement. Write
Z=E+O−g and I=H−2O+2g. A near-sink vertex beats only the root;
every other old vertex beats it. By the inherited rooted-reader law,

    (H,g,E,O) -> (H+2E,g,E+O−g,E+O).

In particular after k≥1 repetitions the profile is
(H,g,E,O)=(I+2^k Z,g,2^(k−1) Z,2^(k−1) Z+g).
These are exact graph operations. They are not arbitrary
arithmetic rewrites of the scalar H.

Read the masks below as outgoing masks in the deletion tournament D,
and append root 8 with outgoing mask m. The direct subset DP verifies
the marked initial profiles.

| digit j | outgoing masks of D | m | H,g,E,O | Z,I | number of appends |
|---:|---|---:|---|---|---:|
|0|0,1,27,51,99,135,47,95|164|723,37,378,415|756,−33|s|
|1|2,72,131,53,39,71,157,59|22|1395,287,799,1000|1512,−31|s−1|
|2|0,33,11,19,39,77,159,63|9|111,49,119,119|189,−29|s+2|

Thus at every s≥2, H=756·2^s−33+2j, E=378·2^s, O=E+g.
The tournaments are strong already after the indicated appends at s=2
(direct reachability), so every later order is strong as a mixed
attachment to a strong core. For an exact target N of this form,
j is unique modulo 756 and then (N+33−2j)/756 dictates s. The
raw core choice and length follow with no search at target height.
Previously recorded common full profiles at individual heights
are distinct from this all-height restricted continuation; the
other earlier matched-even row had intercept −39 rather than −33.

## Which inputs may be safely forgotten?

On normalized roots O=E+g, the scalar near-sink operation reads
only (H,E), with update (H+2E,2E). After a common append the
three outputs therefore recur exactly at the next height, not just
at the next residue. But this projection does not license arbitrary
binary append words. For the fixed-root binary controller of
[the birth-fiber packet](../birth_fibers_binary_growth_20261116/README.md),
let w end in 0 (last vertex a near-sink), x=int_2(w), z the
reverse complement read as an integer, and f(0)=1,
f(2u+1)=f(u), f(2u)=f(u)+u+(u mod 2) for u≥1. Then

    H' = H+2Ez+g(f(z)−1),     E' = 2^|w| E+gx/2.

For two such inputs of equal E but different g a common word
preserves the original H gap *and* the common E iff w is all
zeros. If x>0 then E' splits, and z has a non-leading zero,
so f(z)>1 and H' also reads g. Normalizing at the end of the
word cannot undo that separation. This proves an exact
safe-language boundary within the shared binary controller,
not a no-go for choosing different words or changing source.

## Composite near-sink blocks and target decoding

Instead of one vertex, append arbitrary nonempty B, beaten by
all old nonroot vertices and beating the root uniformly. The
inherited rooted block transfer uses u=H(B) and
v=P_B(2)=sum_{S subset B} H(B[S])H(B−S), empty count 1.
For a normalized root it gives

    H' = u(H−2E)+2vE,    E'=vE,    g'=ug,    O'=E'+g'.

Thus (H,E) is a closed continuation for *all* such near-sink
blocks. Their affine intercepts scale by u and slopes by v;
an old gap 2 scales to 2u, which is why a composite module
does not by itself give a free consecutive-digit move.
Successive uniform near-sink blocks are equivalent to their
ordinal concatenation, so u and v compose multiplicatively.

On the old row, the target equation is

    N = C v + I_j u,    C=756·2^s,    I_j=−33+2j.

Rather than demanding a prime operand, seek a smaller u in the
arithmetic progression I_j u ≡ N (mod C) with
1≤u≤N/(2C+I_j), then ask that v=(N−I_j u)/C occur physically
with H(B)=u. The upper bound follows from v≥2u by the two
one-colour assignments. This is an exact target-decreasing
profile problem. Ignoring the second coordinate would make the
descent circular; the small examples really supply both counts.

- Block masks (2,12,9,49,7,23) have profile (u,v)=(25,304).
  With j=1 and s=6 they give the prime 14,707,961.
- Block masks (2,12,17,37,11,87,31) have profile (49,824).
  With j=1 and s=5 they give 19,932,689 = 7·2,847,527
  with the other factor prime.

The standalone program verifies the block profiles and smaller
literal conjugate placements by subset DP, and primality by
trial division. Together with the block law these certify the
displayed large targets without enumerating all their paths.
The unanswered step is **coverage of that joint profile
progression**, or a changed reader when it is empty. These bands
and the composite options are constructive but do not close
the full Chan–Pak proof.
