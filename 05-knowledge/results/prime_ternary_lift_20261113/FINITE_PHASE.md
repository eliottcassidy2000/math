# A finite changing-core clock cannot renew the +4 digit

2026-11-13, Codex. **PROVED IN NOTES**, in the exact grammar stated here.
The one-frame parity obstruction is inherited from
[clone-clock CLOCKS, section 1](../clone_clocks_negative_collatz_20261113/CLOCKS.md).
This argument extends it to arbitrary finite cycles of different rooted
frames and arbitrary fixed nonnegative length increments. It does not
exclude general core replacement, more than one module, or target-dependent
changes to the module length.

## Statement

For each rooted finite tournament frame (A_i,r_i), put

    f_i(m) = H(A_i[r_i <- TT_m]),   m>=1.

Allow a directed rule i -> j with fixed integer s>=0 and digit d in {2,4}
only when

    f_j(m+s) = 3 f_i(m) + d                         (1)

holds for every m>=1. **Every directed cycle of such rules uses only
digit2.** Therefore a finite grammar of this kind cannot generate
trajectories with arbitrarily many +4 steps: along any trajectory each
+4 edge crosses between distinct strongly connected components, so there
are at most one fewer such edges than there are components.

The theorem concerns identities of whole infinite families. A coincidence
at one m is not an edge. Frames may be strong or nonstrong. There is no
order bound and no matching-number restriction.

## Proof

The run expansion gives

    f_i(m) = sum_{j=1}^R c_{ij} j^m,

with integer c_{ij}, since K_k(TT_m)=sum_{j=1}^k
(-1)^(k-j) binom(k,j) j^m. If R is the greatest nonzero run coefficient,
the highest power coefficient is c_{iR}=a_R>0. Also (R-1)! divides a_R:
permuting the R-1 disjoint internal excursions of a skeleton acts freely
and preserves its validity. This is the inherited excursion-permutation
lemma, with the same nonempty-excursion convention.

Suppose a cycle has p edges, total length increment q, and composed
additive term D. Then

    f_i(m+q) = 3^p f_i(m) + D,
    3^p-1 <= D <= 2(3^p-1).                        (2)

If q=0, positivity makes (2) impossible. If f_i is constant, positivity
again makes (2) impossible. Otherwise q>0 and independence of the
sequences j^m implies j^q=3^p for every nonconstant mode j with nonzero
coefficient. There is at most one such integer j; unique factorization
gives j=3^a with a>=1 and aq=p. Thus

    f_i(m) = A (3^a)^m + B,

where A=a_(3^a) is positive and divisible by (3^a-1)!, hence even.
Hamilton counts are odd, so B is an odd integer. The constant term of
(2) now gives

    B = -D/(3^p-1),     -2 <= B <= -1.

The only odd possibility is B=-1, forcing D=3^p-1. Each digit contributes
its positive power-of-three weight to D. Its minimum occurs precisely
when every digit is2. This proves the claim. The component consequence
follows from the acyclicity of the condensation of a finite directed
graph. No numerical census is used in the proof. QED.

## Boundaries, controls and revised search

- **Equality:** a known frame realizes f(m)=60*3^m-1, hence a digit2 loop
  with s=1. Its excursion vector is (179,180,60,0), independently checked
  in `renewal_probe.py`. The theorem allows it.
- **Arithmetic hostile:** f(m)=3^m-2 is odd and repeats the +4 digit, but
  would require a_3=1, contradicting the free permutation of two internal
  excursions. Positivity of the scalar values is not enough to make a
  tournament frame.
- **Mixed modes:** the new family 6*3^m+6*2^m-5 has a real +4 exit, but
  cannot lie on any cycle of (1), even an all-+2 cycle, because two
  different nonconstant bases cannot satisfy j^q=3^p simultaneously.
- **Scalar growth alone:** with s=0, every edge strictly increases a_1,
  so the entire fixed-module graph is acyclic. Allowing s>0 is the
  substantive strengthening above.
- **Not excluded:** unboundedly many cores; two independently changing
  transitive modules; replacing TT_m by another kind of module; rules
  valid only at target-dependent parameter values; count-preserving
  preparation that changes the family representation. Each needs its
  own entry proof.

This changes the next question. Instead of seeking a finite list of
one-module phases with endlessly renewable digits2 and4, seek a family
with an extra physical resource and a proved rule for renewing that
resource, or a target-specific constructor whose source need not be an
output frame of the preceding step.
