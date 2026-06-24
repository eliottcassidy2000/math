# LRC14 Ramanujan Exact-Period Projectors

The useful Ramanujan-sum move is not another scalar residue fingerprint.  HYP-2978
already showed that scalar divisor, qdiv, and even c_14 speed-profile channels
mix AP/GW, q-witness, K33, petal, and covering routes.

The better object is the primitive shell as an operator on phase functions:

```text
E_q(f) = sum_{a,b mod q} f(a) f(b)c_q(a-b).
```

For LRC14, take `f(a)=N_S(a/q)` or the weak-safe indicator on `Z/qZ`.  This
places the Ramanujan sum between the twist ladder and the Toeplitz/Fejer dual:
it is still exact-period and integer-valued, but it sees a whole phase function
rather than a speed multiset.

External trail:

- Ramanujan's sum is the sum of nth powers of primitive q-th roots, with
  Kluyver's Mobius formula `c_q(n)=sum_{d|gcd(q,n)}d*mu(q/d)`.
- Carmichael-style orthogonality motivates shifted/autocorrelation packet tests.
- Signal-processing work on Ramanujan subspaces treats these sums as integer
  periodicity components for finite signals, matching our "hidden exact period"
  use case.

Sources read in this pass include:

- https://en.wikipedia.org/wiki/Ramanujan%27s_sum
- https://mathtube.org/lecture/video/ramanujan-sums-and-hardy%E2%80%93littlewood-prime-tuple-conjecture
- https://royalsocietypublishing.org/rsta/article/378/2163/20180446/111519/Srinivasa-Ramanujan-and-signal-processing
- https://systems.caltech.edu/dsp/PPVSquare.pdf

The HYP-2979 audit found:

```text
rows audited          21906
no weak q<=42             0
no strict q<=42           2
no strict examples        AP, GW 12->24
```

Named late packets line up with the existing LRC proof geography:

```text
near 12->36            first strict primitive q = 41
petal 10->20           first strict primitive q = 27
petal 13->26           first strict primitive q = 27
P10+GW                 first strict primitive q = 27
P10+K33                first strict primitive q = 27
covering 12->84        first strict primitive q = 41
covering 12->168       first strict primitive q = 41
covering 6->98         first strict primitive q = 25
```

This is not a proof of LRC14.  The AP-neighborhood bank is bounded, and q=14
primitive phase packets still mix many routes.  But it sharpens the shape of a
possible theorem:

```text
Every reduced LRC14 residual has either a primitive exact-period witness,
a harmonic dual certificate, an AP/GW boundary packet, a Ramanujan
danger-energy defect forcing labelled handoff, or the K33/H=7 state lift.
```

Assumption challenge: vertices considered were runners, residues, primitive
phases, q-denominators, Ramanujan modes, endpoint-owner pairs, autocorrelation
shifts, and proof obligations.  The reflection chooses proof carriers and
exact-period modes as tournament vertices.  This preserves weak/strict witness
status and boundary equality, but it destroys endpoint ownership unless C27,
K33, taut-current, or labelled-packet data are explicitly reattached.

The next concrete computation should run the strict primitive-witness audit on
the full HYP-2963 bank and then interval-enclose the late primitive packets.
The key q values to prove exactly are `25,27,34,40,41`.
