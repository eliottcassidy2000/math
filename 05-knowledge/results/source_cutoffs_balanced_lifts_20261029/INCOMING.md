# Final-pull integration: norm charts, named ports and the Gram boundary

Incoming `0ad8de0355` and `51cd593a3e` were read after the main checkpoint.
The first adds the [binary-module and named-port checker](../session110_overnight_20261030/portable_checks/binary_module_counts.js).
The second adds [prime-power-plus-one growth](../prime_power_plus_one_20261029/GROWTH.md)
and [norm-clock arithmetic](../prime_power_plus_one_20261029/ARITHMETIC.md).
The relevant proofs were read; the observations below are scoped deductions,
not a claim that every incoming result was independently re-proved.

## 1. A literal common conic, with different arithmetic fibres

The incoming finite-field clock uses

    z_t=(t+alpha)/(t-alpha), alpha²=D,
    t star s=(ts+D)/(t+s).

Our ramified clock uses

    h(u)=(1+u delta)/(1-u delta), delta²=D,
    u circle v=(u+v)/(1+Duv).

Both are the same Cayley parametrization of the norm-one conic

    X²-DY²=1,
    X=(1+Du²)/(1-Du²), Y=2u/(1-Du²).                       (1)

The coordinate change is projectively [t:1]=[1:u]. The expression t=1/u
is valid only where u is invertible; projective coordinates retain the
other points. Direct multiplication on the conic proves both displayed
operation laws. This is the exact map underlying the resemblance.

The arithmetic specialization is different. For nonsquare D in F_q,
the projective chart gives the whole cyclic norm-one group of order q+1.
For D=p d in our ramified local ring, 1+Duv is always a unit for integral
u,v, and reduction modulo p gives the ADDITIVE law u circle v=u+v.
The special conic is X²=1; its identity component X=1 has the free
Y-coordinate. The integral chart modulo p^k gives our odd group H of
order p^k, proved directly in LOGFREE.md.

Thus the common object preserves norm and multiplication, while the
valuation of D decides which coordinates and group orders survive.
It is not an isomorphism between the two finite groups. The q+1 clock
has even order and acts on a switching representation, whereas our
odd-order clock acts as a literal tournament automorphism. Discarding
the ramification flag would erase precisely the information responsible
for the Gauss-sum parity channels and stabilizer factors.

## 2. Optimal spectral scale is not an exact Gram identity

The incoming deletion theorem for a doubly regular tournament uses

    K^t K=qI-J                                                (2)

to classify every equality case. Our O_p(sqrt N) norm bound supplies
no such identity or lower bound on each nonzero singular value.
Indeed our N=p^(2k) is1 modulo4, so (2) is impossible here: its
common-out-neighbor count would be (N-3)/4, not an integer.

Consequently the incoming sharp half-size deletion margin cannot be
imported into the balanced lift theorem. The statements concern different
levels of information: a spectral upper bound, an exact Gram matrix,
and native modular primeness. L1 proves the last by ideal cuts, not by
pretending the first is the second. A useful next question is a direct
deletion margin for the balanced family with its actual valuation ports.

## 3. Named binary ports are a distinct counting observer

The incoming module checker was replayed: all280 termwise comparisons
pass against the complete palette groups GL3(F2), AGL3(F2), the stabilizer
of one named linear port, and the identity group for three fully named
ports. Their orders are168,1344,24,1. Its uncolored order12 control is
the already-known A000568 value154108311168.

This strengthens the context for COMPOSITION.md's505-class example:
changing the allowable palette operations changes the counted object.
Its weighted regular-composition recurrence retains the actual quotient
action and actual child types. For fixed named binary ports, colors can
be included in those child types and the same canonical decomposition
and Burnside proof applies with the appropriately decorated stabilizer.
Quotienting the bit coordinates by GL3 instead requires its additional
joint action; it is not the same as three named ports. If the bits stand
for new tournament vertices, their mutual arcs must also be retained.

The incoming scalar-product hostile R3[R5] versus R5[R3] reinforces the
same point: equal size, cut imbalance and transitive maximum do not
recover the ordered substitution tree or its automorphism action. Neither
our kernel-size vector nor a spectrum is proposed as that reconstruction.
