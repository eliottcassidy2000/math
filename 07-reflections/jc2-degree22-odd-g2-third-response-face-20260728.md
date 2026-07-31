# The odd-Faber infinity germ is the G2 arrangement, and response is transverse

**Status:** VERIFIED-EXACT FINITE SIDECAR TO PROVED THM-2719.  This records a
typed local calculation for the live THM-2726/third-flux lane.  It is not by
itself a trajectory exclusion or a proof of `JC(2)`.

At the unique infinity point of THM-2719, the second-flux cover tangent is

```text
q^6-15q^4s^2+15q^2s^4-s^6
 =((q+is)^6+(q-is)^6)/2.                              (1)
```

Thus the six lines form the dihedral `I_2(6)` (equivalently `G2` reflection)
arrangement.  The index-two central inversion pairs opposite lines into the
three coarse branches.  This is a literal local realization of the recurring
two/three grammar: six cover directions, a central `C2` quotient, and three
coarse directions.  It is not an imported modular-group action.

Let `R_j=4c_(j,3)+2d c_(j,1)` be THM-2129's response observable.  Exact
Laurent recurrence gives

```text
R_21(1,0,0)/Phi_21(1,0,0)=11/24.                     (2)
```

For `a_21!=0`, solving the first flux for `h` and substituting into the
homogeneous degree-25 response leaves the order-six face

```text
1771/1024 q s(q^2-3s^2)(3q^2-s^2).                  (3)
```

The gcd of `(3)` with `(1)` is one.  Hence response vanishes on none of the
six cover directions.  Since `h` has order six, on each coarse branch the
physical coefficient `q/h^3` has pole order `17`, while the affine response
has pole order `25*6-6=144`.

The preserved conclusion is an exact transverse pole invoice.  Turning it
into a Keller exclusion still requires the source-primitive/dominance input:
one may not infer that an arbitrary rational function with three poles is an
affine source coordinate.  THM-2726 is the intended typed destination for the
physical `q` invoice; the response-`144` sidecar may sharpen its equality
boundary.

Reproduce with

```text
python 04-computation/jc2_degree22_odd_g2_third_response_face_20260728.py
python -O 04-computation/jc2_degree22_odd_g2_third_response_face_20260728.py
```

and compare with

```text
05-knowledge/results/jc2_degree22_odd_g2_third_response_face_20260728.out.
```

The LF-normalized SHA256 hashes are, respectively,
`55ab39ecaa03dc871d3704753d0bf3e0ef983b7d30c7b793855911f583762860`
and
`cfb0122406feab23427f04d24482bc7cb390f4a26deba458b5d4f03d86d68a4d`.
