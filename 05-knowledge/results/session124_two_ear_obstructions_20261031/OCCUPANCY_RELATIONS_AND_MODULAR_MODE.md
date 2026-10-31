> **Curated publication copy.** Executable proofs and certificate data are in
> `code/`; independent audits are in `audit/`. Run `verify_all.py`.
> Historical exploratory scripts mentioned below are intentionally omitted;
> their finite searches are not dependencies of the certified theorems.

# Occupancy identities and a forced arithmetic lower mode

Status: exact identities and independently checked modular obstruction.
These concern the two-fixed-ear/transitive-spine carrier only.

## Eight missing reachable-state identities

Forget good/bad LAST-vertex flags but retain whether each of P,Q,R is empty.
An xy-order partition(P,Q,R), whose Q starts at an x-outneighbor and R starts
at a y-outneighbor, corresponds bijectively to the yx-order partition(P,R,Q).
Consequently, for each of the eight occupancy triples(eP,eQ,eR),

    total_xy(eP,eQ,eR) = total_yx(eP,eR,eQ).

Let E be these eight integer relation rows. They hold for the initial vector
and are preserved by every letter matrix. They are independent; their kernel
has dimension28, equal to the exact reachable linear span found by closure.
The standalone checker does not need that last equality: it proves invariance
on a complete integer kernel basis with112 transition checks.

This explains the false-positive alternating templates in the earlier
36-dimensional nonnegative relaxation: it had omitted cross-order identities.

## Stronger individual-template exclusions

Adding Ev=0 to the Farkas system gives exact certificates for all eight
non-00 eligible length-two templates. No constant-digit affine action
H→9H+d can be iterated forever on any actual seed, for ANY constant d.
Only00 remains possible in this exact classification.

The corresponding length-four calculation gives224 exclusions, leaving only
0000. It strengthens, rather than retracts, the earlier220-certificate result:
the four old survivors were explicitly labelled unsolved by that relaxation,
not claimed to be realizable.

Files:
* `periodic_reachable_length2_certificates.json` (8 identities);
* `periodic_reachable_length4_certificates.json` (224 identities);
* `verify_periodic_reachable_certificates.py L` (standard-library checker).

## A cleaner modular obstruction for complementary pairs

For w∈{03,12,21,30}, let M=M_w, c=output, and define

    q_w = (1/8)c M(M−I)(M−4I)^2(M−9I).

The following exact statements are verified independently:

1. q_w is an INTEGER row.
2. q_w(v_empty)=−4.
3. For every letter a and integer v with Ev=0,
   q_w(M_a v)−q_w(v) is divisible by3.
4. c M(M−I)^2(M−4I)^2(M−9I)=0.

Hence every actual seed v satisfies

    q_w(v) ≡ 2 (mod 3),

so q_w(v) is never zero. This is a genuine arithmetic invariant of reachable
states; no positivity relaxation or numerical solver is used.

By the annihilating polynomial in(4), the scalar count sequence has form

    H(M^k v)=A·9^k+(B+Ck)·4^k+Dk+F,   k≥1.

Applying the polynomial in the definition of q_w annihilates the9- and4-modes
and acts on the linear1-mode with factor−72. Thus

    D = −q_w(v)/9 ≠ 0.

The coefficient is NOT universal as a rational number; its nonvanishing is
universal because its numerator lies in a fixed nonzero residue class.
This forced linear term rules out every pure9^k-plus-constant recurrence,
including all constant-digit radix9 actions for these four templates.

`verify_alternating_mod3.py` independently reconstructs all matrices, proves
integer divisibility in144 coordinates, checks448 kernel transitions modulo3,
and checks four annihilator rows. Normal and optimized receipts agree.
The integer kernel basis is obtained by differences within each signed
occupancy group and one paired representative across orders; it spans the
integer kernel, so modular verification on this basis is sufficient.

## Correction lineage: do not fit away polynomial modes

During exploration a private message claimed
H((12)^k)=(9^k+4·4^k)/5 from only three values. This was FALSE and was retracted
immediately when k=4 was checked: the actual count is1485, not1517.
No theorem, certificate, or published file used that fit. The failure was
assuming away the k·4^k and linear1 modes allowed by the true annihilator.
The exact modular argument above replaces that shortcut and records why a
nonzero lower polynomial mode cannot be discarded.
