# Independent audit of the actual fifth-clock phase family

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT — PASS.**

Auditor: `three_ray_geometry`, separate from the family author,
`certificate_audit`. I read the full
[proof](long_frontier_sep06_phase.md), the complete source, and the
inherited native, dual, actual-decoder and hereditary-profile statements.
No mathematical correction or scope correction was needed.

The family contains exactly 384061 admissible integer scales of the stated
fixed primitive, unitless 5+8 shapes. Each actual thirteen-speed row is in
the required physical box, has full equality W=V_dec, passes the pinned
hereditary profiles, and has an explicit time of clearance greater than
1/6. The named scalar sufficient gates fail, but no assertion is made
that every known safety method fails or that these safe rows are
counterexamples. The scale count is not a speed-height threshold.

## 1. Complete physical entry and profiles

I checked the actual star and cofactor edges against the strict inert atlas
of **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`.
The listed reduced sums have exactly the required prime factors and
exponent bounds. They connect each primitive shape. The uniform cross-ratio
bound excludes every cross atlas edge. The two component shapes have
disjoint prime supports; gcd(t,5P)=1 therefore gives every cross gcd equal
to one. The stated upper-scale physical-sum bound proves the box condition
for the whole interval, and the cross separation proves distinctness.

The first orientation is excluded by the exact minimum distinguished
coefficient t*gcd(v,w)>Q. In the reverse orientation, the distinguished
coefficient is an integer multiple j*gcd(u,w), j>=1. After division, its
target is jtv>349Q, while the entire signed two-term support has absolute
value at most 262Q. Thus every higher multiple is excluded, rather than
only the primitive distinguished coefficient. The same reasoning includes
support-two crossings by allowing one pair coefficient to be zero.

There are exactly 80 and 140 mixed triple supports in the two orientations.
The connected actual edge relations have height at most 355 and span both
internal weighted kernels. Since all bounded support-at-most-three mixed
generators are excluded, this proves full W=V_dec for every family member.
No absence-of-one-orientation or abstract-graph shortcut is used.

For hereditary profiles, every mixed body has gcd one, and a pure smaller
body has fewer than seven labels. The only nontrivial seven-through-twelve
body gcds therefore come from deleting one cofactor P/p from the eight-label
larger shape. Its gcd is exactly p and every complementary label is
coprime to p. These are precisely the six displayed signatures
(p;1,1,1,1,1,1), all present in the pinned joint-shadow bank. The other
signatures have gcd one. This proves the profile assertion uniformly in
t; the direct eight-row subset tests corroborate the symbolic argument.
The profile maxima are exactly (1,1,1,1,1,29).

## 2. Exact phases, scale count and the actual hostile

I verified every row of the word (k,delta) modulo five. It gives
2tk+delta=5 modulo ten and hence the actual half phase tx=1/2 modulo
one. Oddness of every smaller primitive speed then gives clearance 1/2
on all five smaller physical labels.

For each of the 32 pairs consisting of a residue class and a larger speed,
I independently computed the two rational affine endpoints. They lie
strictly between zero and one; thus no integer branch is omitted. The
four uniform bounds are exactly

```text
r=1: 637017991/3311190974
r=2: 1529494468/8277977435
r=3: 2932887917/16555954870
r=4: 1613561814/8277977435
```

Each exceeds 1/6. These endpoint inequalities prove the clearance bound
over the whole continuous scale interval and therefore over every
admissible integer scale.

The count was independently reconstructed by inclusion-exclusion on the
six primes dividing P inside each progression t=r modulo five. For a
squarefree product d of those primes, solve t=d*j=r modulo five and
count the resulting single progression modulo 5d in the exact integer
interval. The four counts are 96015,96015,96014,96017, totaling 384061.
This does not scan the family members.

At t=4912753 I independently computed both literal thirteen-speed
clearances. The nearest-lift phase 7860405/9825506 has minimum
309014/4912753 at the label 6735366, below 1/14. The repaired phase
3930203/9825506 has minimum 1740589/9825506 at the label 374187,
above 1/6. The witness is inside the stated family. It refutes only
the selected nearest-lift rule, not physical safety or the inherited
phase-packet theorem.

## 3. Audit of the stated failed consumers

I matched the native conditions to
[the cross-divisor theorem](open_frontier_sep06_decoder.md) and
[the incoming dual-pair theorem](continuing1_20260906_lrc_dual_pair.md).
The actual cancellation divisor is one in both orientations. The full
larger maximum remains K=9P in every width test.

The 140 forward scores use the exact signed-box radius, whose maximum
here is 262Q at the primitive pair 1:261. The 35 endpoint comparisons
retain the coefficient cap 30D<=Q and fail the required width condition;
their minimum lcm(D,v) is exactly 1041046529754. Every one of the 80
dual widths also fails. The two displayed whole-arc inequalities fail
throughout the scale interval. These are precisely scoped failures of
sufficient criteria, and do not purport to classify every safety mechanism.

The incoming seventh-clock condition max(U)<t is absent. The denominator-23
family has both a different split and a residue hypothesis not satisfied
here. The note correctly describes these as missing hypotheses rather
than counterexamples to the inherited theorems. The stronger leaf-cofactor
cutoff is credited, so an obsolete automatic endpoint cutoff is not
presented as the current comparison.

## 4. Full source inspection and independent replay

I read the complete
[source](../../04-computation/long_frontier_sep06_phase.py), including its
strict atlas factorization, rational phase branches, two independent
scale-count formulas, scalar consumers, pinned profile loader, and exact
signed-box congruence solver. The two interval constraints in that solver
are the complete constraints on its free integer parameter. The independent
coefficient/amplitude inequalities exclude every distinguished multiple.
All gates use explicit exceptions, remain active under optimized execution,
and use integer or rational arithmetic.

The finite physical controls are the first and last eligible scale in
each of the four residue classes. For all eight rows, the source checks
the full graph, both mixed support orientations, all 4095 bodies of
sizes seven through twelve, and the literal physical clearance. The
additional nearest-lift hostile is separate. This is a declared control
universe, not an unreported height census.

Independent normal and optimized execution reproduced the frozen output
byte-for-byte. All **38683 always-active gates PASS**. The source and output
hashes were independently recomputed:

```text
source SHA256 51b8acf006ea19cac1804d06758b24045b0cff0f4b4598c459e032250bbcbe43
output SHA256 c3a77805e6b9d240c93df8c0904295dcfb81c5f77a5244a5d0697363fb6c1191
semantic SHA256 a8c0f182e22fe82240717515e7e6fa32eab9d79cf67a316365af9610523c97a7
profile JSON SHA256 935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f
```

The entire stated family theorem and all explicitly named comparison
claims pass this audit. LRC(14) and general actual-entry phase selection
remain open.
