# Transplanting the recursive arrival law into an arbitrary finite core

**PROVED scoped, same-agent proof review:** every finite positive integer
core has a certified complete top chamber. Using its actual chamber, one
can construct recursive resonant attachments with an exact maximum at every
stage, and re-extract the chamber after a change of binding owner. Half-
clearance cores require an additional spacing guard. The input speeds of an
arbitrary LRC(14) instance are not freely choosable, so this construction
does not prove universal LRC(14).

## 1. Build the valid local model first

Let C be any nonempty finite set of positive integers, and

    f_C(t)=min_(v in C)||vt||,   rho=max_t f_C(t)>0.

The function is continuous and piecewise affine with nonzero slopes between
its genuine corners, so it has a finite set of maximizing times. At each
maximizer t_i there are positive integer descending speeds l_i,r_i such that
f_C(t_i-delta)=rho-l_i delta and f_C(t_i+delta)=rho-r_i delta for sufficiently
small delta>=0. At a tie use the fastest descending active speed; at a
half-integer cusp both directions are allowed with their actual slopes.

Let r_* be the largest value strictly below rho at a genuine corner of
f_C, including0 as needed. Then r_*<rho, and for every r_*<r<=rho the
**entire** safe set is the disjoint union of closed windows

    [t_i-(rho-r)/l_i, t_i+(rho-r)/r_i].                (L1)

This follows by following the actual affine pieces from each maximum until
its first corner. Every component not adjoining a maximum has peak at most
r_*. No other component is omitted, and no local slope is extrapolated past
its first change.

All data are rational and effectively computable. Candidate corner times
are k/(u+v), k/|u-v| and k/(2v), from pairwise distance equality and individual
cusps. Evaluate the lower envelope there, then remove candidates whose two
adjacent envelope slopes agree. Including nonleading corners would give a
valid but unnecessarily short chamber; omitting same-sign owner crossings
would be unsound.

Put Q=lcm of the maximizing-time denominators, p=Q rho (an integer), and
s=min_i(l_i,r_i). Since every maximizing distance is the fractional distance
of an integer speed at t_i, Q rho is indeed integral. We have Q>=2 and
1<=p<=Q/2. Record every side with owner s, not just their count.

## 2. Universal guarded recursive transplant — PROVED

Choose an integer m satisfying

    m(Q-1)>s,
    rho_low := p m(Q-1)/(Q(m(Q-1)+s)) > r_*.           (L2)

There is always such an m. For rho<1/2, choose ANY positive integer sequence
K_0,K_1,... . If rho=1/2, require instead

    K_j >= ceil((Qm+s)/(Qs)) for every j.              (L3)

Start a_0=m,b_0=0 and set

    a_(j+1)=K_j(Q a_j+s),
    b_(j+1)=K_j(Q b_j+p),
    P_0=1, P_(j+1)=Q K_j P_j,
    H_j=sum_(h=0)^j 1/P_h.

For the actual integer speed bank

    V_j=C union {Q a_0,...,Q a_j},

the exact maximum is

    M(V_j)=rho_j=p m/(Qm+s H_j).                      (L4)

Its complete maximizing set is obtained from every original least-owner
side, at displacement x_j/Q in that side's direction, where

    x_j=p H_j/(Qm+s H_j).                             (L5)

Every new speed vanishes at all preceding maximizing points. Their number
at each finite stage is exactly the number of original least-owner sides;
it need not be two. The original core and every older tail remain in V_j.

### Proof: exact joint deadline

At any original center, each speed Q a_h has integral phase. On a side of
owner e, write t=t_i +/- x/Q. Then the core constraint is

    Q r+e x<=p,

and every tail constraint is ||a_h x||>=r. The complete chamber(L1) and the
existence of an s-side imply that some safe point exists iff a joint safe
x>=0 fits Qr+s x<=p. This is the incoming joint-arrival theorem, applied to
certified arbitrary-core windows, not inferred from symmetry.

### Proof: coefficients, all thresholds, and global upper bound

With H_(-1)=0, induction gives

    b_h=(p P_h/Q)H_(h-1),
    a_h=P_h(m+(s/Q)H_(h-1)),
    p a_h-s b_h=p m P_h.

In particular rho_j=(p a_j-s b_j)/(Q a_j+s). Since
H_j<Q/(Q-1), all rho_j>rho_low>r_*; the complete chamber remains valid.

A safe x for tails through h must satisfy x>=(b_h+r)/a_h whenever the
successive first-arrival induction thresholds have been passed. The initial
condition is x>=r/m. At the next tail the threshold is

    theta_h=K_h(p a_h-s b_h)/(K_h(Q a_h+s)+a_h).

Indeed substituting the preceding lower bound places a_(h+1)x strictly
past b_(h+1)-r; safety then forces it at least b_(h+1)+r.
For any later finite stage j>h,

    H_j-H_h <1/((Q-1)K_h P_h),

and a_h(Q-1)>s. These inequalities show rho_j>theta_h. Thus every test
radius r>rho_j passes every preceding induction threshold. Combining
x>=(b_j+r)/a_j with Qr+s x<=p forces r<=rho_j, contradiction. This proves
the global upper bound over all times.

### Proof: old tails remain safe, including the half-core boundary

At x=x_j, an old tail has unwrapped phase

    a_h x_j-b_h
      =p m P_h(H_j-H_(h-1))/(Qm+s H_j).                (L6)

It is at least rho_j, with equality exactly at h=j. If rho<1/2, then
2p<=Q-1 and the geometric tail makes(L6) strictly less than p/(Q-1)<=1/2.
Thus every old tail is safe on the stated branch.

If rho=1/2, put K_min=ceil((Qm+s)/(Qs)) and R=Q K_min. The stronger
spacing gives P_h(H_j-H_(h-1))<R/(R-1), while H_j>=1. Hence(L6) is less
than

    p m R/((R-1)(Qm+s)) <=1/2,

because 2p=Q and Qm<=s(R-1), exactly the guard(L3).
The least-owner core side and newest tail both have clearance rho_j at(L5),
proving attainment. At equality the first-arrival lower bound and the core
deadline force x=x_j and owner s; hence no maximizers are missing. QED.

## 3. A real owner change, and a failed half-core extension

For C={1,6}, the exact data are

    rho=3/7, r_*=2/5, Q=7, p=3, s=1.

Condition(L2) permits m=3. The first added speed21 gives exact maximum9/22.
This directly repairs the THM-724 hostile: its attempted addition7 used m=1,
whose proposed optimizer lies below the certified chamber.

Now add speed23. It preserves the two maxima9/22,13/22 but becomes a binding
wall. For C'={1,6,21,23}, re-extracting the actual geometry gives

    rho'=9/22, r_*'=9/29, Q'=22, p'=9, s'=21.

The cheapest side has changed direction and its speed has changed from1
to21. Taking m'=4 satisfies(L2); the next resonant addition88 has exact
maximum36/109. The method is closed under this change because it recomputes
the chamber, center frame and owners; retaining only the old maximum would
not supply that proof.

The half-core guard is substantive. Starting from C={1}, choose m=3 and
ignore(L3) by using K_0=K_1=1. The predicted third maximum would be12/31
for{1,6,14,30}; the actual maximum is1/3, at7/18 and11/18. An older tail
has folded out of its permitted interval. This refutes the unguarded
half-core extension, not the predecessor theorem, which assumed rho<1/2.

## 4. Composition and a positive infinite margin

For any target r<rho, choose m still larger so rho_low>max(r,r_*).
Then every finite stage and the infinite union have maximum at least r;
the exact infinite value is p m/(Qm+s H_infinity). The finite maximizing
points on each retained side converge, proving attainment for the infinite
bank. We do not infer uniqueness in the infinite limit.

One may also change the core between stages, as above. Given a certified
finite core of value rho_i and a desired fractional loss0<epsilon_i<1,
choose m using its new chamber and

    m > s(1-epsilon_i)/((Q-1)epsilon_i).

Every recursive batch then retains value at least(1-epsilon_i)rho_i.
If such growing finite cores obey sum_i epsilon_i<1, their common safe
sets at radius rho_0(1-sum_i epsilon_i) are nonempty, closed and nested.
Compactness gives a common safe time for their union. Extra attachments
between batches must themselves be certified not to consume unrecorded
margin; there is no permission to append arbitrary prescribed speeds.

This is a constructor for controlled extensions, not an entry theorem for
an arbitrary input set. That quantifier is the remaining LRC(14) gap.

## 5. Scope, recovery and controls

The complete-window transplantation generalizes the preceding
[near-interval echo theorem](../entry_cosets_resonant_echoes_20261030/LRC.md).
The incoming [second-band owner theorem](../odd_eighth_owner_families_20261030/PROOFS.md)
classifies its own wider family uniformly and is retained as a stronger
explicit formula on that family. Its first-arrival-profile hostile shows
why the scalar earliest time is not closed under later additions; here the
full core is retained and the next chamber is re-extracted.

THM-4112 / `01-canon/theorems/THM-4112-antipodal-component-ancestry-chain-and-scale-separated-lrc-families.md`
already gives component-ancestry and scale-separated suppliers. The current
claim is different: exact maxima and marked maximizers for a tailored
recursive construction over any finite core, with explicit complete-chamber
and half-core guards. No new general arbitrary-input LRC closure is claimed.

`core_audit.py` checks all196 primitive triples through12, with588 complete
recursive-bank maximum certificates, five independent pair-clock maxima,
the changed-owner example and the half-core hostile. Its own complete
lower-envelope extractor removes nonleading events; a separate exact safe-
interval implementation checks all maximal points and excludes any larger
rational event maximum. Normal and optimized guards are explicit exceptions.
