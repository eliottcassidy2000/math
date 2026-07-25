---
id: THM-2171
title: "Radix quotient-order pumping and the algebraic finite terminal"
status: >
  PROVED + VERIFIED-EXACT. Adjoining the adjacent quotient-tie mask to the
  relation-carry and owner state repairs the order/distinctness loss in
  THM-2167. Equal augmented states can be pumped while preserving positivity,
  strict order, every tracked relation, and—after harmless common-gcd
  normalization—primitivity. The owner/tie pair is monotone and assumes at
  most 2d values along a d-coordinate radix path, so every feasible bounded
  relation system has an explicit small positive ordered representative. For
  the LRC(14) rank-two carrier the effective path cap is
  26*2729^2=193633466. An exact thirteen-speed witness proves that the
  1/14-safe-set predicate is still mixed on these pump fibres; phase is the
  genuinely missing sidecar.
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2163
related:
  - THM-2163
  - THM-2167
  - THM-2161
script: 04-computation/lrc14_radix_quotient_order_pumping_thm2171.py
output: 05-knowledge/results/lrc14_radix_quotient_order_pumping_thm2171.out
script_sha256: 815da09174e3f102a9561cb6a90f6feaaef1a06937c4532d0334b14cef0f9c42
output_sha256: 7326bc5a7f1968cea68ab978016d6f766c33a0346bf085ea6f6adc5f1764c767
hash_basis: working-tree bytes (LF)
---

# THM-2171 -- radix quotient-order pumping

THM-2167's rank-two carry pump preserves its relations and positivity, but
its owner mask alone can merge coordinates. This theorem identifies a finite
order datum that repairs that loss. It also proves that common
gcd changes are harmless for LRC. The remaining failure is not algebraic:
the safe phase predicate is genuinely nonconstant on the repaired pump
fibres.

## 1. The quotient-order sidecar

Let

```text
0<V_1<...<V_d,                    q>=2,               (1)
```

and use THM-2163's radix notation

```text
V=q^j Z_j+R_j,             0<=R_(j,i)<q^j.            (2)
```

Besides the owner suffix

```text
O_j={i:Z_(j,i)>0},                                (3)
```

define the adjacent quotient-tie mask

```text
T_j={i in {1,...,d-1}:Z_(j,i)<Z_(j,i+1)}.             (4)
```

Thus `T_j` records the cuts between the consecutive equality blocks of the
nondecreasing quotient vector `Z_j`. It uses only `d-1` bits.

Both sidecars are monotone:

```text
O_(j+1) subset O_j,        T_(j+1) subset T_j.         (5)
```

The owner inclusion is immediate. For the tie mask, coordinatewise Euclidean
division gives

```text
Z_(j+1,i)=floor(Z_(j,i)/q).                            (6)
```

If two adjacent entries of `Z_j` are equal, their next quotients are equal.
Consequently a vanished cut can never reappear. This nested-radix fact is
important: for unrelated divisors a quotient inequality could disappear and
reappear, but powers of one radix do not have that defect.

Starting from `j=0`, the owner set loses at most `d` coordinates and the tie
mask loses at most `d-1` cuts. Therefore the pair

```text
P_j=(O_j,T_j)                                         (7)
```

assumes at most

```text
1+d+(d-1)=2d                                         (8)
```

distinct values along any one radix path. Its ambient alphabet can be much
larger; (8) is the effective path count.

## 2. Exact order-preserving pump

Let `A` be any finite family of integer relations on `V`. For every
`a in A`, let

```text
kappa_j^a=-a.Z_j=a.R_j/q^j                            (9)
```

be its THM-2163 carry. Suppose `0<=j<k` and

```text
kappa_j^a=kappa_k^a             for every a in A,
O_j=O_k,
T_j=T_k.                                             (10)
```

Delete radix levels `j,...,k-1`. The reconstructed row is

```text
V'_i=R_(j,i)+q^j Z_(k,i).                             (11)
```

Then:

```text
0<V'_1<...<V'_d,                                      (12)
a.V'=0                         for every a in A.       (13)
```

### Relation proof and converse

For every relation `a.V=0`, equations (9) and (11) give the exact identity

```text
a.V'=q^j(kappa_j^a-kappa_k^a).                        (14)
```

Thus the carry equality in (10) is sufficient for (13). It is also
necessary: for a relation on the original row, the same relation survives
the block deletion if and only if its two boundary carries agree.

### Positivity

If `Z_(k,i)>0`, equation (11) is positive. If `Z_(k,i)=0`, owner equality
forces `Z_(j,i)=0`; hence

```text
V'_i=R_(j,i)=V_i>0.                                   (15)
```

### Strict order

It is enough to compare adjacent coordinates. If

```text
Z_(k,i+1)-Z_(k,i)>=1,
```

then

```text
V'_(i+1)-V'_i
 >=q^j-(q^j-1)
 =1.                                                  (16)
```

If instead the two `k`-quotients are equal, then `i` is absent from `T_k`.
Tie-mask equality makes it absent from `T_j`, so the two `j`-quotients are
also equal. The original strict inequality in (1) now says

```text
R_(j,i)<R_(j,i+1),
```

and (11) again gives `V'_i<V'_(i+1)`. This proves (12).

The equality `T_j=T_k` is load-bearing. Merely knowing that both quotient
vectors are sorted does not control a pair whose highest distinguishing
digits lie inside the deleted block.

## 3. Primitivity is not a missing LRC coordinate

Let `g=gcd(V'_1,...,V'_d)`. The normalized row

```text
V''=V'/g                                               (17)
```

is still a strictly increasing positive integer row and satisfies every
relation in `A`. Moreover multiplication by `g` on the time circle is
surjective and Haar preserving, so

```text
M(V'')=M(V'),          mu(S(V''))=mu(S(V')).          (18)
```

Thus a pump may change the common gcd, but normalization restores
primitivity without changing the lonely-runner target. It can change labelled
residues, so it should be performed only after any residue-dependent
certificate has been discharged.

Combining Sections 2 and 3 repairs all three algebraic losses listed in
THM-2167:

```text
relations:              preserved by equal carries;
positivity:             preserved by equal owners;
order/distinctness:     preserved by equal tie masks;
primitivity:            restored by target-invariant normalization. (19)
```

## 4. Algebraic finite terminal

Let `A={a^(1),...,a^(m)}` be nonzero integer vectors. Put

```text
C=product_(h=1)^m (2||a^(h)||_1-1).                  (20)
```

Suppose the simultaneous relation system

```text
a^(h).V=0,               h=1,...,m,                  (21)
```

has a strictly increasing positive integer solution of length `d`. Then, for
every radix `q>=2`, it has a primitive such solution with

```text
max_i V_i<q^(2dC).                                    (22)
```

Indeed, choose a solution of least maximum and normalize it. Let `J` be least
with `q^J>max_i V_i`. The boundary levels `0,...,J` carry at most `C` carry
tuples and, by (8), at most `2d` sidecar pairs. If

```text
J+1>2dC,                                              (23)
```

two full states repeat. The terminal empty owner occurs only at `J`, since
`O_j` is nonempty for `j<J`; hence the repeated owner is nonempty. Pumping
strictly lowers the largest coordinate while Sections 2--3 keep the row in
the same positive ordered relation class, contradicting minimality. Thus
`J+1<=2dC`, which implies the slightly weaker clean bound (22).

This is a uniform small-representative theorem for the **algebraic** relation
language. It does not bound every solution and, by Section 6 below, it cannot
replace an LRC row by its small representative while retaining the safe-set
predicate.

## 5. The LRC(14) rank-two carrier

For THM-2167's two relations `r,s`,

```text
||r||_infinity,||s||_infinity<=105,
||r||_1,||s||_1<=13*105=1365.                        (24)
```

Hence their carry pair has at most

```text
C<=(2*1365-1)^2
  =2729^2
  =7,447,441                                          (25)
```

values. The owner/tie pair assumes at most `2*13=26` values along the path,
so the effective ordered-path cap is

```text
26*2729^2=193,633,466.                                (26)
```

Whenever a longer path repeats a full state, the pump preserves the two
relations, positivity, strict order, and—after normalization—primitivity.
The adaptive prime from THM-2167 is unchanged because it belongs to the
relation pair, not to the speed magnitudes.

Consequently each feasible rank-two relation code produced by THM-2167 has
an explicitly bounded positive ordered representative. What is not proved is
that this representative remains zero-safe.

## 6. The remaining phase loss is genuine

The following thirteen-speed pump preserves every algebraic datum above but
changes the `1/14` safe-set measure. Take base `q=2`, levels `j=1,k=2`, and

```text
V=(1,4,5,8,9,12,13,16,17,20,21,24,25).              (27)
```

Its two boundary quotients are

```text
Z_1=(0,2,2,4,4,6,6,8,8,10,10,12,12),
Z_2=(0,1,1,2,2,3,3,4,4,5,5,6,6).                   (28)
```

Both have owner suffix `{2,...,13}` and cut mask

```text
{1,3,5,7,9,11}.                                      (29)
```

Let

```text
r=2e_2-e_4,             s=3e_2-e_6.                 (30)
```

They vanish on `V`, have zero carry at both levels, and remain independent
modulo `2`:

```text
r mod 2=e_4,             s mod 2=e_2+e_6.            (31)
```

Deleting level `1` gives

```text
V'=(1,2,...,13).                                      (32)
```

The same two relations vanish on `V'`; both rows are primitive, positive,
and strictly ordered.

Nevertheless, at `t=4/29` the row (27) has

```text
min_i ||V_i t||=3/29>1/14,                            (33)
```

so its safe set has positive measure. For (32), the fourteen points

```text
0,t,2t,...,13t
```

are pairwise at circle distance at least `1/14` exactly when they form the
regular fourteen-grid. Therefore

```text
S(V')={a/14:gcd(a,14)=1},                             (34)
```

a finite set of measure zero. The companion's independent lower-envelope
calculation sharpens (33) to

```text
M(V)=3/29,                 M(V')=1/14.                (35)
```

Thus the repaired finite state is still target-mixed. No carry, owner,
quotient-order, or primitivity argument can silently substitute for a phase
certificate.

## 7. Optional finite-residue sidecar

A fixed modulus can also be preserved at finite cost. Attach

```text
Q_j^(N)=Z_j mod N in (Z/NZ)^d.                        (36)
```

If this residue vector repeats together with the state in (10), then

```text
V-V'=q^j(Z_j-Z_k)=0 mod N                             (37)
```

coordinatewise. Hence any finite CRT bank can be made pump-invariant by
taking one common modulus `N`. Unlike `(O_j,T_j)`, this sidecar is not
monotone and can multiply the state count by `N^d`.

This does not contradict THM-2161 or THM-2163. Those theorems show that fixed
residue banks, even combined with bounded relation data and scalar magnitude,
still mix an adaptive modular or phase certificate. Equation (37) says only
that a *chosen finite* residue predicate can be retained.

## 8. Exact audit

The companion performs, with raising checks active under normal and optimized
Python:

1. `79,086` owner/tie monotonicity transitions over all sorted rows from a
   fixed exact universe, dimensions `2` through `5`, and bases `2` through
   `5`;
2. `5,554` equal-sidecar block deletions, all nontrivial and all preserving
   positivity and strict order;
3. `2,762` relation-splice implications and the same number of converses for
   every coefficient vector in `[-2,2]^4` on the selected exhaustive pump
   universe;
4. the exact witness (27)--(35), including mod-two independence; and
5. common-scale normalization controls.

The exact maximin routine enumerates every breakpoint of each circle-distance
function and every equality point of two affine branches. Normal and
optimized transcripts agree with the stored output. The computation is a
hostile control; Sections 1--7 are the all-row proof. QED.
