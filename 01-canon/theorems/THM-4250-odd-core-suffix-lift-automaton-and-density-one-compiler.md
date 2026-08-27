---
id: THM-4250
title: "Odd-core suffix lift automaton and density-one compiler"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every odd
  core b>=3, the low-binary suffixes that close all THM-4244 complementary
  pairs form a finite carry language with an exact 0/1/2-child lift law.
  Their density among odd residues is monotone and tends to one
  exponentially. A negative collar generalizes the universal -1 suffix;
  b=3 has an exact Fibonacci complement, b=5 recovers THM-4243, and b=7
  gives new modulus-16 prime-power exponent families. Suffix closure is a
  sufficient compiler certificate, not an all-height necessity or factor
  converse.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4244-even-multiplier-odd-core-complementary-pair-factorial-compiler
related:
  - THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families
  - THM-4237-multiplier-six-binary-adjacency-prime-power-factorial-closure
  - THM-4243-multiplier-ten-double-overlap-prime-power-factorial-closure
primary_script: 04-computation/factorial_odd_core_suffix_lift_density_thm4250.py
primary_output: 05-knowledge/results/factorial_odd_core_suffix_lift_density_thm4250.out
independent_audit_script: 04-computation/factorial_odd_core_suffix_lift_density_independent_audit_thm4250.py
independent_audit_output: 05-knowledge/results/factorial_odd_core_suffix_lift_density_independent_audit_thm4250.out
primary_script_sha256: 32770e95c436fec4036833795fbd4bb28192a4f3bb2999cde92a18a96d16f99d
primary_output_sha256: 7006d1e0e4853ccc04c0c0fb317c47eda584a5f8969a6daa3767be03b5b58f79
independent_audit_script_sha256: a86aaef01071b058f4bd7143373950b5c31907ed5222d6e39bd9b159470b283f
independent_audit_output_sha256: 388d3de9cd9f3fb2a9dbf8739b9fc6524cbcc6dddea35bc79931fdcce6a33979
hash_basis: raw LF bytes
audit: >
  PASS. The independent path rederived the lift law from explicit occupied-bit
  sets, rebuilt the schoolbook carry machine, checked every odd 3<=b<=63
  through ell=13, and tested 253,921 lift cells. It confirmed the density bound,
  negative collar, Fibonacci specialization, and b=3,5,7 sequences. The only
  wording guard is that L is sharp for the universal individual-carry reset
  lemma, not asserted minimal for closure at every fixed b.
---

# THM-4250 -- odd-core suffix lift automaton and density-one compiler

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Definition, inheritance, and exact scope

Let `b>=3` be odd. For `ell>=1`, put

```text
M=2^ell,             [x]_ell = the least residue of x modulo M.       (1)
```

Define the odd-core suffix closure set

```text
R_b(ell)={r: 1<=r<M, r odd, and
          [sr]_ell & [(b-s)r]_ell !=0
          for every 1<=s<=(b-1)/2}.                                  (2)
```

Here `&` is bitwise intersection. The pure bit statements below hold for
every `ell>=1`. The THM-4244 universal-suffix setting normally takes

```text
2^(ell-1)>=b,                                                   (3)
```

which in particular makes `-1 mod M` a certificate.

The closest proved mechanism is THM-4244. If an even multiplier is written

```text
a=2^q b,                 q>=1, b odd,                              (4)
```

that theorem proves that its exact THM-3474 candidate barcode is empty iff

```text
(sH)&((b-s)H)!=0             for 1<=s<=(b-1)/2.                   (5)
```

If `H=r mod M` and `r in R_b(ell)`, the common low bit in (2) is also a
common bit of the full integer products in (5). Thus every member of
`R_b(ell)` is a sufficient suffix certificate for THM-4244.

This implication is one-way. A height can acquire its first common bit above
the retained suffix, so failure of (2) is not failure of (5). Likewise,
THM-4244 uses an empty candidate barcode only in the proved direction toward
rational coprimality. No surviving suffix or barcode constructs a factor.

## 2. Exact 0/1/2-child lift law

Fix an odd `r` with `1<=r<M`. For every `1<=u<b`, write

```text
ur = x_u + M q_u,
x_u=[ur]_ell,               h_u=q_u mod 2.                         (6)
```

The two odd lifts of `r` modulo `2M` are

```text
r_epsilon=r+epsilon M,             epsilon in {0,1}.              (7)
```

Their product residues satisfy the exact high-bit formula

```text
[u r_epsilon]_(ell+1)
 =x_u+M*(h_u xor (epsilon*(u mod 2))).                             (8)
```

Let the failed-pair set at level `ell` be

```text
F(r)={s:1<=s<=(b-1)/2 and x_s & x_(b-s)=0}.                       (9)
```

Because `b` is odd, each pair `{s,b-s}` has one even member `e_s` and one
odd member `o_s`. Then:

1. if `F(r)=empty`, both lifts lie in `R_b(ell+1)`;
2. if `F(r)!=empty`, a lift `epsilon` closes iff, for every `s in F(r)`,

   ```text
   h_(e_s)=1,                 h_(o_s) xor epsilon=1;              (10)
   ```

3. equivalently, a failed parent has one closing child precisely when all
   even-member high bits are one and all odd-member high bits agree; its
   unique closing lift is

   ```text
   epsilon=1 xor h_(o_s).                                        (11)
   ```

   In every other failed case it has no closing child.

Thus every parent has exactly zero, one, or two closing children.

### Proof

Equation (8) follows by adding `epsilon*u*M` to (6) and reducing modulo
`2M`. Its low `ell` bits are always `x_u`; its new bit `ell` is toggled by
the lift exactly when `u` is odd.

If a pair already overlaps below bit `ell`, both lifts retain that witness.
For a failed pair, the only possible new witness at the next level is bit
`ell`. The even multiplier's new bit is the fixed value `h_(e_s)`, while the
odd multiplier's new bit is `h_(o_s) xor epsilon`. They are simultaneously
one exactly under (10). One global lift must satisfy every failed pair, which
is possible iff their imposed values of `epsilon` agree. This proves all
three clauses. QED.

### Monotone density recurrence

Let `A_b(ell)` be the number of nonmembers of `R_b(ell)` satisfying the
compatibility condition in clause 3. The lift partition is disjoint, so

```text
|R_b(ell+1)|=2|R_b(ell)|+A_b(ell).                              (12)
```

With

```text
delta_b(ell)=|R_b(ell)|/2^(ell-1),                              (13)
```

equation (12) gives

```text
delta_b(ell+1)=delta_b(ell)+A_b(ell)/2^ell >=delta_b(ell).      (14)
```

This is exact monotonicity, not an empirical trend.

## 3. Finite low-to-high carry automaton

Write the suffix bits as

```text
r=sum_(i=0)^(ell-1) rho_i 2^i,             rho_0=1.             (15)
```

For each multiplier `1<=u<b`, maintain its schoolbook carry `c_u`. On input
bit `rho`, define

```text
y_u=(c_u+u rho) mod 2,
c_u'=floor((c_u+u rho)/2).                                      (16)
```

Starting from `c_u=0`, induction gives

```text
0<=c_u<=u-1.                                                    (17)
```

For every complementary pair maintain a persistent flag

```text
z_s'=z_s or (y_s=1 and y_(b-s)=1),          z_s initially 0.   (18)
```

After exactly `ell` inputs, the emitted word `(y_u)` is the low-`ell` binary
word of `[ur]_ell`. Therefore

```text
r in R_b(ell) iff every z_s is 1 after ell steps.               (19)
```

This is a deterministic finite automaton. A crude state bound is

```text
2^((b-1)/2) * product_(u=1)^(b-1) u
 =2^((b-1)/2)(b-1)!,                                           (20)
```

although the reachable automaton is much smaller in the first cases.

Let `T_0,T_1` be its two transition matrices, let `v` be the state vector
after the forced low bit `rho_0=1`, and let `lambda` select accepting states.
Then

```text
|R_b(ell)|=lambda^T (T_0+T_1)^(ell-1) v.                        (21)
```

Consequently the ordinary generating function

```text
sum_(ell>=1)|R_b(ell)|z^(ell-1)
 =lambda^T(I-z(T_0+T_1))^(-1)v                                (22)
```

is rational, and the count sequence satisfies an integer linear recurrence
with constant coefficients. This is an exact finite-state closed form; it
does not promise a short elementary expression uniformly in `b`.

## 4. A simultaneous carry reset and exponential density one

Put

```text
h=ceil(log_2(b-1)),                L=h+1.                       (23)
```

If the binary word of `r` contains `L` consecutive one-bits, then

```text
r in R_b(ell).                                                    (24)
```

### Proof

Fix `1<=u<b`, and inspect its carry immediately before the run. Put

```text
d=u-1-c_u,                 0<=d<=u-1.                           (25)
```

Under a one input, (16) gives

```text
d'=u-1-floor((2u-1-d)/2)=floor(d/2).                            (26)
```

Since

```text
d<=u-1<=b-2<2^h,                                               (27)
```

after the first `h=L-1` inputs of the run the deficit is zero. On the final
input the output is

```text
(u-1+u) mod 2=1.                                                (28)
```

This holds simultaneously for every `1<=u<b`. Hence all reduced products
have a one at the same final run position, closing every complementary pair.
QED.

The value `L` is sharp for this **universal individual-carry reset lemma**:
for `u=b-1`, a run of only `h` ones from incoming carry zero has, before its
last output, deficit
`floor((u-1)/2^(h-1))=1`, and therefore emits zero there. This does not assert that `L` is
the minimal closure-forcing run length for every particular odd core; pair
interactions can close earlier.

### Density bound

Among the `2^(ell-1)` odd length-`ell` words, bits `1,...,ell-1` are free.
Partition

```text
q=floor((ell-1)/L)                                              (29)
```

disjoint length-`L` blocks from those free positions and ignore leftovers.
Each block is all ones with probability `2^(-L)`, independently across the
chosen blocks. A nonclosing word can contain no such block, so

```text
|R_b(ell)^c|/2^(ell-1)
 <=(1-2^(-L))^floor((ell-1)/L).                                (30)
```

In particular,

```text
lim_(ell->infinity) |R_b(ell)|/2^(ell-1)=1                     (31)
```

with exponential convergence for every fixed odd `b`. Since low-bit closure
implies full-product closure, the actual all-bit overlap predicate also has
density one among odd heights. Equation (31) makes no density assertion about
the sparse subset of prime powers.

## 5. The negative collar

Let `w` be a positive odd integer satisfying

```text
w(b-1)<=2^(ell-1)=M/2.                                        (32)
```

Then

```text
M-w in R_b(ell).                                                (33)
```

Indeed, for every `1<=u<b`,

```text
[u(M-w)]_ell=M-uw in [M/2,M).                                 (34)
```

Every such residue has bit `ell-1` equal to one, so every complementary pair
overlaps there. If

```text
X=floor(M/(2(b-1))),                                           (35)
```

this supplies `(X+1)//2` certified odd residues. The `w=1` case is exactly
THM-4244's universal `-1` suffix under (3). The collar is sufficient and is
not a classification: certified and failed residues both occur immediately
outside it.

## 6. Exact `b=3` Fibonacci law

For `b=3` there is one pair, `{1,2}`. Modulo `2^ell`, multiplication by two
is a left shift with the top bit discarded. Thus

```text
[r]_ell & [2r]_ell !=0
```

iff the linear low-`ell` word of `r` contains adjacent one-bits. There is no
cyclic wrap from the top bit to bit zero.

Let `F_1=F_2=1`. Odd words whose low bit is one and which avoid adjacent ones
are counted by `F_ell`: after the forced initial `10`, the usual final-bit
split gives the Fibonacci recurrence. Hence, for every `ell>=1`,

```text
|R_3(ell)^c|=F_ell,
|R_3(ell)|=2^(ell-1)-F_ell.                                    (36)
```

This recovers THM-4237's adjacent-bit mechanism and upgrades its particular
prime congruence family to the exact suffix count at every length.

## 7. `b=5`: exact recovery of THM-4243

The two predicates are

```text
[r]_ell &[4r]_ell !=0,
[2r]_ell&[3r]_ell !=0.                                        (37)
```

At the least universal length `ell=4`, direct reduction gives

```text
R_5(4)={5,7,13,15}.                                           (38)
```

At `ell=6`, the exact set is

```text
R_5(6)={5,7,13,15,21,23,27,29,31,37,
        39,41,45,47,53,55,57,59,61,63}.                       (39)
```

This is byte-for-byte the twenty-class maximal modulus-64 suffix list in
THM-4243. Thus the general language recovers that theorem's two-coordinate
specialization without changing its factorial scope.

## 8. `b=7`: new modulus-16 and exponent families

The three reduced predicates are

```text
[r]_ell &[6r]_ell !=0,
[2r]_ell&[5r]_ell !=0,
[3r]_ell&[4r]_ell !=0.                                        (40)
```

At `ell=4`, the exact low-suffix set is

```text
R_7(4)={3,5,7,15}.                                            (41)
```

For these four residues the three overlap masks modulo 16 are respectively

```text
r= 3: (2,6,8),        r= 5: (4,8,4),
r= 7: (2,2,4),        r=15: (10,10,12).                      (42)
```

The other odd residues have at least one zero mask, so (41) is exact for the
declared low-16 observer.

Now take multiplier `a=14=2*7`, a prime `p>14`, and `H=p^k`. The odd units
modulo 16 have exponent dividing four. Equation (41) is equivalent to the
following sufficient prime-power table:

```text
k mod 4 =1:       p mod 16 in {3,5,7,15},
k mod 4 =3:       p mod 16 in {7,11,13,15}.                  (43)
```

No even exponent enters this particular modulus-16 suffix certificate.
THM-4244 therefore proves, under (43),

```text
gcd_Q(A_(14p^k-1)^(14p^k+1),A_(14p^k)^(14p^k+1))=1,         (44)
```

and for every exact quadratic `f=alpha+beta*x+gamma*x^2` with all three
coefficients nonzero, at least one of

```text
L(f^(14p^k-1)),       L(f^(14p^k)),       L(f^(14p^k+1))    (45)
```

is nonzero.

THM-4244 already supplied the `p=7 mod 8`, odd-`k` portion of (43). The new
fixed-prime exponent lanes are

```text
p=3 or 5 mod 16,       k=1 mod 4,
p=11 or 13 mod 16,     k=3 mod 4.                            (46)
```

Examples begin with `(p,k)=(19,1),(37,1),(43,3),(29,3)` and repeat every four
exponent levels for the same fixed prime.

## 9. Exact small-core counting recurrences

The reachable automata for `b=3,5,7` have respectively `4,20,29` states. The
primary checker constructs them and certifies the following recurrence
characteristic polynomials, which annihilate the count sequence:

```text
b=3: (x-2)(x^2-x-1),

b=5: (x-2)(x^2+1)(x^2-x-1)
     (x^4-x^3-1)(x^4-x^3-x^2+x-1),

b=7: (x-2)(x^3-x-1)(x^3-x^2-1)
     (x^5-2x^2-2x-1).                                      (47)
```

These are exact finite-matrix certificates, not guessed fits. For a proposed
degree-`d` recurrence on an `S`-state representation, the checker verifies
the first `S` discrepancy values. Each discrepancy is itself an `S`-state
matrix sequence, so Cayley--Hamilton propagates zero to all later indices.

## 10. Hostile controls and destroyed information

1. **Low suffix is not all-height necessity.** For `b=5`, residue `51` is
   not in `R_5(6)`. But

   ```text
   H=1331=51 mod 64,
   H&(4H)=1024,                 (2H)&(3H)=2560,             (48)
   ```

   and `1331 in R_5(11)`. Higher bits repair the failed short observer.

2. **There is no monotonicity in the odd core.** At modulus 16, `r=3` closes
   for `b=3` but not `b=5`, while `r=5` closes for `b=5` but not `b=3`.

3. **Negation and inversion are not symmetries.** At `b=3,ell=3`, `r=3`
   closes while `-3=5` does not. At `b=3,ell=4`, `13` closes while its inverse
   `5` does not.

4. **Every child multiplicity occurs.** At `b=3,ell=3`, parents `1,5,3`
   have respectively zero, one, and two closing children.

5. **A common witness bit is not the full structure.** At the least universal
   levels for odd `b<=63`, the minimum number of bit positions needed
   to hit every pair-overlap mask reaches five; one witness is
   `(b,ell,r)=(33,7,33)`. Negative collars and one-runs have a singleton
   witness, but general closure does not reduce to one common witness bit.

6. **The reset length has a narrow sharpness quantifier.** `L-1` fails as a
   carry-state-independent reset for multiplier `b-1`; it can still suffice
   for a particular core, state, or pair configuration.

7. **Candidate versus factor.** Membership in `R_b(ell)` is used only to
   empty THM-4244's candidate barcode. Nonmembership and nonempty barcodes do
   not construct factors, roots, or bad moment windows.

## 11. Verification universes

The primary companion uses integer modular products and a separately coded
digit/carry transducer. It checks:

```text
odd cores b=3,5,...,63;
20,460 direct/transducer and exact-lift cells over ell0(b)..ell0(b)+3;
2,909 negative-collar cells through ell0(b)+6;
21,328 incoming carry states for the simultaneous-run lemma;
the full least-universal-level atlas and four lift levels for every b<=63;
the b=3 Fibonacci identity through ell=22;
the b=3,5,7 finite-state recurrences and direct counts;
and exact density inequalities at six hostile odd cores.               (49)
```

The independent companion imports no primary code. It represents occupied
bits as explicit sets and separately implements schoolbook carry states. It
checks:

```text
every odd b=3,...,63 and every ell=1,...,13;
253,921 exact lift cells, with child histogram
  {0:53,371, 1:27,619, 2:172,931};
253,921 run-certificate cells and 403 exact density inequalities;
8,235 negative-collar residues;
b=3,5,7 sequences through ell=30;
and explicit zero-, one-, and two-child, short-run, and collar hostiles. (50)
```

Both scripts use explicit runtime gates rather than bare `assert`. Ordinary
and optimized runs are byte-identical to their frozen transcript.

## 12. Type and modular firewalls

- The pure theorem concerns binary suffixes of odd integer heights.
- The factorial transfer is exactly THM-4244's even-multiplier,
  prime-power-height, exact `{0,1,2}` support result.
- The moment statement is a three-slot slice of ambient `SFC(1)`, never
  `SFC(3)` or the full three-variable Factorial Conjecture (MISTAKE-350).
- No moment polynomial is reduced modulo a small prime here. The scripts use
  only integer bit arithmetic, so MISTAKE-363's factorial-denominator failure
  genus does not arise.
- Density one among all odd binary words is not promoted to density one among
  prime powers.

Reproduce from the repository root with

```bash
python3 -B 04-computation/factorial_odd_core_suffix_lift_density_thm4250.py
python3 -B -O 04-computation/factorial_odd_core_suffix_lift_density_thm4250.py
python3 -B 04-computation/factorial_odd_core_suffix_lift_density_independent_audit_thm4250.py
python3 -B -O 04-computation/factorial_odd_core_suffix_lift_density_independent_audit_thm4250.py
```

Each stream must byte-match its declared output.

**QED.**
