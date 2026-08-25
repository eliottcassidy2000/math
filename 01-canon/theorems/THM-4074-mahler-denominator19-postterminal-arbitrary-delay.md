---
id: THM-4074
title: "Mahler denominator-19 post-terminal arbitrary-delay obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.  The denominator-19
  terminal approximants contain a two-parameter family of reachable states
  after their last native one.  LTE forces arbitrarily long consecutive reset
  runs, and a normalized 2-adic exponential quotient is an isometry in the
  odd parameter.  Consequently every finite carry cylinder beginning in one,
  including every such safe cylinder and every greedy-boundary prefix, can be
  programmed after those resets.  Reachable post-terminal follower states are
  unbounded and no uniform finite rejection horizon exists.  This is a finite
  obstruction theorem: it neither produces a nonrejecting infinite orbit nor
  proves universal rejection, so Mahler's Z-number problem remains open.
source: codex-frontier-synthesis-creative-20260825c / Mahler post-terminal lane
audit: >
  PASS.  The primary carry-first replay checks 84 isometry layers for
  0<=s<=6 and 1<=h<=12, all 28,665 odd parameter/residue classes, all
  1,400 safe odd-start target realizations, 2,286 direct denominator-19
  prefix steps, 28 greedy programming controls, and even-parameter, unsafe
  11, and A=1 hostiles.  An independent native-first replay uses direct
  modular exponentiation and bitwise Hensel lifting on all 28,665 odd-start
  words, performs 573,468 lift-candidate gates, independently tests safety by
  every truncated suffix inequality, and directly follows 36,846 ordinary
  orbit steps through s=10.  Its t=1 controls all reject in the displayed
  finite box.  Both scripts use exact integers/Fraction only, contain no
  Python assert nodes or floating constants, and normal and optimized runs
  byte-match their frozen outputs.
depends_on:
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
  - THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation
  - THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction
script: 04-computation/mahler_denominator19_postterminal_arbitrary_delay_thm4074.py
output: 05-knowledge/results/mahler_denominator19_postterminal_arbitrary_delay_thm4074.out
script_sha256: 0611a42923d21f4156550dc35a832b31461976719410a942f7ea9f953de42a13
output_sha256: 42895b87b223222124f41c5b61a75a844086ed0466ab14ace4ce50b34eaafee8
independent_audit_script: 04-computation/mahler_denominator19_postterminal_arbitrary_delay_thm4074_independent_audit.py
independent_audit_output: 05-knowledge/results/mahler_denominator19_postterminal_arbitrary_delay_thm4074_independent_audit.out
independent_audit_script_sha256: b9c86ce801277b0533e2af6f850ba8434b8f00b95cbca2fca02ff31bb58e910d
independent_audit_output_sha256: 2a29e12293dd3d09c9e40a764a8231be8e0efe6932bdb991be17ed24ad0deda8
hash_basis: raw working-tree bytes (LF)
---

# THM-4074 -- denominator 19 defeats every bounded post-terminal horizon

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.**  This theorem
continues the live deterministic task isolated by THM-4072.  It proves that
the task has no uniform finite horizon and no bounded follower-state cap,
even after the last native `1`.  It does not decide any one
infinite orbit.

## 1. Exact terminal launch states

Write

```text
T(u)=ceil(3u/2),
```

and let `d=d_1d_2...` and `Q` be the greedy boundary and safe follower graph
of THM-4072.  In particular,

```text
d_1d_2d_3=101,
```

so the word `100` is a strict safe loop from follower state zero back to
state zero: its last `0` is a reset.

Fix integers

```text
s>=0,                 t>=1 odd,
m=18*2^s*t,
A_(s,t)=9(2^m-1)/19.                                  (1)
```

Both integers in (1) are well defined because the order of `2` modulo `19`
is `18`.  The pure carry period `(100)^infinity` represents the formal
2-adic state `-9/19`.  Since

```text
A_(s,t)=-9/19 mod 2^m,
```

THM-2228's cylinder bijection gives the first `m` ordinary carry digits of
`A_(s,t)` exactly:

```text
c_0...c_(m-1)=(100)^(m/3).                            (2)
```

Thus (2) is safe, its follower state at depth `m` is `q_m=0`, and its state
at depth `m-1` is `q_(m-1)=2`.

The native word is the ordinary binary expansion of `A_(s,t)`.  Since

```text
2^(m-2)<A_(s,t)<2^(m-1)                               (3)
```

for `m>=18`, its last native `1` is exactly at index `m-2`.  Hence the states
at depths `m-1` and `m` are genuinely post-terminal, not merely prefixes
that admit some terminal completion.

Let `u_j=T^j(A_(s,t))`.  Comparing (2) with the period `-9/19`, or directly
iterating the affine recurrence, gives

```text
u_m     =9(3^m-1)/19,
u_(m-1) =6(3^m-1)/19.                                 (4)
```

Equations (1)--(4) are the exact reachability sidecar missing from a search
over arbitrary pairs `(q,u)`.

## 2. LTE forces an unbounded reset runway

Because `m` is even, the two-adic lifting-the-exponent formula gives

```text
v_2(3^m-1)
 =v_2(3-1)+v_2(3+1)+v_2(m)-1
 =s+3.                                                 (5)
```

The other factors in `u_m` are odd.  Therefore

```text
v_2(u_m)=s+3,              v_2(u_(m-1))=s+4.           (6)
```

While an integer state is even, `T(u)=3u/2`; each step removes exactly one
factor of two and emits carry `0`.  Starting immediately after the last
native `1`, at `(q,u)=(2,u_(m-1))`, the next `s+4` carries are consequently
all zero.  The first zero resets because `d_3=1`; every later zero resets at
state zero because `d_1=1`.  Equivalently, from `(0,u_m)` there are `s+3`
consecutive reset edges before the integer state becomes odd.

This already proves arbitrarily long post-terminal nonrejection.  The next
section shows that the odd state after this runway is not confined to a thin
residue family.

## 3. The normalized exponential quotient programs every odd cylinder

Put

```text
k=s+3,
g_s=3^(18*2^s),
[t]_(g_s)=1+g_s+...+g_s^(t-1),
U_s=9*3^k(g_s-1)/(19*2^k).                             (7)
```

The integer `U_s` is odd: `19` divides `g_s-1`, while (5) with `t=1` says
that its exact two-adic valuation is `k`.  After the `k` zero carries from
the depth-`m` state, the resulting odd integer is

```text
W_s(t)
 =3^k u_m/2^k
 =U_s [t]_(g_s).                                      (8)
```

For distinct positive integers `t,t'`, suppose for definiteness that
`t>t'`.  Then

```text
[t]_(g_s)-[t']_(g_s)
 =g_s^t' (g_s^(t-t')-1)/(g_s-1).
```

Since `g_s=1 mod 4`, LTE again gives

```text
v_2(W_s(t)-W_s(t'))=v_2(t-t').                        (9)
```

Thus `t -> W_s(t)` is an isometry modulo every power of two.  In particular,
for every `h>=1` it permutes the residue classes modulo `2^h`, and it maps
odd parameter classes bijectively to odd integer classes.

THM-2228 proves independently that the first `h` carry digits of an integer
are a bijective coordinate on its residue modulo `2^h`.  Odd residues are
exactly the words whose first digit is `1`.  Composing this cylinder
bijection with (9) proves the main programming statement:

> **Post-terminal programming theorem.**  Fix `s>=0` and `h>=1`.  For every
> binary word `v=v_0...v_(h-1)` with `v_0=1`, there is a unique odd class
> `t mod 2^h` such that the reachable terminal launch (1), after its `s+4`
> immediate reset zeros, emits the word `v`.  A positive odd representative
> of that class supplies an actual positive integer launch.

If `v` is in the finite safe language `P_h`, this whole post-terminal segment

```text
0^(s+4) v                                               (10)
```

is nonrejecting.  If `v=d_1...d_h`, its `h` letters are all match edges and
the terminal orbit reaches follower state exactly `h`.

## 4. Exact consequences and boundary

The programming theorem has four all-scale consequences.

1. **No uniform post-terminal rejection horizon.**  For every `H`, choose
   `s` with `s+4>=H`.  The resulting reachable orbit takes at least `H`
   nonrejecting reset edges after its last native `1`.
2. **Reachable follower states are unbounded.**  For every `h`, program the
   greedy prefix `d_1...d_h`; after the reset runway the orbit reaches `q=h`.
3. **Every finite safe odd cylinder occurs.**  If `a_h=|P_h|` is the renewal
   count of THM-3848, exactly `a_h-a_(h-1)` words in `P_h` begin with `1`,
   and every one is realized after the terminal runway.  The first counts are

   ```text
   1,1,2,3,4,6,9,13,20,30,44,67       (1<=h<=12).     (11)
   ```

4. **Finite reset evidence is programmable too.**  Choosing prefixes of
   `(100)^infinity` after the runway supplies any prescribed finite number
   of additional resets while the native word is already permanently zero.

None of these statements passes to an infinite limit with one fixed positive
integer.  Both `t` and `A_(s,t)` generally change when `h` changes.  A single
orbit that never rejects and resets infinitely often would be a Mahler
candidate by THM-4072; this theorem neither exhibits nor excludes one.

## 5. Hostile ledger

The following controls prevent four tempting overstatements.

1. **Oddness is load-bearing.**  If `t` is even, (5) gains `v_2(t)` and the
   displayed zero count changes.  The fixed-`s` odd isometry cannot silently
   include even parameters.
2. **Programming is not safety.**  The map realizes every odd-start word,
   including `11`, which rejects at its second letter because `d_2=0`.
   Safety enters only when the target is explicitly restricted to `P_h`.
3. **Long finite survival is not an infinite survivor.**  In the exact
   `t=1`, `0<=s<=10` hostile box, every displayed ordinary orbit rejects;
   the rejection indices after the last native `1` are

   ```text
   10,6,10,8,15,10,14,12,19,14,21.                    (12)
   ```

   These are zero-based indices after the last native `1`.  The finite row is
   a control, not a universal rejection theorem.
4. **The rational pseudo-orbit is a different completion.**  The number
   `A_(s,t)+9/19=9*2^m/19` has the denominator-19 lower-half phase cycle
   through depth `m` and the hostile next phase `27/38`.  After depth `m`,
   however, its fixed rational tail is not the real tail determined by the
   ordinary carry orbit of `A_(s,t)`.  Its failure cannot be used to declare
   the deterministic integer orbit rejected.

The result therefore refutes bounded-horizon and bounded-state attacks on the
post-terminal problem, not the possibility of a genuinely global argument.

## 6. Exact replays

The primary carry-first replay is

```text
python3 -B 04-computation/mahler_denominator19_postterminal_arbitrary_delay_thm4074.py
python3 -B -O 04-computation/mahler_denominator19_postterminal_arbitrary_delay_thm4074.py
```

It constructs the geometric quotient in (7), exhausts every odd class in the
stated box, then checks direct launch states and hostile controls.  The
independent native-first replay is

```text
python3 -B 04-computation/mahler_denominator19_postterminal_arbitrary_delay_thm4074_independent_audit.py
python3 -B -O 04-computation/mahler_denominator19_postterminal_arbitrary_delay_thm4074_independent_audit.py
```

It instead computes `3^m mod 2^(h+k)`, divides off the exact two-adic
valuation, lifts `t` one bit at a time, tests safety from all truncated suffix
inequalities, and directly follows the ordinary native orbit.  Both outputs
byte-match the frozen artifacts named in the front matter.

## 7. Scope audit

- **THM-2228** supplies the integer/carry cylinder bijection and the exact
  ordinary orbit `T`.
- **THM-3848** supplies the finite safe language `P_h`, its renewal counts,
  and the greedy equality boundary.
- **THM-4072** supplies the follower/reset semantics and proves that a
  nonrejecting, infinitely resetting post-terminal orbit would be exactly the
  unresolved candidate mechanism.
- The new mechanism is the exact reachable family (1)--(4), the LTE runway
  (5)--(6), and the normalized exponential isometry (7)--(9).

Mahler's `3/2` Z-number problem remains open.
