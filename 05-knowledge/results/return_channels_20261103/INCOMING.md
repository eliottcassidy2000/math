# Incoming integration: universal parity and a uniform 3-adic level-raising loop

Read incoming work through `a140b314340c` after the first local checkpoint.
The broader incoming censuses were not rerun. New identities and literal
graph checks here are reproduced by `integration_probe.py/json`.

## I1. The parity barrier survives enlargement of the atom bank

The incoming [ternary carry theorem C1](../ternary_repair_clock_20261103/CARRIES.md)
proves the index-two correction lattice invariant for every supplied
single-bridge block. In our processing-order convention, its universal
calibration is

    row0(P_B)=v(B),
    row2(P_B)=v(a B)-v(B),
    row1(P_B)=v(a a B)-v(a B).

All entries are integers. Each first-coordinate count is the odd Hamilton
count of an actual nonempty tournament, so the first column of P_B is
(1,0,0) modulo 2. Applying THEORY T3's primitive-content argument gives:
**arbitrary sequences of supplied blocks of this join type cannot send an
odd-first primitive error to an even-first primitive error.** This holds
at every block size. A larger atom bank alone cannot return r to a multiple
of t. Re-marking and source replacement remain different operations.

This extension is a corollary of the incoming calibration and our parity
argument; a finite bank scan was not its universal proof. Neither this
barrier nor the incoming integer correction module classifies actual sources.

## I2. A uniform level-raising loop in the original 65 alphabet

The incoming [return-bank theorem R](../adaptive_state_induction_20261103/RETURN_BANKS.md)
uses the saturated bulk quotient of t=(2,-1,0):

    beta=(y,h+2x-y)=vB,     B=[[0,1],[0,2],[1,-1]].

When every return quotient has determinant prime to 3, its 3-adic bulk level
is unchanged. A ternary replacement raises that level, so a finite fixed-seed
bank with those returns cannot close. That theorem is inherited.

Our loops leave its unit-determinant hypothesis explicitly. Put D=L2 b L2,
where L2=`r24 a a a`. This literal processing word costs 19 vertices, and

    t P_D=-4t,
    P_D B=B Q,
    Q=[[715050,1009641],[268050,378483]]=3R,
    R=[[238350,336547],[89350,126161]], det R=-100.    (I1)

For every nonzero integer row beta, let level3(beta) be its common 3-adic
valuation. Then

    level3(beta Q)=level3(beta)+1.                    (I2)

Divide beta by its largest common power of 3. Its nonzero reduction remains
nonzero after R, since det R is a 3-unit; the factor 3 then adds exactly one
level. This proves the identity for every input, not just the tested seeds.
D returns to the SAME even-first error line. It does not reverse I1's
independent parity barrier.

This removes one hypothesis of the inherited obstruction; it does not
negate that theorem or close a source bank. In particular beta Q=3 beta R
is not beta Q=3 beta: the normalized bulk still changes substantially.

There is a sharper integral obstruction. The row lattice Z^2 R has index
100 and is exactly

    (a,b) in Z^2 R iff a+50b=0 mod 100.               (I2a)

Both rows satisfy the congruence; its kernel has index 100, so inclusion
and equal index give equality. The Smith factors are 1,100. Therefore an
exact return `w P_D=3v+m t` from ANY integer state w requires (I2a) for
beta(v), since beta(w)R=beta(v). The U,U2,U4 bulks have residues 25,75,75,
respectively, and fail this gate. Thus this particular exact reset is
excluded even before testing source realizability. Mixed return words,
additional carries and other input seeds are not excluded.

The exact characteristic polynomial is

    char(P_D)=(X+4)(X^2-1093533X-900).

Thus N_k=H(U D^k) obeys

    N_(k+3)=1093529 N_(k+2)+4375032 N_(k+1)+3600 N_k,
    (N0,N1,N2,N3)=
      (31,33868343,37036149733639,40500251957161045607).

For every k>=0 the actual companion seeds give

    H(U2 D^k)=3N_k+2(-4)^k,
    H(U4 D^k)=3N_k+4(-4)^k.                           (I3)

Their vertex clocks are 7+19k, 8+19k, 9+19k. Words and newest ports are
retained. Literal adjacency frontier counting verifies all three families
for k=0,1,2,3, reaching 66 vertices. The bulk levels are k for U and k+1
for U2,U4. No primality or target-coverage conclusion is inferred.

## I3. A separate eight-phase normalized clock

After removing the forced factor 3 at each step,

    R mod 3=[[0,1],[1,2]].

Let eta^2=2eta+1 over F3. This polynomial is irreducible and eta has order
eight in F9*. Encoding (u,v) as u+v eta turns row multiplication by R into
multiplication by eta. Every nonzero normalized bulk therefore has exact
period eight modulo 3. Valuation and unit phase are separate coordinates.

For U, beta0=(25,34) reduces to 1+eta=eta^7. Its relative clock is
alpha_k=eta^(k+7). Frobenius z->z^3 consequently acts by

    k -> 3k+6 mod 8.                                  (I4)

The source phase 6 is essential. The probe checks the field law, all eight
nonzero phases, and the actual source normalization. This is distinct from
the incoming F49 clock and our negative-offset family's 13-adic clock.
The error magnitude is 4^k: lambda is always +1, mu vanishes for k>=1,
and error orientation is (-1)^k. These arithmetic observations forget the
eight-phase bulk even while the full source word retains it.

## I4. Incoming operations sharpen the next test

- [Common two-port observer C1-C3](../adaptive_state_induction_20261103/COMMON_OBSERVER.md)
  supplies a universal 17-query future-H carrier for twelve named birth/discard
  commands. Its port replacement is richer than our one-port append language.
  Dimension 17 and our three integers plus a feasibility bit concern different
  operations and outputs; neither replaces the other without a map.
- [Ternary repairs R3-R5](../ternary_repair_clock_20261103/THEORY.md) supplies
  30 branches through depth four with a larger atom bank and explicit carries.
  The four ideal affine states still fail. These give useful literal endpoints
  for return tests, not an all-size closed induction.
- [Session154 mixed-prefix certificate](../session154_gadgets_20261103/full_prefix/NOTE.md)
  reports the complete primitive-pair orbit modulo 7*13^2 with parity and
  actual prefix order retained. It removes that modular coverage gap. Its
  newly recovered example has composite quotient, so it does not give an
  integer prime/7p inverse or permit discarding the carry.

The next decisive test is to try D or another nonunit return on the supplied
repaired endpoints while retaining normalized bulk phase, integer carry,
port feasibility, and actual graph. Prove a return to a declared source
family. A correct valuation change or a full residue image is insufficient.

**Final incoming update, 2aa09af4b8:** the concurrent
[two-step s78 law](../adaptive_state_induction_20261103/INCOMING.md), I1-I2,
proves an analogous level increase using the enlarged alphabet: Q_s78^2=3U
with det U=49, followed by an exact mod-49 return gate. That update arrived
after our first integration checkpoint. It prompted the independent
mod-100 calculation (I2a) above. The two constructions use different words
and lattices; both demonstrate why removing the valuation obstruction
exposes a remaining integral-image condition rather than source closure.
