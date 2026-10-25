# Independent audit: exact 27-word departure and native refill

2026-10-25. **PASS in the stated shortcut-map convention and finite-prefix
scope.** Both the exact first-bit departure and the all-h native refill
construction are correct. Neither supplies a global Collatz rank, new
descent coverage, a stopping-time bound, or convergence.

Reviewed drafts:

- `collatz_native_refill.md`, SHA-256
  `13595dab7015b4f43b51de86a0ae53de0a73a66bf97365a9b82466e8b0a2ae83`;
- `anchors_collatz_phase.md`, SHA-256
  `80dbc4ff9ec76e204d360de1b1eb2e7ad04242bcc359170f1aac776a37d0040c`.

## First departure

The shortcut is n/2 for even n and (3n+1)/2 for odd n. Independently
iterating n=27 for 59 shortcuts gives endpoint 23 and 37 odd bits. Deriving
the affine numerator alongside that iteration gives exactly

    q=2^59, a=3^37, b=23q-27a, d=q-a.

The comparison rho=b/d has odd denominator, lies strictly between 8 and 9,
and directly replays the same word before returning to itself. Rational
sources with odd denominators have the same parity exactly when their
difference has positive 2-adic valuation. On a common parity branch the
difference is multiplied by 1/2 or 3/2, lowering its valuation by one.
Therefore ell=v2(dn-b) is exactly the first mismatch position, with the
comparison orbit advanced at each step. One must not replace that moving
comparison by the fixed rho at every intermediate step; the draft does not
make that error.

The ordinary clock counts one operation per common even shortcut and two per
common odd shortcut. This proves the stated clock formula. Direct independent
replay verifies n=27 has ell=61, endpoint 53, ordinary clock 100. The first
59-bit block alone ends at 23 at clock 96; its exhaustion count K=1 omits
the next two common bits. The draft distinguishes these stopping points.

## Arbitrary prescribed exact refill

Let c=23d-2b=31a-23q. For h>=1 set M=2^(59h+2) and choose

    t_h=(2^(59h+1)-c)*(ad)^(-1) mod M, 0<=t_h<M.

The inverse exists, and the residue is odd and nonzero because c and ad are
odd. Hence t=t_h+Mu with u>=0 yields infinitely many positive odd integers.
For n=27+qt, the actual guarded word gives m=23+at, which is even. Its next
native step is z=m/2. The three exact identities are

    dn-b=q(4+dt),    dm-b odd,    dz-b=(c+adt)/2.

The first parenthesis is odd; the construction makes
c+adt congruent to 2^(59h+1) modulo 2^(59h+2), giving **exact** valuation
59h+1. Thus ell at shortcut times 0,59,60 is (59,0,59h), and K is
(1,0,h). The corresponding ordinary times are 0,96,97. There is no missing
CRT compatibility or positivity case.

At z, the periodic-word isometry gives precisely h further complete words.
Their endpoint is at shortcut time 60+59h, ordinary time 97+96h, and has
ell=0. No subsequent path or eventual outcome is forced. The positivity and
original-source comparison are retained: n-m=4+dt>0 and m>z>0. These sources
already belonged to the inherited first-block descent cylinder, so this
result adds exact transition information rather than global descent coverage.

The additional t=0 mod7 restriction is compatible because M is a power of
two. It leaves infinitely many lifts and fixes residues (n,m,z)=(6,2,1)
modulo seven for every h. In particular the same current record (ell(m),
m mod7)=(0,2) has next depth 59h for arbitrary h. This refutes sufficiency
of that specific compressed state, not all finite descriptions or all
possible arithmetic sidecars.

## Independent supporting checks

`collatz_independent_audit.py` is standard-library-only, imports no producer
code, and derives the word/carry from native iteration rather than copying
the declared word. Its SHA-256 is
`834f6a0f78defb0ed47f4d6818fbc2f58de9a0302f2aa33e67b9d661ebc9ea36`.

The narrow universe is the rational 59-step cycle, the n=27 departure,
h=1 and h=3 with two positive lifts each, exact suffix exhaustion, and a
mod-seven-matched lift at each h. It passes 917 per-step exact checks. A
separate ordinary-map replay checks clock 100 directly. Normal and optimized
Python output is byte-identical; every gate raises an explicit exception.

The unrestricted quantifier over h is proved by the residue progression,
not inferred from these finite fixtures. No substantive issue was found.
