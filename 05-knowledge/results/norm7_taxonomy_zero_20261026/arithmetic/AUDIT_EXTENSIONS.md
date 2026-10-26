# Independent audit of the prime-power cone transport

2026-10-26. Reviewer: arithmetic_shapes. **Verdict: PASS**, within the
stated raw character-sum scope. No mathematical correction required.
The known relative-link edit is editorial and does not change this verdict.

Reviewed `../motives/EXTENSION_TWIST.md` and `../motives/check_extensions.py`.
Reviewed proof SHA256:
`1e63d3288044edbc389ed98879cddda33f59b3da5125255a4cadd5e819551134`;
checker SHA256:
`cc072719865576bf4b4e0b31b9b989b4c5e4099d7dd06e4ec40d021cb822eae0`.
The full producer replay passes **359 gates** in both normal and optimized
Python, with byte-identical output. A separate standard-library implementation
passes **363 controls**, also byte-identical across modes. It imports no
producer code, uses different polynomial models for each extension field,
and sums literal oriented cubic products rather than using E1 as its oracle.

## Signs and the all-field proof

For each cone edge directed from the cone label toward a larger label,
chi(neighbor-x)=epsilon*chi(x-neighbor). Thus the total correction is
epsilon^k, and the affine change x=a+(b-a)y contributes chi((b-a)^3)=chi(b-a).
This proves E1 and explains precisely why replacing epsilon by -1 cannot
be transported to even extensions of a prime3mod4 field.

The pinned edge lists give the following direct factors. The common cone
factor from neighbors (0,s,t) is chi(s)L(t/s); its orientation contributes
epsilon for H02 and H05 only. The other cones give L(s)^2 for H02/H05,
and L(s)L(t) for H09/H18. Multiplying their base edges gives exactly the
four E2 weights. In H04 the direct weight is
chi(s)chi(t)chi(t-1)chi(t-s)L(t)L(s)^2; in H06 there is the additional
factor chi(s-1). These give the displayed A/B expressions without imposing
epsilon=-1.

The two Legendre identities have the stated signs. For T(u)=u/(u-1),
substituting s=T(v) gives

    T(u)-T(v)=-(u-v)/((u-1)(v-1)).

Consequently both A(T(u)) and B(T(u)) equal epsilon*chi(u-1) times their
original value. The extra factor in B transforms as
chi(T(v)-1)=chi(v-1), so it does not introduce a new epsilon. The same
factor transforms L, making L*A and L*B invariant. On chi(u-1)=-1,
T flips chi(u), and its only allowed fixed point u=2 has chi(u-1)=1.
This validates the cancellation proof of E3 over every odd finite field.
For C05 the two epsilon factors cancel after t=su and s->1/s; this is
the easily missed sign in its transformed B expression.

## Independent literal values

The extension moduli, with coefficients listed from constant term upward,
are [2,1,1] over F3, [1,2,0,1] over F3, and [3,1,1] over F7. They differ
from the producer's chosen models and give fields of sizes9,27,49.

| q | C02 | C04 | C05 | C06 | C09 | C18 |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 16 | 16 | -16 | -16 | -16 | 16 |
| 7 | -128 | 128 | 128 | 128 | 0 | 0 |
| 9 | -208 | -208 | -16 | -16 | 48 | -16 |
| 27 | -8064 | 8064 | 12672 | 12672 | 11520 | 4992 |
| 49 | 2928 | 2928 | 11760 | 11760 | 1904 | -2320 |

All values agree with the producer. The independent controls also verify
the inverse and T Legendre identities and both A/B transformation laws at
every allowed parameter in these fields. Source graph JSON SHA256:
`3036a2093792548afd6113344b46b46f515f5a823cfcf531284bf7d69d58dab9`.

## Type and consequence audit

These are cone sums, allowing independent cone coordinates. They are not
fully injective seven-vertex graph counts. Epsilon=1 makes the quadratic
relation symmetric, so a Paley tournament census is undefined there;
returning null rather than zero is correct. The transitive-set recursion
counts each subset by its unique source when the tournament is defined.

The base-field Legendre recurrence L_(p^2)=2p-L_p^2 has the correct trace
sign. The illustrative Newton expressions6728 and9656 are arithmetically
correct, but neither is presented as a coefficient of a proved pure motive.
Rank, purity, integral model and conductor remain separate obligations.

Reproduce the independent receipt from this directory:

    python3 -B audit_extension_independent.py
    python3 -B -O audit_extension_independent.py

`extension_audit.normal.json` and `extension_audit.optimized.json` have SHA256
`e9cae70c6cfcaf98ab32cc123b3f8a758e11216929b27656432805b2f43e52ce`.
