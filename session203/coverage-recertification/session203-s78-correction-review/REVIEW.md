# Independent review of the S78 endpoint correction

Verdict: **PASS in the stated finite scope; no correction requested.**
The reviewed input is the latest `COVERAGE_AUDIT-frozen.md` and its copied
`COVERAGE_MANIFEST-frozen.sha256`. No shared checkout was modified.

The repaired independent S78 implementation and supplied full-run receipt
support this finite certificate: every odd integer from1 through4,294,967,295,
other than7 and21, is attained by the physical F1/F2 families. The full pass
uses no F0 marks and reconstructs coverage from scratch. This review is not
a third independently authored full-range sieve and does not extend the range,
prove an eventual threshold, or infer universal absence from family misses.

## Endpoint and implementation review

For the legacy one-apex scalar, g(2^k)=2^k+1 for k>=1 follows from its
recurrence (or directly from the sum over zero bits). The old F0 loop includes
v=2^32; its split-table index v>>16 is65536, one past the table. The corrected
fallback to the direct formula makes that exact call safe. Its intended value
4,294,967,297 exceeds the cap, so it cannot legitimately mark a target. This
observation does not retroactively validate undefined behavior; the fresh
NOF0 full run supplies the needed certification.

The split formula is exact: separating the low16 zero-bit contributions gives
k Z(r)+Q(r), while the boundary01 transition gives precisely the subtraction
[k even][top bit of r zero]. The k=0 case remains the direct small table.
The patch changes only the fallback and deterministic endpoint controls.

In the published NOF0 run, both actual coverage families use v<=X-1, since
vmax=floor(X/(1+u))-1 or floor(X/R)-1 with R>=1. Thus no coverage evaluation
uses the problematic endpoint. The corrected self-check intentionally does.
The F1 domain is symmetric in its two one-apex indices, so u<=v is valid.
Since g>=1 and the endpoint indicator is at most1, F1 is at least(1+u)(1+v).
For F2, g h-epsilon e>=h-e>=0, so F2 is at least(1+v)R. These bounds justify
the displayed domains; all values actually marked have physical recipes.

The final scan visits every odd N from1 through X-1. With X=2^32 even,
there is no omitted odd endpoint. The allocated maps exceed every possible
N>>1 index. Concurrent writes are atomic stores of the same value1; implicit
OpenMP loop barriers precede reading. All loop counts and values fit signed
64-bit integers: F1 values are at most2X, and a deliberately loose F2 bound
is190X using the bank's h<=189, R>=1 and g(v)<=v+1. No near-overflow arithmetic
is used in the declared run.

## Independent controls in this review

`check.py` imports no coverage implementation. It validates all20 manifest
entries and independently recomputes:

- F1 pairs47,964,383,242 and work items72,787;
- F2 evaluations134,010,742,236 and work items33,303;
- 2,131,019 labeled tournaments in orders1 through7;
- all22 matching archived dyadic F1 miss counts;
- the2,732-triple bank bounds (maximum h,e,R equal189,45,432).

A separate recursive-g Python sieve reconstructs all odd targets through
65,536 with F1/F2 alone and finds exactly7,21 missing. Sixty-three adversarial
split/endpoint checks, including values on both sides of the old limit,
agree with a separately written bit-sum formula and recurrence. It also builds
256 literal joins of small one-apex tournaments, reverses the second block,
and verifies their Hamilton-path counts by subset DP against F1.

Normal and python -O outputs agree exactly. The supplied sanitizer records,
bank regeneration receipt, corrected source, full output and result JSON are
mutually consistent. I did not repeat the465-second full32 pass; the verdict
combines source/mathematical review, its supplied full execution record, and
the independent smaller tests described above.

The inherited general physical join theorem remains an explicit dependency.
No unrestricted spectrum claim or new all-order theorem is inferred.
