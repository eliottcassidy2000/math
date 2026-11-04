# Exact return versus reusable clock law

2026-11-03. **PROVED coefficient and valuation obstructions**, applied to
the **FINITE-EXACT** named witnesses in `literal/realized_depth3.json` and
the extended `literal/realized_suffix_depth3_moving1_all_prefixes.json`.
None of these obstructions rules out a different strong-suffix transition.

## The coefficient test is stronger than successful re-entry

Let a fixed marked strong suffix have (h,F). Its transitive-prefix family is
H_s=h+2^s F, s>=1. Suppose a fixed output suffix (h',F') and a fixed integer
prefix shift k satisfy

    h'+2^(s+k) F' = 3(h+2^s F)+c                    (C1)

for at least two distinct allowed values of s. Subtracting the two equations
and then substituting back gives, necessarily and sufficiently,

    h'=3h+c,   2^k F'=3F.                           (C2)

Consequently, if F is odd, no k>=1 works with an integral F'. For the
actual source h=2195,F=42153, a same-clock suffix return must have precisely

    (h',F')=(6587,126459) for c=2,
    (h',F')=(6589,126459) for c=4.                  (C3)

This does not forbid growing the total graph, moving old prefix vertices
into the suffix, allowing the suffix transformation to depend on s, or
retaining a richer nested grammar. It excludes a fixed suffix and fixed
positive clock increment as an all-size implementation. In particular,
zero-ready re-entry is a type check, not a proof of the recurrence (C1).

## Complete resynchronization of the three selected witnesses

Let s and t vary independently over all positive integers. The equation

    h'+F'2^t = 3(2195+42153*2^s)+c                  (C4)

has exactly the following solutions for the selected literal outputs.
Indices in this table refer to the original `realized_depth3.json`, whose
indexing is retained for the two-generation lineage:

|Receipt index|c|h'|F'|All positive solutions (s,t)|
|---|---:|---:|---:|---|
|0|2|5655|31673|(2,4)|
|1|2|39331|118273|(2,2)|
|4|4|46073|58294|(2,3)|

Proof for index0. The equation is
31673*2^t-126459*2^s=932, whose right side has 2-adic valuation 2.
Equal exponents make the left side negative. Unequal exponents therefore
force min(s,t)=2. If t=2<s, the left side is negative. If s=2<t,
substitution forces 2^t=16, giving t=4.

Proof for index1. Here
118273*2^t-126459*2^s=-32744, of valuation 3. Equal exponents give
-8186*2^s=-32744 and hence s=t=2. For unequal exponents min(s,t)=3.
If s=3<t, the left side is positive. If t=3<s, it is at most
8*118273-16*126459=-1077160, strictly below -32744.

Proof for index4. Write F'=2*29147. The equation becomes
29147*2^(t+1)-126459*2^s=-39484, of valuation 2.
If s=t+1, the left side is -97312*2^s, whose valuation is s+5>2.
Otherwise min(t+1,s)=2. The possibility t=1 would require
126459*2^s=156072, impossible for s>=1. The other possibility s=2
forces t=3. These are all cases.

Thus neither taking a different prefix length nor independently
resynchronizing these selected input/output clocks reuses the same suffix
pair. The observed return is exact at its original source size.

## A complete small Diophantine compiler

For positive integers a,b and nonzero integer d, every positive solution of

    a*2^t-b*2^s=d

can be found without a bounded exponent search. Write a=2^alpha*a0 and
b=2^beta*b0 with a0,b0 odd; set u=t+alpha,v=s+beta,e=v2(d).
If u!=v, the smaller exponent is e, giving two explicitly solvable cases.
If u=v, solve (a0-b0)*2^u=d. Thus at most three candidates need checking.
When d=0 the boundary is different: solutions exist exactly when a0=b0,
and then u=v is an unbounded line. The script treats that case separately.

`python 05-knowledge/results/zero_ready_strong_switch_20261103/check_clocks.py`
applies this exact integer compiler to all eight literal first-generation
witnesses in the extended all-prefix receipt. The two additional +4 suffixes
have (h',F')=(6377,7907) and (2985,7960); both have the unique positive clock
solution (s,t)=(2,6). The extended receipt reindexes the original six, so its
source filename is recorded explicitly in `clock_check.json`.
The script also checks the compiler against direct bounded arithmetic on
a separate small universe. The all-positive claim follows from the case
proof, not from the bounded control. No finite-field observer is used to
infer an integral transition or a physical tournament automorphism.
