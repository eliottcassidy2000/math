# A physical reset that repairs mismatched paired seeds

The [minimum-seed theorem](PAIRED_AXIOMS.md) says that simultaneously scaling every rooted coordinate from h=5 to h=35 traps a homogeneous recursion in multiples of 5. This construction changes the assumption: preserve only two coordinates on entry, then use a literal attachment to synchronize the entire state. It escapes the seed obstruction rather than contradicting it.

Let a root a have exactly one outneighbor v. Put h=H(T), e=number of full paths ending at a, and eta=number of full paths in T−a ending at v. Let t denote the distinguished-start optional two-slot read from the earlier recursion, and g=H(T−a). Then
```
t=h+eta,        g=e+eta.
```
For the first formula concatenate the other path before the distinguished a-starting path. Every Hamiltonian path splits uniquely this way. The only failures are when the optional path ends at v: then the distinguished path is the singleton a, giving eta failures. For the second formula, a path in the deletion extends by a at the end unless its endpoint is v. In particular **eta is a marked endpoint count**, not a tournament Hamiltonian total; it may equal 7 even though g and h cannot.

Ordinally attach TT_n, n≥3, with sink b and penultimate vertex s. Reverse only the cross arc b→a relative to the ordinal join, then move the root to s. This has sole outneighbor b. The exact map is
```
(h,e,eta) → ( h+2^(n−1)(h+eta), 2^(n−2)(h+eta), h ).
```
The ordinary result follows by the one-arc seam: no old path is lost because b is not a full starting vertex of TT_n, and the distinguished end read at b is 2^(n−1). A path ending s must use b→a. Its path through A is an optional prefix and a distinguished suffix starting at a, and the n−2 earlier tail vertices can each be before b or after the returned suffix and before s. This gives 2^(n−2)t. After deleting s, paths ending b must concatenate A and the transitive remainder, giving eta'=h. This is valid with n≥3, so deleting s leaves at least two tail vertices.

The map ignores e. Hence scaling just (h,eta) by 7 synchronizes all three new readings, and also t' and g' by their displayed sums. It can be used before an arbitrary sequence of the paired transfers. Given h,eta and a prescribed target N, it chooses n by testing that (N−h)/(h+eta) is a power of two at least 4. Conversely on the output slice eta'=h, the target profile gives the smaller h directly. The reset lane cannot itself give 7 or 21: one has N≥5h and N≡h mod4, and for h=1 the only possible multiplier h+eta is 1 or 2. These are statements about this physical lane, not independent all-tournament exclusions.

## Literal entries

Root all four input graphs at 0:

| Lower input | Upper input | Lower reads (h,e,eta) | Upper reads | Exact reset outputs, n=3 |
|---|---|---|---|---|
| (2,4,9,3,15,31,63,127) | (2,4,33,23,39,75,31,127) | (5,2,1) | (35,12,7) | 29 and 203 |
| (16,5,9,19,38,15,63,127) | (16,33,11,19,102,77,143,63) | (29,8,1) | (203,64,7) | 149 and 1,043 |

For the first row, the deletion totals on entry are 3 and 19, not 3 and 21. This is exactly how the reset avoids demanding the forbidden 21 as a physical deletion while making use of its factor-seven endpoint budget. The output three-vectors are (29,12,5) and (203,84,35), giving primitive full states with d=17 and 119. The same input pair and varying n realize 5+6·2^(n−1) and seven times it. The second row gives (149,60,29) and (1043,420,203). All displayed prime claims and full states are checked by direct counting. The replay also checks 330 small literal resets.

This establishes that a paired recursion need not start in perfect lockstep; there are physical ways to assign different auxiliary states, then bring the continuations into exact agreement. Remaining proof task: choose physical input profiles (h,eta) that hit a prescribed target or combine this reset with the other descent lanes. It is not enough that an arithmetic pair fits the power test.

Replay: `python3 sole_exit_reset.py`.
