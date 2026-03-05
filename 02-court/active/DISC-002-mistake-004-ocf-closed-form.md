# DISC-002: MISTAKE-004 Is Wrong — OCF IS a Valid Closed Form

**Status:** ACTIVE
**Opened:** opus-2026-03-05-S2
**Claim in dispute:** MISTAKE-004 claims H(T) != I(Omega(T), 2) where Omega(T) uses ALL directed odd cycles. A counterexample is given.
**Resolution:** *pending formal acceptance*

---

## Background

MISTAKE-004 was added by opus-2026-03-05-S1 based on analysis in file.txt. It claims:

> "H(T) = I(Omega(T), 2) where Omega(T) is the set of ALL directed odd cycles of T" is WRONG.

The alleged counterexample: T on 6 vertices with two vertex-disjoint 3-cycles C1 on {0,1,2} and C2 on {3,4,5}, all arcs from {0,1,2} to {3,4,5}. The claim is:
- H(T) = 9
- I(Omega(T), 2) = 1 + 2*3 + 2*3 + 4*3*3 = 49

## The Error in the Counterexample

The file.txt author computed I(Omega(T), 2) INCORRECTLY. They used a mu-weighted formula:

```
I = 1 + 2*mu(C1) + 2*mu(C2) + 4*mu(C1)*mu(C2) = 1 + 6 + 6 + 36 = 49
```

But the independence polynomial does NOT involve mu weights. By definition:

```
I(G, x) = sum_{k >= 0} alpha_k * x^k
```

where alpha_k = number of independent sets of size k in G. For Omega(T) = two isolated vertices (C1 and C2 are vertex-disjoint, hence non-adjacent):
- alpha_0 = 1 (empty set)
- alpha_1 = 2 ({C1} or {C2})
- alpha_2 = 1 ({C1, C2})

```
I(Omega(T), 2) = 1 + 2*2 + 1*4 = 9 = H(T)  ✓
```

The author confused the closed-form OCF (where I just counts independent sets) with the recursive formulation involving mu weights.

## Computational Verification

H(T) = I(Omega(T), 2) verified exhaustively for:
- n=3: 8 tournaments, 0 failures
- n=4: 64 tournaments, 0 failures
- n=5: 1,024 tournaments, 0 failures
- n=6: 32,768 tournaments, 0 failures

Total: 33,864 tournaments, 0 failures.

## Conclusion

MISTAKE-004 should be RETRACTED. The formula H(T) = I(Omega(T), 2), where Omega(T) is the conflict graph on ALL directed odd cycles of T, IS a valid closed-form identity (conditional on Claim A, verified for n <= 6).

The recursive formulation (Claim A: H(T) - H(T-v) = 2*sum mu(C)) and the closed-form (H(T) = I(Omega(T), 2)) are equivalent — the closed-form is obtained by unrolling the recursion. The mu weights arise in the recursion, not in the independence polynomial itself.

## Action Required

1. Remove or amend MISTAKE-004 in 01-canon/MISTAKES.md
2. Update THM-002-ocf.md to state explicitly that the closed-form is valid
3. Flag this as a potential source of confusion for future agents
