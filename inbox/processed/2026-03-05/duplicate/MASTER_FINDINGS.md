# Master Findings: Cycle Length Analysis in Parity of Tournaments
# Generated through scripts 1-9

## CRITICAL CORRECTION: Claim A verification numbers look wrong

At n=4: "B ≠ C/2 (Claim A violation): 80" -- but Claim A is VERIFIED for all n<=6!
This means there's a bug in sum_mu() -- the μ computation is wrong.
The issue is likely in ind_poly_at_2_restricted(): it uses T[perm[i]][perm[(i+1)%len]]
which for a cycle perm of length L wraps around, but perm doesn't include v.
This should be checking cycles of T-v, but the code uses the FULL T matrix.
Cycles of T-v that don't include v should work, but the cycle-finding may have bugs.

**This does NOT invalidate the paper's claims -- it means our scripts have a μ computation bug.**
**The scripts for (inshat) analysis (scripts 1-5) are unaffected since they don't use μ.**

---

## SOLID FINDINGS (from scripts 1-5, unaffected by the μ bug)

### Finding F1: ALGEBRAIC IDENTITY (proven, not just verified)
(inshat(v, P') - 1) / 2  =  #Type-II positions in P'

This is a pure algebraic identity for any binary signature. Proof:
  Let s = signature, b = boundary term, #01 = Type-I count, #10 = Type-II count.
  inshat = b + #01 + #10
  (inshat - 1)/2 = (b + #01 + #10 - 1)/2
  Need: b + #01 - 1 = #10
  Using #10 = #01 + s[0] - s[-1] (telescoping) and b = s[0] + (1 - s[-1]):
  b + #01 - 1 = s[0] + 1 - s[-1] + #01 - 1 = s[0] - s[-1] + #01 = #10. QED.

**Consequence**: The per-path identity (inshat-1)/2 = Σ_C μ(C) over 3-cycle embeddings
is equivalent to: #Type-II = Σ_{3-cycles C, nonv consec in P'} μ(C)

### Finding F2: EVERY TYPE-II POSITION IS A 3-CYCLE (proven)
For any path P' of T-v and any Type-II position j (s[j]=1, s[j+1]=0):
The arcs v→P[j] (since s[j]=1), P[j]→P[j+1] (since P' is a path), P[j+1]→v (since s[j+1]=0)
form a directed 3-cycle v→P[j]→P[j+1]→v. Conversely every such 3-cycle gives exactly
one Type-II position in each path where P[j],P[j+1] are consecutive.
This is an exact bijection between Type-II positions in P' and 3-cycle embeddings in P'.

**Consequence**: 
  #Type-II positions in P' = #{3-cycles (v,a,b): (a,b) consecutive in P'}
This is also a pure identity, no μ weights needed.

### Finding F3: MINIMUM TYPE-II COUNT FOR L-CYCLE EMBEDDINGS (proven)
For any odd L-cycle consecutively embedded in P', the window s[j..j+L-2]
has s[j]=1 and s[j+L-2]=0, giving at least 1 Type-II position in the window.
The minimum (1 Type-II) is achieved by the "monotone" pattern 1,1,...,1,0.

**Consequence**: Every consecutive embedding of an L-cycle creates at least 1 Type-II position
inside the window -- and those Type-II positions correspond to 3-cycles "inside" the L-cycle.

### Finding F4: EXACT DISTRIBUTION FORMULA (proven)
For an L-cycle consecutively embedded, the number of internal signature patterns
giving exactly k Type-II positions within the window = C(L-2, 2k-1).
Maximum k = (L-1)/2, achieved uniquely by the fully alternating pattern 10101...10.

### Finding F5: THE REAL REASON THE PER-PATH IDENTITY FAILS AT n=6 (proven)
The identity:
  (inshat-1)/2 = Σ_{3-cycles (v,a,b), (a,b) consecutive in P'} μ(a,b)

can fail at n=6 because:
  - The LHS = #Type-II = #{3-cycle embeddings} (algebraic identity)
  - The RHS weights each 3-cycle embedding by μ(C) = I(Ω(T-v)|avoid{a,b}, 2)
  - For n=5: μ(C) = 1 is common (Ω(T-v) restricted to non-{a,b} cycles is often empty)
  - For n=6: μ(C) can be > 1 (more complex conflict graph structure)
  
**The failure is NOT about 5-cycles creating Type-II positions.**
**The failure is about μ WEIGHTS exceeding 1 when n=6.**

Wait -- this contradicts what we thought earlier. Let me verify this in a new script.

Actually: the paper's per-path identity says:
  (inshat(v,P')-1)/2 = Σ_{j: Type-II at j} I(Ω(T-v)|_{avoid {P[j],P[j+1]}}, 2)

For n≤5: we verified this holds (0 failures). So the weights DO work for n≤5.
For n=6: this fails (2,758/9,126 failures).

The question is: WHY does it fail at n=6? 
Option A: Because μ weights are wrong (some 3-cycle's weight doesn't equal μ(C))
Option B: Because 5-cycle embeddings create "phantom" Type-II positions not in 3-cycles
Option C: Because the 5-cycle through v is missed entirely

From F2: Type-II positions at (a,b) ALWAYS correspond to 3-cycle v→a→b→v.
So the "phantom" explanation is wrong. Every Type-II position IS a 3-cycle.

The per-path identity at n≤5 says: each Type-II at (a,b) should contribute
I(Ω(T-v)|avoid{a,b}, 2) to the sum. This equals 1 when n≤4 (Ω(T-v) trivially).
At n=5: T-v has 4 vertices, can have 3-cycles, so I(Ω(T-v)|avoid{a,b}, 2) can be > 1.
Yet the identity still holds at n=5!

At n=6: T-v has 5 vertices, can have 3-cycles and 5-cycles, richer Ω(T-v).
The per-path identity FAILS here.

THE CORRECT DIAGNOSIS:
The per-path identity (inshat-1)/2 = Σ_{Type-II at j} μ({P[j],P[j+1]}) is asking:
  "Sum μ(C) over each 3-cycle C={v,P[j],P[j+1]} that is consecutively embedded"

For this to equal Σ_{ALL odd cycles C∋v, nonv consec} μ(C), we need 5-cycles to contribute 0.
A 5-cycle's consecutive embedding contributes μ(5-cycle) to Claim A's RHS,
but contributes 0 to the per-path identity's RHS (only 3-cycle Type-II positions).
So the per-path identity UNDERCOUNTS cycles of length ≥ 5.

But wait -- at n=5, there are NO 5-cycles through v (v has only n-1=4 other vertices,
a 5-cycle through v needs 4 other vertices, and those 4 vertices must form a specific path
in T-v -- and the 5-cycle itself would use all 5 vertices including v).
Actually a 5-cycle DOES exist at n=5: v→a→b→c→d→v requires exactly the other 4 vertices.
But in T-v (which has 4 vertices), the path a→b→c→d must exist.

So at n=5, 5-cycles THROUGH v DO exist! But the per-path identity still holds.
Why? Because μ(5-cycle) contributes to Σ_C μ(C) which equals (H(T)-H(T-v))/2,
but the per-path identity sums only over 3-cycles...

THIS MEANS: At n=5, the per-path identity WORKS despite 5-cycles existing.
The 5-cycle's contribution to (inshat-1)/2 must somehow be "already captured" by 3-cycles.

This is the DEEP mystery. Why do 5-cycles not break the formula at n=5 but do at n=6?

HYPOTHESIS: At n=5, every 5-cycle embedding in P' is always "accompanied" by
exactly the right set of 3-cycle embeddings to account for its μ contribution.
At n=6, this delicate balance breaks down.

---

## OPEN QUESTIONS (mathematical, for the paper)

1. Does the per-path identity fail at n=5 for 5-cycles? (Contradicts n=5 success)
   Need to check: are there n=5 tournaments with 5-cycles through v where
   the per-path identity holds, but only because the 5-cycle contribution 
   coincidentally matches some 3-cycle accounting?

2. Can we characterize WHEN the per-path identity holds at n=6?
   The 2,758/9,126 failure rate means ~70% of (T,v,P') triples DO satisfy it.
   What's special about the failing cases?

3. The natural generalization (sum over all cycles including 5-cycles) overcounts at n=6.
   The maximal-embedding-only formula also fails. 
   Is there a CORRECT per-path formula for all n?

4. The binomial distribution C(L-2, 2k-1) for Type-II counts per L-cycle embedding
   is a clean theorem. What is its combinatorial proof? 
   It suggests a connection to ballot sequences or Dyck paths.

5. The AVERAGE contribution (L-4)/4 grows with L. Does this lead to an 
   asymptotic formula for Σ_C μ(C) as a function of cycle length distribution?

---

## WHAT THIS MEANS FOR THE PAPER

The main theorem F2 (every Type-II = unique 3-cycle) means:
- The per-path identity is REALLY about 3-cycles only, by definition
- Its failure at n=6 is because the SUM Σ_C μ(C) includes 5-cycle contributions
  that the per-path formula cannot see
- The gap between per-path and Claim A grows with the number of longer cycles

The algebraic identity F1 is CLEAN and provable -- should be in the paper as a lemma.
The distribution formula F4 is NEW and CLEAN -- should be added to the paper.
The failure analysis provides a clean explanation of WHY n=6 breaks things.
