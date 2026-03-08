# Hypothesis Log — Index

Every hypothesis tested in this project, whether confirmed, refuted, or open.
Organized by topic. Each hypothesis has a detail file.

**Status codes:** CONFIRMED | REFUTED | OPEN | PARTIALLY-TRUE

---

## By Status

### CONFIRMED
| ID | Statement | Why it works | Source |
|----|-----------|-------------|--------|
| HYP-001 | M[a,b] = M[b,a] (transfer matrix symmetric) | Unknown mechanism; verified n<=8 | INV-001 |
| HYP-002 | H(T) = H(T^op) | Even-r structure of W; THM-L | opus-S35c11 |
| HYP-003 | M(r) symmetric for ALL r, not just r=1/2 | Unknown; verified n<=5 | opus-S35c11 |
| HYP-004 | W(T,r) = (r-1/2)^{n-1} F(T, (2r+1)/(2r-1)) | Mobius transform; THM-K | opus-S35c11 |
| HYP-005 | Deletion-contraction: H(T) = H(T\e) + H(T/e) | Standard DC on Redei-Berge function | opus-S35c11 |
| HYP-006 | Var(c_0)/Var(H) -> 0 exponentially | Perpendicularity; Fourier hierarchy | opus-S35c11 |
| HYP-007 | Singleton cancellation in M | Algebraic swap proof | THM-055 |
| HYP-008 | GS flip orthogonality Cov(H_sym, H_anti)=0 | Trivial for ANY involution on uniform measure | opus-S35c11 |
| HYP-009 | Blueself impossible at odd n | Grid symmetry + endpoint constraint | INV-039 |
| HYP-010 | Every SC tournament has involution anti-aut | Moon's theorem + Cauchy | THM-024 |
| HYP-011 | F_k(T) = F_{n-1-k}(T^op) complement duality | Path reversal maps ascents to descents | worpitzky_F_at_2.py |
| HYP-012 | Sum_T F(T,x) = A_n(x) * 2^{C(n,2)-(n-1)} | Each perm is HP of 2^{extra} tournaments | worpitzky_restricted_eulerian.py |

### REFUTED
| ID | Statement | Why it fails | First failure | Source |
|----|-----------|-------------|---------------|--------|
| HYP-101 | Per-path identity holds for all n | 3-cycle-only formula misses longer cycles | n=6 | MISTAKE |
| HYP-102 | M(T) = M(T^op) | Fails for generic tournaments | n=4 | reversal_proof_attempt.py |
| HYP-103 | M(T) = M(T^op)^T | Also fails | n=4 | reversal_proof_attempt.py |
| HYP-104 | M is linear in A, A^T, A^2, etc. | No such combination works | n=4 | simple_pattern_search.py |
| HYP-105 | M commutes with A | Fails | n=4 | M_functional_equation.py |
| HYP-106 | M(T\e) symmetric (for deletion-contraction proof) | Deletion breaks tournament structure | n=4 | dc_symmetry_path.py |
| HYP-107 | M(T/e) symmetric (for contraction) | Also fails | n=4 | dc_symmetry_path.py |
| HYP-108 | Even-odd split equivalent to OCF | Split is CONSEQUENCE, not equivalent | - | MISTAKE-008 |
| HYP-109 | Omega(T) always claw-free | Fails at n=9 (90%) | n=9 | INV-032 |
| HYP-110 | I(Omega(T),x) always real-rooted | Fails n=9 (maximally rare) | n=9 | THM-025 |
| HYP-111 | Omega(T) always perfect | Fails n=8 (53.8%) | n=8 | INV-032 |
| HYP-112 | Blueself = SC maximizer | Fails n=6 | n=6 | INV-040 |
| HYP-113 | Fourier degree-4 "two type" decomposition generalizes to n=9 | Dimension >> 200, non-proportional coefficients | n=9 | INV-050 |
| HYP-119 | W(T,r) = (H/2^{n-1})(2r+1)^{n-1} (type B universal) | ARTIFACT of wrong W definition | worpitzky_typeB_verify.py |
| HYP-120 | F(T,x) is palindromic | FALSE: F_k(T) = F_{n-1-k}(T^op), not F_{n-1-k}(T) | worpitzky_ehrhart.py |
| HYP-114 | GS flip orthogonality explains perpendicularity | Trivially true for any involution; no content | - | opus-S35c11 |
| HYP-115 | Contraction approach proves symmetry | Dead end | - | T017 |
| HYP-116 | Contiguous block decomposition | Dead end | - | T035 |
| HYP-117 | Per-vertex decomposition of unmatched counts | Dead end | - | T045 |
| HYP-118 | Cycle bijection under arc reversal | Dead end | - | MISTAKE-005 |

### OPEN
| ID | Statement | Current evidence | Source |
|----|-----------|-----------------|--------|
| HYP-201 | Char poly determines H exactly | 0 ambiguous at n=5; 3 at n=6; 36 at n=7 | char_poly_H.py |
| HYP-202 | Spectral det(I-uA) correlates with H | Corr > 0.94 for optimal u | ihara_deep.py |
| HYP-203 | M has algebraic formula in terms of cycle invariants | Partial: diagonal formula known | THM-053 |
| HYP-204 | F(T,x) ascent poly is ALWAYS unimodal | 100% at n=3-5 (exhaustive), n=6-8 (sampling) | worpitzky_roots.py |
| HYP-205 | F(T,x) is almost always log-concave | 100% except 4/1024 at n=5; 100% at n=6-8 sampling | worpitzky_roots.py |
| HYP-206 | Roots of F(T,x) are real ~89% of the time (stable rate) | n=5: 95%, n=6: 90%, n=7: 89%, n=8: 89% | worpitzky_n6_test.py |

---

## By Topic

### Transfer matrix symmetry
HYP-001, HYP-003, HYP-102, HYP-103, HYP-104, HYP-105, HYP-106, HYP-107, HYP-203

### Complement invariance
HYP-002, HYP-102, HYP-103

### Fourier / W-polynomial
HYP-004, HYP-006, HYP-008, HYP-113, HYP-114

### OCF structure
HYP-005, HYP-101, HYP-108, HYP-115, HYP-116, HYP-117, HYP-118

### Omega(T) properties
HYP-109, HYP-110, HYP-111

### Self-complementary / blueself
HYP-009, HYP-010, HYP-112

### Worpitzky / Forward-edge polynomial
HYP-011, HYP-012, HYP-119, HYP-120, HYP-204, HYP-205, HYP-206

### Spectral
HYP-201, HYP-202

---

## Template for new hypothesis files

```markdown
# HYP-NNN: [Short title]

**Status:** OPEN | CONFIRMED | REFUTED | PARTIALLY-TRUE
**Filed:** [session-id]
**Resolved:** [session-id] (if applicable)

## Statement
[Precise mathematical statement]

## Motivation
[Why we thought this might be true]

## Test method
[How it was tested — script name, parameters, sample size]

## Outcome
[What happened]

## WHY it works/fails
[The crucial insight — what structural/algebraic reason makes this true or false]

## Related
- Variables: [links]
- Hypotheses: [links]
- Theorems: [links]
- Scripts: [links]

## Tags
#tag1 #tag2
```
