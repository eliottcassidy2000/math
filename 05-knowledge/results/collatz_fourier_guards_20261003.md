# Compressed route guards: complete phase sets and the Fourier information they need

2026-10-03 (America/Denver). Continuation of
[the atom and route-memory construction](collatz_atom_route_memory_20261003.md).

**Status:** PROVED elementary statements below; FINITE-EXACT for the explicit
script universes. This compiles legal edits inside an already parameterized
sibling family. Universal coverage of positive integers remains OPEN.

Inheritance: [THM-4501, recursive motif families](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md)
provides common-tail families and finite-prefix realization;
[the trunk chart, Theorem 1.1](collatz_procgen_20260923_trunk_rh.md)
provides the ternary isometry; and
[THM-2218, labelled guard-hole Fourier transform](../../01-canon/theorems/THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy.md)
identifies Fourier phase as information lost by scalar totals and magnitudes.
The latter is a mechanism transferred here, not an LRC consequence.
The hostile controls are the isolated home-cycle padding case and two guards
with equal cosine coefficients and Fourier magnitudes but different legality.
The necessary sidecars are phase origin, designated terminal, sheet, and the
distinction between a compressed head and its actual valuation word.

## 1. Compile every legal phase of a compressed head

Let U(n)=oddpart(3n+1) on positive odd integers. For odd m>0 with 3 not
dividing m, set k0(m)=2 if m=1 mod3, and 1 otherwise. Its odd predecessors are

    h_t = (2^(k0(m)+2t)*m-1)/3,   t>=0.

Fix a compressed head c=(j1,...,jp), ji>=0. Its inverse decoder reads the
indices from right to left. At current target h it rejects if 3 divides h;
otherwise its exponent is a=2ji+k0(h), and its next source is (2^a*h-1)/3.

**Guard theorem.** Exactly 2^p residue classes t modulo M=3^p pass all p
inverse steps. They can be compiled without orbit search as follows:

1. Choose each ai from {2ji+1,2ji+2}.
2. Form A=sum ai and the ordered carry
   S=sum_(i=0)^(p-1) 3^(p-1-i)*2^(a1+...+ai).
3. Find the unique t modulo M such that h_t=2^(-A)*S modulo M.

**Proof.** The identity v3(4^d-1)=1+v3(d) gives
v3(h_t-h_s)=v3(t-s). Thus t modulo M permutes h_t modulo M. For a fixed
exponent word, 3^p divides 2^A*h-S precisely on the displayed residue.
Reverse the final edge using this integrality condition modulo 3 and then
repeat; every reconstructed number is a positive odd integer. Its oddness
also proves that the indicated valuation is exact. Distinct chosen words
cannot share a residue: at each backward step the target modulo 3 fixes
the exponent parity, hence fixes that choice of ai. There are 2^p choices.

This theorem preserves the compressed indices ji. Actual exponents can
change by one between phases. Preserving one specified valuation word is
the narrower single-phase rule proved in the route-memory note.
The first compressed index j1 does not affect legality: changing it only
changes the final reconstructed source, which need not itself admit another
predecessor.

## 2. First-hit validity is not a periodic residue condition alone

Assume the terminal m carries a certified first-hit tail to 1. When h_t>1,
every legal reconstructed earlier source also exceeds 1: a forward route
starting at 1 stays at 1, so it cannot later reach h_t>1. Therefore the
entire decoded route is first-hit valid.

The unique exception is (m,t)=(1,0), where h_t=1. Whenever the arithmetic
head accepts it, the selected block contributes a forbidden extra 1->1.
For example c=(1), t=0 gives 5->1->1. Checking only the initial source>1
would miss this. Reject this isolated parameter or truncate at the first 1.
Do not remove the whole phase-zero class: t=M,2M,... have h_t>1.

For m=1 and head (0,0), the legal phases are {0,3,4,7} modulo 9:

| t | Source and route | Actual head |
|---|---|---|
| 3 | 75 -> 113 -> 85 -> 1 | (1,2) |
| 4 | 151 -> 227 -> 341 -> 1 | (1,1) |
| 7 | 19417 -> 14563 -> 21845 -> 1 | (2,1) |
| 9 | 621377 -> 466033 -> 349525 -> 1 | (2,2) |

These are the first positive parameter representatives of the four phases.
The residue-zero representative t=0 is the padding exception.

## 3. What the Fourier representation retains

For an increment-phase mask g on Z/M, use the unnormalized transform

    ghat(r)=sum_t g(t)*exp(-2*pi*i*r*t/M).

For every length-p head, g has 2^p ones. Parseval gives

    sum_(r!=0)|ghat(r)|^2 = M*2^p-4^p.

Thus all these heads have identical total nonconstant Fourier energy.
That statistic cannot distinguish their legal phase sets.

The smallest stronger hostile is exact. Keep compressed head (0) and root 1.
Starting internal index t=1 gives the route 3->5->1 and permits increment
phases {0,2} modulo 3. Starting t=3 gives 113->85->1 and permits {0,1}.
Their cosine coefficients are both (2,1/2,1/2); their squared Fourier
magnitudes are both (4,1,1). Their transforms are complex conjugates, and
the increment +1 is illegal for the first but legal for the second.
The sine/phase coordinate records precisely that directional distinction.

For signed edges 3N+epsilon=2^k*m, epsilon is s(m)*(-1)^k, where s(m) is
the representative +1 or -1 of m modulo 3. On the finite exponent coordinate
Z/(2M), sheet is therefore the Nyquist character (-1)^k. A same-sheet mask
on k=k0+2t has the exact transform

    Fhat(r)=exp(-pi*i*r*k0/M)*ghat(r mod M),
    Fhat(r+M)=(-1)^k0*Fhat(r).

This is a finite-coordinate identity retaining sheet and ternary phase;
it does not replace the positive-exponent or first-hit guards.

For a library of such masks, cyclic correlation computes overlap and hence
the additional phases covered by a proposed rule. The relative Fourier
phases are needed for that correlation. Implement the masks or integral
group algebra exactly; numerical Fourier magnitudes alone are insufficient.
Even complete phase coverage describes one sibling-index family, not every
positive integer.

## Reproduction and audit boundary

Run `python3 04-computation/experiments/collatz_fourier_guards_20261003.py`.
The saved [output](collatz_fourier_guards_20261003.out) records:

- All 120 compressed heads with lengths 1..4 and indices 0..2; 1554 distinct
  valuation phases; independent backward and forward replay of 14760 index
  cases covering two periods for each head.
- All 45 arithmetic-legal index-zero cases in that universe rejected as
  full first-hit codes. The additional head (1) control catches an early 1
  despite initial source 5.
- Exact Parseval checks in Z[zeta_(3^p)], reduced by
  Phi_(3^p)(X)=X^(2*3^(p-1))+X^(3^(p-1))+1. No floating DFT is used.
- Additional targets 5,7,13 and four stated heads, every phase; these test
  arithmetic legality without assuming any new terminal convergence.
- The four explicit routes, the phase-loss hostile in Q(zeta_3), and 480
  formal-monomial checks of the even/odd Fourier lift.

The all-parameter proofs supply the theorem; finite tests do not supply a
coverage claim. An independent reviewer rederived the phase count and the
isolated first-hit exception without using this implementation.
