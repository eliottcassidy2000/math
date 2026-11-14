# Exact short Collatz certificates for generic tournament windows

**PROVED:** an explicit one-pass feasibility test for a proposed generic *ordinal permutation*, a linear total-valuation bound, and a linear-bit physical positive seed on either Syracuse sheet. This does not say every actual window is generic, or that all odd path counts occur in the generic window family.

## Inheritance and the distinction to keep

The earlier [general Syracuse alphabet](../syracuse_window_alphabet_general_m_20261005.md) reads a window as chronological path edges plus larger-to-smaller chords. For a finite valuation word it uses the coefficient comparison 3^k versus 2^A; exact actual comparisons can differ at a finite set of small starts for a fixed word. It solves generic chord feasibility by difference constraints. Here we retain the entire comparison permutation, including the adjacent ones that the tournament forgets, and solve that lifted problem explicitly. We also bound its physical realization cost uniformly in the order.

## The least level staircase

Let m be the number of values, i=0,...,m−1, and write alpha=log_2(3), b_i=floor(i alpha), phi_i=i alpha−b_i. If a_i≥1 are the prescribed divisions, use q_i=sum_(j=1)^i a_j. The log slopes of the iterates are i alpha−q_i=z_i+phi_i with z_i=b_i−q_i an integer. A sequence of heights represents a valuation word if, up to global translation to make z_0=0,

    z_i−z_(i−1) ≤ eps_i := b_i−b_(i−1)−1,  eps_i in {0,1}.

The phases are all distinct. Given an ascending ordinal permutation pi of the time indices, start at height 0 for pi(0), retaining the current level as you move along pi except for an increase by one at each strict drop in phase. This is the least staircase in sorted order. Check the displayed chronological inequalities on this staircase. The permutation is generic feasible iff this succeeds.

Indeed any two consecutive elements in sorted slope order need a nonnegative integer height separation, at least one when the phase wraps down. Their differences in any feasible assignment are the least staircase gaps plus nonnegative extras. An upward chronological move has minimal height difference in the least staircase. A downward move is automatically legal: its integer height difference is nonpositive. So a violation in the least staircase cannot be fixed, and success directly supplies the divisions. This also gives the exact generic ordinal test without computing with real logarithms: b_i is the bit length of 3^i minus 1, and phase order is the rational order of 3^i/2^b_i. For a target tournament *with a fixed chronological path*, chord comparisons specify only a partial order; generic feasibility means some linear extension passes. This is not a polynomial algorithm for an arbitrary unmarked tournament supplied as input.

## Certificate length and a stronger bound for strong tournaments

Put k=m−1. The height range is at most k by construction, whence total A=q_k=b_k+z_0−z_k≤b_k+k. In a strong window, any integer height boundary strictly between occupied extremes must be crossed upward at least once. Otherwise all visits to the upper side are followed by all visits to the lower side, since both sides are visited and there is no return upward. Every cross arc is forward in time: the chronological seam arc is forced, and all others follow decreasing score. This contradicts strongness. Each upward chronological step crosses at most one integer boundary, and there are at most sum eps_i=b_k−k upward steps. Hence for every generic strong ordinal type our compiler achieves A≤2b_k−k. The singleton is included with A=0.

## Turn the word into actual integers

For sigma=+1 or −1 use S_sigma(x)=(3x+sigma)/2^v2(3x+sigma) on odd positive integers. A word a of total A gives one odd residue class modulo 2^(A+1). Compute it from odd terminal residue modulo 2 by unwinding (2^a y−sigma)/3 in successively larger powers of 2. Induction gives exactly the specified valuations, not just inequalities, because every terminal quotient has been kept odd.

For i≥0 the affine iterate on that cylinder is x_i=3^i x/2^q_i + sigma sum_(r<i) 3^(i−1−r)2^(q_r−q_i). After multiplication by 2^A, every pair has nonzero integral leading-coefficient difference, and the absolute total intercept difference is <3^m 2^A. Thus taking the first seed x>3^m 2^A in the cylinder guarantees that all ranks are the generic ranks (and all values are distinct). This gives a fully explicit O(m)-bit bound, uniformly for both sheets; the cylinder backward pass is the familiar exact-itinerary idea of [Bernstein–Lagarias](https://websites.umich.edu/~lagarias/doc/bernstein.pdf).

A sharper asymptotic existence statement is x=O(2^A) uniformly for all valuation words. For a pair i<j, let l=j−i and s=q_j−q_i. Write x_j=(3^l x_i+sigma C)/2^s. Since at least one division occurs per step, C/2^s ≤ (3/2)^l. Also x_i=(3^i/2^q_i)(x+sigma R_i), with 0≤R_i≤2^A/2. So its potential comparison switch is at absolute x at most
```
2^A/2 + 2^A (3/2)^l / |2^s−3^l|.
```
The last quotient is bounded absolutely. Outside 1/2≤2^s/3^l≤2 it is ≤2^(1−l); inside that range s=O(l), and [Matveev's lower bound for linear forms in logarithms](https://m.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=314&wshow=paper) yields |2^s−3^l|≥c3^l l^(−C_0) with absolute c,C_0>0 after adjusting small cases. A polynomial times 2^(−l) is bounded. Raising x beyond the resulting uniform threshold and choosing its residue class costs an absolute constant times 2^A. Therefore an m-vertex generic strong pattern has positive witnesses on both sheets of binary length at most (2 log_2(3)−1)(m−1)+O(1). The logarithm theorem is not necessary for the effective elementary compiler above.

## What this gets and does not get

The bridge is a deterministic compiler once the ordinal target is supplied, not an instruction for choosing the ordinal target from H. All-order generic realizability transfers to both signs, but a fixed initial integer selects only one itinerary and may not generically shadow a proposed local window after an arbitrary prefix. The projection to H loses adjacencies used by the compiler. A next positive test is to find a graph move that preserves the chronological feasibility inequalities while changing H by a prescribed integer with a bounded carry; the present theorem cannot supply that move.

[Replay](replay_windows.py) independently generates all capped words through order 6 and compares their ordinal support with the test. It also tests prepared positive seeds through order 250 and the strong-bound cases. No claim of literature priority.
