#!/usr/bin/env python3
"""
Bridge the gap: from even-odd split to OCF.

We know:
  sum_S (-1)^|S| Delta(S,R) = 0  [even-odd split, THM-017, proved for all n]
  sum_S Delta(S,R) = delta_H      [definition]
  We want: delta_H = delta_I      [OCF]

Question: can we decompose Delta(S,R) further to get delta_H = delta_I?

Key idea: Delta(S,R) = Li(S)*Rj(R) - Lj(S)*Ri(R).
The even-odd split says the alternating sum vanishes.
OCF says the unsigned sum equals the cycle formula.

Approach 1: Check if Delta(S,R) can be expressed in terms of cycle contributions.
Approach 2: Look at the "size generating function" F(x) = sum_S x^|S| Delta(S,R).
  We know F(-1) = 0 and F(1) = delta_H. What is the structure of F(x)?

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
import random


def h_end(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def h_start(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[0] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def compute_Delta_by_size(n, arc_values):
    """Compute Delta(S,R) for each subset S, grouped by |S|."""
    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)

    def T(a, b):
        if a == b:
            return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    by_size = {}
    for smask in range(1 << m):
        S = [W[bit] for bit in range(m) if smask & (1 << bit)]
        R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
        k = len(S)

        Li_S = h_end(T, S + [I], I)
        Rj_R = h_start(T, [J] + R, J)
        Lj_S = h_end(T, S + [J], J)
        Ri_R = h_start(T, [I] + R, I)

        Delta = Li_S * Rj_R - Lj_S * Ri_R
        by_size[k] = by_size.get(k, 0) + Delta

    return by_size


def compute_delta_H(n, arc_values):
    """Compute delta_H = H(T) - H(T') for arc flip i->j."""
    # Count Ham paths in T using arc 0->1, and in T' using arc 1->0
    V = list(range(n))
    I, J = 0, 1

    def T(a, b):
        if a == b:
            return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

    # adj(i,j) in T: paths with i immediately before j
    adj_T = 0.0
    for p in permutations(V):
        w = 1.0
        for k in range(n - 1):
            w *= T(p[k], p[k + 1])
        # Check if i->j appears consecutively
        for k in range(n - 1):
            if p[k] == I and p[k + 1] == J:
                adj_T += w
                break

    # adj(j,i) in T': paths with j immediately before i
    # In T', the only change is arc i->j becomes j->i
    # But paths not using arc i->j are unchanged
    # Actually, T' has T[j][i]=1 instead of T[i][j]=1. Everything else same.
    arc_values_prime = dict(arc_values)
    arc_values_prime[(I, J)] = 0.0
    arc_values_prime[(J, I)] = 1.0

    def Tp(a, b):
        if a == b:
            return 0
        return arc_values_prime.get((a, b), 1 - arc_values_prime.get((b, a), 0))

    adj_Tp = 0.0
    for p in permutations(V):
        w = 1.0
        for k in range(n - 1):
            w *= Tp(p[k], p[k + 1])
        for k in range(n - 1):
            if p[k] == J and p[k + 1] == I:
                adj_Tp += w
                break

    return adj_T - adj_Tp


def analyze_size_generating_function(n, num_trials=30):
    """Analyze F(x) = sum_{k=0}^m C_k x^k where C_k = sum_{|S|=k} Delta(S,R)."""
    print(f"=== Size generating function at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    for trial in range(min(num_trials, 5)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.randint(0, 1) * 1.0
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        by_size = compute_Delta_by_size(n, arc_values)

        # F(-1) should be 0 (even-odd split)
        F_neg1 = sum((-1) ** k * by_size.get(k, 0) for k in range(m + 1))
        # F(1) should be delta_H
        F_1 = sum(by_size.get(k, 0) for k in range(m + 1))
        delta_H = compute_delta_H(n, arc_values)

        print(f"  trial {trial}: F(-1)={F_neg1:.6f}, F(1)={F_1:.6f}, delta_H={delta_H:.6f}")
        print(f"    Coefficients by |S|: {[f'{by_size.get(k,0):.4f}' for k in range(m+1)]}")

        # Check: is F(x) = (1+x) * Q(x) for some Q?
        # If F(-1)=0, then (1+x) divides F(x).
        # Q(x) = F(x)/(1+x). Compute Q at x=1: Q(1) = F(1)/2 = delta_H/2
        if abs(F_neg1) < 1e-8:
            Q_1 = F_1 / 2
            print(f"    F(1)/2 = {Q_1:.4f} = delta_H/2")

            # Compute Q(x) coefficients by polynomial division
            C = [by_size.get(k, 0) for k in range(m + 1)]
            Q = []
            remainder = C[:]
            for k in range(m):
                q_k = remainder[k]
                Q.append(q_k)
                remainder[k] -= q_k
                remainder[k + 1] -= q_k  # subtract q_k * (1+x) contribution
            print(f"    Q coefficients: {[f'{q:.4f}' for q in Q]}")
            # Check: Q(1) = sum(Q) should equal delta_H/2
            print(f"    Q(1) = {sum(Q):.4f}, delta_H/2 = {delta_H/2:.4f}")

        print()


def analyze_factored_form(n, num_trials=50):
    """
    Since F(-1) = 0, we have F(x) = (1+x)*Q(x).
    Study Q(x): what is Q(1) = delta_H/2? Does Q(1) = delta_I/2?

    At n=4: m=2, F(x) = C_0 + C_1*x + C_2*x^2, Q(x) linear.
    """
    print(f"=== Factored form analysis at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)
    max_err = 0.0

    for trial in range(num_trials):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.randint(0, 1) * 1.0
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        by_size = compute_Delta_by_size(n, arc_values)
        delta_H = compute_delta_H(n, arc_values)

        F_1 = sum(by_size.get(k, 0) for k in range(m + 1))
        err = abs(F_1 - delta_H)
        max_err = max(max_err, err)

    print(f"  Max |F(1) - delta_H|: {max_err:.2e}")
    print(f"  (Confirms F(1) = delta_H = sum_S Delta(S,R))\n")


if __name__ == "__main__":
    for n in range(4, 7):
        analyze_size_generating_function(n, num_trials=5)

    for n in range(4, 7):
        analyze_factored_form(n, num_trials=30)
