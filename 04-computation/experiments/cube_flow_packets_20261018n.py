#!/usr/bin/env python3
"""Exact signed cube flows from the requested sporadic-sequence papers.

Finite universe U(n), n=0..8. Independent tests use a recurrence and the
single-sum formula, and ordered matching decompositions. No real-energy
bound is inferred from a p-adic congruence.
"""
import argparse
from collections import Counter
from itertools import permutations
import json
from math import comb, factorial, prod
from pathlib import Path

EDGES = ((0,1),(0,2),(0,3),(1,0),(1,1),(1,2),
         (2,0),(2,2),(2,3),(3,0),(3,1),(3,3))
CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def boundary(flow):
    return tuple(sum(flow[k] for k,e in enumerate(EDGES) if e[side] == i)
                 for side in (0,1) for i in range(4))


def flows(n):
    triples = [(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)]
    for a in triples:
        for b in triples:
            if a[0]+b[1] > n or a[1]+b[2] > n:
                continue
            for c in triples:
                if a[1]+b[2]+c[1] != n:
                    continue
                d = (n-b[0]-c[0], n-a[0]-b[1], n-a[2]-c[2])
                if min(d) >= 0 and sum(d) == n:
                    yield a+b+c+d


def weight(flow, n):
    return (-1)**(flow[1]+flow[3]+flow[11]) * factorial(n)**4 // prod(factorial(x) for x in flow)


def az(n):
    return sum((-1)**k * 3**(n-3*k) * comb(n,3*k) * comb(n+k,n)
               * comb(3*k,2*k) * comb(2*k,k) for k in range(n//3+1))


def polynomial_product(a, b):
    out = Counter()
    for x, u in a.items():
        for y, v in b.items():
            out[tuple(i+j for i,j in zip(x,y))] += u*v
    return {x:v for x,v in out.items() if v}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    matchings = []
    for p in permutations(range(4)):
        if all((i,p[i]) in EDGES for i in range(4)):
            matchings.append(tuple(int(p[i] == j) for i,j in EDGES))
    check(len(matchings) == 9, "nine matchings of the cube")
    check(Counter(weight(m,1) for m in matchings) == {1:6,-1:3}, "signed initial permanent")
    matrix = [[0]*4 for _ in range(4)]
    for k,(i,j) in enumerate(EDGES):
        matrix[i][j] = -1 if k in (1,3,11) else 1
    check(all(sum(matrix[i][k]*matrix[j][k] for k in range(4)) == 3*(i==j)
              for i in range(4) for j in range(4)), "orthogonal signed cube matrix")
    c = [[(1,1,-1,1)[i]*matrix[i][j] for j in (0,3,1,2)] for i in range(4)]
    check(all(c[i][j] == -c[j][i] for i in range(4) for j in range(4)), "skew conference gauge")
    h4 = [[c[i][j]+(i==j) for j in range(4)] for i in range(4)]
    raw = [row[:] for row in h4]
    h4 = [[raw[i][j]*raw[0][j]*raw[i][0]*raw[0][0] for j in range(4)] for i in range(4)]
    h8 = [[sgn*h4[i][j] for sgn in (1, s) for j in range(4)] for s in (1,-1) for i in range(4)]
    code = {sum(1<<i for i,x in enumerate(row) if x*s == -1) for row in h8 for s in (1,-1)}
    check(len(code) == 16 and all((a^b) in code for a in code for b in code)
          and Counter(v.bit_count() for v in code) == {0:1,4:14,8:1},
          "conference completion and doubling give Hamming code")
    quartic = {(0,)*4:1}
    for row in matrix:
        linear = {tuple(int(i==j) for i in range(4)):row[j] for j in range(4) if row[j]}
        quartic = polynomial_product(quartic,linear)
    power = {(0,)*4:1}
    fock = []
    for n in range(9):
        norm = sum(v*v*prod(factorial(k) for k in monomial) for monomial,v in power.items())
        check(norm == 3**(4*n)*factorial(n)**4, "full Fischer norm identity")
        check(power.get((n,)*4,0) == az(n) and abs(az(n)) <= 9**n,
              "balanced source-target amplitude bound")
        fock.append({"n":n,"monomials":len(power),"weighted_norm":norm,"bound":9**n})
        power = polynomial_product(power,quartic)
    # The adjacent three-channel construction is the Fourier matrix over Z[omega].
    cubic = {(3,0,0):1,(0,3,0):1,(0,0,3):1,(1,1,1):-3}
    power3, fourier3 = {(0,0,0):1}, []
    previous_b, current_b = 0, 1
    for n in range(9):
        value = (-1)**n*power3.get((n,)*3,0)
        check(value == current_b, "Zagier B recurrence against Fourier coefficient")
        norm3 = sum(v*v*prod(factorial(k) for k in a) for a,v in power3.items())
        check(norm3 == 3**(3*n)*factorial(n)**3 and value*value <= 27**n,
              "three-channel Fourier Fischer norm")
        fourier3.append({"n":n,"B_n":value,"weighted_norm":norm3,"squared_bound":27**n})
        numerator = (9*n*(n+1)+3)*current_b-27*n*n*previous_b
        check(numerator % (n+1)**2 == 0, "integral recurrence")
        previous_b, current_b = current_b, numerator//(n+1)**2
        power3 = polynomial_product(power3,cubic)
    # An even clock perturbation is invisible at z=-1, but not as a Laurent identity.
    clocks = [0]*12
    clocks[0] = 2
    gram01 = Counter()
    for j in range(4):
        if matrix[0][j] and matrix[1][j]:
            a, b = EDGES.index((0,j)), EDGES.index((1,j))
            gram01[clocks[a]-clocks[b]] += matrix[0][j]*matrix[1][j]
    check(dict(gram01) == {2:1,0:-1} and sum(v*(-1)**k for k,v in gram01.items()) == 0,
          "integer-clock hostile survives the parity test")
    curl = (1,-1,0,0,-1,1,0,0,0,0,0,0)
    check(boundary(curl) == (0,)*8, "four-cycle keeps all eight boundary totals")
    ordered = Counter({(0,)*12: 1})
    receipts = []
    minimal_weight_loss = None
    for n in range(9):
        universe = list(flows(n))
        check(set(universe) == set(ordered), "integer matching decomposition covers whole fibre")
        signed, unsigned, curl_moves = 0, 0, 0
        primitive = {p: Counter() for p in (5,7) if n and n % p == 0}
        for f in universe:
            check(boundary(f) == (n,)*8, "all source and target margins")
            w = weight(f,n)
            signed += w
            unsigned += abs(w)
            if abs(w) != ordered[f] and minimal_weight_loss is None:
                minimal_weight_loss = {"n":n, "flow":f, "microstate_weight":abs(w),
                                       "ordered_matching_lifts":ordered[f]}
            remaining, chosen = list(f), []
            for _ in range(n):
                m = next((m for m in matchings if all(a>=b for a,b in zip(remaining,m))), None)
                check(m is not None, "constructive Hall peeling")
                chosen.append(m)
                remaining = [a-b for a,b in zip(remaining,m)]
            check(remaining == [0]*12 and len(chosen) == n, "marked flow reconstruction")
            g = tuple(a+b for a,b in zip(f,curl))
            if min(g) >= 0:
                check(boundary(g) == boundary(f) and weight(g,n)*w < 0,
                      "curl reverses sign but preserves margins")
                curl_moves += 1
            for p, ledger in primitive.items():
                number = sum(any(x % p for x in f[3*i:3*i+3]) for i in range(4))
                ledger[number] += w
        check(signed == az(n), "independent single-sum formula")
        if n >= 2:
            check(n**3*signed == (2*n-1)*(7*(n-1)*n+3)*az(n-1)-81*(n-1)**3*az(n-2),
                  "independent three-term recurrence")
        for p, ledger in primitive.items():
            check(ledger[1] == 0 and sum(v for k,v in ledger.items() if k) % p**3 == 0,
                  "primitive contribution divisible by p cubed, finite r=1")
        receipts.append({"n":n, "flows":len(universe), "signed_weight":signed,
                         "unsigned_weight":unsigned, "ordered_matching_lifts":sum(ordered.values()),
                         "signed_ordered_matching_sum":3**n, "feasible_positive_curls":curl_moves,
                         "primitive_row_ledger":{str(p):dict(sorted(v.items())) for p,v in primitive.items()}})
        if n < 8:
            nxt = Counter()
            for f,v in ordered.items():
                for m in matchings:
                    nxt[tuple(a+b for a,b in zip(f,m))] += v
            ordered = nxt
    # Higher exponent checks use the independent exact sum, not sampled flow sets.
    congruences = []
    for p in (5,7,11):
        for r in (1,2):
            n = p**r
            difference = az(n)-az(n//p)
            check(difference % p**(3*r) == 0, "declared finite supercongruence")
            congruences.append({"p":p,"r":r,"modulus":p**(3*r),"remainder":0})
    out = {"status":"FINITE-EXACT; Hall decomposition and curl identities are proved in companion note",
           "edges_in_paper_order":EDGES,"matching_count":9,"cube_cycle_rank":5,
           "curl":curl,"flow_receipts":receipts,"minimal_weight_loss":minimal_weight_loss,
           "signed_matrix":matrix,"skew_conference_gauge":c,"Hamming_code_from_cube":sorted(code),
           "full_Fischer_norm_checks":fock,"clock_hostile_gram_01":dict(sorted(gram01.items())),
           "Fourier3_Zagier_B_norm_checks":fourier3,
           "supercongruence_checks":congruences,"checks":CHECKS}
    encoded = json.dumps(out,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded,end="")


if __name__ == "__main__":
    main()
