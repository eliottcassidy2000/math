#!/usr/bin/env python3
"""Independent finite controls for prime-power Hall growth; no type census.

Universes: full permutations for K_{k,n-k}, 0<=n<=8, 0<=k<=n//2;
complete-bipartite arithmetic for 5<=n<=100; 1<=a,b<=4; and six
explicit balanced prime-power pairs. Exact arithmetic throughout.
"""
from itertools import permutations
from math import comb, factorial, gcd
from pathlib import Path
import json


def v_p(x, p):
    assert x > 0
    e = 0
    while x % p == 0:
        e += 1
        x //= p
    return e


def fact_v(n, p):
    ans = 0
    while n:
        n //= p
        ans += n
    return ans


def model(n, k):
    assert 0 <= k <= n // 2
    if k == 0:
        return factorial(n), True, True
    equal = 2 * k == n
    order = factorial(k) * factorial(n-k) * (2 if equal else 1)
    sign_even = not equal or k % 2 == 0
    eulerian = k % 2 == 0 and (n-k) % 2 == 0
    return order, sign_even, eulerian


def actual_graph(n, k):
    edges = {(u,v) for u in range(k) for v in range(k,n)}
    degree = [0] * n
    for u,v in edges:
        degree[u] += 1
        degree[v] += 1
    aut_order = 0
    odd_sign_automorphisms = 0
    for perm in permutations(range(n)):
        moved = {tuple(sorted((perm[u],perm[v]))) for u,v in edges}
        if moved != edges:
            continue
        aut_order += 1
        odd_sign_automorphisms += sum(perm[u] > perm[v] for u,v in edges) % 2
    return aut_order, odd_sign_automorphisms == 0, all(d % 2 == 0 for d in degree), odd_sign_automorphisms


def digit_leq(k,n,p):
    while k or n:
        if k % p > n % p:
            return False
        k //= p
        n //= p
    return True


def general_indices(n):
    tournament, switching = [], []
    for k in range(n//2 + 1):
        order, sign_even, eulerian = model(n,k)
        coprime = gcd(comb(n,k),15) == 1
        full_sylows = all(v_p(order,p) == fact_v(n,p) for p in (3,5))
        assert coprime == full_sylows
        if sign_even and coprime:
            tournament.append(k)
            if eulerian:
                switching.append(k)
    assert 0 in tournament and 0 in switching
    if n % 2:
        assert switching == [0]
    assert len({k*(n-k) for k in tournament}) == len(tournament)
    return tournament, switching


def sequence_case(a,b,dense=True):
    A,B = 3**a,5**b
    n,M = A*B-1,min(A,B)
    assert n >= 14 and n % 2 == 0 and M % 2 == 1
    assert 2*(M-1) < n
    js = range(1,M) if dense else sorted({1,2,M//2,M-1})
    ks = range(M) if dense else sorted({0,1,2,M//2,M-1})
    for j in js:
        assert v_p(j,3) < a and v_p(j,5) < b
        for p in (3,5):
            assert v_p(n+1-j,p) == v_p(j,p)
    for k in ks:
        if dense:
            assert gcd(comb(n,k),15) == 1
        for p in (3,5):
            assert fact_v(k,p)+fact_v(n-k,p) == fact_v(n,p)
    assert sum(k % 2 == 0 for k in range(M)) == (M+1)//2
    return {"a":a,"b":b,"n":n,"M":M,"tournament_lower_bound":M,"switching_lower_bound":(M+1)//2,"checked_index_mode":"all" if dense else "endpoint samples with exact factorial valuations"}


def main():
    graph_rows = []
    for n in range(9):
        for k in range(n//2+1):
            actual = actual_graph(n,k)
            assert actual[:3] == model(n,k), (n,k,actual,model(n,k))
            graph_rows.append({"n":n,"k":k,"automorphism_order":actual[0],"sign_even":actual[1],"eulerian":actual[2],"odd_sign_automorphisms":actual[3]})
    general = {str(n):dict(zip(("tournament_indices","switching_indices"),general_indices(n))) for n in range(5,101)}
    digit_cases = 0
    for n in range(101):
        for k in range(n+1):
            assert (gcd(comb(n,k),15) == 1) == all(digit_leq(k,n,p) for p in (3,5))
            digit_cases += 1
    extra_examples = {}
    for n in (224,674):
        ti,si = general_indices(n)
        extra_examples[str(n)] = {"tournament_bound":len(ti),"switching_bound":len(si)}
    assert extra_examples == {"224":{"tournament_bound":36,"switching_bound":18},"674":{"tournament_bound":27,"switching_bound":14}}
    sequence = [sequence_case(a,b) for a in range(1,5) for b in range(1,5)]
    balanced=[]
    for b in range(1,7):
        a=1
        while 3**a < 5**b:
            a+=1
        case=sequence_case(a,b,dense=False)
        assert 5**b <= 3**a < 3*5**b
        assert 3*case["M"]**2 > case["n"]+1
        balanced.append(case)
    # Hostile endpoint: the empty graph remains Eulerian for odd n.
    assert general["5"]["switching_indices"] == [0]
    # Hostile equal-part graphs: K_{3,3} is sign-odd; K_{4,4} is sign-even.
    assert model(6,3) == (72,False,False)
    assert model(8,4) == (1152,True,True)
    result={"status":"PASS","graph_permutation_universe":"0 <= n <= 8; 0 <= k <= floor(n/2)","graph_cases":len(graph_rows),"graph_rows":graph_rows,"arithmetic_universe":"5 <= n <= 100; 0 <= k <= floor(n/2)","general_arithmetic_cases":sum(n//2+1 for n in range(5,101)),"selected_general_indices":{str(n):general[str(n)] for n in (5,6,7,8,9,10,14,24,44,74,100)},"digit_condition_cases":digit_cases,"additional_arithmetic_examples":extra_examples,"prime_power_grid":sequence,"balanced_sequence":balanced,"note":"Exact finite controls only. Full-group permutation enumeration is restricted to individual complete bipartite graphs; there is no graph/tournament type census."}
    out=Path(__file__).with_name('independent_controls.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({"status":result["status"],"graph_cases":result["graph_cases"],"general_arithmetic_cases":result["general_arithmetic_cases"],"digit_condition_cases":digit_cases,"additional_arithmetic_examples":extra_examples,"prime_power_grid_cases":len(sequence),"balanced_sequence":balanced},indent=2))

if __name__ == '__main__':
    main()
