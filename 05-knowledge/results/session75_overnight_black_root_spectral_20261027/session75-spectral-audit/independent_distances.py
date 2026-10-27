#!/usr/bin/env python3
"""Recompute maximum pair correlation to join/C3 grammars from scratch.

No producer import or supplied lower bound is used in the optimization.
The optional certificate is compared only after all optima are recomputed.
"""
from functools import lru_cache
import argparse
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def submasks(s):
    t = s
    while t:
        yield t
        t = (t-1) & s


def solve(n):
    require(n in (3,7,11), "fixed audited prime universe")
    squares = {x*x % n for x in range(1,n)}
    k = [[0 if i == j else 1 if (j-i) % n in squares else -1 for j in range(n)] for i in range(n)]
    vertices = [tuple(i for i in range(n) if s >> i & 1) for s in range(1 << n)]
    @lru_cache(None)
    def cross(a,b):
        return sum(k[i][j] for i in vertices[a] for j in vertices[b])
    # reward[s] = maximum sum_{i<j in s} target_ij * baseline_ij.
    reward = [0]*(1 << n)
    counts = [0,0]
    for cardinality in range(2,n+1):
        for s in range(1 << n):
            if len(vertices[s]) != cardinality:
                continue
            first = s & -s
            best = None
            # Each unordered binary partition once: A contains least vertex.
            for b in submasks(s ^ first):
                a = s ^ b
                value = reward[a]+reward[b]+abs(cross(a,b))
                best = value if best is None else max(best,value)
                counts[0] += 1
            # Each unordered ternary partition once: A has least vertex,
            # B has least vertex of its complement. Cycle/reverse exhaust
            # the two quotient orientations; all other 3-quotients join.
            for bc in submasks(s ^ first):
                a = s ^ bc
                second = bc & -bc
                for c in submasks(bc ^ second):
                    b = bc ^ c
                    value = reward[a]+reward[b]+reward[c]+abs(cross(a,b)+cross(b,c)+cross(c,a))
                    best = value if best is None else max(best,value)
                    counts[1] += 1
            require(best is not None, "missing nontrivial partition")
            reward[s] = best
    distance = []
    for s,r in enumerate(reward):
        pairs = len(vertices[s])*(len(vertices[s])-1)//2
        require((pairs-r) % 2 == 0, "correlation-distance parity")
        distance.append((pairs-r)//2)
    require(distance[-1] == {3:0,7:3,11:11}[n], "independent minimum")
    return distance, counts


def witness_distance(node,n):
    matrix = [[0]*n for _ in range(n)]
    seen = set()
    def decode(t):
        if "vertex" in t:
            v=t["vertex"]
            require(type(v) is int and 0<=v<n and v not in seen,"witness vertex labels")
            seen.add(v)
            return [v]
        kind,sgn=t["kind"],t["orientation"]
        require(type(sgn) is int and sgn in (-1,1),"witness quotient orientation")
        children=t["children"]
        q=len(children)
        require((kind=="join" and q==2) or (kind=="cycle" and q==3),"witness quotient type")
        blocks=[decode(c) for c in children]
        for i in range(q):
            for j in range(i+1,q):
                arc=sgn if q==2 or (j-i)%3==1 else -sgn
                for x in blocks[i]:
                    for y in blocks[j]:
                        matrix[x][y]=arc
                        matrix[y][x]=-arc
        return sum(blocks,[])
    decode(node)
    require(seen==set(range(n)),"witness omitted vertices")
    squares={x*x % n for x in range(1,n)}
    exceptions=[]
    for i in range(n):
        for j in range(i+1,n):
            target=1 if (j-i)%n in squares else -1
            require(matrix[i][j] in (-1,1),"missing witness pair")
            if matrix[i][j]!=target:
                exceptions.append([i,j])
    return exceptions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate")
    args = parser.parse_args()
    results = []
    computed = {}
    for n in (3,7,11):
        values, counts = solve(n)
        computed[n] = values
        results.append({"n":n,"minimum_arc_edits":values[-1],"subset_optima_recomputed":len(values),"root_partition_counts":counts})
    if args.certificate:
        with open(args.certificate) as f:
            cases = json.load(f)["cases"]
        require([c["n"] for c in cases] == [3,7,11],"certificate universe")
        for case,result in zip(cases,results):
            n = case["n"]
            require(computed[n] == case["subset_lower_bounds"],"entire optimum vector differs from certificate")
            require(computed[n][-1] == case["minimum"],"certificate minimum")
            exceptions=witness_distance(case["witness"],n)
            require(len(exceptions)==computed[n][-1],"witness does not attain independent optimum")
            result["witness_exceptions"]=exceptions
    print(json.dumps({"status":"PASS","method":"independent maximum-correlation subset optimization, no producer imports","certificate_compared":bool(args.certificate),"cases":results},indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
