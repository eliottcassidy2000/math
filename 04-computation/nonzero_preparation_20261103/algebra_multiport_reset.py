#!/usr/bin/env python3
"""Multi-exit two-birth reset: exact operator and literal lift injections.

Standalone permutation enumeration through order4 for every exit subset/root;
strong-input scalar probe through5. No imported producer or input receipts.
"""
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json

GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def graph(n, code):
    rows = [0] * n
    for i, (a, b) in enumerate(combinations(range(n), 2)):
        if code >> i & 1:
            rows[a] |= 1 << b
        else:
            rows[b] |= 1 << a
    return tuple(rows)


def reset(T, B):
    n, full = len(T), (1 << len(T)) - 1
    first = tuple(row | ((not (B >> i & 1)) << n) for i, row in enumerate(T)) + (B,)
    final = tuple(row | (1 << (n + 1)) for row in first[:n]) + (first[n], 1 << n)
    return first, final


def reach(T):
    R = [r | (1 << i) for i, r in enumerate(T)]
    for k in range(len(T)):
        for i in range(len(T)):
            if R[i] >> k & 1:
                R[i] |= R[k]
    return R


def strong(T):
    return len(T) >= 3 and all(r == (1 << len(T)) - 1 for r in reach(T))


def tables(T, joints=True):
    n, size = len(T), 1 << len(T)
    P = [0] * size
    P[0] = 1
    A = [[0] * size for _ in T]
    E = [[0] * size for _ in T]
    Q = [[[0] * size for _ in T] for _ in T] if joints else None
    words = [[] for _ in range(size)]
    words[0] = [()]
    for k in range(1, n + 1):
        for path in permutations(range(n), k):
            if all(T[a] >> b & 1 for a, b in zip(path, path[1:])):
                s = sum(1 << v for v in path)
                words[s].append(path)
                P[s] += 1
                A[path[0]][s] += 1
                E[path[-1]][s] += 1
                if joints:
                    Q[path[0]][path[-1]][s] += 1
    return P, A, E, Q, words


def mul(A, B):
    C = [0] * len(A)
    for m in range(len(A)):
        s = m
        while True:
            C[m] += A[s] * B[m ^ s]
            if not s:
                break
            s = (s - 1) & m
    return C


def cf(A, B):
    m = len(A) - 1
    return sum(A[s] * B[m ^ s] for s in range(m + 1))


def insert(path, before, vertices):
    i = path.index(before)
    return path[:i] + vertices + path[i:]


def legal(T, path):
    return len(path) == len(set(path)) and all(T[a] >> b & 1 for a, b in zip(path, path[1:]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = ap.parse_args()
    marked = graphcuts = hp_images = cover_images = 0
    finite = []
    for n in range(1, 5):
        nc = nm = 0
        for code in range(1 << (n * (n - 1) // 2)):
            T = graph(n, code)
            P, A, E, Q, words = tables(T)
            oldfull = (1 << n) - 1
            P2 = mul(P, P)
            R = reach(T)
            for B in range(1 << n):
                exits = [b for b in range(n) if B >> b & 1]
                X, U = reset(T, B)
                Pn, An, En, _, _ = tables(U, joints=False)
                AB = [sum(A[b][s] for b in exits) for s in range(1 << n)]
                EB = [sum(E[b][s] for b in exits) for s in range(1 << n)]
                predH = P[oldfull] + 2*cf(P, AB) - cf(EB, AB)
                need(predH == Pn[-1], 'exact multi-exit Hamilton operator')
                need(cf(Pn, An[n])-Pn[-1] == cf(EB, AB), 'new x start gap')
                need(cf(Pn, An[n+1])-Pn[-1] == P[oldfull]-EB[oldfull], 'new y start gap')
                need(all((U[i] & oldfull) == T[i] for i in range(n)), 'old relation retained')
                if strong(T):
                    need(strong(U) == bool(B), 'strong final exactly for nonempty exit set')
                    need(strong(X) == (0 < B < oldfull), 'strong first birth proper nonempty cut')
                if B:
                    image = set()
                    for path in words[oldfull]:
                        first = next(b for b in path if B >> b & 1)
                        lifts = [path + (n+1,n), insert(path, first, (n,)) + (n+1,)]
                        lifts += [insert(path, b, (n+1,n)) for b in exits]
                        for out in lifts:
                            need(legal(U, out) and len(out) == n+2, 'literal Hamilton lift')
                            need(tuple(v for v in out if v<n) == path and out not in image, 'Hamilton lift injection')
                            image.add(out)
                    need(len(image) == (len(exits)+2)*P[oldfull] <= Pn[-1], 'all-order Hamilton lower bound control')
                    hp_images += len(image)
                else:
                    need(Pn[-1] == P[oldfull], 'empty-exit Hamilton boundary')
                for r in range(n):
                    F, Fn = cf(E[r], P), cf(En[r], Pn)
                    QB = [sum(Q[b][r][s] for b in exits) for s in range(1 << n)]
                    predF = F + 2*(cf(mul(E[r], P), AB)+cf(P2,QB)) - cf(mul(E[r], EB), AB) - cf(mul(P, EB), QB)
                    need(predF == Fn, 'exact multi-exit retained-root operator')
                    if B:
                        images, split = set(), 0
                        for s in range(1 << n):
                            for W in words[s]:
                                if not W or W[-1] != r:
                                    continue
                                for V in words[oldfull ^ s]:
                                    BW = [b for b in W if B >> b & 1]
                                    BV = [b for b in V if B >> b & 1]
                                    lifts = [(W,V+(n+1,n))]
                                    if BV:
                                        lifts += [(W,insert(V,BV[0],(n,))+(n+1,)),(W,insert(V,BV[0],(n+1,n)))]
                                    else:
                                        lifts += [(insert(W,BW[0],(n,)),V+(n+1,)),(insert(W,BW[0],(n+1,n)),V)]
                                    if BW and BV:
                                        split += 1
                                        lifts += [(insert(W,BW[0],(n,)),V+(n+1,)),(insert(W,BW[0],(n+1,n)),V)]
                                    for out in lifts:
                                        need(legal(U,out[0]) and legal(U,out[1]) and out[0][-1]==r,
                                             'literal rooted-cover lift')
                                        need(set(out[0]).isdisjoint(out[1]) and set(out[0])|set(out[1])==set(range(n+2)), 'rooted-cover full partition')
                                        need(tuple(v for v in out[0] if v<n)==W and tuple(v for v in out[1] if v<n)==V and out not in images,
                                             'rooted-cover lift injection')
                                        images.add(out)
                        need(len(images)==3*F+2*split <= Fn, 'joint split-cover lower bound')
                        if len(exits)>=2 and any(R[b]>>r&1 for b in exits):
                            need(split>=1 and Fn>=3*F+2, 'reachable multi-exit strictness')
                        if not any(R[b]>>r&1 for b in exits):
                            need(Pn[-1]*F == Fn*P[oldfull], 'unreachable SCC factorization ratio')
                        if Fn == 3*F and Pn[-1] > 3*P[oldfull]:
                            need(len(exits)==1, 'all-tournament positive exact-tripling reduction')
                        cover_images += len(images)
                    else:
                        need(Fn==F, 'empty-exit endpoint boundary')
                    nm += 1
                nc += 1
        graphcuts += nc
        marked += nm
        finite.append({'n':n,'graph_exit_cuts':nc,'marked_cuts':nm})

    probe = []
    for n in range(3,6):
        strong_inputs = exit_sets = hits2 = hits4 = 0
        for code in range(1 << (n*(n-1)//2)):
            T = graph(n,code)
            if not strong(T):
                continue
            strong_inputs += 1
            P,A,E,_,_ = tables(T,joints=False)
            S=[cf(P,a) for a in A]
            EA=[[cf(e,a) for a in A] for e in E]
            for B in range(1 << n):
                V=[b for b in range(n) if B>>b&1]
                if len(V)<2:
                    continue
                exit_sets += 1
                h=P[-1]+2*sum(S[b] for b in V)-sum(EA[b][c] for b in V for c in V)
                need(h >= (len(V)+2)*P[-1], 'complete strong5 scalar probe lower bound')
                hits2 += h-3*P[-1] == 2
                hits4 += h-3*P[-1] == 4
        probe.append({'n':n,'strong_inputs':strong_inputs,'multi_exit_sets':exit_sets,'plus2':hits2,'plus4':hits4})
    need(probe[0]['plus4']==6 and all(r['plus2']==0 for r in probe) and all(r['plus4']==0 for r in probe[1:]), 'literal strong small probe')
    controls=[]
    for name,T,B in [('strong_C3_plus4',(2,4,1),3),('weak_TT2_plus2',(2,0),3),('empty_cut',(2,4,1),0)]:
        X,U=reset(T,B);P,_,E,_,_=tables(T,joints=False);Pn,_,En,_,_=tables(U,joints=False)
        controls.append({'name':name,'source':T,'B':B,'first':X,'target':U,'H':P[-1],'H_new':Pn[-1],
                         'F':[cf(E[r],P) for r in range(len(T))],'F_new':[cf(En[r],Pn) for r in range(len(T))]})
    need(controls[0]['H_new']==13 and controls[1]['H_new']==5, 'sharp scalar boundary examples')
    small=tuple(controls[0]['target'])
    P,A,E,_,_=tables(small,joints=False)
    sigma=[cf(P,a)-P[-1] for a in A]
    need(sigma==[5,5,13,2,1], 'small multi-exit child supplies both scalar letters')
    small_branches=[]
    for a,d in ((3,4),(4,2)):
        _,child=reset(small,1<<a)
        Pc,Ac,Ec,_,_=tables(child,joints=False)
        gaps=[cf(Pc,v)-Pc[-1] for v in Ac]
        need(Pc[-1]==3*P[-1]+d and strong(child), 'literal scalar branch after small multiport positive')
        need(all(g not in (1,2) for g in gaps), 'literal reset-image stopping control')
        small_branches.append({'a':a,'digit':d,'rows':child,'H':Pc[-1],'F_oldroot0':cf(Ec[0],Pc),'start_gaps':gaps})
    result={'scope':'exact all-order formulas/injections have a separate proof; finite controls all labelled inputs<=4 and exits/roots, plus full strong<=5 scalar probe',
            'gates':GATES,'graph_exit_cuts':graphcuts,'marked_cuts':marked,'hamilton_lift_images':hp_images,
            'cover_lift_images':cover_images,'operator_census':finite,'strong_scalar_probe':probe,'controls':controls,
            'small_positive_gaps':sigma,'small_positive_branches':small_branches}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();args.output.write_bytes(raw)
    print(json.dumps({'gates':GATES,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':
    main()
