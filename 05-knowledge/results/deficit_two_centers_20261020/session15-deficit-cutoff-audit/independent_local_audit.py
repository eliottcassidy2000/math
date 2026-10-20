#!/usr/bin/env python3
"""Independent exact audit. No imports from the candidate implementation.

Represent a tournament by an adjacency matrix. Detect strongness by absence
of a dominating cut, independently certify strongness by transitive closure,
and search directly for directed four-cycles with three allowable exits.
"""
from itertools import combinations, permutations
from collections import Counter
from pathlib import Path
import hashlib, json


def matrix(n, code, base=None):
    a = [[False for _ in range(n)] for _ in range(n)]
    if base is not None:
        for i in range(4):
            for j in range(4):
                a[i][j] = bool(base[i][j])
    k = 0
    for i in range(n):
        for j in range(i+1,n):
            if base is not None and j < 4:
                continue
            a[i][j] = bool((code // (2**k)) % 2)
            a[j][i] = not a[i][j]
            k += 1
    return a


def strong_cut(a):
    n = len(a)
    # A tournament is nonstrong iff it has a nonempty proper dominating set.
    for b in range(1,(1<<n)-1):
        if all(a[i][j] for i in range(n) if b>>i&1
               for j in range(n) if not b>>j&1):
            return False
    return True


def strong_closure(a):
    n=len(a)
    r = [[i == j or a[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                r[i][j] = r[i][j] or (r[i][k] and r[k][j])
    return all(all(row) for row in r)


def witnesses(a):
    found=[]
    n=len(a)
    for S in combinations(range(n),4):
        # A directed cycle is itself the strongness certificate.
        C = next(((S[0],)+tail for tail in permutations(S[1:])
                  if all(a[((S[0],)+tail)[k]][((S[0],)+tail)[(k+1)%4]]
                         for k in range(4))), None)
        if C is None:
            continue
        for v in range(n):
            if v not in S:
                indegree=sum(a[u][v] for u in S)
                if indegree >= 3:
                    found.append((C,v,indegree))
    return found


def assert_true(value, description):
    if not value:
        raise RuntimeError(description)


def run_family(name,base):
    n=6
    N=1<<(15 if base is None else 9)
    counts=Counter()
    histogram=Counter()
    digest=hashlib.sha256()
    first=[]
    for code in range(N):
        a=matrix(n,code,base)
        s=strong_cut(a)
        assert_true(s==strong_closure(a),f'strong tests disagree: {name} {code}')
        counts['strong' if s else 'nonstrong']+=1
        if not s:
            continue
        ww=witnesses(a)
        assert_true(bool(ww),f'missing witness: {name} {code}')
        histogram[len(ww)]+=1
        if any(w[2]==4 for w in ww):
            counts['has_four_exits']+=1
        C,v,d=ww[0]
        # Verify the three/four rooted cycle words end at the claimed exits.
        words=[]
        for j,u in enumerate(C):
            if a[u][v]:
                word=C[j+1:]+C[:j+1]
                assert_true(all(a[word[k]][word[k+1]] for k in range(3)), 'bad rooted path')
                assert_true(a[word[-1]][v], 'bad exit')
                words.append(word)
        assert_true(len(words)>=3,'too few paths')
        for P,Q in combinations(words,2):
            rel={P[i]:Q[i] for i in range(4)}
            curr={x:x for x in P}
            order=0
            while True:
                curr={x:rel[curr[x]] for x in P}; order+=1
                if all(curr[x]==x for x in P):
                    break
            assert_true(order in (2,4),'relative map not even')
        line=[code,list(C),v,d]
        digest.update((json.dumps(line,separators=(',',':'))+'\n').encode())
        if len(first)<5:first.append(line)
    return dict(cases=N,counts=dict(counts),witness_count_histogram=dict(sorted(histogram.items())),
                first_witness_records=first,canonical_witness_sha256=digest.hexdigest())


def regular_connectivity_controls():
    ans={}
    for n in (1,3,5,7):
        regular=Counter()
        if n==7:
            # All labelled circulants, not all 2^21 tournaments.
            codes=[]
            for bits in range(1<<((n-1)//2)):
                S={j+1 if bits>>j&1 else n-j-1 for j in range((n-1)//2)}
                a=[[i!=j and (j-i)%n in S for j in range(n)] for i in range(n)]
                codes.append(a)
        else:
            codes=(matrix(n,k) for k in range(1<<(n*(n-1)//2)))
        for a in codes:
            if any(sum(row)!=(n-1)//2 for row in a):continue
            regular['tournaments']+=1
            for s in range((n+2)//3):
                for D in combinations(range(n),s):
                    rem=[i for i in range(n) if i not in D]
                    b=[[a[i][j] for j in rem] for i in rem]
                    assert_true(strong_closure(b),f'regular connectivity failed n={n}')
                    regular['safe_deletions']+=1
        ans[str(n)]=dict(regular)
    return ans


def main():
    # Independent labelled S4 orientation: 0->1->2->3->0, 0->2, 1->3.
    base=[[False,True,True,False],
          [False,False,True,True],
          [False,False,False,True],
          [True,False,False,False]]
    assert_true(strong_closure(base),'base not strong')
    ans={}
    for name,b in [('all_labelled_six',None),('fixed_s4_extensions',base)]:
        ans[name]=run_family(name,b)
        print(name,ans[name]['counts'],flush=True)
    ans['regular_connectivity_controls']=regular_connectivity_controls()
    ans['implementation_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=Path(__file__).with_name('independent_local_results.json')
    out.write_text(json.dumps(ans,indent=2,sort_keys=True)+'\n')
    print('all checks passed')
if __name__=='__main__':main()
