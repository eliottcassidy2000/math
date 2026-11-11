#!/usr/bin/env python3
"""Small fresh hostile controls, not an audit of historical results."""
import itertools,json,math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from census import census,direct_profile,require


def graph(n,mask):
    t=[[0]*n for _ in range(n)]
    for i,(a,b) in enumerate(itertools.combinations(range(n),2)):
        t[a][b]=(mask>>i)&1;t[b][a]=1-t[a][b]
    return t


def induced(t,vs):return [[t[a][b] for b in vs] for a in vs]


def in_family(t):
    n=len(t)
    if n<=3:return True
    for mask in range(1,(1<<n)-1):
        vs=[i for i in range(n) if mask>>i&1]
        if not 1<len(vs)<n:continue
        outside=[i for i in range(n) if not mask>>i&1]
        if all(len({t[v][a] for a in vs})==1 for v in outside):
            # Collapse a module. Substitution closure with hereditary prime
            # restrictions makes any such decomposition a valid recursion.
            return in_family(induced(t,vs)) and in_family(induced(t,outside+[vs[0]]))
    return False


def canonical(t):
    n=len(t)
    return min(sum(t[p[a]][p[b]]<<i for i,(a,b) in enumerate(itertools.combinations(range(n),2))) for p in itertools.permutations(range(n)))


def main():
    A,S,rows,cycle=census(8)
    checks=[]
    for n in range(1,6):
        labels=defaultdict(int);classes=defaultdict(set)
        for mask in range(1<<(n*(n-1)//2)):
            t=graph(n,mask)
            if not in_family(t):continue
            q=direct_profile(t);labels[q]+=1;classes[q].add(canonical(t))
        require(set(classes)==set(A[n]),'independent response support mismatch')
        for q,(m,u) in A[n].items():
            require(len(classes[q])==m,'independent ordinary count mismatch')
            require(Fraction(labels[q],math.factorial(n))==u,'independent labelled mass mismatch')
        checks.append(dict(n=n,orientation_universe=1<<(n*(n-1)//2),family_labeled=sum(labels.values()),classes=sum(map(len,classes.values())),response_types=len(classes)))
    collision=None
    for n in range(1,9):
        pairs=list(A[n])
        for i,a in enumerate(pairs):
            for b in pairs[i+1:]:
                if a[0]==b[0] and a!=b:
                    ca=cycle(a,(1,),(1,));cb=cycle(b,(1,),(1,))
                    if ca[0]!=cb[0]:
                        collision=dict(n=n,H=a[0],profile_A=a,profile_B=b,context='C3[child,1,1]',context_H_A=ca[0],context_H_B=cb[0]);break
            if collision:break
        if collision:break
    require(collision is not None,'missing H-only hostile')
    out=dict(status='FINITE-EXACT',small_labelled_universe=checks,H_only_collision=collision)
    p=Path(__file__).with_name('controls.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
