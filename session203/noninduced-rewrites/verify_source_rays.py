import sympy as S
from itertools import product
import json
M={ 'O':S.Matrix([[1,0,1,0],[1,2,1,2],[0,0,0,0],[0,0,1,1]]), 'Z':S.Matrix([[2,1,1,1],[0,1,0,0],[0,0,1,1],[0,0,0,0]]), 'J':S.Matrix([[0,1,0,1],[0,0,0,0],[1,0,1,1],[0,0,1,0]])}
F=S.Matrix([[2,1,3,3]])
def req(ok,msg):
    if not ok:raise RuntimeError(msg)
for c,A in M.items():req(all(x>=0 for x in 2*F*A-3*F),'growth '+c)
report={}
for q in (3,5,7,9):
    L=0
    while 3**(L+1)<=q*2**(L+1):L+=1
    rays={};words=0
    for n in range(1,L+1):
        for letters in product('OZJ',repeat=n):
            words+=1;A=S.eye(4)
            for c in letters:A=M[c]*A
            N=(A-q*S.eye(4)).nullspace()
            if not N:continue
            req(len(N)==1,'multidimensional kernel needs cone intersection')
            v=N[0];den=S.ilcm(*[x.q for x in v]);v=tuple(int(den*x) for x in v)
            if all(x<=0 for x in v):v=tuple(-x for x in v)
            if not all(x>=0 for x in v):continue
            g=S.igcd(*v);v=tuple(int(x/g) for x in v)
            req(sum(v)>0,'zero ray')
            key=','.join(map(str,v));mat=','.join(map(str,A))
            rays.setdefault(key,{}).setdefault(mat,[]).append(''.join(letters))
    maxcount=max((len(d) for d in rays.values()),default=0)
    req(maxcount<q,'insufficient obstruction')
    report[q]={'length_bound':L,'words_checked':words,'max_distinct_matrices_per_ray':maxcount,'rays':[{'ray':list(map(int,v.split(','))),'matrix_word_classes':list(D.values())} for v,D in sorted(rays.items())]}
print(json.dumps(report,indent=2))
