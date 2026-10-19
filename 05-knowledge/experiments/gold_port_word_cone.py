#!/usr/bin/env python3
"""Exact fixed Gold strategies for the five-seam B8/B9 word family."""
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from functools import reduce
from pathlib import Path
import json
from gold_terminal_chain import count, extensions, require, deficits
from gold_four_collar import table_pred

CORES={'B9':((0,0,1,0,8,24,3,5,189),tuple(range(9))),
       'B8':((0,1,0,4,12,1,3,95),(0,2,3,4,5,6,7,8))}
GRAPHS=[((6,0),(8,3)),((6,3),(8,0)),((8,0),(8,3)),((6,0),(8,0),(8,3)),((6,3),(8,0),(8,3))]
ALPHA=[(),(6,),(8,),(6,8),(8,6)]
BETA=[(),(0,),(3,),(0,3),(3,0)]
PAIRS={'B9':((5,7),(1,2)),'B8':((5,7),(6,7))}

def shuffles(a,b):
    for pos in combinations(range(len(a)+len(b)),len(a)):
        pos=set(pos);i=j=0;w=[]
        for k in range(len(a)+len(b)):
            if k in pos:w.append(a[i]);i+=1
            else:w.append(b[j]);j+=1
        yield tuple(w)

def skeleton(name):
    for s in shuffles((0,2,7),(3,4,5)):
        s=s+(8,)
        if name=='B8':
            for i in range(s.index(0)+1,8):yield s[:i]+(6,)+s[i:]
        else:
            for i,j in combinations(range(9),2):
                w=[];k=0
                for q in range(9):
                    if q==i:w.append(1)
                    elif q==j:w.append(6)
                    else:w.append(s[k]);k+=1
                if w.index(0)<w.index(6):yield tuple(w)

def profiles(w):
    b=();a=()
    for x in w:
        if x not in (0,3):break
        b+=(x,)
    for x in reversed(w):
        if x not in (6,8):break
        a=(x,)+a
    require(a,'empty suffix')
    return BETA.index(b),ALPHA.index(a)-1

def joint(name,words):
    J=[[[0]*4 for _ in range(5)] for _ in range(4)];e,f=PAIRS[name]
    for w in words:
        b,a=profiles(w);i=int(w.index(e[0])>w.index(e[1]));j=int(w.index(f[0])>w.index(f[1]));J[2*i+j][b][a]+=1
    return J

def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def scalar(l,M,r):return sum(l[i]*M[i][j]*r[j] for i in range(len(l)) for j in range(len(r)))
def primitive(x):
    d=reduce(gcd,x);return tuple(v//d for v in x) if d else None

def kernels():
    return [[[sum(all((u,v) in G or w.index(u)<w.index(v) for u in a for v in b) for w in shuffles(a,b)) for b in BETA] for a in ALPHA] for G in GRAPHS]

def cones(K,N):
    rawL=sorted({tuple(x) for k in K for x in k});rawR=sorted({tuple(k[i][j] for i in range(1,5)) for k in K for j in range(5)})
    orbL={(1,)*5};orbR={(1,)*4}
    for k in K:
        orbL.update(primitive(x) for x in (k[1],k[2],[k[3][j]+k[4][j] for j in range(5)]))
        orbR.update(primitive([k[i][j] for i in range(1,5)]) for j in range(3))
        orbR.add(primitive([k[i][3]+k[i][4] for i in range(1,5)]))
    bufL={(1,)*5};bufR={(1,)*4}
    for k in K:
        for M in N.values():
            P=mm(M,k[1:]);Q=mm(k[1:],M)
            bufL.update(primitive(x) for x in P if any(x))
            bufR.update(primitive([row[j] for row in Q]) for j in range(4) if any(row[j] for row in Q))
    return {'raw':(rawL,rawR),'orbit':(sorted(orbL),sorted(orbR)),'buffer':(sorted(bufL),sorted(bufR))}

def assemble(names,gs):
    p=[];offset=[]
    for name in names:
        base,labels=CORES[name];off=len(p);offset.append(off);p += [(x<<off)|((1<<off)-1) for x in base]
    for q,g in enumerate(gs):
        la=CORES[names[q]][1];lb=CORES[names[q+1]][1]
        for u,v in GRAPHS[g]:p[offset[q+1]+lb.index(v)]&=~(1<<(offset[q]+la.index(u)))
    return tuple(p),offset

def main():
    Js={};N={}
    for name,(pred,labels) in CORES.items():
        words=[tuple(labels[x] for x in w) for w in extensions(pred)];sk=list(skeleton(name))
        require(len(sk)==len(set(sk)) and set(sk)==set(words),'independent skeleton count')
        J=joint(name,words);require(J==joint(name,sk),'joint profile count')
        Js[name]=J;N[name]=[[sum(t[i][j] for t in J) for j in range(4)] for i in range(5)]
        aggregate=[[sum(map(sum,J[2*i+j])) for j in range(2)] for i in range(2)]
        e,f=PAIRS[name];got=table_pred(pred,tuple(labels.index(x) for x in e),tuple(labels.index(x) for x in f))
        require(aggregate==got and sum(map(sum,got))==count(pred),'independent core DP')
    K=kernels();C=cones(K,N);results={}
    expected={'raw':{'B9':Fraction(281,2946),'B8':Fraction(57,625)},'orbit':{'B9':Fraction(279,2309),'B8':Fraction(7,60)},'buffer':{'B9':Fraction(23365,173249),'B8':Fraction(791,6015)}}
    for level,(L,R) in C.items():
        common=min(expected[level].values());entry={'left_rays':L,'right_rays':R,'common_margin':str(common),'cores':{}}
        for name,J in Js.items():
            worst=None;checks=[];row_min=[]
            for li,l in enumerate(L):
                local=[]
                for ri,r in enumerate(R):
                    flat=[scalar(l,M,r) for M in J];table=[flat[:2],flat[2:]];Z=sum(flat);D=sum(deficits(table),[])
                    slack=min(Fraction(x,Z) for x in D)
                    if worst is None or slack<worst['margin']:worst=dict(margin=slack,left_index=li,right_index=ri,table=table,deficits=D)
                    scaled=[common.denominator*x-common.numerator*Z for x in D]
                    require(min(scaled)>=0,'cone Gold margin');local+=scaled
                    checks.append(dict(left=li,right=ri,total=Z,deficits=D,scaled_minimum=min(scaled)))
                row_min.append(min(local))
            require(worst['margin']==expected[level][name],'sharp cone minimum')
            worst['margin']=str(worst['margin']);entry['cores'][name]={'worst':worst,'row_minimum_scaled_deficits':row_min,'checks':checks}
        results[level]=entry
    require([len(C[x][0])*len(C[x][1]) for x in ('raw','orbit','buffer')]==[48,48,540],'cone sizes')
    actual=[]
    for names,gs in ((('B9','B9'),(0,)),(('B8','B9'),(2,)),(('B8','B9','B8','B9'),(0,3,1)),(('B9','B8','B9','B8','B9','B9'),(0,1,2,3,4))):
        pred,off=assemble(names,gs);total=count(pred);left=[(1,)*5]
        for i in range(len(names)-1):left.append(tuple(mm(mm([left[-1]],N[names[i]]),K[gs[i]][1:])[0]))
        right=[None]*len(names);right[-1]=(1,)*4
        for i in range(len(names)-2,-1,-1):right[i]=tuple(row[0] for row in mm(mm(K[gs[i]][1:],N[names[i+1]]),[[x] for x in right[i+1]]))
        for i,name in enumerate(names):
            lab=CORES[name][1];e,f=PAIRS[name]
            got=table_pred(pred,tuple(off[i]+lab.index(x) for x in e),tuple(off[i]+lab.index(x) for x in f))
            flat=[scalar(left[i],M,right[i]) for M in Js[name]];expected_table=[flat[:2],flat[2:]]
            require(got==expected_table and sum(flat)==total,'actual word independent DP')
            require(min(Fraction(x,total) for x in sum(deficits(got),[]))>=Fraction(791,6015),'actual word uniform margin')
            actual.append(dict(word=list(names),seams=list(gs),block=i,table=got,total=total))
    pred,off=assemble(('B9','B9'),(0,));old=table_pred(pred,(off[1]+5,off[1]+7),(off[1]+4,off[1]+7));require(min(sum(deficits(old),[]))<0,'old triangular selector obstruction')
    out=dict(status='proved fixed-strategy Gold upgrade via exact port-word cones',strategies=PAIRS,joint_profile_tables=Js,core_totals=N,kernels=K,cone_certificates=results,actual_word_checks=actual,old_triangular_failure=dict(table=old,deficits=deficits(old)))
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
