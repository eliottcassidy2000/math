"""Independent literal word enumerator for one-module excursion coefficients."""
from tournament_kernel import *
from math import factorial,gcd
from itertools import combinations
import json

def word_counts(a,root=0):
    n=len(a);ext=[i for i in range(n) if i!=root];k=len(ext);counts=[0]*(k//2+2)
    for p in permutations(ext):
        for c in range(1,k//2+2):
            for slots in combinations(range(k+1),c):
                # Returning in one exterior vertex would require a two-cycle.
                if any(v-u<2 for u,v in zip(slots,slots[1:])):continue
                selected=set(slots);word=[]
                for j in range(k+1):
                    if j in selected:word.append(root)
                    if j<k:word.append(p[j])
                if all(a[word[i]][word[i+1]] for i in range(len(word)-1)):
                    counts[c]+=1
    return counts

def detector(c):
    s=c-1;n=2*s+1;a=[[0]*n for _ in range(n)]
    for i,j in combinations(range(n),2):a[i][j]=1
    for i in range(s):
        ai=1+i;bi=1+s+i;a[0][ai]=1;a[ai][0]=0;a[bi][0]=1;a[0][bi]=0
        for j in range(s):
            bj=1+s+j;a[ai][bj]=int(i<=j);a[bj][ai]=int(i>j)
    return tuple(map(tuple,a))

def run():
    gcds=[0]*4;cases=0;rows=[]
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=matrix(n,code);W=word_counts(a)
            for c in range(1,len(W)):
                value=factorial(c)*W[c];d=factorial(c)*factorial(c-1)
                require(value%d==0,'excursion divisor');gcds[c]=gcd(gcds[c],value)
            cases+=1
    require(gcds[1:]==[1,2,12],'coefficient gcds')
    for c in range(1,5):
        a=detector(c);W=word_counts(a)
        require(W[c]==factorial(c-1),'sharp detector')
        require(strong(a),'strong detector')
        rows.append({'c':c,'order':len(a),'code':encoding(a),'W':W[1:],'leading_coefficient':factorial(c)*W[c]})
    moduli=[]
    for M in range(2,33):
        c=1
        while factorial(c)*factorial(c-1)%M:c+=1
        moduli.append({'modulus':M,'first_invisible_cover_count':c})
    print(json.dumps({'scope':'FINITE-EXACT all fixed-root quotients order<=5; detector c<=4; universal proof audited separately',
                      'quotients':cases,'coefficient_gcds':gcds[1:],'detectors':rows,'modular_cutoffs':moduli},sort_keys=True,indent=2))
if __name__=='__main__':run()
