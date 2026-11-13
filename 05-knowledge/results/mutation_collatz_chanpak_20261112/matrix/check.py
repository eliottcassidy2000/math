#!/usr/bin/env python3
"""Integer congruence and actual small tournament controls; no assert gates."""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def matrix(n,mask):
    B=[[0]*n for _ in range(n)]
    for bit,(i,j) in enumerate(combinations(range(n),2)):
        B[i][j]=1 if mask>>bit&1 else -1
        B[j][i]=-B[i][j]
    return B

def trans(A):return list(map(list,zip(*A)))
def mul(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def ident(n):return [[int(i==j) for j in range(n)] for i in range(n)]

def mutate(B,k):
    n=len(B)
    return [[-B[i][j] if i==k or j==k else
             B[i][j]+(abs(B[i][k])*B[k][j]+B[i][k]*abs(B[k][j]))//2
             for j in range(n)] for i in range(n)]

def congruence(B,k):
    E=ident(len(B))
    for i in range(len(B)):E[i][k]=-1 if i==k else max(B[i][k],0)
    need(mul(E,E)==ident(len(B)),'involutive integer basis')
    return mul(mul(E,B),trans(E))

def det(B):
    if not B:return 1
    A=[r[:] for r in B]; n=len(A); old=1; sign=1
    for k in range(n-1):
        if A[k][k]==0:
            p=next((i for i in range(k+1,n) if A[i][k]),None)
            if p is None:return 0
            A[k],A[p]=A[p],A[k];sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=A[i][j]*pivot-A[i][k]*A[k][j]
                need(numerator%old==0,'Bareiss exact division')
                A[i][j]=numerator//old
        for i in range(k+1,n):A[i][k]=0
        old=pivot
    return sign*A[-1][-1]

def pf(B):
    if not B:return 1
    if len(B)%2:return 0
    return sum((-1)**(j+1)*B[0][j]*
               pf([[B[a][b] for b in range(len(B)) if b not in (0,j)]
                   for a in range(len(B)) if a not in (0,j)])
               for j in range(1,len(B)))

def rank_mod(B,p):
    A=[[x%p for x in r] for r in B]; rank=0
    for j in range(len(B)):
        i=next((i for i in range(rank,len(B)) if A[i][j]),None)
        if i is None:continue
        A[i],A[rank]=A[rank],A[i]
        inv=pow(A[rank][j],-1,p);A[rank]=[inv*x%p for x in A[rank]]
        for i in range(len(B)):
            if i!=rank:
                f=A[i][j];A[i]=[(x-f*y)%p for x,y in zip(A[i],A[rank])]
        rank+=1
    return rank

def join(A,B):
    a,b=len(A),len(B)
    return [r+[1]*b for r in A]+[[-1]*a+r for r in B]

def hp(B):
    n=len(B)
    if not n:return 1
    @lru_cache(None)
    def end(mask,j):
        rest=mask^(1<<j)
        if not rest:return 1
        return sum(end(rest,i) for i in range(n) if rest>>i&1 and B[i][j]==1)
    return sum(end((1<<n)-1,j) for j in range(n))

def rows_matrix(rows):
    n=len(rows)
    need(all(not (rows[i]>>i&1) for i in range(n)),'loopless')
    need(all(((rows[i]>>j&1)+(rows[j]>>i&1))==1 for i,j in combinations(range(n),2)),'tournament')
    return [[0 if i==j else 1 if rows[i]>>j&1 else -1 for j in range(n)] for i in range(n)]

def cover_end(B,r):
    n=len(B); out=0
    for mask in range(1<<n):
        if not mask>>r&1:continue
        before=[i for i in range(n) if i!=r and mask>>i&1]
        rest=[i for i in range(n) if not mask>>i&1]
        e=0
        for p in permutations(before):
            path=p+(r,)
            e+=all(B[path[i]][path[i+1]]==1 for i in range(len(path)-1))
        out+=e*hp([[B[i][j] for j in rest] for i in rest])
    return out

def main():
    mutations=0
    for n in range(1,5):
        for mask in range(1<<(n*(n-1)//2)):
            B=matrix(n,mask)
            for k in range(n):
                C=B
                for pivot in (k,(k+1)%n,k):
                    D=mutate(C,pivot)
                    need(D==congruence(C,pivot),'mutation congruence')
                    need(det(C)==det(D),'determinant conservation')
                    C=D; mutations+=1
    triangle=matrix(3,5) # 0->1->2->0
    triples=0
    for r in range(5):
        for mask in range(1<<(r*(r-1)//2)):
            R=matrix(r,mask);B=join(triangle,R);n=r+3
            P=ident(n)
            P[0][2]=P[1][2]=1
            for j in range(3,n):P[0][j]=1;P[1][j]=-1
            transformed=mul(mul(trans(P),B),P)
            D=[[0]*n for _ in range(n)];D[0][1]=1;D[1][0]=-1
            for i in range(r):
                D[2][i+3]=3;D[i+3][2]=-3
                for j in range(r):D[i+3][j+3]=R[i][j]
            need(det(P)==1 and transformed==D,'triangle integral normal form')
            need(n-rank_mod(B,3)==1+r-rank_mod(R,3),'triangle ternary nullity')
            need(det(B)==9*det(join([[0]],R)),'triangle determinant law')
            need(hp(B)==3*hp(R),'independent path concatenation')
            triples+=1
    records=[]
    for name,rows in (('two_root',(52,33,34,39,46,0)),
                      ('three_root',(58,60,57,16,32,8))):
        B=rows_matrix(rows)
        literal=sum(all(B[p[i]][p[i+1]]==1 for i in range(5)) for p in permutations(range(6)))
        need(hp(B)==literal==9,'actual H9 controls')
        need(det(B)==pf(B)**2,'independent determinant')
        roots=[j for j in range(6) if sum(B[i][j]==1 for i in range(6))==1]
        records.append({'name':name,'rows':rows,'H':literal,'det':det(B),
                        'pfaffian':pf(B),'rank_mod3':rank_mod(B,3),'indegree_one':roots})
    need([r['det'] for r in records]==[25,81],'different integral classes')
    need([r['indegree_one'] for r in records]==[[3,4],[0,1,2]],'root counts')
    prime=[]
    for name,rows,expected in (
        ('after_5_to_17',(10,28,1,4,45,15),(49,[4,5],[0,3])),
        ('ready_17_chart',(58,60,1,52,36,4),(1,[0,1],[0,1]))):
        B=rows_matrix(rows)
        literal=sum(all(B[p[i]][p[i+1]]==1 for i in range(5)) for p in permutations(range(6)))
        roots=[j for j in range(6) if sum(B[i][j]==1 for i in range(6))==1]
        defects=[cover_end(B,r)-literal for r in roots]
        need(hp(B)==literal==17,'literal prime17 count')
        need((det(B),roots,defects)==expected,'prime readiness separation')
        need(pf(B)**2==det(B),'prime Pfaffian')
        prime.append({'name':name,'rows':rows,'H':literal,'det':det(B),
                      'rank_mod7':rank_mod(B,7),'indegree_one':roots,'defects':defects})
    A=[[0,1,1],[-1,0,1],[-1,-1,0]]
    B=[[0,-2,1],[2,0,1],[-1,-1,0]]
    need([[x%3 for x in r] for r in A]==[[x%3 for x in r] for r in B],'same input residues')
    need((mutate(A,1)[0][2],mutate(B,1)[0][2])==(2,1),'residue-only update hostile')
    source='05-knowledge/results/triples_ports_sieve_20261112/ports/PROOF.md'
    print(json.dumps({'status':'PASS','scope':'explicit finite controls; identities supply all orders',
        'mutation_checks':mutations,'residual_tournaments':triples,'H9_controls':records,
        'prime17_controls':prime,
        'residue_hostile':{'modulus':3,'pivot':1,'new_02_entries':[2,1]},
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256':{source:sha256((ROOT/source).read_bytes()).hexdigest()}},
        indent=2,sort_keys=True))

if __name__=='__main__':main()
