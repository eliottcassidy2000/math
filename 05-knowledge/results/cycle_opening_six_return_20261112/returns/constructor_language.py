"""Exact physical indegree-one L/T language and neutral macro identities."""
from pathlib import Path
from itertools import combinations, product
import json, gzip
import probe
HERE=Path(__file__).resolve().parent
a=probe.a
L=((1,2,0),(0,-1,1),(1,0,1))
T=((2,2,0),(0,0,1),(2,2,1))
I=((1,0,0),(0,1,0),(0,0,1))

def need(ok,label):
    if not ok: raise RuntimeError(label)
def mm(A,B): return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mv(A,v): return tuple(sum(A[i][k]*v[k] for k in range(3)) for i in range(3))
def word_matrix(word):
    ans=I
    for c in word: ans=mm(dict(L=L,T=T)[c],ans)
    return ans

def birth(rows,port):
    n=len(rows); out=list(rows)+[((1<<n)-1)^(1<<port)]
    out[port] |= 1<<n
    return tuple(out)

def construct(rows,root,word):
    for letter in word:
        if letter=='L': rows=birth(rows,root)
        else:
            rows=birth(rows,root)
            rows=birth(rows,root)
        root=len(rows)-1
    return rows,root

def state(rows,root):
    h,s,F,_=probe.profiles(rows)
    full=(1<<len(rows))-1
    d=a.dp(a.induced(rows,[i for i in range(len(rows)) if i!=root]))[0][-1]
    need((full^rows[root]^(1<<root)).bit_count()==1,'indegree-one interface')
    b=F[root]-h
    need(s[root]==d-b,'indegree-one cover identity')
    return h-d,b,d

def strong(rows):
    n=len(rows); reach=1
    while True:
        new=reach
        for i in range(n):
            if reach>>i&1:new |= rows[i]
        if new==reach:break
        reach=new
    if reach!=(1<<n)-1:return False
    full=(1<<n)-1
    rev=tuple(full^row^(1<<i) for i,row in enumerate(rows))
    reach=1
    while True:
        new=reach
        for i in range(n):
            if reach>>i&1:new |= rev[i]
        if new==reach:return reach==full
        reach=new

def main():
    records=[json.loads(x) for x in gzip.open(probe.RESULTS/'metagraph_aspect_sequences_20261111/records.jsonl.gz','rt') if json.loads(x)['n']<=5]
    checks=0; controls=[]
    for rec in records:
        rows=tuple(rec['rows']);n=len(rows);full=(1<<n)-1
        for root in range(n):
            if (full^rows[root]^(1<<root)).bit_count()!=1:continue
            v=state(rows,root)
            for word in ['L','T','LL','TT']:
                out,r=construct(rows,root,word)
                need(state(out,r)==mv(word_matrix(word),v),'literal constructor '+word)
                checks+=1
            if len(controls)<5 or (n==5 and rec['strong'] and len(controls)<8):controls.append((rows,root,v))
    A=word_matrix('TTTLL'); B=word_matrix('TTLTL')
    need(A==((48,48,24),(18,18,9),(54,54,27)),'A matrix')
    need(B==((44,44,22),(24,24,12),(58,58,29)),'B matrix')
    need(mm(A,mm(B,A))==mm(A,mm(A,B)),'ABA=BAA complete state')
    need(mm(B,mm(B,A))==mm(B,mm(A,B)),'ABB=BAB complete state')
    literal=[]
    for rows,root,v in controls:
        left,lr=construct(rows,root,'TTTLL');right,rr=construct(rows,root,'TTLTL')
        lv,rv=state(left,lr),state(right,rr);q=2*v[0]+2*v[1]+v[2]
        need(lv==(24*q,9*q,27*q) and rv==(22*q,12*q,29*q),'literal 51q macros')
        n=len(rows); differences=[(i,j) for i,j in combinations(range(n+8),2) if ((left[i]^right[i])>>j)&1]
        expected=[(n+3,n+5),(n+4,n+5),(n+4,n+6),(n+5,n+6)]
        need(differences==expected,'four precise arc differences')
        need(lr==rr==n+7,'same literal final root')
        need(((1<<(n+8))-1)^left[lr]^(1<<lr)==1<<(n+6),'same final sole predecessor')
        need(((1<<(n+8))-1)^right[rr]^(1<<rr)==1<<(n+6),'same final sole predecessor')
        need(not strong(rows) or (strong(left) and strong(right)),'strong source remains strong')
        literal.append(dict(rows=rows,root=root,input_state=v,q=q,H=51*q,left_state=lv,right_state=rv,changed_arcs=differences,strong_input=strong(rows),strong_outputs=[strong(left),strong(right)]))
    census=[]
    def words(cost):
        if cost==0:yield ''
        if cost>=1:
            for w in words(cost-1):yield w+'L'
        if cost>=2:
            for w in words(cost-2):yield w+'T'
    for cost in range(1,9):
        bank={};collisions=[]
        for word in words(cost):
            M=word_matrix(word);row=tuple(M[0][i]+M[2][i] for i in range(3))
            for prev,PM in bank.get(row,[]):
                collisions.append(dict(words=[prev,word],same_state=PM==M,different_defect=PM[1]!=M[1]))
            bank.setdefault(row,[]).append((word,M))
        census.append(dict(cost=cost,words=sum(map(len,bank.values())),H_rows=len(bank),collisions=collisions))
    need(not any(c['different_defect'] for e in census[:-1] for c in e['collisions']),'minimum cost within L/T language')
    out=dict(status='VERIFIED exact symbolic matrices plus independent literal HP checks',interface_checks=checks,controls=literal,A=A,B=B,census=census)
    (HERE/'constructor_language.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS',checks,'literal interfaces;',len(literal),'literal macro pairs; complete symbolic diamonds')
if __name__=='__main__':main()
