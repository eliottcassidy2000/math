#!/usr/bin/env python3
import itertools,json

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def move(v,c):
    D,E,U,d=v
    if c=='O':return D,0,2*U,d
    if c=='Z':return D+U-E,d,2*U-d,d
    return D+U-E,D-E,D+U,D

words=0; sync=[]; face_count=0
for k in range(11):
    for chars in itertools.product('OZJ',repeat=k):
        w=''.join(chars);v=(1,0,1,0);u=(0,1,0,1)
        for c in w:
            v=move(v,c);u=move(u,c)
            require(v[0]>v[1]>=0 and v[2]>v[3]>=0,'difference cone')
        condition=v[1]==v[3]==0 and v[0]==v[2]
        require(condition==(w=='Z'*k),'common word rigidity')
        if condition:sync.append([w,list(v)])
        face=v[0]==v[2] and v[1]==v[3] and u[0]==u[2] and u[1]==u[3]
        if face:
            face_count+=1
            require((v[1]==0)==(w=='Z'*k),'face zero endpoint coefficient')
        words+=1
# Exact inverse gap identities and necessary parity divisibility.
for N in range(1,10001):
    dh,de=N,6*N
    require((20*dh-3*de,6*dh-de)==(2*N,0),'Fminus inverse')
    require((dh%2==de%2==0)==(N%2==0),'Fminus parity')
    if N%7==0:
        dh,de=N//7,6*N//7
        require((20*dh-de,6*dh-de)==(2*N,0),'Fplus inverse')
        require((dh%2==de%2==0)==(N%14==0),'Fplus parity')
for D in range(1,20,2):
    for E in range(1,20):
        for U in range(D,D+20):
            h=9*D-3*E+11*U;e=3*D-E+3*U;g=2*E
            lower=(20*h-e,6*h-e)
            upper=(20*(h+g)-3*e,6*(h+g)-e)
            delta=(upper[0]-lower[0],upper[1]-lower[1])
            require(delta==(6*(7*E-D-U),12*E),'branch cancellation')
            require(all(v%6==0 for v in delta),'content six')
print(json.dumps({'status':'PASS','all_words_max_length':10,'words_checked':words,'universal_face_words_checked':face_count,'synchronized_derivatives':sync,'inverse_gap_N_checked':10000,'branch_cancellation_cases_checked':3800,'scope':'Exact finite support for infinite proofs in ROUND2.md.'},indent=2))
