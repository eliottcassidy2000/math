#!/usr/bin/env python3
"""Independent full-permutation audit of the order-four-to-six native probe."""
from pathlib import Path
from itertools import permutations
from collections import Counter
from fractions import Fraction
import json, hashlib, argparse

BASE=Path(__file__).resolve().parent
SOURCE=None

def resolve_source(producer_dir=None):
    candidates=([Path(producer_dir).expanduser().resolve()] if producer_dir else
                [BASE.parent/'carrier-transfer', BASE.parent/'session47-carrier-transfer'])
    for directory in candidates:
        path=directory/'native-normal.json'
        if path.is_file():return path
    raise FileNotFoundError('Cannot find native-normal.json; pass --producer-dir with the carrier-transfer directory.')
CHECKS=0

def need(t,msg):
    global CHECKS
    CHECKS+=1
    if not t:raise RuntimeError(msg)

def perm_maps(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    ix={p:k for k,p in enumerate(pairs)}
    return [[(ix[min(p[i],p[j]),max(p[i],p[j])],p[i]>p[j]) for i,j in pairs]
            for p in permutations(range(n))]

MAPS={n:perm_maps(n) for n in (4,6)}
CACHE={}
def canonical(code,n):
    if (code,n) in CACHE:return CACHE[code,n]
    lo=1<<(n*(n-1)//2);hi=-1
    for mapping in MAPS[n]:
        value=0
        for k,flip in mapping:value=(value<<1)|(((code>>k)&1)^flip)
        lo=min(lo,value);hi=max(hi,value)
    converse=((1<<(n*(n-1)//2))-1)^hi
    result=(lo,converse,min(lo,converse));CACHE[code,n]=result
    return result

def from_word(n,w):
    B=[[0]*n for _ in range(n)];chords=iter((w>>k)&1 for k in range((n-1)*(n-2)//2))
    for i in range(n):
        for j in range(i+1,n):
            B[i][j]=1 if j==i+1 else next(chords);B[j][i]=1-B[i][j]
    return B

def encode(B):
    value=0;k=0
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            value|=B[i][j]<<k;k+=1
    return value

def matmul(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]

def matvec(A,x,m):return tuple(sum(a*b for a,b in zip(row,x))%m for row in A)

def audit_graph(n,stored):
    count=1<<((n-1)*(n-2)//2);full=count-1
    chords=[(i,j) for i in range(n) for j in range(i+2,n)];ix={p:k for k,p in enumerate(chords)}
    ordinary=set();merged_to_label={};label_to_merged={};ordinary_to_converse={}
    need(len(stored['owner'])==count,('complete word mass',n))
    for w in range(count):
        B=from_word(n,w);o,c,m=canonical(encode(B),n);label=stored['owner'][w]
        ordinary.add(o);ordinary_to_converse[o]=c
        if m in merged_to_label:need(merged_to_label[m]==label,('merged class split',n,w))
        if label in label_to_merged:need(label_to_merged[label]==m,('merged classes identified',n,w))
        merged_to_label[m]=label;label_to_merged[label]=m
        reflected=sum(((w>>k)&1)<<ix[n-1-j,n-1-i] for k,(i,j) in enumerate(chords))
        need(stored['owner'][reflected]==label,('reflection/converse quotient',n,w))
    need(len(label_to_merged)==len(stored['keys']),('full node support',n))
    need(set(label_to_merged)==set(range(len(stored['keys']))),('dense node labels',n))
    need(stored['keys']==sorted(set(stored['keys'])),('sorted unique representative keys',n))
    bitcount=n*(n-1)//2
    for label,key in enumerate(stored['keys']):
        decoded=sum(((key>>(bitcount-1-k))&1)<<k for k in range(bitcount))
        need(canonical(decoded,n)[2]==label_to_merged[label],('stored key has actual canonical class',n,label))
    colors=Counter();L=[[0]*len(label_to_merged) for _ in label_to_merged]
    # Select closed ports, not the producer's numeric min-complement representative.
    for w in range(count):
        if w>>ix[0,n-1]&1:continue
        v=w^full;u0,v0=stored['owner'][w],stored['owner'][v];u,v=sorted((u0,v0))
        reflected=sum(((w>>k)&1)<<ix[n-1-j,n-1-i] for k,(i,j) in enumerate(chords))
        blue=w==reflected;colors[u,v,'blue' if blue else 'black']+=1
        if u!=v:L[u][u]+=1;L[v][v]+=1;L[u][v]-=1;L[v][u]-=1
    pairs=sorted({(u,v) for u,v,color in colors})
    need([[u,v,colors[u,v,'blue']+colors[u,v,'black']] for u,v in pairs]==stored['edges'],('all edge multiplicities',n))
    need(L==stored['L'],('full literal Laplacian',n))
    need(all(colors[u,v,'black']%2==0 for u,v in pairs),('black reflection pairs',n))
    need(sum(colors.values())==count//2,('one line per complement pair',n))
    reached={0}
    while True:
        new=reached|{j for i in reached for j,w in enumerate(L[i]) if i!=j and w<0}
        if new==reached:break
        reached=new
    need(len(reached)==len(L),('connected full carrier',n))
    result={'order':n,'forward_words':count,'ordinary_classes':len(ordinary),'self_converse_classes':sum(o==c for o,c in ordinary_to_converse.items()),'merged_nodes':len(L),'literal_lines':sum(colors.values()),'blue_lines':sum(v for (a,b,c),v in colors.items() if c=='blue'),'black_lines':sum(v for (a,b,c),v in colors.items() if c=='black'),'blue_loops':sum(v for (a,b,c),v in colors.items() if c=='blue' and a==b),'black_loops':sum(v for (a,b,c),v in colors.items() if c=='black' and a==b),'edge_support':len(pairs),'connected':True,'colored_edges':[[u,v,colors[u,v,'blue'],colors[u,v,'black']] for u,v in pairs]}
    return result,merged_to_label

def main():
    raw=SOURCE.read_bytes();data=json.loads(raw)
    old,oldmap=audit_graph(4,data['old']);new,newmap=audit_graph(6,data['new'])
    need((old['ordinary_classes'],old['merged_nodes'])==(4,3),'complete old classes')
    need((new['ordinary_classes'],new['merged_nodes'])==(56,34),'complete new classes')
    readouts=set();native_records=[]
    for w in range(8):
        B=from_word(4,w);S=[[0]*6 for _ in range(6)]
        for i in range(4):
            for j in range(4):S[i][j]=B[i][j]
            S[4][i]=1;S[i][5]=1
        S[4][5]=1
        path=(4,0,1,2,3,5);kept={frozenset((a,b)) for a,b in zip(path,path[1:])}
        need(all(S[a][b]==1 for a,b in zip(path,path[1:])),'actual suspended Hamiltonian path')
        T=[[0 if i==j else S[i][j] if frozenset((i,j)) in kept else 1-S[i][j] for j in range(6)] for i in range(6)]
        need([r[:4] for r in T[:4]]==from_word(4,w^7),'literal inner complement')
        ns=newmap[canonical(encode(S),6)[2]];nt=newmap[canonical(encode(T),6)[2]]
        os=data['old']['owner'][w];ot=data['old']['owner'][w^7]
        readouts.add((ns,os));readouts.add((nt,ot))
        native_records.append({'word':w,'old_source':os,'old_target':ot,'new_source':ns,'new_target':nt,'path':path})
    need(sorted(readouts)==[tuple(x) for x in data['readouts_new_old']],'all literal endpoint readouts')
    need(all(len({o for x,o in readouts if x==node})==1 for node,_ in readouts),'no readout collisions')
    need((data['old_root'],data['new_root'])==(0,0),'verified compatible roots')
    coords=data['new_reduced_coordinates'];need(coords==list(range(1,34)),'all new nonroot coordinates')
    A=[r[1:] for r in data['old']['L'][1:]];N=[[data['new']['L'][i][j] for j in coords] for i in coords]
    need(A==[[2,-2],[-2,3]],'old grounded matrix')
    # det(A)=2 and gcd(entries)=1 independently give Smith(1,2).
    need(A[0][0]*A[1][1]-A[0][1]*A[1][0]==2,'old determinant')
    Rnew=[[int(j==node) for j in coords] for node,owner in sorted(readouts)]
    Rold=[[int(j==owner) for j in (1,2)] for node,owner in sorted(readouts)]
    M=N+Rnew;J=[[0,0] for _ in N]+Rold
    Z=[[int(node in (3,21)),0] for node in coords]
    MZ=matmul(M,Z);defect=[[v-w for v,w in zip(a,b)] for a,b in zip(MZ,J)]
    inv=[[Fraction(3,2),1],[1,1]];Q=matmul(defect,inv)
    need(all(q.denominator==1 if isinstance(q,Fraction) else True for row in Q for q in row),'integral certificate, not odd-denominator only')
    Q=[[int(x) for x in row] for row in Q]
    need(matmul(Q,A)==defect,'exact integer readout intertwiner MZ-J=QA')
    need(data['integer_certificate']['Z']==Z,'producer common integer lift matches independent construction')
    need(data['integer_certificate']['Q']==Q,'producer full integer relation matches independent construction')
    need(data['integer_certificate']['lift_support_new_nodes']==[3,21],'authenticated native support labels')
    for m in range(2,13):
        xs=((a,b) for a in range(m) for b in range(m))
        for x in xs:
            if not any(matvec(A,x,m)):
                need(matvec(M,matvec(Z,x,m),m)==matvec(J,x,m),'independent fixed lift at every small modulus')
    for test in data['tests']:
        m=test['modulus'];expected={(0,0),(m//2,0)}
        need({tuple(lift['old']) for lift in test['lifts']}==expected,'complete old nontrivial bank and lifts')
        need(test['old_harmonics']==test['extended']==2,'complete finite extension count')
        for lift in test['lifts']:
            x=lift['old'];y=lift['new']
            need(not any(matvec(N,y,m)),'literal full new grounded harmonic')
            need(matvec(Rnew,y,m)==matvec(Rold,x,m),'literal full endpoint readout')
            need(matvec(Z,x,m)==tuple(y),'common integer lift agrees with stored witnesses')
    return {'status':'PASS: independently authenticated complete native order4-to6 readout transfer; no order8-to10 claim','source_native_receipt_sha256':hashlib.sha256(raw).hexdigest(),'old_graph':old,'new_graph':new,'native_records':native_records,'readouts':sorted(readouts),'A':A,'Z':Z,'Q':Q,'identity':'M Z - J = Q A exactly over integers','old_smith':[1,2],'nonzero_old_bank':'(2^(r-1),0) modulo2^r','all_dyadic_extension':True,'all_integer_moduli_extension':True,'checks':CHECKS}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--producer-dir',type=Path,help='Directory containing native-normal.json; defaults to either supported sibling carrier-transfer directory.')
    args=parser.parse_args();SOURCE=resolve_source(args.producer_dir)
    print(json.dumps(main(),sort_keys=True,indent=2))
