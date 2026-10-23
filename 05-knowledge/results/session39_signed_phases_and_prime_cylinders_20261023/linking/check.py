#!/usr/bin/env python3
"""Exact graph discriminant controls; no inherited Python imports.

Rebuilds all order-five fixed-path presentations and all 42 disjoint switches.
The n8 section consumes the frozen edge ledger and terminal potential witness.
"""
import argparse, gzip, hashlib, itertools as it, json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent
CHECKS=0
def need(p,m):
    global CHECKS
    CHECKS+=1
    if not p: raise RuntimeError(m)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def matvec(A,x): return [sum(a*b for a,b in zip(r,x)) for r in A]
def dot(x,y): return sum(a*b for a,b in zip(x,y))
def determinant(A):
    a=[r[:] for r in A]; prev=1; sign=1; n=len(a)
    for k in range(n-1):
        if not a[k][k]:
            j=next((j for j in range(k+1,n) if a[j][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        p=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=p*a[i][j]-a[i][k]*a[k][j]
                need(z%prev==0,'Bareiss exact division');a[i][j]=z//prev
            a[i][k]=0
        prev=p
    return sign*a[-1][-1]
def inverse(A):
    n=len(A);M=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        j=next(j for j in range(k,n) if M[j][k]);M[k],M[j]=M[j],M[k]
        p=M[k][k];M[k]=[x/p for x in M[k]]
        for i in range(n):
            if i!=k:
                p=M[i][k];M[i]=[a-p*b for a,b in zip(M[i],M[k])]
    R=[r[n:] for r in M]
    for i in range(n):
        for j in range(n):need(sum(A[i][k]*R[k][j] for k in range(n))==(i==j),'rational inverse')
    return R
def lap(n,blue,black):
    A=[[0]*n for _ in range(n)]
    for (u,v),w in list(blue.items())+[(e,2*w) for e,w in black.items()]:
        if u!=v:A[u][u]+=w;A[v][v]+=w;A[u][v]-=w;A[v][u]-=w
    return [r[:-1] for r in A[:-1]]
def connected(n,blue,black):
    seen={0};more=True
    while more:
        old=len(seen)
        for (u,v),w in list(blue.items())+list(black.items()):
            if w and (u in seen or v in seen):seen.update([u,v])
        more=len(seen)>old
    return len(seen)==n
def degrees(n,E):
    d=[0]*n
    for (u,v),w in E.items():d[u]+=w;d[v]+=w
    return d
def census(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    pairs=list(it.combinations(range(n),2));top=(1<<len(free))-1
    reflect=[free.index((n-1-j,n-1-i)) for i,j in free]
    def canonical(A):
        ds=[sum(r) for r in A];blocks=[[v for v in range(n) if ds[v]==d] for d in sorted(set(ds))]
        codes=[]
        for parts in it.product(*(it.permutations(b) for b in blocks)):
            p=sum(parts,());codes.append(sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs)))
        return tuple(sorted(ds)),min(codes)
    keys=[];colors=[]
    for x in range(top+1):
        A=[[0]*n for _ in range(n)]
        for i in range(n-1):A[i][i+1]=1
        for k,(i,j) in enumerate(free):
            u,v=(i,j) if x>>k&1 else (j,i);A[u][v]=1
        keys.append(min(canonical(A),canonical(list(map(list,zip(*A))))))
        colors.append(all((x>>i&1)==(x>>j&1) for i,j in enumerate(reflect)))
    nodes=sorted(set(keys));ix={k:i for i,k in enumerate(nodes)};blue=Counter();black=Counter()
    for x in range(top+1):
        if x<(top^x):(blue if colors[x] else black)[tuple(sorted((ix[keys[x]],ix[keys[top^x]])))]+=1
    need(all(w%2==0 for w in black.values()),'literal black twins')
    return nodes,blue,Counter({e:w//2 for e,w in black.items()})
def harmonic(A):
    return [x for x in it.product(range(2),repeat=len(A)) if all(z%2==0 for z in matvec(A,x))]
def forms(A,V):
    q1=[[dot(x,matvec(A,y))//2%2 for y in V] for x in V]
    q=[(F(dot(x,matvec(A,x)),8)-F(dot(x,[A[i][i] for i in range(len(A))]),4))%1 for x in V]
    return q1,q
def histogram(values):return {str(k):v for k,v in sorted(Counter(values).items())}
def edge_rows(E):return [[u,v,w] for (u,v),w in sorted(E.items()) if w]
def check_discriminant(A,V,q1,q):
    R=inverse(A);c=[A[i][i] for i in range(len(A))];X=[];N=len(A)
    for u in V:
        au=matvec(A,u);need(all(z%2==0 for z in au),'two-torsion integer representative')
        X.append([z//2 for z in au])
    for i,x in enumerate(X):
        need((F(dot(x,matvec(R,x))-dot(x,matvec(R,c)),2))%1==q[i],'inverse quadratic equals potential formula')
        for j,y in enumerate(X):
            need(dot(x,matvec(R,y))%1==F(q1[i][j],2),'inverse linking equals half form')
            k=V.index(tuple((a+b)%2 for a,b in zip(V[i],V[j])))
            need((q[k]-q[i]-q[j])%1==F(q1[i][j],2),'quadratic polarization')
        for j in range(N):
            xp=[x[k]+A[k][j] for k in range(N)]
            qp=(dot(xp,matvec(R,xp))-dot(xp,matvec(R,c)))/2
            need(qp%1==q[i],'quadratic representative independence')
    return R
def n5():
    nodes,L,B=census(5);n=len(nodes);need(n==10,'n5 merged nodes');A=lap(n,L,B);V=harmonic(A)
    Q,q=forms(A,V);need(len(V)==8,'three-dimensional harmonic space')
    need(sum(not any(row) for row in Q)==1,'native half form nonsingular')
    tau=determinant(A);need(tau==41096,'native n5 tree count');check_discriminant(A,V,Q,q)
    records=[];positive=[];negative=None;edges=sorted(B)
    for e,f in it.combinations(edges,2):
        a,b=e;c,d=f
        if len({a,b,c,d})<4:continue
        for c,d in [(c,d),(d,c)]:
            C=B.copy();C[e]-=1;C[f]-=1
            adds=[tuple(sorted((a,c))),tuple(sorted((b,d)))]
            for g in adds:C[g]+=1
            if not connected(n,L,C):continue
            need(min(C.values())>=0,'nonnegative switch');need(degrees(n,C)==degrees(n,B),'paired-black degree preservation')
            D=lap(n,L,C);need([D[i][i] for i in range(n-1)]==[A[i][i] for i in range(n-1)],'same characteristic covector')
            W=harmonic(D);need(W==V,'same mod-two harmonic space');Q2,q2=forms(D,V)
            for i,u in enumerate(V):
                uf=list(u)+[0];alpha=uf[a]-uf[d];beta=uf[b]-uf[c]
                need((q2[i]-q[i])%1==F(alpha*beta,2)%1,'quadratic switch law')
                for j,v in enumerate(V):
                    vf=list(v)+[0]
                    need((Q2[i][j]-Q[i][j])%2==(alpha*(vf[b]-vf[c])+beta*(vf[a]-vf[d]))%2,'half-form switch law')
            tree=determinant(D)
            rec={'remove':[list(e),list(f)],'add':[list(g) for g in adds],'tree_count':tree,'same_entire_dyadic_tower':Q2==Q,'same_quadratic_refinement':q2==q,'quadratic_value_counts':histogram(q2)}
            if Q2==Q:
                need(tree%8==0 and tree%16!=0,'same complete two-primary group independently from determinant')
                check_discriminant(D,V,Q2,q2)
                if q2!=q:
                    need(histogram(q2)!=histogram(q),'quadratic nonisometry');positive.append(rec)
                if tree==36680:negative=rec
            records.append(rec)
    need(len(records)==42,'all 42 connected switches');need(sum(z['same_entire_dyadic_tower'] for z in records)==17,'17 entire-tower invisible switches')
    need(len(positive)==4,'four quadratic-visible switches');need(negative is not None and negative['same_quadratic_refinement'],'original hostile also quadratic invisible')
    # Loops contribute neither to A nor to the cut discriminant.
    C=B.copy();C[0,0]+=1;need(lap(n,L,C)==A,'loop blindness');need(degrees(n,C)[0]==degrees(n,B)[0]+2,'loop degree distinction')
    return {'nodes':nodes,'blue_edges':edge_rows(L),'paired_black_edges':edge_rows(B),'root':9,'harmonic_vectors':V,'half_form_on_all_vectors':Q,'native_quadratic_values':[str(x) for x in q],'native_quadratic_value_counts':histogram(q),'native_tree_count':tau,'complete_two_primary':'(Z/2)^3','switch_count':42,'tower_invisible':17,'quadratic_visible_among_invisible':4,'quadratic_invisible_among_invisible':13,'original_hostile':negative,'quadratic_positive_controls':positive,'all_switches':records}
def n8(base):
    p=base/'edges.tsv.gz';need(sha(p)=='aecbcc993eac35fe4c6dca6741b62594131edaa34e157ea78c37859f451a3b4c','frozen n8 edge ledger')
    z=json.loads((HERE/'n8_terminal_lift.json').read_text());v=z['lift'];n=len(v);need(n==3528 and z['k']==8,'terminal potential dimensions');av=[0]*n;energy=0
    with gzip.open(p,'rt') as f:
        need(next(f).strip()=='source\ttarget\tblue\tblack','n8 edge schema')
        for line in f:
            a,b,l,k=map(int,line.split());need(k%2==0,'literal black twins');d=v[a]-v[b];w=l+k
            av[a]+=w*d;av[b]-=w*d;energy+=w*d*d
    need(av==[256*x for x in z['A_lift_divided_by_modulus']],'terminal current replay')
    need(v[-1]==0 and any(x%2 for x in v),'grounded potential has exact order256')
    need(energy==dot(v,av)==14960072960,'exact terminal energy')
    u=(energy//256)%256;need(u==153 and u%8==1,'native linking unit')
    need(13*13*u%256==1,'explicit square normalization');need(3*u%8==3,'odd-scaled linking unit separation')
    # A one-dimensional exact control proves that identical tower need not mean identical linking.
    for k in range(1,6):
        V1={x%2 for x in range(1<<k) if 4*x%(1<<k)==0};V2={x%2 for x in range(1<<k) if 12*x%(1<<k)==0}
        need(V1==V2,'4 and12 exact scalar tower spaces')
        for x in V1:
            a=next(y for y in range(1<<k) if y%2==x and 4*y%(1<<k)==0)
            need((4*a*a//(1<<k))%2==(12*a*a//(1<<k))%2,'4 and12 scalar tower forms')
    need({t*t%4 for t in [1,3]}=={1},'Z4 pairing units1 and3 nonisometric')
    return {'edge_sha256':sha(p),'terminal_lift_sha256':sha(HERE/'n8_terminal_lift.json'),'energy':energy,'energy_divided_by_256':energy//256,'generator_self_linking':'153/256','normalizing_multiplier':13,'normalized_cyclic_pairing':'1/256','odd_scale':3,'scaled_normalized_cyclic_pairing':'3/256','native_highest_unit_mod8':1,'scaled_highest_unit_mod8':3,'inherited_complete_two_primary':{'1':3358,'2':4,'8':1},'inherited_scope':'The full group decomposition is the frozen Session38 theorem; this checker independently replays its terminal integer witness, not the entire 3528-node tower.','counterfeit_scope':'Odd-scaled colored network, not a same-fibre tournament-port counterfeit.'}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--wiring',type=Path,default=HERE.parent/'session38-wiring');ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out={'status':'PASS','n5':n5(),'n8':n8(args.wiring),'checker_sha256':sha(Path(__file__)),'checks':CHECKS}
    out=json.loads(json.dumps(out))
    target=HERE/'receipt.json'
    if args.write:target.write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads(target.read_text())==out,'frozen receipt replay')
    print(json.dumps({'status':'PASS','checks':out['checks'],'n5_switches':42,'same_tower':17,'quadratic_separated':4,'n8_linking_unit_mod8':1,'odd_scaled_unit_mod8':3},sort_keys=True))
if __name__=='__main__':main()
