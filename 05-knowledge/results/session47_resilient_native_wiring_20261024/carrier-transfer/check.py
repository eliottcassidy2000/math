#!/usr/bin/env python3
"""Exact independent finite-ring controls. Standard library; assertions unused."""
import argparse, hashlib, itertools, json, math
from pathlib import Path

CHECKS = 0

def require(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)

def vecs(n,m): return itertools.product(range(m), repeat=n)
def mv(A,x,m): return tuple(sum(a*b for a,b in zip(r,x)) % m for r in A)
def tr(A): return [list(x) for x in zip(*A)]
def mm(A,B,m): return [list(mv(tr(B),r,m)) for r in A]
def kernel(A,m): return {x for x in vecs(len(A[0]),m) if not any(mv(A,x,m))}
def image(A,m): return {mv(A,x,m) for x in vecs(len(A[0]),m)}
def inverse_small(A,m):
    if len(A)==1:
        if math.gcd(A[0][0],m)>1: return None
        return [[pow(A[0][0],-1,m)]]
    if len(A)==2:
        a,b=A[0]; c,d=A[1]; det=a*d-b*c
        if math.gcd(det,m)>1: return None
        u=pow(det,-1,m)
        return [[u*d%m,-u*b%m],[-u*c%m,u*a%m]]
    raise ValueError('inverse control is deliberately restricted to dimension <=2')

def laplacian(n, edges):
    L=[[0]*n for _ in range(n)]
    for i,j,w in edges:
        L[i][i]+=w; L[j][j]+=w; L[i][j]-=w; L[j][i]-=w
    return L

def reduced(n,edges):
    L=laplacian(n,edges)
    return [r[1:] for r in L[1:]]

def connected(n,edges):
    seen={0}
    while True:
        nxt=seen|{j for i,j,w in edges if w and i in seen}|{i for i,j,w in edges if w and j in seen}
        if nxt==seen: return len(seen)==n
        seen=nxt

def audit(A,N,m):
    n=len(A); k=len(N)-n
    B=[r[:n] for r in N[:n]]; C=[r[n:] for r in N[:n]]; D=[r[n:] for r in N[n:]]
    E=[r[n:] for r in N]; F=[r[:n] for r in N]
    H=kernel(A,m); Hnew=kernel(N,m)
    projected={z[:n] for z in Hnew}
    direct=H<=projected
    imE=image(E,m); imA=image(A,m)
    actual_extendible={h for h in H if tuple(-a%m for a in mv(F,h,m)) in imE}
    require(actual_extendible==H&projected,'individual image criterion')
    gate=actual_extendible==H
    annE=kernel(tr(E),m)
    dual=all(mv(tr(F),z,m) in imA for z in annE)
    newrelations={z[:n] for z in image(N,m) if not any(z[n:])}
    contraction=newrelations<=imA
    require(direct==gate==dual==contraction,'equivalent transfer criteria')
    invD=inverse_small(D,m)
    if invD is not None:
        sub=mm(mm(C,invD,m),tr(C),m)
        S=[[(B[i][j]-sub[i][j])%m for j in range(n)] for i in range(n)]
        schur=all(not any(mv(S,h,m)) for h in H)
        require(schur==direct,'unit Schur criterion')
        require(all(len([z for z in Hnew if z[:n]==h])==1 for h in H&projected),'unit D unique extension')
    split=False
    if direct:
        for entries in vecs(n*k,m):
            T=[entries[i*n:(i+1)*n] for i in range(k)]
            P=[[int(i==j) for j in range(n)] for i in range(n)]+T
            defect=mm(N,P,m)
            if all(not any(mv(defect,h,m)) for h in H):
                require(all(tuple(r) in imA for r in defect),'ambient intertwiner row criterion')
                split=True; break
    require(not split or direct,'section implies coverage')
    return 'split' if split else ('nonsplit' if direct else 'loss'), invD is not None

def finite_census():
    rows=[]
    for family in (1,2,3):
        graphs=[]
        if family in (1,2):
            # Old root-v edge weight 1..4; all new edges 0..2, no old edge changed.
            k=family
            pairs=[(0,2),(1,2)] if k==1 else [(0,2),(1,2),(0,3),(1,3),(2,3)]
            for a in range(1,5):
                for ws in itertools.product(range(3),repeat=len(pairs)):
                    edges=[(0,1,a)]+[(u,v,w) for (u,v),w in zip(pairs,ws)]
                    if connected(2+k,edges): graphs.append(([[a]],reduced(2+k,edges)))
        else:
            # Eight positively weighted old triangles, one new vertex, all incident weights 0..2.
            for oldws in itertools.product((1,2),repeat=3):
                old=[(u,v,w) for (u,v),w in zip(((0,1),(0,2),(1,2)),oldws)]
                for ws in itertools.product(range(3),repeat=3):
                    if not any(ws): continue
                    edges=old+[(i,3,w) for i,w in enumerate(ws)]
                    graphs.append((reduced(3,old),reduced(4,edges)))
        status_by_graph=[{} for _ in graphs]
        for m in (2,3,4,5,8):
            count={'split':0,'nonsplit':0,'loss':0,'unit_D':0}
            for index,(A,N) in enumerate(graphs):
                status,unit=audit(A,N,m); count[status]+=1; count['unit_D']+=unit
                status_by_graph[index][m]=status
            rows.append({'family':family,'modulus':m,'graphs':len(graphs),**count})
        for (A,N),statuses in zip(graphs,status_by_graph):
            if len(A)==1: diagonal=[A[0][0]]
            else:
                first=math.gcd(*[a for r in A for a in r])
                determinant=A[0][0]*A[1][1]-A[0][1]*A[1][0]
                diagonal=[first,determinant//first]
            levels=sorted({(d&-d).bit_length()-1 for d in diagonal if d%2==0})
            if levels:
                gate=all(statuses[2**a]!='loss' for a in levels)
                require(gate==(statuses[2**max(levels)]=='split'),'distinct Smith levels iff top additive section')
            else: require(all(statuses[m]=='split' for m in (2,4,8)),'odd determinant trivial dyadic transfer')
    return rows

def hostile_controls():
    out={}
    A=[[2]]; N=[[3,-1],[-1,2]]
    require(kernel(A,8)=={(0,),(4,)},'old double edge has nonzero mod8 guard')
    require(kernel(N,8)=={(0,0)},'retained-endpoint triangle kills guard')
    out['endpoint_loss']={'old':A,'new':N,'modulus':8,'old_kernel':sorted(kernel(A,8)),'new_kernel':sorted(kernel(N,8))}
    N=[[4,-1,-1],[-1,2,0],[-1,0,2]]
    H4=kernel(N,4); lifts=sorted(z for z in H4 if z[0]==2)
    require(lifts==[(2,1,3),(2,3,1)],'nonsplit exact lifts')
    require(all(any(2*a%4 for a in z) for z in lifts),'all lifts have order4')
    require(audit(A,N,4)[0]=='nonsplit','pointwise but no ambient section')
    require(audit(A,N,2)[0]=='loss','mod4 success does not imply mod2 success')
    for m in (4,8,16,32): require({z[0] for z in kernel(N,m)}=={0,m//2},'all larger levels cover old kernel')
    out['nonsplit']={'old':A,'new':N,'modulus':4,'new_kernel':sorted(H4),'old_nonzero_lifts':lifts,'fails_modulus':2}
    A_top=[[2,0],[0,4]]
    N_top=[[4,0,-1,-1],[0,4,0,0],[-1,0,2,0],[-1,0,0,2]]
    require(audit(A_top,N_top,4)[0]=='nonsplit','true saturation top can pass without section')
    require(audit(A_top,N_top,2)[0]=='loss','true saturation top does not suffice for pointwise extension')
    out['true_top_failure']={'old':A_top,'new':N_top,'old_dyadic_exponent':4,'top_modulus_passes':4,'lower_modulus_fails':2}

    A=[[4]]; N=[[6,-2],[-2,4]]
    H=kernel(N,4); projected={z[0] for z in H}
    require(projected=={0,2},'singular Schur exact projection')
    require(20*1%4==0 and 1 not in projected,'cleared determinant false positive')
    require(5*2%4!=0 and 2 in projected,'integral rational Schur false negative')
    out['nonunit_schur']={'old':A,'new':N,'modulus':4,'projected':sorted(projected),'rational_schur':5,'cleared_schur':20}
    # Physical owner split: two old labels, but two actual new owners share v's record.
    old=reduced(2,[(0,1,4)])
    new=reduced(3,[(0,1,4),(0,2,4)])
    P=[[1],[1]]; x=(2,); lifted=mv(P,x,8)
    require(not any(mv(old,x,8)) and not any(mv(new,lifted,8)),'split owner original harmony')
    delta=[[-2,0],[0,2]]
    defect=mv(delta,lifted,8)
    require(defect==(4,4) and sum(defect)%8==0,'collapsed debt hides two physical failures')
    out['owner_split']={'old':old,'new':new,'pullback':P,'old_character':x,'modulus':8,'paired_edge_transfer_delta':delta,'physical_defect':defect}
    return out

def positive_controls():
    # Attach arbitrary unit triangle at old v, pull back its value constantly.
    oldL=laplacian(2,[(0,1,2)])
    newL=laplacian(4,[(0,1,2),(1,2,1),(2,3,1),(3,1,1)])
    P=[[1,0],[0,1],[0,1],[0,1]]; J=[[1,0],[0,1],[0,0],[0,0]]
    product=lambda A,B:[[sum(a*b for a,b in zip(r,c)) for c in zip(*B)] for r in A]
    require(product(newL,P)==product(J,oldL),'exact rooted attachment intertwiner')
    # Uniform two-fold blowup; arbitrary internal fibre edges do not affect pullback.
    oldedges=[(0,1,1),(1,2,2),(0,2,3)]; oldL=laplacian(3,oldedges)
    edges=[(2*u+i,2*v+j,w) for u,v,w in oldedges for i in (0,1) for j in (0,1)]
    edges += [(2*i,2*i+1,1+i) for i in range(3)]
    newL=laplacian(6,edges); P=[[int(i//2==j) for j in range(3)] for i in range(6)]
    require(product(newL,P)==[[2*x for x in r] for r in product(P,oldL)],'exact equitable blowup intertwiner')
    return {'rooted_attachment':'Lnew P = J Lold over integers','uniform_blowup':'Lnew P = 2 P Lold over integers'}

def saturation_controls():
    # Exact Smith-diagonal models with unimodular coordinate changes checked separately.
    rows=[]
    for diagonal in ((1,), (2,), (4,), (6,), (2,4), (2,12), (3,8)):
        s=max((d&-d).bit_length()-1 for d in diagonal)
        if not s:
            require(all(math.gcd(d,2)==1 for d in diagonal),'trivial old dyadic bank')
            continue
        A=[[d*int(i==j) for j in range(len(diagonal))] for i,d in enumerate(diagonal)]
        H=kernel(A,2**s)
        for k in range(s,s+4):
            scaled={tuple(2**(k-s)*a for a in x) for x in H}
            require(scaled==kernel(A,2**k),'dyadic saturation exact set equality')
        rows.append({'smith_diagonal':diagonal,'s':s,'base_kernel_size':len(H),'highest_tested_level':s+3})
    # Explicit non-diagonal integral unimodular conjugacy.
    A=[[2,2],[2,6]] # U diag(2,4) U^t, U=[[1,0],[1,1]].
    H=kernel(A,4)
    for k in range(2,6):
        require({tuple(2**(k-2)*a for a in x) for x in H}==kernel(A,2**k),'non-diagonal saturation')
    return rows

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--out',type=Path); args=parser.parse_args()
    report={'status':'FINITE-EXACT controls; no native enlarged Laplacian claimed','universe':'Connected positive old graph, retained old edge weights, nonnegative added weights; explicit labelled families, not isomorphism census','census':finite_census(),'hostiles':hostile_controls(),'positive':positive_controls(),'saturation':saturation_controls()}
    report['checks']=CHECKS
    payload=json.dumps(report,sort_keys=True,indent=2)+'\n'
    if args.out: args.out.write_text(payload)
    else: print(payload,end='')
if __name__=='__main__': main()
