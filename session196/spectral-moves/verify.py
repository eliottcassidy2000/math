#!/usr/bin/env python3
"""Independent literal quotient words, actual graph DP, and bipartite matching."""
from itertools import permutations, combinations
from math import factorial
import json
from pathlib import Path


def check(c, msg):
    if not c: raise RuntimeError(msg)


def tournament(n, bits):
    rows=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            u,v=(i,j) if bits>>k&1 else (j,i)
            rows[u]|=1<<v; k+=1
    return rows


def validate(rows):
    n=len(rows)
    for i in range(n):
        check(not rows[i]>>i&1, 'loop')
        check(rows[i] < 1<<n, 'bad mask')
        for j in range(i):
            check(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1, 'not a tournament')


def hp(rows):
    n=len(rows)
    if n==0:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    full=(1<<n)-1
    for mask in range(1,1<<n):
        for i in range(n):
            z=dp[mask][i]
            if not z:continue
            choices=rows[i]&~mask&full
            while choices:
                bit=choices&-choices;choices-=bit
                j=bit.bit_length()-1
                dp[mask|bit][j]+=z
    return sum(dp[-1])


def inflate(rows, root, m):
    check(m>=1, 'positive multiplicity required')
    out=rows[:]
    for _ in range(m-1):
        n=len(out); out.append(0)
        # New clone is placed directly after root in the clone order.
        out[root]|=1<<n
        for j in range(n):
            if j==root:continue
            if out[root]>>j&1:out[n]|=1<<j
            else:out[j]|=1<<n
    validate(out)
    return out


def cross_graph(rows, r):
    O=[i for i in range(len(rows)) if rows[r]>>i&1]
    I=[i for i in range(len(rows)) if i!=r and rows[i]>>r&1]
    return O,I,[(x,y) for x in O for y in I if rows[x]>>y&1]


def matching(rows, r):
    O,I,edges=cross_graph(rows,r); neigh={x:[] for x in O}
    for x,y in edges:neigh[x].append(y)
    mates={}
    def augment(x, seen):
        for y in neigh[x]:
            if y in seen:continue
            seen.add(y)
            if y not in mates or augment(mates[y],seen):
                mates[y]=x;return True
        return False
    return sum(augment(x,set()) for x in O)


def quotient_words(rows, root):
    """a[k]: outside vertices once; k indistinguishable, separated root symbols.
    A word is specified uniquely by an outside permutation and a set of gaps.
    This literal counter does not use graph DP, Stirling numbers, or matching.
    """
    n=len(rows); outside=[i for i in range(n) if i!=root]
    a=[0]*(n+1)
    for p in permutations(outside):
        for mask in range(1,1<<n):
            word=[]
            for gap in range(n):
                if mask>>gap&1:word.append(root)
                if gap<n-1:word.append(p[gap])
            if all(rows[x]>>y&1 for x,y in zip(word,word[1:])):
                a[mask.bit_count()]+=1
    return a


def stirling(m,k):
    d=[0]*(k+1);d[0]=1
    for i in range(1,m+1):
        for j in range(k,0,-1):d[j]=j*d[j]+d[j-1]
        d[0]=0
    return d[k]


def evaluate(a,m):
    return sum(z*factorial(k)*stirling(m,k) for k,z in enumerate(a) if z)


def strong(rows):
    n=len(rows)
    for i in range(n):
        seen=1<<i
        while True:
            new=seen
            for j in range(n):
                if seen>>j&1:new|=rows[j]
            if new==seen:break
            seen=new
        if seen!=(1<<n)-1:return False
    return True


def insertion_witness(rows, r):
    """Independent exhaustive maximum matching, then literal path insertion."""
    O,I,edges=cross_graph(rows,r); best=[]
    def visit(k, used, selected):
        nonlocal best
        if len(selected)>len(best):best=selected[:]
        for j in range(k,len(edges)):
            x,y=edges[j]
            if x not in used and y not in used:visit(j+1,used|{x,y},selected+[(x,y)])
    visit(0,set(),[])
    nu=len(best); b=inflate(rows,r,nu+1)
    clones=[r]+list(range(len(rows),len(b)))
    # Their order is irrelevant here: no two are adjacent in the partial path.
    path=[clones[0]]
    for (x,y),c in zip(best,clones[1:]):path.extend([x,y,c])
    for v in range(len(b)):
        if v in path:continue
        for pos in range(len(path)+1):
            if (pos==0 or b[path[pos-1]]>>v&1) and (pos==len(path) or b[v]>>path[pos]&1):
                path.insert(pos,v);break
        else:raise RuntimeError('insertion failed')
    check(all(b[x]>>y&1 for x,y in zip(path,path[1:])), 'bad inserted path')
    check(sorted(path)==list(range(len(b))), 'not Hamiltonian')
    C=set(clones)
    check(sum(v in C and (i==0 or path[i-1] not in C) for i,v in enumerate(path))==nu+1, 'lost a clone run')
    return nu,path,b


def induced(rows, vs):
 return [sum(1<<j for j,v in enumerate(vs) if rows[u]>>v&1) for u in vs]

def profile(rows,r):
 h=e=R=L=0;n=len(rows)
 def legal(p):return all(rows[x]>>y&1 for x,y in zip(p,p[1:]))
 for p in permutations(range(n)):
  if legal(p):h+=1;e+=p[-1]==r
  j=p.index(r)
  R+=legal(p[:j]) and legal(p[j:])
  L+=legal(p[:j+1]) and legal(p[j+1:])
 d=hp(induced(rows,[i for i in range(n) if i!=r]))
 return (h,e,R,L,d)

def glue(x,rx,y,ry,m):
 x=induced(x,[rx]+[i for i in range(len(x)) if i!=rx]);y=induced(y,[ry]+[i for i in range(len(y)) if i!=ry])
 nx=len(x);ny=len(y);n=nx+ny-1+m;rows=[0]*n
 for i in range(nx):
  for j in range(nx):
   if x[i]>>j&1:rows[i]|=1<<j
 mp=[0]+list(range(nx,nx+ny-1))
 for i in range(ny):
  for j in range(ny):
   if y[i]>>j&1:rows[mp[i]]|=1<<mp[j]
 A=list(range(1,nx));B=mp[1:];Z=list(range(nx+ny-1,n))
 for i in A:
  for j in B+Z:rows[i]|=1<<j
 for i in Z:
  rows[i]|=1
  for j in B:rows[i]|=1<<j
 for i in Z:
  for j in Z:
   if i<j:rows[i]|=1<<j
 validate(rows);return rows


def main():
    graphs=roots=dp_checks=star_checks=insertion_checks=0
    amplifier_signs={}; witnesses={}
    # Complete labeled universe, no isomorphism filtering or graph-size inference.
    for n in range(1,6):
        for bits in range(1<<(n*(n-1)//2)):
            rows=tournament(n,bits);validate(rows);graphs+=1
            h=hp(rows)
            for r in range(n):
                roots+=1;a=quotient_words(rows,r);nu=matching(rows,r)
                top=max(k for k,z in enumerate(a) if z)
                check(a[1]==h, 'wrong a1')
                check(top==nu+1, 'wrong matching degree')
                O,I,edges=cross_graph(rows,r)
                has_common=bool(edges) and bool(set.intersection(*[set(e) for e in edges]))
                check((nu==1)==has_common, 'star classification wrong')
                star_checks+=1
                for m in range(1,5):
                    b=inflate(rows,r,m)
                    check(hp(b)==evaluate(a,m), 'Stirling substitution wrong')
                    dp_checks+=1
                    # Every descendant clone has exactly the old cross matching number.
                    for c in [r]+list(range(n,len(b))):
                        check(matching(b,c)==nu, 'clone matching number not retained')
                    if nu==1:
                        alpha=a[2];beta=h-2*alpha
                        check(hp(b)==alpha*(1<<m)+beta, 'affine law wrong')
                if nu==1:
                    beta=h-2*a[2]
                    if beta in (-1,1):
                        amplifier_signs.setdefault(n,{}).setdefault(h,set()).add(-beta)
                if n<=4 or (n==5 and bits%17==0):
                    nu2,path,b=insertion_witness(rows,r)
                    check(nu2==nu,'matching implementations differ');insertion_checks+=1
    # Disjoint-edge hostile: matching number 2, not a star.
    n=5;r=0;rows=[0]*n
    O=[1,2];I=[3,4]
    for x in O:rows[r]|=1<<x
    for y in I:rows[y]|=1<<r
    rows[1]|=1<<2;rows[3]|=1<<4
    for x in O:
        for y in I:
            if (x,y) in [(1,3),(2,4)]:rows[x]|=1<<y
            else:rows[y]|=1<<x
    a=quotient_words(rows,r);nu,path,expanded=insertion_witness(rows,r)
    check(nu==2 and a[3]>0,'disjoint-edge hostile not witnessed')
    witnesses['disjoint_edge_hostile']={'rows':rows,'root':r,'a':a,'matching':nu,'three_clone_rows':expanded,'Hamilton_path':path,'counts':[hp(inflate(rows,r,m)) for m in range(1,7)]}
    # Explicit physical branches and infinite same-sign return families.
    examples=[('H5_plus',[8,1,3,6],0),('H5_minus',[8,1,3,6],1),('H13_minus',[2,12,9,17,7],0),('negative_beta_control',[16,1,3,7,14],0)]
    for label,rows,r in examples:
        a=quotient_words(rows,r);check(matching(rows,r)==1,'bad affine example')
        alpha=a[2];beta=a[1]-2*alpha
        counts=[]
        for m in range(1,9):
            h=hp(inflate(rows,r,m));check(h==alpha*(1<<m)+beta,'family mismatch')
            counts.append(h)
        witnesses[label]={'rows':rows,'root':r,'alpha':alpha,'beta':beta,'counts_m1_to_8':counts,'m2_rows':inflate(rows,r,2)}
    # Check the fixed-template logarithmic bound locally by simultaneous cloning.
    simultaneous=0
    for rows in ([4,1,2],[8,1,3,6],[2,12,9,17,7]):
        n=len(rows);check(strong(rows),'template should be strong')
        for mask in range(1<<n):
            b=rows[:]
            for r in range(n):
                if mask>>r&1:b=inflate(b,r,2)
            h=hp(b)
            check(h>=1+(1<<max(1+(mask>>r&1) for r in range(n))),'triangle growth bound failed')
            simultaneous+=1
    synthesis_checks=0
    small_bank=[]
    for n in range(1,4):
        for bits in range(1<<(n*(n-1)//2)):
            x=tournament(n,bits)
            for r in range(n):small_bank.append((x,r,profile(x,r)))
    for x,rx,(h,e,R,L,d) in small_bank:
        sigma=R-h-d+e; U=R-d
        check(sigma>=0, 'negative rooted merge slack')
        for y,ry,(g,f,Q,V,t) in small_bank:
            alpha=U*V; beta=d*g-sigma*t
            for m in range(1,4):
                b=glue(x,rx,y,ry,m)
                check(hp(b)==alpha*(1<<m)+beta, 'two-piece synthesis mismatch')
                root=len(x)+len(y)-1
                check(matching(b,root)==(1 if U else 0), 'synthesis star mismatch')
                synthesis_checks+=1
    synthesis_examples=[
        ('H29_plus',[18,8,3,5,14],0,[0],0),
        ('H31_minus',[10,0,3,6],0,[4,1,2],0),
    ]
    for label,x,rx,y,ry in synthesis_examples:
        h,e,R,L,d=profile(x,rx);g,f,Q,V,t=profile(y,ry)
        alpha=(R-d)*V;beta=d*g-(R-h-d+e)*t
        vals=[]
        for m in range(1,7):
            b=glue(x,rx,y,ry,m);v=hp(b)
            check(v==alpha*(1<<m)+beta,'published synthesized family failed')
            vals.append(v)
        witnesses[label]={'X':x,'root_X':rx,'Y':y,'root_Y':ry,'alpha':alpha,'beta':beta,'counts_m1_to_6':vals,'m1_rows':glue(x,rx,y,ry,1)}
    out={'status':'PASS','labelled_graphs':graphs,'rooted_graphs':roots,'actual_graph_DP_checks':dp_checks,'star_checks':star_checks,'two_piece_synthesis_DP_checks':synthesis_checks,'constructive_matching_insertion_checks':insertion_checks,'simultaneous_clone_checks':simultaneous,'amplifier_signs':{str(n):{str(h):sorted(s) for h,s in sorted(d.items())} for n,d in amplifier_signs.items()},'witnesses':witnesses}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
