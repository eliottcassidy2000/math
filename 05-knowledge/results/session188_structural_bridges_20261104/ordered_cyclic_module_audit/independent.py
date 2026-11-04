#!/usr/bin/env python3
"""Independent exact positive-kernel, marked-profile and literal C3 audit."""
from pathlib import Path
from itertools import permutations
from math import comb
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def add(*polys):
    d={}
    for p in polys:
        for e,v in p.items():d[e]=d.get(e,0)+v
    return {e:v for e,v in d.items() if v}

def scale(p,k):return {e:v*k for e,v in p.items() if v*k}
def ordinary_mul(*polys):
    z={tuple(0 for _ in range(10)):1}
    for p in polys:
        d={}
        for e,v in z.items():
            for f,w in p.items():
                g=tuple(a+b for a,b in zip(e,f));d[g]=d.get(g,0)+v*w
        z={e:v for e,v in d.items() if v}
    return z

def symbolic():
    one={(0,)*10:1};v=[]
    for i in range(10):e=[0]*10;e[i]=1;v.append({tuple(e):1})
    alpha,beta,a,b,sigma,mu,kappa,tau,x,y=v;mul=ordinary_mul
    h=add(x,y,mul(x,y));D=add(one,scale(mul(sigma,h),-1));U=add(mul(mu,h),mul(x,add(one,y)));V=add(mul(kappa,h),mul(y,add(one,x)))
    Nxx=mul(x,add(one,scale(mul(sigma,y),-1)));Nyy=mul(y,add(one,scale(mul(sigma,x),-1)));Nxy=mul(x,y,add(one,sigma));Nyx=mul(x,y,sigma)
    Cx=x;Cy=mul(y,add(one,x));Rx=mul(x,add(one,y));Ry=y
    AL=add(mul(alpha,D),mul(a,V));BE=add(mul(beta,D),mul(b,U));TA=add(mul(tau,D),Nxy,mul(mu,Cy),mul(kappa,Rx),mul(mu,kappa,h))
    Ex=mul(a,Cx);Ey=mul(a,Cy);Sx=mul(b,Rx);Sy=mul(b,Ry);rx=add(Nxx,mul(mu,Cx));ry=add(Nxy,mul(mu,Cy));lx=add(Nxy,mul(kappa,Rx));ly=add(Nyy,mul(kappa,Ry))
    lhs=add(scale(mul(AL,BE,add(Nxy,scale(Nyx,-1))),2),scale(mul(add(D,scale(TA,2)),add(mul(Ex,Sy),scale(mul(Ey,Sx),-1))),-1),mul(AL,add(mul(rx,Sy),scale(mul(ry,Sx),-1))),mul(BE,add(mul(Ex,ly),scale(mul(Ey,lx),-1))))
    base=add(scale(mul(alpha,beta),2),mul(alpha,b),mul(a,beta))
    mid=add(mul(alpha,b,U),mul(a,beta,V),mul(a,b,add(mul(add(kappa,mu,scale(tau,2)),h),scale(x,2),scale(y,2),mul(x,y))))
    rhs=mul(x,y,add(mul(base,D,D),mul(mid,D),scale(mul(a,b,U,V),2)))
    need(lhs==rhs,'full ten-variable integer numerator identity')
    return [[*e,v] for e,v in sorted(lhs.items())]

def graph(n,code):
    a=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
            k+=1
    return a

def op(a):return [((1<<len(a))-1)^(1<<i)^r for i,r in enumerate(a)]
def endpoints(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        if mask&(mask-1)==0:continue
        for i in range(n):
            if mask>>i&1:d[mask][i]=sum(d[mask^(1<<i)][j] for j in range(n) if mask>>j&1 and a[j]>>i&1)
    return d

def marked(a,mark):
    n=len(a);E=endpoints(a);S=endpoints(op(a));H=list(map(sum,E));q=[[0]*(1<<n) for _ in range(n+1)];q[0][0]=1
    for k in range(1,n+1):
        for mask in range(1,1<<n):
            part=mask
            while part:q[k][mask]+=H[part]*q[k-1][mask^part];part=(part-1)&mask
    profiles=[]
    for counts in (E,S):
        vals=[0]*(n+1)
        for k in range(1,n+1):
            for mask in range(1,1<<n):vals[k]+=counts[mask][mark]*q[k-1][((1<<n)-1)^mask]
        profiles.append(vals)
    need(profiles[0]==profiles[1] and profiles[0][1]>0,'physical marked balance and strictness')
    return profiles[0]

def interfaces(E,ia,im):
    e=len(E);one={(0,0,0):1};out=[dict(one),dict(one),dict(one),dict(one),{},{},{},{}]
    # alpha,beta,a,b,sigma,mu,kappa,tau. Enumerate every physical E path once.
    for length in range(1,e+1):
        for path in permutations(range(e),length):
            if not all(E[path[k]]>>path[k+1]&1 for k in range(length-1)):continue
            first,last=path[0],path[-1];key=(sum(1<<i for i in path),0,0)
            targets=[]
            if ia>>last&1:targets.append(0)
            if not ia>>first&1:targets.append(1)
            if im>>last&1:targets.append(2)
            if not im>>first&1:targets.append(3)
            if not im>>first&1 and im>>last&1:targets.append(4)
            if not ia>>first&1 and im>>last&1:targets.append(5)
            if not im>>first&1 and ia>>last&1:targets.append(6)
            if not ia>>first&1 and ia>>last&1:targets.append(7)
            for k in targets:out[k][key]=out[k].get(key,0)+1
    return out

def positive_value(E,ia,im,X,Y,mx,my):
    nx,ny=len(X),len(Y);e=len(E);bx,by=marked(X,mx),marked(Y,my);one={(0,0,0):1};x={(0,1,0):1};y={(0,0,1):1}
    def mul(*polys):
        out=one
        for p in polys:
            d={}
            for (mask,i,j),v in out.items():
                for (other,k,l),w in p.items():
                    if mask&other or i+k>nx or j+l>ny:continue
                    key=(mask|other,i+k,j+l);d[key]=d.get(key,0)+v*w
            out={key:v for key,v in d.items() if v}
        return out
    alpha,beta,a,b,sigma,mu,kappa,tau=interfaces(E,ia,im);h=add(x,y,mul(x,y));U=add(mul(mu,h),mul(x,add(one,y)));V=add(mul(kappa,h),mul(y,add(one,x)))
    Z=mul(sigma,h);power=one;D1={};D2={};D3={}
    for k in range(e+1):D1=add(D1,power);D2=add(D2,scale(power,k+1));D3=add(D3,scale(power,comb(k+2,2)));power=mul(power,Z)
    need(not power,'nilpotent exterior return termination')
    base=add(scale(mul(alpha,beta),2),mul(alpha,b),mul(a,beta));mid=add(mul(alpha,b,U),mul(a,beta,V),mul(a,b,add(mul(add(kappa,mu,scale(tau,2)),h),scale(x,2),scale(y,2),mul(x,y))))
    K=mul(x,y,add(mul(base,D1),mul(mid,D2),scale(mul(a,b,U,V,D3),2)))
    need(all(v>=0 for v in K.values()),'positive expanded full-support series')
    value=sum(v*bx[i]*by[j] for (mask,i,j),v in K.items() if mask==(1<<e)-1)
    alpha_beta=mul(alpha,beta).get(((1<<e)-1,0,0),0);lower=2*bx[1]*by[1]*alpha_beta
    need(value>=lower>0,'exact strict lower bound')
    return value,lower,bx,by

def build(E,ia,im,X,Y,transitive=False):
    nx,ny=len(X),len(Y);ox=3;oy=3+nx;oe=3+nx+ny;n=oe+len(E);a=[2,4,1]+[r<<ox for r in X]+[r<<oy for r in Y]+[r<<oe for r in E]
    def edge(u,v):a[u]|=1<<v
    for p in range(3):
        for i in range(ox,oy):edge(i,p) if transitive else edge(p,i)
        for j in range(oy,oe):edge(p,j) if transitive else edge(j,p)
    for i in range(ox,oy):
        for j in range(oy,oe):edge(i,j)
    for v in range(len(E)):
        for p in range(3):edge(oe+v,p) if ia>>v&1 else edge(p,oe+v)
        for i in range(ox,oe):edge(oe+v,i) if im>>v&1 else edge(i,oe+v)
    U=a[:]
    for i in range(3):U[i]^=7^(1<<i)
    return a,U,ox,oy

def main():
    algebra=symbolic();cases=[];records=[]
    def add_case(E,ia,im,X,Y,mx=0,my=0,tag='arbitrary_exterior'):
        value,lower,bx,by=positive_value(E,ia,im,X,Y,mx,my);T,U,ox,oy=build(E,ia,im,X,Y);masks=[1,1<<(ox+mx),2,1<<(oy+my)];queries=[]
        for bits in range(16):queries.append((sum(masks[k] for k in (0,1) if bits>>k&1),sum(masks[k] for k in (2,3) if bits>>k&1)))
        cases.append((T,U,queries));records.append(dict(tag=tag,E=E,IA=ia,IM=im,X=X,Y=Y,marked_X=mx,marked_Y=my,balanced_X=bx,balanced_Y=by,predicted_coefficient=value,strict_lower_bound=lower))
    for n in range(4):
        for code in range(1<<(n*(n-1)//2)):
            for ia in range(1<<n):
                for im in range(1<<n):add_case(graph(n,code),ia,im,[0],[0])
    c3=[2,4,1];w5=[14,28,8,16,5]
    for X,Y,mx,my in [(c3,c3,0,1),(w5,[0],4,0),([0],w5,0,4),(w5,c3,4,2)]:
        for n,E in [(0,[]),(1,[0]),(2,[2,0])]:
            for ia,im in [(0,0),((1<<n)-1,0),(0,(1<<n)-1),((1<<n)-1,(1<<n)-1)]:add_case(E,ia,im,X,Y,mx,my,'balanced_marked_children')
    # Excluded transitive attachment with singleton sides has zero strict quartic.
    T,U,ox,oy=build([],0,0,[0],[0],True);queries=[]
    for bits in range(16):queries.append(((1 if bits&1 else 0)|(1<<ox if bits&2 else 0),(2 if bits&4 else 0)|(1<<oy if bits&8 else 0)))
    cases.append((T,U,queries));records.append(dict(tag='excluded_transitive_attachment',predicted_coefficient=0))
    lines=[str(len(cases))]
    for T,U,qs in cases:lines.extend([f'{len(T)} {len(qs)}',' '.join(map(str,T)),' '.join(map(str,U))]+[f'{a} {c}' for a,c in qs])
    with tempfile.TemporaryDirectory() as td:
        exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);raw=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
    need(len(raw)==len(cases)*32,'complete literal output')
    for k,record in enumerate(records):
        vals=[raw[k*32+2*j:k*32+2*j+2] for j in range(16)];co=sum((1 if bits.bit_count()%2==0 else -1)*(f-g) for bits,(f,g) in enumerate(vals));need(co==record['predicted_coefficient'],'literal full quartic against positive functional');record['actual_16_query_pairs']=vals;record['actual_coefficient']=co
    return dict(status='PASS',independent_full_numerator=algebra,symbolic_numerator_monomials=len(algebra),full_exterior_cases=549,marked_child_cases=48,excluded_boundary_cases=1,literal_physical_values=len(raw),records=records,scope='Exact ordered-module C3 quartic, arbitrary exterior and balanced physical marked endpoints; no strong ambient lift.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
