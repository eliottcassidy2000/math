#!/usr/bin/env python3
"""New unary dyadic and odd-digit lifting controls. No inherited census audit.
Run normally and with -O; checks deliberately do not rely on assert.
"""
from collections import deque
from itertools import product
from pathlib import Path
from hashlib import sha256
import argparse,json,random
HERE=Path(__file__).resolve().parent
S0=(-1,3,0,1)
F=[[-1,0,0,0],[0,26,5,-1],[0,8,0,0],[2,5,-1,0]]

def check(x,msg='check failed'):
    if not x:raise RuntimeError(msg)
def mm(A,B,m=None):
    C=[[sum(a*b for a,b in zip(r,c)) for c in zip(*B)] for r in A]
    return C if m is None else [[v%m for v in r] for r in C]
def mv(A,x,m=None):
    y=tuple(sum(a*b for a,b in zip(r,x)) for r in A)
    return y if m is None else tuple(v%m for v in y)
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def power(A,e,m=None):
    B=eye(len(A))
    while e:
        if e&1:B=mm(B,A,m)
        A=mm(A,A,m);e//=2
    return B
G=mm(F,F)
def step(s,a,mod=None):
    r,M,q,E=s
    t=(-r,(a*a+1)*M+a*q-E,(2*a-2)*M,2*r+a*M-q)
    return t if mod is None else (t[0],)+tuple(v%mod for v in t[1:])
def signed_affine(a):
    return [[-(a*a+1),-a,1,0],[-(2*a-2),0,0,0],[-a,1,0,-2],[0,0,0,1]]
def slope_at(k,e,j):
    b=(k+5)//6;return mv(power(F,2*(b+j)+e,1<<k),S0,1<<k)[1]
def compile_dyadic(k,m,sign):
    check(k>=1 and m%2==1 and sign in [-1,1],'invalid target')
    e=int(sign==1);j=0
    for i in range(1,k):
        if slope_at(k,e,j)%(1<<(i+1)) != m%(1<<(i+1)):j+=1<<(i-1)
    run=2*((k+5)//6+j)+e
    check(slope_at(k,e,j)==m%(1<<k),'bit compiler failed')
    check((-1 if run%2==0 else 1)==sign,'exact sign mismatch')
    return run

def physical(run):
    g=[{1},set()];O={0,1}
    for _ in range(run):
        n=len(g);v=n+3;out=[set(z) for z in g]+[set()for _ in range(4)]
        for u in range(n):
            out[u].update(range(n,v))
            if u in O:out[v].add(u)
            else:out[u].add(v)
        for u in range(n,v):out[u].update(range(u+1,v+1))
        g=out;O={n,n+2}
    return g,O

def endpoint_tables(g):
    n=len(g);N=1<<n;incoming=[0]*n
    for u in range(n):
        for v in g[u]:incoming[v]|=1<<u
    D=[[0]*n for _ in range(N)]
    for v in range(n):D[1<<v][v]=1
    for mask in range(1,N):
        bits=mask
        while bits:
            bit=bits&-bits;bits-=bit;v=bit.bit_length()-1;rest=mask^bit
            if not rest:continue
            prev=rest&incoming[v]
            while prev:
                bit2=prev&-prev;prev-=bit2;D[mask][v]+=D[rest][bit2.bit_length()-1]
    return D

def recount(g,O):
    n=len(g);full=(1<<n)-1
    end=endpoint_tables(g);rev=[set()for _ in g]
    for u in range(n):
        for v in g[u]:rev[v].add(u)
    start=endpoint_tables(rev);H=[sum(z)for z in end];H[0]=1
    SS=[sum(z[v]for v in O)for z in start];EE=[sum(z[v]for v in O)for z in end]
    M=sum(SS[a]*H[full^a]for a in range(1,full+1))
    W=sum(EE[a]*SS[full^a]for a in range(1,full))
    r=H[full]-W
    return (r,M,H[full]+r,EE[full])

def full_bfs(n):
    root=(S0[0],)+tuple(x%n for x in S0[1:]);seen={root:0};q=deque([root]);edges=sha256()
    while q:
        s=q.popleft()
        for a in (3,5):
            t=step(s,a,n)
            edges.update(json.dumps([s,a,t],separators=(',',':')).encode())
            if t not in seen:seen[t]=seen[s]+1;q.append(t)
    check(len(seen)==2*n**3,('not full orbit',n,len(seen)))
    return {'modulus':n,'states':len(seen),'edges_checked':2*len(seen),'max_distance':max(seen.values()),'canonical_bfs_edges_sha256':edges.hexdigest()}

def v2(n):
    check(n!=0,'valuation zero input');return (n&-n).bit_length()-1

def main():
    out={'scope':'new constructor controls, not a tournament census'}
    # Polynomial identity independently checks the displayed characteristic factor.
    I=eye(4);J=[[G[i][j]-I[i][j]for j in range(4)]for i in range(4)]
    Q=mm(G,G);Q=[[Q[i][j]-745*G[i][j]+64*I[i][j]for j in range(4)]for i in range(4)]
    check(mm(mm(J,J),Q)==[[0]*4 for _ in range(4)],'characteristic annihilation')
    out['matrix_G']=G
    base=[]
    for e in (0,1):
        x=mv(power(G,4,8),mv(power(F,e,8),S0,8),8)
        d1=tuple((a-b)%4 for a,b in zip(mv(G,x,4),x))
        d2=tuple((a-b)%8 for a,b in zip(mv(power(G,2,8),x,8),x))
        check(d1==(0,2,0,2)and d2==(0,4,0,4),'dyadic increment')
        base.append({'sign':-1 if e==0 else 1,'projected_state_mod8':x,'one_pair_increment_mod4':d1,'two_pair_increment_mod8':d2})
    out['dyadic_base']=base
    universes=[]
    for k in range(1,15):
        n=1<<k;b=(k+5)//6
        for e in (0,1):
            x=mv(power(F,2*b+e,n),S0,n);slopes=[]
            for _ in range(n//2):slopes.append(x[1]);x=mv(G,x,n)
            check(set(slopes)==set(range(1,n,2)),'unary orbit coverage')
            universes.append({'k':k,'sign':-1 if e==0 else 1,'slopes':len(slopes),'burn_in_pairs':b,'checksum':sha256(json.dumps(slopes).encode()).hexdigest()})
    out['complete_dyadic_blocks']=universes
    rng=random.Random(265);bitchecks=0
    for k in range(1,19):
        for sign in (-1,1):
            for _ in range(6):
                m=rng.randrange(1<<(k-1))*2+1;run=compile_dyadic(k,m,sign)
                check(mv(power(F,run,1<<k),S0,1<<k)[1]==m)
                bitchecks+=1
    out['bit_compiler_targets']=bitchecks
    examples=[]
    for m,sign in product((7,21,63),(-1,1)):
        run=compile_dyadic(6,m,sign);s=mv(power(F,run),S0)
        examples.append({'target_mod64':m,'sign':sign,'F5_run':run,'source_order':2+4*run,'actual_state':s,'actual_closing_count_t1':s[0]+2*s[1]})
    out['seven_twentyone_sixtythree_examples']=examples
    literal=[]
    for run in range(4):
        g,O=physical(run);expected=mv(power(F,run),S0);actual=recount(g,O)
        check(actual==expected,'literal graph/source mismatch')
        literal.append({'F5_run':run,'order':len(g),'state':actual})
    out['literal_source_recounts']=literal
    collatz=[]
    for h in range(2,9):
        m=(pow(3,-1,1<<h)*(1+(1<<(h-1))))%(1<<h)
        run=compile_dyadic(h,m,-1);s=mv(power(F,run),S0);H=-1+2*s[1]
        check(v2(3*H+1)==h,'Collatz carry readout')
        collatz.append({'halving_depth':h,'F5_run':run,'order_after_closing':4+4*run,'count_decimal_digits':len(str(H)),'count_mod_2_to_h_plus_2':H%(1<<(h+2))})
    out['collatz_prescribed_depths']=collatz
    # Low-order literal positive and hostile controls for the five carry faces.
    faces=0
    for M,r,t in product(range(1,100,2),(-1,1),range(1,8)):
        H=r+(1<<t)*M
        if H<=0:continue
        expected=(1+v2(3*M-1)if t==1 else 1)if r==-1 else (1 if t==1 else 2+v2(3*M+1)if t==2 else 2)
        check(v2(3*H+1)==expected);faces+=1
    out['carry_face_arithmetic_controls']=faces
    out['hostiles']={'no_integer_inverse':'F5 adds four vertices; neither p-adic projection nor group inverse is a physical contraction','no_scalar_gap_claim':'slope residue7 does not realize scalar H=7','exceptional_primes':[5,761]}
    odd=[]
    for p in (3,7,11,13,43,79):
        A=signed_affine(5);d=p*p-1
        for k in range(1,6):
            n=p**k;T=power(A,d*p**(k-1),n);expected=eye(4)
            for i,v in enumerate((1,-8,-13)):expected[i][3]=(p**(k-1)*v*pow(10,-1,p))%n
            check(T==expected,('odd lift',p,k,T,expected))
            odd.append({'p':p,'k':k,'matrix_checked':True})
    out['prime_power_top_digit_translations']=odd
    mixed=[]
    for a,b in product(range(1,5),repeat=2):
        n=3**a*7**b;A=signed_affine(5)
        if a>=2:
            L=48*3**(a-2)*7**b;T=power(A,L,n)
            check([[z%7**b for z in row]for row in T]==eye(4),'3 lift failed protecting7')
            expected=eye(4)
            for i,v in enumerate((1,-8,-13)):expected[i][3]=(3**(a-1)*2*7**b*pow(10,-1,3)*v)%3**a
            check([[z%3**a for z in row]for row in T]==expected,'3 digit coefficient')
        if b>=2:
            L=48*7**(b-1)*3**a;T=power(A,L,n)
            check([[z%3**a for z in row]for row in T]==eye(4),'7 lift failed protecting3')
            expected=eye(4)
            for i,v in enumerate((1,-8,-13)):expected[i][3]=(7**(b-1)*3**a*pow(10,-1,7)*v)%7**b
            check([[z%7**b for z in row]for row in T]==expected,'7 digit coefficient')
        mixed.append({'a':a,'b':b,'modulus':n})
    out['mixed_independent_digit_controls']=mixed
    out['finite_full_orbits']=[full_bfs(n)for n in (7,21,63)]
    out['passed']=True
    (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'dyadic_blocks':len(universes),'compiler_targets':bitchecks,'finite_orbits':out['finite_full_orbits']},indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--compile',nargs=3,type=int,metavar=('K','ODD_M','SIGN'));ap.add_argument('--literal',action='store_true');args=ap.parse_args()
    if args.compile:
        k,m,sign=args.compile;run=compile_dyadic(k,m,sign)
        out={'k':k,'target':m%(1<<k),'exact_offset':sign,'F5_run':run,'source_order':2+4*run}
        if run<=1000:out['actual_state']=mv(power(F,run),S0)
        if args.literal:
            check(run<=1000,'literal output guard: run>1000');g,O=physical(run);out['out_neighbors']=[sorted(v)for v in g];out['selected']=sorted(O)
        print(json.dumps(out,indent=2))
    else:main()
