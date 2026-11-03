"""Exact projective decoder, full bounded descendant fibres, and literal audit."""
import argparse,hashlib,json,math,subprocess,sys,tempfile
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent
sys.path.insert(0,str(RESULTS/'golden_inverse_defects_20261103'))
import inverse as old
sys.path.insert(0,str(RESULTS/'adaptive_state_induction_20261103'))
import live_port as live
TAU=(2,-1,0)
WORDS={'B':('b',),'Q':('q2',),'S':('s78',),'L2':('r24','a','a','a'),
       'L6':('r24','q4'),'L45':('r19','r32','r48','q4'),
       'L13':('r19','r19','a','a','r18'),'Lp2':('r24','a','a','r21')}
SEEDS=((31,14,25),(95,41,75),(97,40,75))
def need(b,s):
    if not b:raise RuntimeError(s)
def beta(v):h,x,y=v;return y,h+2*x-y
def step2(v,Q):return tuple(sum(v[i]*Q[i][j] for i in range(2)) for j in range(2))
def det(Q):return Q[0][0]*Q[1][1]-Q[0][1]*Q[1][0]
def slope(Q,r):return (Q[0][1]+r*Q[1][1])/(Q[0][0]+r*Q[1][0])
def inverse2(v,Q):
    d=det(Q);return ((v[0]*Q[1][1]-v[1]*Q[1][0])/d,(-v[0]*Q[0][1]+v[1]*Q[0][0])/d)
def digest(p):return hashlib.sha256(Path(p).read_text(encoding='utf-8').encode()).hexdigest()
def run(cmd):return subprocess.run(list(map(str,cmd)),check=True,capture_output=True,text=True).stdout.strip()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    bank={c:tuple(map(tuple,a['matrix'])) for c,a in old.ATLAS.items()}
    P={};Q={};gains={};costs={};intervals={};derivatives={}
    for name,word in WORDS.items():
        A=((1,0,0),(0,1,0),(0,0,1))
        for c in word:A=tuple(old.step(row,bank[c]) for row in A)
        P[name]=A;t=old.step(TAU,A);need(t[2]==0 and t[0]==-2*t[1],('error line',name))
        gains[name]=-t[1];costs[name]=sum(len(old.ATLAS[c]['source_rows']) for c in word)
        row1=beta(A[0]);row0=tuple(a+b for a,b in zip(row1,beta(A[2])));Q[name]=(row0,row1)
        for v in ((1,0,0),(0,1,0),(0,0,1)):need(beta(old.step(v,A))==step2(beta(v),Q[name]),('bulk intertwining',name))
        need(min(x for row in Q[name] for x in row)>0 and det(Q[name]),'positive invertible quotient')
        intervals[name]=tuple(sorted((slope(Q[name],F(1)),slope(Q[name],F(2)))))
        derivatives[name]=F(abs(det(Q[name])),(Q[name][0][0]+Q[name][1][0])**2)
    ordered=sorted(intervals,key=lambda a:intervals[a][0])
    need(all(F(1)<intervals[a][0]<intervals[a][1]<F(2) for a in ordered),'invariant interval')
    need(all(intervals[a][1]<intervals[b][0] for a,b in zip(ordered,ordered[1:])),'disjoint full images')
    need(F(34,25)<intervals[ordered[0]][0] and F(34,25)>1,'seed outside images')
    need(max(derivatives.values())==F(1,45),'uniform contraction bound')
    def decode_bulk(value):
        v=tuple(map(F,value));back=[]
        while v[1]/v[0]!=F(34,25):
            r=v[1]/v[0];hits=[a for a,I in intervals.items() if I[0]<=r<=I[1]]
            if len(hits)!=1:return None
            a=hits[0];u=inverse2(v,Q[a])
            if u[0]<=0 or not 1<=u[1]/u[0]<=2:return None
            need(u[0]<v[0],'decoding height decreases')
            v=u;back.append(a)
            if len(back)>1000:raise RuntimeError('test decoder cap')
        return tuple(reversed(back)),v[0]/25
    # Full membership requires an INTEGER predecessor and the actual seed charge.
    def membership(v):
        b=tuple(map(F,beta(v)));back=[]
        if b[0]<=0 or not 1<=b[1]/b[0]<=2:return None
        while b[1]/b[0]!=F(34,25):
            hits=[a for a,I in intervals.items() if I[0]<=b[1]/b[0]<=I[1]]
            if len(hits)!=1:return None
            a=hits[0];b=inverse2(b,Q[a])
            if b[0]<=0 or not 1<=b[1]/b[0]<=2 or any(z.denominator!=1 for z in b):return None
            back.append(a)
        word=tuple(reversed(back))
        for j,seed in enumerate(SEEDS):
            w=seed
            for a in word:w=old.step(w,P[a])
            if w==v:return j,word
        return None
    ratios={};states=set();checked=0
    for length in range(6):
        for word in product(WORDS,repeat=length):
            for j,seed in enumerate(SEEDS):
                v=seed
                for a in word:v=old.step(v,P[a])
                b=beta(v);r=F(b[1],b[0]);need(ratios.setdefault(r,word)==word,'projective injectivity')
                need(v not in states,'full source injectivity');states.add(v)
                need(decode_bulk(b)==(word,1 if j==0 else 3),'scale and word decode')
                need(membership(v)==(j,word),'literal bank membership decode');checked+=1
    fakes=(tuple(3*x+t for x,t in zip(SEEDS[1],TAU)),
                                         tuple(3*x+2*t for x,t in zip(SEEDS[2],TAU)))
    for v in fakes:need(membership(v) is None,'forbidden second-level source in the macro bank')
    # All five Q4 responses retain their different correction directions.
    phases={c:old.step(TAU,bank[c]) for c in ('q0','q1','q2','q4','q5')}
    need(sum(t==TAU for t in phases.values())==1,'five Q4 channels stay distinct')
    O=[[1,5,73],[0,8,144],[0,6,120]]
    minors=[O[a][c]*O[b][d]-O[a][d]*O[b][c]
            for a in range(3) for b in range(a+1,3)
            for c in range(3) for d in range(c+1,3)]
    need(math.gcd(*[x for row in O for x in row])==1 and math.gcd(*minors)==2
         and 8*120-6*144==96,'observation Smith factors 1,2,48')
    for v in SEEDS:
        hs=[sum(v[i]*O[i][j] for i in range(3)) for j in range(3)]
        h,h1,h2=hs
        need((h,F(20*h1-27*h-h2,16),F(h2-18*h1+17*h,12))==v,'exact 96-index observation inverse')
    # F49 transport is a linearized polynomial A*z+B*z^7, including singular maps.
    def add(x,y):return tuple((a+b)%7 for a,b in zip(x,y))
    def neg(x):return tuple(-a%7 for a in x)
    def mul(x,y):a,b=x;c,d=y;return ((a*c+2*b*d)%7,(a*d+b*c+4*b*d)%7)
    def power(a,k):
        r=(1,0)
        while k:
            if k&1:r=mul(r,a)
            a=mul(a,a);k//=2
        return r
    field={};joint={}
    need(power((0,1),48)==(1,0) and all(power((0,1),k)!=(1,0) for k in range(1,48)),'theta order 48')
    need(power((0,1),38)==(3,4),'literal seed field phase')
    for k in range(48):need(power(power((0,1),k+38),7)==power((0,1),7*k+36+38),'phased Frobenius')
    for name,A in Q.items():
        (a,b),(c,d)=A;f1=(8*d%7,2*c%7);ft=(4*b%7,a%7)
        fb=mul(add(ft,neg(mul(f1,(0,1)))),power((4,5),47));fa=add(f1,neg(fb))
        for y,x in product(range(7),repeat=2):
            alpha=(4*x%7,y);Y,X=step2((y,x),A)
            need(add(mul(fa,alpha),mul(fb,power(alpha,7)))==(4*X%7,Y%7),'full field transport')
        field[name]={'linear':fa,'frobenius':fb}
        split=[]
        for v in ((1,0,1),(0,4,0),(1,3,0)):
            h,x,y=old.step(v,P[name]);split.append(tuple(z%7 for z in (y,h+2*x-y,h-y)))
        need(split[2]==(0,0,gains[name]%7),'joint charge line')
        joint[name]=(split[0][0],split[0][1],split[1][0],split[1][1],split[0][2],split[1][2],split[2][2])
    # Independent Python forward/reverse closure of the unit subgroup with parity.
    def mm(a,b):return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))%7 for i in range(2) for j in range(2))
    G={a:tuple(x%7 for row in Q[a] for x in row) for a in ('B','Q')}
    unit={(1,0,0,1,1):()};todo=list(unit)
    for t in todo:
        for a in G:
            u=mm(t[:4],G[a])+(t[4]*gains[a]%7,)
            if u not in unit:unit[u]=unit[t]+(a,);todo.append(u)
    expected={q+(e,) for q in product(range(7),repeat=4) if (q[0]*q[3]-q[1]*q[2])%7 for e in(1,6)}
    need(set(unit)==expected,'complete GL2(7) with independent error sign')
    for t in unit:
        for a in G:need(mm(t[:4],G[a])+(t[4]*gains[a]%7,) in unit,'independent unit closure')
    need(joint['L2'][4:6]!=(0,0) and gains['L2']%7==5,'nontrivial shear and charge generator')
    frobenius_word=unit[(6,1,0,1,1)]
    frobenius_state=SEEDS[0]
    for a in frobenius_word:frobenius_state=old.step(frobenius_state,P[a])
    zero_state=SEEDS[0]
    for a in ('S','B','S'):zero_state=old.step(zero_state,P[a])
    need(all(x%7==0 for x in beta(zero_state)) and (zero_state[0]-zero_state[2])%7==3,'zero field retains source charge')
    # Exact presentation/marked-state recurrence and its independent word DP.
    N=90;coeff=[0]*(N+1)
    for n in range(N+1):coeff[n]=int(n in(7,8,9))+sum(coeff[n-c] for c in costs.values() if n>=c)
    word_dp=Counter({7:1,8:1,9:1})
    for n in range(N+1):
        for c in costs.values():
            if n+c<=N:word_dp[n+c]+=word_dp[n]
    need(coeff==[word_dp[n] for n in range(N+1)],'rooted source recurrence')
    # Incoming two-live-port language: whole bounded universe, not selected witnesses.
    signatures,matrices,seed,obs=live.compile_carrier();need(len(signatures)==24,'inherited carrier')
    state_list=[seed];index={seed:0};parents=[(-1,-1)];edges=[];H=[]
    for id,v in enumerate(state_list):
        h=live.output(v,obs);H.append(h);row=[]
        for a,M in enumerate(matrices):
            u=live.advance(v,M);hu=live.output(u,obs);need(hu>h,'strict actual growth')
            if hu>295:continue
            if u not in index:index[u]=len(state_list);state_list.append(u);parents.append((id,a))
            row.append((index[u],a))
        edges.append(row)
    def word(id):
        out=[]
        while id: id,a=parents[id];out.append(a)
        return tuple(reversed(out))
    order=sorted(range(len(H)),key=lambda i:H[i],reverse=True);fibres={}
    for h in (31,95,97):
        targets=[3*h+2*j for j in range(3)];tags=[0]*len(H)
        for i in order:
            for j,target in enumerate(targets):
                if H[i]==target:tags[i]|=1<<j
            for k,a in edges[i]:tags[i]|=tags[k]
        ids=[i for i,x in enumerate(H) if x==h];hist=Counter(tags[i] for i in ids)
        need(7 not in hist,'no source has all three descendant targets')
        examples={}
        for tag in sorted(hist):
            i=next(i for i in ids if tags[i]==tag);source=i;paths={}
            for j,target in enumerate(targets):
                if not(tag>>j&1):continue
                here=i;path=[]
                while H[here]!=target:
                    child,a=next((k,a) for k,a in edges[here] if tags[k]>>j&1)
                    path.append(a);here=child
                paths[str(target)]=path
            T,ports,roots=live.decode(word(source))
            examples[str(tag)]={'source_program':word(source),'source_rows':T,'ports':ports,'children':paths}
            for target,path in paths.items():need(live.count(live.decode(word(source)+tuple(path))[0])==int(target),'literal child replay')
        fibres[str(h)]={'fibre':len(ids),'targets':targets,'mask_histogram':dict(sorted(hist.items())),'examples':examples}
    need(fibres['31']['mask_histogram']=={0:6,2:3,4:10,5:18},'complete H31 support')
    need(fibres['95']['mask_histogram']=={0:41,1:37,2:11,3:10,4:8,5:2,6:1},'complete H95 support')
    need(fibres['97']['mask_histogram']=={0:122,1:84,2:8,4:8,5:2},'complete H97 support')
    cache=Path(tempfile.gettempdir())/'math-spectrum31-return-20261103';cache.mkdir(exist_ok=True)
    mpath=cache/'matrices.txt';spath=cache/'states.txt'
    with mpath.open('w',newline='\n') as f:
        f.write(' '.join(map(str,seed+obs))+'\n')
        for M in matrices:f.write(str(len(M))+'\n'+'\n'.join(' '.join(map(str,t)) for t in M)+'\n')
    with spath.open('w',newline='\n') as f:
        for v,p in zip(state_list,parents):f.write(' '.join(map(str,p+v))+'\n')
    run(['g++','-O3','-std=c++17',HERE/'audit.cpp','-o',cache/'audit.exe'])
    audit=run([cache/'audit.exe',mpath,spath])
    jpath=cache/'joint.txt';jpath.write_text('\n'.join(' '.join(map(str,joint[a])) for a in ('B','Q','S','L2'))+'\n')
    run(['g++','-O3','-std=c++17',HERE/'monoid.cpp','-o',cache/'monoid.exe'])
    monoid=run([cache/'monoid.exe',jpath])
    sources=['ternary_repair_clock_20261103/CARRIES.md','golden_inverse_defects_20261103/inverse.py',
             'golden_seam_recurrence_20261103/atom_atlas.json','return_channels_20261103/THEORY.md',
             'return_channels_20261103/INCOMING.md','adaptive_state_induction_20261103/live_port.py',
             'adaptive_state_induction_20261103/LIVE_PORTS.md','adaptive_state_induction_20261103/INCOMING.md',
             'session155_spectrum_20261103/returns/MIXED.md','session155_spectrum_20261103/prefix31/NOTE.md']
    result={'scope':'PROVED projective code / FINITE-EXACT complete descendant fibres',
            'words':WORDS,'matrices':P,'bulk':Q,'gains':gains,'costs':costs,
            'intervals':{a:list(map(str,I)) for a,I in intervals.items()},
            'derivative_bounds':{a:str(x) for a,x in derivatives.items()},
            'projective_addresses_checked':len(ratios),'source_states_checked':checked,
            'rooted_state_coefficients':coeff,'Q4_corrections':phases,'fibre_analysis':fibres,
            'field_transport':field,'joint_mod7':joint,'frobenius_word_with_positive_charge':frobenius_word,
            'frobenius_source_state':frobenius_state,'SBS_source_state':zero_state,
            'action_monoid_audit':monoid,
            'independent_audit':audit,'literal_state_set_sha256':digest(spath),
            'source_sha256':{s:digest(RESULTS/s) for s in sources}}
    receipt=HERE/'receipt.json';serial=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.write:receipt.write_text(serial,encoding='utf-8',newline='\n')
    else:need(json.loads(receipt.read_text())==json.loads(serial),'receipt reproduction')
    print(json.dumps({'status':'PASS','projective_addresses':len(ratios),'source_states':checked,
                      'bounded_states':len(state_list),'audit':audit}))
if __name__=='__main__':main()
