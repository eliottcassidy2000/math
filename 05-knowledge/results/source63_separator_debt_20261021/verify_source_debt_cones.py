#!/usr/bin/env python3
"""Exact staged replay. No orbit census, floating decisions, or disabled asserts."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import base64, hashlib, json, sys, zlib
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
COUNTS={}
def need(v, label):
    if not v: raise RuntimeError(label)
    COUNTS[label]=COUNTS.get(label,0)+1

def frac(x): return str(x)
def affine(w):
    A,Q,B,C,E=1,1,0,0,0
    for a in w:
        need(isinstance(a,int) and a>=1,'positive valuation letter')
        A*=3;B=3*B+Q;C=3*C+(Q*2**a if a==1 else 0)
        E=3*E+(Q*(2**a-4) if a>=3 else 0);Q*=2**a
    return A,Q,B,C,E

def suffix(w):
    lam=F(1);P=D=F(0)
    for a in reversed(w):
        if a==1:P+=lam
        if a>=3:D+=lam*(1-F(4,2**a))
        lam*=F(3,2**a)
    return lam,P,D

def state(w):
    A,Q,B,C,E=affine(w);l,P,D=F(A,Q),F(C,Q),F(E,Q)
    need((l,P,D)==suffix(w),'integer recurrence versus suffix sums')
    need(F(B,Q)==1-l+P-D,'signed affine identity')
    return l,P,D

def cyl(w):
    A,Q,B,C,E=affine(w)
    return ((Q-B)*pow(A,-1,2*Q))%(2*Q),2*Q

def replay(n,w,root_stop=False):
    x=n;clock=0;first_odd=first_T=first_half_T=None;ordinary_clock=0
    for i,a in enumerate(w,1):
        need(x>1 or not root_stop,'ROOT respected at every proper input')
        z=3*x+1
        need((z & -z).bit_length()-1==a,'actual exact valuation')
        ordinary_clock+=1
        for j in range(a):
            z//=2;clock+=1;ordinary_clock+=1
            if first_T is None and z<n:first_T=(clock,ordinary_clock,z)
            if first_half_T is None and 2*z<=n:first_half_T=(clock,ordinary_clock,z)
        x=z
        if first_odd is None and x<n:first_odd=(i,clock,ordinary_clock,x)
    return x,first_odd,first_T,first_half_T

def actual_word(n,stop):
    x=n;w=[]
    for _ in range(10000):
        y=3*x+1;a=(y&-y).bit_length()-1;x=y>>a;w.append(a)
        if stop(x,w):return tuple(w),x
    raise RuntimeError('explicit discovery cap hit')

U=(1,1,1,1,1,4,1,2,1,1,2,1,1,1,2,3,1,1,2,1,2,1,1,1,1,1,3,1,1,1,4,2,2)
V=(1,1,1,3)
W=U+(4,)
ROOT_SUFFIX=(3,1,1,5,4)
l,P,D=state(U);theta=F(2,3);c=D/(l-theta)
need(D>l-1,'source63 outside ordinary unit cone')
need(c<1,'one separator repairs source63 cone')
need(c==F(14060731968178560,14425381885981321),'source63 sharp separator constant')
need(D/(l-1)<=F(20,13),'source63 inherited terminal-packet bound')
need(D/(l-1)>1,'zero separator cannot attain c<=1')
need(state(V)[2]<=c*(state(V)[0]-theta),'changing generator admitted')
for gen in (U,V):
    for i in range(1,len(gen)+1):need(state(gen[:i])[0]>1,'every generator prefix expands')
lw,Pw,Dw=state(W);A,Q,B,C,E=affine(W)
need(cyl(W)==(63,2**55),'source63 complete exact guard')
need(replay(63,W)[0]==61,'source63 original endpoint')
need(Dw==F(215402584767891,140737488355328),'source63 exact reset debt')
need(2-Dw==F(66072391942765,140737488355328),'source63 exact payment margin')
for t in (0,1,2,7):
    n=63+2**55*t;h,fo,ft,fh=replay(n,W)
    need(h==61+2*3**34*t,'source63 cylinder endpoint lift')
    need(fo[:3]==(34,54,88) and ft[:2]==(54,88),'source63 first exits retain actual clocks')
    need(F(C,Q)+F(A,Q)*(n-1)<n-1,'source63 all tested cylinder lifts paid')
root_trace=replay(63,W+ROOT_SUFFIX,True)
need(root_trace[0]==1,'literal63 ROOT grounding')
need((len(W+ROOT_SUFFIX),sum(W+ROOT_SUFFIX))==(39,68),'literal63 ROOT clocks')
# A fixed ROOT word grounds only its one original source, not its other cylinder lifts.
r,mod=cyl(W+ROOT_SUFFIX)
need(r==63 and mod==2**69,'literal63 ROOT complete cylinder')
need(replay(63+mod,W+ROOT_SUFFIX,True)[0]==1+2*3**39,'ROOT cylinder lift is not ROOT')

# An infinite, explicitly defined grounded family through the original63.
# These are NEW sources, not replacements for125-power family members.
H=2*3**33
need(pow(2,H,3**34)==1 and pow(2,H//3,3**34)!=1 and pow(2,H//2,3**34)!=1,'grounded final-valuation exact period')
grounded_modular=[]
for t in (0,1,2,7,31):
    exp=54+H*t
    need((61*pow(2,exp,3**34)-B)%(3**34)==0,'grounded family source integer')
    denom=3**34;num=(61*pow(2,exp,105*denom)-B)%(105*denom)
    need(num%denom==0,'grounded source observer exact division')
    source105=(num//denom)%105
    grounded_modular.append({'t':t,'final_valuation':4+H*t,'source_mod15':source105%15,'source_mod35':source105%35,'first_literal_shortcut_clock':54,'first_literal_ordinary_clock':88,'first_odd_shortcut_clock':54+H*t,'ROOT_shortcut_clock':68+H*t,'ROOT_ordinary_clock':107+H*t})
need(61*(2**54-3**34)-B==2*3**34,'grounded family first literal exit positive margin base')
need(3**34>2**53,'grounded family no earlier final-run literal exit')

# Closure controls: every mix of up to three generators, two separator lengths.
macro_rows=[]
for k in (1,2,3):
    for kinds in product((0,1),repeat=k):
        for gaps in product((1,2),repeat=k-1):
            ww=()
            for i,q in enumerate(kinds):
                if i:ww+=(1,)*gaps[i-1]
                ww+=(U,V)[q]
            ll,pp,dd=state(ww)
            need(dd<=c*(ll-theta),'separator-cone closure')
# Exact equality on the unbounded symbolic one-generator scheme, finite controls.
for k in (0,1,2,3,7,15,31,63):
    prefix=(U+(1,))*k+U
    lp,pp,dp=state(prefix)
    need(dp==c*(lp-theta),'repeated source63 macro exact cone boundary')
    a=2
    while F(3,2**a)*lp>=1:a+=1
    word=prefix+(a,);aa,qq,bb,cc,ee=affine(word);r,mod=cyl(word)
    n=r
    if (qq-aa)*n<=bb:
        n+=mod*((bb-(qq-aa)*n)//((qq-aa)*mod)+1)
    h,fo,ft,fh=replay(n,word,True)
    need(h<n and fo[0]==len(word),'macro first odd exit')
    need(F(ee,qq)<1+c and F(ee,qq)<n-h,'macro positive charge paid')
    need(n-h-F(ee,qq)>1-c,'macro uniform absolute payment margin')
    macro_rows.append({'copies':k+1,'odd_steps':len(word),'shortcut_clock':sum(word),'ordinary_clock':sum(word)+len(word),'last_valuation':a,'debt_resets':sum(x>=3 for x in word),'source':str(n),'endpoint':str(h),'first_literal_shortcut':ft[0],'first_literal_ordinary':ft[1],'mod15_mod35_source':[n%15,n%35],'mod15_mod35_endpoint':[h%15,h%35]})
# Dropping the separator destroys this particular closed cone.
need(state(U+U)[2]>c*(state(U+U)[0]-theta),'no-separator hostile')
# D>=lambda is a sharp obstruction to every finite all-one separator repair.
need(state((3,1,1,1,1))[2]>=state((3,1,1,1,1))[0],'nonpositive-surplus hostile')
# Retain native first-exit versus later-descent hostility.
w55,h55=actual_word(55,lambda x,w:x==53)
need(len(w55)==39 and state(w55)[2]>55-53,'later-descent source55 unpaid')
w4575,h4575=actual_word(4575,lambda x,w:x<4575)
need(state(w4575)[2]>2,'global debt below two is false')

# Frozen inherited source certificates, not newly substituted representatives.
input_path=HERE/'inputs/guarded_children_105_20261018g.json'
raw=input_path.read_bytes();data=json.loads(raw);rows=data['rows']
need(len(rows)==63,'exact inherited 63 template universe')
thresholds=[];root_rank=root_cost=0
for row in rows:
    j=row['j'];n=125**j;v=tuple(row['prefix']);c0=row['c_min'];K=row['guard_bits']
    l0,p0,d0=state(v);a0,q0,b0,_,_=affine(v)
    need(K==sum(v)+c0-2,'variable guard precision')
    need(int(row['P'])==3*a0 and int(row['B'])==3*b0+q0,'template carries inherited exactly')
    lm,pm,dm=state(v+(c0,));R=1+pm/(1-lm)
    need(lm<F(1,2),'template strict coefficient half contraction')
    need(R<9 and R<n,'all-template positive threshold below nine')
    need(dm<2,'all-template minimum-final-valuation debt below two')
    need(max(F(1),dm)<2,'all-final-valuation supremum debt below two')
    for cc in (c0,c0+1,c0+12,c0+64):
        lx,px,dx=state(v+(cc,))
        need(1+px/(1-lx)<=R,'variable final valuation lowers barrier')
        need(dx<2,'variable final valuation preserves debt bound')
        # Compile exact c, including the final odd bit, without materializing 125^J.
        bits=sum(v)+cc-1;J=0;valuation=sum(v)+cc;PP=int(row['P']);BB=int(row['B'])
        for bit in range(bits):
            mm=1<<(bit+3);target=(1<<valuation)%mm
            candidates=[t for t in (J,J+(1<<bit)) if (PP*pow(125,t,mm)+BB)%mm==target]
            need(len(candidates)==1,'unique exact valuation exponent lift')
            J=candidates[0]
        need((J-j)%(1<<K)==0,'exact final valuation retains original exponent guard')
        modulus=1<<(valuation+1);xx=pow(125,J,modulus)
        for a in v+(cc,):
            yy=3*xx+1
            need((yy&-yy).bit_length()-1==a,'modular guard replays literal prefix')
            xx=yy>>a
        observer=[]
        for m in (15,35):
            hh=(PP*pow(125,J,m)+BB)*pow(2,-valuation,m)%m
            # Independent denominator-preserving arithmetic.
            num=PP*pow(125,J,m*(1<<valuation))+BB
            need(num%(1<<valuation)==0,'denominator precision retained')
            need(hh==(num//(1<<valuation))%m,'native register transport')
            observer.append(hh)
        need(observer[0]%5==observer[1]%5,'mod15 mod35 compatibility')
        h105=(observer[0]+15*((observer[1]-observer[0])%7))%105
        need(h105==(PP*pow(125,J,105)+BB)*pow(2,-valuation,105)%105,'CRT decoded native child')
    packed=zlib.decompress(base64.b85decode(row['word_zlib_base85']));word=tuple(packed)
    need(hashlib.sha256(packed).hexdigest()==row['word_sha256'],'retained ROOT word hash')
    need(word[:len(v)]==v and word[len(v)]==row['actual_c'],'retained source and child prefix')
    x=n;child=None;first_half_index=None;clock=0;first_half_T=None
    for i,a in enumerate(word,1):
        need(x>1,'no early ROOT in named source')
        y=3*x+1
        need((y&-y).bit_length()-1==a,'retained ROOT exact native valuation')
        for q in range(a):
            y//=2;clock+=1
            if first_half_T is None and 2*y<=n:first_half_T=clock
        x=y
        if first_half_index is None and 2*x<=n:first_half_index=i
        if i==len(v)+1:child=x
    need(x==1,'named original source reaches ROOT')
    need(child==int(row['child']),'named original child preserved')
    need(first_half_index==len(v)+1==row['first_odd_half_index'],'named first odd half clock')
    need([child%15,child%35]==row['child_residues'],'named child prime registers')
    tail=packed[len(v)+1:]
    need(hashlib.sha256(tail).hexdigest()==row['child_word_sha256'],'retained child ROOT suffix hash')
    need((len(tail),sum(tail))==(row['child_rank'],row['child_cost']),'retained child ROOT clocks')
    # Independent backwards replay uses the exact same word, not a terminating representative.
    x=1
    for a in reversed(word):
        num=(x<<a)-1;need(num%3==0,'ROOT reverse integrality');x=num//3
        need(x%2==1,'ROOT reverse oddness')
    need(x==n,'reverse ROOT word recovers original 125 power')
    need((len(word),sum(word))==(row['word_rank'],row['word_cost']),'named source ROOT clocks')
    root_rank+=len(word);root_cost+=sum(word)
    thresholds.append({'j':j,'prefix_odd_steps':len(v),'prefix_clock':sum(v),'c_min':c0,'actual_seed_c':row['actual_c'],'guard_bits':K,'barrier':str(R),'debt_at_c_min':str(dm),'physical_barrier':str(F(int(row['B']),2**(sum(v)+c0)-int(row['P']))),'first_literal_half_shortcut_seed':first_half_T,'odd_half_shortcut_seed':sum(v)+row['actual_c'],'ROOT_odd_steps':len(word),'ROOT_shortcut_clock':sum(word),'ROOT_ordinary_clock':sum(word)+len(word)})
need((root_rank,root_cost)==(575652,1152705),'all63 inherited total ROOT clocks')
peak=max(thresholds,key=lambda x:F(x['barrier']));peakD=max(thresholds,key=lambda x:F(x['debt_at_c_min']))
need(peak['j']==638 and peakD['j']==638,'sharp common barrier and debt index')
need(next(r for r in rows if r['j']==638)['actual_c']==8,'sharp barrier is physically attained')
# Explicit infinite exponent complement: an all-ones tail after T^2.
L=max(r['prefix_cost']+r['c_min'] for r in rows);J=0
for bit in range(L):
    mm=1<<(bit+3);sol=[s for s in (J,J+(1<<bit)) if (3*pow(125,s,mm)+5)%mm==0]
    need(len(sol)==1,'all-ones exponent complement bit lift');J=sol[0]
need(L==149,'maximum uniform template half clock')
need(all((J-r['j'])%(1<<r['guard_bits']) for r in rows),'infinite exponent progression misses every template guard')
need((3*pow(125,J,1<<(L+2))+5)%(1<<(L+2))==0,'all-ones complement exact congruence')

out={'status':'PROVED scoped separator-cone and variable-terminal implications; FINITE-EXACT controls; universal coverage and ROOT energy OPEN','checks':COUNTS,'total_checks':sum(COUNTS.values()),'source63':{'prefix':U,'first_exit_word':W,'ROOT_suffix':ROOT_SUFFIX,'lambda_prefix':str(l),'P_prefix':str(P),'D_prefix':str(D),'ordinary_cone_min_c':str(D/(l-1)),'separator_theta':str(theta),'sharp_separator_c':str(c),'uniform_family_margin_lower_bound':str(1-c),'first_exit_debt':str(Dw),'first_exit_payment_margin':str(2-Dw),'physical_first_exit_barrier':str(F(B,Q-A)),'positive_first_exit_barrier':str(1+Pw/(1-lw)),'first_exit_odd_clock':34,'first_exit_shortcut_clock':54,'first_exit_ordinary_clock':88,'ROOT_odd_clock':39,'ROOT_shortcut_clock':68,'ROOT_ordinary_clock':107,'source_registers':[3,28],'first_exit_registers':[1,26]},'grounded_source63_family':{'valuation_period':H,'source_formula':'(61*2^(54+H*t)-48215861999407177)/3^34, t>=0','controls':grounded_modular},'macro_controls':macro_rows,'templates':{'input_sha256':hashlib.sha256(raw).hexdigest(),'count':63,'root_rank':root_rank,'root_shortcut_cost':root_cost,'root_ordinary_cost':root_rank+root_cost,'sharp_common_barrier':peak['barrier'],'largest_debt_at_c_min':peakD['debt_at_c_min'],'rows':thresholds},'uncovered_exponent_progression':{'residue':str(J),'bits':L,'modulus':str(1<<L),'guaranteed_growth_tail_shortcut_steps':L,'no_half_descent_through_shortcut_clock':L+2}}
print(json.dumps(out,indent=2,sort_keys=True))
