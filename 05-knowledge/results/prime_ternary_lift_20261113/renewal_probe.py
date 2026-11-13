"""Independent certificates for higher-mode moves and the finite-phase boundary."""
from pathlib import Path
from itertools import product
from math import comb,factorial
import gzip,json
import numpy as np
from core_moves import skeleton,embedding,physical,SOURCE,PLUS2,PLUS4,apex
from frame_probe import ordinary_bank,clone,hp,valid,strong,need
from visit_probe import visits,covers_tt

HERE=Path(__file__).resolve().parent

def power(v):
    return tuple(sum(a*(-1)**(k-j)*comb(k,j) for k,a in enumerate(v,1) if k>=j) for j in range(1,len(v)+1))

def inflate(T,r,M):
    old=[v for v in range(len(T)) if v!=r];labels=old+[r]*len(M);s=len(old)
    return tuple(sum(int(i!=j and ((a==b and M[i-s]>>(j-s)&1) or (a!=b and T[a]>>b&1)))<<j for j,b in enumerate(labels)) for i,a in enumerate(labels))

def run():
    controls=0
    for (n,code),T in sorted(ordinary_bank().items()):
        if n>5:continue
        for r in range(n):
            v=tuple(map(int,visits(np.array(T,dtype=np.int64),r)))
            need(v==tuple(skeleton(T,r,j) for j in range(1,len(v)+1)),'independent permutation-slots certificate')
            controls+=1
    bank=json.loads(gzip.decompress((HERE/'visits.json.gz').read_bytes()))
    chosen=[t for t in bank['transitions'] if (t['source']['profile'][0],t['digit']) in ((25,4),(41,2)) and t['max_runs']==3]
    chosen+= [t for t in bank['transitions'] if t['max_runs']==4]
    records=[]
    for t in chosen:
        A=t['source'];B=t['target'];digit=t['digit']
        for rec in (A,B):
            v=tuple(skeleton(rec['rows'],rec['root'],j) for j in range(1,5))
            need(v==tuple(rec['profile']),'literal higher-mode certificate')
            need(all(a%factorial(j-1)==0 for j,a in enumerate(v,1)),'excursion factorial law')
        E=embedding(A['rows'],A['root'],B['rows'],B['root']);checks=[]
        for m in range(1,6):
            a=clone(A['rows'],A['root'],m);b=physical(A['rows'],A['root'],B['rows'],B['root'],m,E)
            need(valid(a) and valid(b) and strong(a) and strong(b),'physical strong higher-mode graphs')
            x,y=hp(a),hp(b)
            need(y==3*x+digit,'exact higher-mode target')
            need(x==sum(c*j**m for j,c in enumerate(power(A['profile']),1)),'power-basis source identity')
            checks.append([m,x,y])
        M=(2,4,1)
        x=hp(inflate(A['rows'],A['root'],M));y=hp(inflate(B['rows'],B['root'],M))
        need(y==3*x+3*digit and y!=3*x+digit,'hostile cyclic substitution scales digit')
        records.append(dict(digit=digit,source=A,target=B,source_power=power(A['profile']),embedding=E,checks=checks,cyclic_control=[x,y]))
    # Universal identities do not license replacing a transitive module by C3.
    x=hp(inflate(*SOURCE,(2,4,1)))
    ys=[hp(inflate(*B,(2,4,1))) for B in (PLUS2,PLUS4)]
    need((x,ys)==(177,[537,543]),'shared-source cyclic hostile')
    # Check the geometric description of the one-arc +4 construction.
    geometry=[]
    for m in range(1,6):
        A=list(apex('1001'+'0'*m));A[1]^=1<<3;A[3]^=1<<1
        need(valid(A) and hp(tuple(A))==15*2**m-1,'back-arc drawing of source')
        middle=hp(apex('1001'+'0'*m));last=hp(apex('1001'+'0'*m+'110'))
        need(middle==9*2**m-3 and last==45*2**m+1,'straighten then append')
        geometry.append([m,hp(tuple(A)),middle,last])
    # Known positive control for a repeated digit2 clock.
    T=(206,252,56,240,225,193,132,4);r=3
    v=tuple(map(int,visits(np.array(T,dtype=np.int64),r)))
    need(v==(179,180,60,0),'inherited amplitude60 clock')
    clock=[hp(clone(T,r,m)) for m in range(1,6)]
    need(clock==[60*3**m-1 for m in range(1,6)],'repeated +2 control')
    need(power((39,36,9))==(-6,9,9) and 9%factorial(2)!=0,'halved mixed family violates excursion divisibility')
    halving=[]
    for m in range(1,9):
        H=6*3**m+6*2**m-5;r=(H+1)//2;child=3*H+4;Q=(child-1)//2
        need(r%2==1 and (child+1)//2==3*r+1 and Q==9*3**m+9*2**m-6 and Q%2==1,'typed plus branch and first halving')
        halving.append([m,H,child,Q])
    known_witnesses={};wanted={25,79,39,19,9}
    for key,T in sorted(ordinary_bank().items()):
        h=hp(T)
        if h in wanted:known_witnesses.setdefault(h,T)
        if wanted.issubset(known_witnesses):break
    need(wanted.issubset(known_witnesses),'individual small Collatz-chain counts attained with actual witnesses')
    # Arithmetic cycle check: odd constant term plus excursion divisibility.
    cycle_controls=0
    for p in range(1,13):
        for digits in product((2,4),repeat=p):
            D=0
            for d in digits:D=3*D+d
            den=3**p-1
            if D%den==0 and (D//den)%2==1:
                need(set(digits)=={2} and D==den,'finite-phase cycle has only digit2')
            cycle_controls+=1
    # Fixed-module chains have exact entry when the entire output vector matches.
    nexts={}
    for t in bank['transitions']:nexts.setdefault(tuple(t['source']['profile']),[]).append(t)
    best={}
    for v in sorted(nexts,reverse=True):
        cand=[]
        for t in nexts[v]:cand.append([t]+best.get(tuple(t['target']['profile']),[]))
        best[v]=max(cand,key=len)
    chain=max(best.values(),key=len)
    result=dict(status='FINITE-EXACT certificates for stated all-module identities and independent hostile controls',scope='All ordinary frames through5 for independent visit DP; selected higher-mode transitions through8; literal modules m1..5 and C3; all 8190 digit words of lengths1..12. Finite-phase theorem is proved separately, not inferred from words.',independent_root_controls=controls,higher_mode_moves=records,cyclic_hostile=[x,ys],geometry=geometry,known_clock_profile=[179,180,60,0],known_clock_counts=clock,halving_controls=halving,halving_impossible_profile=[39,36,9],small_chain_witnesses=known_witnesses,cycle_word_controls=cycle_controls,longest_same_module_chain_length=len(chain),longest_chain=chain)
    suffix='' if __debug__ else '.optimized'
    (HERE/('renewal'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('INDEPENDENT_ROOTS',controls,'HIGHER_MOVES',len(records),'CYCLE_WORDS',cycle_controls,flush=True)
    print('FAMILIES',[(r['source_power'],r['digit'],r['checks'],r['embedding']) for r in records],flush=True)
    print('LONGEST_CHAIN',len(chain),[(t['source']['profile'],t['digit'],t['target']['profile']) for t in chain],flush=True)

if __name__=='__main__':run()
