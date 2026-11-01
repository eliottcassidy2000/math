"""Literal table replay, prime-integer FRACTRAN, and formula crosschecks."""
from pathlib import Path
from math import prod,comb
from itertools import product
import hashlib,json,sys
import inherited_linear as linear
from machines import *
HERE=Path(__file__).resolve().parent

def read(name):return json.loads((HERE/name).read_text())

def check_table(m):
    seen=set();states=set(m['states']);alphabet=set(m['alphabet']);k=m['tapes']
    for q,s,t,w,mv in m['rows']:
        need(q in states and t in states,'unlisted state')
        need((q,s) not in seen,'nondeterministic table');seen.add((q,s))
        if m['model']=='packed_single':need(s in alphabet and w in alphabet and mv in 'LRS','single table row')
        else:need(len(s)==len(w)==len(mv)==k and set(s+w)<=alphabet and set(mv)<=set('LRS'),'multi table row')

def compare_round(p,tm,rows,v,clock=None):
    data=p['register_primes'];controls=p['control_primes'];q=p['initial_state'];v=list(v)
    sim=LiteralTM(tm,q,v);n=encode(q,v,data);expected_time=0;fsteps=0;isteps=0
    while True:
        next_n,k=fraction_step(n,rows);need(next_n is not None,'unexpected FRACTRAN halt')
        # IR execution until the next FRACTRAN boundary, independent of table.
        while True:
            q,v,c=ir_step(q,v,p['instructions']);expected_time+=c;isteps+=1
            if q in p['boundaries']:break
        nq,nv,dt=sim.next_boundary();need((nq,nv)==(q,v),'TM/IR boundary mismatch')
        need(sim.time==expected_time,'exact unary opcode cost')
        need(encode(nq,nv,data)==next_n,'literal integer/TM mismatch')
        n=next_n;fsteps+=1
        if q==p['initial_state']:break
    return v,fsteps,sim.time,isteps

def blue_cost(v):
    a,b,c,d,e,f,_,n=v
    return 4*a*a+5*e*e+6*a*e+2*(b*b+c*c+d*d+f*f)+14*a+16*e+7*(b+c+d+f)+2*n+24

def matrix_cost(a,x,n):
    y=[sum(c*v for c,v in zip(row,x)) for row in a]
    return sum(v*v for v in x)+3*sum(v*v for v in y)+5*sum(x)+9*sum(y)+6*len(x)+2*n+4

def raw_program(name):return read(name+'.ir.json'),load_tm(HERE/(name+'.tm.json')),load_fractions(HERE/(name+'.fractran'))

def blue_tests():
    p,tm,rows=raw_program('pure_blue');check_table(tm)
    cases=0;max_cost=0
    for vals in product(range(3),repeat=6):
        v=list(vals)+[0,3];out,f,t,i=compare_round(p,tm,rows,v)
        need(out==list(vals[1:])+[vals[0]+vals[4],0,4],'blue update')
        need(f==8+2*(2*vals[0]+sum(vals[1:4])+2*vals[4]+vals[5]),'blue fraction round cost')
        need(t==blue_cost(v),'blue TM round cost');cases+=1;max_cost=max(max_cost,t)
    # Actual stream: initialize once, retain the physical tape throughout.
    sim=LiteralTM(tm,p['initial_state'],p['initial_vector']);n=encode(p['initial_state'],p['initial_vector'],p['register_primes'])
    checkpoints=[];fsteps=0;last_time=0;last_f=0;old=p['initial_vector']
    def narayana(m):return sum(comb(m-2*k,k) for k in range(m//3+1)) if m>=0 else 0
    def value(n):return narayana(n//2) if n%2==0 else narayana(n//2+2)
    for target in range(25):
        if target:
            while True:
                n,k=fraction_step(n,rows);fsteps+=1
                q,v,t=sim.next_boundary()
                need(encode(q,v,p['register_primes'])==n,'blue continuing tape')
                if q==p['initial_state']:break
            need(sim.time-last_time==blue_cost(old),'stream exact TM formula')
        v=sim.values();need(v==[value(target+j) for j in range(6)]+[0,target],'independent binomial readout')
        checkpoints.append({'n':target,'output':v[0],'window':v[:6],'FRACTRAN_steps':fsteps,'TM_transitions':sim.time})
        old=v;last_time=sim.time;last_f=fsteps
    return {'exhaustive_windows':cases,'range':'{0,1,2}^6, scratch=0, clock=3','max_one_round_TM_transitions':max_cost,'checkpoints':checkpoints}

def linear_tests():
    results=[]
    for name,rounds in [('compact_tail21',6),('migration_full',1),('migration_face',2)]:
        p,tm,rows=raw_program(name);check_table(tm);a=p['matrix'];d=len(a);v=p['initial_vector'];rec=[]
        for turn in range(rounds):
            x=v[:d];y=[sum(c*z for c,z in zip(row,x)) for row in a]
            out,f,t,i=compare_round(p,tm,rows,v)
            need(out==y+[0]*d+[v[-1]+1],'matrix output')
            need(f==1+2*d+2*sum(x)+2*sum(y),'matrix fraction cost')
            need(t==matrix_cost(a,x,v[-1]),'matrix TM cost')
            rec.append({'input':x,'output':y,'clock_in':v[-1],'FRACTRAN_steps':f,'IR_steps':i,'TM_transitions':t})
            v=out
        if name=='migration_face':
            for k,r in enumerate(rec,1):need(2*r['output'][0]+r['output'][1]==(52*80**k+27)//79,'radix80 output')
        results.append({'program':name,'rounds':rec})
    # Matrix tests independent of persisted specialized lists, including zero,
    # nilpotent, permuting, dense and coefficient>1 cases.
    cases=[([[0]],[0]),([[0]],[3]),([[1]],[2]),([[3]],[2]),([[0,1],[1,0]],[2,3]),
           ([[2,3],[1,0]],[2,1]),([[0,0],[0,0]],[2,3]),
           ([[0,0,2],[1,0,0],[0,3,0]],[1,2,3])]
    for a,x in cases:
        f,_=linear.compile_linear(a);d=len(x)
        p=to_ir(f['fractions'],f['x']+f['y']+[f['clock']],f['control']+[f['checkpoint']]);p['initial_state']='C'+str(f['checkpoint'])
        tm=to_tm(p,'unary_multi');out,fs,t,is_=compare_round(p,tm,f['fractions'],x+[0]*d+[2])
        need(t==matrix_cost(a,x,2),'generic exact TM cost')
    return {'actual_programs':results,'generic_matrix_cases':len(cases)}

def reverse_tests():
    count=0
    for name in ['pure_blue','compact_tail21','migration_full','migration_face']:
        p=read(name+'.ir.json');b=read(name+'.roundtrip.json');ir=p['instructions'];d=len(b['register_primes'])
        need(b==ir_to_fractran(p),'serialized reverse compiler')
        for q,op in ir.items():
            for level in range(3):
                v=[level]*d;target,w,c=ir_step(q,v,ir)
                n=b['controls'][q]*prod(pr**e for pr,e in zip(b['register_primes'],v))
                n1,k=fraction_step(n,b['fractions']);n2,k2=fraction_step(n1,b['fractions'])
                need(n2==b['controls'][target]*prod(pr**e for pr,e in zip(b['register_primes'],w)),'two-step reverse compiler')
                count+=1
    # Explicit self-loop, zero branch and halt controls: no control cancellation.
    p={'register_primes':[2],'instructions':{'A':['inc',0,'A'],'B':['decjz',0,'B','B'],'H':['halt']}}
    b=ir_to_fractran(p)
    for q in ['A','B']:
        for e in range(4):
            n=b['controls'][q]*b['register_primes'][0]**e
            n1,_=fraction_step(n,b['fractions']);n2,_=fraction_step(n1,b['fractions'])
            expected=e+1 if q=='A' else max(0,e-1)
            need(n2==b['controls'][q]*b['register_primes'][0]**expected,'self-loop trampoline')
    need(fraction_step(b['controls']['H'],b['fractions'])==(None,None),'halt control')
    return {'ordinary_instruction_cases':count,'self_loop_cases':8,'halt_control_checked':True}

def hostile_tests():
    rejected=0
    examples=[([(3,2)],[2],[3]), # missing denominator control
              ([(6,10)],[2],[3,5]), # unreduced
              ([(7,5),(11,10)],[2],[5,7,11]), # fallback preceding guard
              ([(7,20),(11,5)],[2],[5,7,11]), # exponent-two guard
              ([(7,6),(11,3)],[2],[3,7,11,2])] # overlap
    for args in examples:
        try:to_ir(*args)
        except ValueError:rejected+=1
        else:raise ValueError('malformed compiler input accepted')
    # Product integrality after reduction differs from raw-denominator divisibility.
    need(fraction_step(1,[(2,2)])==(1,0),'cancellation hostile')
    return {'compiler_inputs_rejected':rejected,'unreduced_2_over_2_hostile':True}

def general_tests():
    cases=[([(5,6),(7,2)],2), ([(2,2)],1), ([(3,2),(5,2)],2),
           ([(3,2),(5,2)],1), ([(1,12),(5,2)],4), ([(1,1)],7),
           ([],11), ([(9,4),(2,3)],12), ([(9,4),(2,3)],3)]
    records=[]
    for rows,n in cases:
        p=general_fractran_ir(rows);tm=to_tm(p,'unary_multi');check_table(tm)
        v=[];inert=n
        for pr in p['register_primes']:
            e=0
            while inert%pr==0:e+=1;inert//=pr
            v.append(e)
        sim=LiteralTM(tm,'STEP',v);expected,k=fraction_step(n,rows)
        q='STEP';cost=0;steps=0
        if p['instructions'][q][0]!='halt':
            while True:
                q,v,c=ir_step(q,v,p['instructions']);cost+=c;steps+=1
                if q in p['boundaries']:break
            tq,tv,t=sim.next_boundary();need((tq,tv,t)==(q,v,cost),'generic table/cost')
        if expected is None:
            need(p['instructions'][q][0]=='halt' and not sim.step(),'generic halt')
            need(inert*prod(pr**e for pr,e in zip(p['register_primes'],v))==n,'failed guard rollback')
        else:
            need(q=='STEP' and inert*prod(pr**e for pr,e in zip(p['register_primes'],v))==expected,'generic rational semantics')
        records.append({'fractions':rows,'input':n,'output':expected,'chosen_index_0based':k,'IR_steps':steps,'TM_transitions':cost})
    return records

def main():
    output={'status':'PASS','pure_blue':blue_tests(),'linear':linear_tests(),'reverse':reverse_tests(),'hostiles':hostile_tests(),'general_FRACTRAN':general_tests(),
            'scope':'All-input claims are proved in PROOF.md; these finite runs validate emitted tables and receipts.'}
    output['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.fractran') or p.name.endswith(('.tm.json','.ir.json','.roundtrip.json')) or '.tm.rows.' in p.name}
    print(json.dumps(output,sort_keys=True,indent=2))
if __name__=='__main__':main()
