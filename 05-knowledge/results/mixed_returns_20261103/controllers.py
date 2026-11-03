"""Actual graphs, native programs, phase transport and hostile controls."""
from itertools import permutations,combinations
import importlib.util,json,sys,sympy as sp
from support import *

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def chart(n,code):
    rows=[0]*n
    for i in range(n-1):rows[i]|=1<<(i+1)
    for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+2,n)):
        u,v=(i,j) if code>>k&1 else(j,i);rows[u]|=1<<v
    return rows

def pointed(rows,port):
    h=e=d=s=0
    for p in permutations(range(len(rows))):
        bad=[k for k,(u,v) in enumerate(zip(p,p[1:])) if not(rows[u]>>v&1)]
        if not bad:h+=1;e+=p[-1]==port;s+=p[0]==port
        elif len(bad)==1 and p[bad[0]+1]==port:d+=1
    c=h+d;need((c-e)%2==0,'integer pointed coordinates')
    return (h,(c-e)//2,(c+e)//2),s

def phase(v):
    h,x,y=v;return ((4*(h+2*x-y))%7,y%7,(h-y)%7)

def main():
    need(chart(6,78)==ATLAS['s78']['source_rows'],'independent source-chart convention')
    for code in (18413,12083):ATLAS['t'+str(code)]={'source_rows':chart(7,code)}
    input_words=list(SEEDS)+[('a','t18413'),('b','a','a','r21'),('r2','r30'),('a','t12083')]
    states=[vector(w) for w in SEEDS]+list(REPAIRS);inputs=[]
    for w,v in zip(input_words,states):
        rows=actual_target(w);got,s=pointed(rows,len(rows)-1)
        need(got==v and s>0,'every supplied state has an actual positive-start port')
        inputs.append({'word':w,'state':v,'rows':rows,'port':len(rows)-1,'start_count':s})
    native=module('native_compiler',HERE.parent/'ternary_repair_clock_20261103/native.py')
    carry=json.loads((HERE/'carry_complete.json').read_text());programs=[]
    positive=[r for r in carry['results'] if r['answer']]
    for record in positive:
        v=record['input_state'];source=inputs[states.index(tuple(v))];answer=record['answer'];u=answer['output_state']
        a=phase(v);b=phase(u);need(b==((3*a[0])%7,(3*a[1])%7,(3*a[2]+2*record['digit'])%7),'full field/charge branch transport')
        target=tuple(answer['target_rows']);need(all(target[i]>>(i-1)&1 for i in range(1,len(target))),'output port spine')
        program=native.compile_repair(source['rows'],target)
        programs.append({'input_port':source['port'],'output_port':len(target)-1,
            'input_state':v,'output_state':u,'carry':answer['carry'],'old_names_retained':True,'program':program})
    # Choosing the repaired H293 representative removes the next carry entirely.
    new293=positive[0]['answer'];out881=positive[1]['answer']
    need(tuple(out881['output_state'])==tuple(3*x+t for x,t in zip(new293['output_state'],T)),'carry-aware source selection gives an exact next branch')
    program=native.compile_repair(new293['target_rows'],out881['target_rows'])
    programs.append({'input_port':8,'output_port':10,'input_state':new293['output_state'],
        'output_state':out881['output_state'],'carry':[0,0,0],'old_names_retained':True,'program':program})
    # Independent literal replay of every emitted birth/complement/deletion command.
    def complement(rows,path):
        n=len(rows);need(sorted(path)==list(range(n)) and all(rows[a]>>b&1 for a,b in zip(path,path[1:])),'replayed legal path')
        keep={frozenset((a,b)) for a,b in zip(path,path[1:])};out=list(rows)
        for a,b in combinations(range(n),2):
            if frozenset((a,b)) not in keep:out[a]^=1<<b;out[b]^=1<<a
        return tuple(out)
    for rec in programs:
        p=rec['program'];rows=tuple(p['source_rows'])
        for _ in range(p['padding_births']):rows=tuple(r|(1<<len(rows)) for r in rows)+(0,)
        for command in p['commands']:
            n=len(rows);out=command['birth_out']
            rows=tuple(r|(0 if out>>i&1 else 1<<n) for i,r in enumerate(rows))+(out,)
            for path in command['paths']:rows=complement(rows,path)
            rows=tuple(r&((1<<n)-1) for r in rows[:n])
        n=len(p['target_rows']);rows=tuple(r&((1<<n)-1) for r in rows[:n])
        need(rows==tuple(p['target_rows']),'independent full native program replay')
    A=('a','r41');B=('b','r25');need(vector(B)==tuple(3*x-t for x,t in zip(vector(A),T)),'new actual offset pair')
    family=[]
    for k in range(4):
        w=('q2',)*k+('b',);m=vector(A+w)[0];n=vector(B+w)[0]
        need(n==3*m+2,'infinite offset identity sample')
        for prefix,value in ((A,m),(B,n)):need(literal_frontier(actual_target(prefix+w))[0]==value,'independent literal family count')
        family.append({'k':k,'base':m,'companion':n,'words':[A+w,B+w]})
    q=sp.Matrix(MATS['q2']);b=sp.Matrix(MATS['b']);ell=sp.Matrix([1,0,-1])
    need((q*q-18*q+5*sp.eye(3))*b[:,0]==12*ell and q*ell==ell,'all-order new family recurrence')
    for k in range(2):need(family[k+2]['base']==18*family[k+1]['base']-5*family[k]['base']+24,'recurrence controls')
    # Reuse the actual hostile constructors; do not replace their graphs by H.
    prior=HERE.parent/'primegame_port_spectrum_20261103';sys.path.insert(0,str(prior))
    hm=module('actual_hostile',prior/'verify.py');hostile=hm.hostile()
    need(hostile['groups']==[81,27] and hostile['port8_after_sink']==[2390013,1754919],'81/27 hostile')
    old=json.loads((HERE.parent/'inverse65_zeckendorf_20261103/branches.json').read_text())
    counter=module('subset_counter',HERE.parent/'golden_hamilton_transfer_20261103/transfer.py')
    for source in old['H55_sources']:
        rows=tuple(source['rows']);values=[]
        for out in range(1<<len(rows)):
            child=tuple(r|(0 if out>>i&1 else 1<<len(rows)) for i,r in enumerate(rows))+(out,)
            values.append(counter.subset_h(child))
        need(111 not in values and sorted(values)==sorted(source['attachment_H']),'all named one-birth55 controls')
    receipt=json.loads((PRIOR/'birth_verification.json').read_text())
    for r in receipt['positive_extensions']:need(counter.subset_h(tuple(r['rows']))==111,'two/three-birth external positives')
    result={'status':'PASS','input_sources':inputs,'native_programs':programs,'offset_family':family,
        'hostile81_27':hostile,'one_birth55_orientations':768,'external_positive_extensions':len(receipt['positive_extensions']),
        'three_birth_exclusion':'Inherited exact prior certificate, pinned; no redundant full 2^24 rerun.'}
    (HERE/'controllers.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':'PASS','program_moves':[p['program']['native_moves'] for p in programs],
        'program_peaks':[p['program']['peak_order'] for p in programs],'family':[(r['base'],r['companion']) for r in family],'hostile':hostile}))
if __name__=='__main__':main()
