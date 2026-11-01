"""Independent cycle-based audit of pointed one-apex fibres and square shells."""
from collections import Counter
from functools import lru_cache
import itertools
import importlib.util
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).parent))
from audit_compiler import check,ir_run,raw_run
from audit_seven_counters import load_program,load_source
from audit_path_counter import emitted_run

def word(n):return [int(c) for c in bin(n)[2:]]
def cycles(bits):
    return sum((1 if j==i+1 else 2**(j-i-2))
               for i in range(len(bits)) for j in range(i+1,len(bits))
               if bits[i]==1 and bits[j]==0)
def H(bits):return 1+2*cycles(bits)
@lru_cache(None)
def f(n):
    if n<=1:return 1
    q,r=divmod(n,2)
    return f(q) if r else f(q)+q+(q%2)

def main():
    for n in range(4097):check(f(n)==H(word(n)),('cycle/recurrence',n))
    counts=Counter();fixed=Counter()
    for n in range(1,2048):
        w=word(2*n);h=H(w);counts[h]+=1
        if w==[1-b for b in reversed(w)]:fixed[h]+=1
    shells=[]
    for q in range(1,4097):
        u=word(q);w=u+[1-b for b in reversed(u)]
        g=q*q+(q%2)+2*f(q)-1
        check(H(w)==g,('cross-cycle formula',q))
        check(q*q<g<=(q+1)**2,('square shell',q))
        equality=q>=2 and q&(q-1)==0
        check((g==(q+1)**2)==equality,('shell equality',q))
        if g<=2047:shells.append(g)
    for h in range(1,2048,2):
        # n>=1 is essential: the tempting n=0 term would falsely count R_1.
        inverse=sum(f(n)+n+(n%2)==h for n in range(max(1,(h-1)//2),h))
        check(inverse==counts[h],('inverse fibre',h))
        check(fixed[h]==int(h in shells),('fixed fibre',h))
        check(counts[h]%2==fixed[h],('pointed parity',h))
    check(counts[5]==2 and fixed[5]==0,'pointed hostile H5')
    # Extremizer count is an independent SCC-word binomial count, n=t+6.
    source=load_source('extremizer_tail',ROOT/'fractran/extremizer_tail_source.py')
    p,_=load_program(ROOT/'fractran','extremizer_tail')
    counter=[]
    for t in range(25):
        N=t+6;b=N%3;a=(N-4*b)//3
        expected=sum(seq.count(1)==b for seq in itertools.product((0,1),repeat=a+b))
        check(source(t)==expected,('extremizer source',t))
        actual=ir_run(p,t)
        check(actual[0][p['output']]==expected and actual[0][p['error']]==0,('extremizer IR',t))
        if t<=8:check(raw_run(p,t)[:3]==actual,('extremizer literal',t))
        counter.append({'N':N,'count':expected,'IR_steps':actual[2]})
    path=ROOT/'path_counter/square_shell_source.py'
    spec=importlib.util.spec_from_file_location('square_shell_source',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    p=json.loads((ROOT/'path_counter/square_shell.ir.json').read_text())
    square_cases=[]
    for q in range(4097):
        expected=0 if q==0 else q*q+(q%2)+2*f(q)-1
        check(mod.square_shell(q)==expected,('square source',q))
        if q<=31:
            actual=ir_run(p,q)
            check(actual[0][p['output']]==expected and actual[0][p['error']]==0,('square IR',q))
            if q<=10:
                emitted=emitted_run(p,ROOT/'path_counter/square_shell.fractions.00.txt',q)
                check(emitted==actual,('square emitted',q))
            if q<=5:check(raw_run(p,q)[:3]==actual,('square literal',q))
            square_cases.append({'q':q,'G':expected,'IR_steps':actual[2]})
    census=json.loads((ROOT/'path_counter/FIBRES.json').read_text())
    bound=census['bound_H'];counts_full=[0]*(bound+1)
    # Independent word-state tree, not the scalar f recurrence used by producer.
    stack=[(1,0,0,1)]
    while stack:
        q,c,t,last=stack.pop()
        if q%2==0:
            h=1+2*c
            if h<=bound:counts_full[h]+=1
        for bit in (0,1):
            child=2*q+bit
            if child<=2*(bound-1):
                stack.append((child,c+(t+last if bit==0 else 0),2*t+last,bit))
    digest=hashlib.sha256(','.join(map(str,counts_full)).encode()).hexdigest()
    check(digest==census['counts_sha256'],'full bounded fibre digest')
    check(sum(counts_full)==census['pointed_strong_classes_total'],'full fibre total')
    print(json.dumps({'status':'PASS','cycle_recurrence_through_q':4096,
          'inverse_fibres_all_odd_h_through':2047,'square_shells_through_q':4096,
          'endpoint_correction':'R_1=0; inverse formula must restrict n>=1 or h>=3.',
          'marked_hostile':'H5 has two pointed classes100/110 exchanged by converse; the unique unpointed strong4 class is self-converse.',
          'extremizer_counter_cases':counter,'extremizer_literal_inputs':'0..8',
          'square_counter_cases':square_cases,'square_emitted_inputs':'0..10',
          'square_literal_inputs':'0..5','full_fibre_crosscheck_bound':bound,
          'full_fibre_counts_sha256':digest},indent=2))

if __name__=='__main__':main()
