#!/usr/bin/env python3
"""Independent whole-integer checks of priority-aware repeated-word certificates."""
import importlib.util,itertools,json,math
from pathlib import Path
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]/'acceleration'
spec=importlib.util.spec_from_file_location('producer_acceleration',ROOT/'accelerate.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
def check(ok,msg):
    if not ok:raise RuntimeError(msg)
def raw(n,rows):
    for i,f in enumerate(rows):
        if n*f.numerator%f.denominator==0:return n*f.numerator//f.denominator,i
    return n,None
def main():
    guards=0
    for base in itertools.product(range(-3,4),repeat=2):
        for direction in itertools.product(range(-2,3),repeat=2):
            for need in itertools.product(range(3),repeat=2):
                interval=a.eligible_interval(base,direction,need)
                for t in range(31):
                    yes=all(x+t*d>=b for x,d,b in zip(base,direction,need))
                    got=interval is not None and interval[0]<=t and (interval[1] is None or t<=interval[1])
                    check(yes==got,('guard',base,direction,need,t,interval));guards+=1
    # Every two-row program whose 2/3 valuations lie in {-1,0,1}^2.
    fractions=[]
    for d,e in itertools.product(range(-1,2),repeat=2):
        fractions.append((2**max(d,0)*3**max(e,0),2**max(-d,0)*3**max(-e,0)))
    words=finite=infinite=0
    for rows in itertools.product(fractions,repeat=2):
        p=a.Program(rows)
        for aa,bb in itertools.product(range(5),repeat=2):
            initial=2**aa*3**bb;x,inert=p.encode(initial);n=initial;word=[]
            for length in range(1,5):
                n,index=raw(n,p.rows)
                if index is None:break
                word.append(index);c=a.certify_word(p,x,word)
                R=c['repetitions'];v=initial
                rounds=40 if R is None else R
                for repeat in range(rounds):
                    for wanted in word:
                        v,got=raw(v,p.rows);check(got==wanted,('certified word failed',rows,initial,word,R,repeat))
                if R is None:infinite+=1
                else:
                    check(v==p.decode(c['after'],inert),('endpoint',rows,initial,word))
                    check(c['primitive_steps']==R*len(word),'cost')
                    broke=False
                    for wanted in word:
                        v,got=raw(v,p.rows)
                        if got!=wanted:broke=True;break
                    check(broke,('not maximal',rows,initial,word,R));finite+=1
                words+=1
    priority=a.Program([(1,10),(5,2)]);x,_=priority.encode(8);r=a.run(priority,x)
    check(r['status']=='halted' and priority.decode(r['state'])==5 and r['primitive_steps']==3,'priority hostile')
    # Giant acceleration is verified by exact affine/cost arithmetic, not 10^18 replay.
    p=a.Program([(3,2),(1,3)]);x,_=p.encode(2**100)
    # One 3/2 then1/3 would not be executable because3/2 has priority while2 remains;
    # it is deliberately rejected rather than accidentally treated as a drain loop.
    try:a.certify_word(p,x,[0,1])
    except ValueError:pass
    else:raise RuntimeError('unprioritized paired word accepted')
    print(json.dumps({'status':'PASS','guard_membership_cases':guards,'two_row_programs':81,'input_vectors_per_program':25,'certified_words':words,'finite_maximality_cases':finite,'infinite_prefix_controls_40_rounds':infinite,'priority_hostile_terminal':5},indent=2))
if __name__=='__main__':main()
