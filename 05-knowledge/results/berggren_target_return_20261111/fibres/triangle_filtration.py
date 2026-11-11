"""Literal controls for an all-k triangle-chain divisibility theorem.

Random probes are seeded finite searches only; no nonexistence is inferred.
The C++ routine is an independent literal subset path counter from this packet.
"""
from pathlib import Path
from math import gcd
import json, random, subprocess, tempfile
from extension_content import C3,S4,tt,join,extend,h,need

HERE=Path(__file__).resolve().parent

def chain(r,gaps):
    blocks=[tt(gaps[0])]
    for j in range(r): blocks.extend((C3,tt(gaps[j+1])))
    return join(blocks)

def valuation(n):
    v=0
    while n%3==0: v+=1;n//=3
    return v

def main():
    hostile=extend(join((tt(1),C3,tt(1))),1)
    path=(1,4,5,0,2,3)
    orbit=[]
    for shift in range(3):
        rotated=tuple(1+(v-1+shift)%3 if 1<=v<=3 else v for v in path)
        need(all(hostile[a]>>b&1 for a,b in zip(rotated,rotated[1:])),
             'untouched split triangle rotation')
        orbit.append(rotated)
    need(len(set(orbit))==3,'free triangle action')
    pair=extend(join((C3,C3)),9)
    need(h(pair)==41,'one birth on two adjacent triangles')
    explicit=[]
    for r in range(1,6):
        for k in range((r+1)//2+1):
            pairs=min(k,r//2); singles=min(k-pairs,r-2*pairs)
            left=r-2*pairs-singles
            witness=join([pair]*pairs+[S4]*singles+[C3]*left)
            expected=41**pairs*5**singles*3**left
            if len(witness)<=15: need(h(witness)==expected,'literal lower-content witness')
            need(gcd(3**r,expected)==3**max(r-2*k,0),'exact filtration witness gcd')
            explicit.append({'r':r,'k':k,'H':expected,'rows':witness,
                             'content':3**max(r-2*k,0)})
    rng=random.Random(20261111)
    cases=[]; graph_inputs=[]
    # Exhaustive first-step structure is in extension_content.py. Here probe
    # two births where all middle gaps are positive and that formula is silent.
    for r,gaps in [(3,(0,1,1,0)),(3,(0,2,1,0)),(4,(0,1,1,1,0))]:
        source=chain(r,gaps); n=len(source)
        for trial in range(40):
            m1=rng.randrange(1<<n); m2=rng.randrange(1<<(n+1))
            target=extend(extend(source,m1),m2)
            cases.append({'r':r,'gaps':gaps,'trial':trial,'masks':[m1,m2],
                          'rows':target})
            graph_inputs.append(str(n+2)+' 0\n'+' '.join(map(str,target))+'\n')
    with tempfile.TemporaryDirectory(prefix='triangle_filtration_') as td:
        executable=Path(td)/'audit.exe'
        subprocess.run(['g++','-O3','-std=c++17',str(HERE.parent/'arithmetic/audit_state.cpp'),
                        '-o',str(executable)],check=True,capture_output=True,text=True)
        process=subprocess.run([str(executable)],input=''.join(graph_inputs),
                               capture_output=True,text=True,check=True)
        values=[int(line.split()[0]) for line in process.stdout.splitlines()]
    need(len(values)==len(cases),'literal batch size')
    for case,value in zip(cases,values):
        case['H']=value; case['v3']=valuation(value)
        need(case['v3']>=max(case['r']-4,0),'untouched triangle lower bound')
    summary=[]
    for r,gaps in [(3,(0,1,1,0)),(3,(0,2,1,0)),(4,(0,1,1,1,0))]:
        selected=[c for c in cases if c['r']==r and tuple(c['gaps'])==gaps]
        g=3**r
        for c in selected:g=gcd(g,c['H'])
        summary.append({'r':r,'gaps':gaps,'sample_gcd_with_source':g,
                        'min_v3':min(c['v3'] for c in selected),
                        'first_min_witness':min(selected,key=lambda c:c['v3'])})
    out={'status':'PROVED family controls plus seeded finite two-birth probes',
         'noncontiguous_untouched_triangle_hostile':{'rows':hostile,'orbit':orbit},
         'pair_rows':pair,'pair_H':41,'explicit_filtration_witnesses':explicit,
         'probe_summary':summary,'probe_count':len(cases),'probe_seed':20261111,
         'all_probe_counts':[c['H'] for c in cases]}
    HERE.joinpath('triangle_filtration.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'result':'PASS','probe_summary':[
        {k:v for k,v in x.items() if k!='first_min_witness'} for x in summary]}))

if __name__=='__main__':main()
