"""Verify finite bases and the all-order certificate table for six end atoms."""
import json
from pathlib import Path
from inventory import presentations,join,arc,blue,move,need
from five_atom_lift import inflate
from five_atom_probe import atoms
HERE=Path(__file__).resolve().parent

def replay(a,row,p=None):
    m=len(a);p=tuple(range(m)) if p is None else p
    state=join(a,row['atom']);flags=[];states=[];words=[]
    for raw in row['paths']:
        q=tuple(p[v] if v<m else v for v in inflate(raw,m))
        flags.append(blue(state,q));states.append(state);words.append(q)
        state=move(state,q)
    target=join(a,tuple(row['target'][i]>>2 for i in range(2,7)))
    need(state==target,'target or context changed')
    return flags,states,words

def witnesses(i,m,step):
    if i in (0,1):return (0,m-5,m+4,m-6)
    if i==2 and step==1:return (0,m+3,m+2,m)
    if i==3 and step==1:return (0,m+1,m+4,m+3)
    if i==4:return (0,m-1,m+4,m-2)
    return (0,m-3,m+4,m-4)

def main():
    rows=json.loads(HERE.joinpath('five_atom_probe.json').read_text())
    need({tuple(r['atom']) for r in rows}==set(atoms()),'six-atom coverage')
    bases=[];bad=[]
    for i,row in enumerate(rows):
        count=0
        for m in range(1,7):
            for a in presentations(m):
                count+=1;flags,_,_=replay(a,row)
                if any(flags):
                    need(i==0 and a==(30,12,24,16,2) and flags==[False,True,False,False],
                         'new small hostile')
                    p=(0,2,3,4,1)
                    repaired,_,_=replay(a,row,p)
                    need(not any(repaired),'path repair failed')
                    bad.append(dict(atom=i,A=a,original_flags=flags,repair=p))
        need(count==77,'small base incomplete')
        bases.append(dict(atom=i,small_contexts=count))
    checks=0
    for i,row in enumerate(rows):
        for m in range(7,41):
            # TT_m suffices for validating the index/edge identities; the proof
            # uses only source arcs and untouched cross arcs, not its other bits.
            a=tuple(sum(1<<j for j in range(v+1,m)) for v in range(m))
            flags,states,words=replay(a,row)
            need(not any(flags),'large template blue')
            for step in (1,2):
                u,v,r,s=witnesses(i,m,step)
                p=words[step];rho=dict(zip(p,p[::-1]));t=states[step]
                need((rho[v],rho[u])==(r,s),'parametric index identity')
                need(arc(t,u,v) and not arc(t,r,s),'known-arc contradiction')
                checks+=1
    data=dict(status='PASS',finite_bases=bases,sole_repaired_hostile=bad,
              parametric_index_checks=checks,
              scope='six atoms; 77 small context presentations each; proof table m>=7')
    HERE.joinpath('five_atoms_verify.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
