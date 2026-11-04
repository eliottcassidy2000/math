"""Cross-check independent C++ DP against all inherited marked endpoint fields."""
from pathlib import Path
import json, subprocess
P=Path(__file__).resolve().parent

def run(rows,indices):
    n=len(rows); inp=' '.join(map(str,[n,len(indices),*rows,*indices]))
    out=subprocess.run([str(P/'joint')],input=inp,text=True,capture_output=True,check=True).stdout
    h=None;fields={};joints={}
    for line in out.splitlines():
        a=line.split()
        if a[0]=='H':h=int(a[1])
        elif a[0]=='PORT':fields[int(a[1])]=list(map(int,a[2:]))
        elif a[0]=='J':joints[(int(a[1]),int(a[2]))]=int(a[3])
    return h,fields,joints

def main():
    count=0
    for z in json.loads((P/'input_audit.json').read_text())['states']:
        ns=z['names'];inds=[ns.index(x) for x in [50,6,51]]
        h,f,j=run(z['rows'],inds)
        if h!=z['H']:raise RuntimeError('H')
        for pair in z['joint_pairs']:
            for name,v in zip(pair['ports'],pair['port_data']):
                if f[ns.index(name)]!=v:raise RuntimeError((z['stage'],name,f[ns.index(name)],v))
                count+=1
    z=json.loads((P.parent/'returns'/'candidate.json').read_text());ns=z['names'];names=[39,6,51,4,50]
    h,f,j=run(z['rows'],[ns.index(x) for x in names])
    if h!=3*z['source_H']+2:raise RuntimeError('new return')
    out={'H':h,'root':39,'fields_order':['E1','S1','deletion_H','E2','S2'],'fields':{str(x):f[ns.index(x)] for x in names},'joints':[{'end':ns[a],'start':ns[b],'value':v} for (a,b),v in j.items()],'inherited_field_comparisons':count}
    if out['fields']['39'][3]!=h:raise RuntimeError('zero ready')
    (P/'candidate_joint.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS',count,'inherited port fields and candidate named joint fields')
if __name__=='__main__':main()
