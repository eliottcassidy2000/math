#!/usr/bin/env python3
"""Recount every source-bank value by independent predecessor subset DP."""
from pathlib import Path
import hashlib,json,subprocess,tempfile,sys
P=Path(__file__).resolve().parent

def need(x,m):
    if not x:raise RuntimeError(m)
def ports(old,a,c):
    n=len(old);out=list(old)+[1<<(n+1),0]
    for i in range(n):
        if a>>i&1:out[i]|=1<<n
        else:out[n]|=1<<i
        if c>>i&1:out[i]|=1<<(n+1)
        else:out[n+1]|=1<<i
    return out
def subs(mask):
    s=mask
    while True:
        yield s
        if not s:break
        s=(s-1)&mask

def main():
    odd=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent/'session189-odd-transitive';anchor=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else P.parent/'session189-collision-next'
    alljobs=[];sourcepins={};artifacts={};total=0;degrees=[0]*5
    for label,src,files in [('odd',odd,[f'bank.part{i}.json' for i in range(5)]),('anchor',anchor,['bank.json'])]:
        manifest=json.loads((src/'RELEASE_FILES.json').read_text());sourcepins[label]={}
        for f,meta in manifest.items():
            h=hashlib.sha256((src/f).read_bytes()).hexdigest();need(h==meta['sha256'],'frozen source '+f);sourcepins[label][f]=h
        for f in files:
            graphs=json.loads((src/f).read_text())
            for B in graphs:
                n=len(B['rows']);need(all(not(row>>i&1) and all(((row>>j)&1)+((B['rows'][j]>>i)&1)==1 for j in range(n) if j!=i) for i,row in enumerate(B['rows'])),'valid named graph')
                want=[0]+[1<<i for i in range(n)]+[s for s in range(1<<n) if s.bit_count()==2];need(B['masks']==want,'complete named basis');need(len(B['values'])==len(want)**2==len(B['coefficients']),'complete bank lengths')
            alljobs.append((label,src,f,graphs))
    subprocess.run(['g++','-O3','-std=c++17',str(P/'count.cpp'),'-o',str(P/'count')],check=True)
    for label,src,f,graphs in alljobs:
        lines=[]
        for B in graphs:
            for a in B['masks']:
                for c in B['masks']:
                    q=ports(B['rows'],a,c);lines.append(' '.join(map(str,[len(q),*q])))
        out=subprocess.check_output([str(P/'count')],input='\n'.join(lines)+'\n',text=True);numbers=list(map(int,out.split()));need(len(numbers)==len(lines),'complete recount output');offset=0;rebuilt=[]
        for B in graphs:
            bs=B['masks'];w=len(bs);values=numbers[offset:offset+w*w];offset+=w*w;need(values==B['values'],'independent actual named basis values');index={a:i for i,a in enumerate(bs)};coeff=[]
            for a in bs:
                for c in bs:
                    v=sum((-1)**((a^s).bit_count()+(c^t).bit_count())*values[index[s]*w+index[t]] for s in subs(a) for t in subs(c));coeff.append(v);degrees[a.bit_count()+c.bit_count()]+=1
            need(coeff==B['coefficients'],'all named coefficients');total+=len(coeff);R=dict(B);R['values']=values;R['coefficients']=coeff;rebuilt.append(R)
        name=label+'-recount.'+f;(P/name).write_text(json.dumps(rebuilt,separators=(',',':'))+'\n');artifacts[name]=hashlib.sha256((P/name).read_bytes()).hexdigest()
    receipt={'status':'PASS','independently_recounted_graphs':sum(len(g) for _,_,_,g in alljobs),'independently_recounted_named_coefficients':total,'counts_by_total_degree':degrees,'source_pins':sourcepins,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['bank_audit.py','count.cpp']},'artifacts':artifacts}
    (P/'bank_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
