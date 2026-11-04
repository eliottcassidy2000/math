#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations
import hashlib,json,subprocess
import audit as A
P=Path(__file__).resolve().parent

def inflate(q,sizes):
    blocks=[];n=0
    for s in sizes:blocks.append(list(range(n,n+s)));n+=s
    out=[0]*n
    for i,b in enumerate(blocks):
        for j in range(len(b)):
            for k in range(j+1,len(b)):out[b[j]]|=1<<b[k]
        for j,c in enumerate(blocks):
            if q[i]>>j&1:
                for x in b:
                    for y in c:out[x]|=1<<y
    return tuple(out),blocks

def ports(q,a,c):
    n=len(q);out=list(q)+[1<<(n+1),0]
    for i in range(n):
        if a>>i&1:out[i]|=1<<n
        else:out[n]|=1<<i
        if c>>i&1:out[i]|=1<<(n+1)
        else:out[n+1]|=1<<i
    return tuple(out)
def lift(mask,blocks):return sum(1<<v for i,b in enumerate(blocks) if mask>>i&1 for v in b)
def mobius(vals,bs):
    coeff={}
    for a in bs:
        for c in bs:coeff[a,c]=vals[a,c]-sum(v for (u,w),v in coeff.items() if u&a==u and w&c==w)
    return coeff

def main():
    cases=[((2,0),(2,3)),((2,4,1),(2,2,2)),((2,0),(4,4))];req=[];jobs=[]
    for q,sz in cases:
        T,blocks=inflate(q,sz);U,_=inflate(A.opposite(q),sz);bs=A.basis(len(T));start=len(req)
        for graph in [T,U]:
            for a in bs:
                for c in bs:req.append((ports(graph,a,c),1))
        jobs.append((q,sz,T,U,blocks,bs,start))
    # A truly nonzero actual dyadic-profile control: TT4 children satisfy g modulo16.
    q=(2,4,1);sz=(4,4,4);T,b=inflate(q,sz);U,_=inflate(A.opposite(q),sz);control_T=T;control_U=U;cube=[(a,c) for a in [0,1,2,3] for c in [0,1]];ctrlstart=len(req)
    for graph in [T,U]:
        for a,c in cube:req.append((ports(graph,lift(a,b),lift(c,b)),5))
    # Cyclic ambient root->u->v->root; both ports imitate root, same physical core.
    ambstart=len(req)
    for graph in [T,U]:
        n=len(graph);amb=list(graph)+[1<<(n+1),(1<<n)-1]
        for i in range(n):amb[i]|=1<<n
        for a,c in cube:req.append((ports(amb,lift(a,b)|(1<<(n+1)),lift(c,b)|(1<<(n+1))),3))
    text='\n'.join(' '.join(map(str,[len(g),h,*g])) for g,h in req)+'\n';(P/'physical_input.txt').write_text(text)
    subprocess.run(['g++','-O2','-std=c++17',str(P/'physical.cpp'),'-o',str(P/'physical')],check=True)
    output=subprocess.run([str(P/'physical')],input=text,text=True,capture_output=True,check=True).stdout
    (P/'physical_output.txt').write_text(output);v=[list(map(int,s.split())) for s in output.splitlines()];A.require(len(v)==len(req),'physical outputs');banks=[];coord=0;shared=0
    for q,sz,T,U,blocks,bs,start in jobs:
        n=len(T);case={'quotient':q,'sizes':sz,'source_rows':T,'target_rows':U,'basis':bs,'orientations':[]};w=len(bs)**2
        for gi,g in enumerate([T,U]):
            vals={(a,c):v[start+gi*w+i*len(bs)+j][0] for i,a in enumerate(bs) for j,c in enumerate(bs)};cf=mobius(vals,bs);coord+=len(cf)
            # All queried values, all zero entries, cubic and quartic positions retained.
            case['orientations'].append({'values':[vals[a,c] for a in bs for c in bs],'coefficients':[cf[a,c] for a in bs for c in bs]})
            def ev(a,c):return sum(z for (u,w),z in cf.items() if u&a==u and w&c==w)
            A.require(ev(0,0)==vals[0,0],'constant named bank')
        banks.append(case)
    e,o,_,_,_=A.banks(q=(2,4,1),check_inverse=False);ec,oc,_,_,_=A.banks(A.opposite((2,4,1)),False);d={k:e[k]+o[k]-ec[k]-oc[k] for k in e};nonzero=[];controls=[]
    for j,(a,c) in enumerate(cube):
        source=v[ctrlstart+j];target=v[ctrlstart+8+j];delta=[x-y for x,y in zip(source,target)];expected=A.mod(A.evalbank(d,a,c),16)
        A.require([x%16 for x in delta]==[expected,2*expected%16,0,0,0],'literal dyadic profile vector')
        sd=v[ambstart+j];td=v[ambstart+8+j];ad=[x-y for x,y in zip(sd,td)];A.require(ad[0]%16==5*expected%16,'cyclic ambient odd multiplier5')
        if expected:nonzero.append([a,c,expected])
        controls.append({'query':[a,c],'source_profile':source,'target_profile':target,'source_ambient':sd,'target_ambient':td})
    A.require(nonzero,'literal actual nonzero');(P/'physical_banks.json').write_text(json.dumps(banks,separators=(',',':'))+'\n');(P/'literal_controls.json').write_text(json.dumps({'rows':control_T,'target_rows':control_U,'children':b,'modulus':16,'cube':controls,'nonzero':nonzero},sort_keys=True,indent=2)+'\n')
    receipt={'status':'PASS','literal_full_named_banks':2*len(jobs),'named_coefficients':coord,'physical_queries':len(req),'actual_dyadic_profile_queries':16,'actual_fixed_ambient_queries':16,'nonzero_corners':nonzero,'maximum_augmented_order':16,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['audit.py','physical_audit.py','physical.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['physical_banks.json','literal_controls.json','physical_input.txt','physical_output.txt']}}
    (P/'physical_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
