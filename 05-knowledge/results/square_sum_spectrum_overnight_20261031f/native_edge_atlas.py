#!/usr/bin/env python3
"""Transport edge switches through the native two-layer, cut-and-rejoin carrier."""
import json,sys
from pathlib import Path
from edge_atlas import Atlas,Carrier,EdgeCube,edge,edges,need,arity,threshold,code,digest
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031e'))
from native_phase_cycle import OLD,FINAL,inner_data,P

def build(N,word):
    a=Atlas();n=(N-24)//49;m=(n-24)//49;d=N-49*n;chunk=arity(word)
    need(N>=2401*threshold(word)+1200,'protected atlas threshold')
    family='Aplus' if (n+int(d>=57))%2==0 else 'Aminus'
    inner=json.loads((OLD/f'{family}_p{m%2}_d{n-49*m}.json').read_text())
    data=[inner_data(w,m) for w in inner['words']];cuts=[v[0] for v in data];orders=[v[1] for v in data]
    banks={};children={}
    for b in range(-24,25):
        children[b]={c:a.pair(m,word,((b+24)*49+c+24)*chunk) for c in range(-24,25)}
        banks[b]=[a.expand(w,49,children[b]) for w in inner['words']]
    record=next(r for r in json.loads(FINAL.read_text()) if (r['parity'],r['residue'],r['order0'],r['order1'])==(n%2,d,*orders))
    out=[];moves={};used=set()
    for atom in record['word']:
        if atom[0]=='single':out.append(atom[1]);continue
        _,c,i,part,o=atom;source=banks[c][i];p=source.path
        lo,hi=(0,cuts[i]) if part==0 else(cuts[i],len(p));indices=range(lo,hi) if o==1 else range(hi-1,lo-1,-1)
        out.extend(49*p[j]+(-1)**(j+1)*c for j in indices)
        if c in used:continue
        used.add(c);pos={v:j for j,v in enumerate(p)}
        def lift(v):return 49*v+(-1)**(pos[v]+1)*c
        for b,(removed,added) in source.moves.items():
            need(all((pos[u]<cuts[i])==(pos[v]<cuts[i]) for u,v in removed|added),'no varying edge crosses native cut')
            moves[b]=tuple({edge(lift(u),lift(v)) for u,v in part} for part in (removed,added))
    cube=EdgeCube(Carrier(out,moves),anchor=8)
    return a,cube,inner,record,cuts,m,n

def literal(a,N,word,bits,inner,record,cuts,m,n):
    chunk=arity(word);banks={}
    for b in range(-24,25):
        children={}
        for c in range(-24,25):
            start=((b+24)*49+c+24)*chunk
            children[c]=a.t.pair(m,word,bits[start:start+chunk])
        banks[b]=[a.t.expand(w,49,children) for w in inner['words']]
    out=[]
    for atom in record['word']:
        if atom[0]=='single':out.append(atom[1]);continue
        _,c,i,part,o=atom;p=banks[c][i];lo,hi=(0,cuts[i]) if part==0 else(cuts[i],len(p))
        indices=range(lo,hi) if o==1 else range(hi-1,lo-1,-1)
        out.extend(49*p[j]+(-1)**(j+1)*c for j in indices)
    return out

def main():
    records=[]
    for extra in (0,1):
        N=270112+extra;word='';a,cube,inner,r,cuts,m,n=build(N,word)
        import hashlib
        bits=''.join(str(hashlib.sha256(f'native:{extra}:{b}'.encode()).digest()[0]%2) for b in range(cube.bits))
        from_edges=cube.walk(bits);actual=literal(a,N,word,bits,inner,r,cuts,m,n)
        need(edges(from_edges,True)==edges(actual,True),'native literal versus edge-union consumer')
        protected={edge(u,v) for u,v in zip(P,P[1:])}|{edge(482,3362),edge(3362,6242)}
        need(protected<=cube.base and all(e not in cube.controls for e in protected),'native ports constant')
        ports=[min(v for e in cube.moves[b][0]|cube.moves[b][1] for v in e) for b in [0,49,1200,2400]]
        queries=[(e,True) for e in cube.base if any(v in e for v in ports)]
        compiled=cube.compile(queries);need(compiled is not None and len(compiled['pins'])==4,'four arbitrary changed ports cost four bits')
        records.append({'N':N,'bits':cube.bits,'edge_controls':len(cube.controls),'changed_vertex_support':len(cube.owners),
                        'pinned_vertices':ports,'paid_bits':4,'surviving_bits':compiled['free_bits'],'cycle_sha256':digest(from_edges)})
    print(json.dumps({'status':'passed','records':records,'scope':'Native path and seed remain constant; exact query compilation applies to the full edge-cube family.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
