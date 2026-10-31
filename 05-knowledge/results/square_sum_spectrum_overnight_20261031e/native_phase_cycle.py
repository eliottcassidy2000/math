#!/usr/bin/env python3
"""Native protected cycles with2401 phase-fixed tree inputs, decoded from the output itself."""
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
from phase_trees import Trees,threshold,code,prefix,arity,phase
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from independent_ray_channels import P
OLD=HERE.parent/'seeded_eventual_cycle_theorem_20261019'/'certified-pair-transitions'
FINAL=HERE.parent/'square_sum_spectrum_overnight_20261031a'/'n24_forest_table.json'


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def edge(a,b):return tuple(sorted((a,b)))


def inner_data(word,m):
    at=0;ports={}
    for atom in word:
        kind,c=atom[:2]
        if kind=='single':
            if c in (2,7):ports[c]=at
            at+=1
        elif kind=='head':at+=1
        else:at+=m+atom[2]-int(kind=='tail')
    need(abs(ports[2]-ports[7])==1,'fixed A ports')
    return min(ports.values())+1,1 if ports[2]<ports[7] else -1


def decorated_A(trees,pairs,codes,b,inner,n):
    children={c:pairs[codes[c,b]] for c in range(-24,25)}
    out=[trees.expand(w,49,children) for w in inner['words']]
    colors={v:j%2 for j,v in enumerate(out[0])}
    for i,p in enumerate(out):
        need(len(p)==n+i and sorted(p)==list(range(1,n+i+1)), 'signed A ownership')
        need(p[:2]==[69,127] and abs(p.index(2)-p.index(7))==1 and p.index(2)%2==0, 'signed A ports')
        need(all(isqrt(a+z)**2==a+z for a,z in zip(p,p[1:])), 'signed A edges')
        if i:need(all(colors[v]==j%2 for j,v in enumerate(p) if v<=n),'signed A phase')
    return out


def audit(word,extra,flip_zero=False):
    N=2401*threshold(word)+1200+extra
    n=(N-24)//49;m=(n-24)//49;d=N-49*n
    trees=Trees();codes={(a,b):code(word,f'native:{a}:{b}') for a in range(-24,25) for b in range(-24,25)}
    if flip_zero:
        value=codes[0,0];codes[0,0]=('1' if value[0]=='0' else '0')+value[1:]
    pairs={w:trees.pair(m,word,w) for w in sorted(set(codes.values()))}
    family='Aplus' if (n+int(d>=57))%2==0 else 'Aminus'
    inner=json.loads((OLD/f'{family}_p{m%2}_d{n-49*m}.json').read_text())
    data=[inner_data(w,m) for w in inner['words']];cuts=[v[0] for v in data];orders=[v[1] for v in data]
    banks={b:decorated_A(trees,pairs,codes,b,inner,n) for b in range(-24,25)}
    reference=[phase(p) for p in banks[0]]
    need(all([phase(p) for p in bank]==reference for bank in banks.values()),'outer choices share the full source phase')
    final=next(r for r in json.loads(FINAL.read_text()) if (r['parity'],r['residue'],r['order0'],r['order1'])==(n%2,d,*orders))
    output=[]
    for atom in final['word']:
        if atom[0]=='single':output.append(atom[1]);continue
        _,c,i,part,o=atom;p=banks[c][i];lo,hi=(0,cuts[i]) if part==0 else(cuts[i],len(p))
        indices=range(lo,hi) if o==1 else range(hi-1,lo-1,-1)
        output.extend(49*p[j]+(-1)**(j+1)*c for j in indices)
    need(len(output)==N,'native exact length');seen=bytearray(N+1)
    required={edge(a,b) for a,b in zip(P,P[1:])}|{edge(482,3362),edge(3362,6242)}
    ray_prefix=prefix(word)
    for u in range(-24,25):
        if u in (0,2):continue
        for v in range(-24,25):
            t=49*u+v;ray=[2401*x+(-1)**j*t for j,x in enumerate(ray_prefix)]
            required.update(edge(a,b) for a,b in zip(ray,ray[1:]))
    edge_count=len(required)
    for j,v in enumerate(output):
        need(1<=v<=N and not seen[v],'native ownership');seen[v]=1
        u=output[j-1];need(isqrt(u+v)**2==u+v,'native square edge');required.discard(edge(u,v))
    need(not required,'native marked forest')
    # Erase the producer's starting point and direction. Recover its atom clock
    # from the canonical undirected cycle and the fixed native edge 1--8.
    start=output.index(1);step=1 if output[(start+1)%N]==8 else -1
    need(output[(start+step)%N]==8,'native orientation anchor')
    canonical_word=[output[(start+j*step)%N] for j in range(N)]
    clock=0;one_clock=None
    for atom in final['word']:
        if atom==['single',1]:one_clock=clock
        clock+=1 if atom[0]=='single' else (cuts[atom[2]] if atom[3]==0 else n+atom[2]-cuts[atom[2]])
    one_atom=final['word'].index(['single',1])
    direction=1 if final['word'][(one_atom+1)%len(final['word'])]==['single',8] else -1
    recovered_clock=[canonical_word[(direction*(j-one_clock))%N] for j in range(N)]
    need(recovered_clock==output,'canonical cycle restores the integer atom clock')
    output=recovered_clock
    # Consumer: read all signed outer offsets, reassemble split A paths, and
    # independently read all signed inner tails. No original bank is consulted.
    at=0;fragments={};members={}
    for atom in final['word']:
        if atom[0]=='single':
            need(output[at]==atom[1],'fixed final singleton');at+=1;continue
        _,c,i,part,o=atom;length=cuts[i] if part==0 else n+i-cuts[i]
        raw=output[at:at+length]
        if o==-1:raw=raw[::-1]
        fragments[c,part]=[(v+24)//49 for v in raw];members[c]=i
        at+=length
    recovered={}
    for b in range(-24,25):
        path=fragments[b,0]+fragments[b,1];atoms=inner['words'][members[b]];at=0
        for atom in atoms:
            kind,c=atom[:2]
            if kind in ('single','head'):at+=1;continue
            _,c,i,o=atom;length=m+i-int(kind=='tail')
            raw=path[at:at+length]
            if o==-1:raw=raw[::-1]
            child=[(v+24)//49 for v in raw]
            if kind=='tail':child.insert(0,1)
            recovered[c,b]=trees.decode(child,m,word,i)
            at+=length
        need(at==len(path),'reconstructed A length')
    need(recovered==codes,'all2401 source choices decoded from final cycle')
    digest=hashlib.sha256()
    for v in output:digest.update((str(v)+'\n').encode())
    start=output.index(1);step=1 if output[(start+1)%N]==8 else -1
    need(output[(start+step)%N]==8,'canonical native orientation')
    canonical=hashlib.sha256()
    for j in range(N):canonical.update((str(output[(start+j*step)%N])+'\n').encode())
    return {'N':N,'word':word,'leaf_choices':2401*arity(word),'recovered_channels':len(recovered),
            'decoded_from_canonical_undirected_cycle':True,
            'changed_zero_channel':flip_zero,'prescribed_edges':edge_count,'cycle_sha256':digest.hexdigest(),
            'canonical_cycle_sha256':canonical.hexdigest()}


if __name__=='__main__':
    controls=[audit('',0),audit('',1),audit('',0,True)]
    need(controls[0]['canonical_cycle_sha256']!=controls[2]['canonical_cycle_sha256'],'one zero-channel bit changes the physical cycle')
    print(json.dumps({'status':'passed','controls':controls,
                      'scope':'Full native cycles at both boundary parities and a one-zero-channel change with the other2400 choices fixed; general depth uses the complete readable interfaces and the independently audited tree induction.'},indent=2,sort_keys=True))
