#!/usr/bin/env python3
"""Native protected cycles with625 tree inputs, decoded from the output itself."""
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
from tree_pairs import Trees,threshold,code,prefix
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from independent_ray_channels import decorated_A,P
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


def audit(depth,extra,flip_zero=False):
    N=2401*threshold(depth)+1200+extra
    n=(N-24)//49;m=(n-24)//49;d=N-49*n
    trees=Trees();codes={(a,b):code(depth,f'native:{a}:{b}') for a in range(25) for b in range(25)}
    if flip_zero:
        value=codes[0,0];codes[0,0]=('B' if value[0]=='A' else 'A')+value[1:]
    pairs={w:trees.pair(m,depth,w) for w in sorted(set(codes.values()))}
    family='Aplus' if (n+int(d>=57))%2==0 else 'Aminus'
    inner=json.loads((OLD/f'{family}_p{m%2}_d{n-49*m}.json').read_text())
    data=[inner_data(w,m) for w in inner['words']];cuts=[v[0] for v in data];orders=[v[1] for v in data]
    banks=[decorated_A(pairs,codes,b,family,n) for b in range(25)]
    final=next(r for r in json.loads(FINAL.read_text()) if (r['parity'],r['residue'],r['order0'],r['order1'])==(n%2,d,*orders))
    output=[]
    for atom in final['word']:
        if atom[0]=='single':output.append(atom[1]);continue
        _,c,i,part,o=atom;p=banks[abs(c)][i];lo,hi=(0,cuts[i]) if part==0 else(cuts[i],len(p))
        indices=range(lo,hi) if o==1 else range(hi-1,lo-1,-1)
        output.extend(49*p[j]+(-1)**(j+1)*c for j in indices)
    need(len(output)==N,'native exact length');seen=bytearray(N+1)
    required={edge(a,b) for a,b in zip(P,P[1:])}|{edge(482,3362),edge(3362,6242)}
    ray_prefix=prefix('C'*depth)
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
    # Consumer: read positive outer offsets, reassemble split A paths, and
    # independently read positive inner tails. No original bank is consulted.
    at=0;fragments={};members={}
    for atom in final['word']:
        if atom[0]=='single':
            need(output[at]==atom[1],'fixed final singleton');at+=1;continue
        _,c,i,part,o=atom;length=cuts[i] if part==0 else n+i-cuts[i]
        if c>=0:
            raw=output[at:at+length]
            if o==-1:raw=raw[::-1]
            fragments[c,part]=[(v+24)//49 for v in raw];members[c]=i
        at+=length
    recovered={}
    for b in range(25):
        path=fragments[b,0]+fragments[b,1];word=inner['words'][members[b]];at=0
        for atom in word:
            kind,c=atom[:2]
            if kind in ('single','head'):at+=1;continue
            _,c,i,o=atom;length=m+i-int(kind=='tail')
            if c>=0:
                raw=path[at:at+length]
                if o==-1:raw=raw[::-1]
                child=[(v+24)//49 for v in raw]
                if kind=='tail':child.insert(0,1)
                recovered[c,b]=trees.decode(child,m,depth,i)
            at+=length
        need(at==len(path),'reconstructed A length')
    need(recovered==codes,'all625 source choices decoded from final cycle')
    digest=hashlib.sha256()
    for v in output:digest.update((str(v)+'\n').encode())
    start=output.index(1);step=1 if output[(start+1)%N]==8 else -1
    need(output[(start+step)%N]==8,'canonical native orientation')
    canonical=hashlib.sha256()
    for j in range(N):canonical.update((str(output[(start+j*step)%N])+'\n').encode())
    return {'N':N,'depth':depth,'leaf_choices':625*25**depth,'recovered_channels':len(recovered),
            'changed_zero_channel':flip_zero,'prescribed_edges':edge_count,'cycle_sha256':digest.hexdigest(),
            'canonical_cycle_sha256':canonical.hexdigest()}


if __name__=='__main__':
    controls=[audit(0,0),audit(0,1),audit(0,0,True)]
    need(controls[0]['canonical_cycle_sha256']!=controls[2]['canonical_cycle_sha256'],'one zero-channel bit changes the physical cycle')
    print(json.dumps({'status':'passed','controls':controls,
                      'scope':'Full native cycles at both boundary parities and a one-zero-channel change with the other624 choices fixed; general depth uses the complete readable interfaces and the independently audited tree induction.'},indent=2,sort_keys=True))
