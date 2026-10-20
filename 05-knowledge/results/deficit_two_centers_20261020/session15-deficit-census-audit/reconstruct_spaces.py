#!/usr/bin/env python3
"""Independent affine GF(2) reconstruction. Imports no producer code."""
import gzip, hashlib, json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'session15-deficit-centers'
FAMILIES=[(10,2,None),(12,4,None),(14,2,None),(15,3,None),(20,4,None),(21,3,0),(21,3,1),(21,3,2),(28,4,0),(28,4,1),(28,4,2),(28,4,3)]

def require(t, message):
    if not t: raise RuntimeError(message)

def affine(n,q,low):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    idx={e:i for i,e in enumerate(pairs)}; L=len(pairs); rhs=1<<L
    equations=[]
    def equation(a,b,c,d,anti):
        u,v=sorted((a,b)); x,y=sorted((c,d))
        equations.append((1<<idx[u,v])^(1<<idx[x,y])^(rhs if ((a>b)^(c>d)^anti) else 0))
    for i,j in pairs:
        equation(i,j,(i+q)%n,(j+q)%n,False)
        if q in (3,4): equation(i,j,n-1-i,n-1-j,True)
    def fix(i,j,value):
        a,b=sorted((i,j))
        equations.append((1<<idx[a,b])^(rhs if value^(i>j) else 0))
    for i in range(n):fix(i,(i+1)%n,1)
    if low is not None:
        for j in range(n):
            if j%q!=low:fix(low,j,int(j==(low+1)%n))
    pivots={}
    for row in equations:
        while row and row!=rhs:
            p=(row&-row).bit_length()-1
            if p not in pivots:pivots[p]=row;break
            row ^=pivots[p]
        if row==rhs:return None
    free=[j for j in range(L) if j not in pivots]
    def solve(value,homogeneous=False):
        for p in sorted(pivots,reverse=True):
            row=pivots[p]
            bit=((row&value).bit_count()&1)^((row>>L)&1 if not homogeneous else 0)
            if bit:value^=1<<p
        return value
    base=solve(0)
    directions=[solve(1<<j,True) for j in free]
    for row in equations:
        require(((row&base).bit_count()&1)==((row>>L)&1),'base violates equation')
        require(all(not ((row&v).bit_count()&1) for v in directions),'basis violates equation')
    return pairs,equations,base,directions,free

def rows(n,pairs,bits,as_toggle=False):
    out=[0]*n
    for e,(i,j) in enumerate(pairs):
        b=(bits>>e)&1
        if as_toggle:
            if b:out[i]|=1<<j;out[j]|=1<<i
        elif b:out[i]|=1<<j
        else:out[j]|=1<<i
    return out

def file_model(path,n,q,low,pairs,base,dirs):
    with path.open('w') as f:
        f.write(f'{n} {q} {-1 if low is None else low} {len(dirs)}\n')
        f.write(' '.join(map(str,rows(n,pairs,base)))+'\n')
        for d in dirs:f.write(' '.join(map(str,rows(n,pairs,d,True)))+'\n')

def rank(vectors):
    piv={}
    for v in vectors:
        while v:
            p=v.bit_length()-1
            if p not in piv:piv[p]=v;break
            v^=piv[p]
    return len(piv)

def producer_model(spec,pairs,equations,independent_f):
    n=spec['n'];L=len(pairs);fixed={int(k):v for k,v in spec['fixed'].items()}
    free=spec['free'];bits=spec['bits'];K=len(bits)
    require(spec['orbits']==K,'wrong K')
    require(set(fixed).isdisjoint(free) and set(fixed)|set(free)==set(range(K)),'bad partition')
    require(free==sorted(free),'unordered free basis')
    owned={};choices=[[],[]]
    for k,(b0,b1) in enumerate(bits):
        require(len(b0)==n==len(b1),'wrong row size')
        require(all(not (b[i]>>i&1) and not b[i]>>n for b in (b0,b1) for i in range(n)),'invalid arcs')
        for e,(i,j) in enumerate(pairs):
            arc0=(b0[i]>>j&1,b0[j]>>i&1)
            arc1=(b1[i]>>j&1,b1[j]>>i&1)
            if arc0==(0,0) and arc1==(0,0):continue
            require(arc0 in [(0,1),(1,0)] and arc1==arc0[::-1], 'bad complementary pair')
            require(e not in owned,'edge owned twice')
            owned[e]=(k,arc0[0])
    require(len(owned)==L,'missing pairs')
    base=0;dirs=[0]*len(free);fmap={k:t for t,k in enumerate(free)}
    for e,(k,b0) in owned.items():
        base|=(b0^fixed.get(k,0))<<e
        if k in fmap:dirs[fmap[k]]|=1<<e
    for row in equations:
        require(((row&base).bit_count()&1)==((row>>L)&1),'producer base violates full linear equation')
        require(all(not ((row&v).bit_count()&1) for v in dirs),'producer basis violates full linear equation')
    require(rank(dirs)==len(dirs)==independent_f, 'producer affine space not equal to independent space')
    # Both affine spaces now have the same dimension; inclusion proves equality.
    return base,dirs

def main():
    specs=json.loads((SOURCE/'residual_frame_specs.json').read_text())
    specmap={(s['n'],s['q'],s['low_out_orbit']):s for s in specs}
    require(len(specmap)==len(specs),'duplicate specification keys')
    index=json.loads((SOURCE/'witness_certificate_index.json').read_text())
    idx={Path(x['file']).name:x for x in index};require(len(idx)==len(index),'duplicate certificate')
    report=[]
    for n,q,low in FAMILIES:
        name=f'n{n}_q{q}_low{low}'
        model=affine(n,q,low)
        if model is None:
            require((n,q,low) not in specmap,'inconsistent family not rejected')
            report.append(dict(n=n,q=q,low=low,consistent=False));continue
        pairs,eqs,base,dirs,free=model
        spec=specmap[n,q,low]
        require(spec['anti']==(q in [3,4]), 'q2 reflection imposed')
        pbase,pdirs=producer_model(spec,pairs,eqs,len(dirs))
        file_model(ROOT/f'independent_{name}.model',n,q,low,pairs,base,dirs)
        file_model(ROOT/f'producer_{name}.model',n,q,low,pairs,pbase,pdirs)
        cert=f'frames_{name}.txt.witnesses.bin'
        raw=(SOURCE/cert).read_bytes();gz=(SOURCE/(cert+'.gz')).read_bytes()
        entry=idx[cert+'.gz']
        require(gzip.decompress(gz)==raw,'decompression mismatch')
        require(len(raw)==5*(1<<len(dirs)),'certificate length mismatch')
        require(entry['records']==1<<len(dirs),'record count mismatch')
        require(entry['raw_sha256']==hashlib.sha256(raw).hexdigest(),'raw hash mismatch')
        require(entry['gzip_sha256']==hashlib.sha256(gz).hexdigest(),'gzip hash mismatch')
        require(entry['gzip_bytes']==len(gz),'gzip size mismatch')
        report.append(dict(n=n,q=q,low=low,consistent=True,variables=len(pairs),rank=len(pairs)-len(dirs),free=len(dirs),presentations=1<<len(dirs),producer_space_equal=True,certificate_hashes_verified=True,certificate_raw_sha256=entry['raw_sha256']))
    require(len(specs)==sum(x['consistent'] for x in report),'unexpected producer family')
    total=sum(x.get('presentations',0) for x in report)
    require(total==2121984,'wrong total')
    (ROOT/'space_reconstruction.json').write_text(json.dumps(dict(total=total,families=report),indent=2)+'\n')
    print(json.dumps(dict(total=total,families=report),indent=2))
if __name__=='__main__':main()
