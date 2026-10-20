#!/usr/bin/env python3
"""Exact SC class census and odd-orbit toggle audit, without HP/classifier imports.

Two independent presentations are used:
(1) normalized anti-involution frames (k^2 bits) modulo S_k;
(2) all reflection-fixed Hamiltonian-path tilings (k^2 bits).
No full order-(2k+1) tournament or full fixed-path cube is enumerated.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
import argparse, json, math, time


def run(k, outdir, inherited_dir=None):
    n=2*k+1
    pairs=[(i,j) for i in range(k) for j in range(i+1,k)]
    ne=len(pairs); pair_index={ij:z for z,ij in enumerate(pairs)}
    perms=list(permutations(range(k))); actions=[]
    for p in perms:
        dst=list(p); flip=0
        for offset,directed in [(k,True),(k+ne,False)]:
            for z,(i,j) in enumerate(pairs):
                u,v=p[i],p[j]
                q=offset+pair_index[min(u,v),max(u,v)]
                dst.append(q)
                if directed and u>v:flip|=1<<q
        actions.append((dst,flip))

    def act(mask,action):
        dst,flip=action
        ans=flip
        while mask:
            bit=mask&-mask;mask-=bit
            ans^=1<<dst[bit.bit_length()-1]
        return ans

    def adjacency(mask):
        a=[0]*n
        def edge(i,j,v):
            if v:a[i]|=1<<j
            else:a[j]|=1<<i
        for i in range(k):
            edge(i,k,1); edge(k,k+1+i,1)
            edge(i,k+1+i,(mask>>i)&1)
        for z,(i,j) in enumerate(pairs):
            v=(mask>>(k+z))&1
            edge(i,j,v);edge(k+1+j,k+1+i,v)
            v=(mask>>(k+ne+z))&1
            edge(i,k+1+j,v);edge(j,k+1+i,v)
        return a

    def frame_from_adj(a,rho,f):
        ps=[(i,rho[i]) for i in range(n) if i<rho[i]]
        xs=[i if (a[i]>>f)&1 else j for i,j in ps]
        ys=[rho[i] for i in xs]
        mask=sum(((a[xs[i]]>>ys[i])&1)<<i for i in range(k))
        for z,(i,j) in enumerate(pairs):
            mask|=((a[xs[i]]>>xs[j])&1)<<(k+z)
            mask|=((a[xs[i]]>>ys[j])&1)<<(k+ne+z)
        return mask

    @lru_cache(None)
    def hollow_data(hollow):
        images=[act(hollow,ac) for ac in actions]
        h_aut=[p for p,im in zip(perms,images) if im==hollow]
        assert len(h_aut)%2==1
        root_codes=[min(im for p,im in zip(perms,images) if p[i]==0) for i in range(k)]
        least=min(root_codes)
        selected=[i for i in range(k) if root_codes[i]==least]
        assert len(selected)%2==1
        orbit={p[selected[0]] for p in h_aut}
        assert orbit==set(selected)
        chosen=sum(1<<i for i in selected)
        assert all({p[i] for i in selected}==set(selected) for p in h_aut)
        a=adjacency(hollow)
        h=0
        for p in perms:
            for signs in product([0,1],repeat=k):
                half=[i if s==0 else k+1+i for i,s in zip(p,signs)]
                word=half+[k]
                h+=all((a[u]>>v)&1 for u,v in zip(word,word[1:]))
        assert h%2==1
        remaining=set(range(k)); orbit_masks=[]
        while remaining:
            i=min(remaining); orb={p[i] for p in h_aut}; remaining-=orb
            assert len(orb)%2==1
            orbit_masks.append(sum(1<<j for j in orb))
        assert len(orbit_masks)%2==k%2
        return chosen,h,len(h_aut),orbit_masks

    start=time.time(); size=1<<(k*k); class_by_mask=[-1]*size; rows=[]
    for mask in range(size):
        if class_by_mask[mask]>=0:continue
        orbit={act(mask,ac) for ac in actions}
        assert mask==min(orbit)
        ci=len(rows)
        for s in orbit:
            assert class_by_mask[s]<0
            class_by_mask[s]=ci
        stab=math.factorial(k)//len(orbit)
        chosen,h,hsize,orbit_masks=hollow_data(mask & ~((1<<k)-1))
        assert h%stab==0
        assert h%hsize==0 and (h//hsize)%2==1
        a=adjacency(mask);scores=[x.bit_count() for x in a]
        c3=math.comb(n,3)-sum(math.comb(s,2) for s in scores)
        a0=adjacency(mask & ~((1<<k)-1))
        c30=math.comb(n,3)-sum(math.comb(z.bit_count(),2) for z in a0)
        coefficients=[2*k-1-2*a0[i].bit_count() for i in range(k)]
        assert c3==c30+sum(coefficients[i]*((mask>>i)&1) for i in range(k))
        assert c3%2==(k-(mask&((1<<k)-1)).bit_count())%2
        rows.append(dict(id=ci,mask=mask,c3=c3,blue=h//stab,
                         symmetric_HP=h,centralizer=stab,hollow_aut=hsize,
                         hollow_blue_factor=h//hsize,triangle_coefficients=coefficients,
                         selected_pair_mask=chosen,selected_pair_count=chosen.bit_count(),
                         hollow_pair_orbits=orbit_masks))
    # Canonical class involution and all-class assertions.
    for row in rows:
        row['partner']=class_by_mask[row['mask']^row['selected_pair_mask']]
    for row in rows:
        partner=rows[row['partner']]
        assert partner['partner']==row['id']
        assert partner['id']!=row['id']
        assert partner['blue']==row['blue']
        assert (partner['c3']-row['c3'])%2==1
        assert partner['centralizer']==row['centralizer']
        assert partner['symmetric_HP']==row['symmetric_HP']
    # Full independent-orbit toggle cubes in the class set.
    packet_ids={}
    for row in rows:
        generated=[0]
        for om in row['hollow_pair_orbits']:
            generated += [v^om for v in generated]
        ids={class_by_mask[row['mask']^v] for v in generated}
        dim=len(row['hollow_pair_orbits'])
        assert len(ids)==1<<dim
        assert all(rows[i]['blue']==row['blue'] for i in ids)
        assert all(rows[i]['centralizer']==row['centralizer'] for i in ids)
        assert all(rows[i]['symmetric_HP']==row['symmetric_HP'] for i in ids)
        assert sum(rows[i]['c3']%2 for i in ids)==1<<(dim-1)
        root=min(ids); row['cube_packet']=root; row['cube_dimension']=dim
        if root in packet_ids:assert packet_ids[root]==ids
        else:packet_ids[root]=ids
    assert sum(map(len,packet_ids.values()))==len(rows)
    cube_counts=Counter(rows[root]['cube_dimension'] for root in packet_ids)
    # Independent blue HP-cube enumeration. Positional path is 0 -> ... -> n-1.
    coords=[(i,j) for i in range(n) for j in range(i+2,n)]
    ix={ij:q for q,ij in enumerate(coords)}
    refl=[ix[n-1-j,n-1-i] for i,j in coords]
    free=[];used=set()
    for q,r in enumerate(refl):
        if q not in used:
            orbit={q,r};used|=orbit;free.append(sum(1<<z for z in orbit))
    assert len(free)==k*k
    direct=Counter(); rho=list(reversed(range(n)))
    for bits in range(size):
        a=[0]*n
        for i in range(n-1):a[i]|=1<<(i+1)
        t=0
        for z,orb in enumerate(free):
            if (bits>>z)&1:t^=orb
        for q,(i,j) in enumerate(coords):
            if (t>>q)&1:a[j]|=1<<i
            else:a[i]|=1<<j
        frame=frame_from_adj(a,rho,k)
        direct[class_by_mask[frame]]+=1
    assert len(direct)==len(rows)
    assert all(direct[r['id']]==r['blue'] for r in rows)
    # Inherited full tiling census check uses its reported SC classes/blue count,
    # but not its classification algorithm or orbit labels.
    inherited_checked=0
    if inherited_dir:
        path=inherited_dir/f'census_n{n}.json'
        if path.exists():
            old=json.loads(path.read_text());nold=0
            for node in old['nodes']:
                if not node['SC']:continue
                candidates=[t for t,u in enumerate(old['pi']) if u==node['node']]
                t=next(t for t in candidates if all(((t>>i)^(t>>j))&1==0 for i,j in enumerate(old['reflection_permutation'])))
                a=[0]*n
                for i in range(n-1):a[i]|=1<<(i+1)
                for z,(x,y) in enumerate(old['tiles']):
                    i,j=n-x,n-y
                    if (t>>z)&1:a[j]|=1<<i
                    else:a[i]|=1<<j
                mask=frame_from_adj(a,rho,k);row=rows[class_by_mask[mask]]
                assert (row['c3'],row['blue'])==(node['c3'],node['blue'])
                ordinary=old['classes'][node['classes'][0]]
                row['inherited_node']=node['node']
                row['inherited_representative']=node['representative']
                row['mass']=node['mass'];row['type']=node['type']
                row['H_total']=ordinary['H'];row['aut_total']=ordinary['aut']
                nold+=1
            assert nold==len(rows); inherited_checked=nold
    shores=[Counter(r['blue'] for r in rows if r['c3']%2==p) for p in [0,1]]
    assert shores[0]==shores[1]
    result=dict(n=n,k=k,frame_masks=size,classes=len(rows),
                selected_nontrivial_orbit_classes=sum(r['selected_pair_count']>1 for r in rows),
                nontrivial_centralizer_classes=sum(r['centralizer']>1 for r in rows),
                blue_sum=sum(r['blue'] for r in rows),shore_spectrum=dict(sorted(shores[0].items())),
                cube_packet_dimensions=dict(sorted(cube_counts.items())),
                inherited_classes_checked=inherited_checked,
                assertions=['all frame class orbits','hollow group odd','selected orbit odd',
                            'equivariant intrinsic orbit selection','class involution',
                            'no fixed class','blue preserved','c3 parity reversed',
                            'centralizer preserved','symmetric HP count preserved',
                            'independent blue cube equals symmetric-HP/centralizer formula',
                            'free full hollow-orbit toggle cubes','cube B/centralizer/h_rho constant',
                            'cube shores balanced','cube dimension has parity k','integer diagonal triangle affine law',
                            'hollow blue factor odd integer'],
                elapsed_seconds=round(time.time()-start,3),rows=rows)
    (outdir/f'pairing_n{n}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key!='rows'}),flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--max-k',type=int,default=4)
    ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--inherited',type=Path,default=Path('session14-metagraph-census'))
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    reports=[run(k,args.output,args.inherited) for k in range(1,args.max_k+1)]
    (args.output/'verification.json').write_text(json.dumps([{k:v for k,v in r.items() if k!='rows'} for r in reports],indent=2)+'\n')
