"""Independent literal DP, polynomial replay, isomorphisms and hostile controls."""
from collections import Counter
from itertools import combinations
import importlib.util,json,sys
from common import *

def extend(rows,born,mask):
    n=len(rows);out=list(rows)+[0]*born;k=0
    for i,j in combinations(range(n+born),2):
        if j<n:continue
        u,v=(i,j) if mask>>k&1 else (j,i);out[u]|=1<<v;k+=1
    return tuple(out)

def main():
    old=json.loads((OLD/'branches.json').read_text());parents=[tuple(v['rows']) for v in old['H55_sources']]+[tuple(old['positive55to111']['core'])]
    for filename,ids in (('parents55.txt',list(range(7))),('representatives55.txt',[0,2,6])):
        lines=[list(map(int,line.split())) for line in (HERE/filename).read_text().splitlines() if line.strip()]
        need(lines==[[i,len(parents[i]),*parents[i]] for i in ids],'C++ inputs retain the literal pinned parent rows')
    mod=old_graph_module();sampled=0;receipts=[]
    for born,name in ((2,'two_birth.json'),(3,'three_birth.json')):
        data=json.loads((HERE/name).read_text(encoding='utf-8-sig'))
        need([d['id'] for d in data]==(list(range(7)) if born==2 else [0,2,6]),'complete parent selection')
        for d in data:
            i=d['id'];r=parents[i];bits=born*len(r)+born*(born-1)//2;size=1<<bits
            need(d['variables']==bits and d['born']==born and d['transforms_agree_every_mask'],'declared Boolean universe')
            need(sum(d['histogram'].values())==size,'full orientation census')
            need(d['H111_count']==(0 if i<6 else 4 if born==2 else 18),'hostile and positive controls')
            expected=None
            if born==2:
                expected=[0]*size
                for mask,c in d['coefficients']:
                    need(mask.bit_count()<=4,'quadratic-per-new-vertex degree bound');expected[mask]=c
                for bit in range(bits):
                    for mask in range(size):
                        if mask>>bit&1:expected[mask]+=expected[mask^(1<<bit)]
                need(dict(Counter(expected))=={int(h):c for h,c in d['histogram'].items()},'independent Python zeta replay')
            samples={0,size-1,*[1<<b for b in range(bits)],*[size-1-(1<<b) for b in range(bits)],*d['H111_first_masks']}
            random_state=55+100*born+i
            for _ in range(20):
                random_state^=(random_state<<13)&0xffffffff
                random_state^=random_state>>17
                random_state^=(random_state<<5)&0xffffffff
                samples.add(random_state&(size-1))
            need(samples=={int(m) for m in d['probe_H']},'complete independent probe selection')
            for mask in sorted(samples):
                target=extend(r,born,mask);h=mod.subset_h(target);sampled+=1
                if expected is not None:need(h==expected[mask],'literal subset DP versus complete polynomial')
                need(h==d['probe_H'][str(mask)],'literal subset DP versus same-mask transform value')
                if mask in d['H111_first_masks']:
                    need(h==111,'literal positive extension')
                    receipts.append({'parent':i,'born':born,'mask':mask,'H':h,'rows':target})
    # Replay the one-birth obstruction; source labels and each named port remain intact.
    for i,r in enumerate(parents[:6]):
        counts=[mod.subset_h(extend(r,1,m)) for m in range(1<<len(r))]
        need(111 not in counts and sorted(counts)==sorted(old['H55_sources'][i]['attachment_H']),'one-birth hostile unchanged')
    f=HERE.parent/'regular_parent_symmetry_20261029/verify.py';s=importlib.util.spec_from_file_location('isomorphism_control',f)
    iso=importlib.util.module_from_spec(s);s.loader.exec_module(iso)
    maps=[]
    for representative,indices in ((0,(0,1,4)),(2,(2,3,5))):
        for i in indices:
            mapping=iso.iso_maps(parents[representative],parents[i],stop=True)[0]
            need(all(((parents[representative][u]>>v)&1)==((parents[i][mapping[u]]>>mapping[v])&1) for u in range(7) for v in range(7)),'literal isomorphism map')
            maps.append({'representative':representative,'source':i,'permutation':mapping})
    need(not iso.iso_maps(parents[0],parents[2]) and all(not iso.iso_maps(parents[6],parents[i]) for i in (0,2)),'three distinct ordinary parents')
    prior=HERE.parent/'primegame_port_spectrum_20261103';sys.path.insert(0,str(prior))
    s=importlib.util.spec_from_file_location('port_hostile',prior/'verify.py');hm=importlib.util.module_from_spec(s);s.loader.exec_module(hm)
    hostile=hm.hostile();need(hostile['groups']==[81,27] and hostile['port8_after_sink']==[2390013,1754919],'81/27 actual control')
    result={'status':'PASS','literal_sampled_extensions':sampled,'isomorphism_maps':maps,'positive_extensions':receipts,
            'hostile81_27':hostile,'two_birth_hostile_orientations':6*2**15,'three_birth_hostile_orientations':2*2**24,
            'positive_two_birth_count':4,'positive_three_birth_count':18,
            'scope':'No induced H111 extension with at most three added vertices for any of the six bank H55 presentations; four or more remain open.'}
    (HERE/'birth_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','literal_samples':sampled,'hostile_three_birth_orientations':2*2**24,'positive_counts':[4,18],'hostile':hostile}))
if __name__=='__main__':main()
