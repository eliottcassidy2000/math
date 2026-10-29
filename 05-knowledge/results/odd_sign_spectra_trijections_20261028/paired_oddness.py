"""Complete aligned-pair sign modes and converse/mask actions through core four."""
from itertools import combinations, permutations
from collections import Counter
from pathlib import Path
import json
from sign_spectrum import action, linear
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,msg):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(msg)

def run(n):
    edges=list(combinations(range(n),2));e=len(edges);ec=(1<<e)-1;vc=(1<<n)-1
    perms=list(permutations(range(n)))
    transforms=[]
    for p in perms:
        dest,flip=action(edges,p)
        transforms.append(([linear(x,dest) for x in range(1<<e)],flip,
                           [linear(x,p) for x in range(1<<n)]))
    def images(x,signed):
        a=x&ec;b=(x>>e)&ec;c=x>>(2*e)
        return [((am[a]^f) | ((am[b]^f)<<e) | (cm[c]<<(2*e))) if signed
                else (am[a] | (am[b]<<e) | (cm[c]<<(2*e))) for am,f,cm in transforms]
    universe=1<<(2*e+n)
    representatives=[];owner={};mode_seen=set();mode_parities=Counter();joint=0;separate=0;cancel_witness=None
    for x in range(universe):
        if x not in owner:
            orbit=set(images(x,True));q=min(orbit)
            representatives.append(q)
            for y in orbit:owner[y]=q
        if x not in mode_seen:
            orbit=images(x,False);mode_seen.update(orbit)
            a=x&ec;b=(x>>e)&ec;c=x>>(2*e)
            stabilizer=[j for j,y in enumerate(orbit) if y==x]
            ca=[(a&transforms[j][1]).bit_count()%2 for j in stabilizer]
            cb=[(b&transforms[j][1]).bit_count()%2 for j in stabilizer]
            if not any(ca) and not any(cb):separate+=1
            if ca==cb:
                joint+=1;mode_parities[(a.bit_count()%2,c.bit_count()%2)]+=1
                if any(ca) and cancel_witness is None:cancel_witness={'support_I':a,'support_J':b,'mask_support':c,'stabilizer_signs':ca}
    need(joint==len(representatives),'joint sign-even modes equal actual paired classes')
    flips=[((ec if x&1 else 0) | ((ec if x&2 else 0)<<e) | ((vc if x&4 else 0)<<(2*e))) for x in range(8)]
    fixed=[sum(owner[x^f]==x for x in representatives) for f in flips]
    need(all(fixed[j]==0 for j in (1,2,4,5,6)),'all proper core flips and mask flips are free')
    need(len(set(mode_parities.values()))==1 and len(mode_parities)==4,'four Walsh parity sectors equal')
    need(next(iter(mode_parities.values()))==len(representatives)//4,'free Klein action mode division')
    if n%2:need(fixed[7]==0,'odd core anti-automorphism fixed vertex blocks complemented mask')
    else:need(fixed[3]==fixed[7],'even core two anti-mask channels agree')
    seen=set();orbit_sizes=Counter()
    for x in representatives:
        if x in seen:continue
        orbit={owner[x^f] for f in flips};seen.update(orbit);orbit_sizes[len(orbit)]+=1
    need(set(orbit_sizes)<={4,8},'full flip stabilizer is at most diagonal order two')
    need(orbit_sizes[4]==(fixed[3]+fixed[7])//4,'four-orbit formula')
    need(orbit_sizes[8]==(len(representatives)-fixed[3]-fixed[7])//8,'eight-orbit formula')
    return {'core_order':n,'self_converse_target_order':2*n+1,'labelled_inputs':universe,
            'paired_classes':len(representatives),'joint_sign_even_modes':joint,
            'separately_positive_modes':separate,'cancellation_only_modes':joint-separate,
            'cancellation_witness':cancel_witness,
            'walsh_T_and_mask_parity_sector_counts':{str(k):v for k,v in sorted(mode_parities.items())},
            'fixed_classes_by_flip_bits_T_W_mask':fixed,'full_flip_orbit_sizes':dict(orbit_sizes)}

def main():
    rows=[run(n) for n in (2,3,4)]
    need([r['paired_classes'] for r in rows]==[8,88,2752],'classical self-converse controls')
    out={'status':'FINITE-EXACT','rows':rows,'checks':CHECKS}
    (HERE/'paired_oddness.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out))
if __name__=='__main__':main()
