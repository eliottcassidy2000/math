"""Uniform always-ready rebase, with count-preserving boundary exactly k2."""
from pathlib import Path
import argparse,json,sys
from classify_modules import H,strong,need,inherited as a,two_birth as b,converse,canonical_and_aut
from composition import family
HERE=Path(__file__).resolve().parent


def rebase(k):
    source=family(k,1)
    kept=[v for v in range(len(source)) if v not in(k+3,k+5)]
    parent=a.induced(source,kept)
    child=a.birth_beaten_only_by(parent,len(parent)-1)
    return source,kept,parent,child


def full_state(T,r):
    p=a.root_profiles(T)[r];c=a.root_profiles(converse(T))[r]
    W=a.induced(T,[v for v in range(len(T)) if v!=r]);hp,_=a.dp(W);full=len(hp)-1
    K=sum(hp[s]*hp[full^s] for s in range(full+1))
    return [p['h'],c['s'],p['s'],c['F'],p['F'],p['q'],K]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--inverse',action='store_true')
    args=parser.parse_args()
    _,_,W,R=rebase(1)
    wr,rr=b.module_runs(W),b.module_runs(R)
    need(wr[1:]==[31,11,0,0,0,0,0],'complete parent placeholder coefficients')
    need(rr[1:]==[175,125,24,0,0,0,0,0],'complete child placeholder coefficients')
    witnesses=[]
    for k in range(1,7):
        source,kept,W,R=rebase(k);D=1<<k;expected=24*3**k+53*D-3
        need(H(W)==11*D+9 and H(R)==expected,'independent literal all-k response controls')
        need(strong(W) and strong(R),'strong parent and returned source')
        need(a.induced(R,range(len(W)))==a.induced(source,kept),'every persistent arc preserved')
        prof=a.root_profiles(R)[-1]
        need((prof['s'],prof['q'],prof['gap'],prof['delta'])==(11*D+8,11*D+9,1,0),'ready profile')
        pz=a.root_profiles(W)[-1]
        need(pz['s']==1,'unique-start parent')
        # Old labels z,B,A,a,r,prefix,y form the structurally forced path.
        old_path=[k+7,k,k+1,k+2,k+4,*range(k),k+6]
        pos={old:i for i,old in enumerate(kept)}
        path=[pos[v] for v in old_path]
        need(len(path)==len(W) and len(set(path))==len(W),'spanning unique-start witness')
        need(all(W[x]>>y&1 for x,y in zip(path,path[1:])) and W[path[-1]]>>path[0]&1,
             'spanning cycle proving parent strongness')
        first=a.birth_beaten_only_by(R,len(R)-1)
        final=a.birth_beaten_only_by(first,len(first)-1)
        need(H(first)==48*3**k+95*D-13 and H(final)==3*expected+2,'literal ready continuation')
        need((H(R)==H(source))==(k==2),'unique same-count boundary')
        witnesses.append({'k':k,'source_rows':source,'kept_old_labels':kept,
                          'parent_rows':W,'parent_unique_start_path':path,'target_rows':R,
                          'profile':prof,'trace':[H(source),H(W),H(R),H(first),H(final)]})
    source,_,_,returned=rebase(2)
    before=full_state(source,6);after=full_state(returned,8)
    need(before==[425,4,41,2542,449,41,2848] and after==[425,15,52,1491,426,53,1544],
         'independent full rooted response for actual425 rebase')
    delta=[x-y for x,y in zip(before,after)]
    need(all(delta[j]>=0 for j in(0,3,4,6)) and all(delta[j]<=0 for j in(1,2,5)),
         'incoming two-sided amalgamation comparison cone')
    cards=[]
    targets=[canonical_and_aut(T)[0] for T in(returned,converse(returned))]
    for drop in range(len(source)):
        W=a.induced(source,[v for v in range(len(source)) if v!=drop])
        code=canonical_and_aut(W)[0]
        need(code not in targets,'no nine-retained ready425 target in either orientation')
        cards.append({'deleted':drop,'H':H(W),'canonical':code})
    report={'status':'PASS','scope':'Literal k1..6 controls and complete fixed-placeholder coefficient proof',
            'parent_runs':wr,'child_runs':rr,'witnesses':witnesses,
            '425_full_state_before':before,'425_full_state_after':after,
            '425_single_deletion_cards':cards}
    (HERE/'ready_rebase.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS: always1-ready rebase; unique425 same-count boundary; literal1277 continuation')
    if args.inverse:
        sys.path.insert(0,str(HERE.parent.parent/'ready_return_pentagon_20261103'))
        import ready_inverse as inv
        bound,counts,values=inv.ready_values(425)
        need(bound==10 and 425 in values[1] and 425 not in values[2],'complete ready-value boundary425')
        target=values[1][425];P=a.root_profiles(target['T'])[target['root']]
        need(P['h']==425 and P['gap']==1,'inverse witness independently counted')
        all_one=[]
        for n in range(3,bound+1):
            for j in range(inv.class_counts(n)[n]):
                T=inv.unrank(n,j);h,s,q,F=inv.rooted_data(T,0)
                if h+F-1==425:
                    U,r=inv.birth(T,0)
                    all_one.append({'rows':U,'root':r,'parent_order':n,'parent_rank':j})
                inv.tables.cache_clear()
        need(len(all_one)==1,'unique rooted strong1-ready425 class at every order')
        need(canonical_and_aut(all_one[0]['rows'])[0]==canonical_and_aut(returned)[0],
             'constructed return is the unique inherited425 ready class')
        data={'status':'PASS','limit':425,'all_size_parent_bound':bound,'parent_class_counts':counts,
              '425_lanes':[1],'425_witness':target,
              'all_rooted_1_ready425_classes':all_one,
              'lane_values':{str(lane):sorted(bank) for lane,bank in values.items()}}
        (HERE/'ready425_inverse.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
        print('PASS: complete all-order425 ready image has lane1 and excludes lane2')


if __name__=='__main__':main()
