"""Uniform always-ready rebase, with count-preserving boundary exactly k2."""
from pathlib import Path
import argparse,json,sys
from classify_modules import H,strong,need,inherited as a,two_birth as b
from composition import family
HERE=Path(__file__).resolve().parent


def rebase(k):
    source=family(k,1)
    kept=[v for v in range(len(source)) if v not in(k+3,k+5)]
    parent=a.induced(source,kept)
    child=a.birth_beaten_only_by(parent,len(parent)-1)
    return source,kept,parent,child


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
    report={'status':'PASS','scope':'Literal k1..6 controls and complete fixed-placeholder coefficient proof',
            'parent_runs':wr,'child_runs':rr,'witnesses':witnesses}
    (HERE/'ready_rebase.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS: always1-ready rebase; unique425 same-count boundary; literal1277 continuation')
    if args.inverse:
        sys.path.insert(0,str(HERE.parent.parent/'ready_return_pentagon_20261103'))
        import ready_inverse as inv
        bound,counts,values=inv.ready_values(425)
        need(bound==10 and 425 in values[1] and 425 not in values[2],'complete ready-value boundary425')
        target=values[1][425];P=a.root_profiles(target['T'])[target['root']]
        need(P['h']==425 and P['gap']==1,'inverse witness independently counted')
        data={'status':'PASS','limit':425,'all_size_parent_bound':bound,'parent_class_counts':counts,
              '425_lanes':[1],'425_witness':target,
              'lane_values':{str(lane):sorted(bank) for lane,bank in values.items()}}
        (HERE/'ready425_inverse.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
        print('PASS: complete all-order425 ready image has lane1 and excludes lane2')


if __name__=='__main__':main()
