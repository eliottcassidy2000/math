"""Emit an actual strong source and optional legal C,R,C certificate."""
import argparse
import json
from pathlib import Path
from transport import profiles, covers, cyclic, transitive, family_value, family_rewrite, H, require


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',choices=profiles(),default='B11')
    parser.add_argument('--m',type=int,default=2)
    parser.add_argument('--short-arm',type=int,choices=[1,2],default=2)
    parser.add_argument('--rewrite-to',type=int)
    parser.add_argument('--count-direct',action='store_true')
    args=parser.parse_args()
    require(args.m>=1,'m must be positive')
    base=profiles()[args.base]
    rows=cyclic(base,transitive(args.m),transitive(args.short_arm))
    result={'base':args.base,'m':args.m,'short_arm':args.short_arm,
            'source_rows':rows,'H_source':family_value(covers(base),args.m,args.short_arm),
            'strong':'three nonempty cyclic blocks','command_count':0}
    if args.rewrite_to is not None:
        receipt=json.loads((Path(__file__).parent/'verification.json').read_text(encoding='utf-8'))
        key=args.base+'_'+str(args.rewrite_to)
        require(key in receipt['rewrites'],'requested rewrite not in certified atlas')
        step=receipt['rewrites'][key]['certificate']
        step['path_before']=tuple(step['path_before'])
        out=family_rewrite(base,step,args.m,args.short_arm)
        result.update({'command_count':3,'commands':['complement','interior relocation','complement'],
                       'certificate':out,
                       'H_target':family_value(covers(step['rows']),args.m,args.short_arm)})
    if args.count_direct:
        require(len(rows)<=20,'direct subset check is limited to order20')
        require(H(rows)==result['H_source'],'direct source count')
        if 'certificate' in result:
            require(H(result['certificate']['rows'])==result['H_target'],'direct target count')
        result['direct_count_verified']=True
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
