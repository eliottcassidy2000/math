"""Bounded nonminimal next-ready target search after migrated return."""
from itertools import combinations
from collections import Counter
import json,time,argparse
import sparse_repair as s
Q=(14,12,8,64,175,79,535,367,127,447)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--radius',type=int,default=2);args=ap.parse_args()
    X=s.c.birth(Q,9);edges=list(combinations(range(11),2));out=[];stats=Counter();tic=time.time()
    for radius in range(args.radius+1):
        for flips in combinations(edges,radius):
            parent=s.flip(X,flips);h,starts,F,counts=s.c.probe.profiles(parent);stats[radius]+=1
            for port,H in enumerate(counts):
                if H!=881 or starts[port] not in (1,2):continue
                child=s.c.birth(parent,port);mods=[]
                for v in range(11):
                    if v==port:continue
                    hc=s.exp_coeff(s.module_coeff(child,v));dc=s.exp_coeff(s.module_coeff(parent,v,start=port))
                    if hc=={1:71,2:405}:mods.append(dict(module=v,H=hc,defect=dc))
                out.append(dict(radius=radius,flips=flips,port=port,defect=starts[port],strong=s.c.strong(child),rows=child,uniform_modules=mods))
        print('radius',radius,'parents',stats[radius],'small-defect hits',len(out),'uniform',sum(bool(h['uniform_modules']) for h in out),'seconds',round(time.time()-tic,2),flush=True)
    (s.HERE/('next_repair_radius'+str(args.radius)+'.json')).write_text(json.dumps(dict(status='FINITE-EXACT bounded repair search',scope='Q standard advance: every parent arc subset of radius <= requested bound, all latest-root sole predecessors, target881 defects1/2, every target one-vertex module except final root/predecessor',max_radius=args.radius,stats=dict(stats),hits=out),indent=2)+'\n')
if __name__=='__main__':main()
