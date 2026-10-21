#!/usr/bin/env python3
"""Direct exact hostile-ray replay; no finite-interval inference."""
import argparse,hashlib,json,gc
from pathlib import Path
import verify_annulus as a
HERE=Path(__file__).resolve().parent
SCALES=(7,448,28672,8,512,32768)
def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    data=json.loads((HERE/'frontier.json').read_text());results=[]
    for N in SCALES:
        routes=a.Routes(data['frontier'],dict(data['frontier_ranks']))
        e,m,old=a.direct(N,routes);E,M,large=a.direct(64*N,routes)
        a.check(M[-1]-m[-1]==153*N,'ray annulus source budget')
        rows=[[d,e[d],E[d],m[d],M[d],E[d]-e[d],E[d]-e[d]-M[d]+m[d],153*N-E[d]+e[d]] for d in range(a.DEPTH_CAP+1)]
        peak=max(r[5] for r in rows);first=next(r[0] for r in rows if r[5]==peak)
        bad=[r[0] for r in rows if r[7]<0]
        a.check(bool(bad)==(N in (7,8)),'declared hostile-ray outcome')
        dag=routes.audit()
        results.append({'N':N,'old':old,'large':large,'max_increment':peak,'first_peak_depth':first,'budget':153*N,'max_all_cross':max(r[6] for r in rows),'first_full_budget_failure':bad[0] if bad else None,'horizons':rows,'route_dag':dag})
        del routes;gc.collect()
    result={'status':'FINITE-EXACT six stated base scales only; no all-ray inference','scales':SCALES,'horizon_schema':['d','old_E','large_E','arrived_old_norm','arrived_large_norm','increment','all_cross','full_bank'],'results':results,'checks':a.CHECKS,'frontier_sha256':hashlib.sha256((HERE/'frontier.json').read_bytes()).hexdigest()}
    blob=a.encode(result);out=HERE/'ray_receipt.json'
    if args.write:out.write_bytes(blob)
    else:a.check(out.read_bytes()==blob,'saved hostile-ray receipt equality')
    print(json.dumps({'status':'PASS','checks':result['checks'],'receipt_sha256':hashlib.sha256(blob).hexdigest(),'peaks':[[r['N'],r['max_increment'],r['first_peak_depth'],r['budget'],r['max_all_cross']] for r in results]},sort_keys=True))
if __name__=='__main__':main()
