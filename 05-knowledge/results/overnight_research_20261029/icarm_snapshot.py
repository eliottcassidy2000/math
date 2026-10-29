"""Read-only public-source capture and exact point incidence checks."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from datetime import datetime,timezone
import hashlib,json,urllib.request,importlib.util

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
def fetch(url):
    raw=urllib.request.urlopen(url,timeout=75).read(25000000)
    return json.loads(raw),hashlib.sha256(raw).hexdigest()
def main():
    curve,sha=fetch('https://elliptic-rank.icarm.cloud/curve/273.json')
    (ROOT/'curve273.json').write_text(json.dumps(curve,indent=2)+'\n',encoding='utf-8')
    a1,a2,a3,a4,a6=map(int,curve['ainvs'])
    for x,y in curve['points']:
        x,y=F(x),F(y)
        if y*y+a1*x*y+a3*y!=x**3+a2*x*x+a4*x+a6:raise RuntimeError('point incidence')
    oldpath=REPO/'04-computation/elliptic_rank30_optimal_local_observer_certificate_codex_20260821.py'
    spec=importlib.util.spec_from_file_location('inherited_certificate',oldpath)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    if (a1,a2,a3,a4,a6)!=(1,0,0,old.A4,old.A6):raise RuntimeError('changed curve')
    if [tuple(map(F,p)) for p in curve['points']]!=list(old.POINTS):raise RuntimeError('changed ordered points')
    data,dbsha=fetch('https://elliptic-rank.icarm.cloud/database.json')
    curves=data['curves']
    if len(curves)!=data['count']:raise RuntimeError('database count')
    ranks=Counter(c['rank_lower_bound'] for c in curves)
    high=[{k:c.get(k) for k in ('id','curve_key','rank_lower_bound','family','families','created_at','updated_at','commentary')}
          for c in curves if c['rank_lower_bound']>=28]
    out={'retrieved_utc':datetime.now(timezone.utc).isoformat(),
         'source_metadata':{k:v for k,v in data.items() if k!='curves'},
         'curve273_source_sha256':sha,'database_source_sha256':dbsha,
         'curve273_points_checked':len(curve['points']),'matches_august_ordered_certificate':True,
         'rank_lower_bound_histogram':dict(sorted(ranks.items())),
         'at_least_22':sum(v for k,v in ranks.items() if k>=22),'at_least_28':len(high),
         'high_rank_records':high,
         'scope':'Selected leaderboard records, not a height-ordered unbiased population. Incidence and identity checked here; independence replay in the inherited certificate. No rank upper bound.'}
    (ROOT/'icarm_snapshot.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('high_rank_records','source_metadata')},indent=2))
if __name__=='__main__':main()
