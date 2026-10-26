#!/usr/bin/env python3
"""Read exact source blobs into this task's directory; never mutate source repo."""
import hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parent/'session67-repo'
PIN='b4234f728c0516b5fc521c001d69e26f3d2508d5'
PATHS=[
 '04-computation/information_theory_deep.py',
 '05-knowledge/results/information_theory_deep.out',
 '04-computation/ocf_information_theory.py',
 '05-knowledge/results/ocf_information_theory.out',
 '05-knowledge/hypotheses/HYP-3564-relations-not-things-the-relational-reframe.md',
 '05-knowledge/hypotheses/HYP-3595-rib-sc-side-even-and-silent-flip-eulerian.md',
 '05-knowledge/hypotheses/HYP-3602-tournaments-are-intransitivity-among-n-things.md',
 '01-canon/theorems/THM-1382-map-graph-star-flips-nested-cut-cycle.md',
 '01-canon/theorems/THM-1979-tournament-space-is-a-spectrum-single-point-to-continuum.md',
 '01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md',
 '05-knowledge/results/session36_information_curl_contexts_20261023/context/CURL-AND-CLAN-COMPLETION.md',
 '05-knowledge/results/session36_information_curl_contexts_20261023/context/REALIZATION_PHASES.md',
 '05-knowledge/results/session36_information_curl_contexts_20261023/lattice/METRIC_CURL_RECOVERY.md',
]
def git(*args): return subprocess.check_output(['git','-C',str(REPO),*args])
rows=[]
for i,p in enumerate(PATHS,1):
 raw=git('show',f'{PIN}:{p}')
 dst=ROOT/'sources'/f'{i:02d}-{Path(p).name}'
 dst.write_bytes(raw)
 rows.append({'id':i,'path':p,'pin':PIN,'blob':git('rev-parse',f'{PIN}:{p}').decode().strip(),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'local_copy':str(dst.relative_to(ROOT)),'last_change':git('log','-1','--format=%H %aI %s',PIN,'--',p).decode().strip(),'url':f'https://github.com/eliottcassidy2000/math/blob/{PIN}/{p}'})
(ROOT/'SOURCE_MANIFEST.json').write_text(json.dumps({'source_pin':PIN,'source_repo_read_only':str(REPO),'sources':rows},indent=2)+'\n')
print(json.dumps({'pin':PIN,'source_count':len(rows),'sha256_all_verified':all(hashlib.sha256((ROOT/x['local_copy']).read_bytes()).hexdigest()==x['sha256'] for x in rows)},indent=2))
