from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description="Verify intersection and exact provenance pins without workspace fallback")
ap.add_argument('--repo-root',type=Path,help="Repository root used only for provenance labels beginning math/")
args=ap.parse_args()
repo_root=args.repo_root.resolve() if args.repo_root is not None else None
def pinned_path(name):
 label=Path(name)
 if label.is_absolute() or '..' in label.parts:
  raise ValueError('nonrelative provenance label: '+name)
 if label.parts[0]=='math':
  if repo_root is None:
   raise ValueError('--repo-root is required for repository provenance label: '+name)
  return repo_root.joinpath(*label.parts[1:])
 if label.parts[0]=='session22-incoming-center-review':
  return ROOT.joinpath(*label.parts[1:])
 return ROOT.parent/label
checks=0
def ck(ok,msg):
 global checks
 checks+=1
 if not ok: raise ValueError(msg)
pairs=sorted((q*m,q,m) for q in range(2,9) for m in (1,3,5,7) if q*m>=10 and not(q in (6,8) and m>=5) and not(q==7 and m==7))
ck(sorted({n for n,q,m in pairs})==[10,12,14,15,18,20,21,24,25,28,35],'incoming support')
remaining=[p for p in pairs if not(p[1]==2 and p[2]>=5)]
orders=sorted({n for n,q,m in remaining})
ck(orders==[12,15,18,20,21,24,25,28,35],'nine-order support')
ck([p for p in pairs if p not in remaining]==[(10,2,5),(14,2,7)],'exact removed rows')
ck(max(orders)==35 and 35+12+1==48,'unchanged degree and safe onset')
for name,h in json.loads((ROOT/'source_pins.json').read_text()).items():
 ck(hashlib.sha256(pinned_path(name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==h,'source pin '+name)
print(json.dumps({'checks':checks,'remaining_orders':orders,'parameter_rows':remaining,'degree_bound':35,'safe_onset':48,'scope':'intersection bookkeeping, not a center census'},sort_keys=True))
