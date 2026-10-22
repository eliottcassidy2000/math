import json,pathlib,sys
from experiment import frombits
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'session30-decorated-engine'));from predict import predict
cases=[json.loads(x) for x in open(ROOT/'second_order.jsonl')];selected=[cases[0],cases[-1],cases[-4]];out=[]
for r in selected:
 p=predict(frombits(r['n'],r['bits']),derivative_radius=.0125,derivative_points=32);diff={k:p[k]-r['prediction'][k] for k in ['D_A1','D2_A1','y1','y2']};out.append(dict(name=r['name'],differences=diff));print(out[-1],flush=True)
json.dump(out,open(ROOT/'derivative_stability.json','w'),indent=2)
