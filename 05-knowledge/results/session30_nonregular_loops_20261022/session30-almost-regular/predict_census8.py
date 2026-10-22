import pathlib,json,subprocess,sys,math,collections,itertools
from fractions import Fraction
from experiment import frombits,pred
from forest_exact import forest
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'session30-decorated-engine'));from predict import predict
bank=json.load(open(ROOT/'census8_bank.json'));proc=subprocess.Popen([str(ROOT/'count')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True);out=[]
with open(ROOT/'census8_predictions.jsonl','w') as f:
 for i,a in enumerate(bank['almostregular8']):
  b=a['bits'];S=frombits(8,b);proc.stdin.write(f'8 {b}\n');proc.stdin.flush();H=int(proc.stdout.readline());N=forest(S);assert sum(Fraction(v,math.factorial(8)//math.factorial(8-L)) for L,v in enumerate(N))==Fraction(H*128,math.factorial(8))
  target=pred(S,H);pp=predict(S);r=dict(index=i,bits=b,H=H,aut=a['aut'],labelled_weight=math.factorial(8)//a['aut'],N=N,target=target,prediction=pp,error1=target['y']-pp['y1'],error2=target['y']-pp['y2']);out.append(r);f.write(json.dumps(r)+'\n');f.flush();print(i,H,r['error1'],r['error2'],flush=True)
proc.stdin.close();proc.wait()
summary={};summary['classes']=len(out);summary['labelled']=sum(r['labelled_weight'] for r in out);summary['H_histogram_unlabelled']=dict(sorted(collections.Counter(r['H'] for r in out).items()));summary['H_histogram_labelled']={h:sum(r['labelled_weight'] for r in out if r['H']==h) for h in sorted(set(r['H'] for r in out))};summary['second_improves']=sum(abs(r['error2'])<abs(r['error1']) for r in out);summary['second_worsens']=[r['index'] for r in out if abs(r['error2'])>=abs(r['error1'])];summary['worst_absolute_first']=max(out,key=lambda r:abs(r['error1']))['index'];summary['worst_absolute_second']=max(out,key=lambda r:abs(r['error2']))['index'];summary['pairwise_ranking']={}
for name,fn in [('rho',lambda r:r['target']['rho']),('kappa1',lambda r:r['prediction']['kappa1']),('kappa2',lambda r:r['prediction']['kappa2'])]:
 pairs=[(a,b) for a,b in itertools.combinations(out,2) if a['H']!=b['H']];wrong=[(a['index'],b['index']) for a,b in pairs if (a['H']-b['H'])*(fn(a)-fn(b))<0];summary['pairwise_ranking'][name]=dict(strict_H_pairs=len(pairs),reversed=len(wrong),fraction=len(wrong)/len(pairs),pairs=wrong)
json.dump(summary,open(ROOT/'census8_summary.json','w'),indent=2);print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='pairwise_ranking'}),flush=True)
