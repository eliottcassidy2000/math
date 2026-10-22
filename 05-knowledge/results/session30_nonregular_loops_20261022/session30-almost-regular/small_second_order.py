import json,sys,pathlib
from experiment import frombits,pred
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'session30-decorated-engine'));from predict import predict
out=[]
for n,rs in json.load(open(ROOT/'small_classes.json')).items():
 for i,r in enumerate(rs):
  S=frombits(int(n),r['bits']);p=predict(S);target=pred(S,int(r['H']));a=dict(n=int(n),class_index=i,bits=r['bits'],H=r['H'],y=target['y'],prediction=p,err1=target['y']-p['y1'],err2=target['y']-p['y2']);out.append(a);print(n,i,'H',r['H'],'y',target['y'],'y1',p['y1'],'y2',p['y2'],'improve',abs(a['err2'])<abs(a['err1']),flush=True)
json.dump(out,open(ROOT/'small_second_order.json','w'),indent=2)
