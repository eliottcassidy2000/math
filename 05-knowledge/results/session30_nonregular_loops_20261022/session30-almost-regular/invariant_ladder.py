import json,pathlib,sys,sympy as sp,numpy as np
from experiment import frombits,pred
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'session30-decorated-engine'));from predict import predict
r=next(r for r in json.load(open(ROOT/'isospectral_witnesses.json')) if r['n']==8);out=[]
for side in ['first','second']:
 a=r[side];S=frombits(8,a['bits']);M=sp.Matrix(S)/8;Q=M*M;I=sp.eye(8);G=-Q*(I+Q).inv();B=Q*(I-Q).inv();O=M*(I-Q).inv();d=sp.Matrix([B[i,i]-G[i,i] for i in range(8)]);one=sp.ones(8,1);delta=M*one
 f8=sum(2*G[i,i]*B[i,i]-sp.Rational(3,2)*G[i,i]**2 for i in range(8));theta=sum(G[i,j]*(B[i,j]**2-O[i,j]**2-sp.Rational(2,3)*G[i,j]**2) for i in range(8) for j in range(8));dumb=-(d.T*G*d)[0];L=-2*(d.T*G*one)[0];tree=8*(delta.T*(I+Q).inv()*delta)[0]
 exact={k:str(v) for k,v in dict(f8=f8,theta=theta,dumbbell=dumb,lollipop=L,tree=tree,X1=f8+theta+dumb).items()};target=pred(S,int(a['H']));pp=predict(S);out.append(dict(side=side,bits=a['bits'],H=a['H'],rational=exact,target=target,prediction=pp));print(side,exact,'y',target['y'],'ypred1',pp['y1'],'ypred2',pp['y2'],flush=True)
json.dump(out,open(ROOT/'invariant_ladder.json','w'),indent=2)
