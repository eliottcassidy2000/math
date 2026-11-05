"""Physical primitive-interface, decoder, and sparse-language checks."""
import importlib.util,hashlib,itertools,json,math
from pathlib import Path
p=Path(__file__).resolve().parent.parent/'ribbon-terminals'/'verify.py'
if hashlib.sha256(p.read_bytes()).hexdigest()!='9ada1cc1cfaf33030a89122dcf4ed65103aa23f5a2ad01e5dd5c288cde198a3e':raise RuntimeError('dependency hash changed')
spec=importlib.util.spec_from_file_location('terminals',p);t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t)
require=t.require
cells=[dict(rows=[6,4,8,3],b=2,c=3,cycle=[0,1,2,3],word='0013',M=[[5,-1,7],[1,0,1],[7,-1,14]]),dict(rows=[54,60,56,1,8,24],b=2,c=5,cycle=[0,1,2,5,4,3],word='001146',M=[[17,-4,28],[4,-1,6],[28,-6,51]]),dict(rows=[56,33,3,22,6,28],b=0,c=5,cycle=[0,3,4,1,5,2],word='446637',M=[[41,-7,49],[5,-1,6],[75,-12,93]])]
seed=[4,1,2];initial=(3,1,3)

def admissible(h,e):return h>=3 and h%2==1 and 1<=e<=h-2

def decode(h,e):
    out=[]
    while (h,e)!=(3,1):
        if not admissible(h,e):return None
        if 7*e>h and 5*e<h:
            a=0;num=(-h+8*e+2,-3*h+19*e-4);den=(5,5)
        elif 5*e>h and 4*e<h:
            a=1;num=(-7*h+32*e-8,-16*h+73*e-20);den=(1,1)
        elif 9*e>h and 7*e<h:
            a=2;num=(h-8*e+2,17*h-139*e-2);den=(3,21)
        else:return None
        if any(n%d for n,d in zip(num,den)):return None
        hh,ee=(n//d for n,d in zip(num,den))
        if not admissible(hh,ee) or hh>=h:return None
        out.append(a);h,e=hh,ee
    return tuple(reversed(out))

def main():
    for a,z in enumerate(cells):
        G=z['rows'];cy=z['cycle'];M=z['M']
        require(t.matrix(G,z['b'],z['c'])==M,('matrix',a))
        require(t.span(G)<=3 and z['b']<3 and z['c']==len(G)-1,('span premises',a))
        require(sorted(cy)==list(range(len(G))) and all(G[cy[i]]>>cy[(i+1)%len(G)]&1 for i in range(len(G))),('cycle',a))
        require(t.rm([-2,1,1],M)==[-2,1,1],('charge',a))
    physical=[]
    for d in range(3):
        for word in itertools.product(range(3),repeat=d):
            rows=seed;state=initial;mask='044'
            for a in word:
                z=cells[a];rows=t.glue(rows,len(rows)-1,z['rows'],z['b']);state=t.mv(z['M'],state);mask+=z['word']
            require(rows==t.wordgraph(mask),('mask word',word))
            require(t.profile(rows,len(rows)-1)==state,('physical profile',word))
            require(t.span(rows)<=3,('physical span',word))
            physical.append(dict(cells=''.join('ABC'[a] for a in word),order=len(rows),profile=state,word=mask))
    seen={};digest=hashlib.sha256();checks=0
    frontier=[((),initial)]
    for d in range(9):
        nxt=[]
        for word,v in frontier:
            h,e,R=v
            require(decode(h,e)==word,('decode',word,v,decode(h,e)))
            require(v not in seen,('state collision',word,v));seen[v]=word
            require(admissible(h,e) and R==2*h-e-2 and h%4==3,('cone',word,v))
            require(math.gcd(*v)==1,('state content',word,v))
            ha=t.mv(cells[0]['M'],v)[0]
            haa=t.mv(cells[0]['M'],t.mv(cells[0]['M'],v))[0]
            hba=t.mv(cells[0]['M'],t.mv(cells[1]['M'],v))[0]
            require(haa-12*ha-13*h==48*R and hba-61*ha+28*h==64*R,('response identities',word))
            require(math.gcd(h,ha,haa,hba)==1,('observable content',word))
            digest.update(str((word,v)).encode());checks+=1
            if d<8:
                for a,z in enumerate(cells):
                    vv=t.mv(z['M'],v)
                    require(vv[0]>=[11,41,83][a]*h,('growth',word,a))
                    nxt.append((word+(a,),vv))
        frontier=nxt
    pairchecks=0;accepted=[]
    small={(h,e):w for (h,e,R),w in seen.items() if h<=501}
    for h in range(3,502,2):
        for e in range(1,h-1):
            got=decode(h,e);require(got==small.get((h,e)),('bounded pair membership',h,e,got));pairchecks+=1
            if got is not None:accepted.append(dict(h=h,e=e,cells=''.join('ABC'[a] for a in got)))
    print(json.dumps(dict(cells=cells,physical_checks=physical,word_decoder_checks=checks,state_data_sha256=digest.hexdigest(),admissible_pair_checks=pairchecks,accepted_pairs_through501=accepted),sort_keys=True,indent=2))

if __name__=='__main__':main()
