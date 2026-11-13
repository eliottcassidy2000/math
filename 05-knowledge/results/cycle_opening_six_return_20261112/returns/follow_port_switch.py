from probe import *

def child(P,a):
    n=len(P);full=(1<<n)-1
    return tuple(row|((1<<n) if i==a else 0) for i,row in enumerate(P))+(full^(1<<a),)

def main():
    old=json.loads((HERE/'port_collisions.json').read_text())
    seeds={}
    for hit in old['hits']:
        for v in hit['ports']:
            seeds[(tuple(hit['parent']),v)]=hit
    successful=[];summary=Counter();allrows=[]
    for (P,a),base in seeds.items():
        T=child(P,a);r=len(P)
        X=child(T,r);x=len(T)
        parenth,starts,e2,values=profiles(X)
        target=values[x];olddefect=starts[x]
        matches=[i for i in range(len(X)) if i!=x and values[i]==target]
        improved=[i for i in matches if starts[i]<olddefect]
        ready=[i for i in matches if starts[i] in (1,2)]
        summary[(len(P),'seeds')]+=1
        summary[(len(P),'neutral')]+=bool(matches)
        summary[(len(P),'improved')]+=bool(improved)
        summary[(len(P),'ready')]+=bool(ready)
        if matches:
            row=dict(parent=P,port=a,H=base['H'],defect=profiles(P)[1][a],advanced_parent=X,target=target,raw_defect=olddefect,neutral_ports=[dict(port=i,defect=starts[i]) for i in matches],ready=ready,improved=improved)
            allrows.append(row)
            if improved:successful.append(row)
    out=dict(scope='all distinct rootchildren participatingin unequaldefect neutralpairs fromparentbank1..8;exacttwo births thenallneutralentryswitches',summary={str(k):v for k,v in summary.items()},matches=allrows,improvements=successful)
    (HERE/'port_advance.json').write_text(json.dumps(out,indent=2)+'\n')
    print('SUMMARY',out['summary']);print('IMPROVEMENTS',len(successful))
    for row in successful[:20]: print({k:v for k,v in row.items() if k!='advanced_parent'})
if __name__=='__main__':main()
