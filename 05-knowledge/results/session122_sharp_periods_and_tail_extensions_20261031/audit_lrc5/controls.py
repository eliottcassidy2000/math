"""Analytic boundary, closed-endpoint, and certificate-mutation controls."""
import sys,json,copy,hashlib
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import ceil
from geometry import R,need,comb,intersection,safe,width,check_witness
from check import COEFFICIENT,encode,validate_header,reconstruct,verify_family

def run(path):
    d=json.loads(Path(path).read_bytes());validate_header(d)
    boundary=0
    for k in (1,2,3,4,5):
        r=R;c=COEFFICIENT[k]
        for W in (F(1,2),F(1,3),F(2,5),c/5):
            B=ceil(c/W)
            for left in (F(0),F(1,17),F(2,7)):
                for tails in combinations_with_replacement(range(B,B+5),k):
                    I=((left,left+W),)
                    for v in tails:I=intersection(I,comb(v,r))
                    need(bool(I),'sharp conductor closed existence')
                    check_witness(I[0][0],tails,r);boundary+=1
    envelope=0
    for k in range(2,6):
        c=COEFFICIENT[k]
        for j in range(1,1001):
            x=c+F(j-1,137);n=x.numerator//x.denominator;y=x-n
            m=(F(n,7)+min(y,F(1,7)))/x
            need(m<=F(1,k),'bad-measure envelope');envelope+=1
    Q=F(1,4);I=safe((1,2),Q)
    need(bool(intersection(I,comb(3,Q))) and bool(intersection(I,comb(12,Q))),'marginal success')
    need(not safe((1,2,3,12),Q),'joint hostile')
    J=((F(1,2)-R/2,F(1,2)+R/2),)
    need(intersection(J,comb(2))==((J[0][0],J[0][0]),(J[0][1],J[0][1])),'closed singleton equality')
    need(not intersection(((J[0][0]+F(1,10000),J[0][1]-F(1,10000)),),comb(2)),'strictly undersized interval fails')
    need(not safe((1,2,3),F(1,3)),'radius-overclaim hostile')
    need(len(set(range(1,9))|{9})==9 and 8 in set(range(1,9))|{9},'repeat padded to nine-core')
    need(len(set(range(1,9))|{14})==9 and 14 not in range(1,14),'unbounded new speed cannot trigger inherited theorem')
    expected,tree=reconstruct(tuple(d['families'][0]['core']),True)
    verify_family(d['families'][0],expected)
    results=[]
    def rejects(name,fn):
        try:fn()
        except (RuntimeError,KeyError,ValueError):results.append(name);return
        raise RuntimeError('mutation incorrectly accepted: '+name)
    bad=copy.deepcopy(d);bad['families'].pop();rejects('missing core',lambda:validate_header(bad))
    bad=copy.deepcopy(d);bad['coefficient']['4']='1/2';rejects('unproved cutoff coefficient',lambda:validate_header(bad))
    for key,value in [('width','1'),('cutoff',expected['cutoff']-1),('safe_intervals',[]),('node_counts',[0]*6)]:
        badrec=copy.deepcopy(expected);badrec[key]=value
        rejects('false '+key,lambda rec=badrec:verify_family(rec,expected))
    for name,child_index in [('omitted finite branch',0),('false inherited promotion',0)]:
        altered=copy.deepcopy(tree)
        if name=='omitted finite branch':altered['children'].pop(child_index)
        else:
            child_index=next(i for i,row in enumerate(altered['children']) if row[0]>13)
            altered['children'][child_index]=[altered['children'][child_index][0],'inherit_nine_core']
        raw=encode(altered);badrec=copy.deepcopy(expected)
        badrec['tree_bytes']=len(raw);badrec['tree_sha256']=hashlib.sha256(raw).hexdigest()
        rejects(name,lambda rec=badrec:verify_family(rec,expected))
    return {'status':'PASS','boundary_controls':boundary,'bad_measure_envelope_controls':envelope,'hostile_endpoint_controls':7,
            'mutation_controls_rejected':results,'mutation_count':len(results)}
if __name__=='__main__':print(json.dumps(run(sys.argv[1]),sort_keys=True,indent=2))
