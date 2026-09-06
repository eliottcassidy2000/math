#!/usr/bin/env python3
"""Independent finite-grid audit: literal arcs, bounded knapsack, C++ sweep.

No repository producer imports. Every t=1..97096 is inspected in both modes.
The C++ accelerator uses exact integer arithmetic and always-active checks.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path
import json
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
GATES=0
DIGEST=sha256()


def gate(ok,label,data=None):
    global GATES
    if not ok:raise RuntimeError(f'{label}: {data!r}')
    GATES+=1;DIGEST.update((label+':'+repr(data)+'\n').encode())


def literal_components(p,q):
    denominator=14*p*q
    def arcs(u,factor):
        out=[(0,factor)]
        out.extend(((14*k-1)*factor,(14*k+1)*factor) for k in range(1,u))
        out.append((denominator-factor,denominator))
        return out
    aa,bb=arcs(p,q),arcs(q,p)
    i=j=0;pieces=[]
    while i<len(aa) and j<len(bb):
        left=max(aa[i][0],bb[j][0]);right=min(aa[i][1],bb[j][1])
        if left<right:pieces.append((left,right))
        if aa[i][1]<bb[j][1]:i+=1
        elif bb[j][1]<aa[i][1]:j+=1
        else:i+=1;j+=1
    lengths=[b-a for a,b in pieces]
    if pieces[0][0]==0 and pieces[-1][1]==denominator:
        lengths=[lengths[0]+lengths[-1]]+lengths[1:-1]
    return denominator,Counter(lengths)


def make_atlas():
    prime=[True]*357;prime[0]=prime[1]=False
    for p in range(2,19):
        if prime[p]:
            for n in range(p*p,357,p):prime[n]=False
    sums={1}
    for p in range(2,357):
        if prime[p] and p%3==2:
            sums={s*p**e for s in sums for e in range(3) if s*p**e<=356}
    sums.discard(1)
    sums.discard(2)  # No distinct positive pair has primitive sum two.
    gate(len(sums)==94 and 25 in sums and 125 not in sums,'multiplicatively-generated-strict-sums')
    rows=[]
    for total in sorted(sums):
        for p in range(1,(total+1)//2):
            q=total-p
            if gcd(p,q)!=1:continue
            denominator,lengths=literal_components(p,q)
            # Comparison is a control; the accelerator receives the literal
            # interval intersection data, not this resonance formula.
            formula=Counter({2*p:1})
            for k in range(1,(p+q-1)//14+1):formula[min(2*p,p+q-14*k)]+=2
            gate(lengths==formula,'literal-full-arc-multiset',(p,q))
            J=sum(lengths.values());M=sum(n*m for n,m in lengths.items())
            gate(J==2*((p+q-1)//14)+1,'strict-circle-component-count',(p,q))
            gate(J<=51 and 70*M>=denominator
                 and 304500*M>=(4313+37*J)*denominator,
                 'all-real-two-line-envelope-gates',(p,q))
            rows.append((p,q,J,denominator,tuple(sorted(lengths.items())),M))
    gate(len(rows)==5855,'complete-independent-atlas')
    # Preserve a different search order; this affects runtime, not the universe.
    rows.sort(key=lambda r:(-r[2],-r[3],r[0]))
    return rows


CPP=r'''
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using namespace std;
using ll=long long;
ll gates=0;
void need(bool ok,const char* label){++gates;if(!ok)throw runtime_error(label);}
ll up(ll a,ll b){return (a+b-1)/b;}
struct Row{int p,q,J,D,M;vector<pair<int,int>> lengths;};
int knapsack(const vector<int>& values,const array<int,91>& capacity,
             const array<int,91>& cost,int slots,int reserve1=0,int reserve2=0){
    const int neg=-1000000000;
    vector<int> best(slots+1,neg);best[0]=0;
    for(int d:values){
        int count=capacity[d]-(d==reserve1)-(d==reserve2);
        if(count<0)return neg;
        vector<int> next(slots+1,neg);
        for(int n=0;n<=slots;++n)if(best[n]!=neg)
            for(int use=0;use<=count && n+use<=slots;++use)
                next[n+use]=max(next[n+use],best[n]+use*cost[d]);
        best.swap(next);
    }
    return best[slots];
}
void emit(const string& name,const vector<int>& values){
    cout<<name<<" [";
    for(size_t i=0;i<values.size();++i){if(i)cout<<',';cout<<values[i];}
    cout<<"]\n";
}
int main(){
    int n;cin>>n;vector<int> allowed(n);for(int&d:allowed)cin>>d;
    int nr;cin>>nr;vector<Row> rows(nr);
    for(Row&r:rows){int k;cin>>r.p>>r.q>>r.J>>r.D>>r.M>>k;
        for(int i=0;i<k;++i){int a,b;cin>>a>>b;r.lengths.emplace_back(a,b);}}
    // Exact best slope at EACH component count. No imported two-line hull.
    map<int,pair<int,int>> slopes;
    for(const Row&r:rows){auto it=slopes.find(r.J);
        if(it==slopes.end() || ll(r.M)*it->second.second<ll(it->second.first)*r.D)
            slopes[r.J]={r.M,r.D};}
    const int caps[7]={90,30,9,4,2,1,1};
    for(int mode=0;mode<2;++mode){
        vector<int> domain=allowed;int emax=6;string name="profile6";
        if(mode){domain.clear();for(int d=1;d<=90;++d)domain.push_back(d);emax=30;name="scalar30";}
        vector<int> aggregate,components,coupled;
        for(int t=1;t<=97096;++t){
            vector<int> values;array<int,91> capacity{},cost{};
            for(int d:domain)if(t%d==0){
                values.push_back(d);for(int cap:caps)capacity[d]+=(d<=cap);
                cost[d]=7*d*up(t,7*d)-t;
                need(cost[d]==d*((7-(t/d)%7)%7),"exact cost identity");
            }
            int score=knapsack(values,capacity,cost,7);
            need(score>=0 && score%7==0,"seven-slot bounded knapsack");int E=score/7;
            vector<int> es;
            for(int e=1;e<=emax;++e)if(t%e==0){
                ll least=1000000000;
                for(auto item:slopes){int J=item.first,M=item.second.first,D=item.second.second;
                    least=min(least,ll(e)*(up(ll(t/e)*M,D)-J));}
                if(least<=E)es.push_back(e);
            }
            if(es.empty())continue;
            aggregate.push_back(t);bool has_component=false,has_coupled=false;
            map<pair<int,int>,int> conditioned;
            for(int e:es){
                int order=t/e;
                for(const Row&r:rows){
                    ll merged=ll(e)*(up(ll(order)*r.M,r.D)-r.J);
                    if(merged>E)continue;
                    ll credit=-ll(e)*r.J;
                    for(auto item:r.lengths)credit+=ll(e)*item.second*up(ll(order)*item.first,r.D);
                    need(credit>=merged,"component sum versus aggregate");
                    if(credit>E)continue;
                    has_component=true;
                    int dp=e*gcd(order,r.p),dq=e*gcd(order,r.q);
                    need(gcd(dp,dq)==e,"forced pair common gcd");
                    if(dp>90 || dq>90 || !capacity[dp] || !capacity[dq])continue;
                    if(dp>dq)swap(dp,dq);
                    auto key=make_pair(dp,dq);auto it=conditioned.find(key);int EE;
                    if(it==conditioned.end()){
                        int rest=knapsack(values,capacity,cost,5,dp,dq);
                        if(rest<0)EE=-1;
                        else{int s=rest+cost[dp]+cost[dq];need(s%7==0,"forced-pair seven-slot integrality");EE=s/7;}
                        conditioned[key]=EE;
                    }else EE=it->second;
                    need(EE<=E,"conditioning lowers bag maximum");
                    if(EE>=0 && credit<=EE){has_coupled=true;break;}
                }
                if(has_coupled)break;
            }
            if(has_component)components.push_back(t);
            if(has_coupled)coupled.push_back(t);
            need(!has_coupled || has_component,"stage nesting");
        }
        emit(name+"_aggregate",aggregate);emit(name+"_components",components);emit(name+"_coupled",coupled);
    }
    cout<<"CPP_GATES "<<gates<<'\n';
}
'''


def component_credit(t,e,row):
    p,q,J,D,lengths,M=row
    return e*(sum(mult*((t//e*n+D-1)//D) for n,mult in lengths)-J)


def main():
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    gate(sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','pinned-inherited-profiles')
    levels=json.loads(raw)['levels']
    allowed=sorted({c for c,w in levels['6']['profiles']})
    gate(allowed==levels['6']['gcds'] and len(allowed)==42,'all42-exact-core-values')
    rows=make_atlas();lookup={r[:2]:r for r in rows}
    gate(F(lookup[(1,10)][5],lookup[(1,10)][3])==F(1,70)
         and lookup[(1,10)][2]==1
         and F(lookup[(5,348)][5],lookup[(5,348)][3])==F(62,3045)
         and lookup[(5,348)][2]==51,'both-envelope-lines-attained')
    gate(literal_components(1,2)==(28,Counter({2:1})),'contained-origin-component-hostile')
    gate(component_credit(74550,30,lookup[(1,355)])==0,'open-component-one-grid-spacing')
    gate(component_credit(27360,6,lookup[(5,348)])==294,'individual-ceilings-retain-extra-credit')
    gate(component_credit(23760,6,lookup[(25,294)])==246,'last-profile-survivor-component-credit')
    gate(component_credit(23760,6,lookup[(25,294)])<min(component_credit(23760,6,lookup[p]) for p in ((1,10),(5,348))),
         'two-measure-extremizers-do-not-supply-component-minimum')
    aggregate=lambda t,e:min(e*((t*r[5]//e+r[3]-1)//r[3]-r[2]) for r in rows)
    gate([aggregate(57995,e) for e in (1,5,7)]==[828,825,826],'all-e-needed-rounding-hostile')
    data=[str(len(allowed)),' '.join(map(str,allowed)),str(len(rows))]
    for p,q,J,D,lengths,M in rows:
        data.append(' '.join(map(str,(p,q,J,D,M,len(lengths)))))
        data.extend(f'{n} {mult}' for n,mult in lengths)
    with tempfile.TemporaryDirectory(prefix='third-grid-audit-') as directory:
        directory=Path(directory);source=directory/'engine.cpp';binary=directory/'engine'
        source.write_text(CPP)
        subprocess.run(['g++','-std=c++17','-O2',str(source),'-o',str(binary)],check=True,capture_output=True,text=True)
        result=subprocess.run([str(binary)],input='\n'.join(data)+'\n',capture_output=True,text=True,check=True)
    recovered={};cpp_gates=0
    for line in result.stdout.splitlines():
        label,payload=line.split(' ',1)
        if label=='CPP_GATES':cpp_gates=int(payload)
        else:recovered[label]=json.loads(payload)
    producer={}
    for line in (ROOT/'05-knowledge/results/third_20260906_grid_bootstrap.out').read_text().splitlines():
        if line.startswith('SCALE_SET '):
            item=json.loads(line[10:]);producer[item['name']]=item['survivors']
    expected={'profile6_aggregate':(9498,27360),'profile6_components':(8308,23760),'profile6_coupled':(8301,23760),
              'scalar30_aggregate':(34532,88920),'scalar30_components':(32294,74550),'scalar30_coupled':(32272,74520)}
    for name in sorted(expected):
        values=recovered[name]
        gate(values==sorted(set(values)) and min(values)==1,'canonical-complete-clock-set',name)
        gate((len(values),max(values))==expected[name],'exact-count-and-bound',name)
        gate(values==producer[name],'independent-full-survivor-array-comparison',name)
        print('AUDIT_SCALE_SET',json.dumps({'name':name,'count':len(values),'maximum':max(values),
            'sha256':sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest(),'survivors':values},separators=(',',':')))
    exclusions=sorted(set(recovered['profile6_components'])-set(recovered['profile6_coupled']))
    gate(exclusions==[12425,14872,14910,15390,15504,20520,21240],'seven-pair-reservation-exclusions')
    print('PAIR_RESERVATION_EXCLUSIONS',exclusions)
    print('CPP_EXACT_GATES',cpp_gates)
    print('PYTHON_EXACT_GATES',GATES)
    print('SEMANTIC_SHA256',DIGEST.hexdigest())
    print('PASS')


if __name__=='__main__':main()
