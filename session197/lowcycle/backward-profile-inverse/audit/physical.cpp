#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>
#include <map>
#include <array>
#include <stdexcept>
using namespace std;using ll=long long;using G=vector<int>;
void need(bool x,const char*s){if(!x)throw runtime_error(s);}struct F{ll h=0;vector<ll>e,s,R;};
F literal(const G&a){int n=a.size();F f;f.e.resize(n);f.s.resize(n);f.R.resize(n);if(!n){f.h=1;return f;}vector<int>p(n);iota(p.begin(),p.end(),0);do{vector<int>bad(n);for(int i=1;i<n;i++)bad[i]=bad[i-1]+!(a[p[i-1]]>>p[i]&1);if(!bad.back()){f.h++;f.s[p[0]]++;f.e[p.back()]++;}for(int cut=0;cut<n;cut++)if((cut==0||!bad[cut-1])&&bad.back()==bad[cut])f.R[p[cut]]++;}while(next_permutation(p.begin(),p.end()));return f;}
G omit(const G&a,int r){G b(a.size()-1);for(int i=0;i<(int)a.size();i++)if(i!=r)for(int j=0;j<(int)a.size();j++)if(j!=r&&(a[i]>>j&1))b[i-(i>r)]|=1<<(j-(j>r));return b;}
vector<ll>line(string s){replace(s.begin(),s.end(),',',' ');stringstream z(s);vector<ll>v;ll x;while(z>>x)v.push_back(x);return v;}
int main(int argc,char**argv){need(argc>=3,"bank and output");map<array<ll,3>,ll>states;states[{1,1,1}]=1;states[{1,0,2}]=1;long graphs=0,roots=0;auto check=[&](const G&a){F p=literal(a);graphs++;if(p.h>101)return;for(int r=0;r<(int)a.size();r++)if(p.s[r]==1){ll d=literal(omit(a,r)).h;array<ll,3>key={p.h,p.e[r],p.R[r]};if(states.count(key))need(states[key]==d,"profile deletion uniqueness");states[key]=d;roots++;}};check({2,4,1});ifstream in(argv[1]);string s;while(getline(in,s)){auto v=line(s);if(v[0]>9||v[1]!=1)continue;G a(v.begin()+2,v.end());check(a);}ofstream out(argv[2]);for(auto [t,d]:states)out<<t[0]<<","<<t[1]<<","<<t[2]<<","<<d<<"\n";cout<<"PASS canonical_C1_graphs "<<graphs<<" physical_start_one_roots "<<roots<<" distinct_profiles_with_bases "<<states.size()<<" cap 101\n";
if(argc>=4){ifstream rec(argv[3]);int checks=0;while(getline(rec,s)){auto v=line(s);int n=v[0],r=v[1];G a(v.begin()+6,v.end());need((int)a.size()==n,"recipe order");F p=literal(a);need(p.h==v[2]&&p.e[r]==v[3]&&p.R[r]==v[4]&&p.s[r]==1&&literal(omit(a,r)).h==v[5],"literal reconstructed recipe");checks++;}cout<<"PASS reconstructed_recipes "<<checks<<"\n";}
}
