#include <bits/stdc++.h>
using namespace std;using I=long long;using P=array<I,9>;using Key=array<I,9>;struct Rec{int source;string word;};
P mv(P p,char c){auto[h,s,e,L,R,d,C,b,T]=p;if(c=='O')return{h,s,0,L,2*R,d,2*C,2*b,3*T};if(c=='Z')return{h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b};return{h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C};}
int main(int argc,char**argv){ifstream in(argv[1]);int D=argc>2?atoi(argv[2]):7,ix,mask;P p;vector<P>src;vector<pair<int,int>>loc;while(in>>ix>>mask){for(auto &v:p)in>>v;src.push_back(p);loc.push_back({ix,mask});}
map<Key,array<Rec,2>>bank;I visited=0,hits=0;
for(int i=0;i<(int)src.size();i++){
function<void(P,string)>walk=[&](P p,string word){auto[h,s,e,L,R,d,C,b,T]=p;visited++;Key key{h,e,L,R,d,C,b,T,s+2*(s%3)};auto it=bank.find(key);if(it==bank.end())it=bank.emplace(key,array<Rec,2>{Rec{-1,""},Rec{-1,""}}).first;auto&z=it->second;z[s%3-1]={i,word};if(z[0].source>=0&&z[1].source>=0){hits++;cout<<"MATCH";for(auto v:key)cout<<" "<<v;cout<<"\n";for(auto r:z)cout<<loc[r.source].first<<" "<<loc[r.source].second<<" "<<(r.word.empty()?"-":r.word)<<"\n";cout<<flush;}
if((int)word.size()==D)return;for(char c:string("OZJ"))walk(mv(p,c),word+c);};walk(src[i],"");}
cerr<<"sources "<<src.size()<<" depth "<<D<<" visited "<<visited<<" keys "<<bank.size()<<" hits "<<hits<<"\n";
}
