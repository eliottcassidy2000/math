#include <bits/stdc++.h>
using namespace std;using L=long long;
// Key digits1..4 encode the successive endpoint-color pairs00,01,10,11.
map<unsigned,L> tensor(vector<unsigned>a,array<unsigned,2> masks){
 int n=a.size();vector<int>p(n);iota(p.begin(),p.end(),0);map<unsigned,L>out;
 do{
  vector<pair<bool,vector<unsigned>>>gaps;
  for(int g=0;g<=n;g++){
   int left=g?p[g-1]:-1,right=g<n?p[g]:-1;
   bool e=left<0||right<0||((a[left]>>right)&1);
   vector<unsigned>choices;
   for(int f=0;f<2;f++)for(int l=0;l<2;l++)
    if((left<0||!((masks[f]>>left)&1))&&(right<0||((masks[l]>>right)&1)))choices.push_back(1+2*f+l);
   gaps.push_back({e,choices});
  }
  function<void(int,unsigned)>rec=[&](int g,unsigned key){
   if(g==n+1){out[key]++;return;}
   if(gaps[g].first)rec(g+1,key);
   for(unsigned c:gaps[g].second)rec(g+1,key*5+c);
  };rec(0,0);
 }while(next_permutation(p.begin(),p.end()));
 return out;
}
string word(unsigned key){if(!key)return"empty";vector<string>v;while(key){int x=key%5-1;v.push_back(string(1,char('0'+(x>>1)))+char('0'+(x&1)));key/=5;}reverse(v.begin(),v.end());string s;for(auto x:v)s+=x+",";s.pop_back();return s;}
void run(string name,vector<unsigned>p,vector<unsigned>c,array<unsigned,2>pm,array<unsigned,2>cm){
 auto P=tensor(p,pm),C=tensor(c,cm);map<unsigned,L>D=C;
 for(auto [k,v]:P)D[k]-=2*v;D[0]--;for(int x=1;x<=4;x++)D[x]--;
 int pos=0,neg=0;cout<<"FAMILY "<<name<<" Pterms "<<P.size()<<" Cterms "<<C.size()<<endl;
 for(auto [k,v]:D)if(v){if(v>0)pos++;else neg++;cout<<k<<" "<<word(k)<<" "<<v<<endl;}
 cout<<"SIGNS "<<pos<<" "<<neg<<endl;
}
int main(){
 run("A_15_63",{2,24,3,37,13,23},{2,8,67,101,15,87,19},{15,63},{79,127});
 run("B_0_5",{32,5,1,39,47,70,31},{0,5,129,39,47,199,159,27},{0,5},{0,5});
}
