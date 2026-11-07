#include<bits/stdc++.h>
using namespace std;using I=long long;using P=array<I,9>;using Key=array<I,3>;struct Rec{int i,j;};
P mv(P p,char c,bool mod){auto[h,s,e,L,R,d,C,b,T]=p;P q;if(c=='O')q={h,s,0,L,2*R,d,2*C,2*b,3*T};else if(c=='Z')q={h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b};else q={h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C};if(mod)for(auto&x:q){x%=1000000007;if(x<0)x+=1000000007;}return q;}
I X(P p){return p[0]+p[3]-(p[1]+2*(p[1]%3))/3;}
int main(int argc,char**argv){int D=argc>2?atoi(argv[2]):5;if(D>5)throw runtime_error("exact prefix bound");ifstream in(argv[1]);int ix,mask,n;P p;vector<P>src;map<I,array<vector<int>,2>>groups;while(in>>ix>>mask>>n){for(auto&x:p)in>>x;int i=src.size();src.push_back(p);groups[X(p)][p[1]%3-1].push_back(i);}vector<string>w;function<void(string)>gen=[&](string s){if(s.size())w.push_back(s);if((int)s.size()==D)return;for(char c:string("OZJ"))gen(s+c);};gen("");vector<bool>has;for(auto s:w)has.push_back(s.find('J')!=string::npos);I cases=0,pairedkeys=0,tests=0,hits=0,basecount=0;
for(auto &[base,ss]:groups){if(ss[0].empty()||ss[1].empty())continue;basecount++;map<Key,array<vector<Rec>,2>>bank;
for(int r=0;r<2;r++)for(int i:ss[r])for(int j=0;j<(int)w.size();j++){P p=src[i];Key key;for(int k=0;k<3;k++){for(char c:w[j])p=mv(p,c,false);key[k]=X(p);}bank[key][r].push_back({i,j});cases++;}
for(auto &[key,rr]:bank){if(rr[0].empty()||rr[1].empty())continue;bool anyj=false;for(auto r:rr)for(auto z:r)anyj|=has[z.j];if(!anyj)continue;pairedkeys++;
map<array<I,19>,Rec>seen;
for(int r=0;r<2;r++)for(auto z:rr[r]){P p=src[z.i];I c=(p[1]+2*(p[1]%3))/3;array<I,19>sig;for(int k=0;k<19;k++){sig[k]=(p[0]+p[3]-c+1000000007)%1000000007;for(char ch:w[z.j])p=mv(p,ch,true);}tests++;
if(r==0){auto it=seen.find(sig);if(it==seen.end()||has[z.j])seen[sig]=z;}else{auto it=seen.find(sig);if(it!=seen.end()&&(has[z.j]||has[it->second.j])){hits++;cout<<"MATCH "<<base<<"\n"<<it->second.i<<" "<<w[it->second.j]<<"\n"<<z.i<<" "<<w[z.j]<<"\n"<<flush;}}}}
if(basecount%500==0)cerr<<"bases "<<basecount<<" cases "<<cases<<" hits "<<hits<<"\n";
}cerr<<"FINAL bases "<<basecount<<" macros "<<w.size()<<" cases "<<cases<<" three_step_keys "<<pairedkeys<<" full_tests "<<tests<<" hits "<<hits<<"\n";
}
