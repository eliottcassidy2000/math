#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace std;using R=vector<unsigned>;using Z=long long;
void need(bool x,const char*s){if(!x)throw runtime_error(s);}int bit(unsigned x){return __builtin_ctz(x);}
vector<int>digits(int n){vector<int>b;while(n){b.push_back(n&1);n>>=1;}reverse(b.begin(),b.end());return b;}
Z g(int u){auto b=digits(u-1);int n=b.size();Z out=0;for(int s=0;s<(1<<n);s++){int first=-1,last=-1;for(int v=0;v<n;v++)if(s>>v&1)last=v;else if(first<0)first=v;out+=(last<0||!b[last])&&(first<0||b[first]);}return out;}
void birth(R&r,unsigned beaten){int n=r.size();r.push_back(beaten);for(int u=0;u<n;u++)if(!(beaten>>u&1))r[u]|=1<<n;}
R fromword(vector<int>const&w){R r{2,0};for(int c:w)birth(r,c);return r;}
vector<int>prefix(int u,bool B){vector<int>w;for(int b:digits(u-1))w.push_back(b?(B?2:1):3);return w;}
vector<int>qword(int u,int v){auto w=prefix(u,true);w.push_back(1);w.push_back(3);auto b=digits(v-1);for(auto it=b.rbegin();it!=b.rend();it++)w.push_back(2*(*it));return w;}
R directQ(int u,int v){R r{0};for(int b:digits(u-1))birth(r,b?0:1);int z=r.size();birth(r,1);birth(r,1);int root=r.size();birth(r,1u<<z);auto bs=digits(v-1);for(auto it=bs.rbegin();it!=bs.rend();it++)birth(r,*it?1u<<root:0);return r;}
vector<vector<Z>> endpoint(R const&r,bool reverse=false){int n=r.size(),N=1<<n;vector<vector<Z>>d(N,vector<Z>(n));for(int v=0;v<n;v++)d[1<<v][v]=1;for(int s=1;s<N;s++)for(int v=0;v<n;v++)if(d[s][v])for(int u=0;u<n;u++)if(!(s>>u&1)&&(reverse?(r[u]>>v&1):(r[v]>>u&1)))d[s|1<<u][u]+=d[s][v];return d;}
Z H(R const&r){auto d=endpoint(r);return accumulate(d.back().begin(),d.back().end(),Z(0));}
array<Z,2>literal(R const&r){int n=r.size();unsigned full=(1u<<n)-1;Z h=0,c=0;function<void(int,int,unsigned)>go=[&](int first,int last,unsigned used){if(used==full){h++;c+=first==0&&(r[last]&1);return;}unsigned q=r[last]&~used;while(q){unsigned b=q&-q;q-=b;go(first,bit(b),used|b);}};for(int v=0;v<n;v++)go(v,v,1<<v);return{h,c};}
array<Z,5>profile(R const&r,int root){auto e=endpoint(r),s=endpoint(r,true);int F=(1<<r.size())-1;vector<Z>h(F+1);h[0]=1;for(int M=1;M<=F;M++)h[M]=accumulate(e[M].begin(),e[M].end(),Z(0));Z Rv=0;for(int M=0;M<=F;M++)Rv+=s[M][root]*h[F^M];return{h[F],e[F][root],Rv,h[F^(1<<root)],s[F][root]};}
bool transitive(R const&r,unsigned del){int n=r.size();for(int a=0;a<n;a++)if(!(del>>a&1))for(int b=a+1;b<n;b++)if(!(del>>b&1))for(int c=b+1;c<n;c++)if(!(del>>c&1)){bool x=r[a]>>b&1,y=r[b]>>c&1,z=r[c]>>a&1;if((x&&y&&z)||(!x&&!y&&!z))return false;}return true;}
Z scalar(char f,int u,int v,Z x,Z y){int p=u%2,t=v%2;if(f=='N')return x*y+(u-p)*y+(v-t)*x+3LL*u*v-u-v-p-t+p*t;if(f=='M')return x*y+2*u*y+2*v*x+6LL*u*v-u-v-u*t-v*p-p-t;if(f=='P')return x*y+4*v*x+4*u*y+6LL*u*v-u*t-v*p;return y*(x+3*u-p-2)-t*(u-1)+v*(2*x+4*u-p-1);}
int main(int argc,char**argv){need(argc==2,"pair output filename");int maps=0,arcs=0,seeds=0,graphs=0,terminals=0;for(int u=2;u<=128;u++)for(int v=1;v<=8;v++){auto A=directQ(u,v),B=fromword(qword(u,v));int m=digits(u-1).size(),n=A.size();vector<int>f(n);f[0]=0;for(int x=1;x<=m;x++)f[x]=x+1;f[m+1]=m+2;f[m+2]=m+3;f[m+3]=1;for(int x=m+4;x<n;x++)f[x]=x;for(int x=0;x<n;x++)for(int y=x+1;y<n;y++){need((A[x]>>y&1)==(B[f[x]]>>f[y]&1),"independent Q graph map");arcs++;}need(f[m+3]==1,"Q new-root identity");maps++;}
for(int u=2;u<=64;u++){Z D=g(u),p=u%2;auto T=fromword(qword(u,1));need(profile(T,1)==array<Z,5>{3*D+6*u-2*p-2,D+2*u-p-1,3*D+7*u-2*p-3,D+3*u-p-2,2},"independent Q exact profile");seeds++;}
for(int u=2;u<=12;u++)for(int v=1;v<=12;v++){auto r=fromword(qword(u,v));Z h=scalar('Q',u,v,g(u),g(v));auto lc=literal(r);need(H(r)==h&&lc[0]==h&&lc[1]==2,"Q literal paths and cycles");need(transitive(r,3)&&h>=2LL*u+v*(4LL*u+1),"Q physical inverse bound");graphs++;}
for(int u=2;u<=24;u++)for(int k=1;k<=5;k++){auto A=prefix(u,false),B=prefix(u,true);B.push_back(1);for(int t=0;t<k;t++){A.push_back(0);B.push_back(0);}A.push_back(3);B.push_back(3);for(int z=0;z<2;z++){auto r=fromword(z?B:A);Z h=scalar(z?'M':'N',u,1<<k,g(u),1);auto lc=literal(r);need(H(r)==h&&lc[0]==h&&lc[1]==2,"unary terminal physical reduction");terminals++;}}
R boundary{2,4,9,19,39,335,159,319,607,1279,511};need(H(boundary)==533&&literal(boundary)==array<Z,2>{533,2},"533 literal graph");int pairs=0;for(int a=0;a<11;a++)for(int b=a+1;b<11;b++)pairs+=transitive(boundary,(1<<a)|(1<<b));need(!pairs&&transitive(boundary,(1<<2)|(1<<5)|(1<<8)),"533 exact feedback number three");
const int cap=1861;vector<Z>gg(cap+1);for(int u=1;u<=cap;u++)gg[u]=g(u);ofstream out(argv[1]);out<<"15,R,0,0\n";long long found=1,raw=0;array<int,3>gaps{};for(int u=1;u<=cap;u++)for(int v=1;v<=cap;v++)for(char f:{'N','M','P','Q'}){if((f=='N'||f=='M')&&(u<2||v<2||u>v))continue;if(f=='P'&&u>v)continue;if(f=='Q'&&u<2)continue;raw++;Z h=scalar(f,u,v,gg[u],gg[v]);if(h<=cap){out<<h<<','<<f<<','<<u<<','<<v<<'\n';found++;for(int z=0;z<3;z++)if(h==array<int,3>{533,821,1861}[z])gaps[z]++;}}need(gaps==array<int,3>{0,0,0},"independent rectangular exclusions");need(bool(out),"inverse pair file");cout<<"{\"status\":\"PASS\",\"Q_marked_maps\":"<<maps<<",\"map_arc_checks\":"<<arcs<<",\"Q_seed_profiles\":"<<seeds<<",\"Q_literal_graphs\":"<<graphs<<",\"terminal_literal_graphs\":"<<terminals<<",\"inverse_rectangle_bound\":"<<cap<<",\"raw_parameter_pairs\":"<<raw<<",\"inverse_hits\":"<<found<<",\"boundary_H\":533,\"boundary_C\":2,\"feedback_number\":3}\n";}
