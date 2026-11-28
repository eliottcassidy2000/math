#include <bits/stdc++.h>
using namespace std;using U=uint64_t;
constexpr int N=56;
struct O{int h,y;array<int,3>z;};
vector<O> opts;vector<int> slots[3][56];
U allow[3]; vector<int> sel,saved;
mt19937 gen(33);
int P[N][N],Q[N][N],locked[N][N],ip[N][N],iq[N][N];
U sw[N][N][N];
vector<array<int,3>> moves;
int countLeaf=0,countDbl=0,counts[60]={};
bool recolor(){
 countLeaf++;
 vector<int> C[N], ds,holes;
 for(int i:sel)C[opts[i].h].push_back(opts[i].y);
 for(int h=3;h<56;h++)if(h!=28){
   if(C[h].size()>2)return false;
   if(C[h].size()==2)ds.push_back(h);
   if(C[h].empty())holes.push_back(h);
 }
 assert(ds.size()==holes.size());
 countDbl++;
 counts[ds.size()]++;
 vector<pair<int,int>> adj[N];
 for(int h:ds)for(int k:holes){
   auto [a,b]=make_pair(C[h][0],C[h][1]);
   if((sw[h][k][a]>>b&1)==0 && sw[h][k][a]!=0)adj[h].emplace_back(k,a);
   else if((sw[h][k][b]>>a&1)==0 && sw[h][k][b]!=0)adj[h].emplace_back(k,b);
 }
 U seen;int owner[N],w[N];fill(owner,owner+N,-1);
 auto m=[&](auto&&m,int a)->bool{
   for(auto [b,y]:adj[a])if(!(seen>>b&1)){
    seen|=1ULL<<b;
    if(owner[b]<0||m(m,owner[b])){owner[b]=a;w[b]=y;return true;}
   }
   return false;
 };
 for(int a:ds){seen=0;if(!m(m,a))return false;}
 moves.clear();
 for(int b:holes)moves.push_back({owner[b],b,w[b]});
 saved=sel;
 return true;
}
int global_best=53;
bool greedy(){
  int p[N][N],q[N][N],pi[N][N],qi[N][N],m[N][N],d[N];
  U occ[N],ma[N][15];
  vector<int> br;
  for(int h=3;h<N;h++)if(h!=28)br.push_back(h);
  auto fac=[&](int h,int k,int y) {
    U set=1ULL<<y,open=set;
    while(open){
      int v=__builtin_ctzll(open);open&=open-1;
      U add=((1ULL<<pi[k][p[h][v]])|(1ULL<<qi[k][q[h][v]]))&~set;
      set|=add;open|=add;
    }
    return set;
  };
  vector<array<int,3>> track,optsv;
  auto swap_c=[&](int h,int k,int y,U c) {
    track.push_back({h,k,y});
    U dif=(occ[h]^occ[k])&c;occ[h]^=dif;occ[k]^=dif;
    d[h]=__builtin_popcountll(occ[h]);d[k]=__builtin_popcountll(occ[k]);
    for(int i=1;i<=14;i++){
      U df=(ma[h][i]^ma[k][i])&c;ma[h][i]^=df;ma[k][i]^=df;
    }
    for(;c;c&=c-1){
      int v=__builtin_ctzll(c);
      swap(p[h][v],p[k][v]);swap(q[h][v],q[k][v]);swap(m[h][v],m[k][v]);
      pi[h][p[h][v]]=pi[k][p[k][v]]=v;
      qi[h][q[h][v]]=qi[k][q[k][v]]=v;
    }
  };
  for(int rep=0;rep<12;rep++){
    memcpy(p,P,sizeof p);memcpy(q,Q,sizeof q);memcpy(pi,ip,sizeof pi);memcpy(qi,iq,sizeof qi);
    memcpy(m,locked,sizeof m);memset(occ,0,sizeof occ);memset(ma,0,sizeof ma);
    for(int i:sel)occ[opts[i].h]|=1ULL<<opts[i].y;
    for(int h:br){
      d[h]=__builtin_popcountll(occ[h]);
      for(int y=0;y<N;y++)if(m[h][y])ma[h][m[h][y]]|=1ULL<<y;
    }
    track.clear();
    int best=52,last=0;
    for(int step=0;step<300;step++){
      int surplus=0;for(int h:br)surplus+=max(0,d[h]-1);
      if(surplus<global_best){
        global_best=surplus;
        cerr<<"best remaining "<<surplus<<" with "<<track.size()<<" moves\n";
        ofstream out("greedy_best.sol");out<<track.size()<<"\n";
        for(auto [a,b,y]:track)out<<a<<" "<<b<<" "<<y<<"\n";
        for(int i:sel)out<<opts[i].h<<" "<<opts[i].y<<"\n";
      }
      if(!surplus){moves=track;saved=sel;return true;}
      if(surplus<best){best=surplus;last=step;}
      if(step-last>60)break;
      shuffle(br.begin(),br.end(),gen);
      int weightBest=1<<30;
      U good=0; array<int,3> move{};
      int randomWalk=gen()%100;
      for(int h:br)for(int k:br)if(h!=k && d[h]>0 && (d[h]>1 || d[k]==0 || randomWalk<24)){
        if(d[h]<=d[k] && randomWalk>=24)continue;
        U qs=occ[h]&~occ[k];
        if(!qs)continue;
        while(qs){
          int y=__builtin_ctzll(qs);
          U c=fac(h,k,y);qs&=~c;
          int a=__builtin_popcountll(c&occ[h]),b=__builtin_popcountll(c&occ[k]);
          int nh=d[h]+b-a,nk=d[k]+a-b;
          int cost=(nh==0)+(nk==0)-(d[h]==0)-(d[k]==0);
          if(cost>0 && randomWalk>=10)continue;
          bool ok=true;
          for(int i=1;i<=14;i++)if(bool(ma[h][i]&c)!=bool(ma[k][i]&c)){ok=false;break;}
          if(!ok)continue;
          int weight=cost*1024+(gen()%1024);
          if(weight<weightBest){weightBest=weight;good=c;move={h,k,y};}
        }
      }
      if(!good)break;
      swap_c(move[0],move[1],move[2],good);
    }
  }
  return false;
}

int refid[N], owner3[3][N], budget;
unsigned long long nodes=0;
U bmask[1741];
bool anchored(U open,U held,U removed){
 ++nodes;
 if(!open){
   for(U z=held;z;z&=z-1)sel.push_back(refid[__builtin_ctzll(z)]);
   bool ok=recolor() || greedy();
   sel.resize(sel.size()-__builtin_popcountll(held));
   return ok;
 }
 int left=budget-__builtin_popcountll(removed); if(left<0)return false;
 int best=1<<20,j=-1;vector<int> choices;
 for(U z=open;z;z&=z-1){
   int col=__builtin_ctzll(z);
   vector<int> can;
   for(int i:slots[0][col]){
      auto &a=opts[i].z;
      if((allow[1]>>a[1]&allow[2]>>a[2]&1) && __builtin_popcountll(bmask[i]&held)<=left)
         can.push_back(i);
      if((int)can.size()>=best)break;
   }
   if(can.empty())return false;
   if((int)can.size()<best){best=can.size();choices=move(can);j=col;}
 }
 for(int i:choices){
   auto &a=opts[i].z;U rel=bmask[i]&held;
   sel.push_back(i);
   for(int k=1;k<3;k++)allow[k]^=1ULL<<a[k];
   if(anchored((open^(1ULL<<j))|rel,held^rel,removed|rel))return true;
   for(int k=1;k<3;k++)allow[k]^=1ULL<<a[k];
   sel.pop_back();
 }
 return false;
}
int main(int argc,char**argv){
 if(argc>1)gen.seed(atoi(argv[1]));
 int tries=argc>2?atoi(argv[2]):10000;
 ifstream f("transport.in");
 for(int h=0;h<N;h++){
   for(int &v:P[h])f>>v;
   for(int &v:Q[h])f>>v;
   for(int &v:locked[h])f>>v;
   if(P[h][0]>=0)for(int y=0;y<N;y++){ip[h][P[h][y]]=y;iq[h][Q[h][y]]=y;}
 }
 int ct=0;
 for(int h=3;h<N;h++)for(int k=h+1;k<N;k++)if(h!=28&&k!=28){
   U seen=0;
   for(int t=0;t<N;t++)if(!(seen>>t&1)){
     vector<int> bfs{t};U comp=1ULL<<t;bool ok=true; int sigA=0,sigB=0;
     for(size_t i=0;i<bfs.size();i++){
       int y=bfs[i]; sigA|=1<<locked[h][y];sigB|=1<<locked[k][y];
       for(int v:{ip[k][P[h][y]],iq[k][Q[h][y]]}) if(!(comp>>v&1)){
        comp|=1ULL<<v;bfs.push_back(v);
       }
     }
     seen|=comp;
     if((sigA>>1)==(sigB>>1)){ct++;for(int y:bfs)sw[h][k][y]=sw[k][h][y]=comp;}
   }
 }
 cerr<<"safe components "<<ct<<"\n";
 unordered_map<int,int> ref;
 ifstream rf("rigid_51_contacts.json");string line;getline(rf,line);
 for(auto &c:line)if(!isdigit(c))c=' ';
 istringstream rr(line);int a,b;while(rr>>a>>b)ref[a]=b;ref[13]=16;ref[33]=19;
 int x,n;cin>>x>>n;U mask[3]={};vector<int> R;
 for(int i=0;i<n;i++){
   int y,h,c,d,e;cin>>y>>h>>c>>d>>e;opts.push_back({h,y%56,{c,d,e}});
   if(ref[h]==y%56)R.push_back(i);
   for(int k=0;k<3;k++){slots[k][opts[i].z[k]].push_back(i);mask[k]|=1ULL<<opts[i].z[k];}
 }
 U held=0;
 memset(owner3,-1,sizeof(owner3));
 for(int i:R)if(opts[i].h!=13&&opts[i].h!=33){
    const auto &z=opts[i].z;held|=1ULL<<z[0];refid[z[0]]=i;
    for(int k=0;k<3;k++){assert(owner3[k][z[k]]==-1);owner3[k][z[k]]=z[0];}
 }
 assert(__builtin_popcountll(held)==50);
 for(int i=0;i<n;i++)for(int k=1;k<3;k++)
    if(int v=owner3[k][opts[i].z[k]]; v>=0)bmask[i]|=1ULL<<v;
 for(budget=7;budget<=tries;budget++){
   copy(mask,mask+3,allow);sel.clear(); nodes=0;
   if(anchored(mask[0]^held,held,0)){
      cerr<<"found "<<budget<<" "<<nodes<<" "<<countLeaf<<" "<<countDbl<<" "<<moves.size()<<"\n";
      cout<<moves.size()<<"\n";
      for(auto [a,b,y]:moves)cout<<a<<" "<<b<<" "<<y<<"\n";
      for(int i:saved)cout<<opts[i].h<<" "<<opts[i].y<<"\n";
      return 0;
   }
   cerr<<budget<<" "<<nodes<<" "<<countLeaf<<" "<<countDbl<<"\n";
 }
 cerr<<"failed "<<countLeaf<<" "<<countDbl<<"\n";
 for(int k=0;k<=26;k++)if(counts[k])cerr<<k<<":"<<counts[k]<<" ";cerr<<"\n";
 return 1;
}
