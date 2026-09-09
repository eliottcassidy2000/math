#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std; using I=long long; using Word=array<int,7>;
I gates=0,raw_words=0,evals=0;
void need(bool ok,string message){++gates;if(!ok){cerr<<message<<'\n';exit(1);}}
I mod(I a,I b){I r=a%b;return r<0?r+b:r;}
I floorq(I a,I b){return (a-mod(a,b))/b;}
struct Geo{int p,q;I L;vector<pair<I,I>> arcs;};
bool allowed(int s){for(int p=2;p*p<=s;p++){int e=0;while(s%p==0){s/=p;e++;}if(e&&(p%3!=2||e>2))return false;}return s==1||s%3==2;}
Geo geometry(int p,int q){
 Geo G{p,q,14LL*p*q,{}};vector<pair<I,I>> A,B,C;
 for(int k=0;k<=p;k++)A.push_back({max(0LL,(14LL*k-1)*q),min(G.L,(14LL*k+1)*q)});
 for(int k=0;k<=q;k++)B.push_back({max(0LL,(14LL*k-1)*p),min(G.L,(14LL*k+1)*p)});
 size_t i=0,j=0;while(i<A.size()&&j<B.size()){
  I a=max(A[i].first,B[j].first),b=min(A[i].second,B[j].second);if(a<b)C.push_back({a,b});
  if(A[i].second<B[j].second)i++;else if(A[i].second>B[j].second)j++;else{i++;j++;}
 }
 G.arcs.push_back({C.back().first-G.L,C.front().second});G.arcs.insert(G.arcs.end(),C.begin()+1,C.end()-1);return G;
}
pair<I,I> capacity(int n,const Geo&G){
 I base=0,cur=0;map<I,array<int,3>> events;
 for(auto[a,b]:G.arcs){I h=n*(b-a)/G.L,r=n*(b-a)%G.L,A=mod(n*a,G.L),B=mod(n*b,G.L);base+=h;
  if(!r)events[A][2]++;else{events[A][0]++;events[B][1]++;cur+=A>B;}}
 I walls=1000000000,chambers=1000000000;
 for(auto[e,v]:events){walls=min(walls,cur-v[1]-v[2]);cur+=v[0]-v[1];chambers=min(chambers,cur);}
 return {base+min(walls,chambers),base+chambers};
}
I direct(int n,const Geo&G,I phase2){I total=0;for(auto[a,b]:G.arcs)total+=-floorq(-2LL*n*b+phase2,2*G.L)-floorq(2LL*n*a-phase2,2*G.L)-1;return total;}
bool danger(I x,I den,int p){I z=mod(p*x,den);return 14*min(z,den-z)<den;}
I literal(int n,const Geo&G,I phase2){I total=0;for(int j=0;j<n;j++){I x=2*G.L*j+phase2,den=2*G.L*n;total+=danger(x,den,G.p)&&danger(x,den,G.q);}return total;}
void controls(){
 for(auto[p,q]:vector<pair<int,int>>{{1,1},{1,13},{1,112},{1,355},{11,263},{23,323},{33,320},{48,307}}){auto G=geometry(p,q);
  for(int n:{7,14,29,651,1260,1890,1939,1953,2485}){set<I>walls;for(auto[a,b]:G.arcs){walls.insert(mod(n*a,G.L));walls.insert(mod(n*b,G.L));}
   I closed=1000000000,open=1000000000,phase=0;for(I w:walls){I a=direct(n,G,2*w),b=direct(n,G,2*w+1);if(min(a,b)<closed){closed=min(a,b);phase=a<=b?2*w:2*w+1;}open=min(open,b);}
   need(capacity(n,G)==make_pair(closed,open),"all-wall and all-chamber literal integer interval control");need(literal(n,G,phase)==closed,"literal native-grid control");
  }
 }
 need(capacity(1890,geometry(48,307)).second==0,"old7560 absolute-depth hostile survives chamber refinement");
 need(capacity(1953,geometry(33,320)).second==0,"new7812 singleton slot is individually realizable");
 need(gcd(7812,1280)==4&&gcd(7812,132)==12,"new7812 singleton arm has exact actual margins");
}
int gd[91][91];array<unordered_set<unsigned long long>,7> profiles;vector<int> masks;
bool valid(const Word&w){
 int g=0;for(int d:w)g=gd[g][d];if(g!=1)return false;
 int cache[128]={};auto get=[&](auto&&self,int m)->int{if(cache[m])return cache[m];int bit=__builtin_ctz((unsigned)m),r=m&(m-1);return cache[m]=r?gd[w[bit]][self(self,r)]:w[bit];};
 for(int m:masks){int c=get(get,m),k=7-__builtin_popcount((unsigned)m),j=0;array<int,6>A{};
  for(int i=0;i<7;i++)if(!(m&(1<<i)))A[j++]=gd[c][w[i]];sort(A.begin(),A.begin()+k);
  unsigned long long key=c;for(int i=0;i<k;i++)key=(key<<7)|A[i];if(!profiles[k].count(key))return false;
 }return true;
}
void enumerate(const vector<int>&D,Word&w,int pos,int first,vector<Word>&out){
 if(pos==7){raw_words++;if(valid(w))out.push_back(w);return;}
 for(int j=first;j<(int)D.size();j++){w[pos]=D[j];enumerate(D,w,pos+1,j,out);}
}
using Graph=array<array<unsigned char,91>,91>;
bool connected(const Word&w,const Graph&G,int bit){int reached=1;
 while(true){int old=reached;for(int i=0;i<7;i++)if(reached&(1<<i))for(int j=0;j<7;j++)if(i!=j&&(G[w[i]][w[j]]&bit))reached|=1<<j;if(old==reached)return reached==127;}
}
void writeword(ostream&o,const Word&w){for(int i=0;i<7;i++)o<<' '<<w[i];}
int main(int argc,char**argv){
 need(argc==3,"input/output paths required");ifstream in(argv[1]);string dest=argv[2];
 for(int a=0;a<=90;a++)for(int b=0;b<=90;b++)gd[a][b]=gcd(a,b);
 int np;in>>np;for(int i=0;i<np;i++){int k,c;in>>k>>c;unsigned long long key=c;for(int j=0;j<k;j++){int x;in>>x;key=(key<<7)|x;}profiles[k].insert(key);}
 for(int k=6;k>=1;k--)for(int m=1;m<127;m++)if(__builtin_popcount((unsigned)m)==k)masks.push_back(m);need(masks.size()==126,"complete126 masks");
 int nd;in>>nd;vector<vector<int>> domains(nd);for(auto&D:domains){int d;in>>d;D.resize(d);for(int&x:D)in>>x;}
 int nt;in>>nt;vector<pair<int,int>> clocks(nt);set<int> quotients;
 for(auto&[t,d]:clocks){in>>t>>d;need(t%7==0,"zero budget clock");for(int e:domains[d]){need(t%(7*e)==0,"whole alphabet zero budget");quotients.insert(t/e);}}
 need(bool(in),"input complete");controls();vector<Geo> AT;for(int s=3;s<=356;s++)if(allowed(s))for(int p=1;p*2<s;p++)if(gcd(p,s-p)==1)AT.push_back(geometry(p,s-p));need(AT.size()==5855,"full native atlas");
 map<int,Graph> zero;map<int,array<array<int,91>,91>> slots;ofstream zf(dest+"/quotient_zero_pairs.txt");
 ofstream bf(dest+"/critical_banks.txt");set<int>critical={1939,1953};for(int a:{4,8,9,16,18,24})for(int b:{4,8,9,16,18,24})critical.insert(7056/gcd(a,b));
 for(int n:quotients){Graph G{};I local=0;
  array<array<int,91>,91> closed_counts{},open_counts{};
  for(auto&A:AT){int a=gcd(n,A.p),b=gcd(n,A.q);if(a>90||b>90||a%7==0||b%7==0)continue;if(a>b)swap(a,b);
   I lower=0;for(auto[x,y]:A.arcs)lower+=floorq((I)n*(y-x)-1,A.L);if(lower>0)continue;
   auto[c,o]=capacity(n,A);evals++;local++;need(c>=0&&o>=c,"wall/chamber capacity order");
   int bits=(c==0?1:0)|(o==0?2:0);if(bits){G[a][b]|=bits;G[b][a]|=bits;}
   if(bits&&critical.count(n))bf<<n<<' '<<A.p<<' '<<A.q<<' '<<gcd(n,A.p)<<' '<<gcd(n,A.q)<<' '<<bits<<'\n';
   if(c==0)closed_counts[a][b]+=a==b?2:1;
   if(o==0)open_counts[a][b]+=a==b?2:1;
  }
  for(int a=1;a<=90;a++)for(int b=a;b<=90;b++)if(G[a][b]){zf<<n<<' '<<a<<' '<<b<<' '<<(int)G[a][b]<<' '<<closed_counts[a][b]<<' '<<open_counts[a][b]<<'\n';open_counts[b][a]=open_counts[a][b];}zero.emplace(n,G);slots.emplace(n,open_counts);
 }
 cerr<<"QUOTIENTS "<<quotients.size()<<" EVALUATIONS "<<evals<<'\n';
 vector<vector<Word>> words(nd);ofstream df(dest+"/domain_counts.txt");
 for(int d=0;d<nd;d++){I before=raw_words;Word w{};enumerate(domains[d],w,0,0,words[d]);df<<d<<' '<<words[d].size()<<' '<<raw_words-before<<'\n';
  ofstream wf(dest+"/words_"+to_string(d)+".json");wf<<'[';bool first=true;for(auto&v:words[d]){if(!first)wf<<',';first=false;wf<<'[';for(int j=0;j<7;j++){if(j)wf<<',';wf<<v[j];}wf<<']';}wf<<']';
 }
 ofstream out(dest+"/survey_clocks.txt");int closed=0,openclosed=0,slotclosed=0;I totalwords=0;
 ofstream residual(dest+"/residual_7056.txt");
 for(auto[t,d]:clocks){Graph G{};for(int a:domains[d])for(int b:domains[d]){int e=gd[a][b];G[a][b]=zero.at(t/e)[a/e][b/e];}
  I cw=0,ow=0,sw=0;Word fc{},fo{},fs{};
  for(auto&w:words[d]){totalwords++;if(connected(w,G,1)){if(cw==0)fc=w;cw++;}if(connected(w,G,2)){if(ow==0)fo=w;ow++;
    array<int,91> mult{};for(int a:w)mult[a]++;bool possible=true;
    for(int a:domains[d])if(mult[a]&&!(G[a][a]&2)){int bound=0;for(int b:domains[d])if(a!=b&&mult[b]){int e=gd[a][b];bound+=mult[b]*slots.at(t/e)[a/e][b/e];}if(mult[a]>bound)possible=false;}
    if(possible){if(sw==0)fs=w;sw++;if(t==7056){for(int a:w)residual<<a<<' ';residual<<'\n';}}
  }}
  need(ow<=cw,"chamber zero graph contained in closed-wall zero graph");closed+=cw==0;openclosed+=ow==0;
  need(sw<=ow,"distinct-speed slot bound only removes chamber residuals");slotclosed+=sw==0;
  out<<t<<' '<<d<<' '<<words[d].size()<<' '<<cw<<' '<<ow<<' '<<sw;writeword(out,fc);writeword(out,fo);writeword(out,fs);out<<'\n';
 }
 cout<<"CLOCKS "<<nt<<" QUOTIENTS "<<quotients.size()<<" DOMAINS "<<nd<<"\nRAW_WORDS "<<raw_words<<" WORD_CLOCKS "<<totalwords<<"\nCLOSED_WALL_CLOSURES "<<closed<<" CHAMBER_CLOSURES "<<openclosed<<" SLOT_CLOSURES "<<slotclosed<<"\nCAPACITY_EVALUATIONS "<<evals<<" GATES "<<gates<<'\n';
}
