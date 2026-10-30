#include <iostream>
#include <vector>
#include <deque>
#include <limits>
#include <fstream>
#include <algorithm>
#include <cstdlib>
using namespace std;
int main(int argc,char**argv){
 if(argc!=3){cerr<<"usage excursion E output\n";return 2;}
 int E=atoi(argv[1]);if(E<1||E>4)return 2;
 int base=1;for(int i=0;i<E;i++)base*=5;
 int q=base*5,T=4*base,cut=min(64,T+1),D=8*q,N=5*q-5,S=q*q,root=q;
 long long inf=numeric_limits<long long>::max()/4;
 vector<long long> dis(2*S,inf),end(25,inf);vector<int> prev(2*S,-1),letter(2*S),ep(25,-1),el(25),inq(2*S);
 vector<int> inv(cut+1);inv[0]=1;for(int a=1;a<=cut;a++)inv[a]=inv[a-1]*((q+1)/2)%q;
 auto height=[&](int m,int b){int h=0;while(h<E && (m-1)%5==0 && b%5==0){m=(m-1)/5+1;b/=5;h++;}return h;};
 deque<int> que;dis[root]=0;que.push_back(root);inq[root]=1;long long scanned=0,relaxed=0;
 while(!que.empty()){
  int id=que.front();que.pop_front();inq[id]=0;int flag=id/S, native=id%S,m=native/q,b=native%q;
  for(int a=1;a<=cut;a++){
   int mm=(long long)3*m*inv[a]%q,bb=(long long)(3*b+1)*inv[a]%q;
   int ff=flag || a>1;
   long long val=dis[id]+(long long)N*a-(long long)D*height(mm,bb);scanned++;
   if(mm%base==1 && bb%base==0){int x=(mm-1)/base,y=bb/base,j=5*x+y;
    if(ff && val<end[j]){end[j]=val;ep[j]=id;el[j]=a;}continue;}
   int jd=ff*S+mm*q+bb;if(val<dis[jd]){dis[jd]=val;prev[jd]=id;letter[jd]=a;relaxed++;
    if(!inq[jd]){inq[jd]=1;if(!que.empty()&&val<dis[que.front()])que.push_front(jd);else que.push_back(jd);}}
  }
 }
 ofstream out(argv[2]);out<<"{\"cutoff\":"<<cut<<",\"E\":"<<E<<",\"scale\":"<<D<<",\"numerator\":"<<N<<",\"scanned\":"<<scanned<<",\"relaxed\":"<<relaxed<<",\"rows\":[";
 for(int j=0;j<25;j++){
  if(j)out<<",";vector<int>w;int cur=ep[j];if(cur<0){cerr<<"unreachable endpoint\n";return 3;}w.push_back(el[j]);
  while(cur!=root){w.push_back(letter[cur]);cur=prev[cur];if(w.size()>2*S){cerr<<"parent cycle\n";return 4;}}
  reverse(w.begin(),w.end());out<<"{\"x\":"<<j/5<<",\"y\":"<<j%5<<",\"cost\":"<<end[j]<<",\"word\":[";
  for(int k=0;k<(int)w.size();k++){if(k)out<<",";out<<w[k];}out<<"]}";
 }
 out<<"]}\n";cerr<<"E="<<E<<" scanned="<<scanned<<" relaxed="<<relaxed<<"\n";
 // Distances certify every first-return edge; their boundary values alone are not a global proof.
 ofstream cert(string(argv[2])+".dist");for(int f=0;f<2;f++)for(int m=1;m<q;m++)if(m%5)for(int b=0;b<q;b++){int id=f*S+m*q+b;if(dis[id]!=inf && (id==root || !(m%base==1&&b%base==0)))cert<<f<<" "<<m<<" "<<b<<" "<<dis[id]<<"\n";}
}
