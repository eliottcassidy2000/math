#include <bits/stdc++.h>
using namespace std;
using ll=long long;
int m,E; vector<pair<int,int>> edges; vector<int> ctrl,other,stars;
ll countf(int lab){
 bool arc[12][12]={};
 for(int k=0;k<E;k++) {int v=ctrl[k],w=other[k]; int b=(lab>>k&1)?1:-1;
 for(int s=0;s<2;s++)for(int t=0;t<2;t++) { bool a=((s==0?1:-1)==-b);arc[2*v+s][2*w+t]=a;arc[2*w+t][2*v+s]=!a;}}
 ll d[64][12]={};for(int v=0;v<2*m;v++)d[1<<(v/2)][v]=1;
 for(int mask=1;mask<(1<<m);mask++)for(int v=0;v<2*m;v++)if(d[mask][v])for(int w=0;w<2*m;w++)if(!(mask>>(w/2)&1)&&arc[v][w]) d[mask|1<<(w/2)][w]+=d[mask][v];
 ll ans=0;for(int v=0;v<m;v++)ans+=d[(1<<m)-1][2*v];return ans;
}
int centralizer(int lab){vector<int>p(m);iota(p.begin(),p.end(),0);int c=0;do{bool ok=1;for(int k=0;k<E&&ok;k++){int v=p[ctrl[k]],w=p[other[k]];int j=0;while(!((edges[j].first==min(v,w))&&(edges[j].second==max(v,w))))j++;if(ctrl[j]!=v||((lab>>j&1)!=(lab>>k&1)))ok=0;}c+=ok;}while(next_permutation(p.begin(),p.end()));return c;}
int main(){for(m=1;m<=5;m++){E=m*(m-1)/2;edges.clear();for(int i=0;i<m;i++)for(int j=i+1;j<m;j++)edges.push_back({i,j});ll global=0;map<ll,int>hist;vector<tuple<int,int,int>>top;ll checks=0;
for(int q=0;q<(1<<E);q++){ctrl.resize(E);other.resize(E);stars.assign(m,0);for(int k=0;k<E;k++){auto [i,j]=edges[k];ctrl[k]=(q>>k&1)?i:j;other[k]=(q>>k&1)?j:i;stars[ctrl[k]]|=1<<k;}
 int reps=0;for(int s:stars)if(s)reps|=s&-s; vector<int>free;for(int k=0;k<E;k++)if(!(reps>>k&1))free.push_back(k); ll qm=0;
 for(int x=0;x<(1<<(int)free.size());x++){int a=reps;for(int k=0;k<(int)free.size();k++)if(x>>k&1)a|=1<<free[k];ll f=countf(a),A=f,score=0;int opt=a; bool zero=false;checks++;
 for(int v=0;v<m;v++)if(stars[v]){ll f1=countf(a^stars[v]);checks++;assert((f-f1)%2==0);ll B=(f-f1)/2;A-=B;score+=abs(B);if(B==0)zero=true;if(B<0)opt^=stars[v];}
 ll val=A+score;if(m==5&&val==195)assert(!zero);assert(countf(opt)==val);qm=max(qm,val);
 if(val>global){global=val;top.clear();}if(val==global)top.push_back({q,opt,centralizer(opt)});
 }hist[qm]++;}
 cout<<"m "<<m<<" max "<<global<<" DP_calls "<<checks<<" histogram";for(auto[a,b]:hist)cout<<" "<<a<<":"<<b;cout<<"\n";map<int,int>cs;for(auto[q,b,c]:top)cs[c]++;cout<<"top count "<<top.size()<<" centralizers";for(auto[c,n]:cs)cout<<" "<<c<<":"<<n;cout<<"\n";if(m==5){ofstream out("top_m5.txt");for(auto[q,b,c]:top)out<<q<<" "<<b<<" "<<c<<"\n";}
}}
