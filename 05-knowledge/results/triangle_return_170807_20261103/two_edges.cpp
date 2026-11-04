#define main inherited_search_main
#include "search.cpp"
#undef main
// All unordered pairs of distinct old-edge reversals; any resulting mark.
int main(){int n;I target;while(cin>>n>>target){if(n<2||n>18)return 2;T t(n);for(auto&v:t)cin>>v;vector<pair<int,int>>arcs;for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)arcs.push_back({a,b});
 int examined=0,hits=0;I nearest=INT64_MAX;cout<<"{\"n\":"<<n<<",\"target\":"<<target<<",\"hits\":[";
 for(int i=0;i<(int)arcs.size();i++)for(int j=0;j<i;j++){auto[a,b]=arcs[i];auto[c,d]=arcs[j];T u=t;u[a]^=1<<b;u[b]^=1<<a;u[c]^=1<<d;u[d]^=1<<c;auto s=stats(u);examined++;
  for(int r=0;r<n;r++){I z=s.h+2*s.E[r];nearest=min(nearest,z>target?z-target:target-z);if(z==target){if(hits++)cout<<',';cout<<"{\"arcs\":[["<<a<<','<<b<<"],["<<c<<','<<d<<"]],\"root\":"<<r<<",\"H\":"<<s.h<<",\"E2\":"<<s.E[r]<<",\"strong\":"<<(strong(u)?"true":"false")<<",\"rows\":";print(u);cout<<'}'<<flush;}}}
 cout<<"],\"examined\":"<<examined<<",\"nearest\":"<<nearest<<"}\n"<<flush;}}
