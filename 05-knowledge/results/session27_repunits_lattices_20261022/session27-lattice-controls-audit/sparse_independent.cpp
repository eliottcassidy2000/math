#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <climits>
using namespace std;
int main(int argc,char**argv){ifstream f(argc>1?argv[1]:"../session27-papers-review/downloaded_G96_reduced.txt");string s((istreambuf_iterator<char>(f)),{});for(char &c:s)if(c=='['||c==']'||c==';'||c==',')c=' ';istringstream in(s);long long g[96][96];for(auto &r:g)for(auto &x:r)if(!(in>>x))return 2;long long extra;if(in>>extra)return 3;
for(int r=1;r<=4;r++){array<int,4>a={0,1,2,3};long long mn=LLONG_MAX,count=0;while(true){for(int signs=0;signs<(1<<(r-1));signs++){long long v=0;for(int i=0;i<r;i++){int si=i==0?1:((signs>>(i-1)&1)?1:-1);for(int j=0;j<r;j++){int sj=j==0?1:((signs>>(j-1)&1)?1:-1);v+=si*sj*g[a[i]][a[j]];}}mn=min(mn,v);count++;}int k=r-1;while(k>=0&&a[k]==96-r+k)k--;if(k<0)break;a[k]++;for(int j=k+1;j<r;j++)a[j]=a[j-1]+1;}cout<<"support "<<r<<" count "<<count<<" minimum "<<mn<<endl;}}
