#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
using namespace std;
int main(){uint64_t target;string w;int tests=0;while(cin>>target>>w){if(w=="-")w="";int n=w.size()+2;vector<vector<bool>> a(n,vector<bool>(n));a[0][1]=1;for(int i=2;i<n;i++){int m=w[i-2]-'0';for(int j=0;j<2;j++){a[j][i]=(m>>j)&1;a[i][j]=!a[j][i];}for(int j=i+1;j<n;j++)a[i][j]=1;}vector<uint64_t> dp((1<<n)*n);for(int i=0;i<n;i++)dp[(1<<i)*n+i]=1;for(int s=1;s<(1<<n);s++)for(int i=0;i<n;i++)if(s>>i&1)for(int j=0;j<n;j++)if(!(s>>j&1)&&a[i][j])dp[(s|1<<j)*n+j]+=dp[s*n+i];uint64_t got=0;for(int i=0;i<n;i++)got+=dp[((1<<n)-1)*n+i];if(got!=target){cerr<<target<<" "<<got<<" "<<w<<endl;return 1;}tests++;}cout<<"Independent subset-DP checks passed: "<<tests<<endl;}
