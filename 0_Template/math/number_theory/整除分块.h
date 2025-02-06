#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define inf 100000000000000
const int N=2e6+7;
//题目，在1<=a<b<=n的条件下，求gcd(a,b)的和
//##IGNORE##
//这里用到了欧拉函数
//欧拉函数也可用欧拉筛求出
int main()
{
    void solve();
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    solve();
    return 0;
}
bool isprime[N];
vector<ll> p;
ll phi[N]={0,1};//边界条件
void eular(int n)
{
    for(int i=2;i<=n;i++)
    {
        if(!isprime[i])
        {
            p.push_back(i);
            phi[i]=i-1;
        }
        for(auto re:p)
        {
            if(i*re>n)break;
            isprime[i*re]=1;
            if(i%re==0){phi[i*re]=phi[i]*re;break;}
            phi[i*re]=phi[i]*(re-1);
        }
    }
    //到此，欧拉函数就求出来了
    for(int i=1;i<=n;i++)phi[i]+=phi[i-1];  //此处在求函数的前缀和，用于整除分块
}


ll cal(ll n,ll m)
{
    ll l=1,r=0,ans=0;
    while(l<=n)
    {
        r=min((n/(n/l)),(m/(m/l))); 
        ans+=(phi[r]-phi[l-1])*(n/l)*(m/l);
        l=r+1;
    }
    return ans;
}
void solve()
{
    ll n;
    cin>>n;
    eular(n);
    ll ans=(cal(n,n)-n*(n+1)/2)/2;
    cout<<ans<<"\n";
}