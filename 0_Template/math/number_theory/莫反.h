#include<bits/stdc++.h>
using namespace std;
#define ll long long
const int N=5e4+7;
//题目：在1<=a<=n,1<=b<=m的条件下，求满足gcd(a,b)==d的(a,b)的个数

//莫比乌斯反演
//莫比乌斯函数可以用欧拉筛求出
bool isprime[N];
vector<ll> prime;
int mu[N];
void eular()
{
	mu[1]=1;//边界条件
	for(int i=2;i<N;i++)
	{
		if(!isprime[i])
		{
			prime.push_back(i);
			mu[i]=-1;
		}
		for(auto re:prime)
		{
			if(i*re>=N)break;
			isprime[i*re]=1;
			if(i%re==0)
			{
				mu[i*re]=0;
				break;
			}
			mu[i*re]=-mu[i];
		}
	}
	//到此，莫比乌斯函数已经求出

	for(int i=1;i<N;i++)mu[i]+=mu[i-1];  //这里在求函数的前缀和，用于整除分块
}
int main()
{
	ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
	int t=1;
	cin>>t;
	eular();
	void solve();
	while(t--)solve();
	return 0;
}
ll cal(int n,int m,int d)  //n,m,是两个范围的上下界，d为题目所求的gcd(a,b)==d
{
	ll ans=0;
	n/=d;
	m/=d;
	for(int l=1,r=0;l<=min(n,m);)
	{
		r=min(min(n/(n/l),m/(m/l)),min(n,m));
		ans+=(ll)(n/l)*(m/l)*(mu[r]-mu[l-1]);  
		//用莫比乌斯函数解题常常需要求出形如：  个数*个数*mu[i]   的公式
		l=r+1;
	}
	return ans;
}
void solve()
{
	int n,m,d;
	cin>>n>>m>>d;  
	cout<<cal(n,m,d)<<"\n";
}