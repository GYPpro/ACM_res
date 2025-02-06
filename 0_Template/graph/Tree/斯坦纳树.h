/*
n(<=100)个点 m(<=500)条带权无向边G
给定k（<=10) 个节点的点集S，选出G的子图G1
使得S属于G1，G1是联通图  边权和最小
*/
#include<bits/stdc++.h>
using namespace std;
#define int long long
const int N=4e3;
const int inf=2e9;
int n,m,k,p[N],state;
int dp[N][N];
int head[N],to[N],ne[N],w[N],tot;
void add(int x,int y,int z){
	ne[++tot]=head[x];
	to[tot]=y,w[tot]=z;
	head[x]=tot;
}
queue<int>q;
bool vis[N];
void spfa(int s){
	while(!q.empty()){
		int u=q.front();q.pop();
		vis[u]=0;
		for(int i=head[u];i;i=ne[i]){
			int v=to[i];
			if(dp[v][s]>dp[u][s]+w[i]){
				dp[v][s]=dp[u][s]+w[i];
				if(!vis[v]) q.push(v),vis[v]=1;
			}
		}
	}
}
signed main(){
	cin>>n>>m>>k;
	for(int i=1;i<=m;i++){
		int x,y,z;cin>>x>>y>>z;
		add(x,y,z),add(y,x,z);
	}
	state=(1<<k)-1;
	for(int i=1;i<=n;i++) for(int s=0;s<=state;s++) dp[i][s]=inf;
	for(int i=1;i<=k;i++){
		cin>>p[i];
		dp[p[i]][1<<(i-1)]=0;
	}
	for(int s1=1;s1<=state;s1++){
		for(int i=1;i<=n;i++){
			for(int s2=s1&(s1-1);s2;s2=s1&(s2-1)) dp[i][s1]=min(dp[i][s1],dp[i][s2]+dp[i][s1^s2]);//枚举子集 
			if(dp[i][s1]<inf) q.push(i),vis[i]=1;//将这个点看成出发点 
		}
		spfa(s1);
	}
	cout<<dp[p[1]][state];//此时以哪个关键点为根都无所谓，答案是一样的 
}
