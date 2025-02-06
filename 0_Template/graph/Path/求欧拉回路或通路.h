#include<bits/stdc++.h>
using namespace std;
#define ll long long
const ll N=1055;
int main()
{
	ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
	void solve();
	int t=1;
	while(t--)solve();
	return 0;
}
//在已知有欧拉回路或通路的情况下求欧拉回路或欧拉通路
int m;
struct node{int u,v,i;};
vector<node> edge;
vector<pair<int,int>> e[N];
int vis[N]={0};
int ind[N]={0};
stack<int> stk;
bool cmp(pair<int,int> x,pair<int,int> y)
{
    int ix=x.second,iy=y.second;
    int u=x.first,v=y.first;
    return (edge[ix].v^edge[ix].u^u)<(edge[iy].u^edge[iy].v^v);
}
void dfs(int u)
{
    for(auto re:e[u])
    {
        int v=edge[re.second].v^edge[re.second].u^u;
        if(!vis[re.second])continue;
        vis[re.second]--;
        dfs(v);
    }
    stk.push(u);
}
void solve()
{
    memset(vis,0,sizeof(vis));
    int st=550;
    cin>>m;
    edge.clear();
    edge.push_back({-1,-1,0});
    for(int i=1;i<=500;i++)
    {
        ind[i]=0;
        e[i].clear();
    }
    for(int i=1;i<=m;i++)
    {
        int u,v;
        cin>>u>>v;
        edge.push_back({u,v,i});
        e[u].push_back({u,i});  //邻接表，但是存编号
        e[v].push_back({v,i});
        
        ind[u]++;  //记录顶点的度数
        ind[v]++;

        vis[i]++;  //记录边的可遍历次数
        
        st=min(st,min(u,v));
    }
    for(int i=1;i<=500;i++)
    {
        if(!e[i].size())continue;
        sort(e[i].begin(),e[i].end(),cmp);  //要求输出字典序最小的
    }
    for(int i=1;i<=500;i++)
    {
        if(ind[i]&1)  //若有奇数度的点，则需要从奇数度的点开始，找欧拉通路
        {
            st=i;
            break;
        }
    }
    dfs(st);
    while(!stk.empty())  //倒序输出
    {
        int u=stk.top();
        stk.pop();
        cout<<u<<"\n";
    }
}