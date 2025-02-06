
#define	ll long long
#define N 5050
const ll inf=1e17+7;
int n;
struct node  //边
{
	int v;
	ll w;
	int is=0;
};
struct ed  //用于优先队列
{
	int v;
	ll w;
	ll d;
	bool operator < (const ed x) const 
	{
		return x.w<w;
	}
};
queue<ll> ans;  //答案队列
vector<node> e[N];  //邻接表
ll d[N]={0}; // d[i]终点到点i的最短路，也是点i到终点的最短路
void dj(int s)
{
	priority_queue<ed> q;
	vector<int> vis(n+10,0);
	for(int i=1;i<=n;i++)d[i]=inf;//初始化
	d[s]=0;
	q.push((ed){s,d[s],0});  //终点入队
	while(!q.empty())
	{
		auto r=q.top();
		q.pop();
		int u=r.v;

		if(vis[u])continue;
		vis[u]=1;
		d[u]=r.w;

		for(auto re:e[u])
		{
			int v=re.v;
			ll w=re.w;
			if(re.is)continue;//仅搜反边
			if(vis[v]||d[v]<=w+d[u])continue;
			d[v]=w+d[u];  //可松弛
			q.push((ed){v,d[v],0});
		}
	}
}
void astar(int s,int t,ll &k) //A*算法
{
	priority_queue<ed> q;
	vector<int> vis(n+10,0);
	q.push((ed){s,d[s],0});  //起点入队
	while(!q.empty())
	{
		auto temp=q.top();
		q.pop();

		int u=temp.v;

		vis[u]++;//统计出队次数

		if(vis[t]<=k&&u==t)//若只求第k短，只需vis[t]==k，记录并退出
		{
			ans.push(temp.d);
		}
		if(vis[t]==k)return;

		for(auto re:e[u])
		{
			if(!re.is)continue;
			if(vis[u]>k)continue;  //跳过出队次数大于k的

			int v=re.v;
			ll w=re.w;
			ll dis=temp.d;
			dis=dis+w;  //此为从起点到v的距离

			q.push((ed){v,d[v]+dis,dis});  
			//d[v]+dis是估计距离，dis是起点到v的距离
			//按估计距离升序排列
		}
	}
}
void solve()
{
	int s,t; //起点和终点
	int m;  //边数
	ll k;    //第k短
	cin>>n>>m>>k;

	s=n;t=1;   //处理终点和起点

	for(int i=1;i<=m;i++)
	{
		ll c;
		int u,v;
		cin>>u>>v>>c;
		e[u].push_back((node){v,c,1});//建立正边
		e[v].push_back((node){u,c,0});//建立反边
	}

	if(s==t)k++;  //起点与终点重合

	dj(t);//对t求dj，走反边

	astar(s,t,k);  //求s到t的最短的k条路

	for(int i=1;i<=k;i++)
	{
		if(s==t)
		{
			cout<<0<<"\n";
			continue;
		}
		if(!ans.empty())
		{
			cout<<ans.front()<<"\n";
			ans.pop();
		}
		else cout<<-1<<"\n";
	}
}