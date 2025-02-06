void solve(){
	int n,k;
	cin>>n>>k;
	vector<int> a(n+1),ans(n+1);
	deque<int> q;
	for(int i=1;i<=n;i++){
		cin>>a[i];
		while(!q.empty()&&(q.front()<=i-k)) q.pop_front();
		while(!q.empty()&&a[q.back()]<=a[i]) q.pop_back();
		q.push_back(i);
		ans[i]=a[q.front()];
	}
}
