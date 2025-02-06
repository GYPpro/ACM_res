#include<bits/stdc++.h>
using namespace std;
const int N=1e5+100;
#define ll long long
#define int long long

struct node{
	int x,y,z;
}s[N];
vector<int> a(N),b(N);
int find(int x){
	return (a[x]==x) ? x : a[x]=find(a[x]);
}
void merge(int x,int y){
	a[find(x)]=find(y);
}
void solve(){
	int n,m;
	cin>>n>>m;
	vector<int> v(n+1);
	for(int i=0;i<=n;i++) a[i]=i;
	for(int i=0;i<m;i++) cin>>s[i].x>>s[i].y>>s[i].z;
	sort(s,s+m,[&](node a,node b){
		return a.z>b.z;
	});
	for(int i=0;i<m;i++){
		if(find(s[i].x)==find(s[i].y)){
			cout<<s[i].z;
			return;
		}
		if(!b[s[i].x]) b[s[i].x]=s[i].y;
		else merge(s[i].y,b[s[i].x]);
		if(!b[s[i].y]) b[s[i].y]=s[i].x;
		else merge(s[i].x,b[s[i].y]);
	}
	cout<<0;
	return;
}


signed main(){
	ios::sync_with_stdio(false);cin.tie(0);cout.tie(0); 
	int t=1;
//	cin>>t;
	while(t--) solve();
	return 0;
}
