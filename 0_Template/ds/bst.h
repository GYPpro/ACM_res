
const int N=1e5+100;
struct node{
	int s[2];
	int v,p,cnt,sz;
	void init(int p1,int v1){
		p=p1;v=v1;
		cnt=sz=1;
	}
}tr[N];
int root=0,idx=0;
void pushup(int x){
	tr[x].sz=tr[x].cnt+tr[tr[x].s[1]].sz+tr[tr[x].s[0]].sz;
}
void rotate(int x){
	int y=tr[x].p;
	int z=tr[y].p;
	int k=tr[y].s[1]==x;
	tr[y].s[k]=tr[x].s[k^1];
	tr[tr[x].s[k^1]].p=y;
	tr[z].s[tr[z].s[1]==y]=x;
	tr[x].p=z;
	tr[y].p=x;
	tr[x].s[k^1]=y;
	pushup(y);pushup(x);
}
void splay(int x,int k){s
	while(tr[x].p!=k){
		int y=tr[x].p;
		int z=tr[y].p;
		if(z!=k) (tr[y].s[0]==x)^(tr[z].s[0]==y) ? rotate(x) : rotate(y);
		rotate(x);
	}
	if(k==0) root=x;
}
void find(int v){
	int x=root;
	while(tr[x].v!=v && tr[x].s[v>tr[x].v] ) x=tr[x].s[v>tr[x].v];
	splay(x,0);
}
int get_pre(int v){
	find(v);
	int x=root;
	if(tr[x].v<v) return x;
	x=tr[x].s[0];
	while(tr[x].s[1]) x=tr[x].s[1];
	splay(x,0);
	return x;
}
int get_suc(int v){
	find(v);
	int x=root;
	if(tr[x].v>v) return x;
	x=tr[x].s[1];
	while(tr[x].s[0]) x=tr[x].s[0];
	splay(x,0);
	return x;
}
void del(int v){
	int pre=get_pre(v);
	int suc=get_suc(v);
	splay(pre,0);splay(suc,pre);
	int d=tr[suc].s[0];
	if(tr[d].cnt>1){
		tr[d].cnt--;splay(d,0);
	}
	else{
		tr[suc].s[0]=0;splay(suc,0);
	}
}
void insert(int v){
	int x=root;
	int p=0;
	while(x && tr[x].v!=v){
		p=x;x=tr[x].s[v>tr[x].v];
	}
	if(x) tr[x].cnt++;
	else{
		x=++idx;
		tr[p].s[v>tr[p].v]=x;
		tr[x].init(p,v);
	}
	splay(x,0);
}
int get_rank(int v){
	insert(v);
	int res=tr[tr[root].s[0]].sz;
	del(v);
	return res;
}
int get_val(int k){
	int x=root;
	while(1){
		int y=tr[x].s[0];
		if(tr[x].cnt+tr[y].sz<k){
			k-=tr[y].sz+tr[x].cnt;
			x=tr[x].s[1];
		}
		else{
			if(tr[y].sz>=k) x=tr[x].s[0];
			else break;
		}
	}
	splay(x,0);
	return tr[x].v;
}