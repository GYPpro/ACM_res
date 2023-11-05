// #include <basic\subsearch.h>
// #include <math\other\binpow.h>
// #include <string\trie_Tree.h>
#include <ds\segTree.h>
#include <iostream>

int main()
{
    int n;
    cin >> n;
    segTree<int> sgt(n);
    vector<int> arr(n);
    for(int i = 0;i < n;i ++) cin >> arr[i];
    sgt.build(arr);
    int l,r;
    cin >> l >> r;
    cout << sgt.getsum(l,r);
};