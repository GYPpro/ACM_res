#include <template_overAll.h>

template <class T>
int subsearch_1(T l, T r, T tar, vector<T> &arr) // 第一类二分："升序"序列，"找到"就返回元素下标，相等元素返回最"左侧"下标，否则返回-1
{
    while (l <= r)
    {
        T m = (l + r) >> 1;
        if(arr[m] > tar) {
            
        }
    }
}