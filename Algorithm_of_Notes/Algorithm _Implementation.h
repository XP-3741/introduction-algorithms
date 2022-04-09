#pragma once
#include<iostream>
#include<vector>
#include<cstdlib>
#include<ctime>
using std::cout;
using std::cin;
using std::endl;
using std::vector;

#define ARRAY_SIZE 10
typedef int ElemType;
struct array_return_max
{
	int max_left_index = 0;
	int max_right_index = 0;
	ElemType sum = 0;
};
class heap_array
{
private:
	const int array_size;
	int heap_size;
	ElemType* Array;
public:
	heap_array(int length) 
		:array_size(length), heap_size(length) { Array = new ElemType[array_size + 1]; }
	~heap_array()
		{ delete[] Array; }
	ElemType operator[](int i) const
		{ return Array[i]; }
	ElemType& operator[](int i)
		{ return Array[i]; }
	int LENGTH() const { return array_size; }
	int& HEAP_SIZE() { return heap_size; }
};
struct TreeNode
{
	ElemType val;
	TreeNode* left;
	TreeNode* right;
	TreeNode* parent;
	TreeNode(int x) :val(x), left(nullptr), right(nullptr), parent(nullptr) {}
};
#define RED 0
#define BLACK 1
struct RedBlackTreeNode
{
	int color;
	ElemType key;
	RedBlackTreeNode* left;
	RedBlackTreeNode* right;
	RedBlackTreeNode* parent;
	RedBlackTreeNode(int COLOR, int VAL, RedBlackTreeNode* LEFT, RedBlackTreeNode* RIGHT, RedBlackTreeNode* PARENT)
		:color(COLOR), key(VAL), left(LEFT), right(RIGHT), parent(PARENT) {}
};
struct RedBlackTree
{
	RedBlackTreeNode* root;		// 树根
	RedBlackTreeNode* nil;		// 哨兵结点
};

// 归并排序(分治法)
void MERGE(ElemType* A, int p, int q, int r);
void MERGE_SORT(ElemType* A, int p, int r);
// 最大子数组(分治策略)
array_return_max FIND_MAX_CROSSING_SUBARRAY(ElemType* A, int low, int mid, int high);
array_return_max FIND_MAXIMUM_SUBARRAY(ElemType* A, int low, int high);
array_return_max FIND_MAXIMUM_SUBARRAY_LINEAR(ElemType* A, int array_size);	// 非递归、线性时间
// 堆排序
void HEAPSORT(heap_array& A);
// 快速排序
void QUICKSORT(ElemType* A, int p, int r);
// 计数排序
void COUNTING_SORT(ElemType* A, ElemType* B, ElemType k);
// 遗忘比较交换算法
void INSERTION_SORT(ElemType* A);
void INORDER_TREE_WALK(TreeNode* x);
// 二叉搜索树
TreeNode* TREE_SEARCH(TreeNode* x, ElemType k);
TreeNode* TREE_MINIMUM(TreeNode* x);
TreeNode* TREE_MAXIMUM(TreeNode* x);
TreeNode* TREE_SUCCESSOR(TreeNode* x);
void TREE_INSERT(TreeNode* TRoot, TreeNode* z);
void TRANSPLANT(TreeNode* TRoot, TreeNode* u, TreeNode* v);
void TREE_DELETE(TreeNode* TRoot, TreeNode* z);
// 红黑树