#include"Algorithm _Implementation.h"

/*------------------------------------------归并排序(分治法)--------------------------------------------*/
// MERGE(A,p,q,r)
//	n1 = q - p + 1
//	n2 = r - q
//	Let L[1..n1 + 1] and R[1..n2 + 1] be new arrays
//	for i =  1 to n1
//		L[i] = A[p + i -1]
//	for j = 1 to n2
//		R[j] = A[q + j]
//	L[n1 + 1] = ><		// 哨兵牌
//	R[n2 + 1] = ><
//	i = 1
//	j = 1
//	for k = p to r
//		if L[i] <= R[j]
//			A[k] = L[i]
//			i = i + 1
//		else A[k] = R[j]
//			j = j + 1
void MERGE(ElemType*A, int p, int q, int r)
{
	const int n1 = q - p + 1;
	const int n2 = r - q;
	int* L = new int[n1+1];
	int* R = new int[n2+1];
	int i = 0;
	for (i = 0; i < n1; i++)
		L[i] = A[p + i];
	L[i] = INT_MAX;
	int j = 0;
	for (j = 0; j < n2; j++)
		R[j] = A[q + j + 1];
	R[j] = INT_MAX;
	for(int k = p,i=0,j=0;k<r;k++)
		if (L[i] <= R[j]) {
			A[k] = L[i];
			i++;
		}
		else {
			A[k] = R[j];
			j++;
		}
	delete[]L;
	delete[]R;
}
// MERGE-SORT(A,p,r)
//	if p < r
//		q = (p + r)/2
//		MERGE-SORT(A,p,q)
//		MERGE-SORT(A,q+1,r)
//		MERGE(A,p,q,r)
void MERGE_SORT(ElemType* A, int p, int r)
{
	if (p < r) {
		int q = (p + r) / 2;
		MERGE_SORT(A, p, q);
		MERGE_SORT(A, q + 1, r);
		MERGE(A, p, q, r);
	}
}

/*------------------------------------------最大子数组(分治策略)--------------------------------------------*/
// FIND-MAX-CROSSING-SUBARRAY(A,low,mid,high)		/* n */
//	left-sum = -><
//	sum = 0
//	for i = mid downto low
//		sum = sum + A[i]
//		if sum > left-sum
//			left-sum = sum
//			max-left = i
//	right-sum = -><
//	sum = 0
//	for j = mid + 1 to high
//		sum = sum + A[j]
//		if sum > right-sum
//			right-sum = sum
//			max-right = j
//	return (max-left,max-right,left-sum + right-sum)
array_return_max FIND_MAX_CROSSING_SUBARRAY(ElemType* A, int low, int mid, int high)
{
	array_return_max result;
	ElemType left_sum = -INT_MAX;
	ElemType sum = 0;
	for (int i = mid; i >= low; i--) {
		sum += A[i];
		if (left_sum < sum) {
			left_sum = sum;
			result.max_left_index = i;
		}
	}
	result.sum += left_sum;
	ElemType right_sum = -INT_MAX;
	sum = 0;
	for (int i = mid + 1; i <= high; i++) {
		sum += A[i];
		if (right_sum < sum) {
			right_sum = sum;
			result.max_right_index = i;
		}
	}
	result.sum += right_sum;
	return result;
}
// FIND-MAXIMUM-SUBARRAY(A,low,high)		/* nlgn */
//	if high == low
//		return (low,high,A[low])			// base case: only one element
//	else mid = (low + high) / 2 
//		(left-low, left-high, left-sum) = 
//			FIND-MAXIMUM-SUBARRAY(A,low,mid)
//		(right-low, right-high, right-sum) =
//			FIND-MAXIMUM-SUBARRAY(A,mid + 1,high)
//		(cross-low, cross-high, cross-sum) =
//			FIND-MAX-CROSSING-SUBARRAY(A, low, mid, high)
//		if left-sum >= right-sum and left-sum >= cross-sum
//			return (left-low, left-high, left-sum)
//		elseif right-sum >= left-sum and right-sum >= cross-sum
//			return (right-low, right-high, right-sum)
//		else return (cross-low, cross-high, cross-sum)
array_return_max FIND_MAXIMUM_SUBARRAY(ElemType* A, int low, int high)
{
	if (low == high) {
		array_return_max result;
		result.max_left_index = low;
		result.max_right_index = high;
		result.sum = A[low];
		return result;
	}
	else {
		int mid = (low + high) / 2;
		array_return_max left_part = FIND_MAXIMUM_SUBARRAY(A, low, mid);
		array_return_max right_part = FIND_MAXIMUM_SUBARRAY(A, mid + 1, high);
		array_return_max cross_part = FIND_MAX_CROSSING_SUBARRAY(A, low, mid, high);
		if (left_part.sum >= right_part.sum && left_part.sum >= cross_part.sum)
			return left_part;
		else if (right_part.sum >= left_part.sum && right_part.sum >= cross_part.sum)
			return right_part;
		else
			return cross_part;
	}
}
// 非递归、线性时间
// 若已知A[1..j]的最大子数组,基于如下性质将解扩展为A[1..j+1]的最大子数组:
// A[1..j+1]的最大子数组要么是A[1..j]的最大子数组
// 要么是某个子数组A[i..j+1] (1<=i<=j+1)
array_return_max FIND_MAXIMUM_SUBARRAY_LINEAR(ElemType* A, int array_size)
{
	array_return_max result;
	result.sum = A[0];
	for (int i = 1; i < array_size; i++) {
		for (int j = i, sum = 0; j >= 0; j--) {
			sum += A[j];
			if (result.sum < sum) {
				result.sum = sum;
				result.max_left_index = j;
				result.max_right_index = i;
			}
		}
	}
	return result;
}

/*------------------------------------------矩阵乘法的Strassen算法--------------------------------------------*/
// SQUARE-MATRIX-MULTIPLY(A,B)		/* n^3 */
//	n = A.rows
//  let C be new n*n matrix
//	for i = 1 to n
//		for j = 1 to n
//			c[ij] = 0
//			for k = 1 to n
//				c[ij] = c[ij] + a[ik]*b[kj]
//  return C

// SQUARE-MATRIX-MULTIPLY-RECURSIVE(A,B)		/* n^3 */
//	n = A.rows
//  let C be new n*n matrix
//	if n == 1
//		c[11] = a[11] * b[11]
//	else partition A,B and C as in equation (4.9)(P43)
//		C[11] = SQUARE-MATRIX-MULTIPLY-RECURSIVE(A[11],B[11])
//			+ SQUARE-MATRIX-MULTIPLY-RECURSIVE(A[12],B[21])
//		C[12] = SQUARE-MATRIX-MULTIPLY-RECURSIVE(A[11],B[12])
//			+ SQUARE-MATRIX-MULTIPLY-RECURSIVE(A[12],B[22])
//		C[21] = SQUARE-MATRIX-MULTIPLY-RECURSIVE(A[21],B[11])
//			+ SQUARE-MATRIX-MULTIPLY-RECURSIVE(A[22],B[21])
//		C[22] = SQUARE-MATRIX-MULTIPLY-RECURSIVE(A[21],B[12])
//			+ SQUARE-MATRIX-MULTIPLY-RECURSIVE(A[22],B[22])
//  return C

// Strassen方法		/* n^lg7 (n^2.81) */
//	<< Introduction to algorithms >> P45

/*------------------------------------------随机排列数组--------------------------------------------*/
// PERMUTE-BY-SORTING(A)
//	n = A.length
//	let P[1..n] be a new array	// 优先级数组
//  for i = 1 to n
//		P[i] = RANDOM(1, n^3)
//	sort A, using P as sort key 

// RANDOMIZE-IN-PLACE(A)		/* n */
//	n = A.length
//  for i = 1 to n
//		swap A[i] with A[RANDOM(i,n)]

/*------------------------------------------堆排序--------------------------------------------*/
/* 父结点 */
// PARENT(i)		
//	return [i/2]	// 向下取整
int PARENT(int i)
{
	return i / 2;
}

/* 左孩子 */
// LEFT(i)		
//	return 2i
int LEFT(int i)
{
	return 2 * i;
}

/* 右孩子 */
// RIGHT(i)		
//	return 2i+1
int RIGHT(int i)
{
	return 2 * i + 1;
}

/* 维护堆的性质 */
// MAX-HEAPIFY(A, i)
//  l = LEFT(i);
//  r = RIGHT(i);
//	if l <= A.heap_size and A[l] > A[i]
//		largest = l;
//	else
//		largest = i;
//  if r <= A.heap_size and A[r] > A[largest]
//		largest = r;
//	if largest != i
//		exchange A[i] with A[largest]
//		MAX-HEAPIFY(A, largest) 
void MAX_HEAPIFY(heap_array& A, int i)
{
	const int L = LEFT(i);
	const int R = RIGHT(i);
	int largest = 0;
	if (L<=A.HEAP_SIZE() && A[L]>A[i])
		largest = L;
	else
		largest = i;
	if (R<=A.HEAP_SIZE() && A[R]>A[largest])
		largest = R;
	if (i != largest) {
		int temp = A[i];
		A[i] = A[largest];
		A[largest] = temp;
		MAX_HEAPIFY(A, largest);
	}
}

/* 建堆(最大堆) */
// BUILD-MAX-HEAP(A)
//  A.heap_size = A.length
//	for i = [A.length/2] downto 1		// 从第一个非叶子结点开始
//		MAX-HEAPIFY(A, i)
void BUILD_MAX_HEAP(heap_array& A)
{
	for (int i = A.LENGTH() / 2; i >= 1; i--)
		MAX_HEAPIFY(A, i);
}

/* 堆排序算法 */
// HEAPSORT(A)
//	BUILD-MAX-HEAP(A)
//	for i = A.length downto 2
//		exchange A[1] with A[i]
//		A.heap_size = A.heap_size - 1
//		MAX-HEAPIFY(A, 1)
void HEAPSORT(heap_array& A)
{
	BUILD_MAX_HEAP(A);
	for (int i = A.LENGTH(); i >= 2; i--) {
		int temp = A[1];
		A[1] = A[i];
		A[i] = temp;
		A.HEAP_SIZE() = A.HEAP_SIZE() - 1;
		MAX_HEAPIFY(A, 1);
	}
}

/*------------------------------------------快速排序--------------------------------------------*/
// PARTITION(A, p, r)
//	x = A[r]
//	i = p-1
//	for j = p to r-1
//		if A[j] <= x
//			i = i + 1
//			exchange A[i] with A[j]
//	exchange A[i + 1] with A[r]
//	return i + 1
//
int PARTITION(ElemType* A, int p, int r)
{
	int x = A[r];
	int i = p - 1;
	for (int j = p; j <= r - 1; j++)
		if (A[j] < x) {
			i++;
			int temp = A[j];
			A[j] = A[i];
			A[i] = temp;
		}
	int temp = A[r];
	A[r] = A[i + 1];
	A[i + 1] = temp;
	return i + 1;
}
// QUICKSORT(A, p, r)
//	if p < r
//		q = PARTITION(A, p, r)
//		QUICKSORT(A, p, q-1)
//		QUICKSORT(A, q+1, r)
void QUICKSORT(ElemType* A, int p, int r)
{
	if (p < r) {
		int q = PARTITION(A, p, r);
		QUICKSORT(A, p, q - 1);
		QUICKSORT(A, q + 1, r);
	}
}

/*------------------------------------------计数排序--------------------------------------------*/
// COUNTING-SORT(A, B, k)				/* n */
//	let C[0..k] be new array
//	for i = 0 to k
//		C[i] = 0
//	for j = 1 to A.length
//		C[A[j]] = C[A[j]] + 1
//	//C[i] now contians the number of elements equal to i.
//	for i = 1 to k
//		C[i] = C[i] + C[i-1]
//	//C[i] now contians the number of elements less than or equal to i.
//	for j = A.length downto 1
//		B[C[A[j]]] = A[j]
//		C[A[j]] = C[A[j]] - 1
void COUNTING_SORT(ElemType* A, ElemType* B, ElemType k)
{
	ElemType* C = new ElemType[k + 1];
	for (int i = 0; i <= k; i++)
		C[i] = 0;
	for (int j = 0; j < ARRAY_SIZE; j++)
		C[A[j]] = C[A[j]] + 1;
	for (int i = 1; i <= k; i++)
		C[i] = C[i] + C[i - 1];
	for (int j = ARRAY_SIZE - 1; j >= 0; j--) {
		B[C[A[j]]] = A[j];
		C[A[j]] = C[A[j]] - 1;
	}
	delete[]C;
}

/*------------------------------------------桶排序--------------------------------------------*/
// BUCKET-SORT(A)				/* n */
//  n = A.length
//	let B[0...n-1] be a new array
//	for i = 0 to n-1
//		make B[i] an empty list
//	for i = 1 to n
//		insert A[i] into list B[(int)nA[i]]
//	for i = 0 to n-1
//		sort list B[i] with insertion sort
//	concatenate the lists B[0],B[1],...,B[n-1] together in order

/*------------------------------------------遗忘比较交换算法--------------------------------------------*/
// COMPARE-EXCHANGE(A, i, j)
//	if A[i] > A[j]
//		exchange A[i] with A[j]

// INSERTION-SORT(A)
//	for j = 2 to A.length
//		for i = j-1 downto 1
//			COMPARE-EXCHANGE(A, i, i+1)

void COMPARE_EXCHANGE(ElemType* A, int i, int j)
{
	if (A[i] > A[j]) {
		int temp = A[i];
		A[i] = A[j];
		A[j] = temp;
	}
}
void INSERTION_SORT(ElemType* A)
{
	for (int j = 1; j < ARRAY_SIZE; j++)
		for (int i = j - 1; i >= 0; i--)
			COMPARE_EXCHANGE(A, i, i + 1);
}

/*------------------------------------------期望为线性时间的选择算法--------------------------------------------*/
// RANDOMIZED-SELECT(A, p, r, i)
//	if p == r
//		return A[p]
//	q = RANDOMIZED-PARTITION(A, p, r)
//	k = q - p + 1
//	if i == k
//		return A[q]
//	else if i < k
//		return RANDOMIZED-SELECT(A, p, q-1, i)
//	else
//		return RANDOMIZED-SELECT(A, q+1, r, i-k)
//

/*------------------------------------------二叉搜索树--------------------------------------------*/
// 中序遍历
// INORDER-TREE-WALK(x)
//	if x != NIL
//		INORDER-TREE-WALK(x.left)
//		print x.key
//		INORDER-TREE-WALK(x.right)
void INORDER_TREE_WALK(TreeNode* x)
{
	if (x != nullptr) {
		INORDER_TREE_WALK(x->left);
		cout << x->val << " ";
		INORDER_TREE_WALK(x->right);
	}
}

// 查找关键字 k
// TREE-SEARCH(x,k)
//	if x == NIL or k == x.key
//		return x
//	if k < x.key
//		return TREE-SEARCH(x.left,k)
//	else
//		return TREE-SEARCH(x.right,k)
TreeNode* TREE_SEARCH(TreeNode* x, ElemType k)
{
	if (x == nullptr || k == x->val)
		return x;
	if (k < x->val)
		return TREE_SEARCH(x->left, k);
	else
		return TREE_SEARCH(x->right, k);
}

// 最小关键字
// TREE-MINIMUM(x)
//	while x.left != NIL
//		x = x.left
//	return x
TreeNode* TREE_MINIMUM(TreeNode* x)
{
	while (x->left != nullptr)
		x = x->left;
	return x;
}

// 最大关键字
// TREE-MAXIMUM(x)
//	while x.right != NIL
//		x = x.right
//	return x
TreeNode* TREE_MAXIMUM(TreeNode* x)
{
	while (x->right != nullptr)
		x = x->right;
	return x;
}

// 后继
// TREE_SUCCESSOR(x)
//	if x.right != NIL
//		return TREE-MINIMUM(x.right)
//	y = x.p
//	while y != NIL and x == y.right
//		x = y
//		y = y.p
//	return y
TreeNode* TREE_SUCCESSOR(TreeNode* x)
{
	if (x->right != nullptr)
		return TREE_MINIMUM(x->right);
	TreeNode* y = x->parent;
	while (y != nullptr && x == y->right) {
		x = y;
		y = y->parent;
	}
	return y;
}

// 插入
// TREE-INSERT(T,z)
//	y = NIL
//	x = T.root
//	while x != NIL
//		y = x
//		if z.key < x.key
//			x = x.left
//		else
//			x = x.right
//	z.p = y
//  if y == NIL
//		T.root = z		// tree T was empty
//	elseif z.key < y.key
//		y.left = z
//	else
//		y.right = z
void TREE_INSERT(TreeNode* TRoot, TreeNode* z)
{
	TreeNode* y = nullptr;
	TreeNode* x = TRoot;
	while (x != nullptr) {
		y = x;
		if (z->val < x->val)
			x = x->left;
		else
			x = x->right;
	}
	z->parent = y;
	if (y == nullptr)
		TRoot = z;
	else if (z->val < y->val)
		y->left = z;
	else
		y->right = z;
}

// 用另一棵子树替换一棵子树并成为其双亲的孩子
// TRANSPLANT(T,u,v)
//	if u.p == NIL
//		T.root = v
//	elseif u == u.p.left
//		u.p.left = v
//	else
//		u.p.right = v
//	if v != NIL
//		v.p = u.p
void TRANSPLANT(TreeNode* TRoot, TreeNode* u, TreeNode* v)
{
	if (u->parent == nullptr)
		TRoot = v;
	else if (u == u->parent->left)
		u->parent->left = v;
	else
		u->parent->right = v;
	if (v != nullptr)
		v->parent = u->parent;
}

// 删除
// TREE-DELETE(T,z)
//	if z.left == NIL
//		TRANSPLANT(T,z,z.right)
//	elseif z.right == NIL
//		TRANSPLANT(T,z,z.left)
//	else y = TREE-MINiMUM(z.right)
//		if y.p != z
//			TRANSPLANT(T,y,y.right)
//			y.right = z.right
//			y.right.p = y
//		TRANSPLANT(T,z,y)
//		y.left = z.left
//		y.left.p = y
void TREE_DELETE(TreeNode* TRoot, TreeNode* z)
{
	if (z->left == nullptr)
		TRANSPLANT(TRoot, z, z->right);
	else if(z->right==nullptr)
		TRANSPLANT(TRoot, z, z->left);
	else {
		TreeNode* y = TREE_MINIMUM(z->right);
		if (y->parent != z) {
			TRANSPLANT(TRoot, y, y->right);
			y->right = z->right;
			y->right->parent = y;
		}
		TRANSPLANT(TRoot, z, y);
		y->left = z->left;
		y->left->parent = y;
	}
}

/*------------------------------------------红黑树--------------------------------------------*/
// 旋转
// LEFT-ROTATE(T, x)			/* 1 */
//	y = x.right				// set y
//	x.right = y.left		// turn y's left subtree into x's right subtree
//	if y.left != T.nil
//		y.left.p = x
//	y.p = x.p				// link x's parent to y
//	if x.p == T.nil
//		T.root = y
//	elseif x == x.p.left
//		x.p.left = y
//	else x.p.right = y
//	y.left = x				// put x on y's left
//	x.p = y
void LEFT_ROTATE(RedBlackTree* T, RedBlackTreeNode* x)
{
	RedBlackTreeNode* y = x->right;
	x->right = y->left;
	if (y->left != T->nil)
		y->left->parent = x;
	y->parent = x->parent;
	if (x->parent == T->nil)
		T->root = y;
	else if (x == x->parent->left)
		x->parent->left = y;
	else
		x->parent->right = y;
	y->left = x;
	x->parent = y;
}

// RIGHT-ROTATE(T, y)
//	x = y.left				// set x
//	y.left = x.right		// turn x's right subtree into y's left subtree
//	if x.right != T.nil
//		x.right.p = y
//	x.p = y.p				// link y's parent to x
//	if y.p == T.nil
//		t.root = x
//	elseif y == y.p.left
//		y.p.left = x
//	else y.p.right = x
//	x.right = y				// put y on x's right
//	y.p = x
void RIGHT_ROTATE(RedBlackTree* T, RedBlackTreeNode* y)
{
	RedBlackTreeNode* x = y->left;
	y->left = x->right;
	if (x->right != T->nil)
		x->right->parent = y;
	x->parent = y->parent;
	if (y->parent == T->nil)
		T->root = x;
	else if (y == y->parent->left)
		y->parent->left = x;
	else
		y->parent->right = x;
	x->right = y;
	y->parent = x;
}

// 插入
// RB-INSERT-FIXUP(T,z)
//	while z.p.color == RED
//		if z.p == z.p.p.left
//			y = z.p.p.right
//			if y.color == RED
//				z.p.color = BLACK			// case 1
//				y.color = BLACK				// case 1
//				z.p.p.color = RED			// case 1
//				z = z.p.p					// case 1
//			elseif z == z.p.right
//				z = z.p						// case 2
//				LEFT-ROTATE(T, z)			// case 2
//			z.p.color = BLACK				// case 3
//			z.p.p.color = RED				// case 3
//			RIGHT-ROTATE(T, z.p.p)			// case 3
//		else(same as then clause
//				with "right" and "left" exchanged)
//	T.root.color = BLACK
void RB_INSERT_FIXUP(RedBlackTree* T, RedBlackTreeNode* z)
{
	while (z->parent->color == RED) {
		if (z->parent == z->parent->parent->left) {
			RedBlackTreeNode* y = z->parent->parent->right;
			if (y->color == RED) {
				z->parent->color = BLACK;
				y->color = BLACK;
				z->parent->parent->color = RED;
				z = z->parent->parent;
			}
			else {
				if (z == z->parent->right) {
					z = z->parent;
					LEFT_ROTATE(T, z);
				}
				z->parent->color = BLACK;
				z->parent->parent->color = RED;
				RIGHT_ROTATE(T, z->parent->parent);
			}
		}
		else {
			RedBlackTreeNode* y = z->parent->parent->left;
			if (y->color == RED) {
				z->parent->color = BLACK;
				y->color = BLACK;
				z->parent->parent->color = RED;
				z = z->parent->parent;
			}
			else {
				if (z == z->parent->left) {
					z = z->parent;
					RIGHT_ROTATE(T, z);
				}
				z->parent->color = BLACK;
				z->parent->parent->color = RED;
				LEFT_ROTATE(T, z->parent->parent);
			}
		}
	}
	T->root->color = BLACK;
}

// RB-INSERT(T, z)				/* lgn */
//	y = T.nil
//	x = T.root
//  while x != T.nil
//		y = x
//		if z.key < x.key
//		x = x.left
//		else x = x.right
//	z.p = y
//	if y == T.nil
//		T.root = z
//	elseif z.key < y.key
//		y.left = z
//	else y.right = z
//	z.left = T.nil
//	z.right = T.nil
//	z.color = RED
//	RB-INSERT-FIXUP(T,z)
void RB_INSERT(RedBlackTree* T, RedBlackTreeNode* z)
{
	RedBlackTreeNode* y = T->nil;
	RedBlackTreeNode* x = T->root;
	while (x != T->nil) {
		y = x;
		if (z->key < x->key)
			x = x->left;
		else
			x = x->right;
	}
	z->parent = y;
	if (y == T->nil)
		T->root = z;
	else if (z->key < y->key)
		y->left = z;
	else
		y->right = z;
	z->left = T->nil;
	z->right = T->nil;
	z->color = RED;
	RB_INSERT_FIXUP(T, z);
}

// 删除
// RB-TRANSPLANT(T, u, v)
//	if u.p == T.nil
//		T.root = v
//	elseif u == u.p.left
//		u.p.left = v
//	else u.p.right = v
//	v.p = u.p

// RB-DELETE(T, z)
//	y = z
//	y-original-color = y.color
//	if z.left == T.nil
//		x = z.right
//		RB-TRANSPLANT(T, z, z.right)
//	elseif z.right == T.nil
//		x = z.left
//		RB-TRANSPLANT(T, z, z.left)
//	else y = TREE-MINIMUN(z.right)
//		y-original-color = y.color
//		x = y.right
//		if y.p == z
//			x.p = y
//		else RB-TRANSPLANT(T, y, y.right)
//			y.right = z.right
//			y.right.p = y
//		RB-TRANSPLANT(T, z, y)
//		y.left = z.left
//		y.left.p = y
//		y.color = z.color
//	if y-original-color == BLACK
//		RB-DELETE-FIXUP(T, x)

// RB-DELETE-FIXUP(T, x)
//	while x != T.root and x.color == BLACK
//		if x == x.p.left
//			w = x.p.right
//			if w.color == RED
//				w.color = BLACK										// case 1
//				x.p.color = RED										// case 1
//				LEFT-ROTATE(T, x.p)									// case 1
//				w = x.p.right										// case 1
//			if w.left.color == BLACK and w.right.color == BLACK
//				w.color = RED										// case 2
//				x = x.p												// case 2
//			elseif w.right.color == BLACK
//					w.left.color = BLACK							// case 3
//					w.color = RED									// case 3
//					RIGHT-ROTATE(T, w)								// case 3
//					w = x.p.right									// case 3
//				w.color = x.p.color									// case 4
//				x.p.color = BLACK									// case 4
//				w.right.color = BLACK								// case 4
//				LEFT-ROTATE(t, x.p)									// case 4
//				x = T.root											// case 4
//		else (same as then clause with "right" and "left" exchanged)
//	x.color = BLACK
//

/*------------------------------------------动态规划--------------------------------------------*/

/*	设计一个动态规划算法的四个步骤:
*		1. 刻划一个最优解的结构特征
*		2. 递归地定义最优解的值
*		3. 计算最优解的值,通常采用自底向上的方法
*		4. 利用计算出的信息构造一个最优解
*	力扣评论四步骤:
*		1. 定义子问题
*		2. 写出子问题的递推关系
*		3. 确定 DP 数组的计算顺序(DP数组:子问题数组)
*		4. 空间优化
*/


/* ················ 15.1_钢条切割 ················
   如何将长钢条切割成短钢条,使得总价值最高 */
/*
* 自然递归
	法一: Rn = max(Pn, R1+Rn-1, R2+Rn-2, ..., Rn-1+R1)
			Rn - 长度n英尺长钢切割的最大收益
			Pn - 长度为n的钢条的收益
	法二: Rn = max( Pi + Rn-i)	(1<=i<=n)
			左边不切割(Pi),右边递归切割(Rn-i)
*/
// 自顶向下递归实现(法二)
// CUT-ROD(p, n)
//	if n == 0
//		return 0
//	q = -∞
//	for i = 1 to n
//		q = max(q,p[i]+CUT-ROD(p,n-i))
//	return q
/* CUT-ROD 的效率很差,原因在于,CUT-ROD 反复地用相同的参数值对自身进行递归调用
   及反复求解相同的子问题 */

/*
*  动态规划
	法一: 带备忘的自顶向下法
			此方法仍按自然的递归形式编写程序,但过程会保存每个子问题的解
			当需要一个子问题的解时,过程首先检查是否已经保存过此解
	法二: 自底向上法
			这种方法一般需要恰当定义子问题"规模"的概念
			使得任何子问题都只依赖于"更小的"子问题的求解
			因为我们可以将子问题按规模排序,按由小到大的顺序进行求解
			当求解某个子问题时,它所依赖的那些更小的子问题都已经求解完毕,结果已保存
			每个子问题只需求解一次,当我们求解它(第一次遇到它)时
			它的所有前提子问题都已求解完毕
*/
// 自顶向下,加入备忘机制(法一)
// MEMOIZED-CUT-ROD(p, n)
//	let r[0..n] be a new array
//	for i = 0 to n
//		r[i] = -∞
//	return MEMOIZED-CUT-ROD-AUX(p,n,r)

// MEMOIZED-CUT-ROD-AUX(p, n, r)
//	if r[n] >= 0
//		return r[n]
//	if n == 0
//		q = 0
//	else q = -∞
//		for i = 1 to n
//			q = max(q,p[i]+MEMOIZED-CUT-ROD-AUX(p,n-i,r))
//	r[n] = q
//	return q

// 自底向上(法二)
// BOTTOM-UP-CUT-ROD(p, n)
//	let r[0..n] be a new array 
//	r[0] = 0
//	for j = 1 to n
//		q = -∞
//		for i = 1 to j
//			q = max(q,p[i]+r[j-i])
//		r[j] = q
//	return r[n]

/* 重构解
   不仅计算最大收益Ri,还保存最优解对应的第一段钢条的切割长度 */
// EXTENDED-BOTTOM-UP-CUT-ROD(p, n)
//	let r[0..n] and s[0..n] be new arrays
//	r[0] = 0
//	for j = 1 to n
//		q = -∞
//		if q < p[i]+r[j-i]
//			q = p[i]+r[j-i]
//			s[j] = i 
//		r[j] = q
//	return r and s

// PRINT-CUT-ROD-SOLUTION(p, n)
//	(r,s) = EXTENDED-BOTTOM-UP-CUT-ROD(p, n)
//	while n > 0
//		print s[n]
//		n = n - s[n]
/* 接受两个参数,价格表 p 和钢条长度 n
   然后调用 EXTENDED-BOTTOM-UP-CUT-ROD 来计算切割下来的每段钢条长度 s[1..n]
   最后输出长度为 n 的钢条的完整的最优切割方案 */

/* ················ 15.2 矩阵链乘法 ···············
   如何用最少的标量乘法操作完成一个矩阵链相乘的运算
*/
// 两矩阵相乘的标准算法
// MATRIX-MULTIPLY(A,B)
//	if A.columns != B.rows
//		error"incompatible dimensions"
//	else let C be a new A.row*B.columns matrix
//		for i = 1 to A.rows
//			for j = 1 to B.columns
//				Cij = 0
//				for k = 1 to A.columns
//					Cij = Cij + Aik * Bkj
//	return C 

/* 自底向上表格法
*	递推公式:
*		m[i, j] = 0										如果 i = j
*				  min{m[i, k] + m[k+1, j] + Pi-1PkPj}	如果 i < j

	假定矩阵 Ai 的规模 Pi-1 × Pi (i = 1, 2, ..., n)
	它的输入是一个序列 P = <P0, P1, ..., Pn>,其长度为 P.length = n + 1
	此过程用一个辅助表 m[1..n, 1..n] 来保存代价 m[i, j]
	用另一个辅助表 s[1..n-1, 2..n] 记录最优值 m[i, j] 对应的分割点 k
*/
// MATRIX-CHAIN-OEDER(p)
//	n = p.length - 1
//	let m[1..n, 1..n] and s[1..n-1, 2..n] be new tables
//	for i = 1 to n
//		m[i,i] = 0
//	for L = 2 to n				// L is the chain length
//		for i = 1 to n - L + 1
//			j = i + L - 1
//			m[i,j] = ∞
//			for k = i to j -1
//				q = m[i,k] + m[k+1,j] + Pi-1PkPj
//				if q < m[i,j]
//					m[i,j] = q
//					s[i,j] = k
//	return m and s

/*	构造最优解
	输出 <Ai, Ai+1, ..., Aj> 的最有括号方案
	起输入为 MATRIX-CHAIN-OEDER 得到的表 s 及下标 i 和 j
*/
// PRINT-OPTIMAL-PARENS(s, i, j)
//	if i == j
//		print "A"i
//	else print "("
//			PRINT-OPTIMAL-PARENS(s, i, s[i,j])
//			PRINT-OPTIMAL-PARENS(s, s[i,j]+1, j)
//		 print ")"


/* ················ 15.3 动态规划原理 ···············
   适合动态规划的两个关键特征
	- 最优子结构
	- 子问题重叠
*/
/*	最优子结构(P216)
		用动态规划方法求解问题的第一步就是刻画最优解的结构
		如果一个问题的最优解包含其子问题的最优解,就称此问题具有 最优子结构 性质
		使用动态规划时,我们用子问题的最优解来构造原问题的最优解
		因此,要必须小心确保考察了最优解中用到的所有子问题

		在发掘最优子结构性质的过程中,实际上遵循了如下的通用模式:
		1.证明问题最优解的第一个组成部分是做出一个选择,做出这次选择会产生一个或多个待解的子问题
		2.对于一个给定问题,在其可能的第一部选择中,你假定已经知道哪种选择才会得到最优解
		  你现在并不关心这种选择具体如何得到,只是假定已经知道了这种选择
		3.给定可获得最优解选择后,你确定这次选择会产生哪些子问题
	      以及如何最好地刻画子问题空间
		4.利用"剪切-粘贴"技术证明:作为构成原问题的最优解的组成部分,每个子问题的解就是它本身的最优解

		一个刻画子问题空间的好经验是:保持子问题空间尽可能简单,只在必要时才扩展它

		对于不同问题领域,最有子结构的不同体现在两方面:
			- 原问题的最优解中涉及多少个子问题,以及
			- 在确定最优解使用哪些子问题时,我们需要考察多少种选择

		在动态规划中,我们通常自底向上地使用最优子结构
		也就是说,首先求得子问题的最优解,然后求原问题的最优解
		在求解原问题过程中,我们需要在涉及的子问题中做出抉择,选出能够得到原问题最优解的子问题
		原问题最优解的代价通常就是子问题最优解的代价再加上由此次选择直接产生的代价

		子问题之间必须是无关的
		这里,子问题无关的含义是,同一个原问题的一个字问题的解不影响另一个子问题的解
		换个角度来看,我们所面临的困境就是:
		求解一个子问题时用到了某些资源(顶点等),导致这些资源在求解其他子问题时不可用
*/
/*	重叠子问题
		适合用动态规划方法求解的最优化问题应具有的第二个性质是
		子问题空间必须足够"小",及问题的递归算法会反复求解相同的子问题,而不是一直生成新的子问题
		一般来讲,不同子问题总数是输入规模的多项式函数为好
		如果递归算法反复求解相同的子问题,我们就称最优化问题具有 重叠子问题 性质
		与之相对的,适合用分治方法求解的问题通常在递归的每一步都生成全新的子问题
*/
/*	重构最优解
		从实际考虑,我们通常将每个子问题所做的选择存放在一个表中
		这样就不必更具代价来重构这些信息
*/
/*	备忘
		我们可以通过自顶向下策略,同时达到与自底向上动态规划方法相似的效率
		思路就是对自然但低效率的递归算法加入 备忘 机制
		与自底向上方法一样,我们维护一个表记录子问题的解,但仍保持递归算法的控制流程

		带备忘录的递归算法为每个子问题维护一个表项来保存它的解
		每个表项的初值设为一个特殊值,表示尚未填入子问题的解
		当递归调用过程中第一次遇到子问题时,计算其解,并存入对应表项
		随后每次遇到同一个子问题,只是简单的查表,返回其解
*/


/* ················ 15.4  最长公共子序列 ···············
   如何用动态规划方法找到两个序列的最长公共子序列
   给定两个序列 X=<X1, X2, ..., Xm> 和 Y=<Y1, Y2, ..., Yn>
   求 X 和 Y 长度最长的公共子序列
*/
/*
	定义 c[i,j] 表示 Xi 和 Yj 的 LCS 的长度
	如果 i = 0 或 j = 0 , 即一个序列长度为 0,那么LCS的长度为 0
	根据LCS问题的最优子结构性质,可得:
					0						若 i = 0 或 j = 0
		c[i,j] =	c[i-1, j-1] + 1			若 i,j > 0 且 Xi = Yj
					max(c[i,j-1],c[i-1,j])	若 i,j > 0 且 Xi != Yj
*/

// LCS-LENGTH(X,Y)
//	m = X.length
//	n = Y.length
//	let b[1..m, 1..n] and c[0..m, 0..n] be new tables
//	for i = 1 to m
//		c[i,0] = 0
//	for j = 0 to n
//		c[0,j] = 0
//	for i = 1 to m
//		for j = 1 to n
//			if Xi == Yj
//				c[i,j] = c[i-1,j-1] + 1
//				b[i,j] = "↖"
//			elseif c[i-1,j]>=c[i,j-1]
//				c[i,j] = c[i-1,j]
//				b[i,j] = "↑"
//			else c[i,j] = c[i,j-1]
//				b[i,j] = "←"
/*
	LCS-LENGTH 接受两个序列 X = <X1, X2, ..., Xm> 和 Y = <Y1, Y2, ..., Yn> 为输入
	它将 c[i,j] 的值保存在表 c[0..m, 0..n] 中
	过程还维护一个表 b[1..m, 1..n],帮助构造最优解
	b[i,j] 指向的表项对应计算 c[i,j] 时所选择的子问题最优解
*/

// PRINT-LCS(b, X, i, j)		// i = X.length  j = Y.length
//	if i == 0 or j == 0
//		return
//	if b[i,j] == "↖"
//		PRINT-LCS(b, X, i-1, j-1)
//		print Xi
//	elseif b[i,j] == "↑"
//		PRINT-LCS(b, X, i-1, j)
//	else PRINT-LCS(b, X, i, j-1)
/*
	从 b[m,n] 开始,并按箭头方向追踪下去即可
	当在表项 b[i,j] 中遇到一个"↖"时,意味着 Xi = Yj 是一个LCS元素
*/


/* ················ 15.5 最有二叉搜索树 ···············
   用动态规划方法解决在已知关键字分布的前提下,如何构造最有二叉搜索树
*/