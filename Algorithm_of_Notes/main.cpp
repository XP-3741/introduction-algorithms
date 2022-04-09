#include"Algorithm _Implementation.h"

int main()
{
	/*----∑÷÷Œ∑®_≤‚ ‘----*/
	/*ElemType A[ARRAY_SIZE] = { 9,5,7,15,4,9,11,1,3,5 };
	MERGE_SORT(A, 0, ARRAY_SIZE - 1);
	for (int k : A)
		std::cout << k << " ";*/

	/*----◊Ó¥Û◊” ˝◊È_≤‚ ‘----*/
	//ElemType A[8] = { 1,-2, 3, 10, -4, 7 ,2, -5};
	////array_return_max result = FIND_MAXIMUM_SUBARRAY(A, 0, 8 - 1);
	//array_return_max result = FIND_MAXIMUM_SUBARRAY_LINEAR(A, 8);
	//cout << "max_left_index: " << result.max_left_index << endl;
	//cout << "max_right_index: " << result.max_right_index << endl;
	//cout << "sum: " << result.sum << endl;

	/*----∂—≈≈–Ú_≤‚ ‘----*/
	//// 16 14 10 8 7 9 3 2 4 1
	//int test_length = 0;
	//cout << "Enter array size: ";
	//cin >> test_length;
	//heap_array A(test_length);
	//cout << "Enter array elements: " << endl;
	//int index = 1;
	//while (cin >> A[index] && index < 10)
	//	index++;
	//HEAPSORT(A);
	//cout << "HEAPSORT result: " << endl;
	//for (int i = 1; i <= 10; i++)
	//	cout << A[i] << " ";
	//cout << endl << "Done" << endl;

	/*----øÏÀŸ≈≈–Ú_≤‚ ‘----*/
	/*ElemType A[8] = { 1,-2, 3, 10, -4, 7 ,2, -5 };
	QUICKSORT(A, 0, 7);
	for (int k : A)
		cout << k << " ";
	cout << endl;*/

	/*----º∆ ˝≈≈–Ú_≤‚ ‘----*/
	/*ElemType A[ARRAY_SIZE] = { 5, 2, 4, 1, 4, 7 ,2, 5, 9, 1 };
	ElemType* B = new ElemType[ARRAY_SIZE + 1];
	ElemType k = 0;
	for (int i = 0; i < ARRAY_SIZE; i++)
		if (k < A[i])
			k = A[i];
	COUNTING_SORT(A, B, k);
	for (int j = 1; j < ARRAY_SIZE + 1; j++)
		cout << B[j] << " ";
	cout << endl;*/

	/*----“≈Õ¸±»ΩœΩªªªÀ„∑®_≤‚ ‘----*/
	/*ElemType A[ARRAY_SIZE] = { 5, 2, 4, 1, 4, 7 ,2, 5, 9, 1 };
	INSERTION_SORT(A);
	for (int k : A)
		cout << k << " ";
	cout << endl;*/

	/*----∂˛≤ÊÀ—À˜ ˜_≤‚ ‘----*/
	/*srand((int)time(0));
	TreeNode* TRoot = new TreeNode(15);
	for (int i = 0; i < 10; i++) {
		TreeNode* A = new TreeNode(rand()%100);
		TREE_INSERT(TRoot, A);
	}
	cout << "Print all value: ";
	INORDER_TREE_WALK(TRoot);
	cout << endl;
	TreeNode* findNode = TREE_SEARCH(TRoot, 62);
	if (findNode != nullptr)
		cout << "Find " << findNode->val << " success" << endl;
	else
		cout << "Not find 62" << endl;
	cout << "The minimum value: " << TREE_MINIMUM(TRoot)->val << endl;
	cout << "The maximum value: " << TREE_MAXIMUM(TRoot)->val << endl;
	TreeNode* findSuccessor = TRoot->right->left;
	cout << "The " << findSuccessor->val << "'s successor is " << TREE_SUCCESSOR(findSuccessor)->val << endl;
	TreeNode* insertNode = new TreeNode(55);
	TREE_INSERT(TRoot, insertNode);
	cout << "After insert 55 then print all value: ";
	INORDER_TREE_WALK(TRoot);
	cout << endl;
	TREE_DELETE(TRoot, insertNode);
	cout << "After delete 55 then print all value: ";
	INORDER_TREE_WALK(TRoot);
	cout << endl;*/

	return 0;
}