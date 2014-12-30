#ifndef _SUPPORTTREE_H_
#define _SUPPORTTREE_H_
#include <memory>
#include <vector>
#include "../Mathematics/GteVector3.h"
using std::shared_ptr;
using std::vector;


class STreeNode
{
public:
	STreeNode(int lv, const gte::Vector3<double>& pt, bool root);

	bool operator< (const STreeNode &rhs) const;
	void Init();

	bool isRoot;
	int level;
	double length;
	double rad;
	gte::Vector3<double> point;
	//It's Binary Tree
	shared_ptr<STreeNode> children[2];
};





#endif