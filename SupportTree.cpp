#include "SupportTree.h"

STreeNode::STreeNode(int lv, const gte::Vector3<double> &pt, bool root) :level(lv), point(pt), isRoot(root)
{

}

bool STreeNode::operator< (const STreeNode &rhs) const
{
	return level < rhs.level;
}
void
STreeNode::Init()
{
	if (level == 0)
	{
		return;
	}

	if (isRoot)
	{
		if (level == 1)
		{
			rad = 0.8 * children[0]->rad;
			length = children[0]->length;
		}
		else
		{
			auto l = gte::Length(point - children[0]->point);
			length = l + children[0]->length;
			auto ratio = length / children[0]->length;
			ratio = ratio > 1.3 ? ratio : 1.3;
			rad = children[0]->rad * ratio;
		}
		if (rad < 0 || length < 0)
		{
			int x = 1;
		}
		return;

	}

	if (level == 1)
	{
		rad = 0.8 * (children[0]->rad + children[1]->rad);
		length = children[0]->length + children[1]->length;

		if (rad < 0 || length < 0)
		{
			int x = 1;
		}
	}
	else 
	{
		auto l0 = gte::Length(point - children[0]->point);
		auto l1 = gte::Length(point - children[1]->point);
		auto childrad = 0.5 * (children[0]->rad + children[1]->rad);
		auto childlength = children[0]->length + children[1]->length;
		length = l0 + l1 + childlength;
		rad = childrad * (length / childlength);

		if (rad < 0 || length < 0)
		{
			int x = 1;
		}
	}
}