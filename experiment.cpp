#include "experiment.h"
#include "../opengl/glew.h"
#include <chrono>

void experiment::interstellar(std::array<gteVec4, 4>& orig, std::array<gteVec4, 4>& tar, gteVec3& atom)
{
	gte::Matrix4x4<double> M1(orig[0], orig[1], orig[2], orig[3], true);
	gte::Matrix4x4<double> M1Inverse = gte::Inverse<double>(M1);
	gte::Matrix4x4<double> M2(tar[0], tar[1], tar[2], tar[3], true);

	auto _atom = gte::HLift(atom, 1.0);
	_atom = M2 * M1Inverse * _atom;

	for (size_t i = 0; i < 3; ++i)
		atom[i] = _atom[i];
}

std::array<gte::Vector3<double>,4>
experiment::getCrdSystem(const gte::Vector3<double>& Vec, const gte::Vector3<double>& Vertex)
{
	gte::Vector3<double> z(Vec);
	gte::Normalize(z);
	gte::Vector3<double> x(-z[1], z[0], 0.0);
	gte::Normalize(x);
	gte::Vector3<double> y = gte::Cross(z, x);

	std::array<gte::Vector3<double>, 4> crdSys{ { x, y, z, Vertex } };

	return std::move(crdSys);
}

std::array<gte::Vector4<double>,4>
experiment::makeHomogeneous(const std::array<gteVec3, 4>& crdSys)
{
	std::array<gteVec4, 4> ret;
	for (size_t i = 0; i < crdSys.size(); ++i)
	{
		ret[i] = gte::HLift(crdSys[i], 0.0);
	}
	ret[3] = gte::HLift(crdSys[3], 1.0);
	return std::move(ret);
}

void experiment::drawCrdSystem()
{
	glEnable(GL_COLOR_MATERIAL);
	glLineWidth(5.0f);
	//for (size_t i = 0; i < crdSystems.size(); ++i)
	//{

	//}
	for (const auto& csys : crdSystems)
	{
		glColor3d(1.0, 0.0, 0.0);
		glBegin(GL_LINES);
		glVertex3dv(&csys[3][0]);
		gteVec3 xLine = csys[3] + 500.0 * csys[0];
		glVertex3dv(&xLine[0]);

		glColor3d(0.0, 1.0, 0.0);
		glVertex3dv(&csys[3][0]);
		gteVec3 yLine = csys[3] + 500.0 * csys[1];
		glVertex3dv(&yLine[0]);

		glColor3d(0.0, 0.0, 1.0);
		glVertex3dv(&csys[3][0]);
		gteVec3 zLine = csys[3] + 500.0 * csys[2];
		glVertex3dv(&zLine[0]);
		glEnd();
	}
	glDisable(GL_COLOR_MATERIAL);
}

void experiment::drawVectors()
{
	glEnable(GL_COLOR_MATERIAL);
	glLineWidth(5.0f);
	glColor3f(1.0f, 1.0f, 1.0f);
	for (const auto& v : Vectors)
	{
		glBegin(GL_LINES);
		glVertex3dv(&v.first[0]);
		glVertex3dv(&v.second[0]);
		glEnd();
	}
	glDisable(GL_COLOR_MATERIAL);
}

std::vector<gte::Segment3<double>>
experiment::DiscreteCone(const gte::Cone3<double>& c)
{
	std::vector<gte::Vector3<double>> diskPoints(GetConeDisk(c));
	std::vector<gte::Segment3<double>> segments;
	for (const auto &p : diskPoints)
	{
		segments.push_back(gte::Segment3<double>(c.vertex, p));
	}
	return std::move(segments);
}

std::vector<gte::Triangle3<double>>
experiment::DiscreteConeToMesh(const gte::Cone3<double>& c)
{
	std::vector<gte::Vector3<double>> diskPoints(GetConeDisk(c));
	std::vector<gte::Triangle3<double>> triangles;
	auto count = diskPoints.size();
	triangles.reserve(count);
	for (size_t i = 0; i < count - 1; ++i)
	{
		triangles.push_back(gte::Triangle3<double>(c.vertex, diskPoints[i], diskPoints[i + 1]));
	}
	triangles.push_back(gte::Triangle3<double>(c.vertex, diskPoints[count - 1], diskPoints[0]));
	return std::move(triangles);
}

std::vector<gte::Vector3<double>>
experiment::IntrConeToCone(std::vector<gte::Segment3<double>>& c1, const std::vector<gte::Triangle3<double>>& c2)
{
	gte::FIQuery<double, gte::Segment3<double>, gte::Triangle3<double>> fiq;
	std::vector<gte::Vector3<double>> intrp;
	auto c1V = c1[0].p[0];
	for (const auto &seg : c1)
	{
		double minDist = std::numeric_limits<double>::max();
		gte::Vector3<double> nearest;
		for (const auto &tri : c2)
		{
			auto ret = fiq(seg, tri);
			if (ret.intersect)
			{
				auto dist = gte::Length(ret.point - c1V);
				if (dist < minDist)
				{
					minDist = dist;
					nearest = ret.point;
				}
			}
		}
		if (minDist < std::numeric_limits<double>::max())
			intrp.push_back(nearest);
	}
	return std::move(intrp);
}

std::vector<gte::Vector3<double>>
experiment::IntrConeToCone(std::vector<gte::Segment3<double>>& c1, const gte::Cone3<double>& c2)
{
	std::vector<gte::Vector3<double>> ret;
	auto top = c1[0].p[0];
	auto top2 = c2.vertex;
	auto zx = -gte::Vector3<double>::Basis2();
	gte::Vector3<double> edgeVector,edgeVector2;
	double ang1, ang2;
	//we only take nearest intersect point
	gte::FIQuery<double, gte::Segment3<double>, gte::Cone3<double>> fiq;
	for (const auto &seg : c1)
	{
		auto Result = fiq(seg, c2);
		if (Result.intersect)
		{
			auto orig = gte::Angle(seg.p[1] - seg.p[0], zx);
			edgeVector = Result.point[0] - top;
			edgeVector2 = Result.point[0] - top2;
			ang1 = gte::Angle(edgeVector, zx);
			ang2 = gte::Angle(edgeVector2, zx);
			while (fabs(ang2 - ang1) > 0.01)
			{
				double sinTheta = sin(ang2)*c2.cosAngle - cos(ang2)*c2.sinAngle;
				double ll = sinTheta * gte::Length(Result.point[0] - top2);
				gte::Normalize(edgeVector);
				Result.point[0] = Result.point[0] + ll * edgeVector;
				ang2 = gte::Angle(Result.point[0] - top2, zx);

				bool xxyy = true;
			}
			ret.push_back(Result.point[0]);
		}
	}
	return std::move(ret);
}

//in this function , we deafult consider c1 and c2 are not infinite cone.
std::vector<gte::Vector3<double>>
experiment::IntrConeToCone(const gte::Cone3<double>& c1, const gte::Cone3<double>& c2)
{
	std::vector<gte::Vector3<double>> ret;
	if (c1.IsEqualTo(c2) || c1.height == std::numeric_limits<double>::max()
		|| c2.height == std::numeric_limits<double>::max())
		return ret;

	std::vector<gteVec3> diskPoints(GetConeDisk(c1));
	std::vector<gte::Segment3<double>> segments;
	for (const auto &p : diskPoints)
	{
		segments.push_back(gte::Segment3<double>(c1.vertex, p));
	}
	//we only take nearest intersect point
	for (const auto &seg : segments)
	{
		gte::FIQuery<double, gte::Segment3<double>, gte::Cone3<double>> fiq;
		auto Result = fiq(seg, c2);
		if (Result.intersect)
		{
			switch (Result.type)
			{
			case 1:case 2:
				ret.push_back(Result.point[0]);
				break;
			default:
				break;
			}
		}
	}
	return std::move(ret);
}

std::vector<gte::Vector3<double>>
experiment::IntrConeToConeInfinite(const gte::Cone3<double>& c1, const gte::Cone3<double>& c2)
{
	std::vector <gte::Vector3<double>> ret;

	std::vector<gteVec3> diskPoints(GetConeDisk(c1));
	std::vector<gte::Line3<double>> lines;
	for (const auto & p : diskPoints)
	{
		auto dir = p - c1.vertex;
		gte::Normalize(dir);
		lines.push_back(gte::Line3<double>(c1.vertex, dir));
	}
	
	//we only take nearest intersect point;
	for (const auto &line : lines)
	{
		gte::FIQuery<double, gte::Line3<double>, gte::Cone3<double>> fiq;
		auto Result = fiq(line, c2);
		if (Result.intersect)
		{
			switch (Result.type)
			{
			case 1: case 2:
				ret.push_back(Result.point[0]);
				break;
			//case 2:
			//{
			//	auto l1 = Result.point[0] - c1.vertex;
			//	auto l2 = Result.point[1] - c1.vertex;
			//	if (gte::Length(l1) < gte::Length(l2))
			//		ret.push_back(Result.point[0]);
			//	else
			//		ret.push_back(Result.point[1]);
			//}
			//	break;
			//case 3:
			//	ret.push_back(Result.point[0]);
			//	break;
			default:
				break;
			}
		}
	}
	return std::move(ret);
}


std::vector<gte::Vector3<double>>
experiment::GetConeDisk(const gte::Cone3<double>& c1)
{
	double crossProduct = gte::Length(gte::Cross(c1.axis, gteVec3::Basis2()));
	bool parallel = crossProduct < 1e-6;

	//if Cone3 is infinite , make height as 100.0 for convenient
	double height = c1.height == std::numeric_limits<double>::max() ? 100.0 : c1.height;
	double radius = fabs(height * c1.sinAngle / c1.cosAngle);
	gte::Vector3<double> center = c1.vertex + height * c1.axis;

	//how many vertices we need
	size_t times = 18;
	double marchAngle = GTE_C_TWO_PI / times;
	
	//construct transform matrix
	//if parallel , construct a translate matrix
	//else construct a frenet matrix
	//we do not care whether the points is matched precisely
	//because it's a circle.
	gte::Matrix4x4<double> TM;
	if (parallel)
		TM = gte::Matrix4x4<double>(gteVec4::Basis0(), gteVec4::Basis1(), gteVec4::Basis2(),gte::HLift(center,1.0), true);
	else
	{
		gteVec3 zAxis(-c1.axis[0], -c1.axis[1], -c1.axis[2]);
		gteVec3 xAxis(-zAxis[1], zAxis[0], 0.0);
		gte::Normalize(xAxis);
		gteVec3 yAxis = gte::UnitCross<3, double>(zAxis, xAxis);
		TM = gte::Matrix4x4<double>(gte::HLift(xAxis, 0.0), gte::HLift(yAxis, 0.0), gte::HLift(zAxis, 0.0), gte::HLift(center, 1.0), true);
	}

	std::vector<gteVec3> diskPoints;
	double ang = -marchAngle;
	for (size_t i = 0; i < times; ++i)
	{
		ang += marchAngle;
		auto v = gteVec4(radius*cos(ang), radius*sin(ang), 0.0, 1.0);
		v = TM * v;
		diskPoints.push_back(gteVec3(v[0],v[1],v[2]));
	}
	return std::move(diskPoints);
}


void experiment::IntrConeExperiment()
{
	gte::Vector3<double> v1(12.0, 30.5, 9.6), v2(3.0, 10.5, 8.6);
	auto axis = -gte::Vector3<double>(0.0, 1.0, 0.0);
	gte::Normalize(axis);
	gte::Cone3<double> c1(v1, axis, GTE_C_QUARTER_PI, 12.0);
	gte::Cone3<double> c2(v2, axis, GTE_C_QUARTER_PI, 12.0);

	auto center = c1.vertex + c1.height * c1.axis;
	auto center2 = c2.vertex + c2.height * c2.axis;
	Vectors.push_back(std::make_pair(c1.vertex, center));
	Vectors.push_back(std::make_pair(c2.vertex, center2));
	auto disk1 = GetConeDisk(c1);
	auto disk2 = GetConeDisk(c2);

	for (const auto& p : disk1)
	{
		Vectors.push_back(std::make_pair(c1.vertex, p));
	}
	for (const auto& p : disk2)
	{
		Vectors.push_back(std::make_pair(c2.vertex, p));
	}

	auto intrp = IntrConeToCone(c1, c2);
	points.assign(intrp.begin(), intrp.end());
}

void experiment::drawPoints()
{
	glEnable(GL_COLOR_MATERIAL);
	glPointSize(10.0f);
//	glLineWidth(10.0f);
	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_POINTS);
	for (const auto& v : points)
	{
		glVertex3dv(&v[0]);
	}
	glEnd();
	glDisable(GL_COLOR_MATERIAL);
}

//CAN BE ABANDONED!!
std::vector<gte::Vector3<double>>
experiment::IntrConeToTriangle(const std::vector<gte::Segment3<double>>& cone,
const std::vector<gte::Triangle3<double>>& TList)
{
	std::vector<gte::Vector3<double>> ret;
	gte::FIQuery<double, gte::Segment3<double>, gte::Triangle3<double>> fiq;
	double minDist = std::numeric_limits<double>::max();
	auto coneVertex = cone[0].p[0];
	for (const auto &seg : cone)
	{
		for (const auto &tri : TList)
		{
			auto Result = fiq(seg, tri);
			if (Result.intersect)
			{
				double dist = gte::Length(coneVertex - Result.point);
				if (dist < minDist)
				{
					minDist = dist;
					ret.push_back(Result.point);
				}
			}
		}
	}
	return std::move(ret);
}

//CAN BE ABADONED!!
std::vector<gte::Vector3<double>>
experiment::IntrConeToTriangle(const gte::Cone3<double>& cone, const std::vector<gte::Triangle3<double>>& TList)
{
	std::vector<gteVec3> diskPoints(GetConeDisk(cone));
	std::vector<gte::Segment3<double>> segments;
	for (const auto &p : diskPoints)
	{
		segments.push_back(gte::Segment3<double>(p, cone.vertex));
	}

	std::vector<gte::Vector3<double>> ret;
	for (const auto &seg : segments)
	{
		//bug bug bug !! first intersection may not be nearest !!
		for (const auto &tri : TList)
		{
			gte::FIQuery<double, gte::Segment3<double>, gte::Triangle3<double>> fiq;
			auto Result = fiq(seg, tri);
			if (Result.intersect)
			{
				ret.push_back(Result.point);
				break;
			}
		}
	}
	return std::move(ret);
}

std::vector<gte::Triangle3<double>>
experiment::GetTriangleList(TriMesh* mesh)
{
	if (!mesh)
		return std::vector<gte::Triangle3<double>>();

	auto TList = std::vector<gte::Triangle3<double>>();
	std::vector<TriMesh::Face> &faces = mesh->faces;
	std::vector<vec> &vtx = mesh->vertices;
	for (const auto & f : faces)
	{
		std::array<gteVec3, 3> vertex;
		for (size_t i = 0; i < 3; i++)
		{
			vertex[i] = gteVec3(vtx[f[i]][0], vtx[f[i]][1], vtx[f[i]][2]);
		}
		TList.push_back(gte::Triangle3<double>(vertex[0], vertex[1], vertex[2]));
	}
	return std::move(TList);
}

std::vector<gte::Vector3<double>>
experiment::GetSupportPoint(const std::vector<int>& IndexList, 
		const std::vector<gte::Triangle3<double>>& TriList, double resolution)
{
	//Firstly,Projecting triangles to XOY plane.
	typedef std::pair<gte::Triangle2<double>, size_t> Tri2Index;
	std::vector<Tri2Index> xyTriList;
	double xBeg(std::numeric_limits<double>::max()), xEnd(std::numeric_limits<double>::min());
	double ymin(std::numeric_limits<double>::max());
	for (size_t i = 0; i < IndexList.size(); ++i)
	{
		gte::Vector2<double> tri[3];
		for (size_t j = 0; j < 3; j++)
		{
			tri[j] = gte::HProject(TriList[IndexList[i]].v[j]);
		}
		double xbegTemp = (*std::min_element(std::begin(tri), std::end(tri), CompOneDim<2>(0)))[0];
		double xendTemp = (*std::max_element(std::begin(tri), std::end(tri), CompOneDim<2>(0)))[0];
		xBeg = xbegTemp < xBeg ? xbegTemp : xBeg;
		xEnd = xEnd < xendTemp ? xendTemp : xEnd;
		
		double yminTemp = (*std::min_element(std::begin(tri), std::end(tri), CompOneDim<2>(1)))[1];
		ymin = yminTemp < ymin ? yminTemp : ymin;

		xyTriList.push_back(std::make_pair(gte::Triangle2<double>(tri),IndexList[i]));
	}

	//Secondly,Getting intersection pairs from xBeg to xEnd;
	std::map<double, std::vector<IntrPair>> Lines;
	gte::FIQuery<double, gte::Line2<double>, gte::Triangle2<double>> fiq;
	while (xBeg < xEnd)
	{
		std::vector<IntrPair> IPair;
		gte::Line2<double> xScanline(gte::Vector2<double>(xBeg, ymin), gte::Vector2<double>::Basis1());
		for (size_t i = 0; i < xyTriList.size(); ++i)
		{
			auto ret = fiq(xScanline, xyTriList[i].first);
			if (ret.intersect)
			{
				if (ret.numIntersections == 1)
				{
					IPair.push_back(IntrPair(ret.point[0][1], ret.point[0][1], xyTriList[i].second));
				}
				else if (ret.numIntersections == 2)
				{
					IPair.push_back(IntrPair(ret.point[0][1], ret.point[1][1], xyTriList[i].second));
				}
			}
		}
		//It should be sort 
		if (!IPair.empty())
		{
			std::sort(IPair.begin(), IPair.end());
			Lines.insert(std::make_pair(xBeg, IPair));
		}
		xBeg += resolution;
	}

	//Thirdly,Getting sampling points in XOY plane
	typedef std::pair<gte::Vector2<double>, size_t> Vec2Index;
	std::vector<Vec2Index> samples;
	for (auto LIter = Lines.begin(); LIter != Lines.end(); ++LIter)
	{
		std::vector<IntrPair>& line = LIter->second;
		double yBeg = line.front().first;
		const double yEnd = line.back().second;
		const double xNum = LIter->first;
		std::vector<IntrPair>::size_type curr = 0;
		while (yBeg < yEnd)
		{
			if (line[curr].first <= yBeg && yBeg <= line[curr].second)
			{
				samples.push_back(std::make_pair(gte::Vector2<double>(xNum, yBeg), line[curr].TriIndex));
			}
			else if (line[curr].second < yBeg)
			{
				++curr;
				continue;
			}
			yBeg += resolution;
		}
	}

	//Fourthly, Using Line intersection to get final sampling points
	std::vector<gte::Vector3<double>> fsamples;
	gte::FIQuery<double, gte::Line3<double>, gte::Triangle3<double>> l3fiq;
	for (const auto &xoyp : samples)
	{
		auto v3 = gte::HLift(xoyp.first, 0.0);
		gte::Line3<double> l3(v3, gte::Vector3<double>::Basis2());
		auto intrl3 = l3fiq(l3, TriList[xoyp.second]);
		if (intrl3.intersect)
		{
			//for generating the support-connect structure,should minus a small value.
			intrl3.point[2] -= 0.1;
			fsamples.push_back(intrl3.point);
		}
	}
	
	return std::move(fsamples);
}

std::vector<int>
experiment::GetOverhangTriangle(const std::vector<gte::Triangle3<double>>& TList,
double alphcC)
{
	double epsilon = 0.1;
	auto VectorP = gte::Vector3<double>::Basis2();
	std::vector<int> ret;
	int triCount = (int)TList.size();
	for (int i = 0; i < triCount; ++i)
	{
		auto dist = (TList[i].v[0][2] + TList[i].v[1][2] + TList[i].v[2][2]) / 3.0;
		if (dist < epsilon)
			continue;
		auto normal = gte::UnitCross(TList[i].v[1] - TList[i].v[0], TList[i].v[2] - TList[i].v[0]);
		double ang = gte::Angle(VectorP, normal);
		if (ang - GTE_C_HALF_PI > alphcC)
		{
			ret.push_back(i);
		}
	}
	return std::move(ret);
}

typedef std::chrono::time_point<std::chrono::system_clock> sys_clock;
typedef std::chrono::duration<double> used_time;

std::vector<gte::Segment3<double>>
experiment::GenerateSupport(const std::vector<gte::Vector3<double>>& InitialSet,
const std::vector<gte::Triangle3<double>>& TList,
double AlphaC)
{
	using namespace std;
	typedef gte::Vector3<double> GVec3;
	typedef pair<GVec3, multiset<Support>::iterator> GVec3AndPtr;

	//build octree for computing intersection between segment and triangle mesh faster
	TOctree TriOctree(TList, aabb, 9);

	multiset<Support> PointSet;
	size_t Identifier = 0;
	for (const auto & p : InitialSet)
	{
		PointSet.insert(Support(p, Identifier++, AlphaC));
	}

	vector<gte::Segment3<double>> branch;
	while (!PointSet.empty())
	{
		auto top = PointSet.begin();
		auto c1 = DiscreteCone(top->cone());
		
		//compute intersection between c1 to all other cones.
		//IntrPts[i].first is the intersect point.
		//IntrPts[i].second is the iterator of ci.
		vector<GVec3AndPtr> IntrPts;
		for (auto it = ++PointSet.begin(); it != PointSet.end(); ++it)
		{
			//maybe empty because non intersection between two cones!
		//	auto intrs = IntrConeToCone(c1, it->cone());
			auto intrs = IntrConeToCone(c1, DiscreteConeToMesh(it->cone()));
			if (intrs.empty())
				continue;
			auto max = *max_element(intrs.begin(), intrs.end(), CompOneDim<3>(2));
			IntrPts.push_back(make_pair(max, it));
		}
		double dConeCone = std::numeric_limits<double>::max();
		//find the nearest intersection to c1
		auto IntrIter = min_element(IntrPts.begin(), IntrPts.end(), 
			[&top](const GVec3AndPtr& lhs, const GVec3AndPtr& rhs) ->bool
		{
			return gte::Length(lhs.first - top->point()) < gte::Length(rhs.first - top->point());
		});

		if (IntrIter != IntrPts.end())
			dConeCone = gte::Length(IntrIter->first - top->point());

		//compute intersection between c1 to triangle mesh by Octree
		double dConeMesh = std::numeric_limits<double>::max();
		std::vector<GVec3> intrm;
		for (const auto &seg : c1)
		{
			auto ORet = TriOctree.IntrQuery(seg);
			if (ORet.intersect)
				intrm.push_back(ORet.point);
		}
		auto IntrMeshIter = max_element(intrm.begin(), intrm.end(), CompOneDim<3>(2));
		if (IntrMeshIter != intrm.end())
			dConeMesh = gte::Length(*IntrMeshIter - top->point());

		//compute the distance from support to printing ground plane (XOY Plane).
		double dPointPlane = top->point()[2];

		double minDist = std::min(std::min(dConeCone, dConeMesh), dPointPlane);
		if (minDist == dConeCone)
		{
			//c1's point to intersection and ci's point to intersection.
			//then erase Support c1 and ci from PointSet.
			branch.push_back(gte::Segment3<double>(top->point(), IntrIter->first));
			branch.push_back(gte::Segment3<double>(IntrIter->second->point(), IntrIter->first));
			PointSet.erase(IntrIter->second);
			PointSet.insert(Support(IntrIter->first, Identifier++, AlphaC));
		}
		else if (minDist == dConeMesh)
		{
			branch.push_back(gte::Segment3<double>(top->point(), *IntrMeshIter));
		}
		else
		{
			branch.push_back(gte::Segment3<double>(top->point(), GVec3(top->point()[0],top->point()[1],0.0)));
		}
		PointSet.erase(PointSet.begin());
	}

	return std::move(branch);
}

gte::AlignedBox3<double>
experiment::GenerateAABB(TriMesh *mesh)
{
	if (mesh == NULL)
		return gte::AlignedBox3<double>();

	float minX = mesh->vertices[0][0], maxX = minX;
	float minY = mesh->vertices[0][1], maxY = minY;
	float minZ = mesh->vertices[0][2], maxZ = minZ;

	for (size_t i = 1; i < mesh->vertices.size(); ++i)
	{
		if (mesh->vertices[i][0] < minX)
			minX = mesh->vertices[i][0];
		else if (mesh->vertices[i][0] > maxX)
			maxX = mesh->vertices[i][0];

		if (mesh->vertices[i][1] < minY)
			minY = mesh->vertices[i][1];
		else if (mesh->vertices[i][1] > maxY)
			maxY = mesh->vertices[i][1];

		if (mesh->vertices[i][2] < minZ)
			minZ = mesh->vertices[i][2];
		else if (mesh->vertices[i][2] > maxZ)
			maxZ = mesh->vertices[i][2];
	}

	double min = std::min(minX, std::min(minY, minZ));
	double max = std::max(maxX, std::max(maxY, maxZ));

	return gte::AlignedBox3<double>(min, min, min, max, max, max);
}

void test_read_time()
{
	sys_clock start0, start1, start2, end0, end1, end2;

	start0 = std::chrono::system_clock::now();
	TriMesh *mesh0 = TriMesh::load_stl3("dragon.stl");
	end0 = std::chrono::system_clock::now();
	delete mesh0;

	start1 = std::chrono::system_clock::now();
	TriMesh *mesh1 = TriMesh::load_stl2("dragon.stl");
	end1 = std::chrono::system_clock::now();
	delete mesh1;

	start2 = std::chrono::system_clock::now();
	TriMesh *mesh2 = TriMesh::load_stl("dragon.stl");
	end2 = std::chrono::system_clock::now();
	delete mesh2;

	std::chrono::duration<double> time0, time1, time2;
	time0 = end0 - start0;
	time1 = end1 - start1;
	time2 = end2 - start2;

	double t0, t1, t2;
	t0 = time0.count();
	t1 = time1.count();
	t2 = time2.count();

	FILE* f = fopen("time.txt", "w");
	fprintf(f, "%f,%f,%f", t0, t1, t2);
	fclose(f); 
}

std::vector<gte::Vector3<double>>
experiment::GetOverhangPoints(TriMesh *mesh)
{
	std::vector<gte::Vector3<double>> ret;
	mesh->need_neighbors();
	size_t vcount = mesh->vertices.size();
	size_t nbIndex, ncount;
	bool isOverhang;
	for (size_t i = 0; i < vcount; ++i)
	{
		if (mesh->vertices[i][2] < 1)
			continue;
		ncount = mesh->neighbors[i].size();
		isOverhang = true;
		for (size_t j = 0; j < ncount; ++j)
		{
			nbIndex = mesh->neighbors[i][j];
			if (mesh->vertices[i][2] > mesh->vertices[nbIndex][2])
			{
				isOverhang = false;
				break;
			}
		}
		if (isOverhang)
		{
			ret.push_back(gte::Vector3<double>((double)mesh->vertices[i][0],
				(double)mesh->vertices[i][1], (double)mesh->vertices[i][2] - 0.1));
		}
	}
	return std::move(ret);
}

TriMesh*
experiment::GenerateSupportMesh(const std::vector<gte::Segment3<double>>& segs)
{
	TriMesh *mesh = new TriMesh;
	for (const auto &s : segs)
	{
		auto rad = GetStrutRad(s);
		point p1(s.p[0][0], s.p[0][1], s.p[0][2]);
		point p2(s.p[1][0], s.p[1][1], s.p[1][2]);
		GenerateStrut(rad, p2, p1, mesh->vertices, mesh->faces);
	}
	return mesh;
}

void
experiment::GenerateStrut(double rad, const point &p1, const point &p2,
vector<point> &vertices, vector<TriMesh::Face> &faces)
{
	vec medial = p2 - p1;
	float height = len(medial);
	int beginT = vertices.size();
	vector<int> vid;

	//generate base vertices
	vector<vec> basev;
	int numpv = 14;
	for (int i = 0; i<numpv; i++)
	{
		float alpha = (float)i* 2 * M_PIf / numpv;
		basev.push_back(vec(rad*cos(alpha), rad*sin(alpha), 0));
	}
	vertices.insert(vertices.end(), basev.begin(), basev.end());

	float incre = height / 3;
	for (int i = 0; i<3; i++)
	{
		for (unsigned j = 0; j<basev.size(); j++)
		{
			basev[j][2] += incre;
		}
		vertices.insert(vertices.end(), basev.begin(), basev.end());
	}

	int endT = vertices.size();
	for (int i = beginT; i<endT; i++)
	{
		vid.push_back(i);
	}

	//generate side faces
	for (int i = 0; i<3; i++)
	{
		int p = i * numpv;
		int ps = (i + 1) * numpv;
		for (int j = 0; j<numpv; j++)
		{
			TriMesh::Face f;
			f[0] = vid[j + p];
			f[1] = vid[ps + (j + 1) % numpv];
			f[2] = vid[j + ps];
			faces.push_back(f);

			f[0] = vid[j + p];
			f[1] = vid[p + (j + 1) % numpv];
			f[2] = vid[ps + (j + 1) % numpv];
			faces.push_back(f);
		}
	}

	vec bc(0, 0, 0), tc(0, 0, height);
	vertices.push_back(bc);
	int vbc = vertices.size() - 1;
	vertices.push_back(tc);
	int vtc = vertices.size() - 1;
	endT = vertices.size();

	//generate bottom & top faces
	for (int i = 0; i<numpv; i++)
	{
		TriMesh::Face f;
		f[0] = vbc;
		f[1] = vid[(i + 1) % numpv];
		f[2] = vid[i];
		faces.push_back(f);
	}

	int lst = 3 * numpv;
	int end = 4 * numpv;
	for (int i = lst; i<end; i++)
	{
		TriMesh::Face f;
		f[0] = vtc;
		f[1] = vid[i];
		if (i != (end - 1))
			f[2] = vid[i + 1];
		else
			f[2] = vid[lst];

		faces.push_back(f);
	}

	//transform to correct position 
	tc = normalize(tc); medial = normalize(medial);
	gte::Matrix4x4<double> TM;
	TM.MakeIdentity();
	float ang = angle(tc, medial);
	if (fabs(ang - M_PIf) < 1e-6)
	{
		TM.SetCol(3, gte::Vector4<double>(p2[0], p2[1], p2[2], 1.0));
	}
	else if (fabs(ang) < 1e-6)
	{
		TM.SetCol(3, gte::Vector4<double>(p1[0], p1[1], p1[2], 1.0));
	}
	else
	{
		gteVec3 zAxis(medial[0], medial[1], medial[2]);
		gteVec3 xAxis(-zAxis[1], zAxis[0], 0.0);
		gte::Normalize(xAxis);
		gteVec3 yAxis = gte::UnitCross<3, double>(zAxis, xAxis);
		TM = gte::Matrix4x4<double>(gte::HLift(xAxis, 0.0), gte::HLift(yAxis, 0.0), gte::HLift(zAxis, 0.0),
			gteVec4(p1[0],p1[1],p1[2],1.0), true);
	}

	for (int i = beginT; i<endT; i++)
	{
		gteVec4 tv(vertices[i][0], vertices[i][1], vertices[i][2], 1.0);
		auto tmp = TM * tv;
		vertices[i][0] = tmp[0];
		vertices[i][1] = tmp[1];
		vertices[i][2] = tmp[2];
	}

}