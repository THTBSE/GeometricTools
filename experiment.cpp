#include "experiment.h"
#include "../opengl/glew.h"

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
	size_t i = 0;
	for (const auto& v : Vectors)
	{
		glBegin(GL_LINES);
		if (i == 0 || i ==1)
		{
			glColor3f(0.5f, 0.5f, 0.0f);
			++i;
		}
		else
			glColor3f(1.0f, 1.0f, 1.0f);
		glVertex3dv(&v.first[0]);
		glVertex3dv(&v.second[0]);
		glEnd();
	}
	glDisable(GL_COLOR_MATERIAL);
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

	//if Cone3 if infinite , wo make height as 100.0 for convenient
	double height = c1.height == std::numeric_limits<double>::max() ? 100.0 : c1.height;
	double radius = fabs(height * c1.sinAngle / c1.cosAngle);
	gte::Vector3<double> center = c1.vertex + height * c1.axis;

	//how many vertices we need
	size_t times = 36;
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

std::vector<gte::Vector3<double>>
experiment::IntrConeToTriangle(const gte::Cone3<double>& cone, std::vector<gte::Triangle3<double>>& TList)
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
	for (auto i = 0; i < IndexList.size(); ++i)
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
			fsamples.push_back(intrl3.point);
		}
	}
	
	return std::move(fsamples);
}

std::vector<int>
experiment::GetOverhangTriangle(const std::vector<gte::Triangle3<double>>& TList,
double alphcC)
{
	auto VectorP = gte::Vector3<double>::Basis2();
	std::vector<int> ret;
	int triCount = (int)TList.size();
	for (int i = 0; i < triCount; ++i)
	{
		auto normal = gte::UnitCross(TList[i].v[1] - TList[i].v[0], TList[i].v[2] - TList[i].v[0]);
		double ang = gte::Angle(VectorP, normal);
		if (ang - GTE_C_HALF_PI > alphcC)
		{
			ret.push_back(i);
		}
	}
	return std::move(ret);
}