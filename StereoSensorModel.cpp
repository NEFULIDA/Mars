#pragma once
#include "StereoSensorModel.h"

StereoSensorModel::StereoSensorModel(string sAttPathNd, string sEtPathNd, string sEpPathNd, string sInnerParaNd, string sOffsetParaNd, int iCameraIndexNd, bool existCalibND, \
	string sAttPathS1, string sEtPathS1, string sEpPathS1, string sInnerParaS1, string sOffsetParaS1, int iCameraIndexS1, bool existCalibS1, \
	string sAttPathS2, string sEtPathS2, string sEpPathS2, string sInnerParaS2, string sOffsetParaS2, int iCameraIndexS2, bool existCalibS2, \
	string sCamerainstallPath):m_cModelNd(sAttPathNd, sEtPathNd, sEpPathNd, sCamerainstallPath,sInnerParaNd,sOffsetParaNd, iCameraIndexNd, existCalibND), \
	m_cModelS1(sAttPathS1, sEtPathS1, sEpPathS1, sCamerainstallPath,sInnerParaS1,sOffsetParaS1, iCameraIndexS1, existCalibS1),\
	m_cModelS2(sAttPathS2, sEtPathS2, sEpPathS2, sCamerainstallPath,sInnerParaS2,sOffsetParaS2, iCameraIndexS2, existCalibS2){
}

StereoSensorModel::~StereoSensorModel()
{
	
}

bool StereoSensorModel::initImg(string sImgNd, string sImgS1, string sImgS2)
{
	m_cModelNd.initImg(sImgNd);
	m_cModelS1.initImg(sImgS1);
	m_cModelS2.initImg(sImgS2);
	return true;
}

bool StereoSensorModel::initDem(string sDEMfile)
{
	m_dem.open(sDEMfile);
	int iImgWidthDEM = m_dem.getImageWid();
	int iImgHeightDEM = m_dem.getImageHei();
	double trans[4];
	m_dem.getGeoTransform(trans);
	m_dem.readToBuffer_proj(0, trans[0], trans[0] + iImgWidthDEM * trans[1], trans[2], trans[2] + iImgHeightDEM * trans[3]);	//整幅图读入缓存区间，加速后续处理
	m_coord.init(GcsMarsDTMWTK.c_str(), m_dem.getProjectionRef().c_str());	//	初始化转换关系

	//每个相机模型均初始化一次DEM
	m_cModelNd.initDem(sDEMfile);
	m_cModelS1.initDem(sDEMfile);
	m_cModelS2.initDem(sDEMfile);

	m_cModelNd.m_Dem_Coord.init(GcsMarsDTMWTK.c_str(), m_dem.getProjectionRef().c_str());
	m_cModelS1.m_Dem_Coord.init(GcsMarsDTMWTK.c_str(), m_dem.getProjectionRef().c_str());
	m_cModelS2.m_Dem_Coord.init(GcsMarsDTMWTK.c_str(), m_dem.getProjectionRef().c_str());
	
	return true;
}

GeoRange StereoSensorModel::transtoGeoRange(vector<rpcGCP>& vGeoRange)
{
	GeoRange stereoGeoRange = { vGeoRange[0].lat,vGeoRange[0].lon,\
		make_pair(vGeoRange[1].lat - vGeoRange[0].lat, vGeoRange[1].lon - vGeoRange[0].lon),\
		make_pair(vGeoRange[2].lat - vGeoRange[0].lat, vGeoRange[2].lon - vGeoRange[0].lon) };
	return stereoGeoRange;
}

GeoRange StereoSensorModel::calOverlapArea2Img(GeoRange& model1Range, GeoRange& model2Range)
{
	////将角点地理坐标转为起点和向量的表示形式,默认覆盖区域均为平行四边形
	//GeoRange model1Range = { vModel1Range[0].lat,vModel1Range[0].lon,\
	//	make_pair(vModel1Range[1].lat - vModel1Range[0].lat,vModel1Range[1].lon - vModel1Range[0].lon),\
	//	make_pair(vModel1Range[2].lat - vModel1Range[0].lat,vModel1Range[2].lon - vModel1Range[0].lon) };
	//GeoRange model2Range = { vModel2Range[0].lat,vModel2Range[0].lon,\
	//	make_pair(vModel2Range[1].lat - vModel2Range[0].lat, vModel2Range[1].lon - vModel2Range[0].lon), \
	//	make_pair(vModel2Range[2].lat - vModel2Range[0].lat, vModel2Range[2].lon - vModel2Range[0].lon) };

	//模型一到模型二的起点指向
	//pair<double, double> vectormodel1model2(vModel2Range[0].lat - vModel1Range[0].lat, vModel2Range[0].lon - vModel1Range[0].lon);	
	pair<double, double> vectormodel1model2(model2Range.lat - model1Range.lat, model2Range.lon - model1Range.lon);	
	Eigen::Vector2d vector11 = { vectormodel1model2.first,vectormodel1model2.second };

	//两个模型的向量指向
	Eigen::Vector2d vector12 = { model1Range.vector12.first ,model1Range.vector12.second };
	Eigen::Vector2d vector13 = { model1Range.vector13.first ,model1Range.vector13.second };
	
	Eigen::Vector2d vector22 = { model2Range.vector12.first ,model2Range.vector12.second };
	Eigen::Vector2d vector23 = { model2Range.vector13.first ,model2Range.vector13.second };

	//计算 模型一的向量与起点指向向量的夹角
	double cosVal12 = vector11.dot(vector12) / (vector11.norm() * vector12.norm()); //角度cos值
	double angle12 = acos(cosVal12);     //弧度角
	double cosVal13 = vector11.dot(vector13) / (vector11.norm() * vector13.norm());
	double angle13 = acos(cosVal13); 
	double cosVal23 = vector12.dot(vector13) / (vector12.norm() * vector13.norm());
	double angle23 = acos(cosVal23);

	//以模型一起点为原点 根据夹角分为四种情况
	GeoRange overlapArea;
	// 1 模型二的起点在第二象限 在模型一内部 
	if (abs(angle12 + angle13 - angle23) < 0.001)
	{
		Eigen::Vector2d newvector12 = vector12 + vector11.norm() / sin(angle23) * sin(angle12) / vector13.norm() * vector13 - vector11;
		Eigen::Vector2d newvector13 = vector13 + vector11.norm() / sin(angle23) * sin(angle13) / vector12.norm() * vector12 - vector11;
		overlapArea = { model2Range.lat,model2Range.lon,\
			make_pair(newvector12.x(),newvector12.y()),\
			make_pair(newvector13.x(),newvector13.y()) };
	}
	// 2 模型二的起点在第四象限 
	if (abs(angle12 + angle13 + angle23 - 2*M_PI) < 0.001)
	{
		Eigen::Vector2d newvector12 = vector11 + vector22 + vector11.norm() / sin(angle23) * sin(angle12) / vector13.norm() * vector13;
		Eigen::Vector2d newvector13 = vector11 + vector23 + vector11.norm() / sin(angle23) * sin(angle13) / vector12.norm() * vector12;
		overlapArea = { model1Range.lat,model1Range.lon,\
			make_pair(newvector12.x(),newvector12.y()),\
			make_pair(newvector13.x(),newvector13.y()) };
	}
	// 3 模型二的起点在第一象限 
	if (abs(angle13 + angle23 - angle12) < 0.001)
	{
		//重叠区域起点
		double orix, oriy;
		Eigen::Vector2d orixyvector = vector11.norm() / sin(angle23) * sin(angle12) / vector13.norm() * vector13;
		orix = orixyvector.x() + model1Range.lat;
		oriy = orixyvector.y() + model1Range.lon;

		Eigen::Vector2d newvector12 = vector22 - orixyvector + vector11;
		Eigen::Vector2d newvector13 = vector13 - orixyvector;
		overlapArea = { orix, oriy,\
			make_pair(newvector12.x(),newvector12.y()),\
			make_pair(newvector13.x(),newvector13.y()) };
	}
	// 4 模型二的起点在第三象限 
	if (abs(angle12 + angle23 - angle13) < 0.001)
	{
		//重叠区域起点
		double orix, oriy;
		Eigen::Vector2d orixyvector = vector11.norm() / sin(angle23) * sin(angle13) / vector12.norm() * vector12;
		orix = orixyvector.x() + model1Range.lat;
		oriy = orixyvector.y() + model1Range.lon;

		Eigen::Vector2d newvector12 = vector12 - orixyvector;
		Eigen::Vector2d newvector13 = vector23 - orixyvector + vector11;
		overlapArea = { orix, oriy,\
			make_pair(newvector12.x(),newvector12.y()),\
			make_pair(newvector13.x(),newvector13.y()) };
	}

	return overlapArea;
}

GeoRange StereoSensorModel::calOverlapArea3Img()
{
	printf("\n计算重叠区域: ");
	//首先计算S1 S2 模型的重叠区域
	vector<rpcGCP> vS1ModelRange = m_cModelS1.getGeoRange(0, 0, m_cModelS1.m_Img_Width, m_cModelS1.m_Img_Height, true);
	vector<rpcGCP> vS2ModelRange = m_cModelS2.getGeoRange(0, 0, m_cModelS2.m_Img_Width, m_cModelS2.m_Img_Height, true);
	GeoRange S1GeoRange = transtoGeoRange(vS1ModelRange);
	GeoRange S2GeoRange = transtoGeoRange(vS2ModelRange);
	GeoRange overlapAreaS1S2 = calOverlapArea2Img(S1GeoRange, S2GeoRange);
	
	//计算重叠区域与ND模型的重叠区域
	vector<rpcGCP> vNdModelRange = m_cModelNd.getGeoRange(0, 0, m_cModelNd.m_Img_Width, m_cModelNd.m_Img_Height, true);
	GeoRange NdGeoRange = transtoGeoRange(vNdModelRange);
	GeoRange overlapArea = calOverlapArea2Img(overlapAreaS1S2, NdGeoRange);

	return overlapArea;
}

bool StereoSensorModel::segmentGrid(OpticalLinearSensorModel &model1, OpticalLinearSensorModel &model2, GeoRange &range, double resolution)
{
	double lenthVector12 = sqrt(pow(range.vector12.first, 2) + pow(range.vector12.second, 2));
	double lenthVector13 = sqrt(pow(range.vector13.first, 2) + pow(range.vector13.second, 2));
	int numsVector12 = ceil(lenthVector12 / resolution);
	int numsVector13 = ceil(lenthVector13 / resolution);

	for (int i = 0; i < numsVector13-1; i++) //防止超出影像范围，舍弃最边界处
	{
		for (int j = 0; j < numsVector12-1; j++)
		{
			double latOri = range.lat + range.vector12.first / numsVector12 * j + range.vector13.first / numsVector13 * i;
			double lonOri = range.lon + range.vector12.second / numsVector12 * j + range.vector13.second / numsVector13 * i;
			double xOriModel1, yOriModel1;	//	像素起点 模型一
			double xOriModel2, yOriModel2;	//	像素起点 模型二
			model1.fromlatlonh2xy(latOri, lonOri, model1.readHeight(latOri, lonOri), xOriModel1, yOriModel1); //起点有可能定位不准导致超出影像范围
			model2.fromlatlonh2xy(latOri, lonOri, model2.readHeight(latOri, lonOri), xOriModel2, yOriModel2);

			double latEnd = latOri + range.vector12.first / numsVector12 + range.vector13.first / numsVector13; 
			double lonEnd = lonOri + range.vector12.second / numsVector12 + range.vector13.second / numsVector13;
			double xEndModel1, yEndModel1;
			double xEndModel2, yEndModel2;
			model1.fromlatlonh2xy(latEnd, lonEnd, model1.readHeight(latEnd, lonEnd), xEndModel1, yEndModel1);
			model2.fromlatlonh2xy(latEnd, lonEnd, model2.readHeight(latEnd, lonEnd), xEndModel2, yEndModel2);

			double xModel1wid, yModel1Height;
			xModel1wid = xEndModel1 - xOriModel1;
			yModel1Height = yEndModel1 - yOriModel1;
			m_GridsModel1.push_back(DRect(max(0,xOriModel1), max(0,yOriModel1), xModel1wid, yModel1Height)); //起点大于等于0

			double xModel2wid, yModel2Height;
			xModel2wid = min(model2.m_Img.getImageWid(),xEndModel2) - xOriModel2;
			yModel2Height = min(model2.m_Img.getImageHei(),yEndModel2) - yOriModel2;
			m_GridsModel2.push_back(DRect(max(0,xOriModel2), max(0,yOriModel2), xModel2wid, yModel2Height));
		}
	}
	return true;
}

bool StereoSensorModel::segmentGrid3Img(GeoRange &range, double resolution)
{
	//首先分割S1 S2影像 m_GridsModel1 m_GridsModel2
	printf("\n格网化影像: ");
	segmentGrid(m_cModelS1, m_cModelS2, range, resolution);

	//分割ND影像 m_GridsModel3

	int numsVector12 = ceil(sqrt(pow(range.vector12.first, 2) + pow(range.vector12.second, 2)) / resolution);
	int numsVector13 = ceil(sqrt(pow(range.vector13.first, 2) + pow(range.vector13.second, 2)) / resolution);

	for (int i = 0; i < numsVector13-1; i++)
	{
		for (int j = 0; j < numsVector12-1; j++)
		{
			double latOri = range.lat + range.vector12.first / numsVector12 * j + range.vector13.first / numsVector13 * i;
			double lonOri = range.lon + range.vector12.second / numsVector12 * j + range.vector13.second / numsVector13 * i;
			double xOriModel, yOriModel;
			m_cModelNd.fromlatlonh2xy(latOri, lonOri, m_cModelNd.readHeight(latOri, lonOri), xOriModel, yOriModel);

			double latEnd = latOri + range.vector12.first / numsVector12 + range.vector13.first / numsVector13;
			double lonEnd = lonOri + range.vector12.second / numsVector12 + range.vector13.second / numsVector13;
			double xEndModel, yEndModel;
			m_cModelNd.fromlatlonh2xy(latEnd, lonEnd, m_cModelNd.readHeight(latEnd, lonEnd), xEndModel, yEndModel);

			double xModelwid, yModelHeight;
			xModelwid = xEndModel - xOriModel;
			yModelHeight = yEndModel - yOriModel;
			m_GridsModel3.push_back(DRect(max(0, xOriModel), max(0, yOriModel), xModelwid, yModelHeight)); //起点大于等于0
		}
	}

	return true;
}

bool StereoSensorModel::gridMatch(OpticalLinearSensorModel& model1, OpticalLinearSensorModel& model2, string sTiePointsFile)
{
	vector<tiepoint> vTiePoints;
	//逐个格网区间展开匹配
	for (int i = 0; i < m_GridsModel1.size(); i++)
	{
		void* pDataGridImg1 = nullptr;
		void* pDataGridImg2 = nullptr;
		model1.m_Img.readBlock(m_GridsModel1[i].x, m_GridsModel1[i].y, m_GridsModel1[i].width, m_GridsModel1[i].height, 0, pDataGridImg1, GDT_Float32); // 注意 元数据类型不能提取出特征
		model2.m_Img.readBlock(m_GridsModel2[i].x, m_GridsModel2[i].y, m_GridsModel2[i].width, m_GridsModel2[i].height, 0, pDataGridImg2, GDT_Float32);

		matchTest(m_GridsModel1[i].x, m_GridsModel1[i].y, m_GridsModel1[i].width, m_GridsModel1[i].height, pDataGridImg1, \
			m_GridsModel2[i].x, m_GridsModel2[i].y, m_GridsModel2[i].width, m_GridsModel2[i].height, pDataGridImg2, vTiePoints);

		delete[]pDataGridImg1, pDataGridImg1 = nullptr;
		delete[]pDataGridImg2, pDataGridImg2 = nullptr;
	}

	FILE* fp = fopen(sTiePointsFile.c_str(), "w");
	fprintf(fp, "\n");
	for (int i = 0; i < vTiePoints.size(); i++)
	{
		fprintf(fp, "%f\t%f\t%f\t%f\n", vTiePoints[i].x1,vTiePoints[i].y1,vTiePoints[i].x2,vTiePoints[i].y2);
	}
	fclose(fp);

	return true;
}

bool StereoSensorModel::gridMatch(string sTiePointsFile)
{
	vector<tie3point> vTiePoints;
	//逐个格网区间展开匹配
	for (int i = 0; i < m_GridsModel1.size(); i++)
	{
		cout << "\n正在执行第 " << i << "块匹配\n";
		void* pDataGridImg1 = nullptr;
		void* pDataGridImg2 = nullptr;
		void* pDataGridImg3 = nullptr;

		m_cModelS1.m_Img.readBlock(m_GridsModel1[i].x, m_GridsModel1[i].y, m_GridsModel1[i].width, m_GridsModel1[i].height, 0, pDataGridImg1, GDT_Float32); // 注意 元数据类型不能提取出特征
		m_cModelS2.m_Img.readBlock(m_GridsModel2[i].x, m_GridsModel2[i].y, m_GridsModel2[i].width, m_GridsModel2[i].height, 0, pDataGridImg2, GDT_Float32);
		m_cModelNd.m_Img.readBlock(m_GridsModel3[i].x, m_GridsModel3[i].y, m_GridsModel3[i].width, m_GridsModel3[i].height, 0, pDataGridImg3, GDT_Float32);

		threeImgsmatchTest(m_GridsModel1[i].x, m_GridsModel1[i].y, m_GridsModel1[i].width, m_GridsModel1[i].height, pDataGridImg1, \
			m_GridsModel2[i].x, m_GridsModel2[i].y, m_GridsModel2[i].width, m_GridsModel2[i].height, pDataGridImg2, \
			m_GridsModel3[i].x, m_GridsModel3[i].y, m_GridsModel3[i].width, m_GridsModel3[i].height, pDataGridImg3, vTiePoints);

		delete[]pDataGridImg1, pDataGridImg1 = nullptr;
		delete[]pDataGridImg2, pDataGridImg2 = nullptr;
		delete[]pDataGridImg3, pDataGridImg3 = nullptr;
	}

	FILE* fp = fopen(sTiePointsFile.c_str(), "w");
	fprintf(fp, "\n");
	for (int i = 0; i < vTiePoints.size(); i++)
	{
		fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\n", vTiePoints[i].x1, vTiePoints[i].y1, vTiePoints[i].x2, vTiePoints[i].y2, vTiePoints[i].x3, vTiePoints[i].y3);
	}
	fclose(fp);
	return true;
}

bool StereoSensorModel::readtiepoints(string sTiepointsFilepath)
{
	char c_read[1024];
	tie3point point;
	FILE* fp = fopen(sTiepointsFilepath.c_str(), "r");
	while (!feof(fp))
	{
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		if (sscanf(c_read, "%f%f%f%f%f%f", &point.x1, &point.y1, &point.x2, &point.y2,&point.x3,&point.y3) != 6)
			continue;
		m_tiePointsPix.push_back(point);
	}
	return true;
}

bool StereoSensorModel::fromtiepoint2XYZ(OpticalLinearSensorModel& model1, OpticalLinearSensorModel& model2, double x1, double y1, double x2, double y2, double& X, double& Y, double& Z)
{
	//相机1模型
	double et1 = model1.fromxyCalET(x1, y1);
	Eigen::Vector3d XYZs1;//卫星位置
	model1.frometCalpose(et1, XYZs1[0], XYZs1[1], XYZs1[2]);
	Eigen::Matrix3d R_att = model1.frometCalatt(et1);//卫星姿态
	Eigen::Matrix3d R1 = R_att*model1.get_RuMatrix() * model1.m_cameraInstall;
	Eigen::Vector3d XYZc1 = model1.fromxy2Camerapose(x1,y1);//相机坐标系
	Eigen::Vector3d XYZp1 = R1 * XYZc1;

	//相机2模型
	double et2 = model2.fromxyCalET(x2,y2);
	Eigen::Vector3d XYZs2;
	model2.frometCalpose(et2, XYZs2[0], XYZs2[1], XYZs2[2]);
	Eigen::Matrix3d R_att2 = model2.frometCalatt(et2);
	Eigen::Matrix3d R2 = R_att2 *model2.get_RuMatrix()* model2.m_cameraInstall;
	Eigen::Vector3d XYZc2 = model2.fromxy2Camerapose(x2, y2);
	Eigen::Vector3d XYZp2 = R2 * XYZc2;

	double Bx, By, Bz;
	Bx = XYZs2[0] - XYZs1[0];
	By = XYZs2[1] - XYZs1[1];
	Bz = XYZs2[2] - XYZs1[2];
	double Z2, X2, X1, Z1;
	X1 = XYZp1[0];
	X2 = XYZp2[0];
	Z1 = XYZp1[2];
	Z2 = XYZp2[2];
	double N1, N2;

	N1 = (Bx * Z2 - Bz * X2) / (X2 * Z1 - X1 * Z2);
	N2 = (Bx * Z1 - Bz * X1) / (X2 * Z1 - X1 * Z2);

	double Y1, Y2;
	Y1 = XYZp1[1];
	Y2 = XYZp2[1];

	X = XYZs1[0] - N1 * X1;
	Y = (XYZs1[1] - N1 * Y1 + XYZs2[1] - N2 * Y2) / 2.0;
	Z = XYZs1[2] - N1 * Z1;
	return true;
}

tiePoint3Geo StereoSensorModel::fromtiepoint2XYZ(double x1, double y1, double x2, double y2, double x3, double y3)
{
	tiePoint3Geo tie3XYZpoint;
	
	fromtiepoint2XYZ(m_cModelS1, m_cModelS2, x1, y1, x2, y2, tie3XYZpoint.X_s1s2, tie3XYZpoint.Y_s1s2, tie3XYZpoint.Z_s1s2);
	fromXYZ2latlonalt(tie3XYZpoint.X_s1s2, tie3XYZpoint.Y_s1s2, tie3XYZpoint.Z_s1s2, tie3XYZpoint.lat_s1s2, tie3XYZpoint.lon_s1s2, tie3XYZpoint.altitude_s1s2);
	
	fromtiepoint2XYZ(m_cModelS1, m_cModelNd, x1, y1, x3, y3, tie3XYZpoint.X_s1nd, tie3XYZpoint.Y_s1nd, tie3XYZpoint.Z_s1nd);
	fromXYZ2latlonalt(tie3XYZpoint.X_s1nd, tie3XYZpoint.Y_s1nd, tie3XYZpoint.Z_s1nd, tie3XYZpoint.lat_s1nd, tie3XYZpoint.lon_s1nd, tie3XYZpoint.altitude_s1nd);
	
	fromtiepoint2XYZ(m_cModelNd, m_cModelS2, x3, y3, x2, y2, tie3XYZpoint.X_nds2, tie3XYZpoint.Y_nds2, tie3XYZpoint.Z_nds2);
	fromXYZ2latlonalt(tie3XYZpoint.X_nds2, tie3XYZpoint.Y_nds2, tie3XYZpoint.Z_nds2, tie3XYZpoint.lat_nds2, tie3XYZpoint.lon_nds2, tie3XYZpoint.altitude_nds2);

	//将S1和S2交会的值，赋值给最佳值
	tie3XYZpoint.X = tie3XYZpoint.X_s1s2;
	tie3XYZpoint.Y = tie3XYZpoint.Y_s1s2;
	tie3XYZpoint.Z = tie3XYZpoint.Z_s1s2;

	//迭代计算最佳XYZ
	// 旋转矩阵9个参数 *3
	//相机1模型
	double et1 = m_cModelS1.fromxyCalET(x1, y1);
	Eigen::Matrix3d R1 = m_cModelS1.frometCalatt(et1) * m_cModelS1.get_RuMatrix() * m_cModelS1.m_cameraInstall;
	Eigen::Vector3d XYZs1;//卫星位置
	m_cModelS1.frometCalpose(et1, XYZs1[0], XYZs1[1], XYZs1[2]);
	Eigen::Vector3d XYZc1 = m_cModelS1.fromxy2Camerapose(x1, y1);//tanx tany 观测值
	//相机2模型
	double et2 = m_cModelS2.fromxyCalET(x2, y2);
	Eigen::Matrix3d R2 = m_cModelS2.frometCalatt(et2) * m_cModelS2.get_RuMatrix() * m_cModelS2.m_cameraInstall;
	Eigen::Vector3d XYZs2;//卫星位置
	m_cModelS2.frometCalpose(et2, XYZs2[0], XYZs2[1], XYZs2[2]);
	Eigen::Vector3d XYZc2 = m_cModelS2.fromxy2Camerapose(x2, y2);
	//相机3模型
	double et3 = m_cModelNd.fromxyCalET(x3, y3);
	Eigen::Matrix3d R3 = m_cModelNd.frometCalatt(et3) * m_cModelNd.get_RuMatrix() * m_cModelNd.m_cameraInstall;
	Eigen::Vector3d XYZs3;//卫星位置
	m_cModelNd.frometCalpose(et3, XYZs3[0], XYZs3[1], XYZs3[2]);
	Eigen::Vector3d XYZc3 = m_cModelNd.fromxy2Camerapose(x3, y3);

	caltiepointBestXYZ(R1, XYZs1, XYZc1, R2, XYZs2, XYZc2, R3, XYZs3, XYZc3, tie3XYZpoint.X, tie3XYZpoint.Y, tie3XYZpoint.Z);
	fromXYZ2latlonalt(tie3XYZpoint.X, tie3XYZpoint.Y, tie3XYZpoint.Z, tie3XYZpoint.lat, tie3XYZpoint.lon, tie3XYZpoint.altitude);

	return tie3XYZpoint;
}

bool StereoSensorModel::caltiepointBestXYZ(Eigen::Matrix3d R1, Eigen::Vector3d XYZs1, Eigen::Vector3d XYZc1, Eigen::Matrix3d R2, Eigen::Vector3d XYZs2, Eigen::Vector3d XYZc2, Eigen::Matrix3d R3, Eigen::Vector3d XYZs3, Eigen::Vector3d XYZc3, double &X_ori, double &Y_ori, double &Z_ori)
{
	int i = 0;//迭代次数
	do
	{
		// 计算X' Y' Z'3个参数 *3
		Eigen::Vector3d XYZ_pie_s1;
		XYZ_pie_s1[0] = R1(0, 0) * (X_ori - XYZs1[0]) + R1(1, 0) * (Y_ori - XYZs1[1]) + R1(2, 0) * (Z_ori - XYZs1[2]);
		XYZ_pie_s1[1] = R1(0, 1) * (X_ori - XYZs1[0]) + R1(1, 1) * (Y_ori - XYZs1[1]) + R1(2, 1) * (Z_ori - XYZs1[2]);
		XYZ_pie_s1[2] = R1(0, 2) * (X_ori - XYZs1[0]) + R1(1, 2) * (Y_ori - XYZs1[1]) + R1(2, 2) * (Z_ori - XYZs1[2]);
		Eigen::Vector3d XYZ_pie_s2;
		XYZ_pie_s2[0] = R2(0, 0) * (X_ori - XYZs2[0]) + R2(1, 0) * (Y_ori - XYZs2[1]) + R2(2, 0) * (Z_ori - XYZs2[2]);
		XYZ_pie_s2[1] = R2(0, 1) * (X_ori - XYZs2[0]) + R2(1, 1) * (Y_ori - XYZs2[1]) + R2(2, 1) * (Z_ori - XYZs2[2]);
		XYZ_pie_s2[2] = R2(0, 2) * (X_ori - XYZs2[0]) + R2(1, 2) * (Y_ori - XYZs2[1]) + R2(2, 2) * (Z_ori - XYZs2[2]);
		Eigen::Vector3d XYZ_pie_nd;
		XYZ_pie_nd[0] = R3(0, 0) * (X_ori - XYZs3[0]) + R3(1, 0) * (Y_ori - XYZs3[1]) + R3(2, 0) * (Z_ori - XYZs3[2]);
		XYZ_pie_nd[1] = R3(0, 1) * (X_ori - XYZs3[0]) + R3(1, 1) * (Y_ori - XYZs3[1]) + R3(2, 1) * (Z_ori - XYZs3[2]);
		XYZ_pie_nd[2] = R3(0, 2) * (X_ori - XYZs3[0]) + R3(1, 2) * (Y_ori - XYZs3[1]) + R3(2, 2) * (Z_ori - XYZs3[2]);

		//计算矩阵A（6*3）和L（6*1）
		Eigen::Matrix<double, 6, 3> A;
		A(0, 0) = (R1(0, 0) * XYZ_pie_s1[2] - R1(0, 2) * XYZ_pie_s1[0]) / XYZ_pie_s1[2] / XYZ_pie_s1[2];
		A(0, 1) = (R1(1, 0) * XYZ_pie_s1[2] - R1(1, 2) * XYZ_pie_s1[0]) / XYZ_pie_s1[2] / XYZ_pie_s1[2];
		A(0, 2) = (R1(2, 0) * XYZ_pie_s1[2] - R1(2, 2) * XYZ_pie_s1[0]) / XYZ_pie_s1[2] / XYZ_pie_s1[2];
		A(1, 0) = (R1(0, 1) * XYZ_pie_s1[2] - R1(0, 2) * XYZ_pie_s1[1]) / XYZ_pie_s1[2] / XYZ_pie_s1[2];
		A(1, 1) = (R1(1, 1) * XYZ_pie_s1[2] - R1(1, 2) * XYZ_pie_s1[1]) / XYZ_pie_s1[2] / XYZ_pie_s1[2];
		A(1, 2) = (R1(2, 1) * XYZ_pie_s1[2] - R1(2, 2) * XYZ_pie_s1[1]) / XYZ_pie_s1[2] / XYZ_pie_s1[2];

		A(2, 0) = (R2(0, 0) * XYZ_pie_s2[2] - R2(0, 2) * XYZ_pie_s2[0]) / XYZ_pie_s2[2] / XYZ_pie_s2[2];
		A(2, 1) = (R2(1, 0) * XYZ_pie_s2[2] - R2(1, 2) * XYZ_pie_s2[0]) / XYZ_pie_s2[2] / XYZ_pie_s2[2];
		A(2, 2) = (R2(2, 0) * XYZ_pie_s2[2] - R2(2, 2) * XYZ_pie_s2[0]) / XYZ_pie_s2[2] / XYZ_pie_s2[2];
		A(3, 0) = (R2(0, 1) * XYZ_pie_s2[2] - R2(0, 2) * XYZ_pie_s2[1]) / XYZ_pie_s2[2] / XYZ_pie_s2[2];
		A(3, 1) = (R2(1, 1) * XYZ_pie_s2[2] - R2(1, 2) * XYZ_pie_s2[1]) / XYZ_pie_s2[2] / XYZ_pie_s2[2];
		A(3, 2) = (R2(2, 1) * XYZ_pie_s2[2] - R2(2, 2) * XYZ_pie_s2[1]) / XYZ_pie_s2[2] / XYZ_pie_s2[2];

		A(4, 0) = (R3(0, 0) * XYZ_pie_nd[2] - R3(0, 2) * XYZ_pie_nd[0]) / XYZ_pie_nd[2] / XYZ_pie_nd[2];
		A(4, 1) = (R3(1, 0) * XYZ_pie_nd[2] - R3(1, 2) * XYZ_pie_nd[0]) / XYZ_pie_nd[2] / XYZ_pie_nd[2];
		A(4, 2) = (R3(2, 0) * XYZ_pie_nd[2] - R3(2, 2) * XYZ_pie_nd[0]) / XYZ_pie_nd[2] / XYZ_pie_nd[2];
		A(5, 0) = (R3(0, 1) * XYZ_pie_nd[2] - R3(0, 2) * XYZ_pie_nd[1]) / XYZ_pie_nd[2] / XYZ_pie_nd[2];
		A(5, 1) = (R3(1, 1) * XYZ_pie_nd[2] - R3(1, 2) * XYZ_pie_nd[1]) / XYZ_pie_nd[2] / XYZ_pie_nd[2];
		A(5, 2) = (R3(2, 1) * XYZ_pie_nd[2] - R3(2, 2) * XYZ_pie_nd[1]) / XYZ_pie_nd[2] / XYZ_pie_nd[2];

		//Eigen::Matrix<double, 6, 6> P;//单位阵

		Eigen::Matrix<double, 6, 1> L;
		L(0, 0) = XYZc1[0] - XYZ_pie_s1[0] / XYZ_pie_s1[2];
		L(1, 0) = XYZc1[1] - XYZ_pie_s1[1] / XYZ_pie_s1[2];

		L(2, 0) = XYZc2[0] - XYZ_pie_s2[0] / XYZ_pie_s2[2];
		L(3, 0) = XYZc2[1] - XYZ_pie_s2[1] / XYZ_pie_s2[2];

		L(4, 0) = XYZc3[0] - XYZ_pie_nd[0] / XYZ_pie_nd[2];
		L(5, 0) = XYZc3[1] - XYZ_pie_nd[1] / XYZ_pie_nd[2];

		Eigen::Matrix3d N = A.transpose() * A;
		Eigen::Matrix3d N_Inver = N.inverse();
		Eigen::Vector3d W = A.transpose() * L;
		Eigen::Vector3d deita = N_Inver * W;

		//更新地理坐标
		X_ori += deita(0, 0);
		Y_ori += deita(1, 0);
		Z_ori += deita(2, 0);

		double increment = abs(deita(0, 0)) + abs(deita(1, 0)) + abs(deita(2, 0));
		if (increment <= 1)
			break;

		i++;
		if (i > 20)
			return false;

	} while (1);


	return true;
}

bool StereoSensorModel::savetiePointsPixel(string sPixelfilePath)
{
	string s1s2tiePoints = sPixelfilePath + "s1s2.pts";
	string s1ndtiePoints = sPixelfilePath + "s1nd.pts";
	string s2ndtiePoints = sPixelfilePath + "s2nd.pts";

	FILE* fp1 = fopen(s1s2tiePoints.c_str(), "w");
	FILE* fp2 = fopen(s1ndtiePoints.c_str(), "w");
	FILE* fp3 = fopen(s2ndtiePoints.c_str(), "w");

	fprintf(fp1, "\n");
	fprintf(fp2, "\n");
	fprintf(fp3, "\n");

	for (int i = 0; i < m_tiePointsPix.size(); i++)
	{
		fprintf(fp1, "%lf\t %lf\t %lf\t %lf\n", m_tiePointsPix[i].x1, m_tiePointsPix[i].y1, m_tiePointsPix[i].x2, m_tiePointsPix[i].y2);
		fprintf(fp2, "%lf\t %lf\t %lf\t %lf\n", m_tiePointsPix[i].x1, m_tiePointsPix[i].y1, m_tiePointsPix[i].x3, m_tiePointsPix[i].y3);
		fprintf(fp3, "%lf\t %lf\t %lf\t %lf\n", m_tiePointsPix[i].x2, m_tiePointsPix[i].y2, m_tiePointsPix[i].x3, m_tiePointsPix[i].y3);
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	
	return false;
}

bool StereoSensorModel::caltieGeopose(bool eliminateTiepoints)
{
	for (int i = 0; i < m_tiePointsPix.size(); i++)
	{
		printf("%d \n", i);
		tiePoint3Geo point = fromtiepoint2XYZ(m_tiePointsPix[i].x1, m_tiePointsPix[i].y1, m_tiePointsPix[i].x2, m_tiePointsPix[i].y2, m_tiePointsPix[i].x3, m_tiePointsPix[i].y3);
		m_tiePointsGeo.push_back(point);
		tieResidual pixelpoint;
		caltieGeoresidual(m_cModelS1, m_tiePointsPix[i].x1, m_tiePointsPix[i].y1, point.lat, point.lon, point.altitude, pixelpoint.x1_pixel_res, pixelpoint.y1_pixel_res);
		caltieGeoresidual(m_cModelS2, m_tiePointsPix[i].x2, m_tiePointsPix[i].y2, point.lat, point.lon, point.altitude, pixelpoint.x2_pixel_res, pixelpoint.y2_pixel_res);
		caltieGeoresidual(m_cModelNd, m_tiePointsPix[i].x3, m_tiePointsPix[i].y3, point.lat, point.lon, point.altitude, pixelpoint.x3_pixel_res, pixelpoint.y3_pixel_res);
		m_tiePointsRes.push_back(pixelpoint);
	}
	if (eliminateTiepoints)
	{
		eliminat_tiePoints_residual();
	}
	return true;
}

bool StereoSensorModel::updataResidual()
{
	m_tiePointsRes.clear();
	for (int i = 0; i < m_tiePointsPix.size(); i++)
	{
		fromXYZ2latlonalt(m_tiePointsGeo[i].X, m_tiePointsGeo[i].Y, m_tiePointsGeo[i].Z, m_tiePointsGeo[i].lat, m_tiePointsGeo[i].lon, m_tiePointsGeo[i].altitude);
		tieResidual pixelpoint;
		caltieGeoresidual(m_cModelS1, m_tiePointsPix[i].x1, m_tiePointsPix[i].y1, m_tiePointsGeo[i].lat, m_tiePointsGeo[i].lon, m_tiePointsGeo[i].altitude, pixelpoint.x1_pixel_res, pixelpoint.y1_pixel_res);
		caltieGeoresidual(m_cModelS2, m_tiePointsPix[i].x2, m_tiePointsPix[i].y2, m_tiePointsGeo[i].lat, m_tiePointsGeo[i].lon, m_tiePointsGeo[i].altitude, pixelpoint.x2_pixel_res, pixelpoint.y2_pixel_res);
		caltieGeoresidual(m_cModelNd, m_tiePointsPix[i].x3, m_tiePointsPix[i].y3, m_tiePointsGeo[i].lat, m_tiePointsGeo[i].lon, m_tiePointsGeo[i].altitude, pixelpoint.x3_pixel_res, pixelpoint.y3_pixel_res);
		m_tiePointsRes.push_back(pixelpoint);
	}
	return true;
}

bool StereoSensorModel::savetiePointsGeopose(string sfilePath)
{
	//交会点地理坐标
	FILE* fp = fopen((sfilePath + "s1s2nd_Geo.txt").c_str(), "w");
	fprintf(fp, "#S1S2(X Y Z LAT LON ALT) S1ND(X Y Z LAT LON ALT) NDS2(X Y Z LAT LON ALT) BEST(X Y Z LAT LON ALT)\n");
	for (int i = 0; i < m_tiePointsGeo.size(); i++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", \
			m_tiePointsGeo[i].X_s1s2, m_tiePointsGeo[i].Y_s1s2, m_tiePointsGeo[i].Z_s1s2, m_tiePointsGeo[i].lat_s1s2, m_tiePointsGeo[i].lon_s1s2, m_tiePointsGeo[i].altitude_s1s2, \
			m_tiePointsGeo[i].X_s1nd, m_tiePointsGeo[i].Y_s1nd, m_tiePointsGeo[i].Z_s1nd, m_tiePointsGeo[i].lat_s1nd, m_tiePointsGeo[i].lon_s1nd, m_tiePointsGeo[i].altitude_s1nd, \
			m_tiePointsGeo[i].X_nds2, m_tiePointsGeo[i].Y_nds2, m_tiePointsGeo[i].Z_nds2, m_tiePointsGeo[i].lat_nds2, m_tiePointsGeo[i].lon_nds2, m_tiePointsGeo[i].altitude_nds2, \
			m_tiePointsGeo[i].X, m_tiePointsGeo[i].Y, m_tiePointsGeo[i].Z, m_tiePointsGeo[i].lat, m_tiePointsGeo[i].lon, m_tiePointsGeo[i].altitude);
	}
	fclose(fp);

	//交会点高程误差
	FILE* fp_h = fopen((sfilePath + "GeoHight_accuracy.txt").c_str(), "w");
	for (int j = 0; j < m_tiePointsGeo.size(); j++)
	{
		//已知经度纬度，读取dem高程值
		double diffH = m_tiePointsGeo[j].altitude - getHvaluefromDEM(m_tiePointsGeo[j].lat, m_tiePointsGeo[j].lon);
		fprintf(fp_h, "%lf\t %lf\t %lf\n", m_tiePointsGeo[j].lat, m_tiePointsGeo[j].lon, diffH);
	}
	fclose(fp_h);

	return true;
}

bool StereoSensorModel::caltieGeoresidual(OpticalLinearSensorModel& model, double x, double y, double lat, double lon, double alt, double& x_res, double& y_res)
{
	double x_pixel_new = 0.0, y_pixel_new = 0.0;

	model.fromlatlonh2xy(lat, lon, alt, x_pixel_new, y_pixel_new);

	x_res = x - x_pixel_new;
	y_res = y - y_pixel_new;

	return true;
}

bool StereoSensorModel::savetiePointsResival(string sResivalfilePath)
{
	FILE* fp_S1 = fopen((sResivalfilePath + "S1_Residual.txt").c_str(), "w");
	FILE* fp_S2 = fopen((sResivalfilePath + "S2_Residual.txt").c_str(), "w");
	FILE* fp_ND = fopen((sResivalfilePath + "ND_Residual.txt").c_str(), "w");

	fprintf(fp_S1, "#S1(x y x_residual y_residual)\n");
	fprintf(fp_S2, "#S2(x y x_residual y_residual)\n");
	fprintf(fp_ND, "#ND(x y x_residual y_residual)\n");
	for (int i = 0; i < m_tiePointsRes.size(); i++)
	{
		fprintf(fp_S1, "%lf\t%lf\t%lf\t%lf\n", m_tiePointsPix[i].x1, m_tiePointsPix[i].y1, m_tiePointsRes[i].x1_pixel_res, m_tiePointsRes[i].y1_pixel_res);
		fprintf(fp_S2, "%lf\t%lf\t%lf\t%lf\n", m_tiePointsPix[i].x2, m_tiePointsPix[i].y2, m_tiePointsRes[i].x2_pixel_res, m_tiePointsRes[i].y2_pixel_res);
		fprintf(fp_ND, "%lf\t%lf\t%lf\t%lf\n", m_tiePointsPix[i].x3, m_tiePointsPix[i].y3, m_tiePointsRes[i].x3_pixel_res, m_tiePointsRes[i].y3_pixel_res);
	}
	fclose(fp_S1);
	fclose(fp_S2);
	fclose(fp_ND);

	return true;
}

bool StereoSensorModel::eliminat_tiePoints_residual(double x_ResidualThreshold, double y_ResidualThreshold)
{
	//遍历所有连接点
	for (int i = 0; i < m_tiePointsPix.size(); i++)
	{
		//cout <<"保留的连接点数量 " << m_tiePointsPix.size() << endl;
		//判断S1残差是否超出阈值
		if (abs(m_tiePointsRes[i].x1_pixel_res)>=x_ResidualThreshold|| abs(m_tiePointsRes[i].y1_pixel_res) >= y_ResidualThreshold)
		{
			m_tiePointsPix.erase(m_tiePointsPix.begin() + i);
			m_tiePointsGeo.erase(m_tiePointsGeo.begin() + i);
			m_tiePointsRes.erase(m_tiePointsRes.begin() + i);
			i--;
		}
	}

	return true;
}

bool StereoSensorModel::extCalibrate()
{
	int tiepoints_num = m_tiePointsPix.size();

	//迭代计算外定标参数和地物点坐标
	do
	{
		double vError = 0.0;//统计目标函数
		//V=Ax+Bt-L
		Eigen::SparseMatrix<double> spt_A(tiepoints_num * 6, 9);//维度 6*tiepoints_num,9
		vector<T> coefficients_A;

		Eigen::SparseMatrix<double> spt_B(tiepoints_num * 6, tiepoints_num * 3);//纬度 6*tiepoints_num,3*tiepoints_num
		vector<T> coefficients_B;

		Eigen::VectorXd L(6 * tiepoints_num);//维度 6*tiepoints_num,1

		//外参初值 单位/弧度 个数/9
		double phi_s1 = m_cModelS1.m_camera_OffsetAngleMatrix.m_extCalPar_phi * rpd_c();
		double omega_s1 = m_cModelS1.m_camera_OffsetAngleMatrix.m_extCalPar_omega * rpd_c();
		double kappa_s1 = m_cModelS1.m_camera_OffsetAngleMatrix.m_extCalPar_kappa * rpd_c();
		double phi_s2 = m_cModelS2.m_camera_OffsetAngleMatrix.m_extCalPar_phi * rpd_c();
		double omega_s2 = m_cModelS2.m_camera_OffsetAngleMatrix.m_extCalPar_omega * rpd_c();
		double kappa_s2 = m_cModelS2.m_camera_OffsetAngleMatrix.m_extCalPar_kappa * rpd_c();
		double phi_nd = m_cModelNd.m_camera_OffsetAngleMatrix.m_extCalPar_phi * rpd_c();
		double omega_nd = m_cModelNd.m_camera_OffsetAngleMatrix.m_extCalPar_omega * rpd_c();
		double kappa_nd = m_cModelNd.m_camera_OffsetAngleMatrix.m_extCalPar_kappa * rpd_c();
		printf("phi omega kappa:%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", phi_s1, omega_s1, kappa_s1, phi_s2, omega_s2, kappa_s2, phi_nd, omega_nd, kappa_nd);

		//矩阵按行赋值,6行即一组连接点循环一次
		for (int i = 0; i < tiepoints_num; i++)
		{
			//地面点坐标初值 单位/米
			double X_GP = m_tiePointsGeo[i].X;
			double Y_GP = m_tiePointsGeo[i].Y;
			double Z_GP = m_tiePointsGeo[i].Z;

			//3景影像 S1 S2 ND
			//S1的外定标参数系数矩阵A
			double phi = phi_s1, omega = omega_s1, kappa = kappa_s1;

			double Da1Dphi = -sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
			double Da1Domega = sin(phi) * cos(omega) * sin(kappa);
			double Da1Dkappa = -cos(phi) * sin(kappa) + sin(phi) * sin(omega) * cos(kappa);
			double Db1Dphi = 0;
			double Db1Domega = -sin(omega) * sin(kappa);
			double Db1Dkappa = cos(omega) * cos(kappa);
			double Dc1Dphi = -cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
			double Dc1Domega = cos(phi) * cos(omega) * sin(kappa);
			double Dc1Dkappa = sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);

			double Da2Dphi = sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
			double Da2Domega = sin(phi) * cos(omega) * cos(kappa);
			double Da2Dkappa = -cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
			double Db2Dphi = 0;
			double Db2Domega = -sin(omega) * cos(kappa);
			double Db2Dkappa = -cos(omega) * sin(kappa);
			double Dc2Dphi = cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
			double Dc2Domega = cos(phi) * cos(omega) * cos(kappa);
			double Dc2Dkappa = sin(phi) * cos(kappa) - cos(phi) * sin(omega) * sin(kappa);

			double Da3Dphi = cos(phi) * cos(omega);
			double Da3Domega = -sin(phi) * sin(omega);
			double Da3Dkappa = 0;
			double Db3Dphi = 0;
			double Db3Domega = -cos(omega);
			double Db3Dkappa = 0;
			double Dc3Dphi = -sin(phi) * cos(omega);
			double Dc3Domega = -cos(phi) * sin(omega);
			double Dc3Dkappa = 0;

			double et = m_cModelS1.fromxyCalET(m_tiePointsPix[i].x1, m_tiePointsPix[i].y1);
			double ste_pose_X, ste_pose_Y, ste_pose_Z;
			m_cModelS1.frometCalpose(et, ste_pose_X, ste_pose_Y, ste_pose_Z);
			Eigen::Vector3d vec;//卫星指向地面的向量
			vec << X_GP - ste_pose_X, Y_GP - ste_pose_Y, Z_GP - ste_pose_Z;
			Eigen::Matrix3d mat = m_cModelS1.frometCalatt(et).inverse();
			Eigen::Vector3d XYZ_ = mat * vec;
			double X_ = XYZ_[0], Y_ = XYZ_[1], Z_ = XYZ_[2];

			Eigen::Vector3d XX_YY_ZZ = m_cModelS1.get_RuMatrix().inverse() * XYZ_;
			double XX = XX_YY_ZZ[0], YY = XX_YY_ZZ[1], ZZ = XX_YY_ZZ[2];

			Eigen::Vector3d XYZbody = m_cModelS1.get_InstallMatrix() * m_cModelS1.fromxy2Camerapose(m_tiePointsPix[i].x1, m_tiePointsPix[i].y1);
			double V_1 = XX / ZZ - XYZbody[0] / XYZbody[2];
			double V_2 = YY / ZZ - XYZbody[1] / XYZbody[2];
			vError += abs(V_1);
			vError += abs(V_2);

			L(6 * i + 0) = -V_1;
			L(6 * i + 1) = -V_2;

			double dXXdphi = Da1Dphi * X_ + Db1Dphi * Y_ + Dc1Dphi * Z_;
			double dXXdomega = Da1Domega * X_ + Db1Domega * Y_ + Dc1Domega * Z_;
			double dXXdkappa = Da1Dkappa * X_ + Db1Dkappa * Y_ + Dc1Dkappa * Z_;
			double dYYdphi = Da2Dphi * X_ + Db2Dphi * Y_ + Dc2Dphi * Z_;
			double dYYdomega = Da2Domega * X_ + Db2Domega * Y_ + Dc2Domega * Z_;
			double dYYdkappa = Da2Dkappa * X_ + Db2Dkappa * Y_ + Dc2Dkappa * Z_;
			double dZZdphi = Da3Dphi * X_ + Db3Dphi * Y_ + Dc3Dphi * Z_;
			double dZZdomega = Da3Domega * X_ + Db3Domega * Y_ + Dc3Domega * Z_;
			double dZZdkappa = Da3Dkappa * X_ + Db3Dkappa * Y_ + Dc3Dkappa * Z_;

			Eigen::Matrix3d XXYYZZ_XYZ = m_cModelS1.frometCalatt(et) * m_cModelS1.get_RuMatrix();
			double dXXdX = XXYYZZ_XYZ(0, 0), dYYdX = XXYYZZ_XYZ(0, 1), dZZdX = XXYYZZ_XYZ(0, 2);
			double dXXdY = XXYYZZ_XYZ(1, 0), dYYdY = XXYYZZ_XYZ(1, 1), dZZdY = XXYYZZ_XYZ(1, 2);
			double dXXdZ = XXYYZZ_XYZ(2, 0), dYYdZ = XXYYZZ_XYZ(2, 1), dZZdZ = XXYYZZ_XYZ(2, 2);

			coefficients_A.push_back(T(6 * i + 0, 0, (dXXdphi * ZZ - dZZdphi * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 0, 1, (dXXdomega * ZZ - dZZdomega * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 0, 2, (dXXdkappa * ZZ - dZZdkappa * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 1, 0, (dYYdphi * ZZ - dZZdphi * YY) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 1, 1, (dYYdomega * ZZ - dZZdomega * YY) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 1, 2, (dYYdkappa * ZZ - dZZdkappa * YY) / ZZ / ZZ));

			coefficients_B.push_back(T(6 * i + 0, 3 * i + 0, (dXXdX * ZZ - dZZdX * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 0, 3 * i + 1, (dXXdY * ZZ - dZZdY * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 0, 3 * i + 2, (dXXdZ * ZZ - dZZdZ * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 1, 3 * i + 0, (dYYdX * ZZ - dZZdX * YY) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 1, 3 * i + 1, (dYYdY * ZZ - dZZdY * YY) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 1, 3 * i + 2, (dYYdZ * ZZ - dZZdZ * YY) / ZZ / ZZ));

			//S2
			phi = phi_s2, omega = omega_s2, kappa = kappa_s2;
			Da1Dphi = -sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
			Da1Domega = sin(phi) * cos(omega) * sin(kappa);
			Da1Dkappa = -cos(phi) * sin(kappa) + sin(phi) * sin(omega) * cos(kappa);
			Db1Dphi = 0;
			Db1Domega = -sin(omega) * sin(kappa);
			Db1Dkappa = cos(omega) * cos(kappa);
			Dc1Dphi = -cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
			Dc1Domega = cos(phi) * cos(omega) * sin(kappa);
			Dc1Dkappa = sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);

			Da2Dphi = sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
			Da2Domega = sin(phi) * cos(omega) * cos(kappa);
			Da2Dkappa = -cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
			Db2Dphi = 0;
			Db2Domega = -sin(omega) * cos(kappa);
			Db2Dkappa = -cos(omega) * sin(kappa);
			Dc2Dphi = cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
			Dc2Domega = cos(phi) * cos(omega) * cos(kappa);
			Dc2Dkappa = sin(phi) * cos(kappa) - cos(phi) * sin(omega) * sin(kappa);

			Da3Dphi = cos(phi) * cos(omega);
			Da3Domega = -sin(phi) * sin(omega);
			Da3Dkappa = 0;
			Db3Dphi = 0;
			Db3Domega = -cos(omega);
			Db3Dkappa = 0;
			Dc3Dphi = -sin(phi) * cos(omega);
			Dc3Domega = -cos(phi) * sin(omega);
			Dc3Dkappa = 0;

			et = m_cModelS2.fromxyCalET(m_tiePointsPix[i].x2, m_tiePointsPix[i].y2);
			m_cModelS2.frometCalpose(et, ste_pose_X, ste_pose_Y, ste_pose_Z);
			vec << X_GP - ste_pose_X, Y_GP - ste_pose_Y, Z_GP - ste_pose_Z;
			mat = m_cModelS2.frometCalatt(et).inverse();
			XYZ_ = mat * vec;
			X_ = XYZ_[0], Y_ = XYZ_[1], Z_ = XYZ_[2];

			XX_YY_ZZ = m_cModelS2.get_RuMatrix().inverse() * XYZ_;
			XX = XX_YY_ZZ[0], YY = XX_YY_ZZ[1], ZZ = XX_YY_ZZ[2];

			XYZbody = m_cModelS2.get_InstallMatrix() * m_cModelS2.fromxy2Camerapose(m_tiePointsPix[i].x2, m_tiePointsPix[i].y2);
			V_1 = XX / ZZ - XYZbody[0] / XYZbody[2];
			V_2 = YY / ZZ - XYZbody[1] / XYZbody[2];

			vError += abs(V_1);
			vError += abs(V_2);

			L(6 * i + 2) = -V_1;
			L(6 * i + 3) = -V_2;

			dXXdphi = Da1Dphi * X_ + Db1Dphi * Y_ + Dc1Dphi * Z_;
			dXXdomega = Da1Domega * X_ + Db1Domega * Y_ + Dc1Domega * Z_;
			dXXdkappa = Da1Dkappa * X_ + Db1Dkappa * Y_ + Dc1Dkappa * Z_;
			dYYdphi = Da2Dphi * X_ + Db2Dphi * Y_ + Dc2Dphi * Z_;
			dYYdomega = Da2Domega * X_ + Db2Domega * Y_ + Dc2Domega * Z_;
			dYYdkappa = Da2Dkappa * X_ + Db2Dkappa * Y_ + Dc2Dkappa * Z_;
			dZZdphi = Da3Dphi * X_ + Db3Dphi * Y_ + Dc3Dphi * Z_;
			dZZdomega = Da3Domega * X_ + Db3Domega * Y_ + Dc3Domega * Z_;
			dZZdkappa = Da3Dkappa * X_ + Db3Dkappa * Y_ + Dc3Dkappa * Z_;

			XXYYZZ_XYZ = m_cModelS2.frometCalatt(et) * m_cModelS2.get_RuMatrix();
			dXXdX = XXYYZZ_XYZ(0, 0), dYYdX = XXYYZZ_XYZ(0, 1), dZZdX = XXYYZZ_XYZ(0, 2);
			dXXdY = XXYYZZ_XYZ(1, 0), dYYdY = XXYYZZ_XYZ(1, 1), dZZdY = XXYYZZ_XYZ(1, 2);
			dXXdZ = XXYYZZ_XYZ(2, 0), dYYdZ = XXYYZZ_XYZ(2, 1), dZZdZ = XXYYZZ_XYZ(2, 2);

			coefficients_A.push_back(T(6 * i + 2, 3, (dXXdphi * ZZ - dZZdphi * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 2, 4, (dXXdomega * ZZ - dZZdomega * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 2, 5, (dXXdkappa * ZZ - dZZdkappa * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 3, 3, (dYYdphi * ZZ - dZZdphi * YY) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 3, 4, (dYYdomega * ZZ - dZZdomega * YY) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 3, 5, (dYYdkappa * ZZ - dZZdkappa * YY) / ZZ / ZZ));

			coefficients_B.push_back(T(6 * i + 2, 3 * i + 0, (dXXdX * ZZ - dZZdX * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 2, 3 * i + 1, (dXXdY * ZZ - dZZdY * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 2, 3 * i + 2, (dXXdZ * ZZ - dZZdZ * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 3, 3 * i + 0, (dYYdX * ZZ - dZZdX * YY) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 3, 3 * i + 1, (dYYdY * ZZ - dZZdY * YY) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 3, 3 * i + 2, (dYYdZ * ZZ - dZZdZ * YY) / ZZ / ZZ));

			//ND
			phi = phi_nd, omega = omega_nd, kappa = kappa_nd;
			Da1Dphi = -sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
			Da1Domega = sin(phi) * cos(omega) * sin(kappa);
			Da1Dkappa = -cos(phi) * sin(kappa) + sin(phi) * sin(omega) * cos(kappa);
			Db1Dphi = 0;
			Db1Domega = -sin(omega) * sin(kappa);
			Db1Dkappa = cos(omega) * cos(kappa);
			Dc1Dphi = -cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
			Dc1Domega = cos(phi) * cos(omega) * sin(kappa);
			Dc1Dkappa = sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);

			Da2Dphi = sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
			Da2Domega = sin(phi) * cos(omega) * cos(kappa);
			Da2Dkappa = -cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
			Db2Dphi = 0;
			Db2Domega = -sin(omega) * cos(kappa);
			Db2Dkappa = -cos(omega) * sin(kappa);
			Dc2Dphi = cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
			Dc2Domega = cos(phi) * cos(omega) * cos(kappa);
			Dc2Dkappa = sin(phi) * cos(kappa) - cos(phi) * sin(omega) * sin(kappa);

			Da3Dphi = cos(phi) * cos(omega);
			Da3Domega = -sin(phi) * sin(omega);
			Da3Dkappa = 0;
			Db3Dphi = 0;
			Db3Domega = -cos(omega);
			Db3Dkappa = 0;
			Dc3Dphi = -sin(phi) * cos(omega);
			Dc3Domega = -cos(phi) * sin(omega);
			Dc3Dkappa = 0;

			et = m_cModelNd.fromxyCalET(m_tiePointsPix[i].x3, m_tiePointsPix[i].y3);
			m_cModelNd.frometCalpose(et, ste_pose_X, ste_pose_Y, ste_pose_Z);
			vec << X_GP - ste_pose_X, Y_GP - ste_pose_Y, Z_GP - ste_pose_Z;
			mat = m_cModelNd.frometCalatt(et).inverse();
			XYZ_ = mat * vec;
			X_ = XYZ_[0], Y_ = XYZ_[1], Z_ = XYZ_[2];

			XX_YY_ZZ = m_cModelNd.get_RuMatrix().inverse() * XYZ_;
			XX = XX_YY_ZZ[0], YY = XX_YY_ZZ[1], ZZ = XX_YY_ZZ[2];

			XYZbody = m_cModelNd.get_InstallMatrix() * m_cModelNd.fromxy2Camerapose(m_tiePointsPix[i].x3, m_tiePointsPix[i].y3);
			V_1 = XX / ZZ - XYZbody[0] / XYZbody[2];
			V_2 = YY / ZZ - XYZbody[1] / XYZbody[2];

			vError += abs(V_1);
			vError += abs(V_2);

			L(6 * i + 4) = -V_1;
			L(6 * i + 5) = -V_2;

			dXXdphi = Da1Dphi * X_ + Db1Dphi * Y_ + Dc1Dphi * Z_;
			dXXdomega = Da1Domega * X_ + Db1Domega * Y_ + Dc1Domega * Z_;
			dXXdkappa = Da1Dkappa * X_ + Db1Dkappa * Y_ + Dc1Dkappa * Z_;
			dYYdphi = Da2Dphi * X_ + Db2Dphi * Y_ + Dc2Dphi * Z_;
			dYYdomega = Da2Domega * X_ + Db2Domega * Y_ + Dc2Domega * Z_;
			dYYdkappa = Da2Dkappa * X_ + Db2Dkappa * Y_ + Dc2Dkappa * Z_;
			dZZdphi = Da3Dphi * X_ + Db3Dphi * Y_ + Dc3Dphi * Z_;
			dZZdomega = Da3Domega * X_ + Db3Domega * Y_ + Dc3Domega * Z_;
			dZZdkappa = Da3Dkappa * X_ + Db3Dkappa * Y_ + Dc3Dkappa * Z_;

			XXYYZZ_XYZ = m_cModelNd.frometCalatt(et) * m_cModelNd.get_RuMatrix();
			dXXdX = XXYYZZ_XYZ(0, 0), dYYdX = XXYYZZ_XYZ(0, 1), dZZdX = XXYYZZ_XYZ(0, 2);
			dXXdY = XXYYZZ_XYZ(1, 0), dYYdY = XXYYZZ_XYZ(1, 1), dZZdY = XXYYZZ_XYZ(1, 2);
			dXXdZ = XXYYZZ_XYZ(2, 0), dYYdZ = XXYYZZ_XYZ(2, 1), dZZdZ = XXYYZZ_XYZ(2, 2);

			coefficients_A.push_back(T(6 * i + 4, 6, (dXXdphi * ZZ - dZZdphi * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 4, 7, (dXXdomega * ZZ - dZZdomega * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 4, 8, (dXXdkappa * ZZ - dZZdkappa * XX) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 5, 6, (dYYdphi * ZZ - dZZdphi * YY) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 5, 7, (dYYdomega * ZZ - dZZdomega * YY) / ZZ / ZZ));
			coefficients_A.push_back(T(6 * i + 5, 8, (dYYdkappa * ZZ - dZZdkappa * YY) / ZZ / ZZ));

			coefficients_B.push_back(T(6 * i + 4, 3 * i + 0, (dXXdX * ZZ - dZZdX * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 4, 3 * i + 1, (dXXdY * ZZ - dZZdY * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 4, 3 * i + 2, (dXXdZ * ZZ - dZZdZ * XX) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 5, 3 * i + 0, (dYYdX * ZZ - dZZdX * YY) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 5, 3 * i + 1, (dYYdY * ZZ - dZZdY * YY) / ZZ / ZZ));
			coefficients_B.push_back(T(6 * i + 5, 3 * i + 2, (dYYdZ * ZZ - dZZdZ * YY) / ZZ / ZZ));
		}

		//完成稀疏矩阵赋值
		spt_A.setFromTriplets(coefficients_A.begin(), coefficients_A.end());
		spt_B.setFromTriplets(coefficients_B.begin(), coefficients_B.end());

		//解算
		Eigen::SparseMatrix<double> U, V, W;
		U = spt_A.transpose() * spt_A;
		V = spt_B.transpose() * spt_B;//3m*3m
		W = spt_A.transpose() * spt_B;
		Eigen::VectorXd Lx, Lt;
		Lx = spt_A.transpose() * L;
		Lt = spt_B.transpose() * L;

		Eigen::SparseMatrix<double> Vinverse(3 * tiepoints_num, 3 * tiepoints_num); //逆矩阵 维度3m*3m
		Vinverse.reserve(Eigen::VectorXd::Constant(3 * tiepoints_num, 3));//3m列，每列3个空位
		for (int i = 0; i < tiepoints_num; i++)
		{
			Eigen::Matrix3d e3 = V.block(3 * i, 3 * i, 3, 3);
			Eigen::Matrix3d e3inver = e3.inverse();

			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					Vinverse.insert(3 * i + j, 3 * i + k) = e3inver(j, k);
				}
			}
		}
		Vinverse.makeCompressed();

		Eigen::SparseMatrix<double> WVi = W * Vinverse;
		Eigen::SparseMatrix<double> WViWt = WVi * W.transpose();
		Eigen::SparseMatrix<double> M = U - WViWt;//M=U-W*(V逆)*W(转置)
		Eigen::VectorXd N = Lx - W * Vinverse * Lt;//N=Lx-W*(V逆)*Lt
		//解 M*X=N
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_x;
		solver_x.compute(M);
		Eigen::VectorXd X = solver_x.solve(N);//外方位元素改正值
		//解 V*T=S
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_t;
		Eigen::VectorXd S = Lt - W.transpose() * X;
		solver_t.compute(V);
		Eigen::VectorXd T = solver_t.solve(S);//坐标改正值

		//更新外参
		m_cModelS1.m_camera_OffsetAngleMatrix.m_extCalPar_phi += X[0] * dpr_c();
		m_cModelS1.m_camera_OffsetAngleMatrix.m_extCalPar_omega += X[1] * dpr_c();
		m_cModelS1.m_camera_OffsetAngleMatrix.m_extCalPar_kappa += X[2] * dpr_c();
		m_cModelS2.m_camera_OffsetAngleMatrix.m_extCalPar_phi += X[3] * dpr_c();
		m_cModelS2.m_camera_OffsetAngleMatrix.m_extCalPar_omega += X[4] * dpr_c();
		m_cModelS2.m_camera_OffsetAngleMatrix.m_extCalPar_kappa += X[5] * dpr_c();
		m_cModelNd.m_camera_OffsetAngleMatrix.m_extCalPar_phi += X[6] * dpr_c();
		m_cModelNd.m_camera_OffsetAngleMatrix.m_extCalPar_omega += X[7] * dpr_c();
		m_cModelNd.m_camera_OffsetAngleMatrix.m_extCalPar_kappa += X[8] * dpr_c();

		//更新目标点坐标
		for (int i = 0; i < tiepoints_num; i++)
		{
			m_tiePointsGeo[i].X += T[0 + 3 * i];
			m_tiePointsGeo[i].Y += T[1 + 3 * i];
			m_tiePointsGeo[i].Z += T[2 + 3 * i];
		}

		//
		cout <<"误差： " << vError << endl;

	} while (true);

	return true;
}

double StereoSensorModel::getHvaluefromDEM(double lat, double lon)
{
	double north_dem, east_dem;
	m_coord.transCoord(lat, lon, north_dem, east_dem);
	double h_dem;
	m_dem.getBufferValue(east_dem, north_dem, 2, &h_dem);
	return h_dem;
}


