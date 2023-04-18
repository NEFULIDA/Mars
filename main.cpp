#pragma once
#include <iostream>
#include <omp.h>
#include <sciplot/sciplot.hpp>
#include <filesystem>

#include "gdal.h"
#include "spice.h"
#include "OpticalLinearSensorModel.h"
#include "StereoSensorModel.h"
#include "siftgpu.h"
#include "rsgeodll.h"
#include "RpcModel.h"
#include "eigenspice.h"
#include "MOLA.h"

using namespace std;
using namespace sciplot;

// 测试一 HRSC相机安装矩阵
void test_CameraInstall() 
{	
	string filecameraInstall = "E:\\date_keyanrenwu\\Mars\\MEX\\HRSC\\cameraInstall_HEAD.txt";
	MatrixcameraInstall(filecameraInstall);
}

// 测试二 数据解析
void test_Predata_MEX_HRSC()
{
	//行时路径 设定读入路径
	string sEtPath = "../data/i864/it/hi864_0000_s12.txt";
	//姿态 轨道 辅助数据 设定保存路径（我的代码通用格式）
	string spiceatt = "../data/i864/auxiliary/my/hi864_0000_s12.att";
	string spiceep = "../data/i864/auxiliary/my/hi864_0000_s12.eph";

	Predata_MEX_HRSC(spiceatt, spiceep, sEtPath);

	//转换 标准格式输出路径 （蒋哥软件及代码通用格式）
	string stdEPfile = "E:\\leeDa\\水手号峡谷区域—轨道高度500km以内\\辅助数据\\std\\hg733_0000_nd.eph";
	string stdETfile = "E:\\leeDa\\水手号峡谷区域—轨道高度500km以内\\辅助数据\\std\\hg733_0000_nd.it";
	string stdATTfile = "E:\\leeDa\\水手号峡谷区域—轨道高度500km以内\\辅助数据\\std\\hg733_0000_nd.att";
	string stdlineNumfile = "E:\\leeDa\\水手号峡谷区域—轨道高度500km以内\\辅助数据\\std\\hg733_0000_nd.txt";

	int startlineNum = 0; int endlineNum = 0;
	double startET, exposure;

	FILE* fp = fopen(sEtPath.c_str(), "r");
	char c_read[1024];
	fgets(c_read, 1024, fp);
	sscanf(c_read, "%d%lf%lf", &startlineNum, &startET, &exposure);

	while (!feof(fp)) 
	{
		fgets(c_read, 1024, fp);
		endlineNum++;
	}
	fclose(fp);

	fp = fopen(stdlineNumfile.c_str(), "w");
	fprintf(fp, "%d \t%d", startlineNum, endlineNum - 1);
	fclose(fp);

	char utcstr[1024];
	et2utc_c(startET, "C", 3, 35, utcstr);

	calcu_stdEPfile(utcstr, spiceep, stdEPfile);

	calcu_stdATTfile(utcstr, spiceatt, stdATTfile);

	calcu_stdITfile(utcstr, sEtPath, stdETfile);

}

// 测试三 构建线阵相机定位模型 （DEM辅助）输入类型字符串，输入内容 1影像路径 2行时路径 3相机索引号 4姿态数据保存路径 5轨道数据保存路径
void test_GeoImg_DEM() 
{
	// //s1影像
	//string sImgPath = "D:\\Code\\Mars\\data\\i864\\img\\hi864_0000_s12.img";
	// //辅助数据
	//string sEtPath = "D:\\Code\\Mars\\data\\i864\\auxiliary\\hi864_0000_s12.it";
	//string spiceatt = "D:\\Code\\Mars\\data\\i864\\auxiliary\\hi864_0000_s12.att";
	//string spiceep = "D:\\Code\\Mars\\data\\i864\\auxiliary\\hi864_0000_s12.eph";
	// //相机传感器
	//int CameraIndex = 1;

	 //s2影像
	string sImgPath = "../data/i864/img/hi864_0000_s22.img";
	 //辅助数据
	string sEtPath = "../data/i864/auxiliary/hi864_0000_s22.it";
	string spiceatt = "../data/i864/auxiliary/hi864_0000_s22.att";
	string spiceep = "../data/i864/auxiliary/hi864_0000_s22.eph";
	string sInnerParameterS2 = "../data/i864/calibrationParameters/hi864_0000_s22.int";
	string sOffsetParameterS2 = "../data/i864/calibrationParameters/hi864_0000_s22.ext";
	bool bExist_calibrationS2 = false;
	//相机传感器
	int CameraIndex = 2;

	//相机安装
	string sCamerainstallPath = "../data/camera/cameraInstall_HEAD.txt";

	//线阵相机传感器模型实例化
	OpticalLinearSensorModel cModel(spiceatt, sEtPath, spiceep, sCamerainstallPath, sInnerParameterS2, sOffsetParameterS2, CameraIndex, bExist_calibrationS2);

	//初始化
	string sDEMPath = "../data/mola/megr00n270hb.lbl"; //DEM投影数据，采用投影编码
	cModel.initDem(sDEMPath);
	cModel.initImg(sImgPath);
	
	//测试坐标正算
	/*
	vector<rpcGCP> points;
	FILE* fp = fopen("../data/i864/match/Gcp_Img_s2.txt", "r");
	FILE* fp_write = fopen("../data/i864/match/Gcp_Img_s2_toGeo.txt", "w");
	fprintf(fp_write, "\n");

	int time_index = 1;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		rpcGCP point;
		if (fscanf(fp, "%lf %lf %lf %lf %lf", &point.x, &point.y, &point.lat, &point.lon, &point.h) != 5)
			continue;
		cModel.fromxy2latlonaltUsDEM(point.x, point.y, point.lat, point.lon, point.h);
		fprintf(fp_write, "%lf\t%lf\t%lf\t%lf\t%lf\n", point.x, point.y, point.lat, point.lon, point.h);
		points.push_back(point);
		cout << time_index << endl;
		time_index++;
	}
	fclose(fp);
	fclose(fp_write);
	*/

	//测试坐标反算
	/*
	FILE* fp = fopen("../data/i864/match/Gcp_Img_s2_toGeo.txt", "r");
	FILE* fp_write = fopen("../data/i864/match/Gcp_Img_s2_topixel.txt", "w");
	fprintf(fp_write, "\n");
	int time_index = 1;

	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		rpcGCP point;
		if (fscanf(fp, "%lf %lf %lf %lf %lf", &point.x, &point.y, &point.lat, &point.lon, &point.h) != 5)	//sample line 纬度 经度 高程  
			continue;
		cModel.fromlatlonh2xy(point.lat, point.lon, point.h, point.x, point.y);
		fprintf(fp_write, "%lf\t%lf\t%lf\t%lf\t%lf\n", point.x, point.y, point.lat, point.lon, point.h);
		cout << time_index << endl;
		time_index++;
	}
	fclose(fp);
	fclose(fp_write);
	*/
	
	//正射纠正 （可以指定局部范围）
	/*
	vector<rpcGCP> GeoRange = cModel.getGeoRange(500, 500, 400, 400, true);//全景范围
	string sOrthImgfilePath = "../data/i864/orth/hi864_0000_s22.tif";
	cModel.OrthoCorrection(GeoRange, sOrthImgfilePath);
	*/
}

// 测试四 构建立体交会模型
void test_GeoStereo() {
	
	string sAttPathNd = "../data/i864/auxiliary/hi864_0000_nd2.att";
	string sEtPathNd = "../data/i864/auxiliary/hi864_0000_nd2.it";
	string sEpPathNd = "../data/i864/auxiliary/hi864_0000_nd2.eph";
	string sInnerParameterNd = "../data/i864/calibrationParameters/hi864_0000_nd2.int";
	string sOffsetParameterNd = "../data/i864/calibrationParameters/hi864_0000_nd2.ext";
	bool bExist_calibrationNd = true;
	int iCameraIndexNd = 0;

	string sAttPathS1 = "../data/i864/auxiliary/hi864_0000_s12.att";
	string sEtPathS1 = "../data/i864/auxiliary/hi864_0000_s12.it"; 
	string sEpPathS1 = "../data/i864/auxiliary/hi864_0000_s12.eph";
	string sInnerParameterS1 = "../data/i864/calibrationParameters/hi864_0000_s12.int";
	string sOffsetParameterS1 = "../data/i864/calibrationParameters/hi864_0000_s12.ext";
	bool bExist_calibrationS1 = true;
	int iCameraIndexS1 = 1;

	string sAttPathS2 = "../data/i864/auxiliary/hi864_0000_s22.att"; 
	string sEtPathS2 = "../data/i864/auxiliary/hi864_0000_s22.it"; 
	string sEpPathS2 = "../data/i864/auxiliary/hi864_0000_s22.eph";
	string sInnerParameterS2 = "../data/i864/calibrationParameters/hi864_0000_s22.int";
	string sOffsetParameterS2 = "../data/i864/calibrationParameters/hi864_0000_s22.ext";
	bool bExist_calibrationS2 = true;
	int iCameraIndexS2 = 2;

	string sCamerainstallPath = "../data/camera/cameraInstall_HEAD.txt";

	StereoSensorModel cModel(sAttPathNd, sEtPathNd, sEpPathNd, sInnerParameterNd, sOffsetParameterNd, iCameraIndexNd, bExist_calibrationNd, \
		sAttPathS1, sEtPathS1, sEpPathS1, sInnerParameterS1, sOffsetParameterS1,iCameraIndexS1, bExist_calibrationS1, \
		sAttPathS2, sEtPathS2, sEpPathS2, sInnerParameterS2, sOffsetParameterS2, iCameraIndexS2, bExist_calibrationS2, \
		sCamerainstallPath);

	string sDEMPath = "../data/mola/megr00n270hb.lbl";
	cModel.initDem(sDEMPath);

	string sImgNd = "../data/i864/img/hi864_0000_nd2.img"; 
	string sImgS1 = "../data/i864/img/hi864_0000_s12.img"; 
	string sImgS2 = "../data/i864/img/hi864_0000_s22.img";

	cModel.initImg(sImgNd, sImgS1, sImgS2);

	//匹配操作：计算重叠区域 划分网格 分块匹配
	/*
	GeoRange overRange = cModel.calOverlapArea3Img();//重叠区域
	double resolution = 1.0 / 6.0;//分辨率1/6°=10*10km2=300*300pix 
	cModel.segmentGrid3Img(overRange, resolution); //划分网格
	string tiepointfilepath = "../data/i864/match/s1s2nd.pts";
	cModel.gridMatch(tiepointfilepath); //匹配连接点
	*/

	//立体交会+残差剔除
	string tiepointfilepath = "../data/i864/match/s1s2nd.pts";
	cModel.readtiepoints(tiepointfilepath);
	cModel.caltieGeopose(true);
	cModel.savetiePointsPixel("../data/i864/match/");

	//定标前定位结果+残差
	/*string reportSavefile = "../data/i864/calibrationParameters/before/";
	cModel.savetiePointsGeopose(reportSavefile);
	cModel.savetiePointsResival(reportSavefile);*/
	
	//外定标
	cModel.extCalibrate();
	string reportSavefile_after = "../data/i864/calibrationParameters/after/";
	cModel.updataResidual();
	cModel.savetiePointsGeopose(reportSavefile_after);
	cModel.savetiePointsResival(reportSavefile_after);

	//内定标
}

// 测试五 激光条带数据 计算平移量
void test_MOLA() {

	//预测哪些轨道覆盖DEM
	//predictOrbitNum();
	//绘制轨道覆盖曲线
	//drawSatalliteorbit(molafile, "..\\..\\..\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\10819.pdf");
	
	//string demfile = "E:\\IMG_i863_controlPoints\\i863_mola\\run-DEM.tif";//输入 ASP生产的DEM数据
	//string molafile = "E:\\IMG_i863_controlPoints\\i863_mola\\AP10652L.TAB";//输入 MOLA数据

	//string InitGcpfile = "E:\\IMG_i863_controlPoints\\i863_mola\\绘图结果\\AP10652L.txt";//输出 
	//readMolaGpoints(molafile, InitGcpfile, demfile);//函数读取MOLA的经纬度高程，以及对应位置的DEM高程和像素坐标

	//string demfile = "E:\\IMG_i863_controlPoints\\i863_mola\\run-DEM.tif";//输入 ASP生产的DEM数据
	//string InitGcpfile = "E:\\IMG_i863_controlPoints\\i863_mola\\绘图结果\\AP1011.txt";//输入

	//string slopePoints = "E:\\IMG_i863_controlPoints\\i863_mola\\绘图结果\\slopePoints.txt";
	//string Nccfile = "E:\\IMG_i863_controlPoints\\i863_mola\\绘图结果\\offsetH.dat";
	//string BestGcpfile = "E:\\IMG_i863_controlPoints\\i863_mola\\绘图结果\\BestH.dat";

	//calNccOffset(InitGcpfile, demfile, slopePoints, Nccfile);//计算一致相关系数最佳的偏移位置
	//calBestOffset(InitGcpfile, Nccfile, demfile, BestGcpfile);//最终的MOLA数据和DEM数据

	//i864控制点检索
	//string demfile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";//输入 ASP生产的DEM数据
	////string demfile = "E:\\IMG_i864_All\\dem_calibration\\test\\sensorCorrect_result\\run-DEM.tif";//输入 ASP生产的DEM数据(DEM约束过的)
	//
	//string molafile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\AP10228L.TAB";//输入 MOLA数据

	//string InitGcpfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\AP10228L.txt";//输出 
	//readMolaGpoints(molafile, InitGcpfile, demfile);//函数读取MOLA的经纬度高程，以及对应位置的DEM高程和像素坐标

	string InitGcpfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\AP10228L.txt";
	string fileCritePoints = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\critePoints\\AP10228L.txt";
	saveCritePoints(InitGcpfile, fileCritePoints);

	//string demfile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";//输入 ASP生产的DEM数据
	//string demfile = "E:\\IMG_i864_All\\dem_calibration\\test\\sensorCorrect_result\\run-DEM.tif";//输入 ASP生产的DEM数据(DEM约束过的)
	//string InitGcpfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\AP10228L.txt";//输入

	//string slopePoints = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\slopePoints.txt";
	//string BestGcpfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\bestH\\AP10228_best";

	//string Nccfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\AP10228_offset";
	//calNccOffset(InitGcpfile, demfile, slopePoints, Nccfile, BestGcpfile);//计算一致相关系数最佳的偏移位置

	//string Nccfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\offset1.txt";
	//calBestOffset(InitGcpfile, Nccfile, demfile, BestGcpfile);//最终的MOLA数据和DEM数据
	
}

// 测试六 dem格网数据 计算平移量
void test_MOLA_grid() {

	string NCCfile = "E:\\IMG_i864_All\\molaGridMatch\\ncc\\";
	string BestHfile = "E:\\IMG_i864_All\\molaGridMatch\\bestH\\";
	
	//按照DEM范围限定MOLA数据区域
	string MOLA = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\00n270value.tif";//覆盖范围大,投影编码，无效值-32768.0
	string DEM = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	GeoImageIO dem;
	dem.open(DEM.c_str());
	double demTrans[4];
	dem.getGeoTransform(demTrans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = demTrans[2], ori_Dem_lon = demTrans[0];
	double end_Dem_lat = demTrans[2] + demTrans[3] * demHei, end_Dem_lon = demTrans[0] + demTrans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

	GeoImageIO mola;
	mola.open(MOLA.c_str());
	double molaTrans[4];
	mola.getGeoTransform(molaTrans);

	//计算DEM与MOLA重叠区域，在MOLA坐标系下的坐标
	TransCoordAPI coord;
	//coord.init(dem.getProjectionRef().c_str(), GcsMarsDTMWTK.c_str());//验证了二者变换一致
	// 奇怪这一个使用结果却会出现问题
	//coord.init(dem.getProjectionRef().c_str(), mola.getProjectionRef().c_str());
	// 因此只能用这一个
	coord.init(GcsMarsDTMWTK.c_str(), mola.getProjectionRef().c_str());
	cout <<"DEM投影信息：\n" << dem.getProjectionRef() << endl;
	cout <<"地理编码信息：\n" << GcsMarsDTMWTK << endl;
	
	double ori_Mola_west, ori_Mola_north;
	double end_Mola_east, end_Mola_south;
	
	//注意因为要临近搜索平移量，所以mola数据应小于dem区域 dem分辨率*半个=35
	//设置dem搜索范围 160*160 像素
	const int gridnumX = 160;//lon
	const int gridnumY = 40;//lat
	coord.transCoord(ori_Dem_lat + gridnumY * demTrans[3], ori_Dem_lon + gridnumX * demTrans[1], ori_Mola_north, ori_Mola_west);
	coord.transCoord(end_Dem_lat - gridnumY * demTrans[3], end_Dem_lon - gridnumX * demTrans[1], end_Mola_south, end_Mola_east);
	mola.readToBuffer_proj(0, ori_Mola_west, end_Mola_east, ori_Mola_north, end_Mola_south);

	//重叠区域mola数据行数列数
	int molaSample, molaLine;
	molaSample = (end_Mola_east - ori_Mola_west) / molaTrans[1];
	molaLine = (end_Mola_south - ori_Mola_north) / molaTrans[3];

	//计算MOLA坐标到地理坐标
	TransCoordAPI coord_geo;
	coord_geo.init(mola.getProjectionRef().c_str(), GcsMarsDTMWTK.c_str());

	//在MOLA数据中划分格网10*10，搜索范围以格网中心点周围200*200
	int molaGrid = 10;
	int molaGridNums_wid, molaGridNums_hei;
	molaGridNums_wid = molaSample / molaGrid;
	molaGridNums_hei = molaLine / molaGrid;
	cout << "将mola数据划分为： " << molaGridNums_hei << " 行， " << molaGridNums_wid << " 列的格网" << endl;

	//保留高坡度的点集
	FILE* fp_mola = fopen("E:\\IMG_i864_All\\molaGridMatch\\slope\\slope.txt", "w");
	fprintf(fp_mola, "# MOLA经度 MOLA纬度 MOLA高程 DEM高程 DEM_sample DEM_line\n");

	string bestfile_all = BestHfile + "all.txt";
	FILE* fpBestH_all = fopen(bestfile_all.c_str(), "w");
	fprintf(fpBestH_all, "# MOLA经度 MOLA纬度 MOLA高程 DEM经度 DEM纬度 DEM高程 DEM_sample DEM_line DEM初始高程 \n");

	//循环每一块MOLA格网
	for (int i = 0; i < molaGridNums_hei; i++)
	{
		for (int j = 0; j < molaGridNums_wid; j++)
		{
			cout << "执行第： " << i << " 行， 第" << j << " 列的MOLA格网" << endl;
			//本区间内mola点
			vector<GeoPoint> molaGeoPoints;
			for (int m = 0; m < molaGrid; m++)
			{
				for (int n = 0; n < molaGrid; n++)
				{
					GeoPoint molapoint;
					double x = ori_Mola_west + (j * molaGrid + n) * molaTrans[1];
					double y = ori_Mola_north + (i * molaGrid + m) * molaTrans[3];
					mola.getBufferValue(x,y,2,&molapoint.h);
					if (molapoint.h)
					{
						coord_geo.transCoord(y, x, molapoint.lat, molapoint.lon);
						molaGeoPoints.push_back(molapoint);
					}
				}
			}

			//计算高程方差
			double variance = 0.;
			double means = 0.;
			for (int t = 0; t < molaGeoPoints.size(); t++)
			{
				means += molaGeoPoints[t].h / molaGeoPoints.size();
			}
			for (int t = 0; t < molaGeoPoints.size(); t++)
			{
				variance += pow(molaGeoPoints[t].h - means, 2) / molaGeoPoints.size();
			}
			cout << "方差为： " << variance << endl;

			//方差大的区域
			if (variance>1000)
			{
				for (int t = 0; t < molaGeoPoints.size(); t++)
				{
					fprintf(fp_mola, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", molaGeoPoints[t].lon, molaGeoPoints[t].lat, molaGeoPoints[t].h, 0., 0., 0.);
				}

				//计算偏移
				double offsetH[gridnumY][gridnumX];//存储每个偏移点的高程偏差的皮尔逊系数
				double coefficientNcc = 0.0;
				double bestCoefficientNcc = 0.0;
				double sumCoefficientNcc = 0.0;
				int xOffset = 0;
				int yOffset = 0;
				for (int a = -1 * gridnumY / 2; a < gridnumY / 2; a++) //Y轴 纬度偏移 -80至79
				{
					for (int b = -1 * gridnumX / 2; b < gridnumX / 2; b++) //X轴 经度偏移 -80至79
					{
						vector<double> DEMvalue, MOLAvalue;//逐个点计算高程
						for (int k = 0; k < molaGeoPoints.size(); k++)
						{
							double h_value = 0.;
							dem.getBufferValue(molaGeoPoints[k].lon + b * demTrans[1], molaGeoPoints[k].lat + a * demTrans[3], 2, &h_value);
							if (h_value) //非0，则参与计算；
							{
								DEMvalue.push_back(h_value);
								MOLAvalue.push_back(molaGeoPoints[k].h);
							}
						}

						//最佳匹配系数
						coefficientNcc = computerNCC(DEMvalue, MOLAvalue); //返回相关系数(-1到1之间)
						sumCoefficientNcc += (coefficientNcc > 0.0 ? coefficientNcc : 0.0); // 累加系数和
						if (coefficientNcc > bestCoefficientNcc)
						{
							bestCoefficientNcc = coefficientNcc;
							xOffset = b;
							yOffset = a;
						}

						//存储系数值
						offsetH[a + gridnumY / 2][b + gridnumX / 2] = coefficientNcc;

					}
				}
				if (sumCoefficientNcc)
				{
					double threshold = bestCoefficientNcc * gridnumX * gridnumY / sumCoefficientNcc;
					cout << "最佳匹配系数：" << bestCoefficientNcc << endl;
					cout << "阈值：" << threshold << endl;
					int index = j + i * molaGridNums_wid;

					if (bestCoefficientNcc > 0.0 && threshold > 2.0)
					{

						//存储 偏移搜索窗口的高程偏差图 dat文件
						string nccfile = NCCfile + to_string(index) + ".txt";
						FILE* fpoffsetH = fopen(nccfile.c_str(), "w");
						fprintf(fpoffsetH, "# 纬度偏移 经度偏移 NCC一致性\n");
						for (int i = 0; i < gridnumY; i++)
						{
							for (int j = 0; j < gridnumX; j++)
							{
								fprintf(fpoffsetH, "%d\t%d\t%lf\n", i - gridnumY / 2, j - gridnumX / 2, offsetH[i][j]);
							}
							fprintf(fpoffsetH, "\n");
						}
						fclose(fpoffsetH);

						//存储偏移后的数据
						string bestfile = BestHfile + to_string(index) + ".txt";
						FILE* fpBestH = fopen(bestfile.c_str(), "w");
						fprintf(fpBestH, "# MOLA经度 MOLA纬度 MOLA高程 DEM经度 DEM纬度 DEM高程 DEM_sample DEM_line DEM初始高程 \n");
						for (int k = 0; k < molaGeoPoints.size(); k++)
						{
							double h_value = 0.;
							double lon = molaGeoPoints[k].lon + xOffset * demTrans[1];
							double lat = molaGeoPoints[k].lat + yOffset * demTrans[3];
							dem.getBufferValue(lon, lat, 2, &h_value);

							if (h_value)
							{
								double sample = 0.0;
								double line = 0.0;
								dem.getBufferPixel(lon, lat, sample, line);
								fprintf(fpBestH, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", molaGeoPoints[k].lon, molaGeoPoints[k].lat, molaGeoPoints[k].h, lon, lat, h_value, sample, line, 0);
								fprintf(fpBestH_all, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", molaGeoPoints[k].lon, molaGeoPoints[k].lat, molaGeoPoints[k].h, lon, lat, h_value, sample, line, 0);
							}
						}
						fclose(fpBestH);

					}
				}
				else
				{
					cout << "最佳匹配系数：" << bestCoefficientNcc << endl;
				}

			}

		}
	}
	fclose(fp_mola);
	fclose(fpBestH_all);

	//读取MOLA格网内高程值，计算DEM对应的高程值

	//高程点集匹配，保留最高系数的平移量

	//剔除相似性低的点集以及收敛性不明显的点集
}

// 测试七 拼接全球DEM及转换lbl格式为tif（执行一遍即可，不用再执行）
void test_DEM() {
	//拼接全球DEM及转换lbl格式为tif
	//read_Img_Information("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n090hb.lbl");
	GeoImageIO mars_Dem;
	string mars_Dem_file = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\mars_dem_global.tif";
	int mars_Dem_wid, mars_Dem_hei;
	string mars_Dem_wkt;
	GeoImageIO mars_dem_megr88n000hb;
	mars_dem_megr88n000hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr88n000hb.lbl");
	mars_Dem_wid = mars_dem_megr88n000hb.getImageWid() * 4;
	mars_Dem_hei = mars_dem_megr88n000hb.getImageHei() * 4;
	mars_Dem_wkt = mars_dem_megr88n000hb.getProjectionRef();
	double geoTransform_megr88n000hb[4];
	mars_dem_megr88n000hb.getGeoTransform(geoTransform_megr88n000hb);
	double desTransform[] = { geoTransform_megr88n000hb[0],geoTransform_megr88n000hb[1],0,geoTransform_megr88n000hb[2],0.0,geoTransform_megr88n000hb[3] };
	mars_Dem.create(mars_Dem_file.c_str(), "GTiff", GDT_Int16, mars_Dem_wid, mars_Dem_hei, 1, desTransform, mars_Dem_wkt);

	int wid_single = mars_dem_megr88n000hb.getImageWid();
	int hei_single = mars_dem_megr88n000hb.getImageHei();

	//第一块
	void* pData = new int16_t[__int64(wid_single) * hei_single];
	mars_dem_megr88n000hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(0, 0, wid_single, hei_single, 0, pData);

	//第二块
	GeoImageIO mars_dem_megr88n090hb;
	mars_dem_megr88n090hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr88n090hb.lbl");
	mars_dem_megr88n090hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single - 1, 0, wid_single, hei_single, 0, pData);

	//第三块
	GeoImageIO mars_dem_megr88n180hb;
	mars_dem_megr88n180hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr88n180hb.lbl");
	mars_dem_megr88n180hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 2 - 1, 0, wid_single, hei_single, 0, pData);

	//第四块
	GeoImageIO mars_dem_megr88n270hb;
	mars_dem_megr88n270hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr88n270hb.lbl");
	mars_dem_megr88n270hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 3 - 1, 0, wid_single, hei_single, 0, pData);

	//第五块
	GeoImageIO mars_dem_megr44n000hb;
	mars_dem_megr44n000hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44n000hb.lbl");
	mars_dem_megr44n000hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(0, hei_single - 1, wid_single, hei_single, 0, pData);

	//第六块
	GeoImageIO mars_dem_megr44n090hb;
	mars_dem_megr44n090hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44n090hb.lbl");
	mars_dem_megr44n090hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single - 1, hei_single - 1, wid_single, hei_single, 0, pData);

	//第七块
	GeoImageIO mars_dem_megr44n180hb;
	mars_dem_megr44n180hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44n180hb.lbl");
	mars_dem_megr44n180hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 2 - 1, hei_single - 1, wid_single, hei_single, 0, pData);

	//第八块
	GeoImageIO mars_dem_megr44n270hb;
	mars_dem_megr44n270hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44n270hb.lbl");
	mars_dem_megr44n270hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 3 - 1, hei_single - 1, wid_single, hei_single, 0, pData);

	//第九块
	GeoImageIO mars_dem_megr00n000hb;
	mars_dem_megr00n000hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n000hb.lbl");
	mars_dem_megr00n000hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(0, hei_single * 2 - 1, wid_single, hei_single, 0, pData);

	//第十块
	GeoImageIO mars_dem_megr00n090hb;
	mars_dem_megr00n090hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n090hb.lbl");
	mars_dem_megr00n090hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single - 1, hei_single * 2 - 1, wid_single, hei_single, 0, pData);

	//第十一块
	GeoImageIO mars_dem_megr00n180hb;
	mars_dem_megr00n180hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n180hb.lbl");
	mars_dem_megr00n180hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 2 - 1, hei_single * 2 - 1, wid_single, hei_single, 0, pData);

	//第十二块
	GeoImageIO mars_dem_megr00n270hb;
	mars_dem_megr00n270hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n270hb.lbl");
	mars_dem_megr00n270hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 3 - 1, hei_single * 2 - 1, wid_single, hei_single, 0, pData);

	//第十三块
	GeoImageIO mars_dem_megr44s000hb;
	mars_dem_megr44s000hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44s000hb.lbl");
	mars_dem_megr44s000hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(0, hei_single * 3 - 1, wid_single, hei_single, 0, pData);

	//第十四块
	GeoImageIO mars_dem_megr44s090hb;
	mars_dem_megr44s090hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44s090hb.lbl");
	mars_dem_megr44s090hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single - 1, hei_single * 3 - 1, wid_single, hei_single, 0, pData);

	//第十五块
	GeoImageIO mars_dem_megr44s180hb;
	mars_dem_megr44s180hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44s180hb.lbl");
	mars_dem_megr44s180hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 2 - 1, hei_single * 3 - 1, wid_single, hei_single, 0, pData);

	//第十六块
	GeoImageIO mars_dem_megr44s270hb;
	mars_dem_megr44s270hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44s270hb.lbl");
	mars_dem_megr44s270hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 3 - 1, hei_single * 3 - 1, wid_single, hei_single, 0, pData);

	delete[]pData, pData = nullptr;
	mars_Dem.destroy();
}

// 测试八 读取指定位置高程
void test_readGeoH() {
	string demFile = "../data/mola/megr00n270hb.lbl";
	string in_lonlatFile = "../data/mola/s1s2nd_Geo.pts"; 
	string out_lonlataltFile = "../data/mola/i864_alt.txt";

	// 读取文件的经度，纬度
	FILE* fp = fopen(in_lonlatFile.c_str(), "r");
	vector<pair<double, double>> orbit;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2;
		if (fscanf(fp, "%lf%lf", &a1, &a2) != 2)	//1纬度2经度
			continue;
		orbit.push_back(make_pair(a1, a2));
	}
	fclose(fp);

	// 读取DEM高程值
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(demFile);
	//读取影像覆盖范围
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	leftlon = trans[0];
	uplat = trans[2];
	rightlon = trans[0] + trans[1] * dem.getImageWid();
	downlat = trans[2] + trans[3] * dem.getImageHei();
	dem.readToBuffer_proj(0, leftlon, rightlon, uplat, downlat);
	
	TransCoordAPI coord;
	coord.init(GcsMarsDTMWTK.c_str(), dem.getProjectionRef().c_str());

	for (size_t i = 0; i < orbit.size(); i++)
	{
		double north = 0.; double east = 0.;
		coord.transCoord(orbit[i].first, orbit[i].second, north, east);
		double h_value = 0.;
		dem.getBufferValue(east, north, 0, &h_value);

		h_dem.push_back(h_value);
	}

	//存储高程值
	FILE* fpwrite = fopen(out_lonlataltFile.c_str(), "w");
	fprintf(fpwrite, "#纬度 经度 高程\n");
	for (size_t i = 0; i < orbit.size(); i++)
	{
		fprintf(fpwrite, "%lf\t%lf\t%lf\n", orbit[i].first, orbit[i].second, h_dem[i]);
	}
	fclose(fpwrite);
}

void test_IMG_Mola() {
	string demFile = "E:\\IMG_I977\\test2\\sensorCorrect_result\\run-DEM.tif";//DEM影像
	string MolaFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\input\\MolaLists.txt";//MOLA数据列表
	string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\MolaPoints.txt";//MOLA Points （经纬度高程形式）保存路径
	string filePoints_out_xyz = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\MolaPoints_xyz.txt";//MOLA Points （XYZ形式）保存路径
	read_IMG_Mola(demFile, MolaFile, filePoints_out, filePoints_out_xyz);
}

// 测试九 自动提取立体交会生产的DEM中方差较大的区域
void test_dem_slopePoints() {
	// 区域大小限定100*100像素
    // 滑动窗口10个像素
	string DEM = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	GeoImageIO dem;
	dem.open(DEM.c_str());
	double demTrans[4];
	dem.getGeoTransform(demTrans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = demTrans[2], ori_Dem_lon = demTrans[0];
	double end_Dem_lat = demTrans[2] + demTrans[3] * demHei, end_Dem_lon = demTrans[0] + demTrans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

	int windows_size = 100;
	int slip_size = 10;

	vector<rpcGCP> GeoWindow;//保存区域中心
	vector<double> GeoWindow_variance;//保存区域方差
	//全图循环
	for (int line_index = 0; line_index + windows_size + slip_size < demHei; line_index += slip_size)
	{
		for (int sample_index = 0; sample_index + windows_size + slip_size < demWid; sample_index += slip_size)
		{
			//逐块区域统计高程值
			rpcGCP GeoPoints;//区域中心
			GeoPoints.lat = ori_Dem_lat + (line_index + windows_size/2) * demTrans[3];
			GeoPoints.lon = ori_Dem_lon + (sample_index + windows_size/2) * demTrans[1];
			dem.getBufferValue(GeoPoints.lon, GeoPoints.lat, 0, &GeoPoints.h);
			GeoPoints.y = line_index + windows_size / 2;
			GeoPoints.x = sample_index + windows_size / 2;

			vector<double> H_window;
			for (int i = 0; i < windows_size; i++)
			{
				for (int j = 0; j < windows_size; j++)
				{
					double lon = ori_Dem_lon + (sample_index + j) * demTrans[1];
					double lat = ori_Dem_lat + (line_index + i) * demTrans[3];
					double h = 0.;
					dem.getBufferValue(lon, lat, 0, &h);
					if (h > -100000.0)//剔除无效值
					{
						H_window.push_back(h);
					}
				}
			}

			if (H_window.size())//高程值有效则执行
			{
				//计算高程方差
				double variance = 0.;
				double means = 0.;
				for (int t = 0; t < H_window.size(); t++)
				{
					means += H_window[t] / H_window.size();
				}
				for (int t = 0; t < H_window.size(); t++)
				{
					variance += (pow(H_window[t] - means, 2) / H_window.size());
				}
				//cout << "方差为： " << variance << endl;
				if (variance > 3000.0)
				{
					GeoWindow.push_back(GeoPoints);
					GeoWindow_variance.push_back(variance);
				}
			}
			
		}
	}

	//保存区域
	FILE* fp = fopen("E:\\IMG_i864_All\\dem_grid\\slope\\slope3000.txt", "w");
	fprintf(fp, "sample line lon lat variance > 3000.0\n");
	for (int t = 0; t < GeoWindow.size(); t++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", GeoWindow[t].x, GeoWindow[t].y, GeoWindow[t].lon, GeoWindow[t].lat,GeoWindow_variance[t]);
	}
	fclose(fp);
}

// 测试十 读取tab格式的MOLA数据
void read_Tab_Mola() 
{
	// 限定一个经度 纬度 范围
	// 保存在该范围内所有激光点的经度 纬度 高程
	
	//MC-17范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-17\\MOLA17.txt";
	double leftlon = -135.0;
	double rightlon = -90.0;
	double uplat = 0.0;
	double downlat = -30.0;*/

	//MC-18范围
	//点集保存路径
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-18\\MOLA18.txt";
	double leftlon = -90.0;
	double rightlon = -45.0;
	double uplat = 0.0;
	double downlat = -30.0;*/

	//MC-19范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-19\\MOLA19.txt";
	double leftlon = -45.0;
	double rightlon = 0.0;
	double uplat = 0.0;
	double downlat = -30.0;*/

	// MC-20范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-20\\MOLA20.txt";
	double leftlon = 0.0;
	double rightlon = 45.0;
	double uplat = 0.0;
	double downlat = -30.0;*/

	// MC-23范围
	string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-23\\MOLA23.txt";
	double leftlon = 135.0;
	double rightlon = 180.0;
	double uplat = 0.0;
	double downlat = -30.0;

	//MC-24范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-24\\MOLA24.txt";
	double leftlon = -180.0;
	double rightlon = -120.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//MC-25范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-25\\MOLA25.txt";
	double leftlon = -120.0;
	double rightlon = -60.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//MC-26范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-26\\MOLA26.txt";
	double leftlon = -60.0;
	double rightlon = 0.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//MC-27范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-27\\MOLA27.txt";
	double leftlon = 0.0;
	double rightlon = 60.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//MC-28范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-28\\MOLA28.txt";
	double leftlon = 60.0;
	double rightlon = 120.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//mc-29范围
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-29\\MOLA29.txt";
	double leftlon = 120.0;
	double rightlon = 180.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	string tabListFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA\\list.txt";//771个MOLA文件路径索引
	vector<string> tabLists;
	FILE* fp_tablist = fopen(tabListFile.c_str(), "r");
	string tab_file;
	char c_read[1024];
	while (!feof(fp_tablist)) {
		if (fgets(c_read, 1024, fp_tablist) == NULL)
			continue;
		tabLists.push_back(c_read);
	}
	fclose(fp_tablist);

	FILE* fpMolaPoints = fopen(filePoints_out.c_str(), "w");
	fprintf(fpMolaPoints, "# MOLA经度 MOLA纬度 MOLA高程 \n");

	for (int tab_index = 0; tab_index < tabLists.size(); tab_index++)
	{
		// 读取文件的经度 纬度 高程
		string tabfile = tabLists[tab_index];
		cout << tabfile << endl;
		int n = tabfile.find_last_not_of("\n");
		tabfile.erase(n + 1, tabfile.size() - n);
		FILE* fp = fopen(tabfile.c_str(), "r");
		/*vector<pair<double, double>> orbit;
		vector<double> h;*/
		char c_read2[1024];
		while (!feof(fp)) {
			if (fgets(c_read2, 1024, fp) == NULL)
				continue;
			double a1, a2, a3, a4, a5, a6, a7;
			if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &a1, &a2, &a3, &a4, &a5, &a6, &a7) != 7)	//1经度 2纬度 3topgra 4测距 5dem+椭球半径 
				continue;
			//转换经度的表示方法
			if (a1 >= 180.0) {
				a1 -= 360.0;
			}
			if (a1>leftlon&&a1<rightlon&&a2>downlat&&a2<uplat)
			{
				//orbit.push_back(make_pair(a1, a2));
				double h_a5 = a5 - 3396000.0;
				//h.push_back(h_a5);*/
				fprintf(fpMolaPoints, "%lf\t%lf\t%lf\n", a1, a2, h_a5);
			}
			
		}
		fclose(fp);
	}

	fclose(fpMolaPoints);
}

// 测试十二 读取MOLA数据，按照轨道号,一轨轨存储
vector<vector<GeoPoint>> read_orbit_Mola(string fileIMG_MOLA) {

	vector<vector<GeoPoint>> points_IMG_MOLA;//全部的MOLA点集
	vector<GeoPoint> points_orbit_MOLA;//单轨的MOLA点集
	GeoPoint point_single;//MOLA点
	double lat_old = 0.;//上一时刻的MOLA点的纬度值

	FILE* fp = fopen(fileIMG_MOLA.c_str(), "r");
	char c_read[1024];
	while (!feof(fp))
	{
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		if (fscanf(fp, "%lf %lf %lf", &point_single.lon, &point_single.lat, &point_single.h) != 3)
			continue;
		//保存数据
		else
		{
			// 如果连续两点之间纬度差值超过0.1度，存储一轨新的数据
			if (abs(point_single.lat - lat_old) >= 0.1 && points_orbit_MOLA.size() != 0)
			{
				points_IMG_MOLA.push_back(points_orbit_MOLA);//存储一轨
				vector<GeoPoint>().swap(points_orbit_MOLA);//清空
			}

			points_orbit_MOLA.push_back(point_single);
			
			lat_old = point_single.lat;
		}
	}
	fclose(fp);

	cout << "保存了" << points_IMG_MOLA.size() << "轨数据" << endl;
	return points_IMG_MOLA;
}

// 测试十三 分段匹配，计算CCC，保存不一致的点集
void test_ccc_points() {
	string fileIMG_MOLA = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\MolaPoints.txt";
	string demFile = "E:\\IMG_I977\\test2\\sensorCorrect_result\\run-DEM.tif";
	string fileMOLA_filter = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\molapointsfilter.txt";
	saveNone_CCC_points(read_orbit_Mola(fileIMG_MOLA), demFile, fileMOLA_filter);
}

// 测试十四 计算最佳平移量
void test_match_Mola2Dem() {

	// offset save picture
	string script = "C:\\Users\\leeda\\Documents\\offsetsavepicture.plt";
	FILE* fp_script = fopen(script.c_str(), "w");
	fprintf(fp_script, "set term pdfcairo\ncd 'H:\leeda\date_keyanrenwu\Mars\MGS\MOLA\PEDR\i864_mola\output\offset'\nset xlabel 'lat (pixel)'\nset ylabel 'lon (pixel)'\nset zlabel 'offset(m)'\nset size ratio 4\nunset key\nset pm3d map\n");

	string molafile = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\mola_kmeans\\";
	string DEMfile = "D:\\TEST\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";  
	string offsetfile = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\";
	string DEMGeoCoordinatefile = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\偏移后的DEM地理坐标\\";
	string picturefile = "../../../TEST2/IMG_I864/picture/";

	string DEMGeoCoordinatefile_all = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\偏移后的DEM地理坐标\\all.txt";
	FILE* fpDEMGeoCoordinatefile_all = fopen(DEMGeoCoordinatefile_all.c_str(), "w");
	fprintf(fpDEMGeoCoordinatefile_all, "#molalon molalat molah demlon demlat demh no no no\n");
	string ccc = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\ccc.txt"; // i864

	/*string molafile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\mola_kmeans\\";
	string DEMfile = "E:\\IMG_i863\\test_tiecalibration\\sensorCorrect_result\\run-DEM.tif";
	string offsetfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\offset\\";
	string DEMGeoCoordinatefile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\偏移后的DEM地理坐标\\";
	string picturefile = "../../../date_keyanrenwu/Mars/MGS/MOLA/PEDR/i863_mola/output/picture/";

	string DEMGeoCoordinatefile_all = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\偏移后的DEM地理坐标\\all.txt";
	FILE* fpDEMGeoCoordinatefile_all = fopen(DEMGeoCoordinatefile_all.c_str(), "w");
	fprintf(fpDEMGeoCoordinatefile_all, "#molalon molalat molah demlon demlat demh no no no\n");
	string ccc = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\offset\\ccc.txt";*/	//i863

	/*string molafile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\mola_kmeans\\";
	string DEMfile = "E:\\IMG_I818\\test_2\\sensorCorrect_result\\run-DEM.tif";
	string offsetfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\offset\\";
	string DEMGeoCoordinatefile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\偏移后的DEM地理坐标\\";
	string picturefile = "../../../date_keyanrenwu/Mars/MGS/MOLA/PEDR/i818_mola/output/picture/";

	string DEMGeoCoordinatefile_all = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\偏移后的DEM地理坐标\\all.txt";
	FILE* fpDEMGeoCoordinatefile_all = fopen(DEMGeoCoordinatefile_all.c_str(), "w");
	fprintf(fpDEMGeoCoordinatefile_all, "#molalon molalat molah demlon demlat demh no no no\n");
	string ccc = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\offset\\ccc.txt";*/

	/*string molafile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\mola_kmeans\\";
	string DEMfile = "E:\\IMG_I977\\test2\\sensorCorrect_result\\run-DEM.tif";
	string offsetfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\offset\\";
	string DEMGeoCoordinatefile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\偏移后的DEM地理坐标\\";
	string picturefile = "../../../date_keyanrenwu/Mars/MGS/MOLA/PEDR/i977_mola/output/picture/";

	string DEMGeoCoordinatefile_all = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\偏移后的DEM地理坐标\\all.txt";
	FILE* fpDEMGeoCoordinatefile_all = fopen(DEMGeoCoordinatefile_all.c_str(), "w");
	fprintf(fpDEMGeoCoordinatefile_all, "#molalon molalat molah demlon demlat demh no no no\n");
	string ccc = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\offset\\ccc.txt";*/ //i977

	FILE* fpccc = fopen(ccc.c_str(), "w");

	for (int i = 139; i < 140; i++)
	{
		cout << "匹配第" << i << "块" << endl;
		string indexMolafile = molafile + to_string(i) + ".txt";
		string indexOffsetfile = offsetfile + to_string(i) + ".txt"; 
		string indexDEMGeoCoordinatefile = DEMGeoCoordinatefile + to_string(i) + ".txt"; 
		string indexPicturefile = picturefile + to_string(i) + ".pdf";
		double initCoefficientNcc = 0.0;
		double bestCoefficientNcc = 0.0;
		match_Mola2Dem(indexMolafile, DEMfile, indexOffsetfile, indexDEMGeoCoordinatefile, indexPicturefile, initCoefficientNcc, bestCoefficientNcc);
		fprintf(fpccc, "%d\t%lf\t%lf\t%lf\n", i, initCoefficientNcc, bestCoefficientNcc, bestCoefficientNcc / initCoefficientNcc);
		if (bestCoefficientNcc>0.9) // 读取偏移后的DEM地理坐标
		{
			FILE* fp_mola = fopen(indexDEMGeoCoordinatefile.c_str(), "r");
			char c_read1[1024];
			while (!feof(fp_mola))
			{
				if (fgets(c_read1, 1024, fp_mola) == NULL)
					continue;
				double a1, a2, a3, a4, a5, a6, a7, a8, a9;
				if (fscanf(fp_mola, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &a1, &a2, &a3, &a4, &a5, &a6, &a7, &a8, &a9) != 9)
					continue;
				fprintf(fpDEMGeoCoordinatefile_all, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", a1, a2, a3, a4, a5, a6, a7, a8, a9);
			}
			fclose(fp_mola);

			fprintf(fp_script, "set output '%d.pdf'\n", i);
			fprintf(fp_script, "splot '%d.txt'\n", i);
			fprintf(fp_script, "set output\n");
		}
	}
	fclose(fpDEMGeoCoordinatefile_all);
	fclose(fpccc);
	fclose(fp_script);
}

// 测试十五 
void make_gnuplot_script()
{
	// 1 绘制箭头
	/*//读取点集
	string offset = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset.txt";
	vector<rpcGCP> Points;
	FILE* fp = fopen(offset.c_str(), "r");
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		rpcGCP Point;
		if (fscanf(fp, "%lf%lf%lf%lf%lf", &Point.lon, &Point.lat, &Point.x, &Point.y, &Point.h) != 5)	 
			continue;
		Points.push_back(Point);
	}
	fclose(fp);

	//脚本路径
	string script = "C:\\Users\\lee\\Documents\\drawarrow.plt";
	FILE* fp_script = fopen(script.c_str(), "w");
	fprintf(fp_script, "cd 'E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output'\n");
	fprintf(fp_script, "set xlabel 'X-LON (PIXEL)'\n");
	fprintf(fp_script, "set ylabel 'Y-LAT (PIXEL)'\n");

	double ori_lat = -22.999944910792745;
	double ori_lon = -63.687646394017960;
	double resolution_lon = 0.00057747454488099592;
	double resolution_lat = -0.00057747454488099592;

	for (int i = 0; i < Points.size(); i++)
	{
		fprintf(fp_script, "set arrow %d from %lf,%lf rto %lf,%lf nofilled\n", i + 1, (Points[i].lon - ori_lon) / resolution_lon, (Points[i].lat - ori_lat) / resolution_lat, Points[i].x, Points[i].y);
	}
	fprintf(fp_script, "plot 'offset.txt' u (($1+63.687646394017960)/0.00057747454488099592):(($2+22.999944910792745)/-0.00057747454488099592) pt 7 ps 1title 'MOLA Points'\n");
	fclose(fp_script);
	*/
	// 2 绘制所有offset图
	string script = "C:\\Users\\lee\\Documents\\offseti863.plt";
	FILE* fp_script = fopen(script.c_str(), "w");
	fprintf(fp_script, "set term pdfcairo\ncd 'E:\date_keyanrenwu\Mars\MGS\MOLA\PEDR\i863_mola\output\offset'\nset xlabel 'lat (pixel)'\nset ylabel 'lon (pixel)'\nset zlabel 'offset(m)'\nset size ratio 4\nunset key\nset pm3d map\n");
	for (int i = 0; i < 280; i++)
	{
		fprintf(fp_script, "set output '%d.pdf'\n", i);
		fprintf(fp_script, "splot '%d.txt'\n", i);
		fprintf(fp_script, "set output\n");
	}
	fclose(fp_script);
}

// 测试十二 读取STEREO DEM数据,转为XYZ，格式txt
void test_stereo_dem2xyz() {
	//DEM影像
	string demFile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	GeoImageIO dem;
	dem.open(demFile);
	//读取影像覆盖范围
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = trans[2], ori_Dem_lon = trans[0];
	double end_Dem_lat = trans[2] + trans[3] * demHei, end_Dem_lon = trans[0] + trans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

	//分块输出三维点
	//文件夹路径
	string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\stereo_dem\\";

	int index_wid_nums = 2;
	int index_hei_nums = ceil(dem.getImageHei() / dem.getImageWid() * index_wid_nums);
	int index_nums = index_wid_nums * index_hei_nums;

	for (int index = 0; index < index_nums; index++)
	{
		string filePath = filePoints_out + to_string(index) + "I864_stereo_dem_PointsXYZ.txt";
		FILE* fpMolaPoints = fopen(filePath.c_str(), "w");

		//起始行号 列号
		int ori_line = (index / index_wid_nums) * (dem.getImageWid() / index_wid_nums);
		int ori_sample = (index % index_wid_nums) * (dem.getImageWid() / index_wid_nums);

		for (int i=ori_line; i < ori_line + dem.getImageWid() / index_wid_nums ||i > dem.getImageHei(); i++)
		{
			for (int j = ori_sample; j < ori_sample + dem.getImageWid() / index_wid_nums || j > dem.getImageWid(); j++)
			{
				double h_value = 0.;
				dem.getBufferValue(ori_Dem_lon + j * trans[1], ori_Dem_lat + i * trans[3], 0, &h_value);
				if (h_value > -20000.0)
				{
					double X = 0.0;
					double Y = 0.0;
					double Z = 0.0;
					fromlatlonalt2XYZ(ori_Dem_lat + i * trans[3], ori_Dem_lon + j * trans[1], h_value, X, Y, Z);
					fprintf(fpMolaPoints, "%lf\t%lf\t%lf\n", X, Y, Z);
				}

			}
		}
		fclose(fpMolaPoints);
	}
	
}

void test_MOLA_dem2xyz() {
	//MOLA DEM影像
	string demFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\resample\\overrange_i864_136.tif";
	GeoImageIO dem;
	dem.open(demFile);
	//读取影像覆盖范围
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = trans[2], ori_Dem_lon = trans[0];
	double end_Dem_lat = trans[2] + trans[3] * demHei, end_Dem_lon = trans[0] + trans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

	//需要转为地理坐标系
	TransCoordAPI coord;
	coord.init(dem.getProjectionRef().c_str(), GcsMarsDTMWTK.c_str());
	//计算地理范围
	double ori_Dem_lat_Geo, ori_Dem_lon_Geo;
	coord.transCoord(ori_Dem_lat, ori_Dem_lon, ori_Dem_lat_Geo, ori_Dem_lon_Geo);
	double end_Dem_lat_Geo, end_Dem_lon_Geo;
	coord.transCoord(end_Dem_lat, end_Dem_lon, end_Dem_lat_Geo, end_Dem_lon_Geo);
	double resolution_lat = (end_Dem_lat_Geo - ori_Dem_lat_Geo) / demHei;
	double resolution_lon = (end_Dem_lon_Geo - ori_Dem_lon_Geo) / demWid;

	//文件夹路径
	string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\mola_dem\\";

	int index_wid_nums = 2;
	int index_hei_nums = ceil(dem.getImageHei() / dem.getImageWid() * index_wid_nums);
	int index_nums = index_wid_nums * index_hei_nums;

	for (int index = 0; index < index_nums; index++)
	{
		string filePath = filePoints_out + to_string(index) + "I864_mola_dem_PointsXYZ.txt";
		FILE* fpMolaPoints = fopen(filePath.c_str(), "w");

		//起始行号 列号
		int ori_line = (index / index_wid_nums) * (dem.getImageWid() / index_wid_nums);
		int ori_sample = (index % index_wid_nums) * (dem.getImageWid() / index_wid_nums);
		
		for (int i = ori_line; i < ori_line + dem.getImageWid() / index_wid_nums || i > dem.getImageHei(); i++)
		{
			for (int j = ori_sample; j < ori_sample + dem.getImageWid() / index_wid_nums || j > dem.getImageWid(); j++)
			{
				double h_value = 0.;
				dem.getBufferValue(ori_Dem_lon + j * trans[1], ori_Dem_lat + i * trans[3], 0, &h_value);
				if (h_value > -20000.0)
				{
					double X = 0.0, Y = 0.0, Z = 0.0;
					double lat = ori_Dem_lat_Geo + i * resolution_lat; 
					double lon = ori_Dem_lon_Geo + j * resolution_lon;
					fromlatlonalt2XYZ(lat, lon, h_value, X, Y, Z);
					fprintf(fpMolaPoints, "%lf\t%lf\t%lf\n", X, Y, Z);
				}

			}
		}
		fclose(fpMolaPoints);
	}
}

// 测试十五 计算STEREO DEM与MOLA DEM的重叠区域，用于裁剪
void test_overrange() {
	//STEREO DEM影像
	string demFile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	GeoImageIO dem;
	dem.open(demFile);
	//读取影像覆盖范围
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = trans[2], ori_Dem_lon = trans[0];
	double end_Dem_lat = trans[2] + trans[3] * demHei, end_Dem_lon = trans[0] + trans[1] * demWid;

	//MOLA DEM影像
	string MOLAdemFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\mars_dem_global.tif";
	GeoImageIO MOLAdem;
	MOLAdem.open(MOLAdemFile);
	double MOLAtrans[4];
	MOLAdem.getGeoTransform(MOLAtrans);
	
	TransCoordAPI coord;
	// 奇怪这一个使用结果却会出现问题
	//coord.init(dem.getProjectionRef().c_str(), MOLAdem.getProjectionRef().c_str());
	// 因此只能用这一个
	coord.init(GcsMarsDTMWTK.c_str(), MOLAdem.getProjectionRef().c_str());
	double ori_MOLA_lat, ori_MOLA_lon, end_MOLA_lat, end_MOLA_lon;
	coord.transCoord(ori_Dem_lat, ori_Dem_lon, ori_MOLA_lat, ori_MOLA_lon);
	coord.transCoord(end_Dem_lat, end_Dem_lon, end_MOLA_lat, end_MOLA_lon);

	int startx = (ori_MOLA_lon - MOLAtrans[0]) / MOLAtrans[1];
	int starty = (ori_MOLA_lat - MOLAtrans[2]) / MOLAtrans[3];
	int sizex = (end_MOLA_lon - ori_MOLA_lon) / MOLAtrans[1];
	int sizey = (end_MOLA_lat - ori_MOLA_lat) / MOLAtrans[3];

	string cutfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\cutrange\\overrange_i864.tif";
	
	segment_image(MOLAdemFile, cutfile, sizex, sizey, startx, starty);
	
}

// 测试十五 读取VICAR数据
void test15() {
	string fileimg = "../data/i864/img/hi864_0000_nd2.img";
	read_VICAR_img(fileimg);
}


// 测试十六 读取DTM对应位置（一轨激光点）高程值
void test_readDTM_singleMola() {
	string file_mola = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA2\\AP10228L.TAB";
	string file_DTM_beforeCalibration = "D:\\TEST\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	string file_altitude_mola_dtm_beforeCalibration = "D:\\TEST2\\mola\\i864\\10228_before.txt";
	readMolaGpoints(file_mola, file_altitude_mola_dtm_beforeCalibration, file_DTM_beforeCalibration);

	string file_DTM_afterCalibration = "H:\\leeda\\shiyan\\sensorCorrect_result_I864_success\\run-DEM.tif";
	string file_altitude_mola_dtm_afterCalibration = "D:\\TEST2\\mola\\i864\\10228_after.txt";
	readMolaGpoints(file_mola, file_altitude_mola_dtm_afterCalibration, file_DTM_afterCalibration);
}

//804所相机安装+内方位参数生产
void makeInnerPara() {
	FILE* fp = fopen("C:\\Users\\leeda\\Desktop\\inner.txt", "w");
	double ccdlamda = 14e-6;
	double f = 1.6;

	int NumofCCDs = 33;
	fprintf(fp, "NumofCCDs = %d\n", NumofCCDs);
	int pixelsPerCCD = 1024;
	fprintf(fp, "pixelsPerCCD = %d\n", pixelsPerCCD);
	double ori = 219.9528e-3;

	int CCDID = 0;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	double ori_X = -219.9528e-3;
	double end_X = -205.6168e-3;
	double Y = -207.5220e-3;
	vector<cv::Point2f> points_tanx;
	vector<cv::Point2f> points_tany;
	int n_x = 3;	//多项式系数个数
	int n_y = 3;	//多项式系数个数
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	cv::Mat xMat = polyfit(points_tanx, n_x);
	cv::Mat yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 1;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -208.0366e-3;
	end_X = -193.7006e-3;
	Y = -193.7754e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 2;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -195.5563e-3;
	end_X = -181.2203e-3;
	Y = -207.7522e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 3;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -183.0670e-3;
	end_X = -168.7310e-3;
	Y = -193.9835e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 4;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -170.1485e-3;
	end_X = -155.8125e-3;
	Y = -207.9633e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 5;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -157.2622e-3;
	end_X = -142.9262e-3;
	Y = -194.1729e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 6;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -144.2951e-3;
	end_X = -129.9591e-3;
	Y = -208.1500e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 7;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -131.0922e-3;
	end_X = -116.7562e-3;
	Y = -194.3370e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 8;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -117.6013e-3;
	end_X = -103.2653e-3;
	Y = -208.3114e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 9;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -104.1480e-3;
	end_X = -89.8120e-3;
	Y = -194.4750e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 10;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -90.4539e-3;
	end_X = -76.1179e-3;
	Y = -208.4415e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 11;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -76.7763e-3;
	end_X = -62.4403e-3;
	Y = -194.5820e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 12;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -62.9214e-3;
	end_X = -48.5854e-3;
	Y = -208.5371e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 13;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -49.3248e-3;
	end_X = -34.9888e-3;
	Y = -194.6547e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 14;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -35.3225e-3;
	end_X = -20.9865e-3;
	Y = -208.5957e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 15;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -21.2895e-3;
	end_X = -6.9535e-3;
	Y = -194.6926e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 16;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -7.1680e-3;
	end_X = 7.1680e-3;
	Y = -208.6169e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 17;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 6.9535e-3;
	end_X = 21.2895e-3;
	Y = -194.6926e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 18;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 20.9865e-3;
	end_X = 35.3225e-3;
	Y = -208.5957e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 19;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 34.9888e-3;
	end_X = 49.3248e-3;
	Y = -194.6547e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 20;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 48.5854e-3;
	end_X = 62.9214e-3;
	Y = -208.5371e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 21;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 62.4403e-3;
	end_X = 76.7763e-3;
	Y = -194.5820e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 22;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 76.1179e-3;
	end_X = 90.4539e-3;
	Y = -208.4415e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 23;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 89.8120e-3;
	end_X = 104.1480e-3;
	Y = -194.4750e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 24;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 103.2653e-3;
	end_X = 117.6013e-3;
	Y = -208.3114e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 25;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 116.7562e-3;
	end_X = 131.0922e-3;
	Y = -194.3370e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 26;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 129.9591e-3;
	end_X = 144.2951e-3;
	Y = -208.1500e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 27;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 142.9262e-3;
	end_X = 157.2622e-3;
	Y = -194.1729e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 28;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 155.8125e-3;
	end_X = 170.1485e-3;
	Y = -207.9633e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 29;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 168.7310e-3;
	end_X = 183.0670e-3;
	Y = -193.9835e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 30;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 181.2203e-3;
	end_X = 195.5563e-3;
	Y = -207.7522e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 31;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 193.7006e-3;
	end_X = 208.0366e-3;
	Y = -193.7754e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 32;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 205.6168e-3;
	end_X = 219.9528e-3;
	Y = -207.5220e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X垂轨方向
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	fclose(fp);
}
void makeCameraINstall() {
	double seta11 = (3.0 / 60.0 + 53.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta12 = (90.0 + 2.0 / 60.0 + 36.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta13 = (89.0 + 57.0 / 60.0 + 6.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta21 = (89.0 + 57.0 / 60.0 + 24.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta22 = (2.0 / 60.0 + 37.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta23 = (90.0 + 18.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta31 = (90.0 + 2.0 / 60.0 + 54.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta32 = (89.0 + 59.0 / 60.0 + 42.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta33 = (2.0 / 60.0 + 55.0 / 60.0 / 60.0) * M_PI / 180.0;
	
	FILE* fp = fopen("C:\\Users\\leeda\\Desktop\\cameraInstall.txt", "w");//相机到body的
	fprintf(fp, "%lf\t%lf\t%lf\n", cos(seta11), cos(seta21), cos(seta31));
	fprintf(fp, "%lf\t%lf\t%lf\n", cos(seta12), cos(seta22), cos(seta32));
	fprintf(fp, "%lf\t%lf\t%lf\n", cos(seta13), cos(seta23), cos(seta33));

	fclose(fp);
}

//804测试欧拉角闭合
void testeul() {
	//1、轨道六根数转J2000系下的位置速度矢量
	/*Orbit ep_J2000;
	ep_J2000.et = 413593224.356578;
	ep_J2000.X = -6501497.423;
	ep_J2000.Y = -198485.176;
	ep_J2000.Z = 3663475.00;
	ep_J2000.Xv = 3191.535;
	ep_J2000.Yv = -3655.109;
	ep_J2000.Zv = 5467.207;*/
	EpSixElements six;
	six.a = 7471114; six.e = 0.00079; six.i = 63.444; six.AscendingNode = 165.399; six.omega = 27.0; six.M = 6.262329;
	//Orbit ep_J2000;
	//fromElments2PosVel_M(six, &ep_J2000);
	//2、轨道系到J2000的旋转矩阵
	/*SpiceDouble state[6];
	state[0] = ep_J2000.X;
	state[1] = ep_J2000.Y;
	state[2] = ep_J2000.Z;
	state[3] = ep_J2000.Xv;
	state[4] = ep_J2000.Yv;
	state[5] = ep_J2000.Zv;*/
	SpiceDouble rotationOrbit2J2000[3][3];
	//ortationFromOrbit2J2000(state, rotationOrbit2J2000);
	eul2m_c(-1.0 * rpd_c() * six.AscendingNode, -1.0 * rpd_c() * six.i, -1.0 * rpd_c() * six.omega, 3, 1, 3, rotationOrbit2J2000);//顺时针旋转
	//3、Body到J2000的四元数
	SpiceDouble		Qbody2J2000[4];
	Qbody2J2000[1] = 239408848 * pow(2, -31);
	Qbody2J2000[2] = -1737613440 * pow(2, -31);
	Qbody2J2000[3] = 430516512 * pow(2, -31);
	Qbody2J2000[0] = 1161776384 * pow(2, -31);

	SpiceDouble rotationBody2J2000[3][3];
	q2m_c(Qbody2J2000, rotationBody2J2000);
	//4、Body到轨道系 Rbody2orbit=RJ20002orbit*Rbody2J2000
	SpiceDouble rotationBody2Orbit[3][3];
	mtxm_c(rotationOrbit2J2000, rotationBody2J2000, rotationBody2Orbit);
	//5、转为欧拉角
	SpiceDouble roll, pitch, yaw;
	m2eul_c(rotationBody2Orbit, 3, 2, 1, &yaw, &pitch, &roll);//旋转轴顺序 x y z
	cout << "欧拉角： " <<"滚转：" << dpr_c() * roll << "  俯仰：" << dpr_c() * pitch << "  偏航：" << dpr_c() * yaw << endl;
}

//804测试四元数闭合
void testque() {

	//1、 由Body到轨道系的欧拉角，得到旋转矩阵Rbody2orbit；
	SpiceDouble rotationBody2Orbit[3][3];
	SpiceDouble ang3_yaw = -35286e-4, ang2_pitch = -809e-4, ang1_roll = 936e-4;//单位/度
	eul2m_c(rpd_c() * ang3_yaw, rpd_c() * ang2_pitch, rpd_c() * ang1_roll, 3, 2, 1, rotationBody2Orbit);
	//2、轨道系到J2000的旋转矩阵Rorbit2J2000；
	EpSixElements six;
	six.a = 7471114; six.e = 0.00079; six.i = 63.444; six.AscendingNode = 165.399; six.omega = 27.0; six.M = 6.262329;
	SpiceDouble rotationOrbit2J2000[3][3];
	eul2m_c(-1.0 * rpd_c() * six.AscendingNode, -1.0 * rpd_c() * six.i, -1.0 * rpd_c() * six.omega, 3, 1, 3, rotationOrbit2J2000);//顺时针旋转
	//3、Body到J2000的旋转矩阵Rbody2J2000 = Rorbit2J2000 * Rbody2orbit
	SpiceDouble rotationBody2J2000[3][3];
	mxm_c(rotationOrbit2J2000, rotationBody2Orbit, rotationBody2J2000);
	//4、Body到J2000的四元数
	SpiceDouble		Qbody2J2000[4];
	m2q_c(rotationBody2J2000, Qbody2J2000);
	//5、解析四元数
	SpiceDouble q[4];
	q[1] = 239408848 * pow(2, -31);
	q[2] = -1737613440 * pow(2, -31);
	q[3] = 430516512 * pow(2, -31);
	q[0] = -1161776384 * pow(2, -31);
	SpiceDouble Q2M[3][3];
	q2m_c(q, Q2M);
	cout << "解析四元数：" << " q1:" << 239408848 * pow(2, -31) << " q2:" << -1737613440 * pow(2, -31) << " q3:" << 430516512 * pow(2, -31) << " q4:" << -1161776384 * pow(2, -31) << endl;

	cout << "转换四元数：" << " q0:" << Qbody2J2000[0] << " q2:" << Qbody2J2000[1] << " q3:" << Qbody2J2000[2] << " q4:" << Qbody2J2000[3] << endl;
	// -0.1203(pitch) -3.4885(roll) -13.6599(yaw)
}

int main() {
	
	cout << "hello Lee, welcome to Mars!\n";

	//makeCameraINstall();
	//makeInnerPara();
	//testeul();
	//testque();
	// 测试一 HRSC相机安装矩阵
	//test_CameraInstall();

	// 测试二 数据解析
	//test_Predata_MEX_HRSC();

	// 测试三 构建线阵相机定位模型 （DEM辅助）
	//test_GeoImg_DEM();

	// 测试四 立体交会模型
	//test_GeoStereo();

	// 测试五 拟合最佳的激光数据 
	//test_MOLA();
	//test_MOLA_grid();
	//test_dem_slopePoints();

	// 测试六 拼接全球DEM及转换lbl格式为tif（执行一遍即可，不用再执行）
	//test_DEM();

	// 测试七 读取指定位置高程
	//test_readGeoH();
	//read_Img_Information("E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif");
	//read_Img_Information("E:\\date_keyanrenwu\\Mars\\MGS\\MOSAIC\\mc18-256.img");

	// 测试八 读取限定区域的经纬度高程数据
	//read_Tab_Mola();
	// 测试八 检索dem区域的mola点
	//test_IMG_Mola();

	// 测试八 读取，按照轨道存储
	//read_orbit_Mola("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\MolaPoints.txt");

	// 测试八 分段计算一致性
	//test_ccc_points();

	// 测试九 匹配MOLA和DEM的高程值
	//test_match_Mola2Dem();

	// 测试十 生产gnuplot脚本
	//make_gnuplot_script();

	// 测试十一 DEM数据LON LAT H 转为X Y Z 
	//test_stereo_dem2xyz();

	// 测试十二 GDAL影像裁剪,先计算重叠区域
	//test_overrange();

	// 测试十二 GDAL影像重采样
	//ResampleGDAL("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\cutrange\\overrange_i864.tif","E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\resample\\overrange_i864_136.tif",13.6,13.6);
	
	//test_MOLA_dem2xyz();

	//down_sample("E:\\IMG_i871\\img\\hi871_0000_s12_resample.img", "E:\\IMG_i871\\img\\hi871_0000_s12_resample2.img", 0.1);

	// 测试十三 检索激光条带
	/*predictOrbitNum();
	string filereadpath = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA2\\AP10228L.TAB";
	string filewritepath = "D:\\TEST2\\mola\\i864\\10228.txt";
	string Filepath_dem = "D:\\TEST\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	readMolaGpoints(filereadpath, filewritepath, Filepath_dem);*/

	//测试十四 生产绝对定位精度验证，读取dtm高程
	//test_readDTM_singleMola();

	//测试十五
	//test15();

	return 0; 
}
