#include "MOLA.h"

bool predictOrbitNum()
{
	int preOrbitNum;//轨道号预测值

	int initOrbitNum = 13;//初始轨道号 文件AP10013L.TAB
	
	//i864
	double initLon1 = 332.51978;//第5427行 纬度约-28.39°（影像纬度中心）
	double initLon2 = 145.97741;//第28008行 纬度约-28.39°（影像纬度中心）
	//i863
	//double initLon1 = 332.20119;//第5841行 纬度约-30.58°（影像纬度中心）
	//double initLon2 = 146.29699;//第27587行 纬度约-30.58°（影像纬度中心）

	double lonRange = 1.232;//影像经度覆盖范围1.232°
	double lonslope = 0.78;//mola轨道倾斜扫描角度
	double dh = lonRange / 2.0 + 0.78;
	//I864
	double targetAreaLon = 296.928;//影像中心经度 -63.072（+360）
	//I863
	//double targetAreaLon = 38.01;//影像中心经度 -63.072（+360）

	double offsetLon = -28.6305;//相邻轨道经度偏移量
	int num = 1000;//一组ap文件中轨道数
	cout << "开始预测10系列 轨道号：" << endl;

	for (size_t i = 1; i < num; i++)
	{
		initLon1 += offsetLon;
		initLon2 += offsetLon;
		if (initLon1 <= 0)
			initLon1 += 360;
		if (initLon2 <= 0)
			initLon2 += 360;

		if (abs(initLon1 - targetAreaLon) <= dh || abs(initLon2 - targetAreaLon) <= dh)
		{
			preOrbitNum = initOrbitNum + i;
			cout << preOrbitNum << endl;
		}
	}

	initOrbitNum = 1;//初始轨道号 文件AP11001L.TAB

	//i864
	initLon1 = 125.57929;//纬度约-28.39°
	initLon2 = 299.06893;//纬度约-28.39°
	//i863
	//initLon1 = 125.25152;//第5841行 纬度约-30.58°（影像纬度中心）
	//initLon2 = 299.39033;//第27587行 纬度约-30.58°（影像纬度中心）

	offsetLon = -28.62477764;//更新 相邻轨道经度偏移量
	
	cout << "开始预测11系列 轨道号：" << endl;

	for (size_t i = 1; i < num; i++)
	{
		initLon1 += offsetLon;
		initLon2 += offsetLon;
		if (initLon1 <= 0)
			initLon1 += 360;
		if (initLon2 <= 0)
			initLon2 += 360;

		if (abs(initLon1 - targetAreaLon) <= dh || abs(initLon2 - targetAreaLon) <= dh)
		{
			preOrbitNum = initOrbitNum + i;
			cout << preOrbitNum << endl;
		}
	}

	return true;
}

bool drawSatalliteorbit(string filepath,string outFile)
{
	FILE* fp = fopen(filepath.c_str(), "r");
	vector<pair<double, double>> orbit;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6, a7;
		if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &a1, &a2, &a3, &a4, &a5, &a6, &a7) != 7)
			continue;
		orbit.push_back(pair<double,double>(a1, a2));
	}
	fclose(fp);

	//绘制
	int numPoints = orbit.size();
	Vec lon(numPoints), lat(numPoints);
	for (size_t i = 0; i < numPoints; i++)
	{
		lon[i] = orbit[i].first;
		lat[i] = orbit[i].second;
	}

	Plot plot;
	plot.xlabel("lat");
	plot.ylabel("lon");
	plot.legend()
		.atOutsideBottom()
		.displayHorizontal()
		.displayExpandWidthBy(2);
	plot.drawPoints(lat, lon).label("pose").pointSize(0.5);
	plot.save(outFile);

	return true;
}

void readMolaGpoints(string filereadpath, string filewritepath, string Filepath_dem)
{
	// 读取文件的经度，纬度
	FILE* fp = fopen(filereadpath.c_str(), "r");
	vector<pair<double, double>> orbit;//激光范围大于DEM范围
	vector<pair<double, double>> h;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6, a7;
		if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &a1, &a2, &a3, &a4, &a5, &a6, &a7) != 7)	//1经度 2纬度 3topgra 4测距 5dem+椭球半径 
			continue;
		//转换经度的表示方法
		if (a1 >= 180.0) {
			a1 -= 360.0;
		}
		orbit.push_back(make_pair(a1, a2));
		h.push_back(make_pair(a3, a5));
	}
	fclose(fp);

	// 读取DEM高程值
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(Filepath_dem);
	//读取影像覆盖范围
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	leftlon = trans[0];
	uplat = trans[2];
	rightlon = trans[0] + trans[1] * dem.getImageWid();
	downlat = trans[2] + trans[3] * dem.getImageHei();
	dem.readToBuffer_proj(0, leftlon, rightlon, uplat, downlat);

	for (size_t i = 0; i < orbit.size(); i++)
	{
		double h_value = 0.;
		if (orbit[i].first >= leftlon && orbit[i].first <= rightlon && orbit[i].second >= downlat && orbit[i].second <= uplat)
			dem.getBufferValue(orbit[i].first, orbit[i].second, 2, &h_value);

		h_dem.push_back(h_value);
	}

	//存储高程值
	FILE* fpwrite = fopen(filewritepath.c_str(), "w");
	fprintf(fpwrite, "#MOLA经度 MOLA纬度 MOLA高程 DEM高程 DEM_sample DEM_line\n");
	for (size_t i = 0; i < orbit.size(); i++)
	{
		if (h_dem[i]) //非0，则输出；剔除DEM无效值及空洞值；
		{
			//记录DEM的像素坐标
			double sample = 0.0;
			double line = 0.0;
			dem.getBufferPixel(orbit[i].first, orbit[i].second, sample, line);
			fprintf(fpwrite, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", orbit[i].first, orbit[i].second, h[i].second - Mars_Radius, h_dem[i], sample, line);
		}
	}
	fclose(fpwrite);
}

// 测试十一 针对某一景DEM，查找对应的MOLA点对
void read_IMG_Mola(string demFile, string MolaFile, string filePoints_out, string filePoints_out_xyz)
{
	//DEM影像的范围
	//string demFile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif"; // i864 dem
	//string demFile = "E:\\IMG_i863\\test_tiecalibration\\sensorCorrect_result\\run-DEM.tif"; // i863 dem
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

	//MOLA数据列表 
	//string MolaFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\input\\MolaLists.txt"; //  mc-18 mc-25
	//string MolaFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\input\\MolaLists.txt"; //  mc-20 mc-27
	vector<string> MolaLists;
	FILE* fp_molalist = fopen(MolaFile.c_str(), "r");
	string mola_file;
	char c_read[1024];
	while (!feof(fp_molalist)) {
		if (fgets(c_read, 1024, fp_molalist) == NULL)
			continue;
		MolaLists.push_back(c_read);
	}
	fclose(fp_molalist);

	// MOLA Points 保存路径
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\MolaPoints.txt";
	string filePoints_out_xyz = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\MolaPoints_xyz.txt";*/ // i864

	//string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\MolaPoints.txt";
	//string filePoints_out_xyz = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\MolaPoints_xyz.txt"; // i863

	FILE* fpMolaPoints = fopen(filePoints_out.c_str(), "w");
	FILE* fpMolaPoints_xyz = fopen(filePoints_out_xyz.c_str(), "w");
	fprintf(fpMolaPoints, "# MOLA经度 MOLA纬度 MOLA高程 \n");
	fprintf(fpMolaPoints_xyz, "# MOLA_X MOLA_Y MOLA_Z \n");

	for (int file_index = 0; file_index < MolaLists.size(); file_index++)
	{
		// 读取文件的经度 纬度 高程
		string molafile = MolaLists[file_index];

		int n = molafile.find_last_not_of("\n");
		molafile.erase(n + 1, molafile.size() - n);
		cout << molafile << endl;

		FILE* fp = fopen(molafile.c_str(), "r");

		char c_read2[1024];
		while (!feof(fp)) {
			if (fgets(c_read2, 1024, fp) == NULL)
				continue;
			double a1, a2, a3;
			if (fscanf(fp, "%lf%lf%lf", &a1, &a2, &a3) != 3)	//1经度 2纬度 3高程 
				continue;

			if (a1 > leftlon && a1 < rightlon && a2 > downlat && a2 < uplat)
			{
				fprintf(fpMolaPoints, "%lf\t%lf\t%lf\n", a1, a2, a3);
				double x = 0.0;
				double y = 0.0;
				double z = 0.0;
				fromlatlonalt2XYZ(a2, a1, a3, x, y, z);
				fprintf(fpMolaPoints_xyz, "%lf\t%lf\t%lf\n", x, y, z);
			}

		}
		fclose(fp);
	}

	fclose(fpMolaPoints);
	fclose(fpMolaPoints_xyz);

}

void saveCritePoints(string filePoints_in, string filePoints_out)
{
	//读取MOLA控制数据文件
	FILE* fp = fopen(filePoints_in.c_str(), "r");
	vector<pair<double, double>> Geolonlat;//MOLA激光数据的经纬度
	vector<pair<double, double>> GeoHh;//MOLA和dem高程数据
	vector<double> GeoH_mola;//MOLA高程数据
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6;
		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &a1, &a2, &a3, &a4, &a5, &a6) != 6)	//MOLA经度 MOLA纬度 MOLA高程 DEM高程 DEM_sample DEM_line 
			continue;
		Geolonlat.push_back(make_pair(a1, a2));
		GeoHh.push_back(make_pair(a3, a4));
		GeoH_mola.push_back(a3);
	}
	fclose(fp);

	vector<vector<pair<double, double>>> highSlope_lonlat_all;
	vector<vector<pair<double, double>>> highSlope_Hh_all;
	// 自动筛选小撞击坑
	int keepPoints = 15;//一段点集有15个点组成
	int subsection = 150;//分段 搜索局部最小值
	int subsection_num = Geolonlat.size() / subsection;//共分成几段
	auto GeoH_mola_it_begin = GeoH_mola.begin();
	for (size_t i = 0; i < subsection_num; i++)
	{
		cout << "第 " << i + 1 << " 段： ";
		double subsection_min = 10000.;
		auto minPosition = min_element(GeoH_mola_it_begin+i*subsection, GeoH_mola_it_begin + (i+1) * subsection);
		double min_sub = *minPosition;
		int min_index = minPosition - GeoH_mola_it_begin;
		cout << "最小高程值: " << min_sub << " 序号: " << min_index << endl;

		double subsection_max = -10000.;
		auto maxPosition = max_element(GeoH_mola_it_begin + i * subsection, GeoH_mola_it_begin + (i + 1) * subsection);
		double max_sub = *maxPosition;
		int max_index = maxPosition - GeoH_mola_it_begin;
		cout << "最大高程值: " << max_sub << " 序号: " << max_index << endl;
		max_index < keepPoints / 2 ? max_index = keepPoints / 2 : max_index = max_index;

		//计算方差
		double variance_min = 0.;
		double means_min = 0.;
		
		for (int t = 0; t < keepPoints; t++)
		{
			means_min += GeoHh[min_index - keepPoints / 2 + t].first / keepPoints;
		}
		for (int t = 0; t < keepPoints; t++)
		{
			variance_min += pow(GeoHh[min_index - keepPoints / 2 + t].first - means_min, 2) / keepPoints;
		}
		cout << "方差为： " << variance_min << endl;

		double variance_max = 0.;
		double means_max = 0.;

		for (int t = 0; t < keepPoints; t++)
		{
			means_max += GeoHh[max_index - keepPoints / 2 + t].first / keepPoints;
		}
		for (int t = 0; t < keepPoints; t++)
		{
			variance_max += pow(GeoHh[max_index - keepPoints / 2 + t].first - means_max, 2) / keepPoints;
		}
		cout << "方差为： " << variance_max << endl;

		if (variance_max>variance_min)
		{
			if (variance_max > 5000)
			{
				vector<pair<double, double>> highSlope_lonlat;
				vector<pair<double, double>> highSlope_Hh;

				for (int t = 0; t < keepPoints; t++)
				{
					highSlope_lonlat.push_back(Geolonlat[max_index - keepPoints / 2 + t]);
					highSlope_Hh.push_back(GeoHh[max_index - keepPoints / 2 + t]);
				}

				highSlope_Hh_all.push_back(highSlope_Hh);
				highSlope_lonlat_all.push_back(highSlope_lonlat);
			}
		}
		else
		{
			if (variance_min > 5000)
			{
				vector<pair<double, double>> highSlope_lonlat;
				vector<pair<double, double>> highSlope_Hh;

				for (int t = 0; t < keepPoints; t++)
				{
					highSlope_lonlat.push_back(Geolonlat[min_index - keepPoints / 2 + t]);
					highSlope_Hh.push_back(GeoHh[min_index - keepPoints / 2 + t]);
				}

				highSlope_Hh_all.push_back(highSlope_Hh);
				highSlope_lonlat_all.push_back(highSlope_lonlat);
			}
		}
	}

	//保存一份高坡度点集
	FILE* fpSlopePoints = fopen(filePoints_out.c_str(), "w");
	fprintf(fpSlopePoints, "# MOLA经度 MOLA纬度 MOLA高程 DEM初始高程\n");
	for (int i = 0; i < highSlope_lonlat_all.size(); i++)
	{
		for (int j = 0; j < highSlope_lonlat_all[i].size(); j++)
		{
			fprintf(fpSlopePoints, "%lf\t%lf\t%lf\t%lf\n", highSlope_lonlat_all[i][j].first, highSlope_lonlat_all[i][j].second, highSlope_Hh_all[i][j].first, highSlope_Hh_all[i][j].second);
		}
	}
	fclose(fpSlopePoints);
}

void calNccOffset(string file,string demfile,string slopefile, string NCCfile,string BestFile)
{
	//读取MOLA控制数据文件
	FILE* fp = fopen(file.c_str(), "r");
	vector<pair<double, double>> Geolonlat;//MOLA激光数据的经纬度
	vector<pair<double, double>> GeoHh;//MOLA和dem高程数据
	vector<double> GeoH_mola;//MOLA高程数据
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6;
		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &a1, &a2, &a3, &a4, &a5, &a6) != 6)	//MOLA经度 MOLA纬度 MOLA高程 DEM高程 DEM_sample DEM_line 
			continue;
		Geolonlat.push_back(make_pair(a1, a2));
		GeoHh.push_back(make_pair(a3, a4));
		GeoH_mola.push_back(a3);
	}
	fclose(fp);

	vector<vector<pair<double, double>>> highSlope_lonlat_all;
	vector<vector<pair<double, double>>> highSlope_Hh_all;

	//1 手动选择匹配区域
	/*vector<pair<double, double>> highSlope_lonlat;
	vector<pair<double, double>> highSlope_Hh;
	for (size_t i = 0; i < 10; i++)
	{
		highSlope_lonlat.push_back(Geolonlat[1415 + i]);
		highSlope_Hh.push_back(GeoHh[1415 + i]);
	}
	highSlope_lonlat_all.push_back(highSlope_lonlat);
	highSlope_Hh_all.push_back(highSlope_Hh);*/

	/*highSlope_lonlat.resize(0);
	highSlope_Hh.resize(0);
	for (size_t i = 0; i < 8; i++)
	{
		highSlope_lonlat.push_back(Geolonlat[1362 + i]);
		highSlope_Hh.push_back(GeoHh[1362 + i]);
	}
	highSlope_lonlat_all.push_back(highSlope_lonlat);
	highSlope_Hh_all.push_back(highSlope_Hh);*/
	
	////2 自动筛选小撞击坑
	//int keepPoints = 15;//一段点集有15个点组成
	//int subsection = 50;//分段（100点） 搜索局部最小值
	//int subsection_num = Geolonlat.size() / subsection;//共分成几段
	//auto GeoH_mola_it_begin = GeoH_mola.begin();
	//for (size_t i = 0; i < subsection_num; i++)
	//{
	//	cout << "第 " << i + 1 << " 段： ";
	//	double subsection_min = 10000.;
	//	auto minPosition = min_element(GeoH_mola_it_begin+i*subsection, GeoH_mola_it_begin + (i+1) * subsection);
	//	double min_sub = *minPosition;
	//	int min_index = minPosition - GeoH_mola_it_begin;
	//	cout << "最小高程值: " << min_sub << " 序号: " << min_index << endl;

	//	double subsection_max = -10000.;
	//	auto maxPosition = max_element(GeoH_mola_it_begin + i * subsection, GeoH_mola_it_begin + (i + 1) * subsection);
	//	double max_sub = *maxPosition;
	//	int max_index = maxPosition - GeoH_mola_it_begin;
	//	cout << "最大高程值: " << max_sub << " 序号: " << max_index << endl;
	//	max_index < keepPoints / 2 ? max_index = keepPoints / 2 : max_index = max_index;

	//	//计算方差
	//	double variance_min = 0.;
	//	double means_min = 0.;
	//	
	//	for (int t = 0; t < keepPoints; t++)
	//	{
	//		means_min += GeoHh[min_index - keepPoints / 2 + t].first / keepPoints;
	//	}
	//	for (int t = 0; t < keepPoints; t++)
	//	{
	//		variance_min += pow(GeoHh[min_index - keepPoints / 2 + t].first - means_min, 2) / keepPoints;
	//	}
	//	cout << "方差为： " << variance_min << endl;

	//	double variance_max = 0.;
	//	double means_max = 0.;

	//	for (int t = 0; t < keepPoints; t++)
	//	{
	//		means_max += GeoHh[max_index - keepPoints / 2 + t].first / keepPoints;
	//	}
	//	for (int t = 0; t < keepPoints; t++)
	//	{
	//		variance_max += pow(GeoHh[max_index - keepPoints / 2 + t].first - means_max, 2) / keepPoints;
	//	}
	//	cout << "方差为： " << variance_max << endl;

	//	if (variance_max>variance_min)
	//	{
	//		if (variance_max > 5000)
	//		{
	//			vector<pair<double, double>> highSlope_lonlat;
	//			vector<pair<double, double>> highSlope_Hh;

	//			for (int t = 0; t < keepPoints; t++)
	//			{
	//				highSlope_lonlat.push_back(Geolonlat[max_index - keepPoints / 2 + t]);
	//				highSlope_Hh.push_back(GeoHh[max_index - keepPoints / 2 + t]);
	//			}

	//			highSlope_Hh_all.push_back(highSlope_Hh);
	//			highSlope_lonlat_all.push_back(highSlope_lonlat);
	//		}
	//	}
	//	else
	//	{
	//		if (variance_min > 5000)
	//		{
	//			vector<pair<double, double>> highSlope_lonlat;
	//			vector<pair<double, double>> highSlope_Hh;

	//			for (int t = 0; t < keepPoints; t++)
	//			{
	//				highSlope_lonlat.push_back(Geolonlat[min_index - keepPoints / 2 + t]);
	//				highSlope_Hh.push_back(GeoHh[min_index - keepPoints / 2 + t]);
	//			}

	//			highSlope_Hh_all.push_back(highSlope_Hh);
	//			highSlope_lonlat_all.push_back(highSlope_lonlat);
	//		}
	//	}
	//}

	//3 自动筛选坡度大的匹配区域
	//int keepPoints = 15;
	//for (int i = 1; i < Geolonlat.size() - keepPoints; i+=keepPoints)
	//{
	//	//方法一 计算坡度
	//	double Horizontal_distance = sqrt((Geolonlat[i].first - Geolonlat[i - 1].first) * (Geolonlat[i].first - Geolonlat[i - 1].first) + 
	//		(Geolonlat[i].second - Geolonlat[i - 1].second) * (Geolonlat[i].second - Geolonlat[i - 1].second)) * PI * Mars_Radius / 180.0;
	//	double  Vertical_distance = fabs(GeoHh[i].first - GeoHh[i - 1].first);
	//	double slope = Vertical_distance / Horizontal_distance * 100.0;

	//	if (slope > 10)
	//	{
	//		vector<pair<double, double>> highSlope_lonlat;
	//		vector<pair<double, double>> highSlope_Hh;

	//		for (int j = 0; j < keepPoints; j++)
	//		{
	//			highSlope_lonlat.push_back(Geolonlat[i + j]);//向前取5个像素
	//			highSlope_Hh.push_back(GeoHh[i + j]);
	//		}

	//		highSlope_Hh_all.push_back(highSlope_Hh);
	//		highSlope_lonlat_all.push_back(highSlope_lonlat);
	//		i += keepPoints;
	//	}

	//	//方法二 计算方差
	//	/*double variance = 0.;
	//	double means = 0.;
	//	double evaluateslope = 0.;
	//	double begin = GeoHh[i].first;
	//	double end = GeoHh[i + keepPoints].first;
	//	for (int t = 0; t < keepPoints; t++)
	//	{
	//		means += GeoHh[i+t].first / keepPoints;
	//	}
	//	for (int t = 0; t < keepPoints; t++)
	//	{
	//		variance += pow(GeoHh[i+t].first - means, 2) / keepPoints;
	//	}
	//	evaluateslope = variance / abs(begin - end);
	//	cout << "方差为： " << variance << endl;
	//	cout << "评价指标： " << evaluateslope << endl;
	//	if (evaluateslope > 40)
	//	{
	//		vector<pair<double, double>> highSlope_lonlat;
	//		vector<pair<double, double>> highSlope_Hh;

	//		for (int t = 0; t < keepPoints; t++)
	//		{
	//			highSlope_lonlat.push_back(Geolonlat[i + t]);
	//			highSlope_Hh.push_back(GeoHh[i + t]);
	//		}

	//		highSlope_Hh_all.push_back(highSlope_Hh);
	//		highSlope_lonlat_all.push_back(highSlope_lonlat);
	//	}*/
	//}
	
	//4、分段匹配
	int subsection = 500;//分段（50点）
	int subsection_num = Geolonlat.size() / subsection;//共分成几段
	auto GeoH_mola_it_begin = GeoH_mola.begin();
	for (size_t i = 0; i < subsection_num; i++)
	{
		vector<pair<double, double>> highSlope_lonlat;
		vector<pair<double, double>> highSlope_Hh;
		for (int t = 0; t < subsection; t++)
		{
			highSlope_lonlat.push_back(Geolonlat[i * subsection + t]);
			highSlope_Hh.push_back(GeoHh[i * subsection + t]);
		}
		highSlope_Hh_all.push_back(highSlope_Hh);
		highSlope_lonlat_all.push_back(highSlope_lonlat);
	}

	cout << "初始控制点个数: " << Geolonlat.size() << endl;
	cout << "保留高坡度的控制点组数: " << highSlope_lonlat_all.size() << endl;

	//保存一份高坡度点集
	FILE* fpSlopePoints = fopen(slopefile.c_str(), "w");
	fprintf(fpSlopePoints, "# MOLA经度 MOLA纬度 MOLA高程 DEM初始高程\n");
	for (int i = 0; i < highSlope_lonlat_all.size(); i++)
	{
		for (int j = 0; j < highSlope_lonlat_all[i].size(); j++)
		{
			fprintf(fpSlopePoints, "%lf\t%lf\t%lf\t%lf\n", highSlope_lonlat_all[i][j].first, highSlope_lonlat_all[i][j].second, highSlope_Hh_all[i][j].first, highSlope_Hh_all[i][j].second);
		}
	}
	fclose(fpSlopePoints);
	
	// 读取DEM高程值
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(demfile);
	//读取影像覆盖范围
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	leftlon = trans[0];
	uplat = trans[2];
	rightlon = trans[0] + trans[1] * dem.getImageWid();
	downlat = trans[2] + trans[3] * dem.getImageHei();
	dem.readToBuffer_proj(0, leftlon, rightlon, uplat, downlat);
	
	vector<molaOffset> cccOffset_XY;//存储每一块高坡度点集的平移量
	string BestH_all = BestFile + "all.txt";
	FILE* fpBestH_all = fopen(BestH_all.c_str(), "w");
	fprintf(fpBestH_all, "# MOLA经度 MOLA纬度 MOLA高程 DEM经度 DEM纬度 DEM高程 DEM_sample DEM_line DEM初始高程 \n");
	for (int index = 0; index < highSlope_lonlat_all.size(); index++)
	{
		//遍历周围格网内高程差值
		//DEM分辨率0.000607286274678度=35.8米
		const int gridnumX = 140;//lon
		const int gridnumY = 40;//lat 20即可
		double offsetH[gridnumY][gridnumX];//存储每个偏移点的高程偏差的皮尔逊系数
		double coefficientNcc = 0.0;
		double bestCoefficientNcc = 0.0;
		double sumCoefficientNcc = 0.0;
		int xOffset = 0;
		int yOffset = 0;
		for (int i = -1 * gridnumY / 2; i < gridnumY / 2; i++) //Y轴 纬度偏移 -100至99
		{
			for (int j = -1 * gridnumX / 2; j < gridnumX / 2; j++) //X轴 经度偏移 -100至99
			{
				//展示高程值图像
				cv::Mat out(1000, 2000, CV_8UC3, cv::Scalar::all(0));//1000行 800列

				vector<double> DEMvalue, MOLAvalue;

				//弗雷歇距离
				//vector<vector<double>> true_result;
				//vector<vector<double>> test_result;
				//FrechetDistance3D frechetDistance;
				
				double first_h_value = 0.;
				dem.getBufferValue(highSlope_lonlat_all[index][0].first + j * trans[1], highSlope_lonlat_all[index][0].second + i * trans[3], 2, &first_h_value);

				for (int k = 0; k < highSlope_lonlat_all[index].size(); k++)
				{
					double h_value = 0.;
					
					dem.getBufferValue(highSlope_lonlat_all[index][k].first + j * trans[1], highSlope_lonlat_all[index][k].second + i * trans[3], 2, &h_value);
					if (h_value) //非0，则参与计算；
					{
						DEMvalue.push_back(h_value);
						MOLAvalue.push_back(highSlope_Hh_all[index][k].first);

						cv::circle(out, cv::Point(k * 30 + 30, abs(h_value - first_h_value - 500.)), 1, cv::Scalar(0, 0, 255), 1, 4);//列 行
						cv::circle(out, cv::Point(k * 30 + 30, abs(highSlope_Hh_all[index][k].first - highSlope_Hh_all[index][0].first - 500.)), 1, cv::Scalar(0, 255, 0), 1, 4);

						//vector<double> true_point;
						//vector<double> test_point;

						/*double x, y, z;
						fromlatlonalt2XYZ(highSlope_lonlat_all[index][k].second, highSlope_lonlat_all[index][k].first, highSlope_Hh_all[index][k].first, x, y, z);
						true_point.push_back(x);
						true_point.push_back(y);
						true_point.push_back(z);

						fromlatlonalt2XYZ(highSlope_lonlat_all[index][k].second, highSlope_lonlat_all[index][k].first, h_value, x, y, z);
						test_point.push_back(x);
						test_point.push_back(y);
						test_point.push_back(z);

						true_result.push_back(true_point);
						test_result.push_back(test_point);*/
					}
				}
				
				//cv::imshow("elevation",out);
				//cv::waitKey(25);
				
				//coefficientNcc = frechetDistance.frechetDistance(true_result, test_result);

				//最佳匹配系数
				coefficientNcc = computerNCC(DEMvalue, MOLAvalue); //返回相关系数(-1到1之间)
				//cout << "j: " << j << " i: " << i << "参与计算的高程值个数： " << true_result.size() << " 相关系数：" << coefficientNcc << endl;

				sumCoefficientNcc += (coefficientNcc > 0.0 ? coefficientNcc : 0.0); // 累加系数和
				if (coefficientNcc > bestCoefficientNcc && DEMvalue.size() == highSlope_lonlat_all[index].size())
				{
					bestCoefficientNcc = coefficientNcc;
					xOffset = j;
					yOffset = i;
				}

				//存储系数值
				offsetH[i + gridnumY / 2][j + gridnumX / 2] = coefficientNcc;

			}
		}
		cv::destroyWindow("elevation");
		
		//将最佳偏移保存下来
		/*molaOffset mola_off;
		mola_off.lonOffset = xOffset;
		mola_off.latOffset = yOffset;
		mola_off.centralLon = highSlope_lonlat_all[index][keepPoints / 2].first;
		mola_off.centralLat = highSlope_lonlat_all[index][keepPoints / 2].second;
		cccOffset_XY.push_back(mola_off);*/

		// 保存成图像
		//cv::Mat image(gridnumY, gridnumX, CV_16UC1, offsetH);
		//string imagefile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\" + to_string(index) + ".png";
		//cv::imwrite(imagefile,image);

		//cout << "bestCoefficientNcc is : " << bestCoefficientNcc << " xoffset: " << xOffset << " yoffset: " << yOffset << endl;
		/*cccOffset_XY.push_back(make_pair(xOffset, yOffset));*/

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

		//计算平移量质心的方法
		/*IplImage img = cvIplImage(image);
		cv::Point offsetPoint = grayCenter(&img);
		xOffset = offsetPoint.x;
		yOffset = offsetPoint.y;*/

		//计算阈值
		double threshold = bestCoefficientNcc * gridnumX * gridnumY / sumCoefficientNcc;

		cout << "第" << index << "块点集" << endl;
		cout << "经度：" << highSlope_lonlat_all[index][0].first << " 纬度： " << highSlope_lonlat_all[index][0].second << endl;
		cout << "经度方向平移 " << xOffset << "个像素，纬度方向平移" << yOffset << "个像素" << endl;
		cout << "互相关系数最大：" << bestCoefficientNcc << endl;
		cout << "阈值：" << threshold << endl;

		if (bestCoefficientNcc> 0 && threshold > 1.0)
		{
			//存储偏移后的数据
			string bestfile = BestFile + to_string(index) + ".txt";
			FILE* fpBestH = fopen(bestfile.c_str(), "w");
			fprintf(fpBestH, "# MOLA经度 MOLA纬度 MOLA高程 DEM经度 DEM纬度 DEM高程 DEM_sample DEM_line DEM初始高程 \n");
			for (int k = 0; k < highSlope_lonlat_all[index].size(); k++)
			{
				double h_value = 0.;
				double lon = highSlope_lonlat_all[index][k].first + xOffset * trans[1];
				double lat = highSlope_lonlat_all[index][k].second + yOffset * trans[3];
				dem.getBufferValue(lon, lat, 2, &h_value);

				if (h_value)
				{
					double sample = 0.0;
					double line = 0.0;
					dem.getBufferPixel(lon, lat, sample, line);
					fprintf(fpBestH, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", highSlope_lonlat_all[index][k].first, highSlope_lonlat_all[index][k].second, highSlope_Hh_all[index][k].first, lon, lat, h_value, sample, line, highSlope_Hh_all[index][k].second);
					fprintf(fpBestH_all, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", highSlope_lonlat_all[index][k].first, highSlope_lonlat_all[index][k].second, highSlope_Hh_all[index][k].first, lon, lat, h_value, sample, line, highSlope_Hh_all[index][k].second);
				}
			}
			fclose(fpBestH);
		}

	}
	fclose(fpBestH_all);

	//读取点集，拟合曲线
	//vector<cv::Point2f> lon_dem2mola;
	//vector<cv::Point2f> lat_dem2mola;
	//FILE* fps1 = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\s1_1.txt", "r");
	//char c[1024];
	//while (!feof(fps1)) {
	//	if (fgets(c, 1024, fps1) == NULL)
	//		continue;
	//	double a1, a2, a3, a4;
	//	if (fscanf(fps1, "%lf %lf %lf %lf", &a1, &a2, &a3, &a4) != 4)
	//		continue;
	//	lon_dem2mola.push_back(cv::Point2f(a2, a4));
	//	lat_dem2mola.push_back(cv::Point2f(a1, a3));
	//}
	//fclose(fps1);

	//int n = 1;
	//cv::Mat_<float> mat_x = polyfit(lon_dem2mola, n);
	//cv::Mat_<float> mat_y = polyfit(lat_dem2mola, n);
	
	////存储 BestH.dat
	//string bestfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\BestH.dat";
	//FILE* fpBestH = fopen(bestfile.c_str(), "w");
	//fprintf(fpBestH, "# MOLA经度 MOLA纬度 MOLA高程 DEM经度 DEM纬度 DEM高程 DEM_sample DEM_line\n");
	//for (size_t k = 0; k < Geolonlat.size(); k++)
	//{
	//	double h_value = 0.;
	//	double lat, lon, demH, molaH;
	//	double lon_off = mat_x(0, 0) + Geolonlat[k].first * mat_x(1, 0);
	//	double lat_off = mat_y(0, 0) + Geolonlat[k].second * mat_y(1, 0);
	//	lat = Geolonlat[k].second + lat_off * trans[3];
	//	lon = Geolonlat[k].first + lon_off * trans[1];
	//	dem.getBufferValue(lon, lat, 2, &h_value);
	//	if (h_value)
	//	{
	//		//记录平移后DEM的像素坐标
	//		double sample = 0.0;
	//		double line = 0.0;
	//		dem.getBufferPixel(lon, lat, sample, line);
	//		fprintf(fpBestH, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Geolonlat[k].first, Geolonlat[k].second, GeoHh[k].first, lon, lat, h_value, sample, line);
	//	}
	//}
	//fclose(fpBestH);

	////存储 高坡度区域高程对比图
	//FILE* fpH = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\h.txt", "w");
	//fprintf(fpH, "MOLA纬度 MOLA高程 DEM初始高程 DEM最佳高程\n");
	////存储平移后的点集
	//for (int index = 0; index < highSlope_lonlat_all.size(); index++)
	//{
	//	for (int i = 0; i < highSlope_lonlat_all[index].size(); i++)
	//	{
	//		double h_value = 0.;
	//		double lon_off, lat_off;
	//		lon_off = mat_x(0, 0) + (double)highSlope_lonlat_all[index][i].first * mat_x(1, 0);
	//		lat_off = mat_y(0, 0) + highSlope_lonlat_all[index][i].second * mat_y(1, 0);
	//		dem.getBufferValue(highSlope_lonlat_all[index][i].first + lon_off * trans[1], highSlope_lonlat_all[index][i].second + lat_off * trans[3], 2, &h_value);
	//		if (h_value) //非0，则参与计算；
	//		{
	//			fprintf(fpH, "%lf\t%lf\t%lf\t%lf\n", highSlope_lonlat_all[index][i].second, highSlope_Hh_all[index][i].first, highSlope_Hh_all[index][i].second, h_value);
	//		}
	//	}
	//}
	//fclose(fpH);

	////确定 mola经纬度和DEM经纬度差值 与 平移窗口 之间的关系
	////cv::Mat out(500, 500, CV_8UC3, cv::Scalar::all(0));
	////circle(out, cv::Point(100, 100), 5, cv::Scalar(0, 0, 255), CV_FILLED, CV_AA);
	//vector<cv::Point2f> lon_dem2mola;
	//vector<cv::Point2f> lat_dem2mola;
	//FILE* fpoffset = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\off.txt", "w");
	//for (int i = 0; i < cccOffset_XY.size(); i++)
	//{
	//	double x1, x2, y1, y2;
	//	//x2 = highSlope_lonlat_all[i][15].first;//目标值 MOLA 经度 取一组点集的中间点
	//	//y2 = highSlope_lonlat_all[i][15].second;//目标值 MOLA 纬度 取一组点集的中间点
	//	//x1 = highSlope_lonlat_all[i][15].first + cccOffset_XY[i].first * trans[1];
	//	//y1 = highSlope_lonlat_all[i][15].second + cccOffset_XY[i].second * trans[3];
	//	//x1 = cccOffset_XY[i].centralLon - cccOffset_XY[i].lonOffset;//
	//	//y1 = cccOffset_XY[i].centralLat - cccOffset_XY[i].latOffset;
	//	//x2 = cccOffset_XY[i].centralLon;
	//	//y2 = cccOffset_XY[i].centralLat;
	//	x1 = cccOffset_XY[i].lonOffset;
	//	y1 = cccOffset_XY[i].latOffset;
	//	x2 = cccOffset_XY[i].centralLon;
	//	y2 = cccOffset_XY[i].centralLat;
	//	lon_dem2mola.push_back(cv::Point2f(x2, x1));
	//	lat_dem2mola.push_back(cv::Point2f(y2, y1));
	//	fprintf(fpoffset, "%lf\t%lf\t%lf\t%lf\n", x2, x1, y2, y1);
	//}
	//fclose(fpoffset);

	////imshow("散点图", out);
	////cv::waitKey(0);

	////n:多项式阶次
	//int n = 1;
	//cv::Mat_<float> mat_x = polyfit(lon_dem2mola, n);
	//cv::Mat_<float> mat_y = polyfit(lat_dem2mola, n);
	///*for (int i = 0; i < lon_dem2mola.size(); i++)
	//{
	//	double x = mat_x(0, 0) + (double)lon_dem2mola[i].x * mat_x(1, 0);
	//	double y = mat_y(0, 0) + (double)lat_dem2mola[i].x * mat_y(1, 0);
	//	cout << x << "\t" << y << "\n" << lon_dem2mola[i].y << "\t" << lat_dem2mola[i].y << "\n" << abs(lon_dem2mola[i].y - x) << "\t" << abs(lat_dem2mola[i].y - y) << "\n" << endl;
	//}*/

	////高程对比图
	//FILE* fpH = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\h.txt", "w");
	//fprintf(fpH, "MOLA纬度 MOLA高程 DEM初始高程 DEM最佳高程\n");
	//for (int i = 0; i < Geolonlat.size(); i++)
	//{
	//	double dem_value = 0.;
	//	double dem_lat = 0.;
	//	double dem_lon = 0.;
	//	dem_lon = mat_x(0, 0) + Geolonlat[i].first * mat_x(1, 0);
	//	dem_lat = mat_y(0, 0) + Geolonlat[i].second * mat_y(1, 0);
	//	dem.getBufferValue(Geolonlat[i].first + dem_lon * trans[1], Geolonlat[i].second + dem_lat * trans[3], 2, &dem_value);

	//	if (dem_value)
	//	{
	//		fprintf(fpH, "%lf\t%lf\t%lf\t%lf\n", Geolonlat[i].second, GeoHh[i].first, GeoHh[i].second, dem_value);
	//	}
	//}
	//fclose(fpH);
}

void calBestOffset(string MolaFile, string Nccfile, string demfile, string bestfile)
{
	//读取NCC文件
	FILE* fp = fopen(Nccfile.c_str(), "r");
	vector<Nccsampleline> NccPoints;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		Nccsampleline NccPoint;
		if (fscanf(fp, "%d %d %lf", &NccPoint.line, &NccPoint.sample, &NccPoint.coeff) != 3)	//纬度偏移 经度偏移 NCC系数
			continue;
		NccPoints.push_back(NccPoint);
	}
	fclose(fp);

	double bestCoeff = 0.0;
	int bestSample = 0;
	int bestLine = 0;
	for (int i = 0; i < NccPoints.size(); i++)
	{
		if (NccPoints[i].coeff > bestCoeff)
		{
			bestCoeff = NccPoints[i].coeff;
			bestSample = NccPoints[i].sample;
			bestLine = NccPoints[i].line;
		}
	}
	cout << "经度方向最佳偏移: " << bestSample << endl;
	cout << "纬度方向最佳偏移: " << bestLine << endl;
	cout << "最佳一致性: " << bestCoeff << endl;

	// 读取DEM高程值
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(demfile);
	//读取影像覆盖范围
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	leftlon = trans[0];
	uplat = trans[2];
	rightlon = trans[0] + trans[1] * dem.getImageWid();
	downlat = trans[2] + trans[3] * dem.getImageHei();
	dem.readToBuffer_proj(0, leftlon, rightlon, uplat, downlat);

	//读取MOLA控制数据文件
	fp = fopen(MolaFile.c_str(), "r");
	vector<pair<double, double>> Geolonlat;//MOLA激光数据的经纬度
	vector<pair<double, double>> GeoHh;//MOLA和dem高程数据
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6;
		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &a1, &a2, &a3, &a4, &a5, &a6) != 6)	//MOLA经度 MOLA纬度 MOLA高程 DEM高程 DEM_sample DEM_line 
			continue;
		Geolonlat.push_back(make_pair(a1, a2));
		GeoHh.push_back(make_pair(a3, a4));
	}
	fclose(fp);

	//存储 最佳高程结果
	FILE* fpBestH = fopen(bestfile.c_str(), "w");
	fprintf(fpBestH, "# MOLA经度 MOLA纬度 MOLA高程 DEM经度 DEM纬度 DEM高程 DEM_sample DEM_line\n");
	for (size_t k = 0; k < Geolonlat.size(); k++)
	{
		double h_value = 0.;
		double lat, lon, demH, molaH;
		lat = Geolonlat[k].second + bestLine * trans[3];
		lon = Geolonlat[k].first + bestSample * trans[1];
		dem.getBufferValue(lon, lat, 2, &h_value);
		if (h_value)
		{
			//记录平移后DEM的像素坐标
			double sample = 0.0;
			double line = 0.0;
			dem.getBufferPixel(lon, lat, sample, line);
			fprintf(fpBestH, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Geolonlat[k].first, Geolonlat[k].second, GeoHh[k].first, lon, lat, h_value, sample,line);
		}
	}
	fclose(fpBestH);
}

void match_Mola2Dem(string molafile, string DEMfile, string offsetfile, string DEMGeoCoordinatefile, string picturefile, double& initCoefficientNcc, double& bestCoefficientNcc)
{
	// MOLA点集
	vector<GeoPoint> molaPoints_all;
	FILE* fp_mola = fopen(molafile.c_str(), "r");
	char c_read1[1024];
	while (!feof(fp_mola))
	{
		if (fgets(c_read1, 1024, fp_mola) == NULL)
			continue;
		GeoPoint molaPoint;
		if (fscanf(fp_mola, "%lf%lf%lf", &molaPoint.lon, &molaPoint.lat, &molaPoint.h) != 3)
			continue;
		molaPoints_all.push_back(molaPoint);
	}
	fclose(fp_mola);

	if (molaPoints_all.size() > 15)
	{
		// DEM影像
		GeoImageIO dem;
		vector<double> h_dem;
		dem.open(DEMfile);
		double leftlon, rightlon, uplat, downlat;
		double trans[4];
		dem.getGeoTransform(trans);
		int demWid = dem.getImageWid();
		int demHei = dem.getImageHei();
		double ori_Dem_lat = trans[2], ori_Dem_lon = trans[0];
		double end_Dem_lat = trans[2] + trans[3] * demHei, end_Dem_lon = trans[0] + trans[1] * demWid;
		dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

		// 计算平移量
		// 遍历周围格网内高程值
		const int gridnumX = 140; // lon方向
		const int gridnumY = 140; // lat方向 
		double offsetH[gridnumY][gridnumX]; // 存储偏移后的高程一致相关性
		double coefficientNcc = 0.0;

		int xOffset = 0; // 最佳偏移量
		int yOffset = 0;

		for (int i = -1 * gridnumY / 2; i < gridnumY / 2; i++) //Y轴 纬度偏移 -100至99
		{
			for (int j = -1 * gridnumX / 2; j < gridnumX / 2; j++) //X轴 经度偏移 -100至99
			{
				vector<double> DEMvalue, MOLAvalue;

				for (int k = 0; k < molaPoints_all.size(); k++)
				{
					double h_value = 0.;

					dem.getBufferValue(molaPoints_all[k].lon + j * trans[1], molaPoints_all[k].lat + i * trans[3], 2, &h_value);

					if (h_value) //剔除无效值
					{
						DEMvalue.push_back(h_value);
						MOLAvalue.push_back(molaPoints_all[k].h);
					}
				}

				//if (DEMvalue.size() < molaPoints_all.size() - 1) // 若DEM点缺失，则直接相关性为零
				//	coefficientNcc = 0.0;
				//else
				//	coefficientNcc = computerNCC(DEMvalue, MOLAvalue); //返回相关系数(-1到1之间)
				coefficientNcc = computerNCC(DEMvalue, MOLAvalue);

				if (i == 0 && j == 0)
					initCoefficientNcc = coefficientNcc; // 存储初始值

				//initCoefficientNcc += coefficientNcc / gridnumX / gridnumY;//均值ccc

				if (coefficientNcc > bestCoefficientNcc)
				{
					bestCoefficientNcc = coefficientNcc;
					xOffset = j;
					yOffset = i;
				}

				//存储系数值
				offsetH[i + gridnumY / 2][j + gridnumX / 2] = coefficientNcc;
			}
		}
		cout << "一致相关系数均值： " << initCoefficientNcc << endl;
		cout << "最佳一致相关系数： " << bestCoefficientNcc << endl;
		cout << "比值： " << bestCoefficientNcc / initCoefficientNcc << endl;

		// 偏移量文件 (最好直接输出对比图，目前sciplot不能直接绘制)
		FILE* fpoffsetH = fopen(offsetfile.c_str(), "w");
		fprintf(fpoffsetH, "# lat lon 纬度偏移 经度偏移 一致性\n");
		for (int i = 0; i < gridnumY; i++)
		{
			for (int j = 0; j < gridnumX; j++)
			{
				fprintf(fpoffsetH, "%d\t%d\t%lf\n", i - gridnumY / 2, j - gridnumX / 2, offsetH[i][j]);
			}
			fprintf(fpoffsetH, "\n");
		}
		fclose(fpoffsetH);

		// 控制点文件（MOLA数据对应的平移后的DEM）
		FILE* fpcon = fopen(DEMGeoCoordinatefile.c_str(), "w");
		fprintf(fpcon, "#molalon molalat molah demlon demlat demh no no no\n");

		// MOLA和初始DEM高程对比图 
		int pointsNum_ccc = 0; // 点的个数，保证dtm的点有效
		for (int i = 0; i < molaPoints_all.size(); i++)
		{
			double dem_h0 = 0.0;
			double dem_h1 = 0.0;
			dem.getBufferValue(molaPoints_all[i].lon, molaPoints_all[i].lat, 2, &dem_h0);
			double lon = molaPoints_all[i].lon + xOffset * trans[1];
			double lat = molaPoints_all[i].lat + yOffset * trans[3];
			dem.getBufferValue(lon, lat, 2, &dem_h1);
			if (dem_h0&&dem_h1)
			{
				pointsNum_ccc++;
			}
		}
		Vec IndexAltitude(pointsNum_ccc), MOLAAltitude(pointsNum_ccc), DEMAltitude(pointsNum_ccc), DEMAltitude_offset(pointsNum_ccc);
		int index_ccc = 0;
		for (int j = 0;j< molaPoints_all.size();j++)
		{
			double dem_h0 = 0.0;
			double dem_h1 = 0.0;
			dem.getBufferValue(molaPoints_all[j].lon, molaPoints_all[j].lat, 2, &dem_h0);
			double lon = molaPoints_all[j].lon + xOffset * trans[1];
			double lat = molaPoints_all[j].lat + yOffset * trans[3];
			dem.getBufferValue(lon, lat, 2, &dem_h1);
			if (dem_h0 && dem_h1)
			{
				IndexAltitude[index_ccc] = index_ccc + 1;
				MOLAAltitude[index_ccc] = molaPoints_all[j].h;
				DEMAltitude[index_ccc] = dem_h0;
				DEMAltitude_offset[index_ccc] = dem_h1;
				fprintf(fpcon, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", molaPoints_all[j].lon, molaPoints_all[j].lat, molaPoints_all[j].h, lon, lat, DEMAltitude_offset[index_ccc], 0.0, 0.0, 0.0);
				index_ccc++;
			}
			else
				continue;
			
		}
		fclose(fpcon);

		Plot plotAltitude;
		//plotAltitude.fontSize(20);
		//plotAltitude.xlabel("Laser Point Index").fontName("times new roman").fontSize(20);
		//plotAltitude.ylabel("Altitude").fontName("times new roman").fontSize(20);
		plotAltitude.legend().hide();
		plotAltitude.xtics().fontName("times new roman").fontSize(22);
		plotAltitude.ytics().fontName("times new roman").fontSize(22).increment(200);

		/*plotAltitude.legend()
			.displayHorizontal()
			.displayExpandWidthBy(2)
			.atBottom();*/
		/*plotAltitude.drawPoints(IndexAltitude, MOLAAltitude).label("laser").pointType(7);
		plotAltitude.drawPoints(IndexAltitude, DEMAltitude).label("before").pointType(7);
		plotAltitude.drawPoints(IndexAltitude, DEMAltitude_offset).label("after").pointType(7);*/
		plotAltitude.drawCurveWithPoints(IndexAltitude, MOLAAltitude).label("laser").lineColor("red");
		plotAltitude.drawCurveWithPoints(IndexAltitude, DEMAltitude).label("before").lineColor("#04C720");
		plotAltitude.drawCurveWithPoints(IndexAltitude, DEMAltitude_offset).label("after").lineColor("blue");

		plotAltitude.save(picturefile); // 只能支持相对路径

		dem.destroy();
	}
}

// 分段匹配，计算CCC，保存不一致的点集
void saveNone_CCC_points(vector<vector<GeoPoint>> points_IMG_MOLA, string demFile, string fileMOLA_filter)
{
	//保存点集文件
	FILE* fp = fopen(fileMOLA_filter.c_str(), "w");

	//DEM影像
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
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);//读入缓存

	// 读取MOLA数据及对应的DEM数据
	//FILE* fpThreshold = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\threashold.txt", "w");
	int part_index = 0;
	for (int index_oribit = 0; index_oribit < points_IMG_MOLA.size(); index_oribit++) // 逐轨计算
	{
		cout << "计算第" << index_oribit << "轨" << endl;
		int part_points = 15; // 15个点分为一段
		int points_nums = points_IMG_MOLA[index_oribit].size(); // 此轨的点的个数
		int part_nums = points_nums / part_points; // 此轨共分成的段数

		for (int index_part = 0; index_part < part_nums; index_part++)
		{
			vector<double> DEMvalue, MOLAvalue;

			for (int i = 0; i < part_points; i++)
			{
				int index_point = index_part * part_points + i; // 点的编号

				//检索对应DEM的高程
				double h_value = 0.;
				dem.getBufferValue(points_IMG_MOLA[index_oribit][index_point].lon, points_IMG_MOLA[index_oribit][index_point].lat, 2, &h_value);

				// 读取结果为非零有效值
				if (h_value)
				{
					DEMvalue.push_back(h_value); // DEM高程
					MOLAvalue.push_back(points_IMG_MOLA[index_oribit][index_point].h); // MOLA高程
				}

			}

			// 计算方差
			part_index++;
			double variance = computeVariance(MOLAvalue);

			// 计算CCC
			double coefficientNcc = computerNCC(DEMvalue, MOLAvalue); //返回相关系数(-1到1之间)

			//fprintf(fpThreshold, "%d\t%lf\t%lf\n", part_index, variance, coefficientNcc);

			cout << "系数=" << coefficientNcc << " 点个数：" << MOLAvalue.size() << endl;

			//if (coefficientNcc < 0.5 && MOLAvalue.size() >= 10 && variance > 1000.0) // 保存相关系数小、方差大的点集
			if (MOLAvalue.size() >= 13 && variance > 1000.0) // 保存方差大的点集
			{
				for (int j = 0; j < part_points; j++)
				{
					int index_point = index_part * part_points + j; // 点的编号
					fprintf(fp, "%lf\t%lf\t%lf\n", points_IMG_MOLA[index_oribit][index_point].lon, points_IMG_MOLA[index_oribit][index_point].lat, points_IMG_MOLA[index_oribit][index_point].h);
				}
			}
		}
	}
	fclose(fp);
	//fclose(fpThreshold);
}