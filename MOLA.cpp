#include "MOLA.h"

bool predictOrbitNum()
{
	int preOrbitNum;//�����Ԥ��ֵ

	int initOrbitNum = 13;//��ʼ����� �ļ�AP10013L.TAB
	
	//i864
	double initLon1 = 332.51978;//��5427�� γ��Լ-28.39�㣨Ӱ��γ�����ģ�
	double initLon2 = 145.97741;//��28008�� γ��Լ-28.39�㣨Ӱ��γ�����ģ�
	//i863
	//double initLon1 = 332.20119;//��5841�� γ��Լ-30.58�㣨Ӱ��γ�����ģ�
	//double initLon2 = 146.29699;//��27587�� γ��Լ-30.58�㣨Ӱ��γ�����ģ�

	double lonRange = 1.232;//Ӱ�񾭶ȸ��Ƿ�Χ1.232��
	double lonslope = 0.78;//mola�����бɨ��Ƕ�
	double dh = lonRange / 2.0 + 0.78;
	//I864
	double targetAreaLon = 296.928;//Ӱ�����ľ��� -63.072��+360��
	//I863
	//double targetAreaLon = 38.01;//Ӱ�����ľ��� -63.072��+360��

	double offsetLon = -28.6305;//���ڹ������ƫ����
	int num = 1000;//һ��ap�ļ��й����
	cout << "��ʼԤ��10ϵ�� ����ţ�" << endl;

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

	initOrbitNum = 1;//��ʼ����� �ļ�AP11001L.TAB

	//i864
	initLon1 = 125.57929;//γ��Լ-28.39��
	initLon2 = 299.06893;//γ��Լ-28.39��
	//i863
	//initLon1 = 125.25152;//��5841�� γ��Լ-30.58�㣨Ӱ��γ�����ģ�
	//initLon2 = 299.39033;//��27587�� γ��Լ-30.58�㣨Ӱ��γ�����ģ�

	offsetLon = -28.62477764;//���� ���ڹ������ƫ����
	
	cout << "��ʼԤ��11ϵ�� ����ţ�" << endl;

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

	//����
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
	// ��ȡ�ļ��ľ��ȣ�γ��
	FILE* fp = fopen(filereadpath.c_str(), "r");
	vector<pair<double, double>> orbit;//���ⷶΧ����DEM��Χ
	vector<pair<double, double>> h;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6, a7;
		if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &a1, &a2, &a3, &a4, &a5, &a6, &a7) != 7)	//1���� 2γ�� 3topgra 4��� 5dem+����뾶 
			continue;
		//ת�����ȵı�ʾ����
		if (a1 >= 180.0) {
			a1 -= 360.0;
		}
		orbit.push_back(make_pair(a1, a2));
		h.push_back(make_pair(a3, a5));
	}
	fclose(fp);

	// ��ȡDEM�߳�ֵ
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(Filepath_dem);
	//��ȡӰ�񸲸Ƿ�Χ
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

	//�洢�߳�ֵ
	FILE* fpwrite = fopen(filewritepath.c_str(), "w");
	fprintf(fpwrite, "#MOLA���� MOLAγ�� MOLA�߳� DEM�߳� DEM_sample DEM_line\n");
	for (size_t i = 0; i < orbit.size(); i++)
	{
		if (h_dem[i]) //��0����������޳�DEM��Чֵ���ն�ֵ��
		{
			//��¼DEM����������
			double sample = 0.0;
			double line = 0.0;
			dem.getBufferPixel(orbit[i].first, orbit[i].second, sample, line);
			fprintf(fpwrite, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", orbit[i].first, orbit[i].second, h[i].second - Mars_Radius, h_dem[i], sample, line);
		}
	}
	fclose(fpwrite);
}

// ����ʮһ ���ĳһ��DEM�����Ҷ�Ӧ��MOLA���
void read_IMG_Mola(string demFile, string MolaFile, string filePoints_out, string filePoints_out_xyz)
{
	//DEMӰ��ķ�Χ
	//string demFile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif"; // i864 dem
	//string demFile = "E:\\IMG_i863\\test_tiecalibration\\sensorCorrect_result\\run-DEM.tif"; // i863 dem
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(demFile);
	//��ȡӰ�񸲸Ƿ�Χ
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	leftlon = trans[0];
	uplat = trans[2];
	rightlon = trans[0] + trans[1] * dem.getImageWid();
	downlat = trans[2] + trans[3] * dem.getImageHei();

	//MOLA�����б� 
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

	// MOLA Points ����·��
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\MolaPoints.txt";
	string filePoints_out_xyz = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\MolaPoints_xyz.txt";*/ // i864

	//string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\MolaPoints.txt";
	//string filePoints_out_xyz = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\MolaPoints_xyz.txt"; // i863

	FILE* fpMolaPoints = fopen(filePoints_out.c_str(), "w");
	FILE* fpMolaPoints_xyz = fopen(filePoints_out_xyz.c_str(), "w");
	fprintf(fpMolaPoints, "# MOLA���� MOLAγ�� MOLA�߳� \n");
	fprintf(fpMolaPoints_xyz, "# MOLA_X MOLA_Y MOLA_Z \n");

	for (int file_index = 0; file_index < MolaLists.size(); file_index++)
	{
		// ��ȡ�ļ��ľ��� γ�� �߳�
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
			if (fscanf(fp, "%lf%lf%lf", &a1, &a2, &a3) != 3)	//1���� 2γ�� 3�߳� 
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
	//��ȡMOLA���������ļ�
	FILE* fp = fopen(filePoints_in.c_str(), "r");
	vector<pair<double, double>> Geolonlat;//MOLA�������ݵľ�γ��
	vector<pair<double, double>> GeoHh;//MOLA��dem�߳�����
	vector<double> GeoH_mola;//MOLA�߳�����
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6;
		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &a1, &a2, &a3, &a4, &a5, &a6) != 6)	//MOLA���� MOLAγ�� MOLA�߳� DEM�߳� DEM_sample DEM_line 
			continue;
		Geolonlat.push_back(make_pair(a1, a2));
		GeoHh.push_back(make_pair(a3, a4));
		GeoH_mola.push_back(a3);
	}
	fclose(fp);

	vector<vector<pair<double, double>>> highSlope_lonlat_all;
	vector<vector<pair<double, double>>> highSlope_Hh_all;
	// �Զ�ɸѡСײ����
	int keepPoints = 15;//һ�ε㼯��15�������
	int subsection = 150;//�ֶ� �����ֲ���Сֵ
	int subsection_num = Geolonlat.size() / subsection;//���ֳɼ���
	auto GeoH_mola_it_begin = GeoH_mola.begin();
	for (size_t i = 0; i < subsection_num; i++)
	{
		cout << "�� " << i + 1 << " �Σ� ";
		double subsection_min = 10000.;
		auto minPosition = min_element(GeoH_mola_it_begin+i*subsection, GeoH_mola_it_begin + (i+1) * subsection);
		double min_sub = *minPosition;
		int min_index = minPosition - GeoH_mola_it_begin;
		cout << "��С�߳�ֵ: " << min_sub << " ���: " << min_index << endl;

		double subsection_max = -10000.;
		auto maxPosition = max_element(GeoH_mola_it_begin + i * subsection, GeoH_mola_it_begin + (i + 1) * subsection);
		double max_sub = *maxPosition;
		int max_index = maxPosition - GeoH_mola_it_begin;
		cout << "���߳�ֵ: " << max_sub << " ���: " << max_index << endl;
		max_index < keepPoints / 2 ? max_index = keepPoints / 2 : max_index = max_index;

		//���㷽��
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
		cout << "����Ϊ�� " << variance_min << endl;

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
		cout << "����Ϊ�� " << variance_max << endl;

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

	//����һ�ݸ��¶ȵ㼯
	FILE* fpSlopePoints = fopen(filePoints_out.c_str(), "w");
	fprintf(fpSlopePoints, "# MOLA���� MOLAγ�� MOLA�߳� DEM��ʼ�߳�\n");
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
	//��ȡMOLA���������ļ�
	FILE* fp = fopen(file.c_str(), "r");
	vector<pair<double, double>> Geolonlat;//MOLA�������ݵľ�γ��
	vector<pair<double, double>> GeoHh;//MOLA��dem�߳�����
	vector<double> GeoH_mola;//MOLA�߳�����
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6;
		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &a1, &a2, &a3, &a4, &a5, &a6) != 6)	//MOLA���� MOLAγ�� MOLA�߳� DEM�߳� DEM_sample DEM_line 
			continue;
		Geolonlat.push_back(make_pair(a1, a2));
		GeoHh.push_back(make_pair(a3, a4));
		GeoH_mola.push_back(a3);
	}
	fclose(fp);

	vector<vector<pair<double, double>>> highSlope_lonlat_all;
	vector<vector<pair<double, double>>> highSlope_Hh_all;

	//1 �ֶ�ѡ��ƥ������
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
	
	////2 �Զ�ɸѡСײ����
	//int keepPoints = 15;//һ�ε㼯��15�������
	//int subsection = 50;//�ֶΣ�100�㣩 �����ֲ���Сֵ
	//int subsection_num = Geolonlat.size() / subsection;//���ֳɼ���
	//auto GeoH_mola_it_begin = GeoH_mola.begin();
	//for (size_t i = 0; i < subsection_num; i++)
	//{
	//	cout << "�� " << i + 1 << " �Σ� ";
	//	double subsection_min = 10000.;
	//	auto minPosition = min_element(GeoH_mola_it_begin+i*subsection, GeoH_mola_it_begin + (i+1) * subsection);
	//	double min_sub = *minPosition;
	//	int min_index = minPosition - GeoH_mola_it_begin;
	//	cout << "��С�߳�ֵ: " << min_sub << " ���: " << min_index << endl;

	//	double subsection_max = -10000.;
	//	auto maxPosition = max_element(GeoH_mola_it_begin + i * subsection, GeoH_mola_it_begin + (i + 1) * subsection);
	//	double max_sub = *maxPosition;
	//	int max_index = maxPosition - GeoH_mola_it_begin;
	//	cout << "���߳�ֵ: " << max_sub << " ���: " << max_index << endl;
	//	max_index < keepPoints / 2 ? max_index = keepPoints / 2 : max_index = max_index;

	//	//���㷽��
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
	//	cout << "����Ϊ�� " << variance_min << endl;

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
	//	cout << "����Ϊ�� " << variance_max << endl;

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

	//3 �Զ�ɸѡ�¶ȴ��ƥ������
	//int keepPoints = 15;
	//for (int i = 1; i < Geolonlat.size() - keepPoints; i+=keepPoints)
	//{
	//	//����һ �����¶�
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
	//			highSlope_lonlat.push_back(Geolonlat[i + j]);//��ǰȡ5������
	//			highSlope_Hh.push_back(GeoHh[i + j]);
	//		}

	//		highSlope_Hh_all.push_back(highSlope_Hh);
	//		highSlope_lonlat_all.push_back(highSlope_lonlat);
	//		i += keepPoints;
	//	}

	//	//������ ���㷽��
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
	//	cout << "����Ϊ�� " << variance << endl;
	//	cout << "����ָ�꣺ " << evaluateslope << endl;
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
	
	//4���ֶ�ƥ��
	int subsection = 500;//�ֶΣ�50�㣩
	int subsection_num = Geolonlat.size() / subsection;//���ֳɼ���
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

	cout << "��ʼ���Ƶ����: " << Geolonlat.size() << endl;
	cout << "�������¶ȵĿ��Ƶ�����: " << highSlope_lonlat_all.size() << endl;

	//����һ�ݸ��¶ȵ㼯
	FILE* fpSlopePoints = fopen(slopefile.c_str(), "w");
	fprintf(fpSlopePoints, "# MOLA���� MOLAγ�� MOLA�߳� DEM��ʼ�߳�\n");
	for (int i = 0; i < highSlope_lonlat_all.size(); i++)
	{
		for (int j = 0; j < highSlope_lonlat_all[i].size(); j++)
		{
			fprintf(fpSlopePoints, "%lf\t%lf\t%lf\t%lf\n", highSlope_lonlat_all[i][j].first, highSlope_lonlat_all[i][j].second, highSlope_Hh_all[i][j].first, highSlope_Hh_all[i][j].second);
		}
	}
	fclose(fpSlopePoints);
	
	// ��ȡDEM�߳�ֵ
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(demfile);
	//��ȡӰ�񸲸Ƿ�Χ
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	leftlon = trans[0];
	uplat = trans[2];
	rightlon = trans[0] + trans[1] * dem.getImageWid();
	downlat = trans[2] + trans[3] * dem.getImageHei();
	dem.readToBuffer_proj(0, leftlon, rightlon, uplat, downlat);
	
	vector<molaOffset> cccOffset_XY;//�洢ÿһ����¶ȵ㼯��ƽ����
	string BestH_all = BestFile + "all.txt";
	FILE* fpBestH_all = fopen(BestH_all.c_str(), "w");
	fprintf(fpBestH_all, "# MOLA���� MOLAγ�� MOLA�߳� DEM���� DEMγ�� DEM�߳� DEM_sample DEM_line DEM��ʼ�߳� \n");
	for (int index = 0; index < highSlope_lonlat_all.size(); index++)
	{
		//������Χ�����ڸ̲߳�ֵ
		//DEM�ֱ���0.000607286274678��=35.8��
		const int gridnumX = 140;//lon
		const int gridnumY = 40;//lat 20����
		double offsetH[gridnumY][gridnumX];//�洢ÿ��ƫ�Ƶ�ĸ߳�ƫ���Ƥ��ѷϵ��
		double coefficientNcc = 0.0;
		double bestCoefficientNcc = 0.0;
		double sumCoefficientNcc = 0.0;
		int xOffset = 0;
		int yOffset = 0;
		for (int i = -1 * gridnumY / 2; i < gridnumY / 2; i++) //Y�� γ��ƫ�� -100��99
		{
			for (int j = -1 * gridnumX / 2; j < gridnumX / 2; j++) //X�� ����ƫ�� -100��99
			{
				//չʾ�߳�ֵͼ��
				cv::Mat out(1000, 2000, CV_8UC3, cv::Scalar::all(0));//1000�� 800��

				vector<double> DEMvalue, MOLAvalue;

				//����Ъ����
				//vector<vector<double>> true_result;
				//vector<vector<double>> test_result;
				//FrechetDistance3D frechetDistance;
				
				double first_h_value = 0.;
				dem.getBufferValue(highSlope_lonlat_all[index][0].first + j * trans[1], highSlope_lonlat_all[index][0].second + i * trans[3], 2, &first_h_value);

				for (int k = 0; k < highSlope_lonlat_all[index].size(); k++)
				{
					double h_value = 0.;
					
					dem.getBufferValue(highSlope_lonlat_all[index][k].first + j * trans[1], highSlope_lonlat_all[index][k].second + i * trans[3], 2, &h_value);
					if (h_value) //��0���������㣻
					{
						DEMvalue.push_back(h_value);
						MOLAvalue.push_back(highSlope_Hh_all[index][k].first);

						cv::circle(out, cv::Point(k * 30 + 30, abs(h_value - first_h_value - 500.)), 1, cv::Scalar(0, 0, 255), 1, 4);//�� ��
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

				//���ƥ��ϵ��
				coefficientNcc = computerNCC(DEMvalue, MOLAvalue); //�������ϵ��(-1��1֮��)
				//cout << "j: " << j << " i: " << i << "�������ĸ߳�ֵ������ " << true_result.size() << " ���ϵ����" << coefficientNcc << endl;

				sumCoefficientNcc += (coefficientNcc > 0.0 ? coefficientNcc : 0.0); // �ۼ�ϵ����
				if (coefficientNcc > bestCoefficientNcc && DEMvalue.size() == highSlope_lonlat_all[index].size())
				{
					bestCoefficientNcc = coefficientNcc;
					xOffset = j;
					yOffset = i;
				}

				//�洢ϵ��ֵ
				offsetH[i + gridnumY / 2][j + gridnumX / 2] = coefficientNcc;

			}
		}
		cv::destroyWindow("elevation");
		
		//�����ƫ�Ʊ�������
		/*molaOffset mola_off;
		mola_off.lonOffset = xOffset;
		mola_off.latOffset = yOffset;
		mola_off.centralLon = highSlope_lonlat_all[index][keepPoints / 2].first;
		mola_off.centralLat = highSlope_lonlat_all[index][keepPoints / 2].second;
		cccOffset_XY.push_back(mola_off);*/

		// �����ͼ��
		//cv::Mat image(gridnumY, gridnumX, CV_16UC1, offsetH);
		//string imagefile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\" + to_string(index) + ".png";
		//cv::imwrite(imagefile,image);

		//cout << "bestCoefficientNcc is : " << bestCoefficientNcc << " xoffset: " << xOffset << " yoffset: " << yOffset << endl;
		/*cccOffset_XY.push_back(make_pair(xOffset, yOffset));*/

		//�洢 ƫ���������ڵĸ߳�ƫ��ͼ dat�ļ�
		string nccfile = NCCfile + to_string(index) + ".txt";
		FILE* fpoffsetH = fopen(nccfile.c_str(), "w");
		fprintf(fpoffsetH, "# γ��ƫ�� ����ƫ�� NCCһ����\n");
		for (int i = 0; i < gridnumY; i++)
		{
			for (int j = 0; j < gridnumX; j++)
			{
				fprintf(fpoffsetH, "%d\t%d\t%lf\n", i - gridnumY / 2, j - gridnumX / 2, offsetH[i][j]);
			}
			fprintf(fpoffsetH, "\n");
		}
		fclose(fpoffsetH);

		//����ƽ�������ĵķ���
		/*IplImage img = cvIplImage(image);
		cv::Point offsetPoint = grayCenter(&img);
		xOffset = offsetPoint.x;
		yOffset = offsetPoint.y;*/

		//������ֵ
		double threshold = bestCoefficientNcc * gridnumX * gridnumY / sumCoefficientNcc;

		cout << "��" << index << "��㼯" << endl;
		cout << "���ȣ�" << highSlope_lonlat_all[index][0].first << " γ�ȣ� " << highSlope_lonlat_all[index][0].second << endl;
		cout << "���ȷ���ƽ�� " << xOffset << "�����أ�γ�ȷ���ƽ��" << yOffset << "������" << endl;
		cout << "�����ϵ�����" << bestCoefficientNcc << endl;
		cout << "��ֵ��" << threshold << endl;

		if (bestCoefficientNcc> 0 && threshold > 1.0)
		{
			//�洢ƫ�ƺ������
			string bestfile = BestFile + to_string(index) + ".txt";
			FILE* fpBestH = fopen(bestfile.c_str(), "w");
			fprintf(fpBestH, "# MOLA���� MOLAγ�� MOLA�߳� DEM���� DEMγ�� DEM�߳� DEM_sample DEM_line DEM��ʼ�߳� \n");
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

	//��ȡ�㼯���������
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
	
	////�洢 BestH.dat
	//string bestfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\BestH.dat";
	//FILE* fpBestH = fopen(bestfile.c_str(), "w");
	//fprintf(fpBestH, "# MOLA���� MOLAγ�� MOLA�߳� DEM���� DEMγ�� DEM�߳� DEM_sample DEM_line\n");
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
	//		//��¼ƽ�ƺ�DEM����������
	//		double sample = 0.0;
	//		double line = 0.0;
	//		dem.getBufferPixel(lon, lat, sample, line);
	//		fprintf(fpBestH, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Geolonlat[k].first, Geolonlat[k].second, GeoHh[k].first, lon, lat, h_value, sample, line);
	//	}
	//}
	//fclose(fpBestH);

	////�洢 ���¶�����̶߳Ա�ͼ
	//FILE* fpH = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\h.txt", "w");
	//fprintf(fpH, "MOLAγ�� MOLA�߳� DEM��ʼ�߳� DEM��Ѹ߳�\n");
	////�洢ƽ�ƺ�ĵ㼯
	//for (int index = 0; index < highSlope_lonlat_all.size(); index++)
	//{
	//	for (int i = 0; i < highSlope_lonlat_all[index].size(); i++)
	//	{
	//		double h_value = 0.;
	//		double lon_off, lat_off;
	//		lon_off = mat_x(0, 0) + (double)highSlope_lonlat_all[index][i].first * mat_x(1, 0);
	//		lat_off = mat_y(0, 0) + highSlope_lonlat_all[index][i].second * mat_y(1, 0);
	//		dem.getBufferValue(highSlope_lonlat_all[index][i].first + lon_off * trans[1], highSlope_lonlat_all[index][i].second + lat_off * trans[3], 2, &h_value);
	//		if (h_value) //��0���������㣻
	//		{
	//			fprintf(fpH, "%lf\t%lf\t%lf\t%lf\n", highSlope_lonlat_all[index][i].second, highSlope_Hh_all[index][i].first, highSlope_Hh_all[index][i].second, h_value);
	//		}
	//	}
	//}
	//fclose(fpH);

	////ȷ�� mola��γ�Ⱥ�DEM��γ�Ȳ�ֵ �� ƽ�ƴ��� ֮��Ĺ�ϵ
	////cv::Mat out(500, 500, CV_8UC3, cv::Scalar::all(0));
	////circle(out, cv::Point(100, 100), 5, cv::Scalar(0, 0, 255), CV_FILLED, CV_AA);
	//vector<cv::Point2f> lon_dem2mola;
	//vector<cv::Point2f> lat_dem2mola;
	//FILE* fpoffset = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\off.txt", "w");
	//for (int i = 0; i < cccOffset_XY.size(); i++)
	//{
	//	double x1, x2, y1, y2;
	//	//x2 = highSlope_lonlat_all[i][15].first;//Ŀ��ֵ MOLA ���� ȡһ��㼯���м��
	//	//y2 = highSlope_lonlat_all[i][15].second;//Ŀ��ֵ MOLA γ�� ȡһ��㼯���м��
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

	////imshow("ɢ��ͼ", out);
	////cv::waitKey(0);

	////n:����ʽ�״�
	//int n = 1;
	//cv::Mat_<float> mat_x = polyfit(lon_dem2mola, n);
	//cv::Mat_<float> mat_y = polyfit(lat_dem2mola, n);
	///*for (int i = 0; i < lon_dem2mola.size(); i++)
	//{
	//	double x = mat_x(0, 0) + (double)lon_dem2mola[i].x * mat_x(1, 0);
	//	double y = mat_y(0, 0) + (double)lat_dem2mola[i].x * mat_y(1, 0);
	//	cout << x << "\t" << y << "\n" << lon_dem2mola[i].y << "\t" << lat_dem2mola[i].y << "\n" << abs(lon_dem2mola[i].y - x) << "\t" << abs(lat_dem2mola[i].y - y) << "\n" << endl;
	//}*/

	////�̶߳Ա�ͼ
	//FILE* fpH = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\h.txt", "w");
	//fprintf(fpH, "MOLAγ�� MOLA�߳� DEM��ʼ�߳� DEM��Ѹ߳�\n");
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
	//��ȡNCC�ļ�
	FILE* fp = fopen(Nccfile.c_str(), "r");
	vector<Nccsampleline> NccPoints;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		Nccsampleline NccPoint;
		if (fscanf(fp, "%d %d %lf", &NccPoint.line, &NccPoint.sample, &NccPoint.coeff) != 3)	//γ��ƫ�� ����ƫ�� NCCϵ��
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
	cout << "���ȷ������ƫ��: " << bestSample << endl;
	cout << "γ�ȷ������ƫ��: " << bestLine << endl;
	cout << "���һ����: " << bestCoeff << endl;

	// ��ȡDEM�߳�ֵ
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(demfile);
	//��ȡӰ�񸲸Ƿ�Χ
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	leftlon = trans[0];
	uplat = trans[2];
	rightlon = trans[0] + trans[1] * dem.getImageWid();
	downlat = trans[2] + trans[3] * dem.getImageHei();
	dem.readToBuffer_proj(0, leftlon, rightlon, uplat, downlat);

	//��ȡMOLA���������ļ�
	fp = fopen(MolaFile.c_str(), "r");
	vector<pair<double, double>> Geolonlat;//MOLA�������ݵľ�γ��
	vector<pair<double, double>> GeoHh;//MOLA��dem�߳�����
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2, a3, a4, a5, a6;
		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &a1, &a2, &a3, &a4, &a5, &a6) != 6)	//MOLA���� MOLAγ�� MOLA�߳� DEM�߳� DEM_sample DEM_line 
			continue;
		Geolonlat.push_back(make_pair(a1, a2));
		GeoHh.push_back(make_pair(a3, a4));
	}
	fclose(fp);

	//�洢 ��Ѹ߳̽��
	FILE* fpBestH = fopen(bestfile.c_str(), "w");
	fprintf(fpBestH, "# MOLA���� MOLAγ�� MOLA�߳� DEM���� DEMγ�� DEM�߳� DEM_sample DEM_line\n");
	for (size_t k = 0; k < Geolonlat.size(); k++)
	{
		double h_value = 0.;
		double lat, lon, demH, molaH;
		lat = Geolonlat[k].second + bestLine * trans[3];
		lon = Geolonlat[k].first + bestSample * trans[1];
		dem.getBufferValue(lon, lat, 2, &h_value);
		if (h_value)
		{
			//��¼ƽ�ƺ�DEM����������
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
	// MOLA�㼯
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
		// DEMӰ��
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

		// ����ƽ����
		// ������Χ�����ڸ߳�ֵ
		const int gridnumX = 140; // lon����
		const int gridnumY = 140; // lat���� 
		double offsetH[gridnumY][gridnumX]; // �洢ƫ�ƺ�ĸ߳�һ�������
		double coefficientNcc = 0.0;

		int xOffset = 0; // ���ƫ����
		int yOffset = 0;

		for (int i = -1 * gridnumY / 2; i < gridnumY / 2; i++) //Y�� γ��ƫ�� -100��99
		{
			for (int j = -1 * gridnumX / 2; j < gridnumX / 2; j++) //X�� ����ƫ�� -100��99
			{
				vector<double> DEMvalue, MOLAvalue;

				for (int k = 0; k < molaPoints_all.size(); k++)
				{
					double h_value = 0.;

					dem.getBufferValue(molaPoints_all[k].lon + j * trans[1], molaPoints_all[k].lat + i * trans[3], 2, &h_value);

					if (h_value) //�޳���Чֵ
					{
						DEMvalue.push_back(h_value);
						MOLAvalue.push_back(molaPoints_all[k].h);
					}
				}

				//if (DEMvalue.size() < molaPoints_all.size() - 1) // ��DEM��ȱʧ����ֱ�������Ϊ��
				//	coefficientNcc = 0.0;
				//else
				//	coefficientNcc = computerNCC(DEMvalue, MOLAvalue); //�������ϵ��(-1��1֮��)
				coefficientNcc = computerNCC(DEMvalue, MOLAvalue);

				if (i == 0 && j == 0)
					initCoefficientNcc = coefficientNcc; // �洢��ʼֵ

				//initCoefficientNcc += coefficientNcc / gridnumX / gridnumY;//��ֵccc

				if (coefficientNcc > bestCoefficientNcc)
				{
					bestCoefficientNcc = coefficientNcc;
					xOffset = j;
					yOffset = i;
				}

				//�洢ϵ��ֵ
				offsetH[i + gridnumY / 2][j + gridnumX / 2] = coefficientNcc;
			}
		}
		cout << "һ�����ϵ����ֵ�� " << initCoefficientNcc << endl;
		cout << "���һ�����ϵ���� " << bestCoefficientNcc << endl;
		cout << "��ֵ�� " << bestCoefficientNcc / initCoefficientNcc << endl;

		// ƫ�����ļ� (���ֱ������Ա�ͼ��Ŀǰsciplot����ֱ�ӻ���)
		FILE* fpoffsetH = fopen(offsetfile.c_str(), "w");
		fprintf(fpoffsetH, "# lat lon γ��ƫ�� ����ƫ�� һ����\n");
		for (int i = 0; i < gridnumY; i++)
		{
			for (int j = 0; j < gridnumX; j++)
			{
				fprintf(fpoffsetH, "%d\t%d\t%lf\n", i - gridnumY / 2, j - gridnumX / 2, offsetH[i][j]);
			}
			fprintf(fpoffsetH, "\n");
		}
		fclose(fpoffsetH);

		// ���Ƶ��ļ���MOLA���ݶ�Ӧ��ƽ�ƺ��DEM��
		FILE* fpcon = fopen(DEMGeoCoordinatefile.c_str(), "w");
		fprintf(fpcon, "#molalon molalat molah demlon demlat demh no no no\n");

		// MOLA�ͳ�ʼDEM�̶߳Ա�ͼ 
		int pointsNum_ccc = 0; // ��ĸ�������֤dtm�ĵ���Ч
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

		plotAltitude.save(picturefile); // ֻ��֧�����·��

		dem.destroy();
	}
}

// �ֶ�ƥ�䣬����CCC�����治һ�µĵ㼯
void saveNone_CCC_points(vector<vector<GeoPoint>> points_IMG_MOLA, string demFile, string fileMOLA_filter)
{
	//����㼯�ļ�
	FILE* fp = fopen(fileMOLA_filter.c_str(), "w");

	//DEMӰ��
	GeoImageIO dem;
	dem.open(demFile);

	//��ȡӰ�񸲸Ƿ�Χ
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = trans[2], ori_Dem_lon = trans[0];
	double end_Dem_lat = trans[2] + trans[3] * demHei, end_Dem_lon = trans[0] + trans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);//���뻺��

	// ��ȡMOLA���ݼ���Ӧ��DEM����
	//FILE* fpThreshold = fopen("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\threashold.txt", "w");
	int part_index = 0;
	for (int index_oribit = 0; index_oribit < points_IMG_MOLA.size(); index_oribit++) // ������
	{
		cout << "�����" << index_oribit << "��" << endl;
		int part_points = 15; // 15�����Ϊһ��
		int points_nums = points_IMG_MOLA[index_oribit].size(); // �˹�ĵ�ĸ���
		int part_nums = points_nums / part_points; // �˹칲�ֳɵĶ���

		for (int index_part = 0; index_part < part_nums; index_part++)
		{
			vector<double> DEMvalue, MOLAvalue;

			for (int i = 0; i < part_points; i++)
			{
				int index_point = index_part * part_points + i; // ��ı��

				//������ӦDEM�ĸ߳�
				double h_value = 0.;
				dem.getBufferValue(points_IMG_MOLA[index_oribit][index_point].lon, points_IMG_MOLA[index_oribit][index_point].lat, 2, &h_value);

				// ��ȡ���Ϊ������Чֵ
				if (h_value)
				{
					DEMvalue.push_back(h_value); // DEM�߳�
					MOLAvalue.push_back(points_IMG_MOLA[index_oribit][index_point].h); // MOLA�߳�
				}

			}

			// ���㷽��
			part_index++;
			double variance = computeVariance(MOLAvalue);

			// ����CCC
			double coefficientNcc = computerNCC(DEMvalue, MOLAvalue); //�������ϵ��(-1��1֮��)

			//fprintf(fpThreshold, "%d\t%lf\t%lf\n", part_index, variance, coefficientNcc);

			cout << "ϵ��=" << coefficientNcc << " �������" << MOLAvalue.size() << endl;

			//if (coefficientNcc < 0.5 && MOLAvalue.size() >= 10 && variance > 1000.0) // �������ϵ��С�������ĵ㼯
			if (MOLAvalue.size() >= 13 && variance > 1000.0) // ���淽���ĵ㼯
			{
				for (int j = 0; j < part_points; j++)
				{
					int index_point = index_part * part_points + j; // ��ı��
					fprintf(fp, "%lf\t%lf\t%lf\n", points_IMG_MOLA[index_oribit][index_point].lon, points_IMG_MOLA[index_oribit][index_point].lat, points_IMG_MOLA[index_oribit][index_point].h);
				}
			}
		}
	}
	fclose(fp);
	//fclose(fpThreshold);
}