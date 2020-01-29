#pragma once
#include <iostream>
#include <atltime.h>
#include <string>
#include <gdal.h>
#include <gdal_priv.h>


using namespace std;

//Zwraca liczbê dni miêdzy obecnym dniem a dniem zdjêcia (z nazwy pliku jak WV.A2015182.0920.tif)
int daysAgo(string file){
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT]; 
	_splitpath_s(file.c_str(), drive, dir, fname, ext); 
	string fileName = fname;
	int days = 0;
	string sy = fileName.substr(4, 4);
	string sd = fileName.substr(8, 3);
	int y = atoi(sy.c_str());
	int d = atoi(sd.c_str());
	CTime now = CTime::GetCurrentTime();
	int monthsDays[] = {31,28,31,30,31,30,31,31,30,31,30,31};
	if(y % 4 == 0){
		monthsDays[1] = 29;
	}
	int m = 1;
	int md = 0;
	int daysMonthCounter = 0;
	while(m <= 12){
		if(daysMonthCounter + monthsDays[m-1] < d){
			daysMonthCounter += monthsDays[m-1];
			m++;
		}
		else{
			md = d - daysMonthCounter;
			break;
		}
	}

	
	CTime imageTime = CTime(y, m, md, 0,0,0);
	CTimeSpan diff = now - imageTime;
	
	return (int)diff.GetDays();
	
	/*char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT]; 
	_splitpath_s(file.c_str(), drive, dir, fname, ext); 
	string fileName = fname;
	string sy = fileName.substr(4, 4);
	string sm = fileName.substr(12,2);
	string sd = fileName.substr(14,2);
	int y = atoi(sy.c_str());
	int m = atoi(sm.c_str());
	int d = atoi(sd.c_str());
	CTime imageTime = CTime(y, m, d, 0, 0, 0);
	CTime now = CTime::GetCurrentTime();
	CTimeSpan diff = now - imageTime;
	return  (int)diff.GetDays();*/
}

void moisture(string plikWilgotnoscNowy,string plikWilgotnoscObecny, string plikWilgotnoscObecnyDaysAgo, string OSGeo4WShellBat){
		float nodata2 = (float)26.31399917602539;
		string cmdMerge = OSGeo4WShellBat + " gdal_merge ";
		string pusty = "pusty.tif";
		string tempNowyM = "tempNowyM.tif";
		string tempStaryM = "tempStaryM.tif";
		string tempDaysAgo = "tempDaysAgo.tif";
		int daysago = daysAgo(plikWilgotnoscNowy);
		cout << "\nDays ago: " << daysago;
		string cmdDelPusty = "del " + pusty;
		string cmdDelTempNowyM = "del " + tempNowyM;
		string cmdDelTempStaryM = "del " + tempStaryM;
		string cmdDelTempDaysAgo = "del " + tempDaysAgo;
		string delWilgotonoscObecny = "del " + plikWilgotnoscObecny;
		string delDaysAgoObecny = "del " + plikWilgotnoscObecnyDaysAgo;
		string cmdCreatePusty = cmdMerge + "-o " + pusty + " -init 0 -createonly "
			+ plikWilgotnoscObecny + " " + plikWilgotnoscNowy;
		string cmdCreateTempNowyM = cmdMerge + "-o " + tempNowyM + " " + pusty + " " + plikWilgotnoscNowy;
		string cmdCreateTempStaryM = cmdMerge + "-o " + tempStaryM + " " + pusty + " " + plikWilgotnoscObecny;
		string cmdCreateTempDaysAgo = cmdMerge + "-ot Byte -init 255 -createonly -o " + tempDaysAgo + " " + pusty;
		string cmdCreateTempDaysAgo2 = cmdMerge + "-o " + tempDaysAgo + " " + tempDaysAgo + " " + plikWilgotnoscObecnyDaysAgo;

		//usuniecie (jeœli istnieje) pustego pliku tymczasowego
		system(cmdDelPusty.c_str());
		//utworzenie pustego zdjecia (wymiary okalajace z³aczenia obecnego wyniku oraz nowego zdjêcia)
		cout << "\n"<<cmdCreatePusty;
		system(cmdCreatePusty.c_str());
		//wpisanie nowego fragmentu do pustego obrazu
		system(cmdDelTempNowyM.c_str());
		system(cmdCreateTempNowyM.c_str());
		//wpisanie obecnego wyniku do pustego obrazu (w przypadku gdy nowy wychodzi poza obszar zmienia siê rozmiar obrazu wynikowego)
		system(cmdDelTempStaryM.c_str());
		system(cmdCreateTempStaryM.c_str());
		//wpisanie obecnego wyniku DaysAgo do pustego obrazu
		system(cmdDelTempDaysAgo.c_str());
		system(cmdCreateTempDaysAgo.c_str());

		//system(delDaysAgoObecny.c_str());
		system(cmdCreateTempDaysAgo2.c_str());
		//
		//system(cmdCreateDaysAgo.c_str());

		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset *dataSetTempStaryM = (GDALDataset*) GDALOpen(tempStaryM.c_str(), GA_ReadOnly);
		GDALDataset *dataSetTempNowyM = (GDALDataset*) GDALOpen(tempNowyM.c_str(), GA_ReadOnly);
		GDALDataset *dataSetDaysAgo = (GDALDataset*) GDALOpen(tempDaysAgo.c_str(), GA_ReadOnly);
		GDALRasterBand *bandStary = dataSetTempStaryM->GetRasterBand(1);
		GDALRasterBand *bandNowy = dataSetTempNowyM->GetRasterBand(1);
		GDALRasterBand *bandDaysAgo = dataSetDaysAgo->GetRasterBand(1);
		int X = bandStary->GetXSize();
		int Y = bandStary->GetYSize();
		if(X != bandNowy->GetXSize() || Y != bandStary->GetYSize()){
			cout << "\n#Error bad sizes.";
			return;
		}
		float *outData = (float *) CPLMalloc(sizeof(float) * X*Y);
		byte *daysData = (byte *) CPLMalloc(sizeof(byte) * X*Y);
		GDALRasterIO(bandDaysAgo, GF_Read, 0, 0, X, Y, daysData, X, Y, GDT_Byte, 0,0);
		GDALClose(dataSetDaysAgo);
		float *lineStary = (float *) CPLMalloc(sizeof(float)*X);
		float *lineNowy = (float *) CPLMalloc(sizeof(float)*X);
		for(int y=0; y<Y; y++){
			GDALRasterIO(bandStary, GF_Read, 0, y, X, 1, lineStary, X, 1, GDT_Float32, 0,0);
			GDALRasterIO(bandNowy, GF_Read, 0, y, X, 1, lineNowy, X, 1, GDT_Float32, 0,0);
			for(int x=0; x<X; x++){
				if(lineNowy[x] != 0 && lineNowy[x] != nodata2){
					outData[y*X + x] = lineNowy[x];
					daysData[y*X + x] = (byte)(daysago);
				}
				else{
					outData[y*X + x] = lineStary[x];
					int days = daysData[y*X + x];
					daysData[y*X + x] = (byte)(days+1);
				}
			}
		}
		system(delWilgotonoscObecny.c_str());
		GDALDataset *dataSetNowyObecny = poDriver->Create( plikWilgotnoscObecny.c_str(), X, Y, 1, GDALDataType::GDT_Float32, NULL );
		double *geoTransform = (double *) CPLMalloc(sizeof(double)*6);
		dataSetTempNowyM->GetGeoTransform(geoTransform);
		dataSetNowyObecny->SetGeoTransform(geoTransform);
		dataSetNowyObecny->SetProjection(dataSetTempNowyM->GetProjectionRef());
		GDALRasterIO(dataSetNowyObecny->GetRasterBand(1), GF_Write, 0,0,X,Y,outData, X,Y,GDT_Float32,0,0);
		dataSetNowyObecny->FlushCache();

		dataSetDaysAgo = poDriver->Create( plikWilgotnoscObecnyDaysAgo.c_str(), X, Y, 1, GDALDataType::GDT_Byte, NULL );
		dataSetDaysAgo->SetGeoTransform(geoTransform);
		dataSetDaysAgo->SetProjection(dataSetTempNowyM->GetProjectionRef());
		GDALRasterIO(dataSetDaysAgo->GetRasterBand(1), GF_Write, 0,0,X,Y,daysData, X,Y,GDT_Byte,0,0);
		dataSetDaysAgo->FlushCache();
		GDALClose(dataSetDaysAgo);
		GDALClose(dataSetTempStaryM);
		GDALClose(dataSetTempNowyM);
		CPLFree(outData);
		CPLFree(lineStary);
		CPLFree(lineNowy);
}