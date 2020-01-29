#include <gdal.h>
#include <gdal_priv.h>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

class Images{
	vector<GDALDataset*> dataSets;//przechowuje "uchwyty" do plików, kolejnoœæ u³o¿enia to kolejnoœæ kana³ów w wyniku
	vector<int> bandNumbers;//przechowuje numery kana³ów z poszczególnych plików
	int cols, rows;
	double* geoTransform;
	string projection;
public:
	/**
	 * Tworzy pusty zbiór obrazów
	 */
	Images(){
		dataSets = vector<GDALDataset*>();
		bandNumbers = vector<int>();
		geoTransform =  (double *) CPLMalloc(sizeof(double)*6);
		projection = "";
	}
	bool addBand(string fileName, int bandNumber){
		if(ifstream(fileName)) {
			GDALDataset *dataSet = (GDALDataset*) GDALOpen(fileName.c_str(), GA_ReadOnly);
			if(dataSet->GetRasterCount() < bandNumber || bandNumber < 0)
				return false;
			if(dataSets.size() == 0){
				projection = dataSet->GetProjectionRef();
				if( dataSet->GetGeoTransform( geoTransform ) == CE_None )
				{
					// Data structure load ok. Do nothing.
				} else {
					return false;
				}
				cols = dataSet->GetRasterXSize();
				rows = dataSet->GetRasterYSize();
			}
			dataSets.push_back(dataSet);
			bandNumbers.push_back(bandNumber);
			return true;
		}
		return false;
	}
	bool save(string fileName, int cutL, int cutR, int cutU, int cutD, GDALDataType type=GDT_Float32){
		if(cutL < 0) cutL = 0;
		if(cutR < 0) cutR = 0;
		if(cutD < 0) cutD = 0;
		if(cutU < 0) cutU = 0;
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset *poDataset = poDriver->Create( fileName.c_str(), cols-cutL-cutR, rows-cutD-cutU, dataSets.size(), type, NULL );
		geoTransform[0] += cutL*geoTransform[1];
		geoTransform[3] += cutU*geoTransform[5];
		poDataset->SetGeoTransform(geoTransform);
		poDataset->SetProjection(projection.c_str());
		int wiersze = rows-cutD-cutU;
		void *pabyData = (void *) CPLMalloc(sizeof(float) * (cols-cutL-cutR) * (rows-cutD-cutU));
		for(int i=0; i<dataSets.size(); i++){
			GDALRasterBand *band = dataSets[i]->GetRasterBand(bandNumbers[i]);
			for(int l=cutU; l<rows-cutU; l+=wiersze){
				GDALRasterIO(band, GF_Read, cutL, l+(l-cutU)*wiersze, cols-cutL, rows-cutU, pabyData, cols-cutL-cutR, rows-cutD-cutU, GDT_Float32, 0,0);
				GDALRasterIO(poDataset->GetRasterBand(i+1), GF_Write, 0, (l-cutU)*wiersze, cols-cutL-cutR, rows-cutD-cutU, pabyData, cols-cutL-cutR, rows-cutD-cutU, GDT_Float32, 0,0);	
			}
			cout << "\nDodano kanal " << i+1 << "/"<<dataSets.size();
		}
		poDataset->FlushCache();
		return true;
	}
	bool save(string fileName, GDALDataType type=GDT_Float32){
		return save(fileName, 0,0,0,0);
	}
};

