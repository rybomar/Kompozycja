#include <gdal.h>
#include <gdal_priv.h>
#include <iostream>
#include <vector>
#include <fstream>
#include "Tools.h"

using namespace std;

double Log( double n )  
{  
    return log( n ) / log( 2 );  
}

class Images{
	vector<GDALDataset*> dataSets;//przechowuje "uchwyty" do plików, kolejnoœæ u³o¿enia to kolejnoœæ kana³ów w wyniku
	vector<int> bandNumbers;//przechowuje numery kana³ów z poszczególnych plików
	vector<GDALDataType> dataTypes;//przechowuje typy kana³ów, np. GDT_Byte, GDT_UInt16
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
	bool addBand(string fileName, int bandNumber=1){
		if(ifstream(fileName)) {
			GDALDataset *dataSet = (GDALDataset*) GDALOpen(fileName.c_str(), GA_ReadOnly);
			if (dataSet->GetRasterCount() < bandNumber || bandNumber < 0) {
				cout << "\n#dataSet->GetRasterCount()[" << dataSet->GetRasterCount() << "] < bandNumber[" << bandNumber << "] ";
				return false;
			}
			if(dataSets.size() == 0){
				projection = dataSet->GetProjectionRef();
				if( dataSet->GetGeoTransform( geoTransform ) == CE_None )
				{
					// Data structure load ok. Do nothing.
				} else {
					//return false;
				}
				cols = dataSet->GetRasterXSize();
				rows = dataSet->GetRasterYSize();
			}
			dataSets.push_back(dataSet);
			bandNumbers.push_back(bandNumber);
			dataTypes.push_back(dataSet->GetRasterBand(bandNumber)->GetRasterDataType());
			return true;
		}
		cout << "\n#Wrong fileName.";
		return false;
	}
	void clearBands(){
		dataSets.clear();
		bandNumbers.clear();
	}
	bool save(string fileName, int cutL, int cutR, int cutU, int cutD, GDALDataType type= GDT_Unknown){
		type = checkOutputDataType(type);
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
		//void *pabyData = (void *) CPLMalloc(sizeof(float) * (cols-cutL-cutR) * (rows-cutD-cutU));
		int floatSize = sizeof(float);
		unsigned int allocSize = floatSize*cols*rows;
		
		float *pabyData = (float*)VSIMalloc(allocSize);
		float *outData = (float *)VSIMalloc(sizeof(float) * (cols-cutL-cutR) * (rows-cutD-cutU));

		for(int i=0; i<dataSets.size(); i++){
			GDALRasterBand *band = dataSets[i]->GetRasterBand(bandNumbers[i]);
			//for(int l=cutU; l<rows-cutU; l+=wiersze){
				//GDALRasterIO(band, GF_Read, cutL, l+(l-cutU)*wiersze, cols-cutL, rows-cutU, pabyData, cols-cutL-cutR, rows-cutD-cutU, GDT_Float32, 0,0);
				GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);
				for(int r =cutU; r<rows-cutD; r++){
					for(int c=cutL; c<cols-cutR; c++){
						int rOut = r-cutU;
						int cOut = c-cutU;
						double v = pabyData[r*cols + c];
						outData[rOut*(cols-cutL-cutR) + cOut] = pabyData[r*cols + c];
						//float o= (float*)outData[rOut*(cols-cutL-cutR) + cOut];
						// = pabyData[r*cols + c];
					}
				}
				GDALRasterIO(poDataset->GetRasterBand(i+1), GF_Write, 0, 0, cols-cutL-cutR, rows-cutD-cutU, outData, cols-cutL-cutR, rows-cutD-cutU, GDT_Float32, 0,0);	
			//}
			cout << "\nDodano kanal " << i+1 << "/"<<dataSets.size();
		}
		poDataset->FlushCache();
		return true;
	}
	bool save(string fileName, GDALDataType type=GDT_Unknown){
		return save(fileName, 0,0,0,0, type);
	}
	GDALDataType checkOutputDataType(GDALDataType declaredDataType) {
		GDALDataType type = declaredDataType;
		if (type == GDT_Unknown) {
			for (int i = 0; i < dataTypes.size(); i++) {
				if (dataTypes[i] > type) {
					type = dataTypes[i];
				}
			}
		}
		return type;
	}
	bool saveSplit(string fileName, string resultDirName, int cutL, int cutR, int cutU, int cutD, GDALDataType type=GDT_Float32){
		if(ifstream(fileName)) {
			_mkdir(resultDirName.c_str());
			GDALDataset *dataSet = (GDALDataset*) GDALOpen(fileName.c_str(), GA_ReadOnly);
			for(int i=0; i<dataSet->GetRasterCount(); i++){
				string fileNumberName = std::to_string(i+1) + ".tif";
				cout << endl << fileNumberName + "... ";
				clearBands();
				addBand(fileName, i+1);
				string oneBandFileName = resultDirName + "/" + fileNumberName;
				save(oneBandFileName, cutL, cutR, cutU, cutD, type);
				cout << "saved    ";
			}
			return true;
		}
		return false;
	}
	void getMinMax(string type, string srcDir, float *hmin, float *hmax, int bandNumber=1){
		vector<int> R = vector<int>();
		vector<int> I = vector<int>();
		//string srcDir = "J:/MT_SAR_wyniki/";//"E:/MT_SAR/wyniki/vh_p/";
		//string dstDir = "H:/MT_SAR/wynikiLog/";
		int maxv = 13;
		int v = 1;
		for(int i=1; i<=maxv; i++){
			for(int j=1; j<=maxv; j++){
				if(i==j)
					R.push_back(v++);
				else
				if(j>i){
					R.push_back(v++);
					I.push_back(v++);
				}
			}
		}
		if(type == "1"){
			//addBand(srcDir + std::to_string(bandNumber) + ".tif");
		}
		if(type == "R"){//rzeczywiste
			for(int v=0; v<R.size(); v++){
				addBand(srcDir + std::to_string(R[v]) + ".tif");
			}
		}
		if(type == "I"){//urojone
			for(int v=0; v<I.size(); v++){
				addBand(srcDir + std::to_string(I[v]) + ".tif");
			}
		}
		if(type == "E"){//tylko entropia
			addBand(srcDir + std::to_string(1) + ".tif");
		}
		if(type == "A"){//tylko alfa
			addBand(srcDir + std::to_string(2) + ".tif");
		}
		if(type == "V"){//wartoœci w³asne
			for(int v=3; v<=2+maxv; v++){
				addBand(srcDir + std::to_string(v) + ".tif");
			}
		}
		float *pabyData = (float*)CPLMalloc(sizeof(float)*cols*rows);
		float min = 99999999;
		float max = -min;
		int bandNumberMin = -1;
		int bandNumberMax = -1;
		int colMin = -1;
		int colMax = -1;
		int rowMin = -1;
		int rowMax = -1;

		if(type == "1"){
			GDALRasterBand *band = dataSets[bandNumber]->GetRasterBand(bandNumbers[bandNumber]);
			GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);
			for(int r=0; r<rows; r++){
				for(int c=0; c<cols; c++){
					float value = pabyData[r*cols + c];
					if(value < min){
						bandNumberMin = bandNumber;
						colMin = c;
						rowMin = r;
						min = value;
					}
					if(value > max){
						bandNumberMax = bandNumber;
						colMax = c;
						rowMax = r;
						max = value;
					}
				}
			}
		}
		else{
			for(int i=0; i<dataSets.size(); i++){
				cout << "  " << i;
				GDALRasterBand *band = dataSets[i]->GetRasterBand(bandNumbers[i]);
				GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);
				for(int r=0; r<rows; r++){
					for(int c=0; c<cols; c++){
						float value = pabyData[r*cols + c];
						if(value < min){
							bandNumberMin = i+1;
							colMin = c;
							rowMin = r;
							min = value;
						}
						if(value > max){
							bandNumberMax = i+1;
							colMax = c;
							rowMax = r;
							max = value;
						}
					}
				}
			}
		}
		CPLFree(pabyData);
		if(type != "1"){
			clearBands();
		}
		cout << "\n\n\t" << type << ":";
		cout << "\nMIN: " << min<< "         band: " << bandNumberMin << "  row: " << rowMin << "   col: " << colMin;
		cout << "\nMAX: " << max<< "         band: " << bandNumberMax << "  row: " << rowMax << "   col: " << colMax;
		hmin[0] = min;
		hmax[0] = max;
		//getchar();
	}
	void logEA(string srcDir, string destDir, float Vmin, float Vmax){
		cout << "\n\tVmin: " << Vmin << "           Vmax: " << Vmax;
		for(int v=1; v<=15; v++){
				addBand(srcDir + std::to_string(v) + ".tif");
		}
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		double przesE = 0;
		double przesA = 0;
		double przesV = Log(-Vmin + 1);
		double mnoznikE = 254;
		double mnoznikA = 2.83;
		double mnoznikV = (double)(255)/(Log(Vmax + 1) + przesV);
		float *pabyData = (float*)CPLMalloc(sizeof(float)*cols*rows);
		GByte *resData = (GByte*)CPLMalloc(sizeof(GByte)*cols*rows);
		for(int i=0; i<dataSets.size(); i++){
			GDALRasterBand *band = dataSets[i]->GetRasterBand(bandNumbers[i]);
			GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);
			string fileName = destDir + std::to_string(i+1) + ".tif";
			GDALDataset *poDataset = poDriver->Create( fileName.c_str(), cols, rows, 1, GDALDataType::GDT_Byte, NULL );
			poDataset->SetGeoTransform(geoTransform);
			poDataset->SetProjection(projection.c_str());
			double przes = przesV;
			double mnoznik = mnoznikV;
			if(i==0){
				mnoznik = mnoznikE;
				przes = przesE;
			}
			if(i==1){
				mnoznik = mnoznikA;
				przes = przesA;
			}
			if(i>1){
				mnoznik = mnoznikV;
				przes = przesV;
			}
			for(int r=0; r<rows; r++){
				for(int c=0; c<cols; c++){
					double v = pabyData[r*cols + c];
					int res = 0;
					if(v == 0){
						//resData[r*cols + c] = 0;
						res = 0;
					}
					else{
						if(i > 1){
							if(v<0){
								res = (int)((-Log(-v+1) + przes) * mnoznik);
							}
							else{
								res = (int)((Log(v+1) + przes) * mnoznik);
							}
						}
						else{
							res = (int)((v+przes) * mnoznik)+1;
						}
					}
					if(res > 255) res = 255;
					if(res < 0) res = 0;
					resData[r*cols + c] = res;
				}
			}
			GDALRasterIO(poDataset->GetRasterBand(1), GF_Write, 0, 0, cols, rows, resData, cols, rows, GDALDataType::GDT_Byte, 0,0);	
			poDataset->FlushCache();
		}
		CPLFree(pabyData);
		CPLFree(resData);
		clearBands();
	}
	void logEA(string srcDir, string destDir){
		vector<int> V = vector<int>();
		for(int i=1; i<=15; i++){
			V.push_back(i);
		}
		for(int v=0; v<V.size(); v++){
			addBand(srcDir + std::to_string(V[v]) + ".tif");
		}
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		float *pabyData = (float*)CPLMalloc(sizeof(float)*cols*rows);
		GByte *resData = (GByte*)CPLMalloc(sizeof(GByte)*cols*rows);
		for(int i=0; i<dataSets.size(); i++){
			GDALRasterBand *band = dataSets[i]->GetRasterBand(bandNumbers[i]);
			GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);
			string fileName = destDir + std::to_string(V[i]) + ".tif";
			GDALDataset *poDataset = poDriver->Create( fileName.c_str(), cols, rows, 1, GDALDataType::GDT_Byte, NULL );
			poDataset->SetGeoTransform(geoTransform);
			poDataset->SetProjection(projection.c_str());
			float *minA = new float[1];
			float *maxA = new float[1];
			this->getMinMax("1", srcDir, minA, maxA, i);
			float min = minA[0];
			float max = maxA[0];
			//double przes = Log(-min + 1);
			double przes = -min;
			//double mnoznik = (double)(255)/(Log(max + 1) + przes);
			double mnoznik = (double)(255)/(max+przes);
			if(i == 0){
				przes = 0;
				mnoznik = 254;
			}
			if(i==1){
				przes = 0;
				mnoznik = 2.83;
			}

			for(int r=0; r<rows; r++){
				for(int c=0; c<cols; c++){
					double v = pabyData[r*cols + c];
					GByte wynik=0;
					//int wynik = 0;
					if(i==0 || i==1){
						wynik = (int)(v * mnoznik);
					}
					else{
						if(v < 0){
							//wynik = (int)((-Log(-v+1) + przes) * mnoznik);
							wynik = (int)((v+przes) * mnoznik);
						}
						else{
							if(v >0){
								int a=2;
							}
							//wynik = (int)((Log(v+1) + przes) * mnoznik);
							wynik = (int)((v+przes) * mnoznik);
						}
					}
					resData[r*cols + c] = wynik+1;
					if(v == 0){
						resData[r*cols + c] = 0;
					}
				}
			}
			GDALRasterIO(poDataset->GetRasterBand(1), GF_Write, 0, 0, cols, rows, resData, cols, rows, GDALDataType::GDT_Byte, 0,0);	
			poDataset->FlushCache();
		}
		CPLFree(pabyData);
		clearBands();
	}
	void LogKoh(string srcDir, string destDir, double Rmax, double Rmin, double Imax, double Imin){
		double przesR = Log(-Rmin + 1);
		double przesI = Log(-Rmin + 1);
		double mnoznikR = (double)(255)/(Log(Rmax + 1) + przesR);
		double mnoznikI = (double)(255)/(Log(Imax + 1) + przesI);
		vector<int> R = vector<int>();
		vector<int> I = vector<int>();
		int v = 1;
		for(int i=1; i<=13; i++){
			for(int j=1; j<=13; j++){
				if(i==j)
					R.push_back(v++);
				else
				if(j>i){
					R.push_back(v++);
					I.push_back(v++);
				}
			}
		}
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		string type = "R";
		double mnoznik = mnoznikR;
		double przes = przesR;
		if(type == "R"){
			for(int v=0; v<R.size(); v++){
				addBand(srcDir + std::to_string(R[v]) + ".tif");
			}
			float *pabyData = (float*)CPLMalloc(sizeof(float)*cols*rows);
			GByte *resData = (GByte*)CPLMalloc(sizeof(GByte)*cols*rows);
			for(int i=0; i<dataSets.size(); i++){
				GDALRasterBand *band = dataSets[i]->GetRasterBand(bandNumbers[i]);
				GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);
				string fileName = destDir + std::to_string(R[i]) + ".tif";
				GDALDataset *poDataset = poDriver->Create( fileName.c_str(), cols, rows, 1, GDALDataType::GDT_Byte, NULL );
				poDataset->SetGeoTransform(geoTransform);
				poDataset->SetProjection(projection.c_str());
				for(int r=0; r<rows; r++){
					for(int c=0; c<cols; c++){
						double v = pabyData[r*cols + c];
						GByte wynik=0;
						//int wynik = 0;
						if(v < 0){
							wynik = (int)((-Log(-v+1) + przes) * mnoznik);
						}
						else{
							wynik = (int)((Log(v+1) + przes) * mnoznik);
						}
						resData[r*cols + c] = wynik+1;
						if(v == 0){
							resData[r*cols + c] = 0;
						}
					}
				}
				GDALRasterIO(poDataset->GetRasterBand(1), GF_Write, 0, 0, cols, rows, resData, cols, rows, GDALDataType::GDT_Byte, 0,0);	
				poDataset->FlushCache();
			}
			CPLFree(pabyData);
			CPLFree(resData);
		}

		clearBands();
		type = "I";
		if(type == "I"){
			mnoznik = mnoznikI;
			przes = przesI;
			for(int v=0; v<I.size(); v++){
				addBand(srcDir + std::to_string(I[v]) + ".tif");
			}
			float *pabyData = (float*)CPLMalloc(sizeof(float)*cols*rows);
			GByte *resData = (GByte*)CPLMalloc(sizeof(GByte)*cols*rows);
			for(int i=0; i<dataSets.size(); i++){
				GDALRasterBand *band = dataSets[i]->GetRasterBand(bandNumbers[i]);
				GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);
				string fileName = destDir + std::to_string(I[i]) + ".tif";
				GDALDataset *poDataset = poDriver->Create( fileName.c_str(), cols, rows, 1, GDALDataType::GDT_Byte, NULL );
				poDataset->SetGeoTransform(geoTransform);
				poDataset->SetProjection(projection.c_str());
				for(int r=0; r<rows; r++){
					for(int c=0; c<cols; c++){
						double v = pabyData[r*cols + c];
						GByte wynik = 0;
						if(v < 0){
							wynik = (int)((-Log(-v+1) + przes) * mnoznik);
						}
						else{
							wynik = (int)((Log(v+1) + przes) * mnoznik);
						}
						resData[r*cols + c] = wynik;
					}
				}
				GDALRasterIO(poDataset->GetRasterBand(1), GF_Write, 0, 0, cols, rows, resData, cols, rows, GDALDataType::GDT_Byte, 0,0);	
				poDataset->FlushCache();
			}
			CPLFree(pabyData);
			CPLFree(resData);
		}
		

	}

	//Logarytm wed³ug wzorów:
	//TR = 1/2 * ln (R^2 + I^2)
	//TI = atan(I/R);
	//Na razie dzia³a tylko dla macierzy 2x2
	static void LogKohSpecified(string inputFile, string outputFile) {
		if (ifstream(inputFile)) {
			GDALDataset *dataIn = (GDALDataset*)GDALOpen(inputFile.c_str(), GA_ReadOnly);
			GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
			int X = dataIn->GetRasterXSize();
			int Y = dataIn->GetRasterYSize();
			GDALDataset *dataOut = poDriver->Create(outputFile.c_str(), X, Y, dataIn->GetRasterCount(), dataIn->GetRasterBand(1)->GetRasterDataType(), NULL);
			double* geoTransform = (double *)CPLMalloc(sizeof(double) * 6);
			dataIn->GetGeoTransform(geoTransform);
			dataOut->SetGeoTransform(geoTransform);
			dataOut->SetProjection(dataIn->GetProjectionRef());
			int n = sqrt(dataIn->GetRasterCount());
			if (n != 2) {
				throw(new string("Only 4 bands supported at now."));
			}
			//tworzenie listy z warstwami rzeczywistymi (R) oraz urojonymi (I)
			float *inputBand1 = (float*)new float*[(sizeof(float*) *X)];
			float *inputBand2 = (float*)new float*[(sizeof(float*) *X)];
			float *inputBand3 = (float*)new float*[(sizeof(float*) *X)];
			float *inputBand4 = (float*)new float*[(sizeof(float*) *X)];
			for (int y = 0; y < Y; y++) {
				GDALRasterIO(dataIn->GetRasterBand(1), GF_Read, 0, y, X, 1, inputBand1, X, 1, GDT_Float32, 0, 0);
				GDALRasterIO(dataIn->GetRasterBand(2), GF_Read, 0, y, X, 1, inputBand2, X, 1, GDT_Float32, 0, 0);
				GDALRasterIO(dataIn->GetRasterBand(3), GF_Read, 0, y, X, 1, inputBand3, X, 1, GDT_Float32, 0, 0);
				GDALRasterIO(dataIn->GetRasterBand(4), GF_Read, 0, y, X, 1, inputBand4, X, 1, GDT_Float32, 0, 0);
				for (int x = 0; x < X; x++) {
					float r11 = inputBand1[x];
					float i11 = 0;
					float r12 = inputBand2[x];
					float i12 = inputBand3[x];
					float r21 = inputBand4[x];
					float i21 = 0;
					inputBand1[x] = 0.5 * log(r11 * r11 + i11 * i11);
					inputBand2[x] = 0.5 * log(r11 * r12 + i12 * i12);
					inputBand3[x] = atan(i12 / r12);
					inputBand4[x] = 0.5 * log(r21 * r21 + i21 * i21);
				}
				GDALRasterIO(dataOut->GetRasterBand(1), GF_Write, 0, y, X, 1, inputBand1, X, 1, GDT_Float32, 0, 0);
				GDALRasterIO(dataOut->GetRasterBand(2), GF_Write, 0, y, X, 1, inputBand2, X, 1, GDT_Float32, 0, 0);
				GDALRasterIO(dataOut->GetRasterBand(3), GF_Write, 0, y, X, 1, inputBand3, X, 1, GDT_Float32, 0, 0);
				GDALRasterIO(dataOut->GetRasterBand(4), GF_Write, 0, y, X, 1, inputBand4, X, 1, GDT_Float32, 0, 0);
			}
			dataOut->FlushCache();
		}
		else
		{
			throw(new string("Error. Input file not exists."));
		}
	}
	static void normKoh(string inputFile, string outputFile){
		if(ifstream(inputFile)) {
			GDALDataset *dataIn = (GDALDataset*) GDALOpen(inputFile.c_str(), GA_ReadOnly);
			GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
			int X = dataIn->GetRasterXSize();
			int Y = dataIn->GetRasterYSize();
			GDALDataset *dataOut = poDriver->Create( outputFile.c_str(), X, Y, dataIn->GetRasterCount(), dataIn->GetRasterBand(1)->GetRasterDataType(), NULL );
			double* geoTransform = (double *) CPLMalloc(sizeof(double)*6);
			dataIn->GetGeoTransform(geoTransform);
			dataOut->SetGeoTransform(geoTransform);
			dataOut->SetProjection(dataIn->GetProjectionRef());
			//tworzenie listy z warstwami rzeczywistymi (R) oraz urojonymi (I)
			int n = sqrt(dataIn->GetRasterCount());
			vector<int> *R = new vector<int>();
			vector<int> *I = new vector<int>();
			vector<int> *All = new vector<int>();
			int v = 1;
			int vv = 1;
			for(int i=1; i<=n; i++){
				for(int j=1; j<=n; j++){
					All->push_back(vv++);
					if(i==j)
						R->push_back(v++);
					else
					if(j>i){
						R->push_back(v++);
						I->push_back(v++);
					}
				}
			}
			int linesInBufor = 10;
			float *inputBand = (float*) CPLMalloc(sizeof(float*) * Y *X);
			float *res = (float*) CPLMalloc(sizeof(float*) * Y *X);
			for(int type = 2; type <= 2; type++){//type: R lub I
				int nT = 0;
				vector<int> *V = NULL;
				if(type == 0){ nT = R->size(); V = R;}
				if(type == 1){ nT = I->size(); V = I;}
				if(type == 2){ nT = All->size(); V = All;}
				double **Sum = (double**) CPLMalloc(sizeof(double*) * Y);
				for (int y = 0; y < Y; y++) {
					Sum[y] = (double*) CPLMalloc(sizeof(double) * X);
					for(int x=0; x<X; x++){
						Sum[y][x] = 0;
					}
				}
				cout <<"\n\nSuma kwadratow:";
				for(int b=0; b<V->size(); b++){//po wszystkich warstwach
					GDALRasterIO(dataIn->GetRasterBand(V->at(b)), GF_Read, 0, 0, X, Y, inputBand, X, Y, GDT_Float32, 0,0);
					for (int y = 0; y < Y; y++) {
						for (int x = 0; x < X; x++) {
							double inp = inputBand[y * X + x];
							if(y==x && y==0){
								int a = 3;
							}
							Sum[y][x] += inp * inp;
						}
					}
					cout << "\n" << b+1 << " / " << V->size();
				}
				cout <<"\n\nReszta wzorow i zapis:";
				for(int b=0; b<V->size(); b++){//po wszystkich warstwach
					GDALRasterIO(dataIn->GetRasterBand(V->at(b)), GF_Read, 0, 0, X, Y, inputBand, X, Y, GDT_Float32, 0,0);
					for (int y = 0; y < Y; y++) {
						for (int x = 0; x < X; x++) {
							double a = inputBand[y * X + x];
							double b = Sum[y][x];
							double r = a / sqrt(b / V->size());
							res[y * X + x] = inputBand[y * X + x] / sqrt(Sum[y][x]/V->size());
						}
					}
					cout << "\n" << b+1 << " / " << V->size() << "   policzono";
					GDALRasterIO(dataOut->GetRasterBand(V->at(b)), GF_Write, 0, 0, X, Y, res, X, Y, GDT_Float32, 0,0);
					cout << "   zapisano";
				}
				dataOut->FlushCache();
			}
		}
		else
			cout << "\nError. Input file not exists.";
	}
	static void normKohComplex(string inputFile, string outputFile){
		if(ifstream(inputFile)) {
			GDALDataset *dataIn = (GDALDataset*) GDALOpen(inputFile.c_str(), GA_ReadOnly);
			GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
			int X = dataIn->GetRasterXSize();
			int Y = dataIn->GetRasterYSize();
			GDALDataset *dataOut = poDriver->Create( outputFile.c_str(), X, Y, dataIn->GetRasterCount(), dataIn->GetRasterBand(1)->GetRasterDataType(), NULL );
			double* geoTransform = (double *) CPLMalloc(sizeof(double)*6);
			dataIn->GetGeoTransform(geoTransform);
			dataOut->SetGeoTransform(geoTransform);
			dataOut->SetProjection(dataIn->GetProjectionRef());
			//tworzenie listy z warstwami rzeczywistymi (R) oraz urojonymi (I)
			int n = sqrt(dataIn->GetRasterCount());
			vector<int> *R = new vector<int>();
			vector<int> *I = new vector<int>();
			vector<int> *All = new vector<int>();
			int v = 1;
			int vv = 1;
			for(int i=1; i<=n; i++){
				for(int j=1; j<=n; j++){
					All->push_back(vv++);
					if(i==j)
						R->push_back(v++);
					else
					if(j>i){
						R->push_back(v++);
						I->push_back(v++);
					}
				}
			}
			int linesInBufor = 10;
			float *inputBand = (float*) CPLMalloc(sizeof(float*) * Y *X);
			float *res = (float*) CPLMalloc(sizeof(float*) * Y *X);
			for(int type = 2; type <= 2; type++){//type: R lub I
				int nT = 0;
				vector<int> *V = NULL;
				if(type == 0){ nT = R->size(); V = R;}
				if(type == 1){ nT = I->size(); V = I;}
				if(type == 2){ nT = All->size(); V = All;}
				double **Sum = (double**) CPLMalloc(sizeof(double*) * Y);
				for (int y = 0; y < Y; y++) {
					Sum[y] = (double*) CPLMalloc(sizeof(double) * X);
					for(int x=0; x<X; x++){
						Sum[y][x] = 0;
					}
				}
				cout <<"\n\nSuma kwadratow:";
				for(int b=0; b<V->size(); b++){//po wszystkich warstwach
					GDALRasterIO(dataIn->GetRasterBand(V->at(b)), GF_Read, 0, 0, X, Y, inputBand, X, Y, GDT_Float32, 0,0);
					for (int y = 0; y < Y; y++) {
						for (int x = 0; x < X; x++) {
							Sum[y][x] += inputBand[y * X + x] * inputBand[y * X + x];
						}
					}
					cout << "\n" << b+1 << " / " << V->size();
				}
				cout <<"\n\nReszta wzorow i zapis:";
				for(int b=0; b<V->size(); b++){//po wszystkich warstwach
					GDALRasterIO(dataIn->GetRasterBand(V->at(b)), GF_Read, 0, 0, X, Y, inputBand, X, Y, GDT_Float32, 0,0);
					for (int y = 0; y < Y; y++) {
						for (int x = 0; x < X; x++) {
							double a = inputBand[y * X + x];
							double b = Sum[y][x];
							double r = a / sqrt(b / V->size());
							res[y * X + x] = inputBand[y * X + x] / sqrt(Sum[y][x]/(double)V->size());
						}
					}
					cout << "\n" << b+1 << " / " << V->size() << "   policzono";
					GDALRasterIO(dataOut->GetRasterBand(V->at(b)), GF_Write, 0, 0, X, Y, res, X, Y, GDT_Float32, 0,0);
					cout << "   zapisano";
				}
				dataOut->FlushCache();
			}
		}
		else
			cout << "\nError. Input file not exists.";
	}
	void usrednijLinie(string fileLineDefH, string fileLineDefV, string resFileName, GDALDataType type=GDT_Float32){
		FILE* frH = fopen (fileLineDefH.c_str(), "rt");
		FILE* frV = fopen (fileLineDefV.c_str(), "rt");
		vector<int> X1H = vector<int>();
		vector<int> X2H = vector<int>();
		vector<int> Y1H = vector<int>();
		vector<int> Y2H = vector<int>();
		vector<int> X1V = vector<int>();
		vector<int> X2V = vector<int>();
		vector<int> Y1V = vector<int>();
		vector<int> Y2V = vector<int>();
		
		char cline[80];
		while(fgets(cline, 80, frH) != NULL){
			int x1,x2,y1,y2;
			sscanf (cline, "%d %d %d %d", &x1, &y1, &x2, &y2);
			X1H.push_back(x1);
			X2H.push_back(x2);
			Y1H.push_back(y1);
			Y2H.push_back(y2);
		}
		char clineV[80];
		while(fgets(clineV, 80, frV) != NULL){
			int x1,x2,y1,y2;
			sscanf (clineV, "%d %d %d %d", &x1, &y1, &x2, &y2);
			X1V.push_back(x1);
			X2V.push_back(x2);
			Y1V.push_back(y1);
			Y2V.push_back(y2);
		}
		fclose(frV);
		fclose(frH);
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset *poDataset = poDriver->Create( resFileName.c_str(), cols, rows, dataSets.size(), type, NULL );
		poDataset->SetGeoTransform(geoTransform);
		poDataset->SetProjection(projection.c_str());
		int wiersze = rows;
		//void *pabyData = (void *) CPLMalloc(sizeof(float) * (cols) * (rows));
		float *pabyData = (float*)CPLMalloc(sizeof(float)*cols*rows);
		for(int i=0; i<dataSets.size(); i++){
			GDALRasterBand *band = dataSets[i]->GetRasterBand(bandNumbers[i]);
			GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);
			for(int v=0; v<X1H.size(); v++){
				int x1= X1H[v];	int x2= X2H[v];	int y1= Y1H[v];	int y2= Y2H[v];
				int x = x1;
				int y= y1;
				bool ok=true;//czy dalej przegl¹daæ punkty
				while(ok){
					pabyData[y*cols + x] = (pabyData[(y-1)*cols + x] + pabyData[(y+1)*cols + x])/2;
					x++;
					if(x>x2)
						ok = false;
						
				}
			}
			for(int v=0; v<X1V.size(); v++){
				int x1= X1V[v];	int x2= X2V[v];	int y1= Y1V[v];	int y2= Y2V[v];
				int x = x1;
				int y= y1;
				bool ok=true;//czy dalej przegl¹daæ punkty
				while(ok){
						pabyData[y*cols + x] = (pabyData[y*(cols) + x-1] + pabyData[y*(cols) + x+1])/2;
						y++;
						if(y>y2)
							ok = false;
				}
			}
			GDALRasterIO(poDataset->GetRasterBand(i+1), GF_Write, 0, 0, cols, rows, pabyData, cols, rows, GDT_Float32, 0,0);	
		}
		poDataset->FlushCache();
	}
	void convertFuelType2FuelLoad(string fuelLoadOutputFile){
		__int8 *pabyData = (__int8*)CPLMalloc(sizeof(__int8)*cols*rows);
		float *outData = (float*)CPLMalloc(sizeof(float)*cols*rows);
		GDALRasterBand *band = dataSets[0]->GetRasterBand(bandNumbers[0]);
		GDALRasterIO(band, GF_Read, 0, 0, cols, rows, pabyData, cols, rows, GDT_Byte, 0,0);
		for(int r =0; r<rows; r++){
			for(int c=0; c<cols; c++){
				__int8 v = pabyData[r*cols + c];
				double res = 0;
				switch (v)
				{
				case 1:
					res = 0.74;
					break;
				case 2:
					res = 2.0;
					break;
				case 3:
					res = 3.01;
					break;
				case 4:
					res = 5.01;
					break;
				case 5:
					res = 1.0;
					break;
				case 6:
					res = 1.5;
					break;
				case 7:
					res = 1.13;
					break;
				case 8:
					res = 1.5;
					break;
				case 9:
					res = 2.92;
					break;
				case 10:
					res = 3.01;
					break;
				case 11:
					res = 1.5;
					break;
				case 12:
					res = 4.01;
					break;
				case 13:
					res = 7.01;
					break;
				default:
					break;
				}
				outData[r*cols + c] = (float)res;
			}
		}
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		GDALDataset *poDataset = poDriver->Create( fuelLoadOutputFile.c_str(), cols, rows, dataSets.size(), GDT_Float32, NULL );
		poDataset->SetGeoTransform(geoTransform);
		poDataset->SetProjection(projection.c_str());
		GDALRasterIO(poDataset->GetRasterBand(1), GF_Write, 0, 0, cols, rows, outData, cols, rows, GDT_Float32, 0,0);
		poDataset->FlushCache();
	}
};

