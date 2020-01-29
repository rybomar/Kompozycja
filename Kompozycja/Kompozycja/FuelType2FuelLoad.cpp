//#include "Images.h"
//using namespace std;
//
//int main(int argc, char * argv[]){
//	
//	if(argc != 3){
//		cout << "\npathToFuelType.tif resultPathFuelLoad.tif";
//		return 0;
//	}
//
//	GDALAllRegister();
//	Images ims = Images();
//	
//	ims.addBand(argv[1]);
//	ims.convertFuelType2FuelLoad(argv[2]);
//	cout << "\nDone";
//	return 0;
//}