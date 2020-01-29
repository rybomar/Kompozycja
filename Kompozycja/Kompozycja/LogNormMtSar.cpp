//#include "Images.h"
//#include <sstream>
//
//int main(int argc, char * argv[]) {
//	GDALAllRegister();
//	Images ims = Images();
//
//	//MIN MAX
//	float *Rmin = new float[1];
//	float *Rmax = new float[1];
//	float *Imin = new float[1];
//	float *Imax = new float[1];
//	float *Emin = new float[1];
//	float *Emax = new float[1];
//	float *Amin = new float[1];
//	float *Amax = new float[1];
//	float *Vmin = new float[1];
//	float *Vmax = new float[1];
//
//	//for (int i = 0; i < argc; i++) {
//	//	cout << "\n" << i + 1 << ": " << argv[i] << "    ";
//	//}
//	if (argc == 3) {
//		string resultPath = string(argv[2]);
//		string logResultPath = resultPath.substr(0, resultPath.size() - 4) + "_log.tif";
//		Images::LogKohSpecified(argv[1], logResultPath);
//		Images::normKoh(logResultPath, argv[2]);
//	}
//	else {
//		cout << "inputFile outputFile\n\ninputFile: macierz koherencji\noutputFile: wynik normalizacji";
//	}
//	
//	return 0;
//}