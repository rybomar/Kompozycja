#include "Images.h"



int main(int argc, char * argv[]) {
	string manual = " resFilePath firstImgPath bandNumInFirstImg [nextImgPath bandNumInNextImg]*\n Wynik zapisywany z dokaladnoscia zalezna od danych wejsciowych. Obrazy musza miec takie same wymiary w px. Nie nadpisuje istniejacych plikow.";
	if (argc == 2 && (argv[1] == "-man" || argv[1] == "-manual")) {
		cout << endl << manual;
	}
	if (argc < 2 || argc % 2 != 0) {
		cout << "\nZla liczba argumentow.";
		cout << "\n" << manual;
		return 0;
	}
	GDALAllRegister();
	Images ims = Images();
	int cut = 0;
	for (int i = 2; i<argc; i += 2) {
		string fileName = argv[i];
		if (fileName == "-cut") {
			string scut = argv[i + 1];
			cut = stoi(scut);
			cout << "\nObraz zostanie obciety o " << cut << " px przy kazdeja krawedzi.";
			continue;
		}
		string sBandNumber = argv[i + 1];
		int bandNumber = stoi(sBandNumber);
		if (!ims.addBand(fileName, bandNumber)) {
			cout << "\nWARNING: " << "Nie dodano kanalu: " << argv[i] << " (" << argv[i + 1] << ") ";
		}
		else {
			//cout << "\nOdczytano kanal " << i/2 << "/"<<(argc-2)/2;
		}
	}

	if (!ims.save(argv[1], cut, cut, cut, cut)) {
		cout << "\nERROR!: Nie zapisano, sprawdz sciezke pliku. Jesli plik istnieje to usun go.";
	}
	return 0;
}