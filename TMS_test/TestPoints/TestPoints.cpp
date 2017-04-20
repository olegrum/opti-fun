// TestPoints.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

using namespace std;

int main(int argc, char *argv[])
{ 
	if (argc != 6) {
		cout << "testpoints.exe M N A B filename\n";
		return 1;
	}

	int M = stoi(argv[1]); // количестов ячеек (M в степени N - колличество точек)
	int N = stoi(argv[2]); // размерность пространства
	double A = stod(argv[3]); // начало отрезка
	double B = stoi(argv[4]); // конец отрезка
	char* file_name = argv[5]; // имя файла
	double cell_size = (B - A) / M;
	int k = 0; //счетчик

	cout << "Starting a test with the following arguments:\n";
	cout << "M = " << M << "\n";
	cout << "N = " << N << "\n";
	cout << "A = " << A << "\n";
	cout << "B = " << B << "\n";
	cout << "File = " << file_name << "\n";

	// массив matrix с флагами false
	int matrix_size = (int)pow(M, N);
	bool* matrix = new bool[matrix_size];
	for (int i = 0; i < matrix_size; i++) {
		matrix[i] = false;
	}

	// flag for when data points are wrong
	bool is_success = true;

	ifstream file(file_name);

	double* point = new double[N];

	while (is_success) {

		// считываем точки с файла 
		for (int n = 0; n < N; n++) {
			file >> point[n];
		}

		if (file.eof()) {
			break;
		}

		cout << "Testing point (";
		for (int n = 0; n < N; n++) {
			cout << point[n];
			if (n < N - 1) {
				cout << ", ";
			}
		}
		cout << ")\n";
		k = k + 1;

		// вычиление индекса точки в массиве матрикс
		// offset(сдвиг) текущего пространства * размерность текущего пространства
		// + offset для следующего пространства [ + next [ + next .. ]].
		int cell_index = 0;
		for (int i = 0; i < N; i++) {
			int offset = (int)((point[i] - A) / cell_size);
			cell_index += offset * (int)pow(M, N - i - 1);
		}

		if (matrix[cell_index] == false) {
			matrix[cell_index] = true;
		}
		else {
			is_success = false;
			cout << "Error: there is a point already in the cell.\n" << k << "\n";
		}
	}

	file.close();
	delete point;

	for (int i = 0; i < matrix_size; i++) {
		if (!matrix[i]) {
			cout << "Error: missing point. \n"<< k<< "\n";
			is_success = false;
			break;
		}
	}



	if (is_success) {
		cout << "Success: All points are on their places.\n" << k << "\n";
	}

	//for (int i=0; i < matrix_size; ++i)
		//cout << matrix[i] << " ";


	delete matrix;

	system("pause");
	return 0;
}
