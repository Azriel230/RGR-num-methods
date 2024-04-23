#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

void printSystem(const char* msg, std::vector<std::vector<double>> &matrix, std::vector<double> &vec)
{
	std::cout << msg << std::endl;
	for (int i = 0; i < vec.size(); i++)
	{
		for (int j = 0; j < vec.size(); j++)
		{
			if (j + 1 != vec.size())
				std::cout << std::setprecision(4) << matrix[i][j] << std::setw(10);
			else
				std::cout << std::setprecision(4) << matrix[i][j] << std::setw(2);
		}
		std::cout << '|' << std::setw(6) << vec[i] << std::endl;
	}
	std::cout << std::endl;
}

int main()
{
	std::ifstream file("matrix.txt");
	if (!file.is_open())
	{
		std::cout << "Error open file";
		exit(EXIT_FAILURE);
	}

	//Calculating size of system
	std::string str;
	std::getline(file, str);

	int sizeMatrix = 0;
	for (size_t i = 0; i < str.size(); i++)
	{
		if (str[i] == ' ')
			sizeMatrix++;
	}

	if (str.size() == 0)
	{
		std::cout << "Error! You enter empty system!";
		exit(EXIT_FAILURE);
	}

	if (sizeMatrix == 0)
	{
		std::cout << "Error! System size is 0!";
		exit(EXIT_FAILURE);
	}

	file.seekg(0);

	std::vector<std::vector<double>> matrixCoeff(sizeMatrix, std::vector<double>(sizeMatrix));
	std::vector<double> vectorFreeCoeff(sizeMatrix * 1.0);
	std::vector<double> vectorX_indexes(sizeMatrix);

	for (int i = 0; i < sizeMatrix; i++)
	{
		for (int j = 0; j < sizeMatrix; j++)
		{
			file >> matrixCoeff[i][j];  
		}
		file >> vectorFreeCoeff[i];
		vectorX_indexes[i] = i;
	}

	//print entered system
	printSystem("Entered system: ", matrixCoeff, vectorFreeCoeff);

	//forward stroke of the Gaussian method
	for (int k = 0; k < sizeMatrix - 1; k++)
	{
		int index = 0; //number of col to swap
		int max_a = 0; //abs max element a in #k row

		//choose main element in row
		for (int m = k; m < sizeMatrix; m++)
		{
			double a = abs(matrixCoeff[k][m]);
			if (max_a < a)
			{
				max_a = a;
				index = m;
			}
		}

		//swap #index col and #k col 
		for (int s = 0; s < sizeMatrix; s++)
		{
			double temp = matrixCoeff[s][k];
			matrixCoeff[s][k] = matrixCoeff[s][index];
			matrixCoeff[s][index] = temp;
		}

		double x_swap = vectorX_indexes[index];
		vectorX_indexes[index] = vectorX_indexes[k];
		vectorX_indexes[k] = x_swap;

		//do #k step forward stroke of the Gaussian method
		for (int i = k + 1; i < sizeMatrix; i++)
		{
			double a_coeff = matrixCoeff[i][k] / matrixCoeff[k][k];
			for (int j = k; j < sizeMatrix; j++)
			{
				matrixCoeff[i][j] = matrixCoeff[i][j] - a_coeff * matrixCoeff[k][j];
			}
			vectorFreeCoeff[i] = vectorFreeCoeff[i] - a_coeff * vectorFreeCoeff[k];
		}
	}

	//print final system
	printSystem("final system: ", matrixCoeff, vectorFreeCoeff);

	//inverse of the Gaussian method
	std::vector<double> vectorX(sizeMatrix);
	int lastIndex = sizeMatrix - 1;
	double x_index = vectorFreeCoeff[lastIndex] / matrixCoeff[lastIndex][lastIndex];
	vectorX[vectorX_indexes[lastIndex]] = x_index;

	for (int i = lastIndex - 1; i >= 0; i--)
	{
		int index = vectorX_indexes[i];

		double sumCoeff = 0;
		for (int j = i + 1; j < sizeMatrix; j++)
		{
			sumCoeff += matrixCoeff[i][j] * vectorX[vectorX_indexes[j]];
		}
		x_index = (vectorFreeCoeff[i] - sumCoeff) / matrixCoeff[i][i];
		vectorX[index] = x_index;
	}

	//print answer
	for (int i = 0; i < vectorX.size(); i++)
	{
		std::cout << "x" << i + 1 << std::setw(2) << " = " << vectorX[i] << std::endl;
	}

	return 0;
}