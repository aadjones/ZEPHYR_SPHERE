#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

void parseFile(const string &inName, const string &outName) {
	ifstream inFile(inName);
	string line;

	ofstream outFile(outName, ios::binary);
	
	std::getline(inFile, line);
	string xstr = "\"integer nx\"", ystr = "\"integer ny\"", zstr = "\"integer nz\"";
	string p0str = "\"point p0\"", p1str = "\"point p1\"";
	size_t xPos = line.find(xstr);
	size_t yPos = line.find(ystr);
	size_t zPos = line.find(zstr);


	int cell[4];
	string substr = line.substr(xPos + xstr.length(), yPos-xPos-xstr.length());
	cell[0] = stoi(substr);
	substr = line.substr(yPos + ystr.length(), zPos - yPos - ystr.length());
	cell[1] = stoi(substr);
	substr = line.substr(zPos + zstr.length(), line.length()-zPos-zstr.length());
	cell[2] = stoi(substr);
	cell[3] = 1;

	float bounding[6];
	std::getline(inFile, line);
	size_t p0pos = line.find(p0str), p1pos = line.find(p1str);
	substr = line.substr(p0pos + p0str.length(), p1pos - p0pos - p0str.length());
	//sscanf_s(substr.c_str(), " [ %f %f %f", &(bounding[0]), &(bounding[1]), &(bounding[2]));
  // ADJ: sscanf_s is platform-dependent (Microsoft).
	sscanf(substr.c_str(), " [ %f %f %f", &(bounding[0]), &(bounding[1]), &(bounding[2]));
	substr = line.substr(p1pos + p1str.length(), line.length() - p1pos - p1str.length());
	//sscanf_s(substr.c_str(), " [ %f %f %f", &(bounding[3]), &(bounding[4]), &(bounding[5]));
	sscanf(substr.c_str(), " [ %f %f %f", &(bounding[3]), &(bounding[4]), &(bounding[5]));

	std::getline(inFile, line);

	stringstream ss(line);
	string number;

	int channel = 1;
	char str[5] = "VOL\3";
	outFile.write(str, 4 * sizeof(char));
	outFile.write((char*)&channel, sizeof(int));
	outFile.write((char*)cell, 4 * sizeof(int));
	outFile.write((char*)bounding, 6 * sizeof(float));

	try {
		while (std::getline(inFile, line)) {
			stringstream ss(line);
			int sum = 0;
			while (getline(ss, number, ' ')) {
				float tmp = stof(number);
				outFile.write((char*)&tmp, sizeof(float));
				++sum;
			}

			cout << sum << " " << sizeof(float);
		}
	}
	catch (exception e) {
		cout << endl << "non-number, finish parsing" << endl;
	}

	outFile.close();
	inFile.close();
}


int main(int argc, char* argv[]) {
	std::string inName(argv[1]);
  std::string outName(argv[2]);
	parseFile(inName, outName);

	// system("pause");
  system("sleep 0.01");
  return 0;
}
