// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <sstream>
// #include <map>
//
// using std::cout; using std::cerr;
// using std::endl; using std::string;
// using std::ifstream; using std::ostringstream;
// using std::istringstream;
//
// // string readFileIntoString(const string& path) {
// //     auto ss = ostringstream{};
// //     ifstream input_file(path);
// //     if (!input_file.is_open()) {
// //         cerr << "Could not open the file - '"
// //              << path << "'" << endl;
// //         exit(EXIT_FAILURE);
// //     }
// //     ss << input_file.rdbuf();
// //     return ss.str();
// // }
//
// int main()
// {
//     string filename("grades.csv");
//     string file_contents;
//     std::map<int, std::vector<string>> csv_contents;
//     char delimiter = ',';
//
//     file_contents = readFileIntoString(filename);
//
//     istringstream sstream(file_contents);
//     std::vector<string> items;
//     string record;
//
//     int counter = 0;
//     while (std::getline(sstream, record)) {
//         istringstream line(record);
//         while (std::getline(line, record, delimiter)) {
//             items.push_back(record);
//         }
//
//         csv_contents[counter] = items;
//         items.clear();
//         counter += 1;
//     }
//
//     exit(EXIT_SUCCESS);
// }

#include <iostream>
#include <fstream>



double one_column_bin_reader(int idx)
{
  // open file
  std::ifstream file("CSF_TMZ.bin", std::ios::in|std::ios::binary|std::ios::ate);

  if (file.is_open())
  {
    // std::cout<<"file opened"<<std::endl;

    // stream position
    std::streampos size = file.tellg();

    // std::cout<<"Size : "<<size<<std::endl;

		char * memblock = new char [size];

	  file.seekg (0, std::ios::beg);
	  file.read (memblock, size);

		//reinterpret as doubles
		double* double_values = (double*)memblock;

		std::cout<<"conc returning: "<<	double_values[idx]<<std::endl;

		delete[] memblock;
    file.close();
		return double_values[idx];
  }
  else
  {
    std::cout<<"failed"<<std::endl;
  }
  return 0;
}


int main()
{
  one_column_bin_reader(10);
  return 0;
}
