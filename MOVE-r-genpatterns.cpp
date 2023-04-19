#include <climits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <misc/utils.hpp>

void help() {
    std::cout << "MOVE-r-genpatterns: generate patterns from a file." << std::endl << std::endl;
    std::cout << "usage: MOVE-r-genpatterns <file> <length> <number> <patterns file> <forbidden>" << std::endl;
    std::cout << "       randomly extracts <number> substrings of length <length> from <file>," << std::endl;
    std::cout << "       avoiding substrings containing characters in <forbidden>." << std::endl;
    std::cout << "       The output file, <patterns file> has a first line of the form:" << std::endl;
    std::cout << "       # number=<number> length=<length> file=<file> forbidden=<forbidden>" << std::endl;
    std::cout << "       and then the <number> patterns come successively without any separator" << std::endl;
    exit(0);
}

int main(int argc, char *argv[]) {
    if(argc < 5 || 6 < argc) {
        std::cout << "invalid input: wrong number of arguments" << std::endl;
        help();
    }

    std::ifstream input_file(argv[1]);
    if(!input_file.good()) {
        std::cout << "invalid input: could not read <file>" << std::endl;
        help();
    }

    std::cout << std::setprecision(4);
    uint8_t* T;
    input_file.seekg(0,std::ios::end);
	uint64_t n = input_file.tellg();

    uint64_t m = atoi(argv[2]);

    if (m < 0) {
        input_file.close();
        std::cout << "Error: number of patterns must be >= 1" << std::endl;
        help();
    }

    uint64_t N = atoi(argv[3]);

    if (N < 0 || N >= n) {
        input_file.close();
        std::cout << "Error: length must be >= 1 and <= file length" << std::endl;
        help();
    }

    std::ofstream output_file(argv[4]);
    if(!output_file.is_open()) {
        input_file.close();
        std::cout << "invalid input: could not create <patterns file>" << std::endl;
        help();
    }

    std::string forbidden = "";
    
    if (argc == 6) {
        forbidden = argv[5];
    }

    std::string basename = argv[1];
    basename = basename.substr(basename.find_last_of("/\\")+1);
    output_file << "# number=" << N << " length=" << m << " file=" << basename << " forbidden=\n";
	input_file.seekg(0);
    std::cout << "reading text file" << std::flush;
    auto time = now();
    T = (uint8_t*) malloc(n);
	read_from_file(input_file,&T[0],n);
	input_file.close();
    time = log_runtime(time);
    uint64_t pos_rand;
    bool found_forbidden = false;
    std::vector<bool> is_forbidden;

    if (!forbidden.empty()) {
        is_forbidden.resize(256,false);

        for (uint64_t i=0; i<forbidden.size(); i++) {
            is_forbidden[forbidden[i]] = true;
        }
    }

    std::cout << "generating " << N << " petterns of length " << m << std::flush;

    for(uint64_t i=0; i<N; i++) {
        do {
            pos_rand = std::rand()%(n-m);
            found_forbidden = false;

            if (!forbidden.empty()) {
                for (uint64_t i=0; i<m; i++) {
                    if (is_forbidden[T[pos_rand+i]]) {
                        found_forbidden = true;
                        break;
                    }
                }
            }
        } while (found_forbidden);
        
        output_file.write((char*)&T[pos_rand],m);
    }

    output_file.close();
    time = log_runtime(time);
}