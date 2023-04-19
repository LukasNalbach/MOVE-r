#include <iostream>
#include <filesystem>
#include <omp.h>
#include <malloc_count.h>
#include <misc/utils.hpp>
#include <MOVE-r/MOVE-r/build.hpp>

std::string pathprefix_index_file;
std::string name_textfile;
support_mode support = revert_count_locate;
memory_usage_mode memory_usage_mode = automatic;
uint16_t p = omp_get_max_threads();
uint16_t p_r = 256;
uint16_t p_d = 256;
uint16_t a = 8;
double epsilon = 0.125;
std::ofstream mf_idx;
std::ofstream mf_mds;
int ptr = 1;

void help() {
	std::cout << "MOVE-r-build: builds the r-index (extension .MOVE-r is automatically)." << std::endl;
	std::cout << "added to the index file name; the index has size O(r*(1+epsilon)*(a/(a-1)))" << std::endl << std::endl;
	std::cout << "usage: MOVE-r-build [options] <text_file>" << std::endl;
	std::cout << "   -o <basename>     names the index file basename.MOVE-r (default: text_file)" << std::endl;
	std::cout << "   -s <mode>         support mode: revert, revert-count or revert-count-locate" << std::endl;
	std::cout << "                     (default: revert-count-locate)" << std::endl;
	std::cout << "   -mu <mode>         memory usage mode: low, high or auto (low stores some data structures on" << std::endl;
	std::cout << "                     the disk to reduce (-25% max.) the peak memory usage; auto only does this" << std::endl;
	std::cout << "                     if the ram capacity would be exceeded; high does not store anything on the " << std::endl;
	std::cout << "                     disk; default: auto)" << std::endl;
	std::cout << "   -p <integer>      number of threads to use during the construction of the index" << std::endl;
	std::cout << "                     (default: all threads)" << std::endl;
	std::cout << "   -pr <integer>     restricts, how many threads can be used when reverting" << std::endl;
	std::cout << "                     the index afterwards (default: 256)" << std::endl;
	std::cout << "   -pd <integer>     restricts, how many threads can be used when encoding and" << std::endl;
	std::cout << "                     reconstructing the index (default: 256)" << std::endl;
	std::cout << "   -a <integer>      balancing parameter; a must be an integer number and a >= 2 (default: 8)" << std::endl;
	std::cout << "   -e <integer>      epsilon, 0 < epsilon <= 1 (default: 0.125)" << std::endl;
	std::cout << "   -m_idx <m_file>   file to write measurement data of the index construction to" << std::endl;
	std::cout << "   -m_mds <m_file>   file to write measurement data of the construction of the move-datastructures to" << std::endl;
	std::cout << "   <text_file>       input text file" << std::endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
	std::string s = argv[ptr];
	ptr++;

	if (s == "-o") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -o option" << std::endl;
			help();
		}

		pathprefix_index_file = argv[ptr];
		ptr++;
	} else if (s == "-p") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -p option" << std::endl;
			help();
		}

		p = atoi(argv[ptr]);
		ptr++;
	} else if (s == "-mu") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -p option" << std::endl;
			help();
		}

		std::string memory_usage_mode_string = argv[ptr];

		if (memory_usage_mode_string == "low") {
			memory_usage_mode = low;
		} else if (memory_usage_mode_string == "high") {
			memory_usage_mode = high;
		} else if (memory_usage_mode_string == "auto") {
			memory_usage_mode = automatic;
		} else {
			std::cout << "error: unknown mode provided with -mu option" << std::endl;
			help();
		}

		ptr++;
	} else if (s == "-e") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -e option" << std::endl;
			help();
		}

		epsilon = std::stod(argv[ptr]);

		if (epsilon <= 0) {
			std::cout << "error: epsilon <= 0" << std::endl;
			help();
		}

		if (epsilon > 1) {
			std::cout << "error: epsilon > 1" << std::endl;
			help();
		}

		ptr++;
	} else if (s == "-s") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -s option" << std::endl;
			help();
		}

		std::string support_str = argv[ptr];

		if (support_str == "revert") {
			support = revert;
		} else if (support_str == "revert-count") {
			support = revert_count;
		} else if (support_str == "revert-count-locate") {
			support = revert_count_locate;
		} else {
			std::cout << "error: unknown mode provided with -s option" << std::endl;
			help();
		}

		ptr++;
	} else if (s == "-pr") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -pr option" << std::endl;
			help();
		}

		p_r = atoi(argv[ptr]);

		if (p_r < 1) {
			std::cout << "error: p_r < 1" << std::endl;
			help();
		}

		ptr++;
	} else if (s == "-pd") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -pd option" << std::endl;
			help();
		}

		p_d = atoi(argv[ptr]);

		if (p_d < 1) {
			std::cout << "error: pd < 1" << std::endl;
			help();
		}

		ptr++;
	} else if (s == "-a") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -a option" << std::endl;
			help();
		}

		a = atoi(argv[ptr]);

		if (a < 2) {
			std::cout << "error: a < 2" << std::endl;
			help();
		}

		ptr++;
	} else if (s == "-m_idx") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -m_idx option" << std::endl;
			help();
		}

		std::string path_mf_idx = argv[ptr];
		mf_idx.open(path_mf_idx,std::filesystem::exists(path_mf_idx) ? std::ios::app : std::ios::out);

		if (!mf_idx.good()) {
			std::cout << "error: cannot open or create measurement text_file" << std::endl;
			help();
		}

		ptr++;
	} else if (s == "-m_mds") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -m_mds option" << std::endl;
			help();
		}

		std::string path_mf_mds = argv[ptr];
		mf_mds.open(path_mf_mds,std::filesystem::exists(path_mf_mds) ? std::ios::app : std::ios::out);

		if (!mf_mds.good()) {
			std::cout << "error: cannot open or create at least one measurement text_file" << std::endl;
			help();
		}

		ptr++;
	} else {
		std::cout << "error: unrecognized '" << s << "' option" << std::endl;
		help();
	}

	if (!(1 <= p && p <= omp_get_max_threads())) {
		std::cout << "error: incompatible combination of threads and balancing algorithm" << std::endl;
		help();
	}
}

int main(int argc, char** argv) {
	auto time = now();

	if (argc < 2) {
		help();
	}

	while (ptr < argc - 1) {
		parse_args(argv, argc, ptr);
	}

	std::string path_textfile = argv[ptr];

	if (pathprefix_index_file == "") {
		pathprefix_index_file = path_textfile;
	}

	std::cout << std::setprecision(4);
	name_textfile = path_textfile.substr(path_textfile.find_last_of("/\\") + 1);
	std::string path_index_file = pathprefix_index_file.append(".MOVE-r");
	std::cout << "building the r-index of " << path_textfile;
	std::cout << " using " << std::to_string(p) << " threads,";
	std::cout << " a = " << a << " and epsilon = " << epsilon << std::endl;
	std::cout << "the index will be saved to " << path_index_file << std::endl << std::endl;
	std::cout << "reading T" << std::flush;

	malloc_count_reset_peak();
	uint64_t n;
	std::ifstream text_file(path_textfile);

	if (!text_file.good()) {
		std::cout << std::endl << "error: invalid input, could not read textfile" << std::endl;
		help();
	}

	text_file.seekg(0,std::ios::end);
	n = text_file.tellg()+(std::streamsize)+1;
	text_file.seekg(0);
	uint8_t* T = (uint8_t*) malloc(n);
	read_from_file(text_file,&T[0],n-1);
	T[n-1] = 1;
	text_file.close();
	uint64_t time_read_t = time_diff_ns(time,now());
	time = log_runtime(time);
	std::ofstream index_file(path_index_file);

	if (mf_idx.is_open()) {
		mf_idx << "RESULT"
				<< " type=build_index"
				<< " text=" << name_textfile
				<< " index_impl=move_r"
				<< " p=" << p
				<< " a=" << a
				<< " time_read_t=" << time_read_t;
	}

	if (p > n) {
		std::cout << "warning: p is larger than the text length, setting p to 1" << std::endl;
		p = 1;
	}

	if (n+6*254 <= UINT_MAX) {
		if (n+6*254 <= INT_MAX) {
			MOVE_r_build<uint32_t,int32_t>(T,n,index_file,support,memory_usage_mode,p,p_r,p_d,a,epsilon,true,mf_idx,mf_mds,name_textfile);
		} else {
			MOVE_r_build<uint32_t,int64_t>(T,n,index_file,support,memory_usage_mode,p,p_r,p_d,a,epsilon,true,mf_idx,mf_mds,name_textfile);
		}
	} else {
		MOVE_r_build<uint64_t,int64_t>(T,n,index_file,support,memory_usage_mode,p,p_r,p_d,a,epsilon,true,mf_idx,mf_mds,name_textfile);
	}

	if (mf_idx.is_open()) {
		mf_idx.close();
	}
	
	if (mf_mds.is_open()) {
		mf_mds.close();
	}

	index_file.close();
}