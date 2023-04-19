#include <iostream>
#include <filesystem>
#include <malloc_count.h>
#include <misc/utils.hpp>
#include <MOVE-r/MOVE-r.hpp>

int ptr = 1;
uint64_t p_d = 0;
std::ofstream mf;
std::string path_index_file;
std::string name_textfile;
std::ifstream index_file;
std::string path_patternsfile;

void help() {
	std::cout << "MOVE-r-count: count all occurrences of the input patterns." << std::endl << std::endl;
	std::cout << "usage: MOVE-r-count <index_file> <patterns_file>" << std::endl;
	std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
	std::cout << "                              text_name should be the name of the original file" << std::endl;
	std::cout << "   -pd <integer>              the number of threads to use when reconstructing" << std::endl;
	std::cout << "                              the index (default: greatest possible)" << std::endl;
	std::cout << "   <index_file>               index file (with extension .MOVE-r)" << std::endl;
	std::cout << "   <patterns_file>            file in pizza&chili format containing the patterns." << std::endl;
	exit(0);
}

void parse_args(char **argv, int argc, int &ptr) {
	std::string s = argv[ptr];
	ptr++;

	if (s == "-m") {
		if (ptr >= argc - 1) {
			std::cout << "error: missing parameter after -o option." << std::endl;
			help();
		}

		std::string path_mf = argv[ptr];
		mf.open(path_mf,std::filesystem::exists(path_mf) ? std::ios::app : std::ios::out);

		if (!mf.good()) {
			std::cout << "error: cannot open measurement file" << std::endl;
			help();
		}

		ptr++;
		name_textfile = argv[ptr];
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
		} else if (p_d > omp_get_max_threads()) {
			std::cout << "error: the specified number of threads to use during reconstruction" << std::endl;
			std::cout << "       is greater than the number of available threads" << std::endl;
			help();
		}

		ptr++;
	} else {
		std::cout << "error: unknown option " << s << std::endl;
		help();
	}
}

template <typename pos_t, typename idx_t, unsigned int bytes_pos, unsigned int bytes_idx, unsigned int bytes_offs>
void count() {
	uint64_t size_index_file = std::filesystem::file_size(std::filesystem::path(path_index_file));
	std::cout << "index file size: " << format_size(size_index_file) << std::endl;

	if (mf.is_open()) {
		mf << "RESULT"
		   << " type=count"
		   << " text=" << name_textfile
		   << " index_impl=move_r";
	}

	index_file.seekg(6);
	uint16_t p_d_;
	index_file.read((char*)&p_d_,sizeof(uint16_t));
	if (p_d == 0) {
		p_d = std::min(p_d_,(uint16_t)omp_get_max_threads());
	} else if (p_d > p_d_) {
		std::cout << "error: the specified number of threads to use is greater than the maximum" << std::endl;
		std::cout << "       number that has been specified when the index was built" << std::endl;
		help();
	}
	index_file.seekg(0);

	std::cout << "reconstructing the index " << path_index_file << " using " << std::to_string(p_d) << " threads" << std::flush;
	auto t1 = now();
	MOVE_r<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs> move_r;
	move_r.load(index_file,revert_count,true,p_d,mf);
	uint64_t size_index_in_ram = move_r.get_size();
	auto t2 = now();
	uint64_t time_reconstruct_index = time_diff_ns(t1,t2);
	std::cout << "index size in ram: " << format_size(size_index_in_ram) << std::endl;
	move_r.log_data_structure_sizes();
	std::cout << std::endl << "searching patterns ... " << std::endl;
	std::ifstream patterns_file(path_patternsfile);

	/* read header of the pizza&chilli input file, header example:
	# number=7 length=10 file=genome.fasta forbidden=\n\t */
	std::string header;
	std::getline(patterns_file, header);
	uint64_t num_patterns = get_number_of_patterns(header);
	uint64_t m = get_patterns_length(header);
	uint64_t perc;
	uint64_t last_perc = 0;
	uint64_t occ_tot = 0;
	uint64_t time_count = 0;
	uint8_t* P = (uint8_t*) malloc(m);
	std::chrono::steady_clock::time_point t3,t4;

	// extract patterns from file and search them in the index
	for (uint64_t i=0; i<num_patterns; i++) {
		perc = (100*i) / num_patterns;

		if (perc > last_perc) {
			std::cout << perc << "% done .." << std::endl;
			last_perc = perc;
		}

		patterns_file.read((char*)P,m);

		t3 = now();
		occ_tot += move_r.count(P,m);
		t4 = now();
		time_count += time_diff_ns(t3,t4);
	}

	free(P);
	P = NULL;

	patterns_file.close();
	double occ_avg = (double) occ_tot / num_patterns;
	std::cout << "average occurrences per pattern: " << occ_avg << std::endl;
	std::cout << "number of patterns: " << num_patterns << std::endl;
	std::cout << "pattern length: " << m << std::endl;
	std::cout << "total number of occurrences  occ_t: " << occ_tot << std::endl;
	std::cout << "count time: " << format_time(time_count) << std::endl;
	std::cout << "            " << format_time(time_count/num_patterns) << "/pattern" << std::endl;
	std::cout << "            " << format_time(time_count/occ_tot) << "/occurrence" << std::endl << std::endl;

	if (mf.is_open()) {
		mf << " a=" << move_r.get_a();
		mf << " n=" << move_r.get_n();
		mf << " sigma=" << std::to_string(move_r.get_sigma());
		mf << " r=" << move_r.get_r();
		mf << " r_=" << move_r.get_r_();
		mf << " omega_p=" << std::to_string(bytes_pos*8);
		mf << " omega_idx=" << std::to_string(bytes_idx*8);
		mf << " omega_offs=" << std::to_string(bytes_offs*8);
		mf << " size_index_file=" << size_index_file;
		mf << " size_index_in_ram=" << size_index_in_ram;
		move_r.log_data_structure_sizes(mf);
		mf << " m=" << m;
		mf << " num_patterns=" << num_patterns;
		mf << " occ_tot=" << occ_tot;
		mf << " time_reconstruct_index=" << time_reconstruct_index;
		mf << " time_count=" << time_count;
		mf << std::endl;
	}
}

int main(int argc, char **argv) {
	if (argc < 3) {
		help();
	}

	while (ptr < argc - 2) {
		parse_args(argv, argc, ptr);
	}

	std::cout << std::setprecision(4);
	path_index_file = argv[ptr];
	path_patternsfile = argv[ptr+1];
	index_file.open(path_index_file);
	uint8_t support_uint;
	index_file.read((char*)&support_uint,sizeof(uint8_t));

	if (support_uint < 1) {
		std::cout << "error: the index has no count support" << std::endl;
		help();
	}

	uint8_t bytes_pos;
	uint8_t bytes_idx;
	uint8_t bytes_offs;
	index_file.read((char*)&bytes_pos,sizeof(uint8_t));
	index_file.read((char*)&bytes_idx,sizeof(uint8_t));
	index_file.read((char*)&bytes_offs,sizeof(uint8_t));
	index_file.seekg(0);

	switch (bytes_pos) {
		case 1: count<uint8_t,uint8_t,1,1,1>();break;
		case 2: switch (bytes_offs) {
			case 1: switch (bytes_offs) {
				case 1: count<uint16_t,uint8_t,2,1,1>();break;
				case 2: count<uint16_t,uint8_t,2,1,2>();break;}break;
			case 2: switch (bytes_offs) {
				case 1: count<uint16_t,uint16_t,2,2,1>();break;
				case 2: count<uint16_t,uint16_t,2,2,2>();break;}break;}
		case 3: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: count<uint32_t,uint8_t,3,1,1>();break;
				case 2: count<uint32_t,uint8_t,3,1,2>();break;
				case 3: count<uint32_t,uint8_t,3,1,3>();break;}break;
			case 2: switch (bytes_offs) {
				case 1: count<uint32_t,uint16_t,3,2,1>();break;
				case 2: count<uint32_t,uint16_t,3,2,2>();break;
				case 3: count<uint32_t,uint16_t,3,2,3>();break;}break;
			case 3: switch (bytes_offs) {
				case 1: count<uint32_t,uint32_t,3,3,1>();break;
				case 2: count<uint32_t,uint32_t,3,3,2>();break;
				case 3: count<uint32_t,uint32_t,3,3,3>();break;}break;}break;
		case 4: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: count<uint32_t,uint8_t,4,1,1>();break;
				case 2: count<uint32_t,uint8_t,4,1,2>();break;
				case 3: count<uint32_t,uint8_t,4,1,3>();break;
				case 4: count<uint32_t,uint8_t,4,1,4>();break;}break;
			case 2: switch (bytes_offs) {
				case 1: count<uint32_t,uint16_t,4,2,1>();break;
				case 2: count<uint32_t,uint16_t,4,2,2>();break;
				case 3: count<uint32_t,uint16_t,4,2,3>();break;
				case 4: count<uint32_t,uint16_t,4,2,4>();break;}break;
			case 3: switch (bytes_offs) {
				case 1: count<uint32_t,uint32_t,4,3,1>();break;
				case 2: count<uint32_t,uint32_t,4,3,2>();break;
				case 3: count<uint32_t,uint32_t,4,3,3>();break;
				case 4: count<uint32_t,uint32_t,4,3,4>();break;}break;
			case 4: switch (bytes_offs) {
				case 1: count<uint32_t,uint32_t,4,4,1>();break;
				case 2: count<uint32_t,uint32_t,4,4,2>();break;
				case 3: count<uint32_t,uint32_t,4,4,3>();break;
				case 4: count<uint32_t,uint32_t,4,4,4>();break;}break;}break;
		case 5: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: count<uint64_t,uint8_t,5,1,1>();break;
				case 2: count<uint64_t,uint8_t,5,1,2>();break;
				case 3: count<uint64_t,uint8_t,5,1,3>();break;
				case 4: count<uint64_t,uint8_t,5,1,4>();break;
				case 5: count<uint64_t,uint8_t,5,1,5>();break;}break;
			case 2: switch (bytes_offs) {
				case 1: count<uint64_t,uint16_t,5,2,1>();break;
				case 2: count<uint64_t,uint16_t,5,2,2>();break;
				case 3: count<uint64_t,uint16_t,5,2,3>();break;
				case 4: count<uint64_t,uint16_t,5,2,4>();break;
				case 5: count<uint64_t,uint16_t,5,2,5>();break;}break;
			case 3: switch (bytes_offs) {
				case 1: count<uint64_t,uint32_t,5,3,1>();break;
				case 2: count<uint64_t,uint32_t,5,3,2>();break;
				case 3: count<uint64_t,uint32_t,5,3,3>();break;
				case 4: count<uint64_t,uint32_t,5,3,4>();break;
				case 5: count<uint64_t,uint32_t,5,3,5>();break;}break;
			case 4: switch (bytes_offs) {
				case 1: count<uint64_t,uint32_t,5,4,1>();break;
				case 2: count<uint64_t,uint32_t,5,4,2>();break;
				case 3: count<uint64_t,uint32_t,5,4,3>();break;
				case 4: count<uint64_t,uint32_t,5,4,4>();break;
				case 5: count<uint64_t,uint32_t,5,4,5>();break;}break;
			case 5: switch (bytes_offs) {
				case 1: count<uint64_t,uint64_t,5,5,1>();break;
				case 2: count<uint64_t,uint64_t,5,5,2>();break;
				case 3: count<uint64_t,uint64_t,5,5,3>();break;
				case 4: count<uint64_t,uint64_t,5,5,4>();break;
				case 5: count<uint64_t,uint64_t,5,5,5>();break;}break;}break;}

	if (mf.is_open()) {
		mf.close();
	}

	index_file.close();
}