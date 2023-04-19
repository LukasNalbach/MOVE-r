#include <iostream>
#include <filesystem>
#include <climits>
#include <omp.h>
#include <malloc_count.h>
#include <misc/utils.hpp>
#include <MOVE-r/MOVE-r.hpp>

int ptr = 1;
uint16_t p = 0;
uint16_t p_d = 0;
std::string path_index_file;
std::string path_textfile;
std::string name_textfile;
std::ofstream mf;

void help() {
	std::cout << "MOVE-r-revert: reconstruct the original file." << std::endl << std::endl;
	std::cout << "usage: MOVE-r-revert [options] <index_file>" << std::endl;
	std::cout << "   -pr <integer>              number of threads to use while reverting, must not exceed the maximum" << std::endl;
	std::cout << "                              number of threads to use that has been specified when the index was" << std::endl;
	std::cout << "                              built (default: greatest possible)" << std::endl;
	std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
	std::cout << "                              text_name should be the name of the original file" << std::endl;
	std::cout << "   -o <text_file>             output text file" << std::endl;
	std::cout << "   -pd <integer>              the number of threads to use when reconstructing" << std::endl;
	std::cout << "                              the index (default: greatest possible)" << std::endl;
	std::cout << "   <index_file>               index file (with extension .MOVE-r)" << std::endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
	std::string s = argv[ptr];
	ptr++;

	if (s == "-m") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -m option" << std::endl;
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
	} else if (s == "-o") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -o option" << std::endl;
			help();
		}

		path_textfile = argv[ptr];
		ptr++;
	} else if (s == "-pr") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -pr option" << std::endl;
			help();
		}

		p = atoi(argv[ptr]);

		if (p < 1) {
			std::cout << "error: p < 1" << std::endl;
			help();
		}

		if (p > omp_get_max_threads()) {
			std::cout << "error: the specified number of threads to use while reverting" << std::endl;
			std::cout << "       is greater than the number of available threads" << std::endl;
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
void revert_(uint16_t p, std::ifstream& index_file, uint64_t &size_index_file) {
	auto t1 = now();
	MOVE_r<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs> move_r;
	move_r.load(index_file,revert,true,p_d,mf);
	auto t2 = now();
	uint64_t n = move_r.get_n();
	uint64_t size_index_in_ram = move_r.get_size();
	std::cout << "index size in ram: " << format_size(size_index_in_ram) << std::endl;
	move_r.log_data_structure_sizes();
	std::cout << std::endl << "reverting the index" << " using " << std::to_string(p) << " threads" << std::flush;
	uint8_t* T = move_r.revert(p);
	auto t3 = now();
	uint64_t time_reconstruct_index = time_diff_ns(t1,t2);
	uint64_t time_revert = time_diff_ns(t2,t3);
	log_runtime(t2,t3);
	std::cout << "writing the original file" << std::flush;
	std::ofstream text_file(path_textfile);
	write_to_file(text_file,&T[0],n-1);
	text_file.close();
	auto t4 = now();
	free(T);
	T = NULL;
	uint64_t time_write_textfile = time_diff_ns(t3,t4);
	log_runtime(t3,t4);
	uint64_t time_total = time_diff_ns(t1,t4);
	std::cout << "total revert time: " << format_time(time_total) << std::endl << std::endl;

	if (mf.is_open()) {
		mf << " p=" << p;
		mf << " a=" << move_r.get_a();
		mf << " n=" << n;
		mf << " sigma=" << std::to_string(move_r.get_sigma());
		mf << " r=" << move_r.get_r();
		mf << " r_=" << move_r.get_r_();
		mf << " omega_p=" << std::to_string(bytes_pos*8);
		mf << " omega_idx=" << std::to_string(bytes_idx*8);
		mf << " omega_offs=" << std::to_string(bytes_offs*8);
		mf << " size_index_file=" << size_index_file;
		mf << " size_index_in_ram=" << size_index_in_ram;
		move_r.log_data_structure_sizes(mf);
		mf << " time_reconstruct_index=" << time_reconstruct_index;
		mf << " time_revert=" << time_revert;
		mf << " time_write_textfile=" << time_write_textfile;
		mf << " time_total=" << time_total;
		mf << std::endl;
	}
}

int main(int argc, char** argv) {
	if (argc < 3) {
		help();
	}

	while (ptr < argc - 2) {
		parse_args(argv, argc, ptr);
	}

	std::cout << std::setprecision(4);
	path_index_file = argv[ptr];
	std::ifstream index_file(path_index_file);
	index_file.seekg(1);
	uint8_t bytes_pos;
	uint8_t bytes_idx;
	uint8_t bytes_offs;
	index_file.read((char*)&bytes_pos,sizeof(uint8_t));
	index_file.read((char*)&bytes_idx,sizeof(uint8_t));
	index_file.read((char*)&bytes_offs,sizeof(uint8_t));
	uint16_t p_r;
	index_file.read((char*)&p_r,sizeof(uint16_t));

	if (p == 0) {
		p = std::min(p_r,(uint16_t)omp_get_max_threads());
	} else if (p > p_r) {
		std::cout << "error: the specified number of threads to use is greater than the maximum number that has " << std::endl;
		std::cout << "       been specified when the index was built" << std::endl;
		help();
	}

	uint64_t size_index_file = std::filesystem::file_size(std::filesystem::path(path_index_file));
	std::cout << "index file size: " << format_size(size_index_file) << std::endl;

	if (mf.is_open()) {
		mf << "RESULT"
		   << " type=revert"
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

	switch (bytes_pos) {
		case 1: revert_<uint8_t,uint8_t,1,1,1>(p,index_file,size_index_file);break;
		case 2: switch (bytes_offs) {
			case 1: switch (bytes_offs) {
				case 1: revert_<uint16_t,uint8_t,2,1,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint16_t,uint8_t,2,1,2>(p,index_file,size_index_file);break;}break;
			case 2: switch (bytes_offs) {
				case 1: revert_<uint16_t,uint16_t,2,2,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint16_t,uint16_t,2,2,2>(p,index_file,size_index_file);break;}break;}
		case 3: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: revert_<uint32_t,uint8_t,3,1,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint32_t,uint8_t,3,1,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint32_t,uint8_t,3,1,3>(p,index_file,size_index_file);break;}break;
			case 2: switch (bytes_offs) {
				case 1: revert_<uint32_t,uint16_t,3,2,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint32_t,uint16_t,3,2,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint32_t,uint16_t,3,2,3>(p,index_file,size_index_file);break;}break;
			case 3: switch (bytes_offs) {
				case 1: revert_<uint32_t,uint32_t,3,3,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint32_t,uint32_t,3,3,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint32_t,uint32_t,3,3,3>(p,index_file,size_index_file);break;}break;}break;
		case 4: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: revert_<uint32_t,uint8_t,4,1,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint32_t,uint8_t,4,1,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint32_t,uint8_t,4,1,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint32_t,uint8_t,4,1,4>(p,index_file,size_index_file);break;}break;
			case 2: switch (bytes_offs) {
				case 1: revert_<uint32_t,uint16_t,4,2,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint32_t,uint16_t,4,2,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint32_t,uint16_t,4,2,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint32_t,uint16_t,4,2,4>(p,index_file,size_index_file);break;}break;
			case 3: switch (bytes_offs) {
				case 1: revert_<uint32_t,uint32_t,4,3,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint32_t,uint32_t,4,3,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint32_t,uint32_t,4,3,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint32_t,uint32_t,4,3,4>(p,index_file,size_index_file);break;}break;
			case 4: switch (bytes_offs) {
				case 1: revert_<uint32_t,uint32_t,4,4,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint32_t,uint32_t,4,4,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint32_t,uint32_t,4,4,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint32_t,uint32_t,4,4,4>(p,index_file,size_index_file);break;}break;}break;
		case 5: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: revert_<uint64_t,uint8_t,5,1,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint64_t,uint8_t,5,1,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint64_t,uint8_t,5,1,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint64_t,uint8_t,5,1,4>(p,index_file,size_index_file);break;
				case 5: revert_<uint64_t,uint8_t,5,1,5>(p,index_file,size_index_file);break;}break;
			case 2: switch (bytes_offs) {
				case 1: revert_<uint64_t,uint16_t,5,2,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint64_t,uint16_t,5,2,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint64_t,uint16_t,5,2,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint64_t,uint16_t,5,2,4>(p,index_file,size_index_file);break;
				case 5: revert_<uint64_t,uint16_t,5,2,5>(p,index_file,size_index_file);break;}break;
			case 3: switch (bytes_offs) {
				case 1: revert_<uint64_t,uint32_t,5,3,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint64_t,uint32_t,5,3,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint64_t,uint32_t,5,3,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint64_t,uint32_t,5,3,4>(p,index_file,size_index_file);break;
				case 5: revert_<uint64_t,uint32_t,5,3,5>(p,index_file,size_index_file);break;}break;
			case 4: switch (bytes_offs) {
				case 1: revert_<uint64_t,uint32_t,5,4,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint64_t,uint32_t,5,4,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint64_t,uint32_t,5,4,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint64_t,uint32_t,5,4,4>(p,index_file,size_index_file);break;
				case 5: revert_<uint64_t,uint32_t,5,4,5>(p,index_file,size_index_file);break;}break;
			case 5: switch (bytes_offs) {
				case 1: revert_<uint64_t,uint64_t,5,5,1>(p,index_file,size_index_file);break;
				case 2: revert_<uint64_t,uint64_t,5,5,2>(p,index_file,size_index_file);break;
				case 3: revert_<uint64_t,uint64_t,5,5,3>(p,index_file,size_index_file);break;
				case 4: revert_<uint64_t,uint64_t,5,5,4>(p,index_file,size_index_file);break;
				case 5: revert_<uint64_t,uint64_t,5,5,5>(p,index_file,size_index_file);break;}break;}break;}

	if (mf.is_open()) {
		mf.close();
	}

	index_file.close();
}