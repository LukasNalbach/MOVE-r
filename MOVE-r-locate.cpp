#include <iostream>
#include <filesystem>
#include <malloc_count.h>
#include <ips4o.hpp>
#include <misc/utils.hpp>
#include <MOVE-r/MOVE-r.hpp>

int ptr = 1;
uint64_t p_d = 0;
std::string path_index_file;
std::string path_textfile = "";
std::string path_outputfile;
std::string name_textfile;
std::ifstream index_file;
std::string path_patternsfile;
std::ofstream mf;

void help() {
	std::cout << "MOVE-r-locate: locate all occurrences of the input patterns." << std::endl << std::endl;
	std::cout << "usage: MOVE-r-locate [options] <index_file> <patterns>" << std::endl;
	std::cout << "   -c <text_file>             check correctness of each pattern occurrence on" << std::endl;
	std::cout << "                              this text file (must be the indexed text file)" << std::endl;
	std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
	std::cout << "                              text_name should be the name of the original file" << std::endl;
	std::cout << "   -o <output_file>           write pattern occurrences to this file (ASCII)" << std::endl;
	std::cout << "   -pd <integer>              the number of threads to use when reconstructing" << std::endl;
	std::cout << "                              the index (default: greatest possible)" << std::endl;
	std::cout << "   <index_file>               index file (with extension .MOVE-r)" << std::endl;
	std::cout << "   <patterns_file>            file in pizza&chili format containing the patterns" << std::endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
	std::string s = argv[ptr];
	ptr++;

	if (s == "-c") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -c option." << std::endl;
			help();
		}

		path_textfile = argv[ptr];
		ptr++;
	} else if (s == "-m") {
		if (ptr >= argc-1) {
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
	} else if (s == "-o") {
		if (ptr >= argc-1) {
			std::cout << "error: missing parameter after -o option." << std::endl;
			help();
		}

		path_outputfile = argv[ptr];
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
void locate() {
	uint64_t size_index_file = std::filesystem::file_size(std::filesystem::path(path_index_file));
	std::cout << "index file size: " << format_size(size_index_file) << std::endl;

	if (mf.is_open()) {
		mf << "RESULT"
		   << " type=locate"
		   << " text=" << name_textfile
		   << " index_impl=move_r";
	}

	std::cout << std::setprecision(4);
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
    auto t0 = now();
    uint8_t* T;
    bool check_correctness = false;
    std::ofstream output_file(path_outputfile);

    if (path_textfile != "") {
		std::cout << "reading the original text file" << std::flush;	
    	check_correctness = true;
		std::ifstream text_file(path_textfile);
		text_file.seekg(0,std::ios::end);
		uint64_t n = text_file.tellg();
		text_file.seekg(0);
		T = (uint8_t*) malloc(n);
		read_from_file(text_file,&T[0],n);
		text_file.close();
		log_runtime(t0);
    }

	std::cout << "reconstructing the index " << path_index_file << " using " << std::to_string(p_d) << " threads" << std::flush;
    auto t1 = now();
	MOVE_r<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs> move_r;
	move_r.load(index_file,revert_count_locate,true,p_d,mf);
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
	std::getline(patterns_file,header);
	uint64_t num_patterns = get_number_of_patterns(header);
	uint64_t m = get_patterns_length(header);
	uint64_t perc;
	uint64_t last_perc = 0;
	uint64_t occ_tot = 0;
	uint64_t time_locate = 0;
	uint64_t count;
	bool equal;
	uint8_t* P = (uint8_t*) malloc(m);
	std::chrono::steady_clock::time_point t3,t4;
	std::vector<pos_t> Occ(0);

	// extract patterns from file and search them in the index
	for (uint64_t i=0; i<num_patterns; i++) {
		perc = (100*i) / num_patterns;

		if (perc > last_perc) {
			std::cout << perc << "% done .." << std::endl;
			last_perc=perc;
		}

		patterns_file.read((char*)P,m);

		t3 = now();
		move_r.locate(P,m,Occ);
		t4 = now();
		time_locate += time_diff_ns(t3,t4);
		occ_tot += Occ.size();
		bool is_sorted = false;

		// check occurrences
		if (check_correctness) {
			// remove duplicates, if any (there shouldn't be!)
			ips4o::sort(Occ.begin(),Occ.end());
			is_sorted = true;
			Occ.resize(std::distance(Occ.begin(),unique(Occ.begin(),Occ.end())));
			count = move_r.count(P,m);

			if (Occ.size() != count) {
				std::cout << "error: wrong number of located occurrences: " << Occ.size() << "/" << count << std::endl;
				exit(0);
			}

			for (auto o:Occ) {
				equal = true;

				for (pos_t j=0; j<m; j++) {
					if (T[o+j] != P[j]) {
						equal = false;
						break;
					}
				}

				if (!equal) {
					std::cout << "error: wrong occurrence: " << o << " ("  << occ_tot << " occurrences) "<< std::endl;

					for (pos_t i=0; i<m; i++) {
						std::cout << T[o+i];
					}

					std::cout << std::endl << std::endl << "/" << std::endl << std::endl;

					for (pos_t i=0; i<m; i++) {
						std::cout << P[i];
					}

					std::cout << std::endl;
					break;
				}
			}
		}

		if (path_outputfile != "") {
			if (!is_sorted) {
				ips4o::sort(Occ.begin(),Occ.end());
			}
			output_file.write((char*)&Occ[0],Occ.size());
		}

		Occ.resize(0);
	}

	free(P);
	P = NULL;

	if (check_correctness) {
		free(T);
		T = NULL;
	}

	patterns_file.close();
	double occ_avg = (double) occ_tot / num_patterns;

	std::cout << "average occurrences per pattern: " << occ_avg << std::endl;
	std::cout << "number of patterns: " << num_patterns << std::endl;
	std::cout << "pattern length: " << m << std::endl;
	std::cout << "total number of occurrences: " << occ_tot << std::endl;
	std::cout << "locate time: " << format_time(time_locate) << std::endl;
	std::cout << "             " << format_time(time_locate/num_patterns) << "/pattern" << std::endl;
	std::cout << "             " << format_time(time_locate/occ_tot) << "/occurrence" << std::endl << std::endl;

	if (mf.is_open()) {
		mf << " a=" << move_r.get_a();
		mf << " n=" << move_r.get_n();
		mf << " sigma=" << std::to_string(move_r.get_sigma());
		mf << " r=" << move_r.get_r();
		mf << " r_=" << move_r.get_r_();
		mf << " r__=" << move_r.get_r__();
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
		mf << " time_locate=" << time_locate << std::endl;
	}
}

int main(int argc, char** argv) {
	if (argc < 3) {
		help();
	}

	while (ptr < argc - 2) {
		parse_args(argv, argc, ptr);
	}

	path_index_file = argv[ptr];
	path_patternsfile = argv[ptr+1];
	index_file.open(path_index_file);
	uint8_t support_uint;
	index_file.read((char*)&support_uint,sizeof(uint8_t));

	if (support_uint != 2) {
		std::cout << "error: the index has no locate support" << std::endl;
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
		case 1: locate<uint8_t,uint8_t,1,1,1>();break;
		case 2: switch (bytes_offs) {
			case 1: switch (bytes_offs) {
				case 1: locate<uint16_t,uint8_t,2,1,1>();break;
				case 2: locate<uint16_t,uint8_t,2,1,2>();break;}break;
			case 2: switch (bytes_offs) {
				case 1: locate<uint16_t,uint16_t,2,2,1>();break;
				case 2: locate<uint16_t,uint16_t,2,2,2>();break;}break;}
		case 3: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: locate<uint32_t,uint8_t,3,1,1>();break;
				case 2: locate<uint32_t,uint8_t,3,1,2>();break;
				case 3: locate<uint32_t,uint8_t,3,1,3>();break;}break;
			case 2: switch (bytes_offs) {
				case 1: locate<uint32_t,uint16_t,3,2,1>();break;
				case 2: locate<uint32_t,uint16_t,3,2,2>();break;
				case 3: locate<uint32_t,uint16_t,3,2,3>();break;}break;
			case 3: switch (bytes_offs) {
				case 1: locate<uint32_t,uint32_t,3,3,1>();break;
				case 2: locate<uint32_t,uint32_t,3,3,2>();break;
				case 3: locate<uint32_t,uint32_t,3,3,3>();break;}break;}break;
		case 4: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: locate<uint32_t,uint8_t,4,1,1>();break;
				case 2: locate<uint32_t,uint8_t,4,1,2>();break;
				case 3: locate<uint32_t,uint8_t,4,1,3>();break;
				case 4: locate<uint32_t,uint8_t,4,1,4>();break;}break;
			case 2: switch (bytes_offs) {
				case 1: locate<uint32_t,uint16_t,4,2,1>();break;
				case 2: locate<uint32_t,uint16_t,4,2,2>();break;
				case 3: locate<uint32_t,uint16_t,4,2,3>();break;
				case 4: locate<uint32_t,uint16_t,4,2,4>();break;}break;
			case 3: switch (bytes_offs) {
				case 1: locate<uint32_t,uint32_t,4,3,1>();break;
				case 2: locate<uint32_t,uint32_t,4,3,2>();break;
				case 3: locate<uint32_t,uint32_t,4,3,3>();break;
				case 4: locate<uint32_t,uint32_t,4,3,4>();break;}break;
			case 4: switch (bytes_offs) {
				case 1: locate<uint32_t,uint32_t,4,4,1>();break;
				case 2: locate<uint32_t,uint32_t,4,4,2>();break;
				case 3: locate<uint32_t,uint32_t,4,4,3>();break;
				case 4: locate<uint32_t,uint32_t,4,4,4>();break;}break;}break;
		case 5: switch (bytes_idx) {
			case 1: switch (bytes_offs) {
				case 1: locate<uint64_t,uint8_t,5,1,1>();break;
				case 2: locate<uint64_t,uint8_t,5,1,2>();break;
				case 3: locate<uint64_t,uint8_t,5,1,3>();break;
				case 4: locate<uint64_t,uint8_t,5,1,4>();break;
				case 5: locate<uint64_t,uint8_t,5,1,5>();break;}break;
			case 2: switch (bytes_offs) {
				case 1: locate<uint64_t,uint16_t,5,2,1>();break;
				case 2: locate<uint64_t,uint16_t,5,2,2>();break;
				case 3: locate<uint64_t,uint16_t,5,2,3>();break;
				case 4: locate<uint64_t,uint16_t,5,2,4>();break;
				case 5: locate<uint64_t,uint16_t,5,2,5>();break;}break;
			case 3: switch (bytes_offs) {
				case 1: locate<uint64_t,uint32_t,5,3,1>();break;
				case 2: locate<uint64_t,uint32_t,5,3,2>();break;
				case 3: locate<uint64_t,uint32_t,5,3,3>();break;
				case 4: locate<uint64_t,uint32_t,5,3,4>();break;
				case 5: locate<uint64_t,uint32_t,5,3,5>();break;}break;
			case 4: switch (bytes_offs) {
				case 1: locate<uint64_t,uint32_t,5,4,1>();break;
				case 2: locate<uint64_t,uint32_t,5,4,2>();break;
				case 3: locate<uint64_t,uint32_t,5,4,3>();break;
				case 4: locate<uint64_t,uint32_t,5,4,4>();break;
				case 5: locate<uint64_t,uint32_t,5,4,5>();break;}break;
			case 5: switch (bytes_offs) {
				case 1: locate<uint64_t,uint64_t,5,5,1>();break;
				case 2: locate<uint64_t,uint64_t,5,5,2>();break;
				case 3: locate<uint64_t,uint64_t,5,5,3>();break;
				case 4: locate<uint64_t,uint64_t,5,5,4>();break;
				case 5: locate<uint64_t,uint64_t,5,5,5>();break;}break;}break;}

	if (mf.is_open()) {
		mf.close();
	}

	index_file.close();
}