#pragma once

#include <chrono>
#include <fstream>
#include <climits>
#include <unistd.h>

uint64_t ram_size() {
    uint64_t pages = sysconf(_SC_PHYS_PAGES);
    uint64_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages*page_size;
}

std::chrono::steady_clock::time_point now() {
    return std::chrono::steady_clock::now();
}

std::string format_time(uint64_t ns) {
    std::string time_str;

    if (ns > 10000000000) {
        time_str = std::to_string(ns/1000000000) + " s";
    } else if (ns > 10000000) {
        time_str = std::to_string(ns/1000000) + " ms";
    } else if (ns > 10000) {
        time_str = std::to_string(ns/1000) + " us";
    } else {
        time_str = std::to_string(ns) + " ns";
    }

    return time_str;
}

std::string format_size(uint64_t B) {
    std::string size_str;

    if (B > 10000000000) {
        size_str = std::to_string(B/1000000000) + " GB";
    } else if (B > 10000000) {
        size_str = std::to_string(B/1000000) + " MB";
    } else if (B > 10000) {
        size_str = std::to_string(B/1000) + " KB";
    } else {
        size_str = std::to_string(B) + " B";
    }
    
    return size_str;
}

uint64_t time_diff_ns(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
}

uint64_t time_diff_ns(std::chrono::steady_clock::time_point t) {
    return time_diff_ns(t,std::chrono::steady_clock::now());
}

std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    std::cout << ", in ~ " << format_time(time_diff_ns(t1,t2)) << std::endl;
    return std::chrono::steady_clock::now();
}

std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t) {
    return log_runtime(t,std::chrono::steady_clock::now());
}

void log_message(std::string message) {
    std::cout << message << std::flush;
}

void header_error() {
	std::cout << "Error: malformed header in patterns file" << std::endl;
	std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
	exit(0);
}

uint64_t get_number_of_patterns(std::string header) {
	uint64_t start_pos = header.find("number=");

	if (start_pos == std::string::npos or start_pos + 7 >= header.size()) {
		header_error();
	}

	start_pos += 7;
	uint64_t end_pos = header.substr(start_pos).find(" ");

	if (end_pos == std::string::npos) {
		header_error();
	}

	uint64_t n = std::atoi(header.substr(start_pos).substr(0, end_pos).c_str());
	return n;
}

uint64_t get_patterns_length(std::string header) {
	uint64_t start_pos = header.find("length=");

	if (start_pos == std::string::npos or start_pos + 7 >= header.size()) {
		header_error();
	}

	start_pos += 7;
	uint64_t end_pos = header.substr(start_pos).find(" ");

	if (end_pos == std::string::npos) {
		header_error();
	}

	uint64_t n = std::atoi(header.substr(start_pos).substr(0, end_pos).c_str());
	return n;
}

void read_from_file(std::ifstream& in, uint8_t* data, uint64_t size) {
    uint64_t size_left = size;
    uint64_t bytes_to_read;

    while (size_left > 0) {
        bytes_to_read = std::min(size_left,(uint64_t)INT_MAX);
        in.read((char*)&data[size-size_left],bytes_to_read);
        size_left -= bytes_to_read;
    }
}

void write_to_file(std::ofstream& out, uint8_t* data, uint64_t size) {
    uint64_t size_left = size;
    uint64_t bytes_to_write;

    while (size_left > 0) {
        bytes_to_write = std::min(size_left,(uint64_t)INT_MAX);
        out.write((char*)&data[size-size_left],bytes_to_write);
        size_left -= bytes_to_write;
    }
}

template <typename uint_t>
inline static uint_t log2_clz(uint_t x) {
    if constexpr (sizeof (uint_t) < 8) {
        return 31 - __builtin_clz(x);
    } else {
        return 63 - __builtin_clzl(x);
    }
}

template <typename T, typename Alloc = std::allocator<T>>
class default_init_allocator : public Alloc {
    using a_t = std::allocator_traits<Alloc>;
    using Alloc::Alloc;

    public:
    template<typename U> struct rebind {};
    template<typename U> void construct (U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value) {::new(static_cast<void*>(ptr)) U;}
    template<typename U, typename... Args> void construct (U* ptr, Args&&... args) {}
};