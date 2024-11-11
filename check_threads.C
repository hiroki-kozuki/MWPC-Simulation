#include <iostream>
#include <omp.h>

int main() {
	// Get number of available threads.
	int num_threads = omp_get_max_threads();
	
	// Print number of threads.
	std::cout << "Number of available threads: " << num_threads << std::endl;
	return 0;
	}
