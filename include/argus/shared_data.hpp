#ifndef SHARED_DATA_H
#define SHARED_DATA_H

#include <Eigen/Dense> // Assuming you include your Eigen types here
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <cstring> // For std::memcpy

// --- CONSTANTS ---
// Define maximum buffer sizes for pre-allocation in shared memory
constexpr size_t MAX_FLOAT_BUFFER_SIZE = 10 * 1024 * 1024; // 10MB
constexpr size_t MAX_INT_BUFFER_SIZE = 5 * 1024 * 1024; // 5MB
// Max sizes for shared memory buffers (adjust as needed for your simulation)
constexpr size_t MAX_DATA_SIZE = 1 * 1024 * 1024; // 1MB for each float buffer
// --- DATA PAYLOAD STRUCTS (Flattened for Shared Memory) ---

// Data structure for Process 1 (Simulation) -> Process 2 (Processor)
struct SimToProcData
{
	int mat_rows;
	int mat_cols;
	int vec_size;
	// Flat buffer to store the Eigen::MatrixXf and Eigen::VectorXf data
	double float_matrix_buffer[MAX_FLOAT_BUFFER_SIZE / sizeof(double)];
	double float_vector_buffer[MAX_FLOAT_BUFFER_SIZE / sizeof(double)];
};

// Data structure for Process 2 (Processor) -> Process 1 (Simulation/Visualization)
struct ProcToSimData
{
	int mat_i_rows;
	int mat_i_cols;
	int mat_f_rows;
	int mat_f_cols;
	// Flat buffer for Eigen::MatrixXi
	int int_matrix_buffer[MAX_INT_BUFFER_SIZE / sizeof(int)];
	// Flat buffer for Eigen::MatrixXf
	double float_matrix_buffer[MAX_FLOAT_BUFFER_SIZE / sizeof(float)];
};

// --- MAIN SHARED MEMORY STRUCTURE ---

using namespace boost::interprocess;

struct SharedMemoryData
{
	// ------------------------------------
	// 1. SIM -> PROC CHANNEL (For Sim Data)
	// ------------------------------------
	SimToProcData sim_to_proc;

	// Mutex to protect 'sim_to_proc' and the 'sim_ready' flag
	interprocess_mutex sim_mutex;
	// Condition variable for Process 2 to wait for new simulation data
	interprocess_condition sim_condition;
	// Flag to track if new data is available
	bool sim_ready = false;

	// ------------------------------------
	// 2. PROC -> SIM CHANNEL (For Processed Data)
	// ------------------------------------
	ProcToSimData proc_to_sim;

	// Mutex to protect 'proc_to_sim' and the 'proc_ready' flag
	interprocess_mutex proc_mutex;
	// Condition variable for Process 1 to wait for new processed data
	interprocess_condition proc_condition;
	// Flag to track if new data is available
	bool proc_ready = false;

	// Constructor required for construction in managed_shared_memory
	SharedMemoryData() = default;
};

#endif // SHARED_DATA_H
