std::ofstream debugFile;
//int rank;

// debug output to stderr
#define DEBUG_PRINT_ENABLED 1
#define DEBUG_FILE "/user/kmarcus2/Elemental-code/debug/%i.debug"
#define DEBUG_FORMAT "DEBUG " << rank << " [" << __FILE__ << ":" << __LINE__ << " " << __DATE__ << " " << __TIME__ << "] "
#define DEBUG(x) do { if (DEBUG_PRINT_ENABLED) debugFile << DEBUG_FORMAT << x << std::endl; } while (0)
#define DEBUG_VAR(x) do { if (DEBUG_PRINT_ENABLED) debugFile << DEBUG_FORMAT << #x << ": " << x << std::endl; } while (0)

// print macros for stdout
#define PRINT(x) do { std::cout << x << std::endl; } while (0)
#define PRINT_VAR(x) do { std::cout << #x << ": " << x << std::endl; } while (0)

// ends the mpi program and flushes debug output
#define END_PROGRAM \
    if (DEBUG_PRINT_ENABLED) { \
        debugFile.flush(); \
        debugFile.close(); \
    } \
    MPI_Barrier(MPI_COMM_WORLD); \
    MPI_Finalize(); \
    return 1;
    


void libMPIdebugInit() {
    
    if (DEBUG_PRINT_ENABLED) {
        char buff[128];
        std::cout << rank << std::endl;
        sprintf(buff, DEBUG_FILE, rank);
        debugFile.open(buff);
    }
    
}