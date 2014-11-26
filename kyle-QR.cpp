#include "elemental.hpp"
using namespace elem;

int main( int argc, char* argv[] )
{
Initialize( argc, argv );

    DistMatrix<double> A;

    Read( A, "matrix.txt", ASCII );

    // forme an explicit QR decomposition,
    // which returns the matrices Q and R directly.
    DistMatrix<double> R;
    qr::Explicit( A, R );
    // NOTE: A has been overwritten with Q

    Finalize();
    return 0;

}
