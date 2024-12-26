Компиляция OMP и OMP+MPI:

mpixlC <file.omp> -std=c++11 -qsmp=omp -O3 -o <file.cpp>


Компиляция MPI:

module load SpectrumMPI

mpixlC <file.omp> -std=c++11 -O3 -o <file.cpp>


Запуск OMP:

bsub -n 1 -o omp.out -e omp.err -R "affinity[core(N)]" -m "polus-c3-ib polus-c4-ib" -q normal "OMP_NUM_THREADS=N mpiexec ./omp"

N - число нитей


Запуск MPI:

bsub -n N -o mpi.out -e omp160180-16.err -m "polus-c3-ib polus-c4-ib" -q normal "mpiexec ./mpi"

N - число процессов

Запуск MPI+OMP:

bsub -n M -o mpi+omp.out -e mpi+omp.err -R "affinity[core(N)]" -m "polus-c3-ib polus-c4-ib" -q normal "OMP_NUM_THREADS=N mpiexec ./mpi+omp"

M - число процессов
N - число нитей
