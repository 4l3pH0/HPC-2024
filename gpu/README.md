## Подгрузка модулей:
module load openmpi
module load pgi

## Компиляция с gpu:

mpic++ -O3 -lm -acc -std=c++11 -ta=tesla:cc60,time -Minfo=accel -o a.out gpu.cpp -Msafeptr

## Запуск

bsub -n 1 -m "polus-c2-ib polus-c3-ib polus-c4-ib" -gpu "num=1:mode=exclusive_process" -q short -o n160_180-1-gpu "mpiexec -n 1 ./a.out"
bsub -n 2 -m "polus-c2-ib polus-c3-ib polus-c4-ib" -R "affinity[core(4,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]" -gpu "num=2:mode=exclusive_process" -q short -o n160_180-2-gpu "mpiexec -n 2 ./a.out"



