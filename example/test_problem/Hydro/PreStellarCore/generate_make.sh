# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --nlevel=18 \
                       --flu_scheme=MHM_RP --slope=PPM --flux=HLLC --mhd=False \
                       --gravity=True --pot_scheme=SOR --unsplit_gravity=True --fftw=FFTW3 \
                       --eos=MULTIGAMMA --barotropic=True \
                       --hdf5=True --mpi=True --gpu=True \
                       --bitwise_reproducibility=True --debug=False --double=False \
                       "$@"
