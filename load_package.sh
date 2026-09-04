#!/bin/bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "This script should be sourced, not executed. Use 'source ${0}'"
    exit 1
fi

# An option given as an argument is used directly, so this can be sourced from a
# script or a batch job: `source load_package.sh 1`. Without one it asks, as
# before. Reading from a closed stdin previously selected nothing, printed
# "Invalid option selected" and returned successfully, leaving ROOT unloaded and
# the failure to surface later as a confusing import error.
# Chosen through the environment rather than an argument, deliberately. ROOT's
# thisroot.sh locates itself with ${BASH_ARGV[0]}, so a script sourced with any
# argument poisons that lookup and ROOTSYS is derived from the wrong place.
option="${TIDAL_LOAD_OPTION:-}"

if [ -z "$option" ]; then
    echo "Select an option to activate:"
    echo "1: Activate ROOT"
    echo "2: Activate SVFIT"
    echo "3: Activate CP-TOOLS"
    echo "4: Activate OPTION4"

    if ! read -p "Enter the option number: " option; then
        echo "No option given and no terminal to ask. Set one, for example:" >&2
        echo "   TIDAL_LOAD_OPTION=1 source load_package.sh   # ROOT" >&2
        return 1
    fi
fi

case $option in
    1)
        echo "Load ROOT..."
        source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh
        ;;
    2)
        echo "Load SVFIT..."
        LD_LIBRARY_PATH=/usr/local/lib
        CPPYY_BACKEND_LIBRARY=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/libcppyy_backend3_9.so
        source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/bin/thisroot.sh
        export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/
        python -c "import TauAnalysis.ClassicSVfit.wrapper.pybind_wrapper" && echo "Module is working" || echo "Module import failed"
        ;;
    3)
        echo "Load CP-TOOLS..."
        source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh
        ;;
    4)
        echo "Activating OPTION4..."
        ;;
    *)
        echo "Invalid option selected: '$option'" >&2
        return 1
        ;;
esac
