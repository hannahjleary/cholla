#!/bin/bash

#-- This script needs to be sourced in the terminal from the login node, e.g.
#   source ./setup.crc.a100.sh

module purge
export MODULEPATH=/software/rhel9/spack-2026/modules/Core
module load gcc/13.4.0 openmpi/4.1.8-cuda hdf5/1.14.6