/*!
* \file hlld_cuda_tests.cpp
* \author Robert 'Bob' Caddy (rvc@pitt.edu)
* \brief Test the code units within hlld_cuda.cu
*
*/

// STL Includes
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <valarray>

// External Includes
#include <gtest/gtest.h>    // Include GoogleTest and related libraries/headers

// Local Includes
#include "../global/global_cuda.h"
#include "../utils/gpu.hpp"
#include "../utils/testing_utilities.h"
#include "../utils/mhd_utilities.h"
#include "../riemann_solvers/hlld_cuda.h"   // Include code to test

#if defined(CUDA) && defined(HLLD)
    // =========================================================================
    // Integration tests for the entire HLLD solver. Unit tests are below
    // =========================================================================

    // =========================================================================
    /*!
    * \brief Test fixture for simple testing of the HLLD Riemann Solver.
    Effectively takes the left state, right state, fiducial fluxes, and
    custom user output then performs all the required running and testing
    *
    */
    class tMHDCalculateHLLDFluxesCUDA : public ::testing::Test
    {
    protected:
        // =====================================================================
        /*!
        * \brief Compute and return the HLLD fluxes
        *
        * \param[in] leftState The state on the left side in conserved
        * variables. In order the elements are: density, x-momentum,
        * y-momentum, z-momentum, energy, passive scalars, x-magnetic field,
        * y-magnetic field, z-magnetic field.
        * \param[in] rightState The state on the right side in conserved
        * variables. In order the elements are: density, x-momentum,
        * y-momentum, z-momentum, energy, passive scalars, x-magnetic field,
        * y-magnetic field, z-magnetic field.
        * \param[in] gamma The adiabatic index
        * \param[in] direction Which plane the interface is. 0 = plane normal to
        * X, 1 = plane normal to Y, 2 = plane normal to Z. Defaults to 0.
        * \return std::vector<double>
        */
        std::vector<Real> computeFluxes(std::vector<Real> stateLeft,
                                        std::vector<Real> stateRight,
                                        Real const &gamma,
                                        int const &direction=0)
        {

            // Rearrange X, Y, and Z values if a different direction is chosen
            // besides default
            stateLeft  = _cycleXYZ(stateLeft, direction);
            stateRight = _cycleXYZ(stateRight, direction);

            // Simulation Paramters
            int const nx        = 1;  // Number of cells in the x-direction?
            int const ny        = 1;  // Number of cells in the y-direction?
            int const nz        = 1;  // Number of cells in the z-direction?
            int const nGhost    = 0;  // Isn't actually used it appears
            int nFields         = 8;  // Total number of conserved fields
            #ifdef  SCALAR
                nFields += NSCALARS;
            #endif  // SCALAR
            #ifdef  DE
                nFields++;
            #endif  //DE

            // Launch Parameters
            dim3 const dimGrid (1,1,1);  // How many blocks in the grid
            dim3 const dimBlock(1,1,1);  // How many threads per block

            // Create the std::vector to store the fluxes and declare the device
            // pointers
            std::vector<Real> testFlux(nFields);
            Real *devConservedLeft;
            Real *devConservedRight;
            Real *devTestFlux;

            // Allocate device arrays and copy data
            CudaSafeCall(cudaMalloc(&devConservedLeft,  nFields*sizeof(Real)));
            CudaSafeCall(cudaMalloc(&devConservedRight, nFields*sizeof(Real)));
            CudaSafeCall(cudaMalloc(&devTestFlux,       nFields*sizeof(Real)));

            CudaSafeCall(cudaMemcpy(devConservedLeft,
                         stateLeft.data(),
                         nFields*sizeof(Real),
                         cudaMemcpyHostToDevice));
            CudaSafeCall(cudaMemcpy(devConservedRight,
                         stateRight.data(),
                         nFields*sizeof(Real),
                         cudaMemcpyHostToDevice));

            // Run kernel
            hipLaunchKernelGGL(Calculate_HLLD_Fluxes_CUDA,
                               dimGrid,
                               dimBlock,
                               0,
                               0,
                               devConservedLeft,   // the "left" interface
                               devConservedRight,  // the "right" interface
                               devTestFlux,
                               nx,
                               ny,
                               nz,
                               nGhost,
                               gamma,
                               direction,
                               nFields);

            CudaCheckError();
            CudaSafeCall(cudaMemcpy(testFlux.data(),
                                    devTestFlux,
                                    nFields*sizeof(Real),
                                    cudaMemcpyDeviceToHost));

            // Make sure to sync with the device so we have the results
            cudaDeviceSynchronize();
            CudaCheckError();

            return testFlux;
        }
        // =====================================================================

        // =====================================================================
        /*!
        * \brief Check if the fluxes are correct
        *
        * \param[in] fiducialFlux The fiducial flux in conserved variables. In
        * order the elements are: density, x-momentum,
        * y-momentum, z-momentum, energy, passive scalars, x-magnetic field,
        * y-magnetic field, z-magnetic field.
        * \param[in] scalarFlux The fiducial flux in the passive scalars
        * \param[in] thermalEnergyFlux The fiducial flux in the dual energy
        * thermal energy
        * \param[in] testFlux The test flux in conserved variables. In order the
        * elements are: density, x-momentum,
        * y-momentum, z-momentum, energy, passive scalars, x-magnetic field,
        * y-magnetic field, z-magnetic field.
        * \param[in] customOutput Any custom output the user would like to
        * print. It will print after the default GTest output but before the
        * values that failed are printed
        * \param[in] direction Which plane the interface is. 0 = plane normal to
        * X, 1 = plane normal to Y, 2 = plane normal to Z. Defaults to 0.
        */
        void checkResults(std::vector<Real> fiducialFlux,
                          std::vector<Real> scalarFlux,
                          Real thermalEnergyFlux,
                          std::vector<Real> const &testFlux,
                          std::string const &customOutput = "",
                          int const &direction=0)
        {
            // Field names
            std::vector<std::string> fieldNames{"Densities",
                                                "X Momentum",
                                                "Y Momentum",
                                                "Z Momentum",
                                                "Energies",
                                                "X Magnetic Field",
                                                "Y Magnetic Field",
                                                "Z Magnetic Field"};
            #ifdef  DE
                fieldNames.push_back("Thermal energy (dual energy)");
                fiducialFlux.push_back(thermalEnergyFlux);
            #endif  //DE
            #ifdef  SCALAR
                std::vector<std::string> scalarNames{"Scalar 1", "Scalar 2", "Scalar 3"};
                fieldNames.insert(fieldNames.begin()+5,
                                  scalarNames.begin(),
                                  scalarNames.begin() + NSCALARS);

                fiducialFlux.insert(fiducialFlux.begin()+5,
                                    scalarFlux.begin(),
                                    scalarFlux.begin() + NSCALARS);
            #endif  //SCALAR

            // Rearrange X, Y, and Z values if a different direction is chosen
            // besides default
            fiducialFlux = _cycleXYZ(fiducialFlux, direction);

            ASSERT_TRUE(    (fiducialFlux.size() == testFlux.size())
                        and (fiducialFlux.size() == fieldNames.size()))
                << "The fiducial flux, test flux, and field name vectors are not all the same length" << std::endl
                << "fiducialFlux.size() = " << fiducialFlux.size() << std::endl
                << "testFlux.size() = "     << testFlux.size()     << std::endl
                << "fieldNames.size() = "   << fieldNames.size()   << std::endl;

            // Check for equality
            for (size_t i = 0; i < fieldNames.size(); i++)
            {
                // Check for equality and if not equal return difference
                double absoluteDiff;
                int64_t ulpsDiff;

                bool areEqual = testingUtilities::nearlyEqualDbl(fiducialFlux[i],
                                                                 testFlux[i],
                                                                 absoluteDiff,
                                                                 ulpsDiff);
                EXPECT_TRUE(areEqual)
                    << std::endl << customOutput << std::endl
                    << "There's a difference in "      << fieldNames[i]   << " Flux" << std::endl
                    << "The direction is:       "      << direction       << " (0=X, 1=Y, 2=Z)" << std::endl
                    << "The fiducial value is:       " << fiducialFlux[i] << std::endl
                    << "The test value is:           " << testFlux[i]     << std::endl
                    << "The absolute difference is:  " << absoluteDiff    << std::endl
                    << "The ULP difference is:       " << ulpsDiff        << std::endl;
            }
        }
        // =====================================================================

        // =====================================================================
        /*!
         * \brief Convert a vector of quantities in primitive variables  to
         * conserved variables
         *
         * \param[in] input The state in primitive variables. In order the
         * elements are: density, x-momentum,
         * y-momentum, z-momentum, energy, passive scalars, x-magnetic field,
         * y-magnetic field, z-magnetic field.
         * \return std::vector<Real> The state in conserved variables. In order
         * the elements are: density, x-momentum,
         * y-momentum, z-momentum, energy, passive scalars, x-magnetic field,
         * y-magnetic field, z-magnetic field.
         */
        std::vector<Real> primitive2Conserved(std::vector<Real> const &input,
                                              double const &gamma,
                                              std::vector<Real> const &primitiveScalars)
        {
            std::vector<Real> output(input.size());
            output.at(0) = input.at(0);  // Density
            output.at(1) = input.at(1) * input.at(0);  // X Velocity to momentum
            output.at(2) = input.at(2) * input.at(0);  // Y Velocity to momentum
            output.at(3) = input.at(3) * input.at(0);  // Z Velocity to momentum
            output.at(4) = mhdUtils::computeEnergy(input.at(4),
                                                   input.at(0),
                                                   input.at(1),
                                                   input.at(2),
                                                   input.at(3),
                                                   input.at(5),
                                                   input.at(6),
                                                   input.at(7),
                                                   gamma);  // Pressure to Energy
            output.at(5) = input.at(5);  // X Magnetic Field
            output.at(6) = input.at(6);  // Y Magnetic Field
            output.at(7) = input.at(7);  // Z Magnetic Field

            #ifdef SCALAR
                std::vector<Real> conservedScalar(primitiveScalars.size());
                std::transform(primitiveScalars.begin(),
                               primitiveScalars.end(),
                               conservedScalar.begin(),
                               [&](Real const &c){ return c*output.at(0); });
                output.insert(output.begin()+5,
                              conservedScalar.begin(),
                              conservedScalar.begin() + NSCALARS);
            #endif //SCALAR
            #ifdef  DE
                output.push_back(mhdUtils::computeThermalEnergy(output.at(4),
                                                                output.at(0),
                                                                output.at(1),
                                                                output.at(2),
                                                                output.at(3),
                                                                output.at(5 + NSCALARS),
                                                                output.at(6 + NSCALARS),
                                                                output.at(7 + NSCALARS),
                                                                gamma));
            #endif  //DE
            return output;
        }
        // =====================================================================

        // =====================================================================
        /*!
         * \brief On test start make sure that the number of NSCALARS is allowed
         *
         */
        void SetUp()
        {
            #ifdef  SCALAR
                ASSERT_LE(NSCALARS, 3) << "Only up to 3 passive scalars are currently supported in HLLD tests. NSCALARS = " << NSCALARS;
                ASSERT_GE(NSCALARS, 1) << "There must be at least 1 passive scalar to test with passive scalars. NSCALARS = " << NSCALARS;
            #endif  //SCALAR
        }
        // =====================================================================
    private:
        // =====================================================================
        /*!
         * \brief Cyclically permute the vector quantities in the list of
         * conserved variables so that the same interfaces and fluxes can be
         * used to test the HLLD solver in all 3 directions.
         *
         * \param[in,out] conservedVec The std::vector of conserved variables to
         * be cyclically permutated
         * \param[in] direction Which plane the interface is. 0 = plane normal
         * to X, 1 = plane normal to Y, 2 = plane normal to Z
         *
         * \return std::vector<Real> The cyclically permutated list of conserved
         * variables
         */
        std::vector<Real> inline _cycleXYZ(std::vector<Real> conservedVec,
                                           int const &direction)
        {
            switch (direction)
            {
            case 0:  // Plane normal to X. Default case, do nothing
                ;
                break;
            case 1:  // Plane normal to Y
            case 2:  // Plane normal to Z
                // Fall through for both Y and Z normal planes
                {
                    size_t shift = 3 - direction;
                    auto momentumBegin = conservedVec.begin()+1;
                    auto magneticBegin = conservedVec.begin()+5;
                    #ifdef  SCALAR
                        magneticBegin += NSCALARS;
                    #endif  //SCALAR

                    std::rotate(momentumBegin, momentumBegin+shift, momentumBegin+3);
                    std::rotate(magneticBegin, magneticBegin+shift, magneticBegin+3);
                }
                break;
            default:
                throw std::invalid_argument(("Invalid Value of `direction`"
                    " passed to `_cycleXYZ`. Value passed was "
                     + std::to_string(direction) + ", should be 0, 1, or 2."));
                break;
            }
            return conservedVec;
        }
        // =====================================================================
    };
    // =========================================================================

    // =========================================================================
    /*!
    * \brief Test the HLLD Riemann Solver using various states and waves from
    * the Brio & Wu Shock tube
    *
    */
    TEST_F(tMHDCalculateHLLDFluxesCUDA,
           BrioAndWuShockTubeCorrectInputExpectCorrectOutput)
    {
        // Constant Values
        Real const gamma = 2.;
        Real const Vz    = 0.0;
        Real const Bx    = 0.75;
        Real const Bz    = 0.0;
        std::vector<Real> const primitiveScalar{1.1069975296, 2.2286185018, 3.3155141875};

        // States
        std::vector<Real> const                       // | Density | X-Velocity | Y-Velocity | Z-Velocity | Pressure | X-Magnetic Field | Y-Magnetic Field | Z-Magnetic Field | Adiabatic Index | Passive Scalars |
            leftICs                = primitive2Conserved({1.0,       0.0,         0.0,        Vz,           1.0,       Bx,                1.0     ,          Bz},               gamma,            primitiveScalar),
            leftFastRareLeftSide   = primitive2Conserved({0.978576,  0.038603,   -0.011074,   Vz,           0.957621,  Bx,                0.970288,          Bz},               gamma,            primitiveScalar),
            leftFastRareRightSide  = primitive2Conserved({0.671655,  0.647082,   -0.238291,   Vz,           0.451115,  Bx,                0.578240,          Bz},               gamma,            primitiveScalar),
            compoundLeftSide       = primitive2Conserved({0.814306,  0.506792,   -0.911794,   Vz,           0.706578,  Bx,               -0.108819,          Bz},               gamma,            primitiveScalar),
            compoundPeak           = primitive2Conserved({0.765841,  0.523701,   -1.383720,   Vz,           0.624742,  Bx,               -0.400787,          Bz},               gamma,            primitiveScalar),
            compoundRightSide      = primitive2Conserved({0.695211,  0.601089,   -1.583720,   Vz,           0.515237,  Bx,               -0.537027,          Bz},               gamma,            primitiveScalar),
            contactLeftSide        = primitive2Conserved({0.680453,  0.598922,   -1.584490,   Vz,           0.515856,  Bx,               -0.533616,          Bz},               gamma,            primitiveScalar),
            contactRightSide       = primitive2Conserved({0.231160,  0.599261,   -1.584820,   Vz,           0.516212,  Bx,               -0.533327,          Bz},               gamma,            primitiveScalar),
            slowShockLeftSide      = primitive2Conserved({0.153125,  0.086170,   -0.683303,   Vz,           0.191168,  Bx,               -0.850815,          Bz},               gamma,            primitiveScalar),
            slowShockRightSide     = primitive2Conserved({0.117046, -0.238196,   -0.165561,   Vz,           0.087684,  Bx,               -0.903407,          Bz},               gamma,            primitiveScalar),
            rightFastRareLeftSide  = primitive2Conserved({0.117358, -0.228756,   -0.158845,   Vz,           0.088148,  Bx,               -0.908335,          Bz},               gamma,            primitiveScalar),
            rightFastRareRightSide = primitive2Conserved({0.124894, -0.003132,   -0.002074,   Vz,           0.099830,  Bx,               -0.999018,          Bz},               gamma,            primitiveScalar),
            rightICs               = primitive2Conserved({0.128,     0.0,         0.0,        Vz,           0.1,       Bx,               -1.0,               Bz},               gamma,            primitiveScalar);

        for (size_t direction = 0; direction < 3; direction++)
        {
            // Initial Condition Checks
            {
                std::string const outputString {"Left State:  Left Brio & Wu state\n"
                                                "Right State: Left Brio & Wu state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0, 1.21875, -0.75, 0, 0, 0.0, 0, 0};
                std::vector<Real> const scalarFlux{0, 0, 0};
                Real thermalEnergyFlux = 0.0;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Brio & Wu state\n"
                                                "Right State: Right Brio & Wu state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0, 0.31874999999999998, 0.75, 0, 0, 0.0, 0, 0};
                std::vector<Real> const scalarFlux{0, 0, 0};
                Real thermalEnergyFlux = 0.0;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left Brio & Wu state\n"
                                                "Right State: Right Brio & Wu state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.20673357746080057, 0.4661897584603672, 0.061170028480309613, 0, 0.064707291981509041, 0.0, 1.0074980455427278, 0};
                std::vector<Real> const scalarFlux{0.22885355953447648, 0.46073027567244362, 0.6854281091039145};
                Real thermalEnergyFlux = 0.20673357746080046;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Brio & Wu state\n"
                                                "Right State: Left Brio & Wu state\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.20673357746080057, 0.4661897584603672, 0.061170028480309613, 0, -0.064707291981509041, 0.0, -1.0074980455427278, 0};
                std::vector<Real> const scalarFlux{-0.22885355953447648, -0.46073027567244362, -0.6854281091039145};
                Real thermalEnergyFlux = -0.20673357746080046;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }

            // Cross wave checks
            {
                std::string const outputString {"Left State:  Left of left fast rarefaction\n"
                                                "Right State: Right of left fast rarefaction\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.4253304970883941, 0.47729308161522394, -0.55321646324583107, 0, 0.92496835095531071, 0.0, 0.53128887284876058, 0};
                std::vector<Real> const scalarFlux{0.47083980954039228, 0.94789941519098619, 1.4101892974729979};
                Real thermalEnergyFlux = 0.41622256825457099;
                std::vector<Real> const testFluxes = computeFluxes(leftFastRareLeftSide,
                                                                   leftFastRareRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of left fast rarefaction\n"
                                                "Right State: Left of left fast rarefaction\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.070492123816403796, 1.2489600267034342, -0.71031457071286608, 0, 0.21008080091470105, 0.0, 0.058615131833681167, 0};
                std::vector<Real> const scalarFlux{0.078034606921016325, 0.15710005136841393, 0.23371763662029341};
                Real thermalEnergyFlux = 0.047345816580591255;
                std::vector<Real> const testFluxes = computeFluxes(leftFastRareRightSide,
                                                                   leftFastRareLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of compound wave\n"
                                                "Right State: Right of compound wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.4470171023231666, 0.60747660800918468, -0.20506357956052623, 0, 0.72655525704800772, 0.0, 0.76278089951123285, 0};
                std::vector<Real> const scalarFlux{0.4948468279606959, 0.99623058485843297, 1.482091544807598};
                Real thermalEnergyFlux = 0.38787931087981475;
                std::vector<Real> const testFluxes = computeFluxes(compoundLeftSide,
                                                                   compoundRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of compound wave\n"
                                                "Right State: Left of compound wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.38496850292724116, 0.66092864409611585, -0.3473204105316457, 0, 0.89888639514227009, 0.0, 0.71658566275120927, 0};
                std::vector<Real> const scalarFlux{0.42615918171426637, 0.85794792823389721, 1.2763685331959034};
                Real thermalEnergyFlux = 0.28530908823756074;
                std::vector<Real> const testFluxes = computeFluxes(compoundRightSide,
                                                                   compoundLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of Compound Wave\n"
                                                "Right State: Peak of Compound Wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.41864266180405574, 0.63505764056357727, -0.1991008813536404, 0, 0.73707474818824525, 0.0, 0.74058225030218761, 0};
                std::vector<Real> const scalarFlux{0.46343639240225803, 0.93299478173931882, 1.388015684704111};
                Real thermalEnergyFlux = 0.36325864563467081;
                std::vector<Real> const testFluxes = computeFluxes(compoundLeftSide,
                                                                   compoundPeak,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Peak of Compound Wave\n"
                                                "Right State: Left of Compound Wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.39520761138156862, 0.6390998385557225, -0.35132701297727598, 0, 0.89945171879176522, 0.0, 0.71026545717401468, 0};
                std::vector<Real> const scalarFlux{0.43749384947851333, 0.88076699477714815, 1.3103164425435772};
                Real thermalEnergyFlux = 0.32239432669410983;
                std::vector<Real> const testFluxes = computeFluxes(compoundPeak,
                                                                   compoundLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Peak of Compound Wave\n"
                                                "Right State: Right of Compound Wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.4285899590904928, 0.6079309920345296, -0.26055320217638239, 0, 0.75090757444649436, 0.0, 0.85591904930227747, 0};
                std::vector<Real> const scalarFlux{0.47444802592454061, 0.95516351251477749, 1.4209960899845735};
                Real thermalEnergyFlux = 0.34962629086469987;
                std::vector<Real> const testFluxes = computeFluxes(compoundPeak,
                                                                   compoundRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of Compound Wave\n"
                                                "Right State: Peak of Compound Wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.39102247793946454, 0.65467021266207581, -0.25227691377588229, 0, 0.76271525822813691, 0.0, 0.83594460438033491, 0};
                std::vector<Real> const scalarFlux{0.43286091709705776, 0.8714399289555731, 1.2964405732397004};
                Real thermalEnergyFlux = 0.28979582956267347;
                std::vector<Real> const testFluxes = computeFluxes(compoundRightSide,
                                                                   compoundPeak,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of contact discontinuity\n"
                                                "Right State: Right of contact discontinuity\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.40753761783585118, 0.62106392255463172, -0.2455554035355339, 0, 0.73906344777217226, 0.0, 0.8687394222350926, 0};
                std::vector<Real> const scalarFlux{0.45114313616335622, 0.90824587528847567, 1.3511967538747176};
                Real thermalEnergyFlux = 0.30895701155896288;
                std::vector<Real> const testFluxes = computeFluxes(contactLeftSide,
                                                                   contactRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of contact discontinuity\n"
                                                "Right State: Left of contact discontinuity\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.13849588572126192, 0.46025037934770729, 0.18052412687974539, 0, 0.35385590617992224, 0.0, 0.86909622543144227, 0};
                std::vector<Real> const scalarFlux{0.15331460335320088, 0.30865449334158279, 0.45918507401922254};
                Real thermalEnergyFlux = 0.30928031735570188;
                std::vector<Real> const testFluxes = computeFluxes(contactRightSide,
                                                                   contactLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Slow shock left side\n"
                                                "Right State: Slow shock right side\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{3.5274134848883865e-05, 0.32304849716274459, 0.60579784881286636, 0, -0.32813070621836449, 0.0, 0.40636483121437972, 0};
                std::vector<Real> const scalarFlux{3.9048380136491711e-05, 7.8612589559210735e-05, 0.00011695189454326261};
                Real thermalEnergyFlux = 4.4037784886918126e-05;
                std::vector<Real> const testFluxes = computeFluxes(slowShockLeftSide,
                                                                   slowShockRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Slow shock right side\n"
                                                "Right State: Slow shock left side\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.016514307834939734, 0.16452009375678914, 0.71622171077118635, 0, -0.37262428139914472, 0.0, 0.37204015363322052, 0};
                std::vector<Real> const scalarFlux{-0.018281297976332211, -0.036804091985367396, -0.054753421923485097};
                Real thermalEnergyFlux = -0.020617189878790236;
                std::vector<Real> const testFluxes = computeFluxes(slowShockRightSide,
                                                                   slowShockLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right fast rarefaction left side\n"
                                                "Right State: Right fast rarefaction right side\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.026222824218991747, 0.22254903570732654, 0.68544334213642255, 0, -0.33339172106895454, 0.0, 0.32319665359522443, 0};
                std::vector<Real> const scalarFlux{-0.029028601629558917, -0.058440671223894146, -0.086942145734385745};
                Real thermalEnergyFlux = -0.020960370728633469;
                std::vector<Real> const testFluxes = computeFluxes(rightFastRareLeftSide,
                                                                   rightFastRareRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right fast rarefaction right side\n"
                                                "Right State: Right fast rarefaction left side\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.001088867226159973, 0.32035322820305906, 0.74922357263343131, 0, -0.0099746892805345766, 0.0, 0.0082135595470345102, 0};
                std::vector<Real> const scalarFlux{-0.0012053733294214947, -0.0024266696462237609, -0.0036101547366371614};
                Real thermalEnergyFlux = -0.00081785194236053073;
                std::vector<Real> const testFluxes = computeFluxes(rightFastRareRightSide,
                                                                   rightFastRareLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
    * \brief Test the HLLD Riemann Solver using various states and waves from
    * the Dai & Woodward Shock tube
    *
    */
    TEST_F(tMHDCalculateHLLDFluxesCUDA,
           DaiAndWoodwardShockTubeCorrectInputExpectCorrectOutput)
    {
        // Constant Values
        Real const gamma = 5./3.;
        Real const coef = 1. / (std::sqrt(4. * M_PI));
        Real const Bx = 4. * coef;
        std::vector<Real> const primitiveScalar{1.1069975296, 2.2286185018, 3.3155141875};

        // States
        std::vector<Real> const                       // | Density | X-Velocity | Y-Velocity | Z-Velocity | Pressure | X-Magnetic Field | Y-Magnetic Field | Z-Magnetic Field | Adiabatic Index | Passive Scalars |
            leftICs                 = primitive2Conserved({1.08,     0.0,         0.0,         0.0,         1.0,       Bx,                3.6*coef,          2*coef},           gamma,            primitiveScalar),
            leftFastShockLeftSide   = primitive2Conserved({1.09406,  1.176560,    0.021003,    0.506113,    0.970815,  1.12838,           1.105355,          0.614087},         gamma,            primitiveScalar),
            leftFastShockRightSide  = primitive2Conserved({1.40577,  0.693255,    0.210562,    0.611423,    1.494290,  1.12838,           1.457700,          0.809831},         gamma,            primitiveScalar),
            leftRotationLeftSide    = primitive2Conserved({1.40086,  0.687774,    0.215124,    0.609161,    1.485660,  1.12838,           1.458735,          0.789960},         gamma,            primitiveScalar),
            leftRotationRightSide   = primitive2Conserved({1.40119,  0.687504,    0.330268,    0.334140,    1.486570,  1.12838,           1.588975,          0.475782},         gamma,            primitiveScalar),
            leftSlowShockLeftSide   = primitive2Conserved({1.40519,  0.685492,    0.326265,    0.333664,    1.493710,  1.12838,           1.575785,          0.472390},         gamma,            primitiveScalar),
            leftSlowShockRightSide  = primitive2Conserved({1.66488,  0.578545,    0.050746,    0.250260,    1.984720,  1.12838,           1.344490,          0.402407},         gamma,            primitiveScalar),
            contactLeftSide         = primitive2Conserved({1.65220,  0.578296,    0.049683,    0.249962,    1.981250,  1.12838,           1.346155,          0.402868},         gamma,            primitiveScalar),
            contactRightSide        = primitive2Conserved({1.49279,  0.578276,    0.049650,    0.249924,    1.981160,  1.12838,           1.346180,          0.402897},         gamma,            primitiveScalar),
            rightSlowShockLeftSide  = primitive2Conserved({1.48581,  0.573195,    0.035338,    0.245592,    1.956320,  1.12838,           1.370395,          0.410220},         gamma,            primitiveScalar),
            rightSlowShockRightSide = primitive2Conserved({1.23813,  0.450361,   -0.275532,    0.151746,    1.439000,  1.12838,           1.609775,          0.482762},         gamma,            primitiveScalar),
            rightRotationLeftSide   = primitive2Conserved({1.23762,  0.450102,   -0.274410,    0.145585,    1.437950,  1.12838,           1.606945,          0.493879},         gamma,            primitiveScalar),
            rightRotationRightSide  = primitive2Conserved({1.23747,  0.449993,   -0.180766,   -0.090238,    1.437350,  1.12838,           1.503855,          0.752090},         gamma,            primitiveScalar),
            rightFastShockLeftSide  = primitive2Conserved({1.22305,  0.424403,   -0.171402,   -0.085701,    1.409660,  1.12838,           1.447730,          0.723864},         gamma,            primitiveScalar),
            rightFastShockRightSide = primitive2Conserved({1.00006,  0.000121,   -0.000057,   -0.000028,    1.000100,  1.12838,           1.128435,          0.564217},         gamma,            primitiveScalar),
            rightICs                = primitive2Conserved({1.0,      0.0,         0.0,         1.0,         0.2,       Bx,                4*coef,            2*coef},           gamma,            primitiveScalar);

        for (size_t direction = 0; direction < 3; direction++)
        {
            // Initial Condition Checks
            {
                std::string const outputString {"Left State:  Left Dai & Woodward state\n"
                                                "Right State: Left Dai & Woodward state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0, 1.0381971863420549, -1.1459155902616465, -0.63661977236758127, 0, 0.0, 0, -1.1102230246251565e-16};
                std::vector<Real> const scalarFlux{0,0,0};
                Real thermalEnergyFlux = 0.0;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Dai & Woodward state\n"
                                                "Right State: Right Dai & Woodward state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0, 0.35915494309189522, -1.2732395447351625, -0.63661977236758127, -0.63661977236758172, 0.0, 2.2204460492503131e-16, -1.1283791670955123};
                std::vector<Real> const scalarFlux{0,0,0};
                Real thermalEnergyFlux = 0.0;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left Dai & Woodward state\n"
                                                "Right State: Right Dai & Woodward state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.17354924587196074, 0.71614983677687327, -1.1940929411768009, -1.1194725181819352, -0.11432087006939984, 0.0, 0.056156000248263505, -0.42800560867873094};
                std::vector<Real> const scalarFlux{0.19211858644420357, 0.38677506032368902, 0.57540498691841158};
                Real thermalEnergyFlux = 0.24104061926661174;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Dai & Woodward state\n"
                                                "Right State: Left Dai & Woodward state\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.17354924587196074, 0.71614983677687327, -1.1940929411768009, -0.14549552299758384, -0.47242308031148195, 0.0, -0.056156000248263505, -0.55262526758377528};
                std::vector<Real> const scalarFlux{-0.19211858644420357, -0.38677506032368902, -0.57540498691841158};
                Real thermalEnergyFlux = -0.24104061926661174;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }

            // Cross wave checks
            {
                std::string const outputString {"Left State:  Left of left fast shock\n"
                                                "Right State: Right of left fast shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.96813688187727132, 3.0871217875403394, -1.4687093290523414, -0.33726008721080036, 4.2986213406773457, 0.0, 0.84684181393860269, -0.087452560407274671};
                std::vector<Real> const scalarFlux{1.0717251365527865, 2.157607767226648, 3.2098715673061045};
                Real thermalEnergyFlux = 1.2886155333980993;
                std::vector<Real> const testFluxes = computeFluxes(leftFastShockLeftSide,
                                                                   leftFastShockRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of left fast shock\n"
                                                "Right State: Left of left fast shock\n"
                                                "HLLD State: Left Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{1.3053938862274184, 2.4685129176021858, -1.181892850065283, -0.011160487372167127, 5.1797404608257249, 0.0, 1.1889903073770265, 0.10262704114294516};
                std::vector<Real> const scalarFlux{1.4450678072086958, 2.9092249669830292, 4.3280519500627666};
                Real thermalEnergyFlux = 2.081389946702628;
                std::vector<Real> const testFluxes = computeFluxes(leftFastShockRightSide,
                                                                   leftFastShockLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of left rotation/Alfven wave\n"
                                                "Right State: Right of left rotation/Alfven wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.96326128304298586, 2.8879592118317445, -1.4808188010794987, -0.20403672861184916, 4.014027751838869, 0.0, 0.7248753989305099, -0.059178137562467162};
                std::vector<Real> const scalarFlux{1.0663278606879119, 2.1467419174572049, 3.1937064501984724};
                Real thermalEnergyFlux = 1.5323573637968553;
                std::vector<Real> const testFluxes = computeFluxes(leftRotationLeftSide,
                                                                   leftRotationRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of left rotation/Alfven wave\n"
                                                "Right State: Left of left rotation/Alfven wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.96353754504060063, 2.8875487093397085, -1.4327309336053695, -0.31541343522923493, 3.9739842521208342, 0.0, 0.75541746728406312, -0.13479771672887678};
                std::vector<Real> const scalarFlux{1.0666336820367937, 2.1473576000564334, 3.1946224007710313};
                Real thermalEnergyFlux = 1.5333744977458499;
                std::vector<Real> const testFluxes = computeFluxes(leftRotationRightSide,
                                                                   leftRotationLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of left slow shock\n"
                                                "Right State: Right of left slow shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.88716095730727451, 2.9828594399125663, -1.417062582518549, -0.21524331343191233, 3.863474778369334, 0.0, 0.71242370728996041, -0.05229712416644372};
                std::vector<Real> const scalarFlux{0.98208498809672407, 1.9771433235295921, 2.9413947405483505};
                Real thermalEnergyFlux = 1.4145715457049737;
                std::vector<Real> const testFluxes = computeFluxes(leftSlowShockLeftSide,
                                                                   leftSlowShockRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of left slow shock\n"
                                                "Right State: Left of left slow shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{1.042385440439527, 2.7732383399777376, -1.5199872074603551, -0.21019362664841068, 4.1322001036232585, 0.0, 0.72170937317481543, -0.049474715634396704};
                std::vector<Real> const scalarFlux{1.1539181074575644, 2.323079478570472, 3.4560437166206879};
                Real thermalEnergyFlux = 1.8639570701934713;
                std::vector<Real> const testFluxes = computeFluxes(leftSlowShockRightSide,
                                                                   leftSlowShockLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of contact discontinuity\n"
                                                "Right State: Right of contact discontinuity\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.95545795601418737, 2.8843900822429749, -1.4715039715239722, -0.21575736014726318, 4.0078718055059257, 0.0, 0.72241353110189066, -0.049073560388753337};
                std::vector<Real> const scalarFlux{1.0576895969443709, 2.1293512784652289, 3.1678344087247892};
                Real thermalEnergyFlux = 1.7186185770667382;
                std::vector<Real> const testFluxes = computeFluxes(contactLeftSide,
                                                                   contactRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of contact discontinuity\n"
                                                "Right State: Left of contact discontinuity\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.86324813554422819, 2.8309913324581251, -1.4761428591480787, -0.23887765947428419, 3.9892942559102793, 0.0, 0.72244123046603836, -0.049025527032060034};
                std::vector<Real> const scalarFlux{0.95561355347926669, 1.9238507665182214, 2.8621114407298114};
                Real thermalEnergyFlux = 1.7184928987481187;
                std::vector<Real> const testFluxes = computeFluxes(contactRightSide,
                                                                   contactLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of right slow shock\n"
                                                "Right State: Right of right slow shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.81125524370350677, 2.901639500435365, -1.5141545346789429, -0.262600896007809, 3.8479660419540087, 0.0, 0.7218977970017596, -0.049091614519593846};
                std::vector<Real> const scalarFlux{0.89805755065482806, 1.8079784457999033, 2.6897282701827465};
                Real thermalEnergyFlux = 1.6022319728249694;
                std::vector<Real> const testFluxes = computeFluxes(rightSlowShockLeftSide,
                                                                   rightSlowShockRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of right slow shock\n"
                                                "Right State: Left of right slow shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.60157947557836688, 2.3888357198399746, -1.9910500022202977, -0.45610948442354332, 3.5359430988850069, 0.0, 1.0670963294022622, 0.05554893654378229};
                std::vector<Real> const scalarFlux{0.66594699332331575, 1.3406911495770899, 1.994545286188885};
                Real thermalEnergyFlux = 1.0487665253534804;
                std::vector<Real> const testFluxes = computeFluxes(rightSlowShockRightSide,
                                                                   rightSlowShockLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of right rotation/Alfven wave\n"
                                                "Right State: Right of right rotation/Alfven wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.55701691287884714, 2.4652223621237814, -1.9664615862227277, -0.47490477894092042, 3.3900659850690529, 0.0, 1.0325648885587542, 0.059165409025635551};
                std::vector<Real> const scalarFlux{0.61661634650230224, 1.2413781978573175, 1.8467974773272691};
                Real thermalEnergyFlux = 0.9707694646266285;
                std::vector<Real> const testFluxes = computeFluxes(rightRotationLeftSide,
                                                                   rightRotationRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of right rotation/Alfven wave\n"
                                                "Right State: Left of right rotation/Alfven wave\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.55689116371132596, 2.4648517303940851, -1.7972202655166787, -0.90018282739798461, 3.3401033852664566, 0.0, 0.88105841856465605, 0.43911718823267476};
                std::vector<Real> const scalarFlux{0.61647714248450702, 1.2410979509359938, 1.8463805541782863};
                Real thermalEnergyFlux = 0.9702629326292449;
                std::vector<Real> const testFluxes = computeFluxes(rightRotationRightSide,
                                                                   rightRotationLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of right fast shock\n"
                                                "Right State: Right of right fast shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.48777637414577313, 2.3709438477809708, -1.7282900552525988, -0.86414423547773778, 2.8885015704245069, 0.0, 0.77133731061645838, 0.38566794697432505};
                std::vector<Real> const scalarFlux{0.53996724117661621, 1.0870674521621893, 1.6172294888076189};
                Real thermalEnergyFlux = 0.84330016382608752;
                std::vector<Real> const testFluxes = computeFluxes(rightFastShockLeftSide,
                                                                   rightFastShockRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of right fast shock\n"
                                                "Right State: Left of right fast shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.040639426423817904, 1.0717156491947966, -1.2612066401572222, -0.63060225433149875, 0.15803727234007203, 0.0, 0.042555541396817498, 0.021277678888288909};
                std::vector<Real> const scalarFlux{0.044987744655527385, 0.090569777630660403, 0.13474059488003065};
                Real thermalEnergyFlux = 0.060961577855018087;
                std::vector<Real> const testFluxes = computeFluxes(rightFastShockRightSide,
                                                                   rightFastShockLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
    * \brief Test the HLLD Riemann Solver using various states and waves from
    * the Ryu & Jones 4d Shock tube
    *
    */
    TEST_F(tMHDCalculateHLLDFluxesCUDA,
           RyuAndJones4dShockTubeCorrectInputExpectCorrectOutput)
    {
        // Constant Values
        Real const gamma = 5./3.;
        Real const Bx = 0.7;
        std::vector<Real> const primitiveScalar{1.1069975296, 2.2286185018, 3.3155141875};

        // States
        std::vector<Real> const                           // | Density | X-Velocity | Y-Velocity |  Z-Velocity |  Pressure | X-Magnetic Field | Y-Magnetic Field | Z-Magnetic Field | Adiabatic Index | Passive Scalars |
            leftICs                    = primitive2Conserved({1.0,       0.0,         0.0,          0.0,          1.0,       Bx,                0.0,               0.0},              gamma,            primitiveScalar),
            hydroRareLeftSide          = primitive2Conserved({0.990414,  0.012415,    1.458910e-58, 6.294360e-59, 0.984076,  Bx,                1.252355e-57,      5.366795e-58},     gamma,            primitiveScalar),
            hydroRareRightSide         = primitive2Conserved({0.939477,  0.079800,    1.557120e-41, 7.505190e-42, 0.901182,  Bx,                1.823624e-40,      8.712177e-41},     gamma,            primitiveScalar),
            switchOnSlowShockLeftSide  = primitive2Conserved({0.939863,  0.079142,    1.415730e-02, 7.134030e-03, 0.901820,  Bx,                2.519650e-02,      1.290082e-02},     gamma,            primitiveScalar),
            switchOnSlowShockRightSide = primitive2Conserved({0.651753,  0.322362,    8.070540e-01, 4.425110e-01, 0.490103,  Bx,                6.598380e-01,      3.618000e-01},     gamma,            primitiveScalar),
            contactLeftSide            = primitive2Conserved({0.648553,  0.322525,    8.072970e-01, 4.426950e-01, 0.489951,  Bx,                6.599295e-01,      3.618910e-01},     gamma,            primitiveScalar),
            contactRightSide           = primitive2Conserved({0.489933,  0.322518,    8.073090e-01, 4.426960e-01, 0.489980,  Bx,                6.599195e-01,      3.618850e-01},     gamma,            primitiveScalar),
            slowShockLeftSide          = primitive2Conserved({0.496478,  0.308418,    8.060830e-01, 4.420150e-01, 0.489823,  Bx,                6.686695e-01,      3.666915e-01},     gamma,            primitiveScalar),
            slowShockRightSide         = primitive2Conserved({0.298260, -0.016740,    2.372870e-01, 1.287780e-01, 0.198864,  Bx,                8.662095e-01,      4.757390e-01},     gamma,            primitiveScalar),
            rotationLeftSide           = primitive2Conserved({0.298001, -0.017358,    2.364790e-01, 1.278540e-01, 0.198448,  Bx,                8.669425e-01,      4.750845e-01},     gamma,            primitiveScalar),
            rotationRightSide          = primitive2Conserved({0.297673, -0.018657,    1.059540e-02, 9.996860e-01, 0.197421,  Bx,                9.891580e-01,      1.024949e-04},     gamma,            primitiveScalar),
            fastRareLeftSide           = primitive2Conserved({0.297504, -0.020018,    1.137420e-02, 1.000000e+00, 0.197234,  Bx,                9.883860e-01, -    4.981931e-17},     gamma,            primitiveScalar),
            fastRareRightSide          = primitive2Conserved({0.299996, -0.000033,    1.855120e-05, 1.000000e+00, 0.199995,  Bx,                9.999865e-01,      1.737190e-16},     gamma,            primitiveScalar),
            rightICs                   = primitive2Conserved({0.3,       0.0,         0.0,          1.0,          0.2,       Bx,                1.0,               0.0},              gamma,            primitiveScalar);

        for (size_t direction = 0; direction < 3; direction++)
        {
            // Initial Condition Checks
            {
                std::string const outputString {"Left State:  Left Ryu & Jones 4d state\n"
                                                "Right State: Left Ryu & Jones 4d state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0, 0.75499999999999989, 0, 0, 2.2204460492503131e-16, 0.0, 0, 0};
                std::vector<Real> const scalarFlux{0,0,0};
                Real thermalEnergyFlux = 0.0;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Ryu & Jones 4d state\n"
                                                "Right State: Right Ryu & Jones 4d state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-5.5511151231257827e-17, 0.45500000000000013, -0.69999999999999996, -5.5511151231257827e-17, 0, 0.0, 0, -0.69999999999999996};
                std::vector<Real> const scalarFlux{-6.1450707278254418e-17, -1.2371317869019906e-16, -1.8404800947169341e-16};
                Real thermalEnergyFlux = 0.0;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left Ryu & Jones 4d state\n"
                                                "Right State: Right Ryu & Jones 4d state\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.092428729855986602, 0.53311593977445149, -0.39622049648437296, -0.21566989083797167, -0.13287876964320211, 0.0, -0.40407579574102892, -0.21994567048141428};
                std::vector<Real> const scalarFlux{0.10231837561464294, 0.20598837745492582, 0.30644876517012837};
                Real thermalEnergyFlux = 0.13864309478397996;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Ryu & Jones 4d state\n"
                                                "Right State: Left Ryu & Jones 4d state\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.092428729855986602, 0.53311593977445149, -0.39622049648437296, 0.21566989083797167, 0.13287876964320211, 0.0, 0.40407579574102892, -0.21994567048141428};
                std::vector<Real> const scalarFlux{-0.10231837561464294, -0.20598837745492582, -0.30644876517012837};
                Real thermalEnergyFlux = -0.13864309478397996;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }

            // Cross wave checks
            {
                std::string const outputString {"Left State:  Left side of pure hydrodynamic rarefaction\n"
                                                "Right State: Right side of pure hydrodynamic rarefaction\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.074035256375659553, 0.66054553664209648, -6.1597070943493028e-41, -2.9447391900433873e-41, 0.1776649658235645, 0.0, -6.3466063324344113e-41, -3.0340891384335242e-41};
                std::vector<Real> const scalarFlux{0.081956845911157775, 0.16499634214430131, 0.24546494288869905};
                Real thermalEnergyFlux = 0.11034221894046368;
                std::vector<Real> const testFluxes = computeFluxes(hydroRareLeftSide,
                                                                   hydroRareRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right side of pure hydrodynamic rarefaction\n"
                                                "Right State: Left side of pure hydrodynamic rarefaction\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.013336890338886076, 0.74071279157971992, -6.1745213352160876e-41, -2.9474651270630147e-41, 0.033152482405470307, 0.0, 6.2022392844946449e-41, 2.9606965476795895e-41};
                std::vector<Real> const scalarFlux{0.014763904657692993, 0.029722840565719184, 0.044218649135708464};
                Real thermalEnergyFlux = 0.019189877201961154;
                std::vector<Real> const testFluxes = computeFluxes(hydroRareRightSide,
                                                                   hydroRareLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of switch on slow shock\n"
                                                "Right State: Right of switch on slow shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.19734622040826083, 0.47855039640569758, -0.3392293209655618, -0.18588204716255491, 0.10695446263054809, 0.0, -0.3558357543098733, -0.19525093130352045};
                std::vector<Real> const scalarFlux{0.21846177846784187, 0.43980943806215089, 0.65430419361309078};
                Real thermalEnergyFlux = 0.2840373040888583;
                std::vector<Real> const testFluxes = computeFluxes(switchOnSlowShockLeftSide,
                                                                   switchOnSlowShockRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of switch on slow shock\n"
                                                "Right State: Left of switch on slow shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.097593254768855386, 0.76483698872352757, -0.02036438492698419, -0.010747481940703562, 0.25327551496496836, 0.0, -0.002520109973016129, -0.00088262199017708799};
                std::vector<Real> const scalarFlux{0.10803549193474633, 0.21749813322875222, 0.32357182079044206};
                Real thermalEnergyFlux = 0.1100817647375162;
                std::vector<Real> const testFluxes = computeFluxes(switchOnSlowShockRightSide,
                                                                   switchOnSlowShockLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of contact discontinuity\n"
                                                "Right State: Right of contact discontinuity\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.2091677440314007, 0.5956612619664029, -0.29309091669513981, -0.16072556008504282, 0.19220050968424285, 0.0, -0.35226977371803297, -0.19316940226499904};
                std::vector<Real> const scalarFlux{0.23154817591476573, 0.46615510432814616, 0.69349862290347741};
                Real thermalEnergyFlux = 0.23702444986592192;
                std::vector<Real> const testFluxes = computeFluxes(contactLeftSide,
                                                                   contactRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of contact discontinuity\n"
                                                "Right State: Left of contact discontinuity\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.15801775068597168, 0.57916072367837657, -0.33437339604094024, -0.18336617461176744, 0.16789791355547545, 0.0, -0.3522739911439669, -0.19317084712861482};
                std::vector<Real> const scalarFlux{0.17492525964231936, 0.35216128279157616, 0.52391009427617696};
                Real thermalEnergyFlux = 0.23704936434506069;
                std::vector<Real> const testFluxes = computeFluxes(contactRightSide,
                                                                   contactLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of slow shock\n"
                                                "Right State: Right of slow shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.11744487326715558, 0.66868230621718128, -0.35832022960458892, -0.19650694834641164, 0.057880816021092185, 0.0, -0.37198011453582402, -0.20397277844271294};
                std::vector<Real> const scalarFlux{0.13001118457092631, 0.26173981750473918, 0.38939014356639379};
                Real thermalEnergyFlux = 0.1738058891582446;
                std::vector<Real> const testFluxes = computeFluxes(slowShockLeftSide,
                                                                   slowShockRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of slow shock\n"
                                                "Right State: Left of slow shock\n"
                                                "HLLD State: Left Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.038440990187426027, 0.33776683678923869, -0.62583241538732792, -0.3437911783906169, -0.13471828103488348, 0.0, -0.15165427985881363, -0.082233932588833825};
                std::vector<Real> const scalarFlux{0.042554081172858457, 0.085670301959209896, 0.12745164834795927};
                Real thermalEnergyFlux = 0.038445630017261548;
                std::vector<Real> const testFluxes = computeFluxes(slowShockRightSide,
                                                                   slowShockLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of rotation/Alfven wave\n"
                                                "Right State: Right of rotation/Alfven wave\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.0052668366104996478, 0.44242247672452317, -0.60785196341731951, -0.33352435102145184, -0.21197843894720192, 0.0, -0.18030635192654354, -0.098381113757603278};
                std::vector<Real> const scalarFlux{-0.0058303751166299484, -0.011737769516117116, -0.017462271505355991};
                Real thermalEnergyFlux = -0.0052395622905745485;
                std::vector<Real> const testFluxes = computeFluxes(rotationLeftSide,
                                                                   rotationRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of rotation/Alfven wave\n"
                                                "Right State: Left of rotation/Alfven wave\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.005459628948343731, 0.4415038084184626, -0.69273580053867279, -0.0051834737482743809, -0.037389286119015486, 0.0, -0.026148289294373184, -0.69914753968916865};
                std::vector<Real> const scalarFlux{-0.0060437957583491572, -0.012167430087241717, -0.018101477236719343};
                Real thermalEnergyFlux = -0.0054536013916442853;
                std::vector<Real> const testFluxes = computeFluxes(rotationRightSide,
                                                                   rotationLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left of fast rarefaction\n"
                                                "Right State: Right of fast rarefaction\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.0059354802028144249, 0.44075681881443612, -0.69194176811725872, -0.0059354802028144804, -0.040194357552219451, 0.0, -0.027710302430178135, -0.70000000000000007};
                std::vector<Real> const scalarFlux{-0.0065705619215052757, -0.013227920997059845, -0.019679168822056604};
                Real thermalEnergyFlux = -0.0059354109546219782;
                std::vector<Real> const testFluxes = computeFluxes(fastRareLeftSide,
                                                                   fastRareRightSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right of fast rarefaction\n"
                                                "Right State: Left of fast rarefaction\n"
                                                "HLLD State: Right Double Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-3.0171858819483255e-05, 0.45503057873272706, -0.69998654276213712, -3.0171858819427744e-05, -0.00014827469339251387, 0.0, -8.2898844654399895e-05, -0.69999999999999984};
                std::vector<Real> const scalarFlux{-3.340017317660794e-05, -6.7241562798797897e-05, -0.00010003522597924373};
                Real thermalEnergyFlux = -3.000421709818028e-05;
                std::vector<Real> const testFluxes = computeFluxes(fastRareRightSide,
                                                                   fastRareLeftSide,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
    * \brief Test the HLLD Riemann Solver using various states and waves from
    * the Einfeldt Strong Rarefaction (EFR)
    *
    */
    TEST_F(tMHDCalculateHLLDFluxesCUDA,
           EinfeldtStrongRarefactionCorrectInputExpectCorrectOutput)
    {
        // Constant Values
        Real const gamma = 5./3.;
        Real const V0 = 2.;
        Real const Vy = 0.0;
        Real const Vz = 0.0;
        Real const Bx = 0.0;
        Real const Bz = 0.0;

        std::vector<Real> const primitiveScalar{1.1069975296, 2.2286185018, 3.3155141875};

        // States
        std::vector<Real> const                       // | Density | X-Velocity | Y-Velocity | Z-Velocity | Pressure | X-Magnetic Field | Y-Magnetic Field | Z-Magnetic Field | Adiabatic Index | Passive Scalars |
            leftICs                = primitive2Conserved({1.0,      -V0,          Vy,          Vz,          0.45,      Bx,                0.5,               Bz},               gamma,            primitiveScalar),
            leftRarefactionCenter  = primitive2Conserved({0.368580, -1.180830,    Vy,          Vz,          0.111253,  Bx,                0.183044,          Bz},               gamma,            primitiveScalar),
            leftVxTurnOver         = primitive2Conserved({0.058814, -0.125475,    Vy,          Vz,          0.008819,  Bx,                0.029215,          Bz},               gamma,            primitiveScalar),
            midPoint               = primitive2Conserved({0.034658,  0.000778,    Vy,          Vz,          0.006776,  Bx,                0.017333,          Bz},               gamma,            primitiveScalar),
            rightVxTurnOver        = primitive2Conserved({0.062587,  0.152160,    Vy,          Vz,          0.009521,  Bx,                0.031576,          Bz},               gamma,            primitiveScalar),
            rightRarefactionCenter = primitive2Conserved({0.316485,  1.073560,    Vy,          Vz,          0.089875,  Bx,                0.159366,          Bz},               gamma,            primitiveScalar),
            rightICs               = primitive2Conserved({1.0,       V0,          Vy,          Vz,          0.45,      Bx,                0.5,               Bz},               gamma,            primitiveScalar);

        for (size_t direction = 0; direction < 3; direction++)
        {
            // Initial Condition Checks
            {
                std::string const outputString {"Left State:  Left Einfeldt Strong Rarefaction state\n"
                                                "Right State: Left Einfeldt Strong Rarefaction state\n"
                                                "HLLD State: Right"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-2, 4.5750000000000002, -0, -0, -6.75, 0.0, -1, -0};
                std::vector<Real> const scalarFlux{-2.2139950592000002, -4.4572370036000004, -6.6310283749999996};
                Real thermalEnergyFlux = -1.3499999999999996;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Einfeldt Strong Rarefaction state\n"
                                                "Right State: Right Einfeldt Strong Rarefaction state\n"
                                                "HLLD State: Left"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{2, 4.5750000000000002, 0, 0, 6.75, 0.0, 1, 0};
                std::vector<Real> const scalarFlux{2.2139950592000002, 4.4572370036000004, 6.6310283749999996};
                Real thermalEnergyFlux = 1.3499999999999996;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left Einfeldt Strong Rarefaction state\n"
                                                "Right State: Right Einfeldt Strong Rarefaction state\n"
                                                "HLLD State: Left Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0, -1.4249999999999998, -0, -0, 0, 0.0, 0, -0};
                std::vector<Real> const scalarFlux{0,0,0};
                Real thermalEnergyFlux = 0.0;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Einfeldt Strong Rarefaction state\n"
                                                "Right State: Left Einfeldt Strong Rarefaction state\n"
                                                "HLLD State: Left Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0, 10.574999999999999, 0, 0, 0, 0.0, 0, 0};
                std::vector<Real> const scalarFlux{0,0,0};
                Real thermalEnergyFlux = 0.0;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }

            // Intermediate state checks
            {
                std::string const outputString {"Left State:  Left Einfeldt Strong Rarefaction state\n"
                                                "Right State: Left rarefaction center\n"
                                                "HLLD State: Right"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.43523032140000006, 0.64193857338676208, -0, -0, -0.67142479846795033, 0.0, -0.21614384652000002, -0};
                std::vector<Real> const scalarFlux{-0.48179889059681413, -0.9699623468164007, -1.4430123054318851};
                Real thermalEnergyFlux = -0.19705631998499995;
                std::vector<Real> const testFluxes = computeFluxes(leftICs,
                                                                   leftRarefactionCenter,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left rarefaction center\n"
                                                "Right State: Left Einfeldt Strong Rarefaction state\n"
                                                "HLLD State: Right"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-2, 4.5750000000000002, -0, -0, -6.75, 0.0, -1, -0};
                std::vector<Real> const scalarFlux{-2.2139950592000002, -4.4572370036000004, -6.6310283749999996};
                Real thermalEnergyFlux = -1.3499999999999996;
                std::vector<Real> const testFluxes = computeFluxes(leftRarefactionCenter,
                                                                   leftICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left rarefaction center\n"
                                                "Right State: Left Vx turnover point\n"
                                                "HLLD State: Right Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.023176056428381629, -2.0437812714100764e-05, 0, 0, -0.00098843768795337005, 0.0, -0.011512369309265979, 0};
                std::vector<Real> const scalarFlux{-0.025655837212088663, -0.051650588155052128, -0.076840543898599858};
                Real thermalEnergyFlux = -0.0052127803322822184;
                std::vector<Real> const testFluxes = computeFluxes(leftRarefactionCenter,
                                                                   leftVxTurnOver,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left Vx turnover point\n"
                                                "Right State: Left rarefaction center\n"
                                                "HLLD State: Right Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.43613091609689758, 0.64135749005731213, 0, 0, -0.67086080671260462, 0.0, -0.21659109937066717, 0};
                std::vector<Real> const scalarFlux{-0.48279584670145054, -0.9719694288205295, -1.445998239926636};
                Real thermalEnergyFlux = -0.19746407621898149;
                std::vector<Real> const testFluxes = computeFluxes(leftVxTurnOver,
                                                                   leftRarefactionCenter,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Left Vx turnover point\n"
                                                "Right State: Midpoint\n"
                                                "HLLD State: Right Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.0011656375857387598, 0.0062355370788444902, 0, 0, -0.00055517615333601446, 0.0, -0.0005829533231464588, 0};
                std::vector<Real> const scalarFlux{-0.0012903579278217153, -0.0025977614899708843, -0.0038646879530001054};
                Real thermalEnergyFlux = -0.00034184143405415065;
                std::vector<Real> const testFluxes = computeFluxes(leftVxTurnOver,
                                                                   midPoint,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Midpoint\n"
                                                "Right State: Left Vx turnover point\n"
                                                "HLLD State: Right Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-0.0068097924351817191, 0.010501781004354172, 0, 0, -0.0027509360975397175, 0.0, -0.0033826654536986789, 0};
                std::vector<Real> const scalarFlux{-0.0075384234028349319, -0.015176429414463658, -0.022577963432775162};
                Real thermalEnergyFlux = -0.001531664896602873;
                std::vector<Real> const testFluxes = computeFluxes(midPoint,
                                                                   leftVxTurnOver,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Midpoint\n"
                                                "Right State: Right Vx turnover point\n"
                                                "HLLD State: Left Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.0013952100758668729, 0.0061359407125797273, 0, 0, 0.00065984543596031629, 0.0, 0.00069776606396793105, 0};
                std::vector<Real> const scalarFlux{ 0.001544494107257657, 0.0031093909889746947, 0.0046258388010795683};
                Real thermalEnergyFlux = 0.00040916715364737997;
                std::vector<Real> const testFluxes = computeFluxes(midPoint,
                                                                   rightVxTurnOver,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Vx turnover point\n"
                                                "Right State: Midpoint\n"
                                                "HLLD State: Left Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.0090024688079190333, 0.011769373146023688, 0, 0, 0.003725251767222792, 0.0, 0.0045418689996141555, 0};
                std::vector<Real> const scalarFlux{0.0099657107306674268, 0.020063068547205749, 0.029847813055181766};
                Real thermalEnergyFlux = 0.0020542406295284269;
                std::vector<Real> const testFluxes = computeFluxes(rightVxTurnOver,
                                                                   midPoint,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Vx turnover point\n"
                                                "Right State: Right rarefaction center\n"
                                                "HLLD State: Left Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.023310393229073981, 0.0033086897645311728, 0, 0, 0.0034208520409618887, 0.0, 0.011760413130542123, 0};
                std::vector<Real> const scalarFlux{0.025804547718589466, 0.051949973634547723, 0.077285939467198722};
                Real thermalEnergyFlux = 0.0053191138878843835;
                std::vector<Real> const testFluxes = computeFluxes(rightVxTurnOver,
                                                                   rightRarefactionCenter,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right rarefaction center\n"
                                                "Right State: Right Vx turnover point\n"
                                                "HLLD State: Left Star"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.33914253809565298, 0.46770133685446141, 0, 0, 0.46453338019960133, 0.0, 0.17077520175095764, 0};
                std::vector<Real> const scalarFlux{0.37542995185416178, 0.75581933514738364, 1.1244318966408966};
                Real thermalEnergyFlux = 0.1444638874418068;
                std::vector<Real> const testFluxes = computeFluxes(rightRarefactionCenter,
                                                                   rightVxTurnOver,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right rarefaction center\n"
                                                "Right State: Right Einfeldt Strong Rarefaction state\n"
                                                "HLLD State: Left"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{0.33976563660000003, 0.46733255780629601, 0, 0, 0.46427650313257612, 0.0, 0.17108896296000001, 0};
                std::vector<Real> const scalarFlux{0.37611972035917141, 0.75720798400261535, 1.1264977885722693};
                Real thermalEnergyFlux = 0.14472930749999999;
                std::vector<Real> const testFluxes = computeFluxes(rightRarefactionCenter,
                                                                   rightICs,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Right Einfeldt Strong Rarefaction state\n"
                                                "Right State: Right rarefaction center\n"
                                                "HLLD State: Left"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{2, 4.5750000000000002, 0, 0, 6.75, 0.0, 1, 0};
                std::vector<Real> const scalarFlux{2.2139950592000002, 4.4572370036000004, 6.6310283749999996};
                Real thermalEnergyFlux = 1.3499999999999996;
                std::vector<Real> const testFluxes = computeFluxes(rightICs,
                                                                   rightRarefactionCenter,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
    * \brief Test the HLLD Riemann Solver with the degenerate state
    *
    */
    TEST_F(tMHDCalculateHLLDFluxesCUDA,
           DegenerateStateCorrectInputExpectCorrectOutput)
    {
        // Constant Values
        Real const gamma = 5./3.;
        std::vector<Real> const primitiveScalar{1.1069975296, 2.2286185018, 3.3155141875};

        // State
        std::vector<Real> const      // | Density | X-Velocity | Y-Velocity | Z-Velocity | Pressure | X-Magnetic Field | Y-Magnetic Field | Z-Magnetic Field | Adiabatic Index | Passive Scalars |
            state = primitive2Conserved({1.0,       1.0,         1.0,         1.0,         1.0,       3.0E4,             1.0,               1.0},              gamma,            primitiveScalar);

        std::vector<Real> const fiducialFlux{1, -449999997, -29999, -29999, -59994, 0.0, -29999, -29999};
        std::vector<Real> const scalarFlux{1.1069975296000001, 2.2286185018000002, 3.3155141874999998};
        Real thermalEnergyFlux = 1.5;
        std::string const outputString {"Left State:  Degenerate state\n"
                                        "Right State: Degenerate state\n"
                                        "HLLD State: Left Double Star State"};

        // Compute the fluxes and check for correctness
        // Order of Fluxes is rho, vec(V), E, vec(B)
        // If you run into issues with the energy try 0.001953125 instead.
        // That's what I got when running the Athena solver on its own. Running
        // the Athena solver with theses tests gave me -0.00080700946455175148
        // though
        for (size_t direction = 0; direction < 3; direction++)
        {
            std::vector<Real> const testFluxes = computeFluxes(state,
                                                               state,
                                                               gamma,
                                                               direction);
            checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
    * \brief Test the HLLD Riemann Solver with all zeroes
    *
    */
    TEST_F(tMHDCalculateHLLDFluxesCUDA,
           AllZeroesExpectAllZeroes)
    {
        // Constant Values
        Real const gamma = 5./3.;

        // State
        size_t numElements = 8;
        #ifdef SCALAR
            numElements += 3;
        #endif // SCALAR

        std::vector<Real> const state(numElements, 0.0);
        std::vector<Real> const fiducialFlux(8,0.0);
        std::vector<Real> const scalarFlux(3,0.0);
        Real thermalEnergyFlux = 0.0;

        std::string const outputString {"Left State:  All zeroes\n"
                                        "Right State: All zeroes\n"
                                        "HLLD State: Right Star State"};

        for (size_t direction = 0; direction < 3; direction++)
        {
            // Compute the fluxes and check for correctness
            // Order of Fluxes is rho, vec(V), E, vec(B)
            std::vector<Real> const testFluxes = computeFluxes(state,
                                                               state,
                                                               gamma,
                                                               direction);
            checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
    * \brief Test the HLLD Riemann Solver with negative pressure, energy, and
      density.
    *
    */
    TEST_F(tMHDCalculateHLLDFluxesCUDA,
           UnphysicalValuesExpectAutomaticFix)
    {
        // Constant Values
        Real const gamma                  = 5./3.;

        // States
        std::vector<Real>                // | Density | X-Momentum | Y-Momentum | Z-Momentum | Energy   | X-Magnetic Field | Y-Magnetic Field | Z-Magnetic Field | Adiabatic Index | Passive Scalars |
            negativePressure              = { 1.0,      1.0,         1.0,         1.0,         1.5,       1.0,               1.0,               1.0},
            negativeEnergy                = { 1.0,      1.0,         1.0,         1.0,        -(5-gamma), 1.0,               1.0,               1.0},
            negativeDensity               = {-1.0,      1.0,         1.0,         1.0,         1.0,       1.0,               1.0,               1.0},
            negativeDensityEnergyPressure = {-1.0,     -1.0,        -1.0,        -1.0,         -gamma,    1.0,               1.0,               1.0},
            negativeDensityPressure       = {-1.0,      1.0,         1.0,         1.0,        -1.0,       1.0,               1.0,               1.0};

        #ifdef SCALAR
            std::vector<Real> const conservedScalar{1.1069975296, 2.2286185018, 3.3155141875};
            negativePressure.insert(negativePressure.begin()+5, conservedScalar.begin(), conservedScalar.begin() + NSCALARS);
            negativeEnergy.insert(negativeEnergy.begin()+5, conservedScalar.begin(), conservedScalar.begin() + NSCALARS);
            negativeDensity.insert(negativeDensity.begin()+5, conservedScalar.begin(), conservedScalar.begin() + NSCALARS);
            negativeDensityEnergyPressure.insert(negativeDensityEnergyPressure.begin()+5, conservedScalar.begin(), conservedScalar.begin() + NSCALARS);
            negativeDensityPressure.insert(negativeDensityPressure.begin()+5, conservedScalar.begin(), conservedScalar.begin() + NSCALARS);
        #endif  // SCALAR
        #ifdef  DE
            negativePressure.push_back(mhdUtils::computeThermalEnergy(negativePressure.at(4),negativePressure.at(0),negativePressure.at(1),negativePressure.at(2),negativePressure.at(3),negativePressure.at(5 + NSCALARS),negativePressure.at(6 + NSCALARS),negativePressure.at(7 + NSCALARS),gamma));
            negativeEnergy.push_back(mhdUtils::computeThermalEnergy(negativeEnergy.at(4),negativeEnergy.at(0),negativeEnergy.at(1),negativeEnergy.at(2),negativeEnergy.at(3),negativeEnergy.at(5 + NSCALARS),negativeEnergy.at(6 + NSCALARS),negativeEnergy.at(7 + NSCALARS),gamma));
            negativeDensity.push_back(mhdUtils::computeThermalEnergy(negativeDensity.at(4),negativeDensity.at(0),negativeDensity.at(1),negativeDensity.at(2),negativeDensity.at(3),negativeDensity.at(5 + NSCALARS),negativeDensity.at(6 + NSCALARS),negativeDensity.at(7 + NSCALARS),gamma));
            negativeDensityEnergyPressure.push_back(mhdUtils::computeThermalEnergy(negativeDensityEnergyPressure.at(4),negativeDensityEnergyPressure.at(0),negativeDensityEnergyPressure.at(1),negativeDensityEnergyPressure.at(2),negativeDensityEnergyPressure.at(3),negativeDensityEnergyPressure.at(5 + NSCALARS),negativeDensityEnergyPressure.at(6 + NSCALARS),negativeDensityEnergyPressure.at(7 + NSCALARS),gamma));
            negativeDensityPressure.push_back(mhdUtils::computeThermalEnergy(negativeDensityPressure.at(4),negativeDensityPressure.at(0),negativeDensityPressure.at(1),negativeDensityPressure.at(2),negativeDensityPressure.at(3),negativeDensityPressure.at(5 + NSCALARS),negativeDensityPressure.at(6 + NSCALARS),negativeDensityPressure.at(7 + NSCALARS),gamma));
        #endif  //DE

        for (size_t direction = 0; direction < 3; direction++)
        {
            {
                std::string const outputString {"Left State:  Negative Pressure\n"
                                                "Right State: Negative Pressure\n"
                                                "HLLD State: Left Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{1, 1.5, 0, 0, -1.6254793235168146e-16, 0, 0, 0};
                std::vector<Real> const scalarFlux{1.1069975296000001, 2.2286185018000002, 3.3155141874999998};
                Real thermalEnergyFlux = -1.5;
                std::vector<Real> const testFluxes = computeFluxes(negativePressure,
                                                                   negativePressure,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Negative Energy\n"
                                                "Right State: Negative Energy\n"
                                                "HLLD State: Left Star State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{1, 1.5, 0, 0, -1.5, 0, 0, 0};
                std::vector<Real> const scalarFlux{1.1069975296000001, 2.2286185018000002, 3.3155141874999998};
                Real thermalEnergyFlux = -6.333333333333333;
                std::vector<Real> const testFluxes = computeFluxes(negativeEnergy,
                                                                   negativeEnergy,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Negative Density\n"
                                                "Right State: Negative Density\n"
                                                "HLLD State: Left State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{1, 1E+20, 1e+20, 1e+20, -5e+19, 0, 0, 0};
                std::vector<Real> const scalarFlux{1.1069975296000002e+20, 2.2286185018000002e+20, 3.3155141874999997e+20};
                Real thermalEnergyFlux = -1.5000000000000001e+40;
                std::vector<Real> const testFluxes = computeFluxes(negativeDensity,
                                                                   negativeDensity,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Negative Density, Energy, and Pressure\n"
                                                "Right State: Negative Density, Energy, and Pressure\n"
                                                "HLLD State: Right State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{-1, 1E+20, 1E+20, 1E+20, 1.5E+20, 0, 0, 0};
                std::vector<Real> const scalarFlux{-1.1069975296000002e+20, -2.2286185018000002e+20, -3.3155141874999997e+20};
                Real thermalEnergyFlux = 1.5000000000000001e+40;
                std::vector<Real> const testFluxes = computeFluxes(negativeDensityEnergyPressure,
                                                                   negativeDensityEnergyPressure,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
            {
                std::string const outputString {"Left State:  Negative Density and Pressure\n"
                                                "Right State: Negative Density and Pressure\n"
                                                "HLLD State: Left State"};
                // Compute the fluxes and check for correctness
                // Order of Fluxes is rho, vec(V), E, vec(B)
                std::vector<Real> const fiducialFlux{1, 1e+20, 1e+20, 1e+20, -1.5e+20, 0, 0, 0};
                std::vector<Real> const scalarFlux{1.1069975296000002e+20, 2.2286185018000002e+20, 3.3155141874999997e+20};
                Real thermalEnergyFlux = -1.5000000000000001e+40;
                std::vector<Real> const testFluxes = computeFluxes(negativeDensityPressure,
                                                                   negativeDensityPressure,
                                                                   gamma,
                                                                   direction);
                checkResults(fiducialFlux, scalarFlux, thermalEnergyFlux, testFluxes, outputString, direction);
            }
        }
    }
    // =========================================================================

    // =========================================================================
    // End of integration tests for the entire HLLD solver. Unit tests are below
    // =========================================================================

    // =========================================================================
    // Unit tests for the contents of the _hlldInternal namespace
    // =========================================================================
    /*!
     * \brief A struct to hold some basic test values
     *
     */
    namespace
    {
        struct testParams
        {
            // List of cases
            std::vector<std::string> names{"Case 1", "Case 2"};

            // Conserved Variables
            double gamma = 5./3.;
            std::valarray<double> densityL  {21.50306776645775  , 48.316634031589935};
            std::valarray<double> densityR  {81.1217731762265   , 91.02955738853635};
            std::valarray<double> momentumXL{38.504606872151484 , 18.984145880030045};
            std::valarray<double> momentumXR{ 8.201811315045326 , 85.24863367778745};
            std::valarray<double> momentumYL{ 7.1046427940455015, 33.76182584816693};
            std::valarray<double> momentumYR{13.874767484202021 , 33.023492551299974};
            std::valarray<double> momentumZL{32.25700338919422  , 89.52561861038686};
            std::valarray<double> momentumZR{33.85305318830181  ,  8.664313303796256};
            std::valarray<double> energyL   {65.75120838109942  , 38.461354599479826};
            std::valarray<double> energyR   {18.88982523270516  , 83.65639784178894};
            std::valarray<double> magneticXL{92.75101068883114  , 31.588767769990532};
            std::valarray<double> magneticXR{93.66196246448985  , 84.3529879134052};
            std::valarray<double> magneticYL{12.297499156516622 , 63.74471969570406};
            std::valarray<double> magneticYR{84.9919141787549   , 35.910258841630984};
            std::valarray<double> magneticZL{46.224045698787776 , 37.70326455170754};
            std::valarray<double> magneticZR{34.852095153095384 , 24.052685003977757};
            // Star States
            std::valarray<double> densityStarL  {28.520995251761526 , 54.721668215064945};
            std::valarray<double> densityStarR  {49.09069570738605  , 72.68000504460609};
            std::valarray<double> momentumStarXL{48.96082367518151  , 97.15439466280228};
            std::valarray<double> momentumStarXR{65.74705433463932  , 94.5689655974538};
            std::valarray<double> momentumStarYL{44.910034185328996 , 78.60179936059853};
            std::valarray<double> momentumStarYR{51.642522487399276 , 44.63864007208728};
            std::valarray<double> momentumStarZL{39.78163555990428  , 63.01612978428839};
            std::valarray<double> momentumStarZR{33.47900698769427  , 52.19410653341197};
            std::valarray<double> energyStarL   { 6.579867455284738 , 30.45043664908369};
            std::valarray<double> energyStarR   {90.44484278669114  , 61.33664731346812};
            std::valarray<double> magneticStarXL{49.81491527582234  , 62.379765828560906};
            std::valarray<double> magneticStarXR{67.77402751903804  , 64.62226739788758};
            std::valarray<double> magneticStarYL{62.09348829143065  , 54.27916744403672};
            std::valarray<double> magneticStarYR{26.835645069149873 , 98.97444628327318};
            std::valarray<double> magneticStarZL{62.765890944643196 , 93.26765455509641};
            std::valarray<double> magneticStarZR{ 7.430231695917344 , 10.696380763901459};
            // Double Star State
            std::valarray<double> momentumDoubleStarXL{75.42525315887075  , 83.87480678359029};
            std::valarray<double> momentumDoubleStarYL{22.56132540660678  , 76.11074421934487};
            std::valarray<double> momentumDoubleStarZL{27.83908778933224  , 28.577101567661465};
            std::valarray<double> energyDoubleStar    {45.83202455707669  , 55.4553014145573};
            std::valarray<double> magneticDoubleStarY {20.943239839455895 , 83.8514810487021};
            std::valarray<double> magneticDoubleStarZ {83.3802438268807   , 80.36671251730783};
            // Fluxes
            std::valarray<double> densityFluxL     {12.939239309626116 , 81.71524586517073};
            std::valarray<double> momentumFluxXL   {65.05481464917627  , 56.09885069707803};
            std::valarray<double> momentumFluxYL   {73.67692845586782  ,  2.717246983403787};
            std::valarray<double> momentumFluxZL   {16.873647595664387 , 39.70132983192873};
            std::valarray<double> energyFluxL      {52.71888731972469  , 81.63926176158796};
            std::valarray<double> magneticFluxXL   {67.7412464028116   , 42.85301340921149};
            std::valarray<double> magneticFluxYL   {58.98928445415967  , 57.04344459221359};
            std::valarray<double> magneticFluxZL   {29.976925743532302 , 97.73329827141359};
            std::valarray<double> momentumStarFluxX{74.90125547448865  , 26.812722601652684};
            std::valarray<double> momentumStarFluxY{16.989138610622945 , 48.349566649914976};
            std::valarray<double> momentumStarFluxZ{38.541822734846185 , 61.22843961052538};
            std::valarray<double> energyStarFlux   {19.095105176247017 , 45.43224973313112};
            std::valarray<double> magneticStarFluxY{96.23964526624277  , 33.05337536594796};
            std::valarray<double> magneticStarFluxZ{86.22516928268347  , 15.62102082410738};

            // Derived/Primitive variables
            std::valarray<double> velocityXL     = momentumXL / densityL;
            std::valarray<double> velocityXR     = momentumXR / densityR;
            std::valarray<double> velocityYL     = momentumYL / densityL;
            std::valarray<double> velocityYR     = momentumYR / densityR;
            std::valarray<double> velocityZL     = momentumZL / densityL;
            std::valarray<double> velocityZR     = momentumZR / densityR;
            std::valarray<double> totalPressureStarL{66.80958736783934  , 72.29644038317676};
            std::vector<double> gasPressureL;
            std::vector<double> gasPressureR;
            std::vector<double> totalPressureL;
            std::vector<double> totalPressureR;
            // Star State
            std::valarray<double> velocityStarXL = momentumStarXL / densityStarL;
            std::valarray<double> velocityStarXR = momentumStarXR / densityStarR;
            std::valarray<double> velocityStarYL = momentumStarYL / densityStarL;
            std::valarray<double> velocityStarYR = momentumStarYR / densityStarR;
            std::valarray<double> velocityStarZL = momentumStarZL / densityStarL;
            std::valarray<double> velocityStarZR = momentumStarZR / densityStarR;
            // Double Star State
            std::valarray<double> velocityDoubleStarXL = momentumDoubleStarXL / densityStarL;
            std::valarray<double> velocityDoubleStarYL = momentumDoubleStarYL / densityStarL;
            std::valarray<double> velocityDoubleStarZL = momentumDoubleStarZL / densityStarL;
            // Other
            std::valarray<double> speedM            {68.68021569453585  , 70.08236749169825};
            std::valarray<double> speedSide         {70.37512772923496  ,  3.6579130085113265};
            testParams()
            {
                for (size_t i = 0; i < names.size(); i++)
                {
                    gasPressureL.push_back(mhdUtils::computeGasPressure(energyL[i], densityL[i], momentumXL[i], momentumYL[i], momentumZL[i], magneticXL[i], magneticYL[i], magneticZL[i], gamma));
                    gasPressureR.push_back(mhdUtils::computeGasPressure(energyR[i], densityR[i], momentumXR[i], momentumYR[i], momentumZR[i], magneticXR[i], magneticYR[i], magneticZR[i], gamma));
                    totalPressureL.push_back(mhdUtils::computeTotalPressure(gasPressureL.back(), magneticXL[i], magneticYL[i], magneticZL[i]));
                    totalPressureR.push_back(mhdUtils::computeTotalPressure(gasPressureL.back(), magneticXR[i], magneticYR[i], magneticZR[i]));
                }
            }
        };
    }
    // =========================================================================

    // =========================================================================
    /*!
     * \brief Test the _hlldInternal::_approximateWaveSpeeds function
     *
     */
    TEST(tMHDHlldInternalApproximateWaveSpeeds,
         CorrectInputExpectCorrectOutput)
    {
        testParams const parameters;
        std::vector<double> const fiducialSpeedL      {-22.40376497145191, -11.190385012513822};
        std::vector<double> const fiducialSpeedR      {24.295526347371595, 12.519790189404299};
        std::vector<double> const fiducialSpeedM      {-0.81760587897407833, -0.026643804611559244};
        std::vector<double> const fiducialSpeedStarL  {-19.710500632936679, -4.4880642018724357};
        std::vector<double> const fiducialSpeedStarR  {9.777062240423124, 9.17474383484066};
        std::vector<double> const fiducialDensityStarL{24.101290139122913, 50.132466596958501};
        std::vector<double> const fiducialDensityStarR{78.154104734671265, 84.041595114910123};

        double testSpeedL = 0;
        double testSpeedR = 0;
        double testSpeedM = 0;
        double testSpeedStarL = 0;
        double testSpeedStarR = 0;
        double testDensityStarL = 0;
        double testDensityStarR = 0;

        for (size_t i = 0; i < parameters.names.size(); i++)
        {
            _hlldInternal::_approximateWaveSpeeds(parameters.densityL[i],
                                                  parameters.momentumXL[i],
                                                  parameters.momentumYL[i],
                                                  parameters.momentumZL[i],
                                                  parameters.velocityXL[i],
                                                  parameters.velocityYL[i],
                                                  parameters.velocityZL[i],
                                                  parameters.gasPressureL[i],
                                                  parameters.totalPressureL[i],
                                                  parameters.magneticXL[i],
                                                  parameters.magneticYL[i],
                                                  parameters.magneticZL[i],
                                                  parameters.densityR[i],
                                                  parameters.momentumXR[i],
                                                  parameters.momentumYR[i],
                                                  parameters.momentumZR[i],
                                                  parameters.velocityXR[i],
                                                  parameters.velocityYR[i],
                                                  parameters.velocityZR[i],
                                                  parameters.gasPressureR[i],
                                                  parameters.totalPressureR[i],
                                                  parameters.magneticXR[i],
                                                  parameters.magneticYR[i],
                                                  parameters.magneticZR[i],
                                                  parameters.gamma,
                                                  testSpeedL,
                                                  testSpeedR,
                                                  testSpeedM,
                                                  testSpeedStarL,
                                                  testSpeedStarR,
                                                  testDensityStarL,
                                                  testDensityStarR);
            // Now check results
            testingUtilities::checkResults(fiducialSpeedL[i],
                                           testSpeedL,
                                           parameters.names.at(i) + ", SpeedL");
            testingUtilities::checkResults(fiducialSpeedR.at(i),
                                           testSpeedR,
                                           parameters.names.at(i) + ", SpeedR");
            testingUtilities::checkResults(fiducialSpeedM.at(i),
                                           testSpeedM,
                                           parameters.names.at(i) + ", SpeedM");
            testingUtilities::checkResults(fiducialSpeedStarL.at(i),
                                           testSpeedStarL,
                                           parameters.names.at(i) + ", SpeedStarL");
            testingUtilities::checkResults(fiducialSpeedStarR.at(i),
                                           testSpeedStarR,
                                           parameters.names.at(i) + ", SpeedStarR");
            testingUtilities::checkResults(fiducialDensityStarL.at(i),
                                           testDensityStarL,
                                           parameters.names.at(i) + ", DensityStarL");
            testingUtilities::checkResults(fiducialDensityStarR.at(i),
                                           testDensityStarR,
                                           parameters.names.at(i) + ", DensityStarR");
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
     * \brief Test the _hlldInternal::_starFluxes function in the non-degenerate
     * case
     *
     */
     TEST(tMHDHlldInternalStarFluxes,
          CorrectInputNonDegenerateExpectCorrectOutput)
    {
        testParams const parameters;

        std::vector<double> const fiducialVelocityStarY    {12.831290892281075, 12.92610185957192};
        std::vector<double> const fiducialVelocityStarZ    {48.488664548015286, 9.0850326944201107};
        std::vector<double> const fiducialEnergyStar       {1654897.6912410262, 956.83439334487116};
        std::vector<double> const fiducialMagneticStarY    {-186.47142421374559, 2.6815421494204679};
        std::vector<double> const fiducialMagneticStarZ    {-700.91191100481922, 1.5860591049546646};
        std::vector<double> const fiducialDensityStarFlux  {506.82678248238807, 105.14430372486369};
        std::vector<double> const fiducialMomentumStarFluxX{135208.06632708258, 14014.840899433098};
        std::vector<double> const fiducialMomentumStarFluxY{25328.25203616685, 2466.5997745560339};
        std::vector<double> const fiducialMomentumStarFluxZ{95071.711914347878, 1530.7490710422007};
        std::vector<double> const fiducialEnergyStarFlux   {116459061.8691024, 3440.9679468544314};
        std::vector<double> const fiducialMagneticStarFluxY{-13929.399086330559, -166.32034689537392};
        std::vector<double> const fiducialMagneticStarFluxZ{-52549.811458376971, -34.380297363339892};

        double testVelocityStarY;
        double testVelocityStarZ;
        double testEnergyStar;
        double testMagneticStarY;
        double testMagneticStarZ;
        double testDensityStarFlux;
        double testMomentumStarFluxX;
        double testMomentumStarFluxY;
        double testMomentumStarFluxZ;
        double testEnergyStarFlux;
        double testMagneticStarFluxY;
        double testMagneticStarFluxZ;

        for (size_t i = 0; i < parameters.names.size(); i++)
        {
            _hlldInternal::_starFluxes(parameters.speedM[i],
                                       parameters.speedSide[i],
                                       parameters.densityL[i],
                                       parameters.velocityXL[i],
                                       parameters.velocityYL[i],
                                       parameters.velocityZL[i],
                                       parameters.momentumXL[i],
                                       parameters.momentumYL[i],
                                       parameters.momentumZL[i],
                                       parameters.energyL[i],
                                       parameters.totalPressureL[i],
                                       parameters.magneticXL[i],
                                       parameters.magneticYL[i],
                                       parameters.magneticZL[i],
                                       parameters.densityStarL[i],
                                       parameters.totalPressureStarL[i],
                                       parameters.densityFluxL[i],
                                       parameters.momentumFluxXL[i],
                                       parameters.momentumFluxYL[i],
                                       parameters.momentumFluxZL[i],
                                       parameters.energyFluxL[i],
                                       parameters.magneticFluxYL[i],
                                       parameters.magneticFluxZL[i],
                                       testVelocityStarY,
                                       testVelocityStarZ,
                                       testEnergyStar,
                                       testMagneticStarY,
                                       testMagneticStarZ,
                                       testDensityStarFlux,
                                       testMomentumStarFluxX,
                                       testMomentumStarFluxY,
                                       testMomentumStarFluxZ,
                                       testEnergyStarFlux,
                                       testMagneticStarFluxY,
                                       testMagneticStarFluxZ);

            // Now check results
            testingUtilities::checkResults(fiducialVelocityStarY[i],
                                            testVelocityStarY,
                                            parameters.names.at(i) + ", VelocityStarY");
            testingUtilities::checkResults(fiducialVelocityStarZ[i],
                                            testVelocityStarZ,
                                            parameters.names.at(i) + ", VelocityStarZ");
            testingUtilities::checkResults(fiducialEnergyStar[i],
                                           testEnergyStar,
                                           parameters.names.at(i) + ", EnergyStar");
            testingUtilities::checkResults(fiducialMagneticStarY[i],
                                            testMagneticStarY,
                                            parameters.names.at(i) + ", MagneticStarY");
            testingUtilities::checkResults(fiducialMagneticStarZ[i],
                                            testMagneticStarZ,
                                            parameters.names.at(i) + ", MagneticStarZ");
            testingUtilities::checkResults(fiducialDensityStarFlux[i],
                                            testDensityStarFlux,
                                            parameters.names.at(i) + ", DensityStarFlux");
            testingUtilities::checkResults(fiducialMomentumStarFluxX[i],
                                            testMomentumStarFluxX,
                                            parameters.names.at(i) + ", MomentumStarFluxX");
            testingUtilities::checkResults(fiducialMomentumStarFluxY[i],
                                            testMomentumStarFluxY,
                                            parameters.names.at(i) + ", MomentumStarFluxY");
            testingUtilities::checkResults(fiducialMomentumStarFluxZ[i],
                                            testMomentumStarFluxZ,
                                            parameters.names.at(i) + ", MomentumStarFluxZ");
            testingUtilities::checkResults(fiducialEnergyStarFlux[i],
                                            testEnergyStarFlux,
                                            parameters.names.at(i) + ", EnergyStarFlux");
            testingUtilities::checkResults(fiducialMagneticStarFluxY[i],
                                            testMagneticStarFluxY,
                                            parameters.names.at(i) + ", MagneticStarFluxY");
            testingUtilities::checkResults(fiducialMagneticStarFluxZ[i],
                                            testMagneticStarFluxZ,
                                            parameters.names.at(i) + ", MagneticStarFluxZ");
        }
    }

    /*!
     * \brief Test the _hlldInternal::_starFluxes function in the degenerate
     * case
     *
     */
     TEST(tMHDHlldInternalStarFluxes,
          CorrectInputDegenerateExpectCorrectOutput)
    {
        testParams const parameters;

        // Used to get us into the degenerate case
        double const totalPressureStarMultiplier = 1E15;

        std::vector<double> const fiducialVelocityStarY    {0.33040135813215948, 0.69876195899931859};
        std::vector<double> const fiducialVelocityStarZ    {1.500111692877206, 1.8528943583250035};
        std::vector<double> const fiducialEnergyStar       {2.7072182962581443e+18, -76277716432851392};
        std::vector<double> const fiducialMagneticStarY    {12.297499156516622, 63.744719695704063};
        std::vector<double> const fiducialMagneticStarZ    {46.224045698787776, 37.703264551707541};
        std::vector<double> const fiducialDensityStarFlux  {506.82678248238807, 105.14430372486369};
        std::vector<double> const fiducialMomentumStarFluxX{135208.06632708258, 14014.840899433098};
        std::vector<double> const fiducialMomentumStarFluxY{236.85804348470396, 19.08858135095122};
        std::vector<double> const fiducialMomentumStarFluxZ{757.76012607552047, 83.112898961023902};
        std::vector<double> const fiducialEnergyStarFlux   {1.9052083339008875e+20, -2.7901725119926531e+17};
        std::vector<double> const fiducialMagneticStarFluxY{58.989284454159673, 57.043444592213589};
        std::vector<double> const fiducialMagneticStarFluxZ{29.976925743532302, 97.733298271413588};

        double testVelocityStarY;
        double testVelocityStarZ;
        double testEnergyStar;
        double testMagneticStarY;
        double testMagneticStarZ;
        double testDensityStarFlux;
        double testMomentumStarFluxX;
        double testMomentumStarFluxY;
        double testMomentumStarFluxZ;
        double testEnergyStarFlux;
        double testMagneticStarFluxY;
        double testMagneticStarFluxZ;

        for (size_t i = 0; i < parameters.names.size(); i++)
        {
            _hlldInternal::_starFluxes(parameters.speedM[i],
                                        parameters.speedSide[i],
                                        parameters.densityL[i],
                                        parameters.velocityXL[i],
                                        parameters.velocityYL[i],
                                        parameters.velocityZL[i],
                                        parameters.momentumXL[i],
                                        parameters.momentumYL[i],
                                        parameters.momentumZL[i],
                                        parameters.energyL[i],
                                        parameters.totalPressureL[i],
                                        parameters.magneticXL[i],
                                        parameters.magneticYL[i],
                                        parameters.magneticZL[i],
                                        parameters.densityStarL[i],
                                        parameters.totalPressureStarL[i] * totalPressureStarMultiplier,
                                        parameters.densityFluxL[i],
                                        parameters.momentumFluxXL[i],
                                        parameters.momentumFluxYL[i],
                                        parameters.momentumFluxZL[i],
                                        parameters.energyFluxL[i],
                                        parameters.magneticFluxYL[i],
                                        parameters.magneticFluxZL[i],
                                        testVelocityStarY,
                                        testVelocityStarZ,
                                        testEnergyStar,
                                        testMagneticStarY,
                                        testMagneticStarZ,
                                        testDensityStarFlux,
                                        testMomentumStarFluxX,
                                        testMomentumStarFluxY,
                                        testMomentumStarFluxZ,
                                        testEnergyStarFlux,
                                        testMagneticStarFluxY,
                                        testMagneticStarFluxZ);

            // Now check results
            testingUtilities::checkResults(fiducialVelocityStarY[i],
                                            testVelocityStarY,
                                            parameters.names.at(i) + ", VelocityStarY");
            testingUtilities::checkResults(fiducialVelocityStarZ[i],
                                            testVelocityStarZ,
                                            parameters.names.at(i) + ", VelocityStarZ");
            testingUtilities::checkResults(fiducialEnergyStar[i],
                                            testEnergyStar,
                                            parameters.names.at(i) + ", EnergyStar");
            testingUtilities::checkResults(fiducialMagneticStarY[i],
                                            testMagneticStarY,
                                            parameters.names.at(i) + ", MagneticStarY");
            testingUtilities::checkResults(fiducialMagneticStarZ[i],
                                            testMagneticStarZ,
                                            parameters.names.at(i) + ", MagneticStarZ");
            testingUtilities::checkResults(fiducialDensityStarFlux[i],
                                            testDensityStarFlux,
                                            parameters.names.at(i) + ", DensityStarFlux");
            testingUtilities::checkResults(fiducialMomentumStarFluxX[i],
                                            testMomentumStarFluxX,
                                            parameters.names.at(i) + ", MomentumStarFluxX");
            testingUtilities::checkResults(fiducialMomentumStarFluxY[i],
                                            testMomentumStarFluxY,
                                            parameters.names.at(i) + ", MomentumStarFluxY");
            testingUtilities::checkResults(fiducialMomentumStarFluxZ[i],
                                            testMomentumStarFluxZ,
                                            parameters.names.at(i) + ", MomentumStarFluxZ");
            testingUtilities::checkResults(fiducialEnergyStarFlux[i],
                                            testEnergyStarFlux,
                                            parameters.names.at(i) + ", EnergyStarFlux");
            testingUtilities::checkResults(fiducialMagneticStarFluxY[i],
                                            testMagneticStarFluxY,
                                            parameters.names.at(i) + ", MagneticStarFluxY");
            testingUtilities::checkResults(fiducialMagneticStarFluxZ[i],
                                            testMagneticStarFluxZ,
                                            parameters.names.at(i) + ", MagneticStarFluxZ");
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
     * \brief Test the _hlldInternal::_nonStarFluxes function
     *
     */
    TEST(tMHDHlldInternalNonStarFluxes,
         CorrectInputExpectCorrectOutput)
    {
        testParams const parameters;

        std::vector<double> const fiducialDensityFlux  {38.504606872151484, 18.984145880030045};
        std::vector<double> const fiducialMomentumFluxX{-3088.4810263278778, 2250.9966820900618};
        std::vector<double> const fiducialMomentumFluxY{-1127.8835013070616, -2000.3517480656785};
        std::vector<double> const fiducialMomentumFluxZ{-4229.5657456907293, -1155.8240512956793};
        std::vector<double> const fiducialMagneticFluxY{-8.6244637840856555, 2.9729840344910059};
        std::vector<double> const fiducialMagneticFluxZ{-56.365490339906408, -43.716615275067923};
        std::vector<double> const fiducialEnergyFlux   {-12344.460641662206, -2717.2127176227905};

        double testDensityFlux;
        double testMomentumFluxX;
        double testMomentumFluxY;
        double testMomentumFluxZ;
        double testMagneticFluxY;
        double testMagneticFluxZ;
        double testEnergyFlux;

        for (size_t i = 0; i < parameters.names.size(); i++)
        {
            _hlldInternal::_nonStarFluxes(parameters.momentumXL[i],
                                          parameters.velocityXL[i],
                                          parameters.velocityYL[i],
                                          parameters.velocityZL[i],
                                          parameters.totalPressureL[i],
                                          parameters.energyL[i],
                                          parameters.magneticXL[i],
                                          parameters.magneticYL[i],
                                          parameters.magneticZL[i],
                                          testDensityFlux,
                                          testMomentumFluxX,
                                          testMomentumFluxY,
                                          testMomentumFluxZ,
                                          testMagneticFluxY,
                                          testMagneticFluxZ,
                                          testEnergyFlux);

            // Now check results
            testingUtilities::checkResults(fiducialDensityFlux[i],
                                            testDensityFlux,
                                            parameters.names.at(i) + ", DensityFlux");
            testingUtilities::checkResults(fiducialMomentumFluxX[i],
                                           testMomentumFluxX,
                                           parameters.names.at(i) + ", MomentumFluxX");
            testingUtilities::checkResults(fiducialMomentumFluxY[i],
                                           testMomentumFluxY,
                                           parameters.names.at(i) + ", MomentumFluxY");
            testingUtilities::checkResults(fiducialMomentumFluxZ[i],
                                           testMomentumFluxZ,
                                           parameters.names.at(i) + ", MomentumFluxZ");
            testingUtilities::checkResults(fiducialMagneticFluxY[i],
                                           testMagneticFluxY,
                                           parameters.names.at(i) + ", MagneticFluxY");
            testingUtilities::checkResults(fiducialMagneticFluxZ[i],
                                           testMagneticFluxZ,
                                           parameters.names.at(i) + ", MagneticFluxZ");
            testingUtilities::checkResults(fiducialEnergyFlux[i],
                                           testEnergyFlux,
                                           parameters.names.at(i) + ", EnergyFlux");
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
     * \brief Test the _hlldInternal::_dotProduct function
     *
     */
    TEST(tMHDHlldInternalDotProduct,
         CorrectInputExpectCorrectOutput)
    {
        testParams const parameters;

        std::vector<double> const fiducialDotProduct{5149.7597411033557,6127.2319832451567};

        double testDotProduct;

        for (size_t i = 0; i < parameters.names.size(); i++)
        {
            testDotProduct = _hlldInternal::_dotProduct(parameters.momentumXL[i],
                                                        parameters.momentumYL[i],
                                                        parameters.momentumZL[i],
                                                        parameters.magneticXL[i],
                                                        parameters.magneticYL[i],
                                                        parameters.magneticZL[i]);

            // Now check results
            testingUtilities::checkResults(fiducialDotProduct[i],
                                           testDotProduct,
                                           parameters.names.at(i) + ", DotProduct");
            }
    }
    // =========================================================================

    // =========================================================================
    /*!
     * \brief Test the _hlldInternal::_doubleStarState function. Non-degenerate
     * state
     *
    */
    TEST(tMHDHlldInternalDoubleStarState,
         CorrectInputNonDegenerateExpectCorrectOutput)
    {
        testParams const parameters;

        double const fixedEpsilon = 7E-12;

        std::vector<double> const fiducialVelocityDoubleStarY{-1.5775383335759607, 3.803188977150934};
        std::vector<double> const fiducialVelocityDoubleStarZ{-3.4914062207842482, -4.2662645349592765};
        std::vector<double> const fiducialMagneticDoubleStarY{45.259313435283325, 71.787329583230417};
        std::vector<double> const fiducialMagneticDoubleStarZ{36.670978215630669, 53.189673238238178};
        std::vector<double> const fiducialEnergyDoubleStarL  {-2048.1953674500514, -999.79694164635089};
        std::vector<double> const fiducialEnergyDoubleStarR  {1721.0582276783764, 252.04716752257781};

        double testVelocityDoubleStarY;
        double testVelocityDoubleStarZ;
        double testMagneticDoubleStarY;
        double testMagneticDoubleStarZ;
        double testEnergyDoubleStarL;
        double testEnergyDoubleStarR;

        for (size_t i = 0; i < parameters.names.size(); i++)
        {
            _hlldInternal::_doubleStarState(parameters.speedM[i],
                                            parameters.magneticXL[i],
                                            parameters.totalPressureStarL[i],
                                            parameters.densityStarL[i],
                                            parameters.velocityStarYL[i],
                                            parameters.velocityStarZL[i],
                                            parameters.energyStarL[i],
                                            parameters.magneticStarYL[i],
                                            parameters.magneticStarZL[i],
                                            parameters.densityStarR[i],
                                            parameters.velocityStarYR[i],
                                            parameters.velocityStarZR[i],
                                            parameters.energyStarR[i],
                                            parameters.magneticStarYR[i],
                                            parameters.magneticStarZR[i],
                                            testVelocityDoubleStarY,
                                            testVelocityDoubleStarZ,
                                            testMagneticDoubleStarY,
                                            testMagneticDoubleStarZ,
                                            testEnergyDoubleStarL,
                                            testEnergyDoubleStarR);

            // Now check results
            testingUtilities::checkResults(fiducialVelocityDoubleStarY[i],
                                           testVelocityDoubleStarY,
                                           parameters.names.at(i) + ", VelocityDoubleStarY");
            testingUtilities::checkResults(fiducialVelocityDoubleStarZ[i],
                                           testVelocityDoubleStarZ,
                                           parameters.names.at(i) + ", VelocityDoubleStarZ");
            testingUtilities::checkResults(fiducialMagneticDoubleStarY[i],
                                           testMagneticDoubleStarY,
                                           parameters.names.at(i) + ", MagneticDoubleStarY");
            testingUtilities::checkResults(fiducialMagneticDoubleStarZ[i],
                                           testMagneticDoubleStarZ,
                                           parameters.names.at(i) + ", MagneticDoubleStarZ");
            testingUtilities::checkResults(fiducialEnergyDoubleStarL[i],
                                           testEnergyDoubleStarL,
                                           parameters.names.at(i) + ", EnergyDoubleStarL");
            testingUtilities::checkResults(fiducialEnergyDoubleStarR[i],
                                           testEnergyDoubleStarR,
                                           parameters.names.at(i) + ", EnergyDoubleStarR",
                                           fixedEpsilon);
        }
    }

    /*!
     * \brief Test the _hlldInternal::_doubleStarState function in the
     * degenerate state.
     *
    */
    TEST(tMHDHlldInternalDoubleStarState,
         CorrectInputDegenerateExpectCorrectOutput)
    {
        testParams const parameters;

        std::vector<double> const fiducialVelocityDoubleStarY{1.5746306813243216, 1.4363926014039052};
        std::vector<double> const fiducialVelocityDoubleStarZ{1.3948193325212686, 1.1515754515491903};
        std::vector<double> const fiducialMagneticDoubleStarY{62.093488291430653, 54.279167444036723};
        std::vector<double> const fiducialMagneticDoubleStarZ{62.765890944643196, 93.267654555096414};
        std::vector<double> const fiducialEnergyDoubleStarL  {6.579867455284738, 30.450436649083692};
        std::vector<double> const fiducialEnergyDoubleStarR  {90.44484278669114, 61.33664731346812};

        double testVelocityDoubleStarY;
        double testVelocityDoubleStarZ;
        double testMagneticDoubleStarY;
        double testMagneticDoubleStarZ;
        double testEnergyDoubleStarL;
        double testEnergyDoubleStarR;

        for (size_t i = 0; i < parameters.names.size(); i++)
        {
            _hlldInternal::_doubleStarState(parameters.speedM[i],
                                            0.0,
                                            parameters.totalPressureStarL[i],
                                            parameters.densityStarL[i],
                                            parameters.velocityStarYL[i],
                                            parameters.velocityStarZL[i],
                                            parameters.energyStarL[i],
                                            parameters.magneticStarYL[i],
                                            parameters.magneticStarZL[i],
                                            parameters.densityStarR[i],
                                            parameters.velocityStarYR[i],
                                            parameters.velocityStarZR[i],
                                            parameters.energyStarR[i],
                                            parameters.magneticStarYR[i],
                                            parameters.magneticStarZR[i],
                                            testVelocityDoubleStarY,
                                            testVelocityDoubleStarZ,
                                            testMagneticDoubleStarY,
                                            testMagneticDoubleStarZ,
                                            testEnergyDoubleStarL,
                                            testEnergyDoubleStarR);
            // Now check results
            testingUtilities::checkResults(fiducialVelocityDoubleStarY[i],
                                            testVelocityDoubleStarY,
                                            parameters.names.at(i) + ", VelocityDoubleStarY");
            testingUtilities::checkResults(fiducialVelocityDoubleStarZ[i],
                                            testVelocityDoubleStarZ,
                                            parameters.names.at(i) + ", VelocityDoubleStarZ");
            testingUtilities::checkResults(fiducialMagneticDoubleStarY[i],
                                            testMagneticDoubleStarY,
                                            parameters.names.at(i) + ", MagneticDoubleStarY");
            testingUtilities::checkResults(fiducialMagneticDoubleStarZ[i],
                                            testMagneticDoubleStarZ,
                                            parameters.names.at(i) + ", MagneticDoubleStarZ");
            testingUtilities::checkResults(fiducialEnergyDoubleStarL[i],
                                            testEnergyDoubleStarL,
                                            parameters.names.at(i) + ", EnergyDoubleStarL");
            testingUtilities::checkResults(fiducialEnergyDoubleStarR[i],
                                            testEnergyDoubleStarR,
                                            parameters.names.at(i) + ", EnergyDoubleStarR");
        }
    }
    // =========================================================================

    // =========================================================================
    /*!
     * \brief Test the _hlldInternal::_doubleStarFluxes function
     *
     */
    TEST(tMHDHlldInternalDoubleStarFluxes,
         CorrectInputExpectCorrectOutput)
    {
        testParams const parameters;

        std::vector<double> const fiducialMomentumDoubleStarFluxX{1937.3388606704509, -21.762854649386174};
        std::vector<double> const fiducialMomentumDoubleStarFluxY{-1555.8040962754276, 39.237503643804175};
        std::vector<double> const fiducialMomentumDoubleStarFluxZ{-801.91650203165148, -64.746529703562871};
        std::vector<double> const fiducialEnergyDoubleStarFlux   {2781.4706748628528, 136.89786983482355};
        std::vector<double> const fiducialMagneticDoubleStarFluxY{-2799.7143456312342, 141.2263259922299};
        std::vector<double> const fiducialMagneticDoubleStarFluxZ{1536.9628864256708, -31.569502877970095};


        double testMomentumDoubleStarFluxX;
        double testMomentumDoubleStarFluxY;
        double testMomentumDoubleStarFluxZ;
        double testEnergyDoubleStarFlux;
        double testMagneticDoubleStarFluxY;
        double testMagneticDoubleStarFluxZ;

        for (size_t i = 0; i < parameters.names.size(); i++)
        {
            _hlldInternal::_doubleStarFluxes(parameters.speedSide[i],
                                             parameters.momentumStarFluxX[i],
                                             parameters.momentumStarFluxY[i],
                                             parameters.momentumStarFluxZ[i],
                                             parameters.energyStarFlux[i],
                                             parameters.magneticStarFluxY[i],
                                             parameters.magneticStarFluxZ[i],
                                             parameters.densityStarL[i],
                                             parameters.velocityStarXL[i],
                                             parameters.velocityStarYL[i],
                                             parameters.velocityStarZL[i],
                                             parameters.energyStarL[i],
                                             parameters.magneticStarYL[i],
                                             parameters.magneticStarZL[i],
                                             parameters.velocityDoubleStarXL[i],
                                             parameters.velocityDoubleStarYL[i],
                                             parameters.velocityDoubleStarZL[i],
                                             parameters.energyDoubleStar[i],
                                             parameters.magneticDoubleStarY[i],
                                             parameters.magneticDoubleStarZ[i],
                                             testMomentumDoubleStarFluxX,
                                             testMomentumDoubleStarFluxY,
                                             testMomentumDoubleStarFluxZ,
                                             testEnergyDoubleStarFlux,
                                             testMagneticDoubleStarFluxY,
                                             testMagneticDoubleStarFluxZ);

            // Now check results
            testingUtilities::checkResults(fiducialMomentumDoubleStarFluxX[i],
                                           testMomentumDoubleStarFluxX,
                                           parameters.names.at(i) + ", MomentumDoubleStarFluxX");
            testingUtilities::checkResults(fiducialMomentumDoubleStarFluxY[i],
                                           testMomentumDoubleStarFluxY,
                                           parameters.names.at(i) + ", MomentumDoubleStarFluxY");
            testingUtilities::checkResults(fiducialMomentumDoubleStarFluxZ[i],
                                           testMomentumDoubleStarFluxZ,
                                           parameters.names.at(i) + ", MomentumDoubleStarFluxZ");
            testingUtilities::checkResults(fiducialEnergyDoubleStarFlux[i],
                                           testEnergyDoubleStarFlux,
                                           parameters.names.at(i) + ", EnergyDoubleStarFlux");
            testingUtilities::checkResults(fiducialMagneticDoubleStarFluxY[i],
                                           testMagneticDoubleStarFluxY,
                                           parameters.names.at(i) + ", MagneticDoubleStarFluxY");
            testingUtilities::checkResults(fiducialMagneticDoubleStarFluxZ[i],
                                           testMagneticDoubleStarFluxZ,
                                           parameters.names.at(i) + ", MagneticDoubleStarFluxZ");
            }
    }
    // =========================================================================

    // =========================================================================
    /*!
     * \brief Test the _hlldInternal::_returnFluxes function
     *
     */
    TEST(tMHDHlldInternalReturnFluxes,
         CorrectInputExpectCorrectOutput)
    {
        double const dummyValue    = 999;
        double const densityFlux   = 1;
        double const momentumFluxX = 2;
        double const momentumFluxY = 3;
        double const momentumFluxZ = 4;
        double const energyFlux    = 5;
        double const magneticFluxY = 6;
        double const magneticFluxZ = 7;

        int threadId = 0;
        int n_cells = 10;
        int nFields = 8;  // Total number of conserved fields
        #ifdef  SCALAR
            nFields += NSCALARS;
        #endif  // SCALAR
        #ifdef  DE
            nFields++;
        #endif  //DE

        // Lambda for finding indices and check if they're correct
        auto findIndex = [](std::vector<double> const &vec,
                            double const &num,
                            int const &fidIndex,
                            std::string const &name)
        {
            int index = std::distance(vec.begin(), std::find(vec.begin(), vec.end(), num));
            // EXPECT_EQ(fidIndex, index) << "Error in " << name << " index" << std::endl;

            return index;
        };

        for (size_t direction = 0; direction < 3; direction++)
        {
            int o1, o2, o3;
            if (direction==0) {o1 = 1; o2 = 2; o3 = 3;}
            if (direction==1) {o1 = 2; o2 = 3; o3 = 1;}
            if (direction==2) {o1 = 3; o2 = 1; o3 = 2;}

            std::vector<double> testFluxArray(nFields*n_cells, dummyValue);

            // Fiducial Indices
            int const fiducialDensityIndex   = threadId;
            int const fiducialMomentumIndexX = threadId + n_cells * o1;
            int const fiducialMomentumIndexY = threadId + n_cells * o2;
            int const fiducialMomentumIndexZ = threadId + n_cells * o3;
            int const fiducialEnergyIndex    = threadId + n_cells * 4;
            int const fiducialMagneticYIndex = threadId + n_cells * (o2 + 4 + NSCALARS);
            int const fiducialMagneticZIndex = threadId + n_cells * (o3 + 4 + NSCALARS);

            _hlldInternal::_returnFluxes(threadId,
                                         o1,
                                         o2,
                                         o3,
                                         n_cells,
                                         testFluxArray.data(),
                                         densityFlux,
                                         momentumFluxX,
                                         momentumFluxY,
                                         momentumFluxZ,
                                         energyFlux,
                                         magneticFluxY,
                                         magneticFluxZ);

            // Find the indices for the various fields
            int densityLoc    = findIndex(testFluxArray, densityFlux,   fiducialDensityIndex,   "density");
            int momentumXLocX = findIndex(testFluxArray, momentumFluxX, fiducialMomentumIndexX,  "momentum X");
            int momentumYLocY = findIndex(testFluxArray, momentumFluxY, fiducialMomentumIndexY,  "momentum Y");
            int momentumZLocZ = findIndex(testFluxArray, momentumFluxZ, fiducialMomentumIndexZ,  "momentum Z");
            int energyLoc     = findIndex(testFluxArray, energyFlux,    fiducialEnergyIndex,    "energy");
            int magneticYLoc  = findIndex(testFluxArray, magneticFluxY, fiducialMagneticYIndex, "magnetic Y");
            int magneticZLoc  = findIndex(testFluxArray, magneticFluxZ, fiducialMagneticZIndex, "magnetic Z");

            for (size_t i = 0; i < testFluxArray.size(); i++)
            {
                // Skip the already checked indices
                if ((i != densityLoc)    and
                    (i != momentumXLocX) and
                    (i != momentumYLocY) and
                    (i != momentumZLocZ) and
                    (i != energyLoc)     and
                    (i != magneticYLoc)  and
                    (i != magneticZLoc))
                {
                    EXPECT_EQ(dummyValue, testFluxArray.at(i))
                        << "Unexpected value at index that _returnFluxes shouldn't be touching" << std::endl
                        << "Index     = " << i         << std::endl
                        << "Direction = " << direction << std::endl;
                }
            }
        }
    }
    // =========================================================================
#endif  // CUDA & HLLD