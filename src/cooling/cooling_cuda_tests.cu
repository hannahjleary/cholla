/*!
* \file cooling_cuda_tests.cpp
* \author Helena Richie (helenarichie@pitt.edu)
* \brief Testing the CIE cooling rate function in cooling_cuda.cu
*
*/

// STL Includes
#include <iostream>
#include <string>

// External Includes
#include <gtest/gtest.h> // Include GoogleTest and related libraries/headers

// Local Includes
#include "../global/global_cuda.h"
#include "../utils/gpu.hpp"
#include "../utils/testing_utilities.h"
#include "../cooling/cooling_cuda.h" // Include code to test

#ifdef COOLING_GPU

TEST(tCOOLINGPracticeTest, PracticeTestExpectCorrectOutput) // test suite name, test name
{
    Real const testn = 5;
    Real const testT = 5;
    Real const testNumber = CIE_cool(testn, testT);

    Real const fiducialNumber = 100;

    Real absoluteDiff;
    Real ulpsDiff;

    bool isTrue;

    isTrue = testingUtilities::nearlyEqualDbl(fiducialNumber, testNumber, absoluteDiff, ulpsDiff);
    
    EXPECT_TRUE(isTrue)
            << "The fiducial value is:       " << fiducialNumber  << std::endl
            << "The test value is:           " << testNumber      << std::endl
            << "The absolute difference is:  " << absoluteDiff    << std::endl
            << "The ULP difference is:       " << ulpsDiff        << std::endl;
}

#endif // COOLING