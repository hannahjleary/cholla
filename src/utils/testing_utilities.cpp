/*!
 * \file testing_utilites.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Implementation file for various testing related utility functions
 *
 */

// STL includes
#include <limits>
#include <cmath>
#include <cstring>
#include <iostream>

// External Includes
#include <gtest/gtest.h>    // Include GoogleTest and related libraries/headers

// Local includes
#include "../utils/testing_utilities.h" // Include the header file
#include "../system_tests/system_tester.h" // provide systemTest class

namespace testingUtilities
{
    // =========================================================================
    int64_t ulpsDistanceDbl(double const &a, double const &b)
    {
        // Save work if the floats are equal.
        // Also handles +0 == -0
        if (a == b) return 0;

        const auto maxInt = std::numeric_limits<int64_t>::max();

        // If either one is NaN then they are not equal, max distance.
        if (std::isnan(a) || std::isnan(b)) return maxInt;

        // If one's infinite and they're not equal, max distance.
        if (std::isinf(a) || std::isinf(b)) return maxInt;

        int64_t ia, ib;
        std::memcpy(&ia, &a, sizeof(double));
        std::memcpy(&ib, &b, sizeof(double));

        // Don't compare differently-signed floats.
        if ((ia < 0) != (ib < 0)) return maxInt;

        // Return the absolute value of the distance in ULPs.
        int64_t distance = ia - ib;
        if (distance < 0) distance = -distance;

        return distance;
    }
    // =========================================================================

    // =========================================================================
    bool nearlyEqualDbl(double  const &a,
                        double  const &b,
                        double  &absoluteDiff,
                        int64_t &ulpsDiff,
                        double  const &fixedEpsilon, // = 1E-14 by default
                        int     const &ulpsEpsilon)  // = 4 by default
    {
        // Compute differences
        ulpsDiff = ulpsDistanceDbl(a, b);
        absoluteDiff = std::abs(a - b);

        // Perform the ULP check which is for numbers far from zero
        if (ulpsDiff <= ulpsEpsilon)
        {
            return true;
        }
        // Perform the absolute check which is for numbers near zero
        else if (absoluteDiff <= fixedEpsilon)
        {
            return true;
        }
        // if none of the checks have passed indicate test failure
        else
        {
            return false;
        }
    }
    // =========================================================================

    // =========================================================================
    void checkResults(double fiducialNumber,
                      double testNumber,
                      std::string outString,
                      double fixedEpsilon,
                      int ulpsEpsilon)
    {
        // Check for equality and if not equal return difference
        double absoluteDiff;
        int64_t ulpsDiff;
        bool areEqual;

        if ((fixedEpsilon < 0) and (ulpsEpsilon < 0))
        {
            areEqual = testingUtilities::nearlyEqualDbl(fiducialNumber,
                                                        testNumber,
                                                        absoluteDiff,
                                                        ulpsDiff);
        }
        else if ((fixedEpsilon > 0) and (ulpsEpsilon < 0))
        {
            areEqual = testingUtilities::nearlyEqualDbl(fiducialNumber,
                                                        testNumber,
                                                        absoluteDiff,
                                                        ulpsDiff,
                                                        fixedEpsilon);
        }
        else
        {
            areEqual = testingUtilities::nearlyEqualDbl(fiducialNumber,
                                                        testNumber,
                                                        absoluteDiff,
                                                        ulpsDiff,
                                                        fixedEpsilon,
                                                        ulpsEpsilon);
        }

        EXPECT_TRUE(areEqual)
            << "Difference in "                << outString       << std::endl
            << "The fiducial value is:       " << fiducialNumber  << std::endl
            << "The test value is:           " << testNumber      << std::endl
            << "The absolute difference is:  " << absoluteDiff    << std::endl
            << "The ULP difference is:       " << ulpsDiff        << std::endl;
    }
    // =========================================================================

  void wrapperEqual(int i, int j, int k, std::string dataSetName,
		    double test_value, double fid_value, double fixedEpsilon=5.0E-12) {

    std::string outString;
    outString += dataSetName;
    outString += " dataset at [";
    outString += i;
    outString += ",";
    outString += j;
    outString += ",";
    outString += k;
    outString += "]";

    checkResults(fid_value,test_value,outString,fixedEpsilon);
  }

  void analyticConstant(systemTest::SystemTestRunner testObject, std::string dataSetName, double value) {
    std::vector<size_t> testDims(3,1);
    std::vector<double> testData = testObject.loadTestFieldData(dataSetName,testDims);
    for (size_t i = 0; i < testDims[0]; i++)
      {
	for (size_t j = 0; j < testDims[1]; j++)
	  {
	    for (size_t k = 0; k < testDims[2]; k++)
	      {
		size_t index = (i * testDims[1] * testDims[2]) + (j * testDims[2]) + k;

		wrapperEqual(i,j,k,dataSetName,testData.at(index),value);
	      }
	  }
      }
  }

  void analyticSine(systemTest::SystemTestRunner testObject, std::string dataSetName,
		    double constant, double amplitude,
		    double kx, double ky, double kz, double phase, double tolerance)
  {
    std::vector<size_t> testDims(3,1);
    std::vector<double> testData = testObject.loadTestFieldData(dataSetName,testDims);
    for (size_t i = 0; i < testDims[0]; i++)
      {
	for (size_t j = 0; j < testDims[1]; j++)
	  {
	    for (size_t k = 0; k < testDims[2]; k++)
	      {
		double value = constant + amplitude*std::sin(kx*i+ky*j+kz*k+phase);
		size_t index = (i * testDims[1] * testDims[2]) + (j * testDims[2]) + k;
		wrapperEqual(i,j,k,dataSetName,testData.at(index),value,tolerance);
	      }
	  }
      }
  }




}
