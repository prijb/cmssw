#include <cppunit/extensions/HelperMacros.h>

#include "PhysicsTools/NearestNeighbors/interface/NearestNeighbors.h"

#include <iostream>

using namespace cms::nanoflann;

class testNanoFlann : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testNanoFlann);
  CPPUNIT_TEST(checkAll);
  CPPUNIT_TEST_SUITE_END();

public:
  void checkAll();
};

CPPUNIT_TEST_SUITE_REGISTRATION(testNanoFlann);

void testNanoFlann::checkAll() {
  const size_t num_neighbors = 3;

  PointCloud<float>::Points points{
      {{0., 0., 0.}},
      {{0., 1., 0.}},
      {{0., 0., 1.}},
      {{0., 1., 1.}},
      {{0., 0., 10.}},
  };

  PointCloud<float>::Points queries{
      {{0., 1., 0.}},
      {{0., 0., 1.}},
  };

  std::vector<float> expected{1, 0, 3, 2, 0, 3};

  std::vector<float> output;
  CPPUNIT_ASSERT_NO_THROW(output = PointCloud<float>::knn<float>(points, queries, num_neighbors));
  CPPUNIT_ASSERT(output.size() == expected.size());
  for (unsigned i = 0; i < output.size(); ++i) {
    CPPUNIT_ASSERT(output[i] == expected[i]);
    if (i % num_neighbors == 0)
      std::cout << std::endl;
    std::cout << output[i] << ", ";
  }
  std::cout << std::endl;
}
