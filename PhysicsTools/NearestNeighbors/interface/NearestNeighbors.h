#ifndef PhysicsTools_NearestNeighbors_NearestNeighbors_H__
#define PhysicsTools_NearestNeighbors_NearestNeighbors_H__

#include <memory>
#include <vector>

#include "nanoflann.hpp"

namespace cms::nanoflann {

  template <typename CoordinateType, int DIM = 3, typename IndexType = size_t, typename Distance = ::nanoflann::metric_L2>
  class PointCloud {
    typedef PointCloud<CoordinateType, DIM, IndexType, Distance> self_type;
    typedef typename Distance::template traits<CoordinateType, self_type>::distance_t metric_type;
    typedef ::nanoflann::KDTreeSingleIndexAdaptor<metric_type, self_type, DIM, IndexType> index_type;

  public:
    typedef std::array<CoordinateType, DIM> Point;
    typedef std::vector<Point> Points;

    PointCloud() {}
    PointCloud(const Points& points) : pts_(points) {}
    ~PointCloud() {}

    const Points& points() const { return pts_; }
    Points& points() { return pts_; }

    size_t size() const { return pts_.size(); }

    void add(const Point& point) { pts_.push_back(point); }

    template <typename OutputType>
    static std::vector<OutputType> knn(const Points& supports,
                                       const Points& queries,
                                       size_t num_neighbors,
                                       size_t leaf_max_size = 10) {
      PointCloud cloud(supports);
      index_type index(DIM, cloud, ::nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size));
      index.buildIndex();

      std::vector<OutputType> out_indices(queries.size() * num_neighbors);
      for (size_t i = 0; i < queries.size(); i++) {
        std::vector<IndexType> ids(num_neighbors);
        std::vector<CoordinateType> dists_sqr(num_neighbors);
        auto num_found = index.knnSearch(queries[i].data(), num_neighbors, &ids[0], &dists_sqr[0]);
        assert(num_found == num_neighbors);
        std::copy(ids.begin(), ids.end(), out_indices.begin() + (i * num_neighbors));
      }

      return out_indices;
    }

    // === interface needed by nanoflann ===
    inline size_t kdtree_get_point_count() const { return pts_.size(); }

    inline CoordinateType kdtree_get_pt(const size_t idx, const size_t dim) const {
      assert(dim >= 0 && dim < DIM);
      return pts_[idx][dim];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const {
      return false;
    }
    // =====================================

  private:
    Points pts_;
  };

}  // namespace cms::nanoflann

#endif  // PhysicsTools_NearestNeighbors_NearestNeighbors_H__
