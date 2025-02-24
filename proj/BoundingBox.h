#ifndef _BOUNDINGBOX_H

#define _BOUNDINGBOX_H

#include <algorithm>
#include "vector3d.h"
#include "ray.h"
#include "constant.h"

namespace prender {

class BoundingBox {
public:
  inline BoundingBox(const Vector3d &minPos = Vector3d(0,0,0), const Vector3d &maxPos = Vector3d(0,0,0))
    : m_min(minPos)
    , m_max(maxPos)
    , m_centerPos((minPos+maxPos)/2.0)
  {
  }
  inline BoundingBox(const double minPos[3], const double maxPos[3]) 
  {
    SetBox(minPos, maxPos);
  }

  inline double CalcVolume() const 
  {
    return (m_max.x-m_min.x)*(m_max.y-m_min.y)*(m_max.z-m_min.z);
  }
  inline double CalcSurfaceArea() const 
  {
    double diff_x = m_max.x-m_min.x, diff_y = m_max.y-m_min.y, diff_z = m_max.z-m_min.z;
    return diff_x*diff_y + diff_x*diff_z + diff_y*diff_z;
  }

  template <typename FLOATING> static FLOATING CalcSurfaceArea(FLOATING min[3], FLOATING max[3]) 
  {
    FLOATING diff[3] = {max[0]-min[0], max[1]-min[1], max[2]-min[2]};
    return diff[0]*diff[1] + diff[1]*diff[2] + diff[0]*diff[2];
  }

  template <typename FLOATING> void SetBox(const FLOATING min[3], const FLOATING max[3]) 
  {
    SetBox(Vector3d(static_cast<double>(min[0]), static_cast<double>(min[1]), static_cast<double>(min[2])),
      Vector3d(static_cast<double>(max[0]), static_cast<double>(max[1]), static_cast<double>(max[2])));
  }

  void SetBox(const Vector3d &min_, const Vector3d &max_)
  {
    this->m_min = min_; this->m_max = max_;
    m_centerPos = (m_min+m_max)/2.0;
  }

  inline bool CheckIntersection(const Ray &ray, double &distance) const 
  {
    Vector3d t_min(PS_INF, PS_INF, PS_INF), t_max(-PS_INF, -PS_INF, -PS_INF);
    double fastest_out_t = PS_INF;
    double latest_in_t = -PS_INF;

    double min_array[3] = {m_min.x, m_min.y, m_min.z};
    double max_array[3] = {m_max.x, m_max.y, m_max.z};
    double ray_dir_array[3] = {ray.dir.x, ray.dir.y, ray.dir.z};
    double ray_start_array[3] = {ray.org.x, ray.org.y, ray.org.z};

    return CheckIntersection(ray_dir_array, ray_start_array, min_array, max_array, distance);
  }

  template <typename FLOATING> inline static bool CheckIntersection(const FLOATING rayDir[3], const FLOATING rayOrig[3], const FLOATING minbox[3], const FLOATING maxbox[3], FLOATING &distance) 
  {
    FLOATING fastest_out_t = PS_INF;
    FLOATING latest_in_t = -PS_INF;

    for (int i=0; i<3; i++) 
	{
      FLOATING t_min = PS_INF, t_max = -PS_INF;
      const FLOATING inv = rayDir[i] != 0 ? 1.0/rayDir[i] : static_cast<FLOATING>(PS_INF);
      const FLOATING t1 = (minbox[i] - rayOrig[i])*inv;
      const FLOATING t2 = (maxbox[i] - rayOrig[i])*inv;
      t_min = std::min(t1, t2);
      t_max = std::max(t1, t2);
      fastest_out_t = std::min(fastest_out_t, t_max);
      latest_in_t = std::max(latest_in_t, t_min);
      if (latest_in_t > fastest_out_t) return false;
    }

    if (latest_in_t > 0)
      distance = latest_in_t;// = ray.begin + ray.dir * latest_in_t;
    else
      distance = 0;

    return true;

  }

 inline bool CheckIntersection(const Ray& ray) const
 {
    double fastest_out_t = PS_INF;
    double latest_in_t = -PS_INF;

    {
      double t_min = PS_INF, t_max = -PS_INF;
      const double inv = ray.dir.x != 0 ? 1.0/ray.dir.x : PS_INF;
	  const double t1 = (m_min.x - ray.org.x)*inv;
      const double t2 = (m_max.x - ray.org.x)*inv;
      t_min = std::min(t1, t2);
      t_max = std::max(t1, t2);
      fastest_out_t = std::min(fastest_out_t, t_max);
      latest_in_t = std::max(latest_in_t, t_min);
      if (latest_in_t > fastest_out_t) return false;
    }
    {
      double t_min = PS_INF, t_max = -PS_INF;
      const double inv = ray.dir.y != 0 ? 1.0/ray.dir.y : PS_INF;
	  const double t1 = (m_min.y - ray.org.y)*inv;
      const double t2 = (m_max.y - ray.org.y)*inv;
      t_min = std::min(t1, t2);
      t_max = std::max(t1, t2);
      fastest_out_t = std::min(fastest_out_t, t_max);
      latest_in_t = std::max(latest_in_t, t_min);
      if (latest_in_t > fastest_out_t) return false;
    }
    {
      double t_min = PS_INF, t_max = -PS_INF;
      const double inv = ray.dir.z != 0 ? 1.0/ray.dir.z : PS_INF;
	  const double t1 = (m_min.z - ray.org.z)*inv;
      const double t2 = (m_max.z - ray.org.z)*inv;
      t_min = std::min(t1, t2);
      t_max = std::max(t1, t2);
      fastest_out_t = std::min(fastest_out_t, t_max);
      latest_in_t = std::max(latest_in_t, t_min);
      if (latest_in_t > fastest_out_t) return false;
    }
    return true;
  }

#if 10
  // Reference: http://d.hatena.ne.jp/ototoi/20090925/p1
  inline static bool CheckIntersection4floatAABB(
	  const __m128 bboxes[2][3], // 4boxes: min-max[2] * xyz[3] * boxes[4](__m128)
	  const __m128 rayOrig[3],   // ray origin
	  const __m128 rayInverseDir[3], // ray inversed dir
	  const int raySign[3],         // ray xyz direction => +:0, -:1
	  __m128 tmin, __m128 tmax,     // ray range tmin-tmax
	  bool results[4]               // intersection results
	  ) {

    // tmin = max(tmin, (box[min or max]-rayOrig)/rayDir)
    // tmax = min(tmax, (box[max or min]-rayOrig)/rayDir)
    // if tmax > tmin, intersects. Otherwise, no intersections

    // x
    tmin = _mm_max_ps(
      tmin, _mm_mul_ps(_mm_sub_ps(bboxes[raySign[0]][0], rayOrig[0]), rayInverseDir[0])
    );
    tmax = _mm_min_ps(
      tmax, _mm_mul_ps(_mm_sub_ps(bboxes[1 - raySign[0]][0], rayOrig[0]), rayInverseDir[0])
    );

    // y
    tmin = _mm_max_ps(
      tmin, _mm_mul_ps(_mm_sub_ps(bboxes[raySign[1]][1], rayOrig[1]), rayInverseDir[1])
      );
    tmax = _mm_min_ps(
      tmax, _mm_mul_ps(_mm_sub_ps(bboxes[1 - raySign[1]][1], rayOrig[1]), rayInverseDir[1])
      );

    // z
    tmin = _mm_max_ps(
      tmin, _mm_mul_ps(_mm_sub_ps(bboxes[raySign[2]][2], rayOrig[2]), rayInverseDir[2])
      );
    tmax = _mm_min_ps(
      tmax, _mm_mul_ps(_mm_sub_ps(bboxes[1 - raySign[2]][2], rayOrig[2]), rayInverseDir[2])
      );

    int ret = _mm_movemask_ps(_mm_cmpge_ps(tmax, tmin));
    if (ret == 0) return false;
    for (int i = 0; i < 4; i++) results[i] = ((ret >> i) & 0x1) != 0;
    return true;
  }

  inline static bool CheckIntersection2doubleAABB(
    const __m128d bboxes[2][3], // 4boxes: min-max[2] * xyz[3] * boxes[2](__m128)
    const __m128d rayOrig[3],   // ray origin
    const __m128d rayInverseDir[3], // ray inversed dir
    const int raySign[3],         // ray xyz direction => +:0, -:1
    __m128d tmin, __m128d tmax,     // ray range tmin-tmax
    bool results[2]               // intersection results
    ) {

    // tmin = max(tmin, (box[min or max]-rayOrig)/rayDir)
    // tmax = min(tmax, (box[max or min]-rayOrig)/rayDir)
    // if tmax > tmin, intersects. Otherwise, no intersections

    // x
    tmin = _mm_max_pd(
      tmin, _mm_mul_pd(_mm_sub_pd(bboxes[raySign[0]][0], rayOrig[0]), rayInverseDir[0])
      );
    tmax = _mm_min_pd(
      tmax, _mm_mul_pd(_mm_sub_pd(bboxes[1 - raySign[0]][0], rayOrig[0]), rayInverseDir[0])
      );

    // y
    tmin = _mm_max_pd(
      tmin, _mm_mul_pd(_mm_sub_pd(bboxes[raySign[1]][1], rayOrig[1]), rayInverseDir[1])
      );
    tmax = _mm_min_pd(
      tmax, _mm_mul_pd(_mm_sub_pd(bboxes[1 - raySign[1]][1], rayOrig[1]), rayInverseDir[1])
      );

    // z
    tmin = _mm_max_pd(
      tmin, _mm_mul_pd(_mm_sub_pd(bboxes[raySign[2]][2], rayOrig[2]), rayInverseDir[2])
      );
    tmax = _mm_min_pd(
      tmax, _mm_mul_pd(_mm_sub_pd(bboxes[1 - raySign[2]][2], rayOrig[2]), rayInverseDir[2])
      );

    int ret = _mm_movemask_pd(_mm_cmpge_pd(tmax, tmin));
    if (ret == 0) return false;
    for (int i = 0; i < 2; i++) results[i] = ((ret >> i) & 0x1) != 0;
    return true;
  }
#endif

  static BoundingBox CompoundBoxes(const BoundingBox &b1, const BoundingBox &b2) {
    return BoundingBox(Vector3d(std::min(b1.min().x, b2.min().x), std::min(b1.min().y, b2.min().y), std::min(b1.min().z, b2.min().z)), 
        Vector3d(std::max(b1.max().x, b2.max().x), std::max(b1.max().y, b2.max().y), std::max(b1.max().z, b2.max().z)));
  }

  inline void MergeAnotherBox(const BoundingBox &b2) 
  {
    m_min.x = std::min(m_min.x, b2.m_min.x);
    m_min.y = std::min(m_min.y, b2.m_min.y);
    m_min.z = std::min(m_min.z, b2.m_min.z);
    m_max.x = std::max(m_max.x, b2.m_max.x);
    m_max.y = std::max(m_max.y, b2.m_max.y);
    m_max.z = std::max(m_max.z, b2.m_max.z);
    m_centerPos = (m_min+m_max)/2;
  }

  const Vector3d &min() const {return m_min;}
  const Vector3d &max() const {return m_max;}
  const Vector3d &position() const {return m_centerPos;}

private:
  Vector3d m_min;
  Vector3d m_max;
  Vector3d m_centerPos;
};

}

#endif
