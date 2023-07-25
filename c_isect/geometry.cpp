#include "geometry.h"

#include <cmath>
#include <utility>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <vector>

namespace geometry {

    void VoxelMap2D::addLineVoxel(std::vector<size_t>& lineVoxels, int x, int y) const {
        if (x >= 0 && y >= 0 && x + y * x_bins < x_bins * y_bins)
            lineVoxels.push_back((x + y * x_bins));
    }

    // get voxels from line using Bresenham's algorithm
    std::vector<size_t> VoxelMap2D::getVoxelsFromLine(const std::pair<Point, Point> line) const {
        std::vector<size_t> voxels;

        // if x-y line is too long, the angle is too big meaning bad sample
        if(std::sqrt(std::pow(line.first.x - line.second.x, 2) +
                  std::pow(line.first.y - line.second.y, 2)) > 0.1) {
            return voxels;
        }

        double x1 = line.first.x / precision + x_bins / 2.0;
        double y1 = line.first.y / precision + y_bins / 2.0;
        double x2 = line.second.x / precision + x_bins / 2.0;
        double y2 = line.second.y / precision + y_bins / 2.0;
        int incx = (x2 > x1) ? 1 : -1;
        int incy = (y2 > y1) ? 1 : -1;
        double dx = x2 - x1;
        double dy = y2 - y1;
        dx=std::abs(dx);
        dy=std::abs(dy);
        int x = static_cast <int> (std::floor(x1));
        int y = static_cast <int> (std::floor(y1));
        int x_final = static_cast <int> (std::floor(x2));
        int y_final = static_cast <int> (std::floor(y2));
        double e;
        double inc1, inc2;
        addLineVoxel(voxels, x, y);
        if (dx > dy) {
            e = 2 * dy - dx;
            inc1 = 2 * (dy - dx);
            inc2 = 2 * dy;
            while (x != x_final) {
                if (e >= 0.5) {
                    y += incy;
                    e += inc1;
                } else
                    e += inc2;
                x += incx;
                addLineVoxel(voxels, x, y);
            }
        } else {
            e = 2 * dx - dy;
            inc1 = 2 * (dx - dy);
            inc2 = 2 * dx;
            while (y != y_final)  {
                if (e >= 0.5) {
                    x += incx;
                    e += inc1;
                } else
                    e += inc2;
                y += incy;
                addLineVoxel(voxels, x, y);
            }
        }
        return voxels;
    }

    void VoxelMap2D::addInitialToGrid(const geometry::Muon& muon) {
        std::vector<size_t> voxels = muon.voxels;
        double theta = muon.theta;

        size_t voxels_count = voxels.size();
        double weight = 1.0 / voxels_count;
        for (const auto& v : voxels) {
            grid[v].total_weight += weight;
            grid[v].weight.push_back(weight);
            grid[v].theta.push_back(theta);
        }
    }


    void VoxelMap2D::addWeightedToGrid(const geometry::Muon& muon) {
        std::vector<size_t> voxels = muon.voxels;
        double theta = muon.theta;

        double av_theta_sum = 0;
        for (const auto& v : voxels) {
            av_theta_sum += grid[v].av_theta;
        }

        size_t voxels_count = voxels.size();
        for (const auto& v : voxels) {
            // if many voxels then line goes through many materials, so adjust theta based on av_theta
            if (voxels_count > 3) {
                double theta_adj = theta * (grid[v].av_theta * voxels_count / av_theta_sum);
                grid[v].theta.push_back(theta_adj);
            } else {
                grid[v].theta.push_back(theta);
            }
        }
    }

    void VoxelMap2D::computeInitialAvgGridTheta(const std::vector<geometry::Muon> &muons) {
        for (const geometry::Muon& muon : muons) {
            addInitialToGrid(muon);
        }

        for (auto& v : grid) {
            if (v.total_weight == 0){
                v.av_theta = 0;
                continue;
            }
            double sum_theta = 0;
            for (size_t i{}; i < v.theta.size(); ++i) {
                sum_theta += v.theta[i] * v.weight[i];
            }
            v.av_theta = sum_theta / v.total_weight;
        }
    }

    void VoxelMap2D::computeWeightedAvgGridTheta(const std::vector<geometry::Muon> &muons) {
        for (auto& v : grid) {
            v.theta.clear();
        }

        for (const geometry::Muon& muon : muons) {
            addWeightedToGrid(muon);
        }

        for (auto& v : grid) {
            if (v.total_weight == 0){
                v.av_theta = 0;
                continue;
            }
            double sum_theta = 0;
            for (size_t i{}; i < v.theta.size(); ++i) {
                sum_theta += v.theta[i] * v.weight[i];
            }
            v.av_theta = sum_theta / v.total_weight;
        }
    }

    std::vector<double> VoxelMap2D::getAvgGridThetas() const {
        std::vector<double> avg_grid_thetas;
        avg_grid_thetas.reserve(grid.size());
        for (const auto& v : grid) {
            avg_grid_thetas.push_back(v.av_theta);
        }
        return avg_grid_thetas;
    }

    void VoxelMap2D::computeConfInterval(){
        for (auto& v : grid) {
            if (v.total_weight == 0) {
                v.conf_interval = 0; // later set to max*10
                continue;
            }
            double sum_variance = 0;
            for (size_t i{}; i < v.theta.size(); ++i) {
                sum_variance += (v.av_theta - v.theta[i]) * (v.av_theta - v.theta[i]) * v.weight[i];
            }
            v.conf_interval = sum_variance / v.total_weight;
        }

        double max_conf_interval = 0;
        for (auto& v : grid) {
            if (v.conf_interval > max_conf_interval)
                max_conf_interval = v.conf_interval;
        }

        for (auto& v : grid) {
            if (v.total_weight == 0) // if no events in voxel set conf interval to max*10
                v.conf_interval = max_conf_interval*10;
            else if (v.total_weight <= 1) // if one or fewer events in voxel set conf interval to max
                v.conf_interval = max_conf_interval;
        }
    }


    void VoxelMap2D::output_avg_theta(const std::string& out_file) const{
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t y{}; y < y_bins; ++y) {
            for (size_t x{}; x < x_bins; ++x) {
                myfile << grid[x + y * x_bins].av_theta << "," << grid[x + y * x_bins].conf_interval << std::endl;
            }
        }

        myfile.close();
    }

    double theta(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {
        double v_p21_x, v_p21_y, v_p21_z;
        double v_p43_x, v_p43_y, v_p43_z;

        v_p21_x = p2.x - p1.x;
        v_p21_y = p2.y - p1.y;
        v_p21_z = p2.z - p1.z;

        v_p43_x = p4.x - p3.x;
        v_p43_y = p4.y - p3.y;
        v_p43_z = p4.z - p3.z;

        double dot = v_p21_x * v_p43_x + v_p21_y * v_p43_y + v_p21_z * v_p43_z;
        double mag_p21 = std::sqrt(v_p21_x * v_p21_x + v_p21_y * v_p21_y + v_p21_z * v_p21_z);
        double mag_p43 = std::sqrt(v_p43_x * v_p43_x + v_p43_y * v_p43_y + v_p43_z * v_p43_z);

        double theta = std::acos(dot / (mag_p21 * mag_p43));

        return theta;
    }

    void sub_v3_v3v3(double r[3], const double a[3], const double b[3]) {
        r[0] = a[0] - b[0];
        r[1] = a[1] - b[1];
        r[2] = a[2] - b[2];
    }

    void madd_v3_v3v3fl(double r[3], const double a[3], const double b[3], double f) {
        r[0] = a[0] + b[0] * f;
        r[1] = a[1] + b[1] * f;
        r[2] = a[2] + b[2] * f;
    }

    double dot_v3v3(const double a[3], const double b[3]) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    bool isect_line_plane_v3(double r_isect_co[3],
                             const double l1[3],
                             const double l2[3],
                             const double plane_co[3],
                             const double plane_no[3]) {
        double u[3], h[3];
        double dot;

        sub_v3_v3v3(u, l2, l1);
        sub_v3_v3v3(h, l1, plane_co);
        dot = dot_v3v3(plane_no, u);

        if (fabsf(dot) > 1e-6) {
            double lambda = -dot_v3v3(plane_no, h) / dot;
            madd_v3_v3v3fl(r_isect_co, l1, u, lambda);
            return true;
        }

        /* The segment is parallel to plane */
        return false;
    }

    // get intersections between first line and top of the plane with materials
    // and second line and bottom of the plane with materials
    std::pair<Point, Point> lineIntersect2D(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {


        const double plane_co1[3] = {0, 0, 0.05}; // coords of top of the plane
        const double plane_co2[3] = {0, 0, -0.05}; // coords of bottom of the plane

        const double plane_no[3] = {0, 0, 1};
        double r_isect_co1[3];
        double r_isect_co2[3];

        const double l1[3] = {p1.x, p1.y, p1.z};
        const double l2[3] = {p2.x, p2.y, p2.z};
        const double l3[3] = {p3.x, p3.y, p3.z};
        const double l4[3] = {p4.x, p4.y, p4.z};

        isect_line_plane_v3(r_isect_co1, l1, l2, plane_co1, plane_no);
        isect_line_plane_v3(r_isect_co2, l3, l4, plane_co2, plane_no);

        return std::make_pair(Point(r_isect_co1[0], r_isect_co1[1], r_isect_co1[2]),
                              Point(r_isect_co2[0], r_isect_co2[1], r_isect_co2[2]));
    }


    double distance(const Point &p1, const Point &p2) {
        return std::sqrt(std::pow(p1.x - p2.x, 2) +
                         std::pow(p1.y - p2.y, 2) +
                         std::pow(p1.z - p2.z, 2));
    }

    const Point middlePoint(const std::pair<Point, Point> &points) {
        Point p_res;
        p_res.x = (std::get<0>(points).x + std::get<1>(points).x) / 2;
        p_res.y = (std::get<0>(points).y + std::get<1>(points).y) / 2;
        p_res.z = (std::get<0>(points).z + std::get<1>(points).z) / 2;

        return p_res;
    }

    const size_t VoxelMap2D::getVoxelFromPoint(const Point &p) const {
        size_t x_idx = std::floor(p.x / precision + x_bins / 2.0);
        size_t y_idx = std::floor(p.y / precision + y_bins / 2.0);

        return (x_idx + y_idx * x_bins);
    }

}