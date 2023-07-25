#include "Geometry.h"

#include <cmath>
#include <utility>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <sstream>


namespace geometry {

    const size_t VoxelMap2D::getVoxelFromPoint(const Point &p) const {
        size_t x_idx = std::floor(p.x / precision + x_bins / 2.0);
        size_t y_idx = std::floor(p.y / precision + y_bins / 2.0);

        return (x_idx + y_idx * x_bins);
    }

    void VoxelMap2D::event(const Point &hit, const double th, const double e) {
//        if (hit.z < -0.2 || hit.z > 0.2)
//            return;

        const auto v_idx = getVoxelFromPoint(hit);
        if (v_idx >= x_bins * y_bins)
            return;

        ++grid[v_idx].count;
        grid[v_idx].theta.push_back(th);
        grid[v_idx].energy.push_back(e);
    }

    void VoxelMap2D::computeGridTheta() {
        std::for_each(grid.begin(), grid.end(), [](auto &cell) {

            // Density hist
            // cell.density_hist = std::vector<size_t>(100);
            // const double hist_bin = 0.05/100.;
            // for(size_t i{}; i<cell.count; ++i) {
            // 	double density = cell.theta[i] * cell.energy[i];
            // 	if (density >= 0.01 && density < 0.05) {
            // 		const size_t hist_idx = std::floor(density / hist_bin);
            // 		++cell.density_hist[hist_idx];
            // 	}
            // }

            // Theta hist
            // cell.theta_hist = std::vector<size_t>(200);
            // const double hist_bin_size = 0.3/200.;
            // std::for_each(cell.theta.cbegin(), cell.theta.cend(), [&](const auto& th){
            // 	if (th >= 0.3) {
            // 		return;
            // 	}
            // 	const size_t hist_idx = std::floor(th / hist_bin_size);
            // 	++cell.theta_hist[hist_idx];
            // });


            if (!cell.count) {
                // 	cell.rms_theta = 0.0;
                // 	cell.rms_density = 0.0;

                cell.av_theta = 0.0;
                // 	cell.av_density = 0.0;

                return;
            }

            // if (cell.count > 20) {
            // for(size_t i{}; i<cell.count; ++i) {
            // 	cell.rms_theta += std::pow(cell.theta[i], 2);
            // 	cell.rms_density += std::pow(cell.energy[i] * cell.theta[i] * 1000. * 10., 2);
            // }

            // cell.rms_theta /= cell.count;
            // cell.rms_theta = std::sqrt(cell.rms_theta);

            // cell.rms_density /= cell.count;
            // cell.rms_density = std::sqrt(cell.rms_density);
            // }


            const double th_sum = std::accumulate(cell.theta.cbegin(), cell.theta.cend(), 0.0);
            cell.av_theta = th_sum / double(cell.count);

            // double density = 0.0;
            // for(size_t i{}; i<cell.count; ++i)
            // 	density += cell.energy[i] * cell.theta[i];
            // cell.av_density = density / double(cell.count);

            // const double e_sum = std::accumulate(cell.energy.cbegin(), cell.energy.cend(), 0.0);
            // cell.av_density *= (e_sum / double(cell.energy.size())) * 1000. * 10.;


            // if (cell.count > 20) {
            // 	const double e_sum = std::accumulate(cell.energy.cbegin(), cell.energy.cend(), 0.0);
            // 	cell.av_density = (e_sum / double(cell.energy.size())) * cell.av_theta * 1000. * 10.;
            // }
            // else {
            // 	double d_max {};
            // 	for (size_t i{}; i < cell.theta.size(); ++i) {
            // 		if (double new_d_max = cell.theta[i] * cell.energy[i] * 1000; d_max < new_d_max)
            // 			d_max = new_d_max;
            // 		// d_sum += cell.theta[i] * cell.energy[i] * 1000;
            // 	}
            // 	// cell.density = d_sum / double(cell.theta.size());
            // 	cell.density = d_max;
            // }
        });
    }

    void VoxelMap2D::print_count() const {
        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j)
                std::cout << grid[j + i * y_bins].count << " ";
            std::cout << std::endl;
        }
    }

    void VoxelMap2D::print_theta() const {
        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j)
                std::cout << grid[j + i * y_bins].av_theta << " ";
            std::cout << std::endl;
        }
    }

    void VoxelMap2D::print_density() const {
        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j)
                std::cout << grid[j + i * y_bins].av_density << " ";
            std::cout << std::endl;
        }
    }

    void VoxelMap2D::output_map_count(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j)
                myfile << grid[j + i * y_bins].count << " ";
            myfile << std::endl;
        }

        myfile.close();
    }

    void VoxelMap2D::output_map_theta(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j) {
                // if (grid[j + i*y_bins].count < 20)
                // 	myfile << 0.0 << std::endl;
                // else
                myfile << grid[j + i * y_bins].av_theta << std::endl;
            }
        }

        myfile.close();
    }

    void VoxelMap2D::output_map_rms_theta(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j) {
                if (grid[j + i * y_bins].count < 20)
                    myfile << 0.0 << std::endl;
                else
                    myfile << grid[j + i * y_bins].rms_theta << std::endl;
            }
        }

        myfile.close();
    }

    void VoxelMap2D::output_map_density(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j)
                myfile << grid[j + i * y_bins].av_density << " ";
            myfile << std::endl;
        }

        myfile.close();
    }

    void VoxelMap2D::output_map_rms_density(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j) {
                if (grid[j + i * y_bins].count < 20)
                    myfile << 0.0 << std::endl;
                else
                    myfile << grid[j + i * y_bins].rms_density << std::endl;
            }
        }

        myfile.close();
    }

    void VoxelMap2D::output_rms(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i)
            for (size_t j{}; j < x_bins; ++j)
                myfile << grid[j + i * y_bins].count << "," << grid[j + i * y_bins].rms_theta << ","
                       << grid[j + i * y_bins].rms_density << std::endl;

        myfile.close();
    }

    void VoxelMap2D::output_av(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i)
            for (size_t j{}; j < x_bins; ++j)
                myfile << grid[j + i * y_bins].count << "," << grid[j + i * y_bins].av_theta << ","
                       << grid[j + i * y_bins].av_density << std::endl;

        myfile.close();
    }

    void VoxelMap2D::output_density(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i)
            for (size_t j{}; j < x_bins; ++j)
                myfile << grid[j + i * y_bins].count << "," << grid[j + i * y_bins].av_density << std::endl;

        myfile.close();
    }

    void VoxelMap2D::output_hist_theta(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j) {
                std::stringstream out;
                std::for_each(grid[j + i * y_bins].theta_hist.cbegin(), grid[j + i * y_bins].theta_hist.cend(),
                              [&](const auto &bin) {
                                  out << bin << ",";
                              });
                out.seekp(-1, out.cur);
                out << std::endl;
                myfile << out.rdbuf();
            }
        }

        myfile.close();
    }

    void VoxelMap2D::output_hist_density(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j) {
                std::stringstream out;
                std::for_each(grid[j + i * y_bins].density_hist.cbegin(), grid[j + i * y_bins].density_hist.cend(),
                              [&](const auto &bin) {
                                  out << bin << ",";
                              });
                out.seekp(-1, out.cur);
                out << std::endl;
                myfile << out.rdbuf();
            }
        }

        myfile.close();
    }

    void VoxelMap2D::output_theta_vals(std::string out_file) {
        std::ofstream myfile;
        myfile.open(out_file);

        size_t max_theta_size = 0;

        for (size_t y{}; y < y_bins; ++y) {
            for (size_t x{}; x < x_bins; ++x) {
                size_t theta_size = grid[x + y * x_bins].count;
                if (theta_size > max_theta_size) {
                    max_theta_size = theta_size;
                }
            }
        }

        const size_t COUNT_SIZE = max_theta_size;

        for (size_t y{}; y < y_bins; ++y) {
            for (size_t x{}; x < x_bins; ++x) {
                std::stringstream out;
                out << grid[x + y * x_bins].count << ",";


                if (grid[x + y * y_bins].count < COUNT_SIZE)
                    grid[x + y * y_bins].theta.resize(COUNT_SIZE);

                std::sort(grid[x + y * x_bins].theta.begin(), grid[x + y * x_bins].theta.end(), std::greater<double>());
                for (size_t th{}; th < COUNT_SIZE; ++th)
                    out << grid[x + y * x_bins].theta[th] << ",";

                out.seekp(-1, out.cur);
                out << std::endl;
                myfile << out.rdbuf();
            }
        }

        myfile.close();
    }

    void VoxelMap2D::output_max_theta(std::string out_file) const {
        std::ofstream myfile;
        myfile.open(out_file);

        for (size_t i{}; i < y_bins; ++i) {
            for (size_t j{}; j < x_bins; ++j) {
                if (!grid[j + i * y_bins].count)
                    myfile << 0.000001 << std::endl;
                else
                    myfile << *std::max_element(grid[j + i * y_bins].theta.cbegin(), grid[j + i * y_bins].theta.cend())
                           << std::endl;
            }
        }

        myfile.close();
    }


    const double theta(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {
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

    const std::pair<Point, Point> lineIntersect3D(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {
        Point p13, p43, p21;

        p13.x = p1.x - p3.x;
        p13.y = p1.y - p3.y;
        p13.z = p1.z - p3.z;
        p43.x = p4.x - p3.x;
        p43.y = p4.y - p3.y;
        p43.z = p4.z - p3.z;

        p21.x = p2.x - p1.x;
        p21.y = p2.y - p1.y;
        p21.z = p2.z - p1.z;

        const double d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
        const double d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
        const double d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
        const double d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
        const double d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

        const double denom = d2121 * d4343 - d4321 * d4321;
        const double numer = d1343 * d4321 - d1321 * d4343;

        const double mua = numer / denom;
        const double mub = (d1343 + d4321 * (mua)) / d4343;

        return std::make_pair(Point(p1.x + mua * p21.x, p1.y + mua * p21.y, p1.z + mua * p21.z),
                              Point(p3.x + mub * p43.x, p3.y + mub * p43.y, p3.z + mub * p43.z));
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
    const std::pair<Point, Point> lineIntersect2D(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {


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

}