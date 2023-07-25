#ifndef __GEOMETRY_H_
#define __GEOMETRY_H_

#include <vector>
#include <iostream>


/*
Works only for meters in all dimensions.
*/

namespace geometry
{

    struct Point {
        Point() = default;
        Point(double x, double y, double z) : x{x}, y{y}, z{z} { }

        double x, y, z;

        inline void printPoint() const { std::cout << "(" << x << "," << y << "," << z << ")" << std::endl; }
    };

    struct Voxel2D {
        size_t count;
        double av_theta, av_density;
        double rms_theta, rms_density;
        std::vector<double> theta, energy;
        std::vector<size_t> theta_hist, density_hist;
    };

    class VoxelMap2D {
    private:
        uint16_t x_bins, y_bins;
        double precision;
        std::vector<Voxel2D> grid;

        const size_t getVoxelFromPoint(const Point& p) const;

    public:
        VoxelMap2D(uint16_t x_bins, uint16_t y_bins, double precision) :	x_bins{x_bins},
                                                                            y_bins{y_bins},
                                                                            grid{std::vector<Voxel2D>(x_bins*y_bins)},
                                                                            precision{precision} { }

        void event(const Point& hit, const double th, const double energy);

        void computeGridTheta();

        void print_count() const;
        void print_theta() const;
        void print_density() const;

        void output_map_count(std::string out_file) const;
        void output_map_theta(std::string out_file) const;
        void output_map_rms_theta(std::string out_file) const;
        void output_map_density(std::string out_file) const;
        void output_map_rms_density(std::string out_file) const;
        void output_rms(std::string out_file) const;
        void output_av(std::string out_file) const;
        void output_density(std::string out_file) const;
        void output_hist_theta(std::string out_file) const;
        void output_hist_density(std::string out_file) const;
        void output_theta_vals(std::string out_file);
        void output_max_theta(std::string out_file) const;
    };


    const double theta(const Point& p1, const Point& p2, const Point& p3, const Point& p4);
    const std::pair<Point,Point> lineIntersect3D(const Point&, const Point&, const Point&, const Point&);
    const std::pair<Point,Point> lineIntersect2D(const Point&, const Point&, const Point&, const Point&);
    double distance(const Point&, const Point&);
    const Point middlePoint(const std::pair<Point,Point>&);

}


#endif // __GEOMETRY_H_