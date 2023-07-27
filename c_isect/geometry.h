#ifndef C_TEST_GEOMETRY_H
#define C_TEST_GEOMETRY_H

#include <vector>
#include <iostream>
#include <cmath>

/*
Works only for meters in all dimensions.
*/

namespace geometry
{

    struct Point {
        Point() = default;
        Point(double x, double y, double z) : x{x}, y{y}, z{z} { }

        double x, y, z;
    };

    struct Voxel2D {
        double total_weight = 0;
        std::vector<double> theta, weight;
        double av_theta = std::nan("1");
        double conf_interval = std::nan("1");
    };

    struct Muon {
        std::vector<size_t> voxels;
        double theta;
    };

    class VoxelMap2D {
    private:
        uint16_t x_bins, y_bins;
        double precision;
        std::vector<Voxel2D> grid;

        //const size_t getVoxelFromPoint(const Point& p) const;
        void addWeightedToGrid(const Muon &muon);
        void addInitialToGrid(const Muon &muon);

    public:
        VoxelMap2D(uint16_t x_bins, uint16_t y_bins, double precision) :	x_bins{x_bins},
                                                                            y_bins{y_bins},
                                                                            grid{std::vector<Voxel2D>(x_bins*y_bins)},
                                                                            precision{precision} {}

        [[nodiscard]] std::vector<size_t> getVoxelsFromLine(std::pair<Point, Point>) const;
        void computeWeightedAvgGridTheta(const std::vector<geometry::Muon> &muons);
        void output_avg_theta(const std::string& out_file) const;
        [[nodiscard]] std::vector<double> getAvgGridThetas() const;
        void computeConfInterval();

        void computeInitialAvgGridTheta(const std::vector<geometry::Muon> &muons);

        const size_t getVoxelFromPoint(const Point &p) const;

        void addLineVoxel(std::vector<size_t> &lineVoxels, int x, int y) const;
    };


    double theta(const Point& p1, const Point& p2, const Point& p3, const Point& p4);
    const std::pair<Point,Point> lineIntersect3D(const Point&, const Point&, const Point&, const Point&);
    std::pair<Point,Point> lineIntersect2D(const Point&, const Point&, const Point&, const Point&);
    double distance(const Point&, const Point&);

    const Point middlePoint(const std::pair<Point,Point>&);

}

#endif //C_TEST_GEOMETRY_H
