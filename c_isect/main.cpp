#include <cmath>

#include "csv.h"
#include "geometry.h"


int main(int argc, char** argv)
{
    io::CSVReader<14> in(argv[1]);

    double 	scintillator_above_0_x,
            scintillator_above_0_y,
            scintillator_above_0_z,
            scintillator_above_1_x,
            scintillator_above_1_y,
            scintillator_above_1_z,
            scintillator_below_0_x,
            scintillator_below_0_y,
            scintillator_below_0_z,
            scintillator_below_1_x,
            scintillator_below_1_y,
            scintillator_below_1_z,
            scintillator_above_1_E,
            scintillator_below_0_E;

    geometry::VoxelMap2D grid_map(40, 40, 0.01);
    std::vector<geometry::Muon> muons;

    while(in.read_row(	scintillator_above_0_x,
                          scintillator_above_0_y,
                          scintillator_above_0_z,
                          scintillator_above_1_x,
                          scintillator_above_1_y,
                          scintillator_above_1_z,
                          scintillator_below_0_x,
                          scintillator_below_0_y,
                          scintillator_below_0_z,
                          scintillator_below_1_x,
                          scintillator_below_1_y,
                          scintillator_below_1_z,
                          scintillator_above_1_E,
                          scintillator_below_0_E)) {

        geometry::Point p0_ab(scintillator_above_0_x, scintillator_above_0_y, scintillator_above_0_z);
        geometry::Point p1_ab(scintillator_above_1_x, scintillator_above_1_y, scintillator_above_1_z);
        geometry::Point p0_bw(scintillator_below_0_x, scintillator_below_0_y, scintillator_below_0_z);
        geometry::Point p1_bw(scintillator_below_1_x, scintillator_below_1_y, scintillator_below_1_z);

        const std::pair<geometry::Point, geometry::Point> isect_pts =
                geometry::lineIntersect2D(p0_ab, p1_ab, p0_bw, p1_bw);

        const double p_th = geometry::theta(p0_ab, p1_ab, p0_bw, p1_bw);
        if (!p_th || std::isnan(p_th)) // TODO: why is this necessary?
            continue;

        const double p_th_adj = p_th / geometry::distance(isect_pts.first, isect_pts.second);

        std::vector<size_t> voxels = grid_map.getVoxelsFromLine(isect_pts);
//        printf("%lu\n", voxels.size());
//        std::vector<size_t> voxels = {grid_map.getVoxelFromPoint(geometry::middlePoint(isect_pts))};
        muons.push_back({voxels, p_th_adj});
    }

    grid_map.computeInitialAvgGridTheta(muons);
    std::vector<double> old_average_grid_thetas = grid_map.getAvgGridThetas();

    for (int i = 0; i < 1000; ++i) {
        grid_map.computeWeightedAvgGridTheta(muons);
        std::vector<double> new_average_grid_thetas = grid_map.getAvgGridThetas();
        // check if converged
        bool converged = true;
        for (int j = 0; j < old_average_grid_thetas.size(); ++j) {
            if (std::abs(old_average_grid_thetas[j] - new_average_grid_thetas[j]) > 0.0001) {
                converged = false;
                break;
            }
        }
        if (converged) {
            std::cout << "Converged after " << i+1 << " iterations." << std::endl;
            break;
        }
        old_average_grid_thetas = new_average_grid_thetas;
    }
    grid_map.computeConfInterval();
//    std::vector<double> avg_grid_thetas = grid_map.getAvgGridThetas();
//    double max_theta = 0;
//    // get max theta from avg_grid_thetas
//    for (int i = 0; i < avg_grid_thetas.size(); ++i) {
//        if (avg_grid_thetas[i] > max_theta)
//            max_theta = avg_grid_thetas[i];
//    }
//    std::cout << "Max theta: " << max_theta << std::endl;
    grid_map.output_avg_theta(argv[2]);
}
