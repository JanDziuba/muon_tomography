#include <cmath>

#include "csv.h"
#include "Geometry.h"


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

	size_t counter {1};
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

		const auto p_res = geometry::middlePoint(geometry::lineIntersect3D(p0_ab, p1_ab, p0_bw, p1_bw));
		const auto p_th = geometry::theta(p0_ab, p1_ab, p0_bw, p1_bw);

		if (!p_th || std::isnan(p_th))
			continue;
		grid_map.event(p_res, p_th, scintillator_above_1_E);


		++counter;
	}

	grid_map.computeGridTheta();

	// grid_map.output_map_count("0_grid_count.txt");
	// grid_map.output_map_rms_theta(argv[2]);
	// grid_map.output_map_rms_theta("0_grid_rms_theta.txt");
	// grid_map.output_map_rms_density("0_grid_rms_density.txt");
	// grid_map.output_rms(argv[2]);
	// grid_map.output_av(argv[2]);
	// grid_map.output_density(argv[2]);
	// grid_map.output_hist_theta(argv[2]);
	// grid_map.output_hist_density(argv[2]);
	grid_map.output_theta_vals(argv[2]);
	// grid_map.output_map_theta("1_grid_av_theta.txt");
	// grid_map.output_max_theta(argv[2]);

}
