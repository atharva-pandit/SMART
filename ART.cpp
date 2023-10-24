#include "basic_func.hpp"
#include "imageOpt.hpp"

constexpr int maxvox {500};		//number of voxels in each dimension
auto voxgrid = new int[maxvox][maxvox][maxvox];	//3D voxel grid maxvox*maxvox*maxvox
auto labels = new int[maxvox][maxvox][maxvox]; //3D voxel grid to store local maxima
float scale {maxvox/100.0};	// voxels per mm; maxvox/scale = extent of grid in mm

#include "ART.hpp"
#include "boxhit.hpp"
#include "label3D.hpp"

int main() {

Timer time;

constexpr int iterations{10}, perc{2};
constexpr bool do_ART {0}, do_MART {0}, do_SBP {1};
constexpr float tracer_rad {4.35327};
int labelcount {};

std::vector<CoordPixel>	pts0,pts1,pts2,pts3,pts4,pts5,pts6;
std::vector<CoordFloat> angle_deg, angle_rad;
angle_deg = readCoord<CoordFloat>("datafiles/angles.dat");
angle_rad = DegToRad(angle_deg);

Eigen::RowVector3f rs, l;
rs << (Eigen::RowVector3f() << -690,0,0).finished();

CoordIntGrey ptline;
CoordFloat Hit1, Hit2, B1{-50.0,-50.0,-50.0}, B2{49.9,49.9,49.9}, lit;

for (int t=0; t<1; t++) {

	std::string outfile {"datafiles/centres/StaticPacking/SB800_p"+std::to_string(perc)+"_t"+std::to_string(t)+".dat"};
	clearfile(outfile);

	std::vector<std::string> imagenames = generate_fnames(angle_deg, t, perc);

	pts0 = get_global(binarise(readImage(imagenames[0].c_str()), "minimum"), angle_rad[0]);
	pts1 = get_global(binarise(readImage(imagenames[1].c_str()), "minimum"), angle_rad[1]);
	pts2 = get_global(binarise(readImage(imagenames[2].c_str()), "minimum"), angle_rad[2]);
	pts3 = get_global(binarise(readImage(imagenames[3].c_str()), "minimum"), angle_rad[3]);
	pts4 = get_global(binarise(readImage(imagenames[4].c_str()), "minimum") ,angle_rad[4]);
	pts5 = get_global(binarise(readImage(imagenames[5].c_str()), "minimum"), angle_rad[5]);
	pts6 = get_global(binarise(readImage(imagenames[6].c_str()), "minimum"), angle_rad[6]);

	std::vector<CoordIntGrey> Ray0, Ray1, Ray2, Ray3, Ray4, Ray5, Ray6;

	l = rs*rotmat(angle_rad[0].x,angle_rad[0].y,angle_rad[0].z,false);
	lit.x = l[0]; lit.y = l[1]; lit.z = l[2];
	for (CoordPixel& pt : pts0) {
		if (pt.label != -1) {
			if ( CheckLineBox(B1, B2, lit, pt, Hit1, Hit2) ) {
				ptline.raypts = bresenham(Hit1,Hit2);
				ptline.greyval = pt.greyval;
				Ray0.emplace_back(ptline);
			}
			else {std::cout << "0 missed" << '\n';}
		}
	}

	l = rs*rotmat(angle_rad[1].x,angle_rad[1].y,angle_rad[1].z,false);
	lit.x = l[0]; lit.y = l[1]; lit.z = l[2];
	for (CoordPixel& pt : pts1) {
		if (pt.label != -1) {
			if ( CheckLineBox(B1, B2, lit, pt, Hit1, Hit2) ) {
				ptline.raypts = bresenham(Hit1,Hit2);
				ptline.greyval = pt.greyval;
				Ray1.emplace_back(ptline);
			}
			else {std::cout << "1 missed" << '\n';}
		}
	}
	
	l = rs*rotmat(angle_rad[2].x,angle_rad[2].y,angle_rad[2].z,false);
	lit.x = l[0]; lit.y = l[1]; lit.z = l[2];
	for (CoordPixel& pt : pts2) {
		if (pt.label != -1) {
			if ( CheckLineBox(B1, B2, lit, pt, Hit1, Hit2) ) {
				ptline.raypts = bresenham(Hit1,Hit2);
				ptline.greyval = pt.greyval;
				Ray2.emplace_back(ptline);
			}
			else {std::cout << "2 missed" << '\n';}
		}
	}

	l = rs*rotmat(angle_rad[3].x,angle_rad[3].y,angle_rad[3].z,false);
	lit.x = l[0]; lit.y = l[1]; lit.z = l[2];
	for (CoordPixel& pt : pts3) {
		if (pt.label != -1) {
			if ( CheckLineBox(B1, B2, lit, pt, Hit1, Hit2) ) {
				ptline.raypts = bresenham(Hit1,Hit2);
				ptline.greyval = pt.greyval;
				Ray3.emplace_back(ptline);
			}
			else {std::cout << "3 missed" << '\n';}
		}
	}

	l = rs*rotmat(angle_rad[4].x,angle_rad[4].y,angle_rad[4].z,false);
	lit.x = l[0]; lit.y = l[1]; lit.z = l[2];
	for (CoordPixel& pt : pts4) {
		if (pt.label != -1) {
			if ( CheckLineBox(B1, B2, lit, pt, Hit1, Hit2) ) {
				ptline.raypts = bresenham(Hit1,Hit2);
				ptline.greyval = pt.greyval;
				Ray4.emplace_back(ptline);
			}
			else {std::cout << "4 missed" << '\n';}
		}
	}

	l = rs*rotmat(angle_rad[5].x,angle_rad[5].y,angle_rad[5].z,false);
	lit.x = l[0]; lit.y = l[1]; lit.z = l[2];
	for (CoordPixel& pt : pts5) {
		if (pt.label != -1) {
			if ( CheckLineBox(B1, B2, lit, pt, Hit1, Hit2) ) {
				ptline.raypts = bresenham(Hit1,Hit2);
				ptline.greyval = pt.greyval;
				Ray5.emplace_back(ptline);
			}
			else {std::cout << "5 missed" << '\n';}
		}
	}

	l = rs*rotmat(angle_rad[6].x,angle_rad[6].y,angle_rad[6].z,false);
	lit.x = l[0]; lit.y = l[1]; lit.z = l[2];
	for (CoordPixel& pt : pts6) {
		if (pt.label != -1) {
			if ( CheckLineBox(B1, B2, lit, pt, Hit1, Hit2) ) {
				ptline.raypts = bresenham(Hit1,Hit2);
				ptline.greyval = pt.greyval;
				Ray6.emplace_back(ptline);
			}
			else {std::cout << "6 missed" << '\n';}
		}
	}

	// set all voxel values to 0
	clear_voxgrid();

	if (do_ART) {
		for (int i=0; i<iterations; i++) {
			for (CoordIntGrey& rayNum : Ray0) ART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray1) ART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray2) ART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray3) ART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray4) ART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray5) ART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray6) ART(rayNum.raypts, rayNum.greyval);
		}
	}
	else if (do_MART) {
		// set all voxel values to 1
		setupVoxgrid_MART();
		for (int i=0; i<iterations; i++) {
			for (CoordIntGrey& rayNum : Ray0) MART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray1) MART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray2) MART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray3) MART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray4) MART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray5) MART(rayNum.raypts, rayNum.greyval);
			for (CoordIntGrey& rayNum : Ray6) MART(rayNum.raypts, rayNum.greyval);
		}
	}

	else if (do_SBP) {
		for (CoordIntGrey& rayNum : Ray0) simpback(rayNum.raypts);
		for (CoordIntGrey& rayNum : Ray1) simpback(rayNum.raypts);
		for (CoordIntGrey& rayNum : Ray2) simpback(rayNum.raypts);
		for (CoordIntGrey& rayNum : Ray3) simpback(rayNum.raypts);
		for (CoordIntGrey& rayNum : Ray4) simpback(rayNum.raypts);
		for (CoordIntGrey& rayNum : Ray5) simpback(rayNum.raypts);
		for (CoordIntGrey& rayNum : Ray6) simpback(rayNum.raypts);
	}

	std::cout << "Frame " << t << " done." << '\n';

initialise_labelgrid( voxgrid_max() );
std::map<int,int> label_cnt, label_list, big_blobs;
std::map<int,int>::iterator it;

int max_blob_size = static_cast<int>(1.3333*3.1415*std::pow(tracer_rad*scale,3));
int min_blob_size = static_cast<int>(max_blob_size/10.0);

std::cout << "here 1" << '\n';
label_cnt = add_labels(labelcount);

for (it = label_cnt.begin(); it != label_cnt.end(); it++) {
	//if (it->second > max_blob_size) big_blobs.insert({it->first, it->second});
	if (it->second > min_blob_size) label_list.insert({it->first, it->second});
}

//CoordInt pointpoint;
//std::vector<CoordId> distmap;
//std::string filefile {"points_dist.dat"};
//clearfile(filefile);

//std::cout << "here 2" << '\n';
////while (!big_blobs.empty()) {
//	for (it = big_blobs.begin(); it != big_blobs.end(); it++) {
//		std::vector<CoordInt> pointList;
//		for (int i=0; i<maxvox; i++) {
//			for (int j=0; j<maxvox; j++) {
//				for (int k=0; k<maxvox; k++) {
//					if ( labels[i][j][k] == it->first ) {
//						pointpoint.x = i;
//						pointpoint.y = j;
//						pointpoint.z = k;
//						pointList.emplace_back(pointpoint);
//					}
//				}
//			}
//		}
//		distmap = distanceMap(pointList);
//		for (CoordId& pt : distmap) label_cnt = add_labels(labelcount);
//		//std::map<int,int> big_blobs.clear();
//		for (it = label_cnt.begin(); it != label_cnt.end(); it++) {
//			//if (it->second > max_blob_size/2) {big_blobs.insert({it->first, it->second}); std::cout << "help\n";}
//			if (it->second > min_blob_size/5) label_list.insert({it->first, it->second});
//		}
//	}
////}

std::cout << "here 3" << '\n';
//writeCoordId(filefile, distmap);
std::vector<CoordFloat> centres = find_cluster_centres(label_list);
writeCoord(outfile, centres);
delete_voxgrid();

std::cout << "Time taken: " << time.elapsed() << " seconds\n";
}

return 0;
}
