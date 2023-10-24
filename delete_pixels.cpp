#include "basic_func.hpp"
#include "imageOpt.hpp"
#include <map>

/* Check if a 3D line segment intersects a sphere ----------------------------- */
bool lineptdistcheck(CoordFloat& lit, CoordPixel& px, CoordFloat& cent, const float& trac_rad)
{
	CoordFloat x0x1 {cent.x - lit.x, cent.y - lit.y, cent.z - lit.z};
	CoordFloat x0x2 {cent.x - px.x, cent.y - px.y, cent.z - px.z};
	CoordFloat x2x1 {px.x - lit.x, px.y - lit.y, px.z - lit.z};

	if ( (normSq(cross_pdt(x0x1, x0x2)) / normSq(x2x1)) <= std::pow(trac_rad,2) ) return 1;
	else return 0;
}

/* Label connected pixels --------------------------------------------- */
void label_connected_components(int px, std::vector<CoordPixel>& image, int width)
{
	if (image[px+1].label == 0) {
		image[px+1].label = image[px].label;
		label_connected_components(px+1, image, width);
	}

	if (image[px-1].label == 0) {
		image[px-1].label = image[px].label;
		label_connected_components(px-1, image, width);
	}

	if (image[px+width].label == 0) {
		image[px+width].label = image[px].label;
		label_connected_components(px+width, image, width);
	}

	if (image[px-width].label == 0) {
		image[px-width].label = image[px].label;
		label_connected_components(px-width, image, width);
	}

	if (image[px+width+1].label == 0) {
		image[px+width+1].label = image[px].label;
		label_connected_components(px+width+1, image, width);
	}

	if (image[px+width-1].label == 0) {
		image[px+width-1].label = image[px].label;
		label_connected_components(px+width-1, image, width);
	}

	if (image[px-width+1].label == 0) {
		image[px-width+1].label = image[px].label;
		label_connected_components(px-width+1, image, width);
	}

	if (image[px-width-1].label == 0) {
		image[px-width-1].label = image[px].label;
		label_connected_components(px-width-1, image, width);
	}
}

/* Restore pixels to image ---------------------------------------- */
void restore_pixels(int px, std::vector<CoordPixel>& pts, std::vector<uint8_t>& image, std::map<int,int>& grey_lim, int width)
{
	if ( (pts[px+1].label == -1) && ((image[px+1] < image[px]) || image[px+1] < grey_lim[pts[px].label]) ) {
		pts[px+1].label = pts[px].label;
		restore_pixels( px+1, pts, image, grey_lim, width );
	}
	if ( (pts[px-1].label == -1) && ((image[px-1] < image[px]) || image[px-1] < grey_lim[pts[px].label]) ) {
		pts[px-1].label = pts[px].label;
		restore_pixels( px-1, pts, image, grey_lim, width );
	}
	if ( (pts[px+width].label == -1) && ((image[px+width] < image[px]) || image[px+width] < grey_lim[pts[px].label]) ) {
		pts[px+width].label = pts[px].label;
		restore_pixels( px+width, pts, image, grey_lim, width );
	}
	if ( (pts[px-width].label == -1) && ((image[px-width] < image[px]) || image[px-width] < grey_lim[pts[px].label]) ) {
		pts[px-width].label = pts[px].label;
		restore_pixels( px-width, pts, image, grey_lim, width );
	}
	if ( (pts[px+width+1].label == -1) && ((image[px+width+1] < image[px]) || image[px+width+1] < grey_lim[pts[px].label]) ) {
		pts[px+width+1].label = pts[px].label;
		restore_pixels( px+width+1, pts, image, grey_lim, width );
	}
	if ( (pts[px+width-1].label == -1) && ((image[px+width-1] < image[px]) || image[px+width-1] < grey_lim[pts[px].label]) ) {
		pts[px+width-1].label = pts[px].label;
		restore_pixels( px+width-1, pts, image, grey_lim, width );
	}
	if ( (pts[px-width+1].label == -1) && ((image[px-width+1] < image[px]) || image[px-width+1] < grey_lim[pts[px].label]) ) {
		pts[px-width+1].label = pts[px].label;
		restore_pixels( px-width+1, pts, image, grey_lim, width );
	}	
	if ( (pts[px-width-1].label == -1) && ((image[px-width-1] < image[px]) || image[px-width-1] < grey_lim[pts[px].label]) ) {
		pts[px-width-1].label = pts[px].label;
		restore_pixels( px-width-1, pts, image, grey_lim, width );
	}
}


int main() {

Timer time;

constexpr int perc{2}, maxvox{500};
//constexpr float trac_rad {4.35327}, scale {5.0};
constexpr float trac_rad {1.5}, scale {5.0};
int width, height, channels, labelcount;

std::vector<std::string> imagenames, outimagenames;
std::vector<CoordPixel>	pts0,pts1,pts2,pts3,pts4,pts5,pts6;
std::vector<uint8_t> img0, img1, img2, img3, img4, img5, img6;
std::vector<CoordFloat> angle_deg, angle_rad;
angle_deg = readCoord<CoordFloat>("datafiles/angles.dat");
angle_rad = DegToRad(angle_deg);

Eigen::RowVector3f rs, l;
rs << (Eigen::RowVector3f() << -690,0,0).finished();

CoordFloat lit;

for (int t=0; t<1; t++) {

	//std::vector<CoordFloat> centres = readCoord<CoordFloat>("datafiles/centres/StaticPacking/SB800_p"+std::to_string(perc)+"_t"+std::to_string(t)+".dat");
	std::vector<CoordFloat> centres = readCoord<CoordFloat>("centres.dat");

//	for(CoordFloat& cent : centres) {
//		cent.x = (cent.x-maxvox/2)/scale;
//		cent.y = (cent.y-maxvox/2)/scale;
//		cent.z = (cent.z-maxvox/2)/scale;
//	}

	imagenames = generate_fnames(angle_deg, t, perc);
	stbi_info(imagenames[0].c_str(), &width, &height, &channels);
	outimagenames = generate_outfnames(angle_deg, t, perc);
	std::vector<uint8_t> img0 = readImage(imagenames[0].c_str());

	pts0 = get_global(binarise(img0, "minimum"), angle_rad[0]);

	l = rs*rotmat(angle_rad[0].x,angle_rad[0].y,angle_rad[0].z,false);
	lit.x = l[0]; lit.y = l[1]; lit.z = l[2];
	for (CoordFloat& cent : centres) {
		for (CoordPixel& pt : pts0) {
			if (pt.label != -1) {
				if (lineptdistcheck(lit, pt, cent, trac_rad)) {
					pt.greyval = 0;
					pt.label = -1;
				}
				else pt.label = 0;
			}
		}
	}

//	labelcount = 0;
//	for (CoordPixel& pt : pts0) {
//		if (pt.label == 0) {
//			pt.label = ++labelcount;
//			label_connected_components(pt.pxnum, pts0, width);
//		}
//	}

//	std::map<int,int> label_cnt, grey_lim;
//	std::map<int, float> grey_avg, grey_stddev;

//	for (int i=-1; i<=labelcount; i++) {
//		label_cnt.insert({i,0});
//		grey_avg.insert({i,0.0});
//		grey_stddev.insert({i,0.0});
//	}

//	for (CoordPixel& pt : pts0) if (pt.label != -1) label_cnt[pt.label] += 1;
//	for (CoordPixel& pt : pts0) {
//		if (label_cnt[pt.label] < 1000) pt.label = -1;
//		else grey_avg[pt.label] += img0[pt.pxnum];

//	}

//	for (int i=0; i<=labelcount; i++) if (grey_avg[i] != 0) grey_avg[i] /= label_cnt[i];

//	for (CoordPixel& pt : pts0)	{
//		if (pt.label != -1) grey_stddev[pt.label] += std::pow((img0[pt.pxnum] - grey_avg[pt.label]),2);
//	}

//	for (int i=0; i<=labelcount; i++) {
//		if (label_cnt[i] != 0) grey_stddev[i] = std::pow((grey_stddev[i] / static_cast<float>(label_cnt[i])),0.5);
//	}

//	for (int i=0; i<=labelcount; i++) grey_lim[i] = static_cast<int>(grey_avg[i] - 0.0*grey_stddev[i]);

//	for (CoordPixel& pt : pts0) if (pt.label != -1) restore_pixels(pt.pxnum, pts0, img0, grey_lim, width);

	for (CoordPixel& pt : pts0) if (pt.label == -1) img0[pt.pxnum] = 0;
	stbi_write_png(outimagenames[0].c_str(), width, height, 1, img0.data(), width);
}

return 0;
}

