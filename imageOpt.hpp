#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <sstream>
#include <boost/histogram.hpp>

//// Read image : convert rgba to grey ----------------------------------------------
std::vector<uint8_t> readImage(const char* filename) {
	int width, height, num_channels;
	uint8_t* img_data = stbi_load(filename, &width, &height, &num_channels, 4);
	std::vector<uint8_t> greyimage(width * height);

	// Greyscale from RGB channels (ignore A) = 0.299*r + 0.587*g + 0.114*b
	for (int i = 0; i < width*height; ++i) {
		greyimage[i] = static_cast<uint8_t>(0.299*img_data[4*i]+0.587*img_data[4*i+1]+0.114*img_data[4*i+2]);
	}

	stbi_image_free(img_data); // free memory

	return greyimage;
}

//// Binirisation loop with global threshold -----------------------------------------------
std::vector<uint8_t> run_bin( std::vector<uint8_t>& greyimage, int& threshold ) {

	for (int i = 0; i < greyimage.size(); ++i) {
		if (greyimage[i] >= threshold) greyimage[i] = 0;
		else greyimage[i] = 255;
	}

	return greyimage;
}

//// Test if histogram is bimodal ---------------------------------------------------
bool bimodal(boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>& imghist) {
	bool b = false;
	int modes {0}, k;

	for (k=1; k<255; k++) {
		if ((imghist[k-1] < imghist[k]) && (imghist[k+1] < imghist[k])) {
			modes++;
			if (modes > 2) return false;
		}
	}
	if (modes==2) b = true;

	return b;
}

//// Minimum method -----------------------------------------------------------
int minimum ( boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>& hist ) {
	auto temp_hist = hist;
	int i {}, iter {}, max {}, threshold {};

	for (auto&& x : boost::histogram::indexed(hist)) {
		if (*x > 0) max = x.index();
	}

	while( !bimodal(hist) ) {

		for (i = 1; i < 255; ++i) {
			temp_hist[i] = (hist[i-1] + hist[i] + hist[i+1])/3;
		}
		temp_hist[0] = (hist[0] + hist[1])/3;
		temp_hist[255] = (hist[254] + hist[255])/3;

		hist = temp_hist;
		iter++;
	
		if (iter > 10000) {
			std::cout << "Threshold not found, try different method." << '\n';
			std::exit(0);
		}
	}

	for (i = 1; i < max; ++i) {
		if ((hist[i-1] > hist[i]) && (hist[i+1] >= hist[i])) {
			threshold = i;
			break;
		}
	}

	return threshold;
}

//// Default method from ImageJ ---------------------------------------------------
int ijdefault (boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>& hist ) {
	int i{}, level {}, min {}, max {255}, movingIndex, inc;
	double result {}, s1 {}, s2 {}, s3 {}, s4 {};

	while ( (hist[min]==0) && (min<255) ) min++;
	while ( (hist[max]==0) && (max>0) ) max--;

	if ( min >= max ) {
		level = 128;
		return level;
	}

	movingIndex = min;
	inc = std::max( max/40, 1 );

	do {
		s1 = s2 = s3 = s4 = 0.0;
		for (i = min; i<=movingIndex; i++) {
			s1 += i*hist[i];
			s2 += hist[i];
		}
		for (i = movingIndex+1; i<=max; i++) {
			s3 += i*hist[i];
			s4 += hist[i];
		}
		result = (s1/s2 + s3/s4)/2.0;
		movingIndex++;
	} while ( (movingIndex+1) <= result && movingIndex < (max-1) );

	level = static_cast<int>(result);
	return level;
}

//// Intermodes method -----------------------------------------------------------
int intermodes ( boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>& hist ) {

	int iter{}, threshold{-1}, i{}, tt{};
	double prev {}, curr {}, nxt {hist[0]};
	
	while( !bimodal(hist) ) {
		prev = 0; curr = 0; nxt = hist[0];
		for (i=0; i<255; i++) {
			prev = curr;
			curr = nxt;
			nxt = hist[i+1];
			hist[i] = (prev + curr + nxt)/3;
		}
		hist[255] = (curr + nxt)/3;
		iter++;
		if (iter > 10000) {
			threshold = -1;
			std::cout << "Threshold not found, try different method." << '\n';
			std::exit(0);
		}
	}

	for (i=1; i<255; i++) {
		if ( (hist[i-1] < hist[i]) && (hist[i+1] < hist[i]) ) tt += i;
	}

	threshold = static_cast<int> (tt/2.0);

	return threshold;
}

//// ISOData method ------------------------------------------------------------------------
int isodata ( boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>& hist ) {
	int i, l, toth, totl, h, g {};

	for (i=1; i<256; i++) {
		if (hist[i] > 0) {
			g = i+1;
			break;
		}
	}

	while (true) {
		l=0; totl=0;
		for (i=0; i<g; i++) {
			totl += hist[i];
			l += hist[i] * i;
		}

		h =0; toth = 0;
		for (i=g+1; i<256; i++) {
			toth += hist[i];
			h += hist[i]*i;
		}

		if (totl > 0 && toth >0) {
			l /= totl; h/= toth;
			if (g == static_cast<int>((l+h)/2.0)) break;
		}
		g++;
		if (g > 254) {
			std::cout << "Threshold not found, try different method." << '\n';
			std::exit(0);
		}
	}

	return g;
}

//// OTSU method ------------------------------------------------------------------------
int otsu ( boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double>>>& hist ) {

	int ih, threshold{-1}, num_pixels{0};
	double total_mean, bcv, term, max_bcv {};
	std::vector<double> normhist(256), meanhist(256), cumhist(256);

	// Calculate total number of pixels
	for (ih=0; ih<256; ih++) {
		num_pixels += hist[ih];
	}

	term = 1.0/static_cast<double>(num_pixels);

	// Calculate normalised histogram
	for (ih=0; ih<256; ih++) {
		normhist[ih] = term * hist[ih];
	}

	// Calculate cumulative normalised histogram
	cumhist[0] = normhist[0];
	for (ih=1; ih<256; ih++) {
		cumhist[ih] = cumhist[ih-1] + normhist[ih];
	}

	// Calculate mean grey-level
	meanhist[0] = 0.0;
	for (ih=1; ih<256; ih++) {
		meanhist[ih] = meanhist[ih-1] + ih*normhist[ih];
	}

	total_mean = meanhist[255];

	// Calculate background class variance at each greylevel and find threshold that maximizes it
	for (ih=0; ih<256; ih++) {
		bcv = total_mean * cumhist[ih] - meanhist[ih];
		bcv *= bcv / (cumhist[ih] * (1.0 - cumhist[ih]));

		if ( max_bcv < bcv ) {
			max_bcv = bcv;
			threshold = ih;
		}
	}

	return threshold;
}

//// Create histogram and call methods ---------------------------------------------------
std::vector<uint8_t> binarise( std::vector<uint8_t> greyimg, const std::string& bintype ) {

	int threshold {};

	auto hist = boost::histogram::make_histogram( boost::histogram::axis::regular<>(256, 0, 255) );
	for (uint8_t& value : greyimg) hist(value);
	//for (auto&& x : boost::histogram::indexed(hist)) x.index(), x.bin().lower(), x.bin().upper(), *x;

	if (bintype == "manual") {
		std::cout << "Please enter threshold for binirisation (0-255) : ";
		std::cin >> threshold;
	}

	else if (bintype == "ijdefault") threshold = ijdefault( hist );
	else if (bintype == "minimum") threshold = minimum( hist );
	else if (bintype == "intermodes") threshold = intermodes( hist );
	else if (bintype == "isodata") threshold = isodata( hist );
	else if (bintype == "otsu") threshold = otsu( hist );

	else {
		std::cout << "Please choose valid binarisation method." << '\n';
		std::exit(0);
	}

	greyimg = run_bin(greyimg, threshold);
	return greyimg;
}

std::vector<CoordPixel> get_global(std::vector<uint8_t> greyimage, CoordFloat& angle) {

	int width {1936}, height {1536}, i, j;
	float conv {75e-3};	// in mm
	CoordPixel px;
	std::vector<CoordPixel> imageGlobal;
	Eigen::Vector3f pt;
	Eigen::Matrix3f rot = rotmat(angle.x,angle.y,angle.z,true);

	for (i=0; i<height; ++i) {
		for (j=0; j<width; ++j) {
			if (greyimage[i*width+j] != 0) {
				pt << 110.0, (968.0-j)*conv, (768.0-i)*conv;
				pt = rot*pt;
				px.x = pt(0); px.y = pt(1); px.z = pt(2);
				px.greyval = greyimage[i*width+j];
				px.pxnum = i*width+j;
				px.label = 0;
				imageGlobal.emplace_back(px);
			}
			else {
				px.x = 0; px.y = 0; px.z = 0;
				px.greyval = 0;
				px.pxnum = i*width+j;
				px.label=-1;
				imageGlobal.emplace_back(px);
			}
		}
	}

	return imageGlobal;
}

std::vector<std::string> generate_fnames(std::vector<CoordFloat>& angles, int& t, const int& p) {

	std::vector<std::string> filenames;
	for (CoordFloat& angle : angles) {
		std::ostringstream theta, gamma, phi, time, perc;
		theta.precision(1);	gamma.precision(1); phi.precision(1), time.precision(0);
		theta << std::fixed << angle.x;
		gamma << std::fixed << angle.y;
		phi << std::fixed << angle.z;
		time << std::fixed << t;
		perc << std::fixed << p;
		//std::string imagename {"datafiles/projections/StaticPacking/p"+perc.str()+"_t"+time.str()+"_the"+theta.str()+"_phi"+gamma.str()+"_gam"+phi.str()+".png"};
		std::string imagename {"p"+perc.str()+"_t"+time.str()+"_the"+theta.str()+"_phi"+gamma.str()+"_gam"+phi.str()+".png"};
		filenames.emplace_back(imagename);
	}
	return filenames;
}

std::vector<std::string> generate_outfnames(std::vector<CoordFloat>& angles, int& t, const int& p) {

	std::vector<std::string> filenames;
	for (CoordFloat& angle : angles) {
		std::ostringstream theta, gamma, phi, time, perc;
		theta.precision(1);	gamma.precision(1); phi.precision(1), time.precision(0);
		theta << std::fixed << angle.x;
		gamma << std::fixed << angle.y;
		phi << std::fixed << angle.z;
		time << std::fixed << t;
		perc << std::fixed << p;
		std::string imagename {"greydel_p"+perc.str()+"_t"+time.str()+"_the"+theta.str()+"_phi"+gamma.str()+"_gam"+phi.str()+".png"};
		filenames.emplace_back(imagename);
	}
	return filenames;
}
