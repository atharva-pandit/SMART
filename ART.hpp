std::vector<CoordInt> bresenham(CoordFloat& strt, CoordFloat& end) {

	CoordInt stvox {static_cast<int>(strt.x * scale + maxvox/2),
					static_cast<int>(strt.y * scale + maxvox/2),
					static_cast<int>(strt.z * scale + maxvox/2)};
	CoordInt endvox {static_cast<int>(end.x * scale + maxvox/2),
					static_cast<int>(end.y * scale + maxvox/2),
					static_cast<int>(end.z * scale + maxvox/2)};

	std::vector<CoordInt> voxlist;
	int i, p1, p2;

	const int dx = endvox.x - stvox.x;
	const int dy = endvox.y - stvox.y;
	const int dz = endvox.z - stvox.z;
	const int x_inc = (dx < 0) ? -1 : 1;
	const int y_inc = (dy < 0) ? -1 : 1;
	const int z_inc = (dz < 0) ? -1 : 1;
	const int l = abs(dx); const int m = abs(dy); const int n = abs(dz);
	const int dx2 = l << 1; const int dy2 = m << 1; const int dz2 = n << 1;

	if ((l >= m) && (l >= n)) {
		p1 = dy2 - l;
		p2 = dz2 - l;
		for (i = 0; i < l; i++) {
			stvox.x += x_inc;
			if (p1 > 0) {
				stvox.y += y_inc;
				p1 -= dx2;
			}
			if (p2 > 0) {
				stvox.z += z_inc;
				p2 -= dx2;
			}
			p1 += dy2;
			p2 += dz2;
			voxlist.emplace_back(stvox);
		}
	}
	else if ((m >= l) && (m >= n)) {
		p1 = dx2 - m;
		p2 = dz2 - m;
		for (i = 0; i < m; i++){
			stvox.y += y_inc;
			if (p1 > 0){
				stvox.x += x_inc;
				p1 -= dy2;
			}
			if (p2 > 0){
				stvox.z += z_inc;
				p2 -= dy2;
			}
			p1 += dx2;
			p2 += dz2;
			voxlist.emplace_back(stvox);
		}
	}
	else {
		p1 = dy2 - n;
		p2 = dx2 - n;
		for (i = 0; i < n; i++){
			stvox.z += z_inc;
			if (p1 > 0){
				stvox.y += y_inc;
				p1 -= dz2;
			}
			if (p2 > 0){
				stvox.x += x_inc;
				p2 -= dz2;
			}
			p1 += dy2;
			p2 += dx2;
			voxlist.emplace_back(stvox);
		}
	}
	return voxlist;
}

void ART(std::vector<CoordInt>& raypts, int& greyval) {

	float sum_w {}, sum_E {}, mu {1.0};

	// summation of E along ray, because w_ij=1
	for (CoordInt& ray : raypts) {
		if ((ray.x < 0) || (ray.y < 0) || (ray.z < 0) || (ray.x >= maxvox) || (ray.y >= maxvox) || (ray.z >= maxvox)) {
			std::cout << "Out of bounds !" << '\n';
			std::cout << ray.x << " " << ray.y << " " << ray.z << '\n';
			std::exit(EXIT_FAILURE);
		}
		sum_E += voxgrid[ray.x][ray.y][ray.z];
		sum_w += 1;
	}

	// ART formula
	for (CoordInt& ray : raypts) {
		voxgrid[ray.x][ray.y][ray.z] += mu*(greyval-sum_E)/sum_w;
	}
}

void MART(std::vector<CoordInt>& raypts, int& greyval) {

	float sum_w {}, sum_E {}, mu {1.0};

	// summation of E along ray, because w_ij=1
	for (CoordInt& ray : raypts) {
		if ((ray.x < 0) || (ray.y < 0) || (ray.z < 0) || (ray.x >= maxvox) || (ray.y >= maxvox) || (ray.z >= maxvox)) {
			std::cout << "Out of bounds !" << '\n';
			std::cout << ray.x << " " << ray.y << " " << ray.z << '\n';
			std::exit(EXIT_FAILURE);
		}
		sum_E += voxgrid[ray.x][ray.y][ray.z];
		sum_w += 1;
	}

	// MART formula
	for (CoordInt& ray : raypts) {
		voxgrid[ray.x][ray.y][ray.z] *= std::pow((greyval/sum_E),mu);
	}
}

void simpback(std::vector<CoordInt>& raypts) {
	for (CoordInt& ray : raypts) voxgrid[ray.x][ray.y][ray.z] += 1;
}

void clear_voxgrid() {
	for (int i=0; i<maxvox; ++i) {
		for (int j=0; j<maxvox; ++j) {
			for (int k=0; k<maxvox; ++k) {
				voxgrid[i][j][k] = 0.0;
			}
		}
	}
}

void setupVoxgrid_MART() {
	for (int i=0; i<maxvox; i++) {
		for (int j=0; j<maxvox; j++) {
			for (int k=0; k<maxvox; k++) {
				voxgrid[i][j][k] = 1.0;
			}
		}
	}
}

void delete_voxgrid() { delete [] voxgrid; }
