/* Get max value from reconstructed grid ------- */
float voxgrid_max ()
{
	float max_greyvalue {};

	for (int i=0; i<maxvox; ++i) {
		for (int j=0; j<maxvox; ++j) {
			for (int k=0; k<maxvox; ++k) {
				if ( voxgrid[i][j][k] > max_greyvalue ) {
					max_greyvalue = voxgrid[i][j][k];
				}
			}
		}
	}
	return max_greyvalue;
}

/* Label background and foreground voxels to -1 & 0 ------ */
void initialise_labelgrid (float maxgrey)
{
	for (int i=0; i<maxvox; ++i) {
		for (int j=0; j<maxvox; ++j) {
			for (int k=0; k<maxvox; ++k) {
				if ( voxgrid[i][j][k] < maxgrey*0.5 ) labels[i][j][k] = -1;
				else labels[i][j][k] = 0;
			}
		}
	}
}

/* Label connected components --------------- */
void update_labels( int i, int j, int k, int labelcount, std::map<int,int>& lab_cnt )
{
	if (labels[i+1][j][k] == 0) {
		labels[i+1][j][k] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j,k, labelcount, lab_cnt);
	}

	if (labels[i][j+1][k] == 0) {
		labels[i][j+1][k] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i,j+1,k, labelcount, lab_cnt);
	}

	if (labels[i][j][k+1] == 0) {
		labels[i][j][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i,j,k+1, labelcount, lab_cnt);
	}

	if (labels[i-1][j][k] == 0) {
		labels[i-1][j][k] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j,k, labelcount, lab_cnt);
	}

	if (labels[i][j-1][k] == 0) {
		labels[i][j-1][k] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i,j-1,k, labelcount, lab_cnt);
	}

	if (labels[i][j][k-1] == 0) {
		labels[i][j][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i,j,k-1, labelcount, lab_cnt);
	}

	if (labels[i+1][j+1][k] == 0) {
		labels[i+1][j+1][k] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j+1,k, labelcount, lab_cnt);
	}

	if (labels[i-1][j-1][k] == 0) {
		labels[i-1][j-1][k] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j-1,k, labelcount, lab_cnt);
	}

	if (labels[i+1][j][k+1] == 0) {
		labels[i+1][j][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j,k+1, labelcount, lab_cnt);
	}

	if (labels[i-1][j][k-1] == 0) {
		labels[i-1][j][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j,k-1, labelcount, lab_cnt);
	}

	if (labels[i][j+1][k+1] == 0) {
		labels[i][j+1][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i,j+1,k+1, labelcount, lab_cnt);
	}

	if (labels[i][j-1][k-1] == 0) {
		labels[i][j-1][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i,j-1,k-1, labelcount, lab_cnt);
	}

	if (labels[i+1][j-1][k] == 0) {
		labels[i+1][j-1][k] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j-1,k, labelcount, lab_cnt);
	}

	if (labels[i-1][j+1][k] == 0) {
		labels[i-1][j+1][k] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j+1,k, labelcount, lab_cnt);
	}

	if (labels[i+1][j][k-1] == 0) {
		labels[i+1][j][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j,k-1, labelcount, lab_cnt);
	}

	if (labels[i-1][j][k+1] == 0) {
		labels[i-1][j][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j,k+1, labelcount, lab_cnt);
	}

	if (labels[i][j+1][k-1] == 0) {
		labels[i][j+1][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i,j+1,k-1, labelcount, lab_cnt);
	}

	if (labels[i][j-1][k+1] == 0) {
		labels[i][j-1][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i,j-1,k+1, labelcount, lab_cnt);
	}

	if (labels[i+1][j+1][k+1] == 0) {
		labels[i+1][j+1][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j+1,k+1, labelcount, lab_cnt);
	}

	if (labels[i-1][j-1][k-1] == 0) {
		labels[i-1][j-1][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j-1,k-1, labelcount, lab_cnt);
	}

	if (labels[i+1][j+1][k-1] == 0) {
		labels[i+1][j+1][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j+1,k-1, labelcount, lab_cnt);
	}

	if (labels[i+1][j-1][k+1] == 0) {
		labels[i+1][j-1][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j-1,k+1, labelcount, lab_cnt);
	}

	if (labels[i-1][j+1][k+1] == 0) {
		labels[i-1][j+1][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j+1,k+1, labelcount, lab_cnt);
	}

	if (labels[i+1][j-1][k-1] == 0) {
		labels[i+1][j-1][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i+1,j-1,k-1, labelcount, lab_cnt);
	}

	if (labels[i-1][j+1][k-1] == 0) {
		labels[i-1][j+1][k-1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j+1,k-1, labelcount, lab_cnt);
	}

	if (labels[i-1][j-1][k+1] == 0) {
		labels[i-1][j-1][k+1] = labelcount;
		lab_cnt[labelcount] += 1;
		update_labels(i-1,j-1,k+1, labelcount, lab_cnt);
	}
}

/* Add labels to all voxels, connected components have the same label -------- */
std::map<int,int> add_labels (int& labelcount)
{
	std::map<int, int> lab_cnt;

	for (int i=0; i<maxvox; i++) {
		for (int j=0; j<maxvox; j++) {
			for (int k=0; k<maxvox; k++) {
				if (labels[i][j][k] == 0) {
					labels[i][j][k] = ++labelcount;
					lab_cnt.insert({labelcount,1});
					update_labels(i,j,k, labelcount, lab_cnt);
				}
			}
		}
	}

	return lab_cnt;
}

/* Geometric mean of each cluster = particle centre ----------------- */
std::vector<CoordFloat> find_cluster_centres(std::map<int,int>& label_list)
{
	std::vector<CoordFloat> centres(label_list.size(),{0.0,0.0,0.0});
	std::map<int,int>::iterator it;
	std::map<int,int>::iterator it_begin = label_list.begin();
	std::map<int,int>::iterator it_end = label_list.end();
	int ind {}, i {};

	for (int i=0; i<maxvox; i++) {
		for (int j=0; j<maxvox; j++) {
			for (int k=0; k<maxvox; k++) {
				it = label_list.find(labels[i][j][k]);
				if ( it != it_end ) {
					ind = std::distance(it_begin, it);
					centres[ind].x += i;
					centres[ind].y += j;
					centres[ind].z += k;
				}
			}
		}
	}

	for (it=label_list.begin(); it != label_list.end(); it++) {
		centres[i].x /= it->second;
		centres[i].y /= it->second;
		centres[i].z /= it->second;
		i++;
	}

	return centres;
}

std::vector<CoordId> distanceMap(std::vector<CoordInt>& data)
{
	const int mxsq {maxvox*maxvox}, mx {maxvox};
	int ptsmarked {}, dist_now {1};
	std::map<int,CoordId>::iterator pt;
	CoordId newpoint;
	std::map<int,CoordId> points;
	std::vector<CoordId> finalpts;

	for (CoordInt& dd : data) {
		newpoint.x = dd.x;
		newpoint.y = dd.y;
		newpoint.z = dd.z;
		newpoint.t = 0;
		newpoint.id = 0;
		points.insert({dd.x*mxsq + dd.y*mx + dd.z, newpoint});
	}

	for(pt = points.begin(); pt != points.end(); pt++) {
		if(pt->second.id == 0){
			if ( (labels[pt->second.x+1][pt->second.y][pt->second.z] == -1) ||
			(labels[pt->second.x][pt->second.y+1][pt->second.z] == -1) ||
			(labels[pt->second.x][pt->second.y][pt->second.z+1] == -1) ||
			(labels[pt->second.x-1][pt->second.y][pt->second.z] == -1) ||
			(labels[pt->second.x][pt->second.y-1][pt->second.z] == -1) ||
			(labels[pt->second.x][pt->second.y][pt->second.z-1] == -1) ||
			(labels[pt->second.x+1][pt->second.y+1][pt->second.z] == -1) ||
			(labels[pt->second.x-1][pt->second.y-1][pt->second.z] == -1) ||
			(labels[pt->second.x+1][pt->second.y][pt->second.z+1] == -1) ||
			(labels[pt->second.x-1][pt->second.y][pt->second.z-1] == -1) ||
			(labels[pt->second.x][pt->second.y+1][pt->second.z+1] == -1) ||
			(labels[pt->second.x][pt->second.y-1][pt->second.z-1] == -1) ||
			(labels[pt->second.x+1][pt->second.y-1][pt->second.z] == -1) ||
			(labels[pt->second.x-1][pt->second.y+1][pt->second.z] == -1) ||
			(labels[pt->second.x+1][pt->second.y][pt->second.z-1] == -1) ||
			(labels[pt->second.x-1][pt->second.y][pt->second.z+1] == -1) ||
			(labels[pt->second.x][pt->second.y+1][pt->second.z-1] == -1) ||
			(labels[pt->second.x][pt->second.y-1][pt->second.z+1] == -1) ||
			(labels[pt->second.x+1][pt->second.y+1][pt->second.z+1] == -1) ||
			(labels[pt->second.x-1][pt->second.y-1][pt->second.z-1] == -1) ||
			(labels[pt->second.x+1][pt->second.y+1][pt->second.z-1] == -1) ||
			(labels[pt->second.x+1][pt->second.y-1][pt->second.z+1] == -1) ||
			(labels[pt->second.x-1][pt->second.y+1][pt->second.z+1] == -1) ||
			(labels[pt->second.x+1][pt->second.y-1][pt->second.z-1] == -1) ||
			(labels[pt->second.x-1][pt->second.y+1][pt->second.z-1] == -1) ||
			(labels[pt->second.x-1][pt->second.y-1][pt->second.z+1] == -1) )
				{ pt->second.id = dist_now; ptsmarked++; }
		}
	}

	while(ptsmarked != points.size()) {
		for(pt = points.begin(); pt != points.end(); pt++) {
			if (pt->second.id == 0) {
				if ( (points[(pt->second.x+1)*mxsq + pt->second.y*mx + pt->second.z].id == dist_now) ||
					(points[pt->second.x*mxsq + (pt->second.y+1)*mx + pt->second.z].id == dist_now) ||
					(points[pt->second.x*mxsq + pt->second.y*mx + pt->second.z+1].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + pt->second.y*mx + pt->second.z].id == dist_now) ||
					(points[pt->second.x*mxsq + (pt->second.y-1)*mx + pt->second.z].id == dist_now) ||
					(points[pt->second.x*mxsq + pt->second.y*mx + pt->second.z-1].id == dist_now) ||
					(points[(pt->second.x+1)*mxsq + (pt->second.y+1)*mx + pt->second.z].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + (pt->second.y-1)*mx + pt->second.z].id == dist_now) ||
					(points[(pt->second.x+1)*mxsq + pt->second.y*mx + pt->second.z+1].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + pt->second.y*mx + pt->second.z-1].id == dist_now) ||
					(points[pt->second.x*mxsq + (pt->second.y+1)*mx + pt->second.z+1].id == dist_now) ||
					(points[pt->second.x*mxsq + (pt->second.y+1)*mx + pt->second.z-1].id == dist_now) ||
					(points[(pt->second.x+1)*mxsq + (pt->second.y-1)*mx + pt->second.z].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + (pt->second.y+1)*mx + pt->second.z].id == dist_now) ||
					(points[(pt->second.x+1)*mxsq + pt->second.y*mx + pt->second.z-1].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + pt->second.y*mx + pt->second.z+1].id == dist_now) ||
					(points[pt->second.x*mxsq + (pt->second.y+1)*mx + pt->second.z-1].id == dist_now) ||
					(points[pt->second.x*mxsq + (pt->second.y-1)*mx + pt->second.z+1].id == dist_now) ||
					(points[(pt->second.x+1)*mxsq + (pt->second.y+1)*mx + pt->second.z+1].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + (pt->second.y-1)*mx + pt->second.z-1].id == dist_now) ||
					(points[(pt->second.x+1)*mxsq + (pt->second.y+1)*mx + pt->second.z-1].id == dist_now) ||
					(points[(pt->second.x+1)*mxsq + (pt->second.y-1)*mx + pt->second.z+1].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + (pt->second.y+1)*mx + pt->second.z+1].id == dist_now) ||
					(points[(pt->second.x+1)*mxsq + (pt->second.y-1)*mx + pt->second.z-1].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + (pt->second.y+1)*mx + pt->second.z-1].id == dist_now) ||
					(points[(pt->second.x-1)*mxsq + (pt->second.y-1)*mx + pt->second.z+1].id == dist_now) )
						{ pt->second.id = dist_now+1; ptsmarked++; }
			}
		}
		dist_now++;
	}

	for(pt = points.begin(); pt != points.end(); pt++) {
		if (pt->second.id >= dist_now-1) {
			finalpts.emplace_back(pt->second);
			labels[pt->second.x][pt->second.y][pt->second.z] = 0;
		}
		else labels[pt->second.x][pt->second.y][pt->second.z] = -1;
	}

	return finalpts;
}
