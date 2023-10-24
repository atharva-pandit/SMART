#include "basic_func.h"
#include "nanoflann.hpp"

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float,PointCloud>,PointCloud,3>;

int main() {

Timer time;
// Output file
int perc {5};
const std::string ptsfile {"datafiles/centres/SB5_traj"+std::to_string(perc)+".dat"};
clearfile(ptsfile);

// Setup variables
float threshSq {5}, best_dist, cost;	// Max dist moved by a particle in 1 timestep
size_t i,j,k,l,t,id10,id21,id32,nmatches10,nmatches21,nmatches32, tmin {0}, tmax {19}, new_id {1};
PointCloud data0, data1, data2, data3, memorylist;
std::vector<PointCloud> fwd_data, bck_data;
CoordId surePt; std::vector<CoordId> comm_data;

// Forward direction starting from tmin -----------------------------------------
data0 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(tmin)+ "_cent.dat");
data1 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(tmin+1)+ "_cent.dat");
data2 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(tmin+2)+ "_cent.dat");
data3 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(tmin+3)+ "_cent.dat");

KDTree index0(3,data0,{10});	KDTree index1(3,data1,{10});
KDTree index2(3,data2,{10});	KDTree index3(3,data3,{10});

std::vector<nanoflann::ResultItem<uint32_t, float>> resultSet10, resultSet21, resultSet32;

// Initialise using 4FBE - forward
for (i=0; i<data0.kdtree_get_point_count(); i++) {
	data0.pts[i].id = new_id++;
	float queryPt10[3];
	queryPt10[0] = data0.pts[i].x;
	queryPt10[1] = data0.pts[i].y;
	queryPt10[2] = data0.pts[i].z;
	nmatches10 = index1.radiusSearch(queryPt10, threshSq , resultSet10);
	best_dist = threshSq, cost = 0;		// reset cost_func

	// go through all neighbours
	for (j=0; j<nmatches10; j++) {
		id10 = resultSet10[j].first;
		float neib10[3] {data1.pts[id10].x, data1.pts[id10].y, data1.pts[id10].z};	// nn in t=1
		float vel10[3] {neib10[0]-queryPt10[0],neib10[1]-queryPt10[1],neib10[2]-queryPt10[2]};	// velocity = x1-x0
		float queryPt21[3] {neib10[0]+vel10[0],neib10[1]+vel10[1],neib10[2]+vel10[2]};	// prediction for t=2 ; x+vt
		nmatches21 = index2.radiusSearch(queryPt21, threshSq , resultSet21);	// search nn of predicted point in t=2

		for (k=0; k<nmatches21; k++) {
			id21 = resultSet21[k].first;
			float neib21[3] {data2.pts[id21].x, data2.pts[id21].y, data2.pts[id21].z};	// nn of predicted point in t=2
			float vel21[3] {neib21[0]-neib10[0],neib21[1]-neib10[1],neib21[2]-neib10[2]};	// velocity = x2-x1
			float acc20[3] {vel21[0]-vel10[0],vel21[1]-vel10[1],vel21[2]-vel10[2]}; // acceleration = v2-v1
			// prediction for t=3; qPt = x + v1*(2t) + 0.5*a*(2t)^2
			float queryPt32[3] {queryPt10[0] + 2*vel10[0] + 2*acc20[0],
								 queryPt10[1] + 2*vel10[1] + 2*acc20[1],
								 queryPt10[2] + 2*vel10[2] + 2*acc20[2]};
			nmatches32 = index3.radiusSearch(queryPt32, threshSq , resultSet32);	// search nn of predicted point in t=3

			for (l=0; l<nmatches32; l++) {
				id32 = resultSet32[l].first; cost = resultSet32[l].second;	// nn of predicted point in t=3; id & dist
				if ( cost < best_dist ) {				// smallest dist is the best path
					best_dist = cost;
					data1.pts[id10].id = data0.pts[i].id;
					data2.pts[id21].id = data0.pts[i].id;
					data3.pts[id32].id = data0.pts[i].id;

					float vel32[3] {data3.pts[id32].x-neib21[0],data3.pts[id32].y-neib21[1],data3.pts[id32].z-neib21[2]};
					data1.pts[id10].vx = vel10[0]; data1.pts[id10].vy = vel10[1]; data1.pts[id10].vz = vel10[2];
					data2.pts[id21].vx = vel21[0]; data2.pts[id21].vy = vel21[1]; data2.pts[id21].vz = vel21[2];
					data3.pts[id32].vx = vel32[0]; data3.pts[id32].vy = vel32[1]; data3.pts[id32].vz = vel32[2];

					data2.pts[id21].ax = acc20[0]; data2.pts[id21].ay = acc20[1]; data2.pts[id21].az = acc20[2];
					float acc31[3] {vel32[0]-vel21[0],vel32[1]-vel21[1],vel32[2]-vel21[2]};
					data3.pts[id32].ax = acc31[0]; data3.pts[id32].ay = acc31[1]; data3.pts[id32].az = acc31[2];
				}
			}
		}
	}
}

fwd_data.emplace_back(data0);
fwd_data.emplace_back(data1);
fwd_data.emplace_back(data2);
fwd_data.emplace_back(data3);

for (t=(tmin+3); t<tmax; t++) {
	if (t==(tmin+3)) data0 = data3;
	data1 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(t+1)+ "_cent.dat");

	KDTree index0(3,data0,{10}); KDTree index1(3,data1,{10});

	// First check if "remembered points" from t-1 have a match in t+1
	for (i=0; i<memorylist.kdtree_get_point_count(); i++) {
		float queryPt10[3];
		queryPt10[0] = data0.pts[i].x + 2*data0.pts[i].vx + 2*data0.pts[i].ax;
		queryPt10[1] = data0.pts[i].y + 2*data0.pts[i].vy + 2*data0.pts[i].ay;
		queryPt10[2] = data0.pts[i].z + 2*data0.pts[i].vz + 2*data0.pts[i].az;
		nmatches10 = index1.radiusSearch(queryPt10, threshSq , resultSet10);	// search nn in t=1
		if (nmatches10 > 0) {
			id10 = resultSet10[0].first;	// select closest point
			data1.pts[id10].id = memorylist.pts[i].id;
			data1.pts[id10].vx = (data1.pts[id10].x - memorylist.pts[i].x) / 2.0;	//dt = 2
			data1.pts[id10].vy = (data1.pts[id10].y - memorylist.pts[i].y) / 2.0;
			data1.pts[id10].vz = (data1.pts[id10].z - memorylist.pts[i].z) / 2.0;
			data1.pts[id10].ax = (data1.pts[id10].vx - memorylist.pts[i].vx) / 2.0;
			data1.pts[id10].ay = (data1.pts[id10].vy - memorylist.pts[i].vy) / 2.0;
			data1.pts[id10].az = (data1.pts[id10].vz - memorylist.pts[i].vz) / 2.0;
		}
	}
	memorylist.pts.clear();

	for (i=0; i<data0.kdtree_get_point_count(); i++) {
		if (data0.pts[i].id != 0) {
			float queryPt10[3];
			queryPt10[0] = data0.pts[i].x + data0.pts[i].vx + 0.5*data0.pts[i].ax;
			queryPt10[1] = data0.pts[i].y + data0.pts[i].vy + 0.5*data0.pts[i].ay;
			queryPt10[2] = data0.pts[i].z + data0.pts[i].vz + 0.5*data0.pts[i].az;
			nmatches10 = index1.radiusSearch(queryPt10, threshSq , resultSet10);	// search nn in t=1

			if (nmatches10 > 0) {
				id10 = resultSet10[0].first;	// select closest point
				data1.pts[id10].id = data0.pts[i].id;
				data1.pts[id10].vx = data1.pts[id10].x - data0.pts[i].x;
				data1.pts[id10].vy = data1.pts[id10].y - data0.pts[i].y;
				data1.pts[id10].vz = data1.pts[id10].z - data0.pts[i].z;
				data1.pts[id10].ax = data1.pts[id10].vx - data0.pts[i].vx;
				data1.pts[id10].ay = data1.pts[id10].vy - data0.pts[i].vy;
				data1.pts[id10].az = data1.pts[id10].vz - data0.pts[i].vz;
			}
			
			else {	// "remember" the particle for 1 more frame
				memorylist.pts.emplace_back(data0.pts[i]);
			}
		}
	}

	fwd_data.emplace_back(data1);
	data0 = data1;
}

// REVERSE DIRECTION starting from tmax ---------------------------------------------------------------------

data0 = fwd_data[tmax];
data1 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(tmax-1)+ "_cent.dat");
data2 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(tmax-2)+ "_cent.dat");
data3 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(tmax-3)+ "_cent.dat");

KDTree bindex0(3,data0,{10});	KDTree bindex1(3,data1,{10});
KDTree bindex2(3,data2,{10});	KDTree bindex3(3,data3,{10});

std::vector<nanoflann::ResultItem<uint32_t, float>> resultSet01, resultSet12, resultSet23;

// Initialise using 4FBE - backward
for (i=0; i<data0.kdtree_get_point_count(); i++) {
	if (data0.pts[i].id != 0) {	// only established trajectories
		float queryPt10[3];
		queryPt10[0] = data0.pts[i].x;
		queryPt10[1] = data0.pts[i].y;
		queryPt10[2] = data0.pts[i].z;
		nmatches10 = bindex1.radiusSearch(queryPt10, threshSq , resultSet01);
		best_dist = threshSq, cost = 0;		// reset cost_func

		// go through all neighbours
		for (j=0; j<nmatches10; j++) {
			id10 = resultSet01[j].first;
			float neib10[3] {data1.pts[id10].x, data1.pts[id10].y, data1.pts[id10].z};	// nn in t=1
			float vel10[3] {neib10[0]-queryPt10[0],neib10[1]-queryPt10[1],neib10[2]-queryPt10[2]};	// velocity = x1-x0
			float queryPt21[3] {neib10[0]+vel10[0],neib10[1]+vel10[1],neib10[2]+vel10[2]};	// prediction for t=2 ; x+vt
			nmatches21 = bindex2.radiusSearch(queryPt21, threshSq , resultSet12);	// search nn of predicted point in t=2

			for (k=0; k<nmatches21; k++) {
				id21 = resultSet12[k].first;
				float neib21[3] {data2.pts[id21].x, data2.pts[id21].y, data2.pts[id21].z};	// nn of predicted point in t=2
				float vel21[3] {neib21[0]-neib10[0],neib21[1]-neib10[1],neib21[2]-neib10[2]};	// velocity = x2-x1
				float acc20[3] {vel21[0]-vel10[0],vel21[1]-vel10[1],vel21[2]-vel10[2]}; // acceleration = v2-v1
				// prediction for t=3; qPt = x + v1*(2t) + 0.5*a*(2t)^2
				float queryPt32[3] {queryPt10[0] + 2*vel10[0] + 2*acc20[0],
									 queryPt10[1] + 2*vel10[1] + 2*acc20[1],
									 queryPt10[2] + 2*vel10[2] + 2*acc20[2]};
				nmatches32 = bindex3.radiusSearch(queryPt32, threshSq , resultSet23);	// search nn of predicted point in t=3

				for (l=0; l<nmatches32; l++) {
					id32 = resultSet23[l].first; cost = resultSet23[l].second;	// nn of predicted point in t=3; id & dist
					if ( cost < best_dist ) {				// smallest dist is the best path
						best_dist = cost;
						data1.pts[id10].id = data0.pts[i].id;
						data2.pts[id21].id = data0.pts[i].id;
						data3.pts[id32].id = data0.pts[i].id;

						float vel32[3] {data3.pts[id32].x-neib21[0],data3.pts[id32].y-neib21[1],data3.pts[id32].z-neib21[2]};
						data1.pts[id10].vx = vel10[0]; data1.pts[id10].vy = vel10[1]; data1.pts[id10].vz = vel10[2];
						data2.pts[id21].vx = vel21[0]; data2.pts[id21].vy = vel21[1]; data2.pts[id21].vz = vel21[2];
						data3.pts[id32].vx = vel32[0]; data3.pts[id32].vy = vel32[1]; data3.pts[id32].vz = vel32[2];

						data2.pts[id21].ax = acc20[0]; data2.pts[id21].ay = acc20[1]; data2.pts[id21].az = acc20[2];
						float acc31[3] {vel32[0]-vel21[0],vel32[1]-vel21[1],vel32[2]-vel21[2]};
						data3.pts[id32].ax = acc31[0]; data3.pts[id32].ay = acc31[1]; data3.pts[id32].az = acc31[2];
					}
				}
			}
		}
	}
}

bck_data.insert(bck_data.begin(), data0);
bck_data.insert(bck_data.begin(), data1);
bck_data.insert(bck_data.begin(), data2);
bck_data.insert(bck_data.begin(), data3);

for (t=(tmax-3); t>tmin; t--) {
	if (t==(tmax-3)) data0 = data3;
	data1 = readPtCloud("datafiles/centres/SB5_p"+std::to_string(perc)+"_t" +std::to_string(t-1)+ "_cent.dat");

	KDTree bindex0(3,data0,{10});	KDTree bindex1(3,data1,{10});

	// First check if "remembered points" from t-1 have a match in t+1
	for (i=0; i<memorylist.kdtree_get_point_count(); i++) {
		float queryPt10[3];
		queryPt10[0] = data0.pts[i].x + 2*data0.pts[i].vx + 2*data0.pts[i].ax;
		queryPt10[1] = data0.pts[i].y + 2*data0.pts[i].vy + 2*data0.pts[i].ay;
		queryPt10[2] = data0.pts[i].z + 2*data0.pts[i].vz + 2*data0.pts[i].az;
		nmatches10 = bindex1.radiusSearch(queryPt10, threshSq , resultSet01);	// search nn in t=1
		if (nmatches10 > 0) {
			id10 = resultSet01[0].first;	// select closest point
			data1.pts[id10].id = memorylist.pts[i].id;
			data1.pts[id10].vx = (data1.pts[id10].x - memorylist.pts[i].x) / 2.0;	//dt = 2
			data1.pts[id10].vy = (data1.pts[id10].y - memorylist.pts[i].y) / 2.0;
			data1.pts[id10].vz = (data1.pts[id10].z - memorylist.pts[i].z) / 2.0;
			data1.pts[id10].ax = (data1.pts[id10].vx - memorylist.pts[i].vx) / 2.0;
			data1.pts[id10].ay = (data1.pts[id10].vy - memorylist.pts[i].vy) / 2.0;
			data1.pts[id10].az = (data1.pts[id10].vz - memorylist.pts[i].vz) / 2.0;
		}
	}
	memorylist.pts.clear();

	for (i=0; i<data0.kdtree_get_point_count(); i++) {
		if (data0.pts[i].id != 0) {
			float queryPt10[3];
			queryPt10[0] = data0.pts[i].x + data0.pts[i].vx + 0.5*data0.pts[i].ax;
			queryPt10[1] = data0.pts[i].y + data0.pts[i].vy + 0.5*data0.pts[i].ay;
			queryPt10[2] = data0.pts[i].z + data0.pts[i].vz + 0.5*data0.pts[i].az;
			nmatches10 = bindex1.radiusSearch(queryPt10, threshSq , resultSet01);	// search nn in t=1

			if (nmatches10 > 0) {
				id10 = resultSet01[0].first;	// select closest point
				data1.pts[id10].id = data0.pts[i].id;
				data1.pts[id10].vx = data1.pts[id10].x - data0.pts[i].x;
				data1.pts[id10].vy = data1.pts[id10].y - data0.pts[i].y;
				data1.pts[id10].vz = data1.pts[id10].z - data0.pts[i].z;
				data1.pts[id10].ax = data1.pts[id10].vx - data0.pts[i].vx;
				data1.pts[id10].ay = data1.pts[id10].vy - data0.pts[i].vy;
				data1.pts[id10].az = data1.pts[id10].vz - data0.pts[i].vz;
			}

			else {	// "remember" the particle for 1 more frame
				memorylist.pts.emplace_back(data0.pts[i]);
			}
		}
	}

	bck_data.insert(bck_data.begin(), data1);
	data0 = data1;
}

// Go through both data sets to find common trajectories
for (t=0; t<tmax; t++) {
	data0 = fwd_data[t]; data1 = bck_data[t];
	for (i=0; i<data0.kdtree_get_point_count(); i++) {
		if (data0.pts[i].id == data1.pts[i].id) {
			if ((data0.pts[i].x == data1.pts[i].x) && (data0.pts[i].y == data1.pts[i].y) && (data0.pts[i].z == data1.pts[i].z)) {
				surePt.x = data0.pts[i].x;
				surePt.y = data0.pts[i].y;
				surePt.z = data0.pts[i].z;
				surePt.t = t;
				surePt.id = data0.pts[i].id;
				comm_data.emplace_back(surePt);
			}
		}
	}
}

writeCoordId(ptsfile, comm_data);

return 0;
}
