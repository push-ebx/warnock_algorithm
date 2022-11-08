#define _USE_MATH_DEFINES
#include "plane.h"

#define f 50.0
#define n 1.0

#define a -(f + n) / (f - n)
#define b -2*f*n / (f - n)

vector<vector<int>> sides = {{0, 1, 2}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {0, 2, 5, 3}};
int colors[] = {CYAN, RED, YELLOW, BLUE, BROWN, GREEN, RED, BLUE, GREEN};
int color = 1;

struct centroid { float x, y, z; };

centroid get_centroid(const vector<array<float, 3>> &coords) {
	float x = 0, y = 0, z = 0;
	for (size_t i = 0; i < vertex; i++) x += coords[i][0], y += coords[i][1], z += coords[i][2];
	return {x /= vertex, y /= vertex, z /= vertex};
}

void rotate(vector<array<float, 3>> &coords, float ang[3]) {
	centroid cent = get_centroid(coords);
	Matrix<float> rotate_x(vector<vector<float>>({{1, 0, 0}, {0, cos(ang[0]), sin(ang[0])}, {0, -sin(ang[0]), cos(ang[0])}}));
	Matrix<float> rotate_y(vector<vector<float>>({{cos(ang[1]), 0, -sin(ang[1])}, {0, 1, 0}, {sin(ang[1]), 0, cos(ang[1])}}));
	Matrix<float> rotate_z(vector<vector<float>>({{cos(ang[2]), sin(ang[2]), 0}, {-sin(ang[2]), cos(ang[2]), 0}, {0, 0, 1}}));
	Matrix<float> res = rotate_x*rotate_y*rotate_z;
	
	for (size_t i = 0; i < vertex; i++) {
		Matrix<float> temp(vector<vector<float>>({{coords[i][0]-cent.x, coords[i][1]-cent.y, coords[i][2]-cent.z}}));
		Matrix<float> new_coords = temp*res;
		coords[i][0] = new_coords(0, 0)+cent.x; coords[i][1] = new_coords(0, 1)+cent.y; coords[i][2] = new_coords(0, 2)+cent.z;
	}
	delete ang;
}

void scale(vector<array<float, 3>> &coords, float coef) {
	centroid cent = get_centroid(coords);
	coef = coef > 0 ? coef : -1 / coef;
	Matrix<float> scale_mat(vector<vector<float>>({{coef, 0, 0}, {0, coef, 0}, {0, 0, coef}}));

	for (size_t i = 0; i < vertex; i++) {
		Matrix<float> temp(vector<vector<float>>({{coords[i][0]-cent.x, coords[i][1]-cent.y, coords[i][2]-cent.z}}));
		Matrix<float> new_coords = temp*scale_mat;
		coords[i][0] = new_coords(0, 0)+cent.x; coords[i][1] = new_coords(0, 1)+cent.y; coords[i][2] = new_coords(0, 2)+cent.z;
	}
}

void translate(vector<array<float, 3>> &coords, float dx, float dy, float dz) { // можно переписать с offset[3]
	for (size_t i = 0; i < vertex; i++) coords[i][0] += dx, coords[i][1] += dy, coords[i][2] += dz;
}

bool intersection(float line1[2][2], float line2[2][2], float *i_x, float *i_y)  {
	float c1_x, c1_y, c2_x, c2_y;
	c1_x = line1[1][0] - line1[0][0];
	c1_y = line1[1][1] - line1[0][1];
	c2_x = line2[1][0] - line2[0][0];
	c2_y = line2[1][1] - line2[0][1];

	float s, t;
	s = (-c1_y * (line1[0][0] - line2[0][0]) + c1_x * (line1[0][1] - line2[0][1])) / (-c2_x * c1_y + c1_x * c2_y);
	t = ( c2_x * (line1[0][1] - line2[0][1]) - c2_y * (line1[0][0] - line2[0][0])) / (-c2_x * c1_y + c1_x * c2_y);
	
	if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
		if (i_x != NULL)
			*i_x = line1[0][0] + (t * c1_x);
		if (i_y != NULL)
			*i_y = line1[0][1] + (t * c1_y);
		return 1;
	}

	return 0;
}

void fill(const vector<array<float, 3>> &coords) {
	setcolor((color++) % 10);
	int count_vertex = coords.size();
	float bounds[2][2] = {{coords[0][0], coords[0][1]}, {coords[0][0], coords[0][1]}};
	for (int i = 1; i < count_vertex; i++) {
		int x = coords[i][0], y = coords[i][1];
		if (x < bounds[0][0]) bounds[0][0] = x;
		else if (x > bounds[1][0]) bounds[1][0] = x;
		if (y < bounds[0][1]) bounds[0][1] = y;
		else if (y > bounds[1][1]) bounds[1][1] = y;
	}
	for (float i = bounds[0][1]+ 0.01; i <= bounds[1][1]; i++) {
		float pts[20] = {0};
		int count = 0;

		float pt_vec[2][2] = {{bounds[0][0] - 1, i}, {bounds[1][0] + 1, i}};
		for (int j = 0; j < count_vertex; j++)
		{
			float line[2][2] = {{coords[j][0], coords[j][1]},
													{coords[(j + 1) % count_vertex][0], coords[(j + 1) % count_vertex][1]}};
			float x, y;
			if (intersection(pt_vec, line, &x, &y)) pts[count++] = x;
		}
		if (count) {
			qsort(pts, count, sizeof(float), [](const void* x, const void* y) { return (int)(*(float*)x - *(float*)y); });
			for (int j = 0; j < count; j++) {
				if (pts[j+1]) line(pts[j], i, pts[j + 1], i);
			}
		}
	}
}

void sort_sides(vector<array<float, 3>> &coords) {
	int side_count = sides.size();
	double min_z[side_count];

	for (size_t i = 0; i < side_count; i++) {
		min_z[i] = numeric_limits<double>::max();
		for (size_t j = 0; j < sides[i].size(); j++) {
			double z = coords[sides[i][j]][2];
			min_z[i] = z < min_z[i] ? z : min_z[i];
		}
	}
	
	for (size_t i = 0; i < side_count; i++)
		for(size_t j = side_count - 1; j > i; j--)
			if (min_z[j - 1] > min_z[j]) {
				swap(min_z[j - 1], min_z[j]);
				swap(sides[j - 1], sides[j]);
				swap(colors[j - 1], colors[j]);
			}
}

float border_line[2][2] = {{30, 30}, {700, 30}};

void show(vector<array<float, 3>> &coords) {
	vector<array<float, 3>> proj_coords(vertex);
	proj_coords = coords;

	// float fNear = 0.1f;
	// float fFar = 500.0f;
	// float fFov = 80.0; // угол обзора
	// float fAspectRatio = (float)width / (float)height;
	// float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159f);

	// Matrix<float> matProj(vector<vector<float>>({
	// 		{fAspectRatio * fFovRad, 0, 0, 0}, 
	// 		{0, fFovRad, 0, 0},
	// 		{0, 0, fFar / (fFar - fNear), 1.0f},
	// 		{0, 0, (-fFar * fNear) / (fFar - fNear), 0.0f}
	// }));
	
	// for (size_t i = 0; i < vertex; i++) {
	// 	Matrix<float> temp(vector<vector<float>>({{coords[i][0], coords[i][1], coords[i][2], 1}}));
	// 	Matrix<float> new_coords = temp*matProj;
	// 	proj_coords[i][0] = new_coords(0, 0); proj_coords[i][1] = new_coords(0, 1); proj_coords[i][2] = new_coords(0, 2);
	// 	cout << proj_coords[i][0] << " x " << coords[i][0] << "\n";
	// 	cout << proj_coords[i][1] << " y " << coords[i][1] << "\n";
	// 	cout << proj_coords[i][2] << " z " << coords[i][2] << "\n";
	// }

	//////////////
	// for (int i = 0; i < vertex; i++) {
	// 	float x = coords[i][0], y = coords[i][1], z = coords[i][2];
	// 	x = n * x / (-z);
	// 	y = n * y / (-z);
	// 	z = (a * z + b) / -z;
	// 	proj_coords[i][0] = x, proj_coords[i][1] = y, proj_coords[i][2] = z;
	// }
	//////////////

	// sort_sides(coords);
	float x, y;
	for (size_t i = 0; i < sides.size(); i++) {
		int count_vertex = sides[i].size();
		vector<array<float, 3>> plane(count_vertex);

		for (size_t j = 0; j < count_vertex; j++) {
			plane[j] = proj_coords[sides[i][j]];
			plane[j][0] = (plane[j][0] * 0.5 + 0.5) * width;
      plane[j][1] = (plane[j][1] * 0.5 + 0.5) * height;
		}
		
		// fill(plane);
		setcolor(WHITE);
		setlinestyle(SOLID_LINE, 1, 3);
		for (size_t j = 0; j < count_vertex; j++) {
    	line(plane[j][0], plane[j][1], plane[(j+1)%count_vertex][0], plane[(j+1)%count_vertex][1]);
		}

		
    // }
		// setcolor(RED);
		// line(border_line[0][0], border_line[0][1], border_line[1][0], border_line[1][1]);
		// for (size_t j = 0; j < count_vertex; j++) {
		// 	float line1[2][2] = {{plane[j][0], plane[j][1]}, {plane[(j+1)%count_vertex][0], plane[(j+1)%count_vertex][1]}};
		// 	if (intersection(line1, border_line, &x, &y)) {
		// 		if (line1[0][1] < border_line[0][1]) {
		// 			cout << x << " " << y << "\n";
		// 			line1[0][0] = x;
		// 			line1[0][1] = y;
		// 		} else {
		// 			line1[1][0] = x;
		// 			line1[1][1] = y;
		// 		}
		// 	}
		// 	setcolor(WHITE);
		// 	line(line1[0][0], line1[0][1], line1[1][0], line1[1][1]);	
		// }
		
		// int poly[8] = {
		// 	(int)round(triangle[0][0]), (int)round(triangle[0][1]),
		// 	(int)round(triangle[1][0]), (int)round(triangle[1][1]),
		// 	(int)round(triangle[2][0]), (int)round(triangle[2][1]),
		// 	(int)round(triangle[0][0]), (int)round(triangle[0][1]),
		// };
		// setfillstyle(SOLID_FILL, i);
		// fillpoly(4, poly);
	}
	vector<array<float, 3>> poly;
		
		// float ground = -1.5;
		// float light_dir[3] = {0.0, -1.0, 0.0};

		// // for (int i = 0; i < vertex; i++) {
			
		// for(int j = 0; j < coords.size(); j++){
		// 		array<float, 3> pt = coords[j];

		// 		float t = (ground - pt[1]) / light_dir[1];

		// 		pt[0] += t*light_dir[0];
		// 		pt[1] = ground; 
		// 		pt[2] += t*light_dir[2];   
		// 		poly.push_back(pt);
		// }
		
		// poly.color = RGB(0, 0, 0);
		setcolor(RED);
		for (size_t j = 0; j < poly.size(); j++) {
			cout << poly[j][2] << "\n"; 
			line(poly[j][0], poly[j][1], poly[(j+1)%vertex][0], poly[(j+1)%vertex][1]);
		}
	color = 1;
}

int main() {
	int win = initwindow(width, height, "ddd"), key;

	// vector<array<float, 3>> prism_1 = {{100, 100, 10}, {500, 100, 10}, {300, 100, 200}, {100, 400, 10}, {500, 400, 10}, {300, 400, 200}};
	vector<array<float, 3>> prism_1 = {{-0.5, -0.5, 2.0},
																			{0.5, -0.5, 2.0},
																			{0.0, -0.5, 2.8},
																			{-0.5,  0.5, 2.0},
																			{0.5,  0.5, 2.0},
																			{0.0,  0.5, 2.8}};
	// vector<array<float, 3>> prism_2 = {{100+300, 100, 0}, {500+300, 100, 0}, {300+300, 100, 200}, {100+300, 400, 0}, {500+300, 400, 0}, {300+300, 400, 200}};

	// show(prism_1);

	// rotate(prism_1, new float[3]{deg2rad(-angle*5), deg2rad(0), deg2rad(0)});
  setbkcolor(RGB(55, 55, 55));
	cleardevice();
	show(prism_1);

	while (1) {
		setlinestyle(SOLID_LINE, 1, 3);
		// setcolor(RED); line(300,50,300,350);
		// setcolor(GREEN); line(300,350,700,350);
		// setcolor(YELLOW); line(300,350,50,600);
		setlinestyle(SOLID_LINE, 1, 1);
		key = getch();
		if (key == S_KEY) rotate(prism_1, new float[3]{deg2rad(angle), deg2rad(0), deg2rad(0)});
		else if (key == W_KEY) rotate(prism_1, new float[3]{deg2rad(-angle), deg2rad(0), deg2rad(0)});
		else if (key == D_KEY) rotate(prism_1, new float[3]{deg2rad(0), deg2rad(angle), deg2rad(0)});
		else if (key == A_KEY) rotate(prism_1, new float[3]{deg2rad(0), deg2rad(-angle), deg2rad(0)});
		else if (key == X_KEY) rotate(prism_1, new float[3]{deg2rad(0), deg2rad(0), deg2rad(angle)});
		else if (key == Z_KEY) rotate(prism_1, new float[3]{deg2rad(0), deg2rad(0), deg2rad(-angle)});
		else if (key == PLUS_KEY) scale(prism_1, inc_coef);
		else if (key == MINUS_KEY) scale(prism_1, -inc_coef);
		else if (key == UP_KEY) translate(prism_1, 0, -speed, 0);
		else if (key == DOWN_KEY) translate(prism_1, 0, speed, 0);
		else if (key == LEFT_KEY) translate(prism_1, -speed, 0, 0);
		else if (key == RIGHT_KEY) translate(prism_1, speed, 0, 0);
		else if (key == 91) translate(prism_1, 0, 0, speed);
		else if (key == 93) translate(prism_1, 0, 0, -speed);
		show(prism_1); swapbuffers(); cleardevice(); 
		// Sleep(delay);
	}
	return 0;
}