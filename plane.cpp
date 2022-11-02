#define _USE_MATH_DEFINES
#include "plane.h"

int colors[] = {CYAN, RED, YELLOW, BLUE, BROWN, GREEN, RED, BLUE, GREEN};
int color = 1;
bool is_1 = true;

struct centroid { float x, y, z; };

centroid get_centroid(const vector<array<float, 3>> coords) {
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

void translate(vector<array<float, 3>> &coords, float dx, float dy, float dz) {
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
		if (i_x != NULL) *i_x = line1[0][0] + (t * c1_x);
		if (i_y != NULL) *i_y = line1[0][1] + (t * c1_y);
		return 1;
	}

	return 0;
}

void fill(const vector<array<float, 3>> coords) {
	setcolor((color++) % 10);
	float bounds[2][2] = {{coords[0][0], coords[0][1]}, {coords[0][0], coords[0][1]}};
	for (int i = 1; i < vertex; i++) {
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
		for (int j = 0; j < vertex; j++)
		{
			float line[2][2] = {{coords[j][0], coords[j][1]},
													{coords[(j + 1) % vertex][0], coords[(j + 1) % vertex][1]}};
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

void show(const vector<array<float, 3>> coords) {
  for (size_t j = 0; j < vertex; j++) {
    line(coords[j][0], coords[j][1], coords[(j+1)%vertex][0], coords[(j+1)%vertex][1]);	
  }
  fill(coords);
  if (color > 2) color = 1;
}



int main() {
	int win = initwindow(width, height, "cg"), key;
	
	// float plane1[vertex][3] = {{100, 100, -100}, {200, 100, -100}, {200, 200, 0}, {50, 200, 0}};
  // float plane2[vertex][3] = {{300, 100, -100}, {400, 100, -100}, {400, 200, 0}, {250, 200, 0}};
	
  vector<array<float, 3>> plane1 = {{100, 100, 100}, {200, 100, -100}, {200, 200, 0}, {50, 200, 0}};
  vector<array<float, 3>> plane2 = {{300, 100, 100}, {400, 100, -100}, {400, 200, 0}, {250, 200, 0}};

  setbkcolor(RGB(55, 55, 55));
	cleardevice();
	show(plane1);
  show(plane2);
  
  float window[4][2][2] = {{{0, 0}, {width, 0}}, {{width, 0}, {width, height}}, {{width, height}, {0, height}}, {{0, height}, {0, 0}}};
  float x,y;
	while (1) {
		key = getch();

    // for (size_t j = 0; j < vertex; j++) {
    //   float line1[2][2] = {{plane1[j][0], plane1[j][1]}, {plane1[(j+1)%vertex][0], plane1[(j+1)%vertex][1]}};
    //   for (size_t i = 0; i < 4; i++) {
    //     cout << (intersection(line1, window[i], &x, &y) ? "true" : "") << "\n";
    //   }
    // }

    // for (size_t j = 0; j < vertex; j++) {
    //   cout << (plane1[j][0] >= 0 && plane1[j][0] <= width && plane1[j][1] >= 0 && plane1[j][1] <= height ? "in" : "") << "\n";
    // }

		if (key == S_KEY) rotate((is_1 ? plane1 : plane2), new float[3]{deg2rad(angle), deg2rad(0), deg2rad(0)});
		else if (key == W_KEY) rotate((is_1 ? plane1 : plane2), new float[3]{deg2rad(-angle), deg2rad(0), deg2rad(0)});
		else if (key == D_KEY) rotate((is_1 ? plane1 : plane2), new float[3]{deg2rad(0), deg2rad(angle), deg2rad(0)});
		else if (key == A_KEY) rotate((is_1 ? plane1 : plane2), new float[3]{deg2rad(0), deg2rad(-angle), deg2rad(0)});
		else if (key == X_KEY) rotate((is_1 ? plane1 : plane2), new float[3]{deg2rad(0), deg2rad(0), deg2rad(angle)});
		else if (key == Z_KEY) rotate((is_1 ? plane1 : plane2), new float[3]{deg2rad(0), deg2rad(0), deg2rad(-angle)});
		else if (key == PLUS_KEY) scale((is_1 ? plane1 : plane2), inc_coef);
		else if (key == MINUS_KEY) scale((is_1 ? plane1 : plane2), -inc_coef);
		else if (key == UP_KEY) translate((is_1 ? plane1 : plane2), 0, -speed, 0);
		else if (key == DOWN_KEY) translate((is_1 ? plane1 : plane2), 0, speed, 0);
		else if (key == LEFT_KEY) translate((is_1 ? plane1 : plane2), -speed, 0, 0);
		else if (key == RIGHT_KEY) translate((is_1 ? plane1 : plane2), speed, 0, 0);
    else if(key == ESC_KEY) is_1 = !is_1;
		show(plane1); show(plane2); swapbuffers(); cleardevice(); 
		Sleep(delay);
	}
	return 0;
}